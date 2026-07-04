#!/usr/bin/env python3
"""
maf2con.py

Find constrained regions from MAF/MAF.GZ alignments and write BED intervals.

Coverage-aware definition (per reference base):
  - N_cov = number of aligned A/C/G/T among the target assemblies at that column
    (gap/N/missing are excluded from BOTH numerator and denominator).
  - The column is callable only when it clears the coverage floor:
      N_cov >= --min-depth
      AND (when a call-rate floor applies) N_cov >= ceil(--min-call-rate * N_exp),
      where N_exp is the expected aligned depth of the reference chromosome.
  - It is "constrained" when the unique major A/C/G/T allele fraction (n_maj / N_cov)
    is strictly greater than --min-major-similarity (default 0.99); with --fixed-only
    only 100% monomorphic columns (zero mismatch) are reported.

N_exp is derived per reference chromosome automatically by a cheap pre-scan when a
call-rate floor is in effect (--min-call-rate > 0) and --expected-depth is not given,
so autosomes, chrX and chrY are handled correctly in a single run. The engine lives
in realignpro/maf2con_cov.py; see ConstraintVariantAnalysisPlan_v0.3.md for rationale.

Typical usage:
  realignpro maf2con --input merged.maf.gz --ref-id hg38 --target-ids all
  realignpro maf2con --input merged.maf.gz --ref-id hg38 --target-ids all --fixed-only
  realignpro maf2con --input merged.maf.gz --ref-id hg38 --target-ids all \
      --expected-depth 464 --emit-depth
  # disable the relative coverage floor (keep only the absolute --min-depth):
  realignpro maf2con --input merged.maf.gz --ref-id hg38 --target-ids all --min-call-rate 0
"""

from __future__ import annotations

import argparse
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, List, Optional, Set

from .maf2bed import (
    DEFAULT_OUT_QSIZE,
    DEFAULT_START_METHOD,
    DEFAULT_WORK_QSIZE,
    _split_csv_strict,
    _split_tokens,
    list_maf_ids,
)
from .maf2con_cov import maf2con_cov_multiprocessing


@dataclass(frozen=True)
class Maf2ConConfig:
    input_maf: str
    output_bed: str
    total_threads: int
    ref_id: str
    target_ids: List[str]
    outgroup_ids: List[str]
    min_major_similarity: float = 0.99
    min_target_count: int = 2
    min_depth: int = 2
    min_call_rate: float = 0.99
    expected_depth: Optional[int] = None
    fixed_only: bool = False
    emit_depth: bool = False
    work_qsize: int = DEFAULT_WORK_QSIZE
    out_qsize: int = DEFAULT_OUT_QSIZE
    start_method: str = DEFAULT_START_METHOD


def _derive_default_output_bed(input_maf: str) -> str:
    """
    Derive output BED path from input MAF path.
      foo.maf.gz -> foo.con.bed
      foo.maf    -> foo.con.bed
      otherwise  -> foo.con.bed
    """
    p = Path(input_maf)
    name = p.name

    if name.endswith(".maf.gz"):
        out_name = name[:-7] + ".con.bed"
    elif name.endswith(".maf"):
        out_name = name[:-4] + ".con.bed"
    else:
        out_name = name + ".con.bed"

    return str(p.with_name(out_name))


def _dedupe_preserve_order(values: Iterable[str]) -> List[str]:
    seen: Set[str] = set()
    out: List[str] = []
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return out


def _resolve_target_ids(input_maf: str, target_ids_arg: str, outgroup_ids: List[str]) -> List[str]:
    raw = (target_ids_arg or "").strip()
    if not raw:
        raise ValueError("--target-ids is required in run mode.")

    outgroup = set(outgroup_ids)
    if raw.lower() == "all":
        parsed = list_maf_ids(input_maf)
    else:
        parsed = _split_csv_strict(raw, "--target-ids")

    target_ids = [sid for sid in _dedupe_preserve_order(parsed) if sid not in outgroup]
    if not target_ids:
        raise ValueError("--target-ids is empty after parsing and outgroup exclusion.")
    return target_ids


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="realignpro maf2con",
        description="Coverage-aware constrained-region caller: per reference column, keep "
        "positions where the aligned assemblies clear a coverage floor and their major "
        "allele fraction exceeds a threshold (or are 100% fixed). Writes BED intervals.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-ids", "--ids", dest="ids_maf", metavar="MAF",
        help="List unique species/assembly IDs in the given MAF/MAF.GZ and exit.",
    )
    parser.add_argument(
        "-i", "--input", dest="input_maf", metavar="MAF",
        help="Input MAF file path (.maf or .maf.gz).",
    )
    parser.add_argument(
        "-o", "--output", dest="output_bed", metavar="BED", default=None,
        help="Output BED path. If omitted, derived from input by replacing .maf(.gz) with .con.bed.",
    )
    parser.add_argument(
        "-t", "--threads", dest="total_threads", type=int, default=4,
        help="Total threads/processes INCLUDING 1 reader + 1 writer. Must be >= 3.",
    )
    parser.add_argument(
        "--ref-id", dest="ref_id", required=False,
        help="Reference species ID (coordinate system for BED output). REQUIRED in run mode.",
    )
    parser.add_argument(
        "--target-ids", dest="target_ids", metavar="all|ID1,ID2,...", required=False,
        help="Target species IDs as comma-separated values (no spaces), or 'all' for every ID.",
    )
    parser.add_argument(
        "--outgroup-ids", dest="outgroup_ids", nargs="*", default=None,
        help="Optional species IDs to exclude from the calculation (space- and/or comma-separated).",
    )

    grp_con = parser.add_argument_group("constraint threshold")
    grp_con.add_argument(
        "--min-major-similarity", dest="min_major_similarity", type=float, default=0.99,
        help="Strict lower bound (>) for the major-allele fraction among aligned bases.",
    )
    grp_con.add_argument(
        "--fixed-only", dest="fixed_only", action="store_true",
        help="Strictest mode: report only 100%% fixed (monomorphic) columns -- a single distinct "
        "aligned allele, zero mismatch. Ignores --min-major-similarity.",
    )
    grp_con.add_argument(
        "--min-target-count", dest="min_target_count", type=int, default=2,
        help="Minimum number of present target assemblies required to assess a block.",
    )

    grp_cov = parser.add_argument_group("coverage floor")
    grp_cov.add_argument(
        "--min-depth", dest="min_depth", type=int, default=2,
        help="Absolute minimum aligned A/C/G/T count (N_cov) required to call a column.",
    )
    grp_cov.add_argument(
        "--min-call-rate", dest="min_call_rate", type=float, default=0.99,
        help="Relative coverage floor: require N_cov >= ceil(rate * N_exp) per column, where "
        "N_exp is the reference chromosome's expected depth (auto-derived, or --expected-depth). "
        "Set 0 to disable and use only --min-depth.",
    )
    grp_cov.add_argument(
        "--expected-depth", dest="expected_depth", type=int, default=None,
        help="Expected aligned depth (single global N_exp). If given, the per-chromosome pre-scan "
        "is skipped (saves one read pass). Otherwise N_exp is auto-derived per reference chromosome.",
    )

    grp_out = parser.add_argument_group("output / performance")
    grp_out.add_argument(
        "--emit-depth", dest="emit_depth", action="store_true",
        help="Append N_cov and major-allele frequency columns to each BED interval.",
    )
    grp_out.add_argument(
        "--work-qsize", dest="work_qsize", type=int, default=DEFAULT_WORK_QSIZE,
        help="Max number of buffered work items (blocks) between reader and workers.",
    )
    grp_out.add_argument(
        "--out-qsize", dest="out_qsize", type=int, default=DEFAULT_OUT_QSIZE,
        help="Max number of buffered result items between workers and writer.",
    )
    grp_out.add_argument(
        "--start-method", dest="start_method", default=DEFAULT_START_METHOD,
        choices=["spawn", "fork", "forkserver"], help="Multiprocessing start method.",
    )

    return parser


def _validate_and_build_config(args: argparse.Namespace) -> Maf2ConConfig:
    if not args.input_maf:
        raise ValueError("--input is required (or use --ids).")
    input_maf = args.input_maf.strip()
    if not input_maf:
        raise ValueError("--input is empty.")

    output_bed = args.output_bed.strip() if args.output_bed else _derive_default_output_bed(input_maf)

    total_threads = int(args.total_threads)
    if total_threads < 3:
        raise ValueError("--threads must be >= 3.")

    if not args.ref_id or not args.ref_id.strip():
        raise ValueError("--ref-id is required in run mode.")
    ref_id = args.ref_id.strip()

    outgroup_ids = _split_tokens(args.outgroup_ids or [])
    target_ids = _resolve_target_ids(input_maf, args.target_ids or "", outgroup_ids)

    min_major_similarity = float(args.min_major_similarity)
    if not (0.0 <= min_major_similarity < 1.0):
        raise ValueError("--min-major-similarity must be in [0, 1) because comparison is strict greater-than.")

    min_target_count = int(args.min_target_count)
    if min_target_count < 1:
        raise ValueError("--min-target-count must be >= 1.")

    min_depth = int(args.min_depth)
    if min_depth < 1:
        raise ValueError("--min-depth must be >= 1.")

    min_call_rate = float(args.min_call_rate)
    if not (0.0 <= min_call_rate <= 1.0):
        raise ValueError("--min-call-rate must be in [0, 1].")

    expected_depth = int(args.expected_depth) if args.expected_depth is not None else None
    if expected_depth is not None and expected_depth < 1:
        raise ValueError("--expected-depth must be >= 1.")

    work_qsize = int(args.work_qsize)
    out_qsize = int(args.out_qsize)
    if work_qsize < 1:
        raise ValueError("--work-qsize must be >= 1.")
    if out_qsize < 1:
        raise ValueError("--out-qsize must be >= 1.")

    return Maf2ConConfig(
        input_maf=input_maf,
        output_bed=output_bed,
        total_threads=total_threads,
        ref_id=ref_id,
        target_ids=target_ids,
        outgroup_ids=outgroup_ids,
        min_major_similarity=min_major_similarity,
        min_target_count=min_target_count,
        min_depth=min_depth,
        min_call_rate=min_call_rate,
        expected_depth=expected_depth,
        fixed_only=bool(args.fixed_only),
        emit_depth=bool(args.emit_depth),
        work_qsize=work_qsize,
        out_qsize=out_qsize,
        start_method=args.start_method,
    )


def main(argv: Optional[List[str]] = None) -> int:
    """
    Entry point.

    Modes:
      1) ID listing:  maf2con --ids input.maf
      2) Run:         maf2con --input <maf> --ref-id <id> --target-ids <all|id...> [options]
    """
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    if args.ids_maf:
        try:
            ids = list_maf_ids(args.ids_maf)
        except Exception as e:
            print(f"[ERROR] failed to read MAF for IDs: {e}", file=sys.stderr)
            return 1
        print(", ".join(ids))
        return 0

    try:
        cfg = _validate_and_build_config(args)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        parser.print_help(sys.stderr)
        return 1

    # Derive N_exp per reference chromosome only when a relative floor is in effect
    # and the user did not pin it explicitly.
    auto_depth = cfg.min_call_rate > 0.0 and cfg.expected_depth is None

    maf2con_cov_multiprocessing(
        maf_path=cfg.input_maf,
        bed_path=cfg.output_bed,
        ref_id=cfg.ref_id,
        target_ids=cfg.target_ids,
        outgroup_ids=cfg.outgroup_ids,
        total_threads=cfg.total_threads,
        work_qsize=cfg.work_qsize,
        out_qsize=cfg.out_qsize,
        min_major_similarity=cfg.min_major_similarity,
        min_target_count=cfg.min_target_count,
        min_depth=cfg.min_depth,
        min_call_rate=cfg.min_call_rate,
        expected_depth=cfg.expected_depth,
        auto_depth=auto_depth,
        emit_depth=cfg.emit_depth,
        fixed_only=cfg.fixed_only,
        start_method=cfg.start_method,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
