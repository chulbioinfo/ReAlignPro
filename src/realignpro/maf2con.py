#!/usr/bin/env python3
"""
maf2con.py

Find constrained regions from MAF/MAF.GZ alignments and write BED3 intervals.

Constrained definition (per reference base):
  - All target species must exist in the block.
  - The most common A/C/G/T allele in the target group must be at least
    --min-major-similarity (default: 0.99).
  - Gap, N, and other ambiguous bases are not major-allele candidates, but they
    remain in the denominator so missing/ambiguous columns are penalized.

Typical usage:
  realignpro maf2con --input merged.maf.gz --ref-id hg38 --target-ids all
  realignpro maf2con --input merged.maf --output constrained.bed --ref-id hg38 --target-ids hg38,hs1,hs2
"""

from __future__ import annotations

import argparse
import multiprocessing as mp
import sys
import threading
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Dict, Iterable, List, Optional, Set, Tuple

from .maf2bed import (
    DEFAULT_OUT_QSIZE,
    DEFAULT_START_METHOD,
    DEFAULT_WORK_QSIZE,
    OUT_STOP_IDX,
    WORK_STOP,
    _split_csv_strict,
    _split_tokens,
    list_maf_ids,
    open_text_maybe_gzip,
    parse_maf_block_lines,
    writer_thread,
)


DNA_BASES = {"A", "C", "G", "T"}


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
        description="Find target-group constrained MAF positions by major-allele similarity and write BED3 intervals.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "-ids",
        "--ids",
        dest="ids_maf",
        metavar="MAF",
        help="List unique species/assembly IDs in the given MAF/MAF.GZ and exit.",
    )
    parser.add_argument(
        "-i",
        "--input",
        dest="input_maf",
        metavar="MAF",
        help="Input MAF file path (.maf or .maf.gz).",
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output_bed",
        metavar="BED",
        default=None,
        help="Output BED file path. If omitted, derived from input by replacing .maf(.gz) with .con.bed.",
    )
    parser.add_argument(
        "-t",
        "--threads",
        dest="total_threads",
        type=int,
        default=4,
        help="Total threads/processes INCLUDING 1 reader + 1 writer. Must be >= 3.",
    )
    parser.add_argument(
        "--ref-id",
        dest="ref_id",
        required=False,
        help="Reference species ID (used for coordinate system in BED output). REQUIRED in run mode.",
    )
    parser.add_argument(
        "--target-ids",
        dest="target_ids",
        metavar="all|ID1,ID2,...",
        required=False,
        help="Target species IDs as comma-separated values (no spaces), or 'all' for every ID in the MAF.",
    )
    parser.add_argument(
        "--outgroup-ids",
        dest="outgroup_ids",
        nargs="*",
        default=None,
        help="Optional species IDs to exclude from constrained-region calculation (space-separated and/or comma-separated).",
    )
    parser.add_argument(
        "--min-major-similarity",
        dest="min_major_similarity",
        type=float,
        default=0.99,
        help="Minimum major allele frequency within the target group.",
    )
    parser.add_argument(
        "--min-target-count",
        dest="min_target_count",
        type=int,
        default=2,
        help="Minimum number of target species required after outgroup exclusion.",
    )
    parser.add_argument(
        "--work-qsize",
        dest="work_qsize",
        type=int,
        default=DEFAULT_WORK_QSIZE,
        help="Max number of buffered work items (blocks) between reader and workers.",
    )
    parser.add_argument(
        "--out-qsize",
        dest="out_qsize",
        type=int,
        default=DEFAULT_OUT_QSIZE,
        help="Max number of buffered result items between workers and writer.",
    )
    parser.add_argument(
        "--start-method",
        dest="start_method",
        default=DEFAULT_START_METHOD,
        choices=["spawn", "fork", "forkserver"],
        help="Multiprocessing start method.",
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
    if not (0.0 <= min_major_similarity <= 1.0):
        raise ValueError("--min-major-similarity must be in [0, 1].")

    min_target_count = int(args.min_target_count)
    if min_target_count < 1:
        raise ValueError("--min-target-count must be >= 1.")

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
        work_qsize=work_qsize,
        out_qsize=out_qsize,
        start_method=args.start_method,
    )


def call_major_base(
    target_nts: Iterable[str],
    min_major_similarity: float = 0.99,
) -> Tuple[Optional[str], float]:
    """
    Return the unique major A/C/G/T base and its frequency when it meets the threshold.

    Gap, N, and other ambiguous values are excluded from candidate counts but included
    in the denominator.
    """
    nts = [nt.upper() for nt in target_nts]
    denom = len(nts)
    if denom == 0:
        return None, 0.0

    counts = Counter(nt for nt in nts if nt in DNA_BASES)
    if not counts:
        return None, 0.0

    top_count = max(counts.values())
    top_bases = [base for base, count in counts.items() if count == top_count]
    ratio = top_count / float(denom)

    if len(top_bases) != 1:
        return None, ratio
    if ratio < min_major_similarity:
        return None, ratio

    return top_bases[0], ratio


def matrix2con(
    block: Dict[str, Dict[str, Any]],
    ref_id: str,
    target_ids: List[str],
    outgroup_ids: List[str],
    min_major_similarity: float = 0.99,
    min_target_count: int = 2,
) -> List[str]:
    """
    Identify target-group constrained positions on the reference, then merge adjacent hits
    into BED3 intervals within the current MAF block.
    """
    outgroup = set(outgroup_ids)
    targets = [sid for sid in target_ids if sid not in outgroup]
    required = set(targets)

    if len(targets) < min_target_count:
        return []
    if not required.issubset(block.keys()):
        return []

    ref = block[ref_id]
    ref_seq = ref["seq"]
    ref_chr = ref["chr"]
    ref_start = ref["start"]
    ref_strand = ref["strand"]
    ref_src_size = ref["srcSize"]

    if ref_strand == "+":
        pos_ref = ref_start - 1
        step = +1
    else:
        forward_end = ref_src_size - ref_start
        pos_ref = forward_end
        step = -1

    merged: List[Tuple[str, int, int]] = []
    cur_start: Optional[int] = None
    cur_end: Optional[int] = None

    def flush_current() -> None:
        nonlocal cur_start, cur_end
        if cur_start is not None and cur_end is not None:
            merged.append((ref_chr, cur_start, cur_end))
        cur_start, cur_end = None, None

    aln_idx = -1
    for ref_nt in ref_seq:
        aln_idx += 1
        if ref_nt == "-":
            continue

        pos_ref += step

        target_nts: List[str] = []
        for sid in targets:
            sid_seq = block[sid]["seq"]
            nt = sid_seq[aln_idx].upper() if aln_idx < len(sid_seq) else "-"
            target_nts.append(nt)

        major_base, _ratio = call_major_base(target_nts, min_major_similarity)
        is_hit = major_base is not None

        if is_hit:
            base_start = pos_ref
            base_end = pos_ref + 1

            if cur_start is None:
                cur_start, cur_end = base_start, base_end
            else:
                if ref_strand == "+":
                    if base_start == cur_end:
                        cur_end = base_end
                    else:
                        flush_current()
                        cur_start, cur_end = base_start, base_end
                else:
                    if base_end == cur_start:
                        cur_start = base_start
                    else:
                        flush_current()
                        cur_start, cur_end = base_start, base_end
        else:
            flush_current()

    flush_current()
    return [f"{chrom}\t{start}\t{end}\n" for chrom, start, end in merged]


def worker_proc(
    work_q: "mp.JoinableQueue",
    out_q: "mp.Queue",
    ref_id: str,
    target_ids: List[str],
    outgroup_ids: List[str],
    min_major_similarity: float,
    min_target_count: int,
) -> None:
    """
    Worker process:
      - pull (block_index, block_dict)
      - compute merged constrained BED lines for that block
      - push (block_index, bed_list) to out_q
    """
    while True:
        item = work_q.get()
        try:
            if item is WORK_STOP:
                break
            idx, block = item
            bed_list = matrix2con(
                block,
                ref_id,
                target_ids,
                outgroup_ids,
                min_major_similarity=min_major_similarity,
                min_target_count=min_target_count,
            )
            out_q.put((idx, bed_list))
        finally:
            work_q.task_done()

    out_q.put((OUT_STOP_IDX, None))


def maf2con_multiprocessing(
    maf_path: str,
    bed_path: str,
    ref_id: str,
    target_ids: List[str],
    outgroup_ids: List[str],
    total_threads: int,
    work_qsize: int,
    out_qsize: int,
    min_major_similarity: float = 0.99,
    min_target_count: int = 2,
    start_method: str = DEFAULT_START_METHOD,
) -> None:
    """
    Orchestrate reader, worker processes, and ordered writer for maf2con.
    """
    n_workers = max(1, total_threads - 2)
    ctx = mp.get_context(start_method)

    work_q: "mp.JoinableQueue" = ctx.JoinableQueue(maxsize=work_qsize)
    out_q: "mp.Queue" = ctx.Queue(maxsize=out_qsize)

    procs = [
        ctx.Process(
            target=worker_proc,
            args=(
                work_q,
                out_q,
                ref_id,
                target_ids,
                outgroup_ids,
                min_major_similarity,
                min_target_count,
            ),
        )
        for _ in range(n_workers)
    ]
    for p in procs:
        p.start()

    t_writer = threading.Thread(target=writer_thread, args=(out_q, bed_path, n_workers), daemon=True)
    t_writer.start()

    idx = 0
    in_block = False
    block_lines: List[str] = []

    with open_text_maybe_gzip(maf_path) as f:
        for raw in f:
            line = raw.rstrip("\n")

            if not in_block:
                if line.startswith("a"):
                    in_block = True
                    block_lines = [line]
                continue

            if line.strip() == "":
                block = parse_maf_block_lines(block_lines, ref_id)
                if block is not None:
                    work_q.put((idx, block))
                    idx += 1
                in_block = False
                block_lines = []
                continue

            block_lines.append(line)

        if in_block and block_lines:
            block = parse_maf_block_lines(block_lines, ref_id)
            if block is not None:
                work_q.put((idx, block))
                idx += 1

    for _ in range(n_workers):
        work_q.put(WORK_STOP)

    work_q.join()

    for p in procs:
        p.join()

    t_writer.join()


def main(argv: Optional[List[str]] = None) -> int:
    """
    Entry point.

    Modes:
      1) ID listing mode:
         maf2con.py --ids input.maf
      2) Run mode:
         maf2con.py --input <maf> --ref-id <id> --target-ids <all|id...> [options]
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

    maf2con_multiprocessing(
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
        start_method=cfg.start_method,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
