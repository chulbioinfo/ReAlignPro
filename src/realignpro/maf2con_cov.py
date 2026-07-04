#!/usr/bin/env python3
"""
maf2con_cov.py — coverage-aware constrained-region engine for `realignpro maf2con`.

This is the engine behind `maf2con`. It calls a reference base as constrained using a
per-column, coverage-aware tally, which stays correct when the set of aligned assemblies
varies from block to block (e.g., sex chromosomes):

  1. No whole-block gate. A block is assessed on whatever target assemblies are
     actually present, instead of being skipped unless every target is present.
  2. Per-column coverage-aware denominator. The major-allele fraction is
     ``n_maj / N_cov`` where ``N_cov`` is the number of aligned A/C/G/T at that
     column (gap/N/missing excluded from both numerator and denominator).
  3. Per-column coverage floor. A column is only callable when
     ``N_cov >= min_depth`` and, when a relative floor is in effect,
     ``N_cov >= ceil(min_call_rate * N_exp)``, where ``N_exp`` is the expected
     aligned depth of the reference chromosome. ``N_exp`` is derived per reference
     chromosome by a cheap pre-scan when the relative floor is used and
     ``expected_depth`` is not given, so autosomes, chrX and chrY are handled in one run.

``fixed_only`` restricts hits to 100%-fixed (monomorphic) columns.

Homologous cross-alignment between X and Y (PAR / XTR / gametolog regions) is not
separated here; pooling homologous lineages only dilutes the major-allele fraction,
so its dominant effect is a conservative false negative. Use ``emit_depth`` to keep
N_cov / major-allele frequency per interval for downstream filtering. See
ConstraintVariantAnalysisPlan_v0.3.md for the rationale.
"""

from __future__ import annotations

import math
import multiprocessing as mp
import threading
from collections import Counter
from typing import Any, Dict, List, Optional, Tuple

from .maf2bed import (
    OUT_STOP_IDX,
    WORK_STOP,
    open_text_maybe_gzip,
    parse_maf_block_lines,
    writer_thread,
)

DNA_BASES = {"A", "C", "G", "T"}


def call_major_base_cov(
    target_nts,
    min_major_similarity: float = 0.99,
    min_depth: int = 1,
    fixed_only: bool = False,
) -> Tuple[Optional[str], float, int]:
    """
    Coverage-aware major-allele test for one column.

    Denominator is ``N_cov`` = number of aligned A/C/G/T (gap/N/other excluded from
    BOTH numerator and denominator). Returns ``(major_base, maj_freq, n_cov)`` when
    a unique major allele strictly exceeds ``min_major_similarity`` and
    ``N_cov >= min_depth``; otherwise ``(None, maj_freq_or_0, n_cov)``.

    When ``fixed_only`` is set, the column is a hit only if it is 100% fixed, i.e.
    a single distinct A/C/G/T allele across ALL aligned bases (``maj_freq == 1.0``);
    ``min_major_similarity`` is then ignored.
    """
    counts = Counter(nt.upper() for nt in target_nts if nt.upper() in DNA_BASES)
    n_cov = sum(counts.values())
    if n_cov < min_depth or not counts:
        return None, 0.0, n_cov

    top_count = max(counts.values())
    top_bases = [base for base, count in counts.items() if count == top_count]
    ratio = top_count / float(n_cov)

    if len(top_bases) != 1:
        return None, ratio, n_cov
    if fixed_only:
        # 100% fixed: exactly one distinct aligned allele (no mismatch at all).
        if len(counts) != 1:
            return None, ratio, n_cov
        return top_bases[0], ratio, n_cov
    if ratio <= min_major_similarity:
        return None, ratio, n_cov
    return top_bases[0], ratio, n_cov


def matrix2con_cov(
    block: Dict[str, Dict[str, Any]],
    ref_id: str,
    target_ids: List[str],
    outgroup_ids: List[str],
    min_major_similarity: float = 0.99,
    min_target_count: int = 2,
    min_depth: int = 2,
    min_call_rate: float = 0.0,
    expected_depth_map: Optional[Dict[str, int]] = None,
    expected_depth_global: Optional[int] = None,
    emit_depth: bool = False,
    fixed_only: bool = False,
) -> List[str]:
    """
    Coverage-aware constrained positions on the reference for one MAF block,
    merged into BED intervals. No whole-block gate: the tally uses whichever
    targets are present in this block. When ``fixed_only`` is set, only 100%
    fixed (monomorphic) columns are reported.
    """
    if ref_id not in block:
        return []
    outgroup = set(outgroup_ids)
    # Assessed set = targets that are actually present here (order preserved).
    targets = [sid for sid in target_ids if sid in block and sid not in outgroup]
    if len(targets) < min_target_count:
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
        pos_ref = ref_src_size - ref_start
        step = -1

    # Resolve the expected-depth call-rate floor for this reference chromosome.
    exp = None
    if expected_depth_map is not None:
        exp = expected_depth_map.get(ref_chr)
    if exp is None:
        exp = expected_depth_global
    # ceil so that "call rate >= r" is exact: N_cov >= ceil(r * N_exp)
    # (e.g. r=0.99, N_exp=118 -> require N_cov >= 117, i.e. <= 1% missing).
    call_rate_floor = math.ceil(min_call_rate * exp) if (exp and min_call_rate > 0.0) else 0
    depth_floor = max(min_depth, call_rate_floor)

    seqs = [block[sid]["seq"] for sid in targets]

    merged: List[Tuple[str, int, int, int, float]] = []
    cur: Optional[Dict[str, Any]] = None

    def flush() -> None:
        nonlocal cur
        if cur is not None:
            merged.append((ref_chr, cur["s"], cur["e"], cur["n"], cur["f"]))
        cur = None

    aln_idx = -1
    for ref_nt in ref_seq:
        aln_idx += 1
        if ref_nt == "-":
            continue
        pos_ref += step

        col = [s[aln_idx] for s in seqs if aln_idx < len(s)]
        base, ratio, n_cov = call_major_base_cov(col, min_major_similarity, depth_floor, fixed_only)
        is_hit = base is not None

        if is_hit:
            base_start = pos_ref
            base_end = pos_ref + 1
            if cur is None:
                cur = {"s": base_start, "e": base_end, "n": n_cov, "f": ratio}
            else:
                if ref_strand == "+":
                    contiguous = base_start == cur["e"]
                else:
                    contiguous = base_end == cur["s"]
                if contiguous:
                    if ref_strand == "+":
                        cur["e"] = base_end
                    else:
                        cur["s"] = base_start
                    cur["n"] = min(cur["n"], n_cov)
                    cur["f"] = min(cur["f"], ratio)
                else:
                    flush()
                    cur = {"s": base_start, "e": base_end, "n": n_cov, "f": ratio}
        else:
            flush()

    flush()

    if emit_depth:
        return [f"{c}\t{s}\t{e}\t{n}\t{f:.4f}\n" for (c, s, e, n, f) in merged]
    return [f"{c}\t{s}\t{e}\n" for (c, s, e, _n, _f) in merged]


def scan_expected_depth(
    maf_path: str,
    ref_id: str,
    target_set: set,
    outgroup_set: set,
) -> Dict[str, int]:
    """
    Cheap single pass over the MAF: for each reference chromosome, return the
    maximum per-block assessed depth (number of aligned target assemblies).
    Used as ``N_exp(chrom)`` for the optional call-rate floor. No sequence work.
    """
    n_exp: Dict[str, int] = {}
    cur_chr: Optional[str] = None
    cur_depth = 0

    def commit() -> None:
        nonlocal cur_chr, cur_depth
        if cur_chr is not None and cur_depth > n_exp.get(cur_chr, 0):
            n_exp[cur_chr] = cur_depth
        cur_chr = None
        cur_depth = 0

    with open_text_maybe_gzip(maf_path) as f:
        for raw in f:
            if not raw:
                continue
            c0 = raw[0]
            if c0 == "a":
                commit()
                continue
            if c0 == "s":
                parts = raw.split()
                if len(parts) < 2:
                    continue
                src = parts[1]
                if "." in src:
                    sid, chrom = src.split(".", 1)
                else:
                    sid, chrom = src, src
                if sid == ref_id:
                    cur_chr = chrom
                if sid in target_set and sid not in outgroup_set:
                    cur_depth += 1
            elif raw.strip() == "":
                commit()
        commit()

    return n_exp


def _worker_proc_cov(
    work_q: "mp.JoinableQueue",
    out_q: "mp.Queue",
    ref_id: str,
    target_ids: List[str],
    outgroup_ids: List[str],
    params: Dict[str, Any],
) -> None:
    while True:
        item = work_q.get()
        try:
            if item is WORK_STOP:
                break
            idx, block = item
            bed_list = matrix2con_cov(
                block,
                ref_id,
                target_ids,
                outgroup_ids,
                min_major_similarity=params["min_major_similarity"],
                min_target_count=params["min_target_count"],
                min_depth=params["min_depth"],
                min_call_rate=params["min_call_rate"],
                expected_depth_map=params["expected_depth_map"],
                expected_depth_global=params["expected_depth_global"],
                emit_depth=params["emit_depth"],
                fixed_only=params["fixed_only"],
            )
            out_q.put((idx, bed_list))
        finally:
            work_q.task_done()

    out_q.put((OUT_STOP_IDX, None))


def maf2con_cov_multiprocessing(
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
    min_depth: int = 2,
    min_call_rate: float = 0.0,
    expected_depth: Optional[int] = None,
    auto_depth: bool = False,
    emit_depth: bool = False,
    fixed_only: bool = False,
    start_method: str = "spawn",
) -> None:
    """
    Reader + coverage-aware workers + ordered writer. When ``auto_depth`` is set,
    a cheap pre-scan derives ``N_exp`` per reference chromosome first.
    """
    expected_depth_map: Optional[Dict[str, int]] = None
    if auto_depth:
        expected_depth_map = scan_expected_depth(
            maf_path, ref_id, set(target_ids), set(outgroup_ids)
        )

    params: Dict[str, Any] = {
        "min_major_similarity": min_major_similarity,
        "min_target_count": min_target_count,
        "min_depth": min_depth,
        "min_call_rate": min_call_rate,
        "expected_depth_map": expected_depth_map,
        "expected_depth_global": expected_depth,
        "emit_depth": emit_depth,
        "fixed_only": fixed_only,
    }

    n_workers = max(1, total_threads - 2)
    ctx = mp.get_context(start_method)

    work_q: "mp.JoinableQueue" = ctx.JoinableQueue(maxsize=work_qsize)
    out_q: "mp.Queue" = ctx.Queue(maxsize=out_qsize)

    procs = [
        ctx.Process(
            target=_worker_proc_cov,
            args=(work_q, out_q, ref_id, target_ids, outgroup_ids, params),
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
