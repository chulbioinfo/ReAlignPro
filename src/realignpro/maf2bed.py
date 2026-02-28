#!/usr/bin/env python3
"""
maf2bed.py (formerly maf2var.py)

Extract variant-like positions from MAF/MAF.GZ and write BED3 intervals.

This version replaces the CTL control-file mode with command-line options.

Extra utility:
  -ids / --ids <input.maf[.gz]>
    Print unique assembly/species IDs found in the MAF file.
    IDs are extracted from 's' lines: src field before the first dot (.)
    Output format: id1, id2, id3

Variant definition (per reference base):
  - All target species must exist in the block
  - Target species share exactly one allele at that column
  - That allele must NOT appear among other species in the block
  - Species in OUTGROUP_IDS are skipped entirely when building allele sets

Output:
  - BED3 lines (chrom, start, end), UCSC 0-based half-open
  - Adjacent hits are merged into blocks within each MAF block to reduce output size

Typical usage:
  realignpro maf2bed --input merged.maf.gz --ref-id hg38 --target-ids hg38,m11_Homo_sapiens
  python3 maf2bed.py --input merged.maf --output custom.bed --threads 8 --ref-id hg38 --target-ids hg38,m11_Homo_sapiens --outgroup-ids mm10

Notes:
  - --input is required
  - --output is optional; if omitted, it is derived from --input by replacing:
      *.maf.gz -> *.bed
      *.maf    -> *.bed
    otherwise, appends ".bed"
  - --threads defaults to 4 (must be >= 3). total_threads includes 1 reader + 1 writer.
  - --ref-id and --target-ids are required (target-ids must be comma-separated)
  - --outgroup-ids is optional
  - --work-qsize and --out-qsize default to 200
"""

from __future__ import annotations

import argparse
import gzip
import sys
import threading
import multiprocessing as mp
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple, Set, Iterable


# Default multiprocessing start method
# - "spawn" is safest across platforms
# - on Linux HPC, "fork" can be faster (use with care if other threaded libs are involved)
DEFAULT_START_METHOD = "spawn"

# Queue defaults
DEFAULT_WORK_QSIZE = 200
DEFAULT_OUT_QSIZE = 200

# Sentinel values
WORK_STOP = None
OUT_STOP_IDX = -1


# ----------------------------
# CLI config
# ----------------------------
@dataclass(frozen=True)
class Maf2BedConfig:
    input_maf: str
    output_bed: str
    total_threads: int
    ref_id: str
    target_ids: List[str]
    outgroup_ids: List[str]
    work_qsize: int = DEFAULT_WORK_QSIZE
    out_qsize: int = DEFAULT_OUT_QSIZE
    start_method: str = DEFAULT_START_METHOD


def _split_csv_strict(value: str, arg_name: str = "--target-ids") -> List[str]:
    """Split a comma-separated ID list. Spaces are not allowed as separators.

    Valid:
      --target-ids hg38,m11_Homo_sapiens

    Invalid:
      --target-ids hg38 m11_Homo_sapiens
    """
    v = (value or "").strip()
    if not v:
        return []
    # Disallow whitespace to avoid silent mis-parsing.
    if re.search(r"\s", v):
        raise ValueError(
            f"{arg_name} must be comma-separated (no spaces). Example: hg38,m11_Homo_sapiens"
        )
    parts = [p.strip() for p in v.split(",")]
    return [p for p in parts if p]


def _split_tokens(values: Iterable[str]) -> List[str]:
    """
    Split tokens that may contain commas, returning cleaned non-empty items.
    Intended for optional multi-value arguments such as --outgroup-ids.
    """
    out: List[str] = []
    for v in values:
        if v is None:
            continue
        v = v.strip()
        if not v:
            continue
        parts = [p.strip() for p in v.split(",")]
        out.extend([p for p in parts if p])
    return out


def _derive_default_output_bed(input_maf: str) -> str:
    """
    Derive output BED path from input MAF path.
      foo.maf.gz -> foo.bed
      foo.maf    -> foo.bed
      otherwise  -> foo.bed (append)
    """
    p = Path(input_maf)
    name = p.name

    if name.endswith(".maf.gz"):
        out_name = name[:-7] + ".bed"  # remove ".maf.gz"
    elif name.endswith(".maf"):
        out_name = name[:-4] + ".bed"  # remove ".maf"
    else:
        out_name = name + ".bed"

    return str(p.with_name(out_name))


def build_arg_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        prog="realignpro maf2bed",
        description="Extract target-shared, others-absent positions from MAF and write BED3 intervals.",
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
        help="Output BED file path. If omitted, derived from input by replacing .maf(.gz) with .bed.",
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
        metavar="ID1,ID2,...",
        required=False,
        help="Target species IDs (comma-separated (no spaces)). REQUIRED in run mode.",
    )
    parser.add_argument(
        "--outgroup-ids",
        dest="outgroup_ids",
        nargs="*",
        default=None,
        help="Optional outgroup species IDs to exclude (space-separated and/or comma-separated).",
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


def _validate_and_build_config(args: argparse.Namespace) -> Maf2BedConfig:
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

    if not args.target_ids:
        raise ValueError("--target-ids is required in run mode.")
    target_ids = _split_csv_strict(args.target_ids, "--target-ids")
    if not target_ids:
        raise ValueError("--target-ids is empty after parsing.")
    outgroup_ids = _split_tokens(args.outgroup_ids or [])

    work_qsize = int(args.work_qsize)
    out_qsize = int(args.out_qsize)
    if work_qsize < 1:
        raise ValueError("--work-qsize must be >= 1.")
    if out_qsize < 1:
        raise ValueError("--out-qsize must be >= 1.")

    return Maf2BedConfig(
        input_maf=input_maf,
        output_bed=output_bed,
        total_threads=total_threads,
        ref_id=ref_id,
        target_ids=target_ids,
        outgroup_ids=outgroup_ids,
        work_qsize=work_qsize,
        out_qsize=out_qsize,
        start_method=args.start_method,
    )


# ----------------------------
# I/O helpers
# ----------------------------
def open_text_maybe_gzip(path: str):
    """
    Open .maf or .maf.gz in text mode.
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="replace", newline="")
    return open(path, "rt", encoding="utf-8", errors="replace", newline="")


# ----------------------------
# Utility: list assembly/species IDs in MAF
# ----------------------------
def list_maf_ids(maf_path: str) -> List[str]:
    """
    Return sorted unique species/assembly IDs found in the MAF file.

    Extraction rule:
      - From lines starting with 's'
      - Take the 2nd field (src)
      - Species/assembly ID is the substring before the first dot '.'
        e.g., 'hg38.chr1' -> 'hg38'
              'Taeniopygia_guttata.chr1' -> 'Taeniopygia_guttata'
              'chr1' (no dot) -> 'chr1' (kept as-is)
    """
    ids: Set[str] = set()
    with open_text_maybe_gzip(maf_path) as f:
        for raw in f:
            line = raw.strip()
            if not line or line[0] != "s":
                continue
            parts = line.split()
            if len(parts) < 2:
                continue
            src = parts[1]
            sid = src.split(".", 1)[0] if "." in src else src
            if sid:
                ids.add(sid)
    return sorted(ids)


# ----------------------------
# MAF parsing
# ----------------------------
def parse_maf_block_lines(block_lines: List[str], ref_id: str) -> Optional[Dict[str, Dict[str, Any]]]:
    """
    Parse one MAF block (paragraph ended by a blank line).

    Expected MAF 's' fields:
      s src start size strand srcSize text

    Returns:
      species -> {
        "seq": aligned sequence text,
        "chr": chromosome name (part after first '.'),
        "start": int, 0-based MAF start,
        "size": int, ungapped size,
        "strand": '+' or '-',
        "srcSize": int, source sequence length
      }

    Returns None if ref_id is missing in the block.
    """
    d: Dict[str, Dict[str, Any]] = {}

    for line in block_lines:
        if not line or line[0] != "s":
            continue
        parts = line.strip().split()
        if len(parts) < 7:
            continue

        src = parts[1]
        start = int(parts[2])
        size = int(parts[3])
        strand = parts[4]
        src_size = int(parts[5])
        text = parts[6]

        if "." in src:
            species, chrom = src.split(".", 1)
        else:
            species, chrom = src, src

        d[species] = {
            "seq": text,
            "chr": chrom,
            "start": start,
            "size": size,
            "strand": strand,
            "srcSize": src_size,
        }

    if ref_id not in d:
        return None
    return d


# ----------------------------
# Variant calling logic (with in-block merging)
# ----------------------------
def matrix2var(
    block: Dict[str, Dict[str, Any]],
    ref_id: str,
    target_ids: List[str],
    outgroup_ids: List[str],
) -> List[str]:
    """
    Identify target-shared, others-absent positions on the reference, then merge adjacent hits
    into blocks within the current MAF block to reduce output size.

    Output:
      - BED3 lines, UCSC 0-based half-open (chromStart inclusive, chromEnd exclusive)

    Strand handling:
      - '+' strand: walk forward
      - '-' strand: MAF start is relative to reverse-complemented source; convert to forward-strand
        coordinates while scanning the alignment text.
    """
    required = set(target_ids)
    outgroup = set(outgroup_ids)

    # Gate: process the block only if ALL targets exist in the block
    if not required.issubset(block.keys()):
        return []

    ref = block[ref_id]
    ref_seq = ref["seq"]
    ref_chr = ref["chr"]
    ref_start = ref["start"]      # 0-based MAF
    ref_strand = ref["strand"]
    ref_src_size = ref["srcSize"]

    # Initialize reference coordinate cursor
    if ref_strand == "+":
        pos_ref = ref_start - 1
        step = +1
    else:
        # Forward-strand interval is:
        #   [srcSize - start - size, srcSize - start)
        forward_end = ref_src_size - ref_start
        pos_ref = forward_end
        step = -1

    species_in_block = list(block.keys())

    # Merge adjacent hit bases on the fly
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
        other_nts: List[str] = []

        # Build allele lists for this alignment column
        for sid in species_in_block:
            # Skip outgroups entirely
            if sid in outgroup:
                continue

            sid_seq = block[sid]["seq"]
            nt = sid_seq[aln_idx].upper() if aln_idx < len(sid_seq) else "-"

            if sid in required:
                target_nts.append(nt)
            else:
                other_nts.append(nt)

        tset = set(target_nts)
        oset = set(other_nts)

        is_hit = False
        if len(tset) == 1:
            allele = next(iter(tset))
            if allele not in oset:
                is_hit = True

        if is_hit:
            base_start = pos_ref
            base_end = pos_ref + 1

            if cur_start is None:
                cur_start, cur_end = base_start, base_end
            else:
                if ref_strand == "+":
                    # Adjacent if new base starts exactly at current end
                    if base_start == cur_end:
                        cur_end = base_end
                    else:
                        flush_current()
                        cur_start, cur_end = base_start, base_end
                else:
                    # On '-' strand, scan goes high->low. Adjacent if base_end == cur_start.
                    if base_end == cur_start:
                        cur_start = base_start
                    else:
                        flush_current()
                        cur_start, cur_end = base_start, base_end
        else:
            flush_current()

    flush_current()

    return [f"{chrom}\t{start}\t{end}\n" for chrom, start, end in merged]


# ----------------------------
# Worker / writer
# ----------------------------
def worker_proc(
    work_q: "mp.JoinableQueue",
    out_q: "mp.Queue",
    ref_id: str,
    target_ids: List[str],
    outgroup_ids: List[str],
) -> None:
    """
    Worker process:
      - pull (block_index, block_dict)
      - compute merged BED lines for that block
      - push (block_index, bed_list) to out_q

    JoinableQueue contract:
      - every get() must be paired with exactly one task_done()
    """
    while True:
        item = work_q.get()
        try:
            if item is WORK_STOP:
                break
            idx, block = item
            bed_list = matrix2var(block, ref_id, target_ids, outgroup_ids)
            out_q.put((idx, bed_list))
        finally:
            work_q.task_done()

    out_q.put((OUT_STOP_IDX, None))


def writer_thread(out_q: "mp.Queue", out_path: str, n_workers: int) -> None:
    """
    Writer thread:
      - receive (idx, bed_list) from out_q
      - write in increasing idx order
      - buffer out-of-order results
    """
    next_idx = 0
    buffer: Dict[int, List[str]] = {}
    stop_count = 0

    with open(out_path, "w", encoding="utf-8", newline="") as fp:
        while True:
            idx, bed_list = out_q.get()

            if idx == OUT_STOP_IDX:
                stop_count += 1
                if stop_count == n_workers:
                    # Flush any remaining contiguous buffered results
                    while next_idx in buffer:
                        bl = buffer.pop(next_idx)
                        if bl:
                            fp.write("".join(bl))
                        next_idx += 1
                    break
                continue

            if idx == next_idx:
                if bed_list:
                    fp.write("".join(bed_list))
                next_idx += 1
                while next_idx in buffer:
                    bl = buffer.pop(next_idx)
                    if bl:
                        fp.write("".join(bl))
                    next_idx += 1
            else:
                buffer[idx] = bed_list or []


# ----------------------------
# Orchestration
# ----------------------------
def maf2bed_multiprocessing(
    maf_path: str,
    bed_path: str,
    ref_id: str,
    target_ids: List[str],
    outgroup_ids: List[str],
    total_threads: int,
    work_qsize: int,
    out_qsize: int,
    start_method: str = DEFAULT_START_METHOD,
) -> None:
    """
    Orchestrate the pipeline.

    total_threads includes:
      - 1 reader (main process reading/parsing MAF)
      - 1 writer (writer thread)
      - remaining are worker processes

    n_workers = max(1, total_threads - 2)
    """
    n_workers = max(1, total_threads - 2)

    ctx = mp.get_context(start_method)

    work_q: "mp.JoinableQueue" = ctx.JoinableQueue(maxsize=work_qsize)
    out_q: "mp.Queue" = ctx.Queue(maxsize=out_qsize)

    # Start workers
    procs = [
        ctx.Process(target=worker_proc, args=(work_q, out_q, ref_id, target_ids, outgroup_ids))
        for _ in range(n_workers)
    ]
    for p in procs:
        p.start()

    # Start writer
    t_writer = threading.Thread(target=writer_thread, args=(out_q, bed_path, n_workers), daemon=True)
    t_writer.start()

    # Stream blocks from MAF
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

            # End of block (blank line)
            if line.strip() == "":
                block = parse_maf_block_lines(block_lines, ref_id)
                if block is not None:
                    work_q.put((idx, block))
                    idx += 1
                in_block = False
                block_lines = []
                continue

            block_lines.append(line)

        # Flush last block if file doesn't end with blank line
        if in_block and block_lines:
            block = parse_maf_block_lines(block_lines, ref_id)
            if block is not None:
                work_q.put((idx, block))
                idx += 1

    # Stop workers
    for _ in range(n_workers):
        work_q.put(WORK_STOP)

    # Wait until all tasks are done
    work_q.join()

    # Join workers
    for p in procs:
        p.join()

    # Writer exits after receiving n_workers stop signals
    t_writer.join()


def main(argv: Optional[List[str]] = None) -> int:
    """
    Entry point.

    Modes:
      1) ID listing mode:
         maf2bed.py --ids input.maf
      2) Run mode:
         maf2bed.py --input <maf> --ref-id <id> --target-ids <id...> [options]
    """
    parser = build_arg_parser()
    args = parser.parse_args(argv)

    # Mode 1: --ids
    if args.ids_maf:
        try:
            ids = list_maf_ids(args.ids_maf)
        except Exception as e:
            print(f"[ERROR] failed to read MAF for IDs: {e}", file=sys.stderr)
            return 1
        print(", ".join(ids))
        return 0

    # Mode 2: run
    try:
        cfg = _validate_and_build_config(args)
    except Exception as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        parser.print_help(sys.stderr)
        return 1

    maf2bed_multiprocessing(
        maf_path=cfg.input_maf,
        bed_path=cfg.output_bed,
        ref_id=cfg.ref_id,
        target_ids=cfg.target_ids,
        outgroup_ids=cfg.outgroup_ids,
        total_threads=cfg.total_threads,
        work_qsize=cfg.work_qsize,
        out_qsize=cfg.out_qsize,
        start_method=cfg.start_method,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
