#!/usr/bin/env python3
from __future__ import annotations

"""
ReAlignPro_fa2maf.py
"""

import argparse
import concurrent.futures
import dataclasses
import datetime as _dt
import os
import re
import shlex
import shutil
import subprocess
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

DNA_BASES = set("ACGT")

# Optional helpers for common UCSC/RefSeq naming for GRCh38/hg38.
# These are NOT required for non-human references, but improve auto-resolution.
UCSC_TO_REFSEQ_GRCH38 = {
    "chr1": "NC_000001.11", "chr2": "NC_000002.12", "chr3": "NC_000003.12",
    "chr4": "NC_000004.12", "chr5": "NC_000005.10", "chr6": "NC_000006.12",
    "chr7": "NC_000007.14", "chr8": "NC_000008.11", "chr9": "NC_000009.12",
    "chr10": "NC_000010.11", "chr11": "NC_000011.10", "chr12": "NC_000012.12",
    "chr13": "NC_000013.11", "chr14": "NC_000014.9", "chr15": "NC_000015.10",
    "chr16": "NC_000016.10", "chr17": "NC_000017.11", "chr18": "NC_000018.10",
    "chr19": "NC_000019.10", "chr20": "NC_000020.11", "chr21": "NC_000021.9",
    "chr22": "NC_000022.11", "chrX": "NC_000023.11", "chrY": "NC_000024.10",
    "chrM": "NC_012920.1",
}
REFSEQ_TO_UCSC_GRCH38 = {v: k for k, v in UCSC_TO_REFSEQ_GRCH38.items()}


def tss_tag_2xwindow(window: int) -> str:
    return f"{2 * int(window)}bp"


# -----------------------------
# FASTA I/O
# -----------------------------
@dataclasses.dataclass
class FastaRecord:
    rid: str
    header: str
    seq: str


class FastaIO:
    @staticmethod
    def read_fasta(path: Path) -> List[FastaRecord]:
        records: List[FastaRecord] = []
        cur_id: Optional[str] = None
        cur_header: Optional[str] = None
        seq_parts: List[str] = []
        with path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.rstrip("\n")
                if not line:
                    continue
                if line.startswith(">"):
                    if cur_id is not None:
                        records.append(FastaRecord(rid=cur_id, header=cur_header or cur_id, seq="".join(seq_parts)))
                    cur_header = line[1:].strip()
                    cur_id = cur_header.split()[0]
                    seq_parts = []
                else:
                    seq_parts.append(line.strip())
            if cur_id is not None:
                records.append(FastaRecord(rid=cur_id, header=cur_header or cur_id, seq="".join(seq_parts)))
        if not records:
            raise RuntimeError(f"No FASTA records found in: {path}")
        return records

    @staticmethod
    def write_fasta(records: List[Tuple[str, str]], out_path: Path, wrap: int = 60) -> None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as out:
            for rid, seq in records:
                out.write(f">{rid}\n")
                for i in range(0, len(seq), wrap):
                    out.write(seq[i:i + wrap] + "\n")


# -----------------------------
# Sorting helpers
# -----------------------------
class OrthologSorter:
    _re_id_cov = re.compile(r"\|id%=(\d+(?:\.\d+)?)\|cov%=(\d+(?:\.\d+)?)\b")

    @staticmethod
    def id_cov_product(record_id: str) -> float:
        m = OrthologSorter._re_id_cov.search(record_id)
        if not m:
            return -1.0
        try:
            return float(m.group(1)) * float(m.group(2))
        except ValueError:
            return -1.0

    @staticmethod
    def sort_records(query_id: str, records: List[Tuple[str, str]]) -> List[Tuple[str, str]]:
        query_seq = None
        others: List[Tuple[str, str]] = []
        for rid, seq in records:
            if rid == query_id:
                query_seq = seq
            else:
                others.append((rid, seq))
        if query_seq is None:
            raise RuntimeError("Query record not found when sorting merged FASTA.")
        others_sorted = sorted(others, key=lambda x: (OrthologSorter.id_cov_product(x[0]), x[0]), reverse=True)
        return [(query_id, query_seq)] + others_sorted


class RecordIdParser:
    @staticmethod
    def assembly_from_merged_record_id(record_id: str) -> Optional[str]:
        parts = record_id.split("|")
        if len(parts) >= 2 and parts[1]:
            return parts[1]
        return None


# -----------------------------
# Reference contig resolver
# -----------------------------
class RefContigResolver:
    """
    Resolve which contig name in the reference genome (--ref_genome) corresponds to the contig used in query_id.

    Goals:
      - Work out-of-the-box when query contig names match the reference FASTA.
      - Handle common UCSC <-> RefSeq naming for GRCh38/hg38.
      - Provide safe heuristics (strip/add 'chr', chrM<->MT), otherwise fail fast if ambiguous.

    This resolver returns the reference FASTA contig NAME to use for srcSize lookups.
    """
    @staticmethod
    def _strip_chr(x: str) -> str:
        return x[3:] if x.startswith("chr") else x

    @staticmethod
    def _toggle_chr(x: str) -> str:
        return RefContigResolver._strip_chr(x) if x.startswith("chr") else ("chr" + x)

    @staticmethod
    def resolve(query_contig: str, ref_contigs: List[str]) -> str:
        ref_set = set(ref_contigs)

        # 1) Exact match
        if query_contig in ref_set:
            return query_contig

        # 2) Common UCSC<->RefSeq for GRCh38
        if query_contig in UCSC_TO_REFSEQ_GRCH38 and UCSC_TO_REFSEQ_GRCH38[query_contig] in ref_set:
            return UCSC_TO_REFSEQ_GRCH38[query_contig]
        if query_contig in REFSEQ_TO_UCSC_GRCH38 and REFSEQ_TO_UCSC_GRCH38[query_contig] in ref_set:
            return REFSEQ_TO_UCSC_GRCH38[query_contig]

        # 3) Toggle chr prefix
        toggled = RefContigResolver._toggle_chr(query_contig)
        if toggled in ref_set:
            return toggled

        # 4) chrM/MT special cases
        if query_contig == "chrM" and "MT" in ref_set:
            return "MT"
        if query_contig == "MT" and "chrM" in ref_set:
            return "chrM"

        # 5) Heuristic: unique suffix match (safe only if unique)
        q0 = RefContigResolver._strip_chr(query_contig)
        candidates = []
        for c in ref_contigs:
            c0 = RefContigResolver._strip_chr(c)
            if c0 == q0:
                candidates.append(c)
        if len(candidates) == 1:
            return candidates[0]

        # 6) Heuristic: exact end-match (rare, but sometimes ref contigs are like "gi|...|chr1")
        candidates = [c for c in ref_contigs if c.endswith(query_contig)]
        if len(candidates) == 1:
            return candidates[0]

        raise RuntimeError(
            "Cannot resolve reference contig for query contig "
            f"'{query_contig}'. Provide a reference FASTA that contains this contig name, "
            "or rename query coordinates/headers to match the reference FASTA contigs."
        )

    @staticmethod
    def acceptable_ref_contig_names(expected_query_contig: str) -> set:
        """
        Build a conservative set of acceptable contig names for reciprocal-hit checks.
        This is independent of reading the full ref FASTA index (kept lightweight for workers).
        """
        s = {expected_query_contig, RefContigResolver._toggle_chr(expected_query_contig)}
        if expected_query_contig in UCSC_TO_REFSEQ_GRCH38:
            s.add(UCSC_TO_REFSEQ_GRCH38[expected_query_contig])
        if expected_query_contig in REFSEQ_TO_UCSC_GRCH38:
            s.add(REFSEQ_TO_UCSC_GRCH38[expected_query_contig])
        if expected_query_contig == "chrM":
            s.add("MT")
        if expected_query_contig == "MT":
            s.add("chrM")
        return s


# -----------------------------
# Tree utilities
# -----------------------------
class NewickTipOrder:
    _tip_pat = re.compile(r"(m\d+_[A-Za-z0-9_.]+)")

    @staticmethod
    def read_tip_order(nwk_path: Path) -> List[str]:
        txt = nwk_path.read_text(encoding="utf-8", errors="replace")
        compact = re.sub(r"\s+", "", txt)
        seen = set()
        order: List[str] = []
        for m in NewickTipOrder._tip_pat.finditer(compact):
            tip = m.group(1)
            if tip not in seen:
                seen.add(tip)
                order.append(tip)
        if not order:
            raise RuntimeError(f"No tips parsed from Newick: {nwk_path}")
        return order

    @staticmethod
    def write_tree_ordered_fasta(
        merged_records: List[Tuple[str, str]],
        query_id: str,
        tree_tips: List[str],
        out_path: Path,
    ) -> None:
        query_seq: Optional[str] = None
        asm_to_records: Dict[str, List[Tuple[str, str]]] = {}
        for rid, seq in merged_records:
            if rid == query_id:
                query_seq = seq
                continue
            asm = RecordIdParser.assembly_from_merged_record_id(rid)
            if asm is None:
                continue
            asm_to_records.setdefault(asm, []).append((rid, seq))

        if query_seq is None:
            raise RuntimeError("Query record not found when creating tree-ordered FASTA.")

        out_records: List[Tuple[str, str]] = [(query_id, query_seq)]
        for asm in tree_tips:
            if asm not in asm_to_records:
                continue
            recs_sorted = sorted(
                asm_to_records[asm],
                key=lambda x: (OrthologSorter.id_cov_product(x[0]), x[0]),
                reverse=True,
            )
            out_records.extend(recs_sorted)

        FastaIO.write_fasta(out_records, out_path)


# -----------------------------
# Progress
# -----------------------------
class ProgressReporter:
    def __init__(self, total: int, enabled: bool = True):
        self.total = max(1, total)
        self.enabled = enabled
        self.start_time = time.time()
        self.done = 0

    def start(self, workers: int, stage: str, extra: str = "") -> None:
        if not self.enabled:
            return
        msg = f"[START:{stage}] total={self.total} workers={workers}"
        if extra:
            msg += f" | {extra}"
        print(msg, flush=True)

    def update(self, label: str, forward_pass: int, rbh_pass: int, accepted: int) -> None:
        if not self.enabled:
            return
        self.done += 1
        pct = 100.0 * self.done / float(self.total)
        elapsed = time.time() - self.start_time
        print(
            f"[{self.done:3d}/{self.total:<3d} {pct:5.1f}%] {label} | "
            f"forward_pass={forward_pass} rbh_pass={rbh_pass} accepted={accepted} | elapsed={elapsed:.0f}s",
            flush=True,
        )

    def finish(self, stage: str) -> None:
        if not self.enabled:
            return
        elapsed = time.time() - self.start_time
        print(f"[FINISH:{stage}] elapsed={elapsed:.0f}s", flush=True)


# -----------------------------
# External command runner
# -----------------------------
class CommandRunner:
    def __init__(self, dry_run: bool = False):
        self.dry_run = dry_run

    def run(self, cmd: List[str], log_path: Optional[Path] = None, cwd: Optional[Path] = None) -> None:
        if self.dry_run:
            print("[DRY-RUN]", " ".join(cmd), flush=True)
            return

        if log_path is not None:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open("w", encoding="utf-8") as log:
                log.write(f"# CMD: {' '.join(cmd)}\n")
                log.write(f"# PWD: {str(cwd) if cwd else os.getcwd()}\n")
                log.write(f"# TIME: {_dt.datetime.now().isoformat()}\n\n")
                proc = subprocess.run(
                    cmd,
                    cwd=str(cwd) if cwd else None,
                    stdout=log,
                    stderr=subprocess.STDOUT,
                    check=False,
                    text=True,
                )
                if proc.returncode != 0:
                    raise RuntimeError(f"Command failed (exit {proc.returncode}). See log: {log_path}")
        else:
            proc = subprocess.run(cmd, cwd=str(cwd) if cwd else None, check=False)
            if proc.returncode != 0:
                raise RuntimeError(f"Command failed (exit {proc.returncode}): {' '.join(cmd)}")


# -----------------------------
# FASTA index (.fai) helper
# -----------------------------
class FaiIndex:
    def __init__(self):
        self.cache: Dict[Path, Dict[str, int]] = {}

    @staticmethod
    def read_fai(fai_path: Path) -> Dict[str, int]:
        d: Dict[str, int] = {}
        with fai_path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                if not line.strip():
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 2:
                    continue
                d[parts[0]] = int(parts[1])
        return d

    def ensure_and_load(self, fasta_path: Path, samtools: str, runner: CommandRunner, log_path: Path) -> Dict[str, int]:
        fai_path = fasta_path.with_suffix(fasta_path.suffix + ".fai")
        if fasta_path not in self.cache:
            if not fai_path.exists():
                runner.run([samtools, "faidx", str(fasta_path)], log_path=log_path)
            self.cache[fasta_path] = self.read_fai(fai_path)
        return self.cache[fasta_path]

    def get_len(self, fasta_path: Path, contig: str, samtools: str, runner: CommandRunner, log_path: Path) -> int:
        d = self.ensure_and_load(fasta_path, samtools, runner, log_path)
        if contig not in d:
            raise RuntimeError(f"Contig '{contig}' not found in fai for {fasta_path}")
        return d[contig]

    def contigs(self, fasta_path: Path, samtools: str, runner: CommandRunner, log_path: Path) -> List[str]:
        d = self.ensure_and_load(fasta_path, samtools, runner, log_path)
        return list(d.keys())


# -----------------------------
# LASTZ
# -----------------------------
@dataclasses.dataclass
class AlignmentHit:
    score: int
    target_name: str
    target_strand: str
    target_zstart: int
    target_end: int
    query_name: str
    query_strand: str
    query_zstart_plus: int
    query_end_plus: int
    identity_pct: float
    coverage_pct: float

    @property
    def target_start_1based(self) -> int:
        return self.target_zstart + 1

    @property
    def target_end_1based(self) -> int:
        return self.target_end


@dataclasses.dataclass
class Best2:
    best: AlignmentHit
    second: Optional[AlignmentHit] = None

    def ratio(self) -> float:
        if self.second is None or self.second.score <= 0:
            return float("inf")
        return self.best.score / float(self.second.score)


class LastzRunner:
    def __init__(self, lastz_path: Path, runner: CommandRunner, extra_params: str):
        # Accept either an explicit filesystem path or a binary name available in PATH.
        resolved = lastz_path
        if not resolved.exists():
            found = shutil.which(str(lastz_path))
            if found:
                resolved = Path(found)
        self.lastz_path = resolved
        self.runner = runner
        self.extra_params = extra_params
        if not self.lastz_path.exists():
            raise FileNotFoundError(f"LASTZ not found: {lastz_path} (resolved: {self.lastz_path})")

    @staticmethod
    def format_fields() -> str:
        return "score,name1,strand1,zstart1,end1,name2,strand2,zstart2+,end2+,id%,cov%"

    @staticmethod
    def _parse_pct(x: str) -> float:
        x = x.strip()
        if x.endswith("%"):
            x = x[:-1]
        return float(x)

    def run_lastz(self, target_fa: Path, query_fa: Path, out_tsv: Path, log_path: Path) -> None:
        """Run LASTZ and write results to out_tsv.

        Implementation note:
          - LASTZ writes its tabular output directly to a file via --output.
          - To ensure provenance is captured inside the TSV itself (not only in a separate .log),
            we first write to a temporary file and then prepend metadata lines (e.g., '# CMD: ...')
            to the final TSV (method A).
        """
        out_tsv.parent.mkdir(parents=True, exist_ok=True)

        # Prepare LASTZ command
        target_arg = f"{target_fa}[multiple,unmask,nameparse=darkspace]"
        query_arg = f"{query_fa}[multiple,unmask,nameparse=darkspace]"
        cmd: List[str] = [str(self.lastz_path), target_arg, query_arg]
        if self.extra_params.strip():
            cmd.extend(shlex.split(self.extra_params.strip()))

        # Write output to a temp file first, then prepend provenance to the final TSV.
        tmp_tsv = out_tsv.with_suffix(out_tsv.suffix + ".tmp")
        if tmp_tsv.exists():
            tmp_tsv.unlink()

        cmd.extend([f"--format=general:{self.format_fields()}", f"--output={str(tmp_tsv)}"])

        # Run and capture stderr/stdout to log_path (CommandRunner already writes # CMD/#CWD/#TIME there).
        self.runner.run(cmd, log_path=log_path)

        # Prepend the executed command (and context) to the final TSV for reproducibility.
        # NOTE: we use the same cmd list (including --output=<tmp>) to reflect the actual call.
        executed_cmd_str = " ".join(shlex.quote(x) for x in cmd)
        cwd_str = str(Path.cwd())
        time_str = _dt.datetime.now().isoformat()

        with out_tsv.open("w", encoding="utf-8") as out:
            out.write(f"# CMD: {executed_cmd_str}\n")
            out.write(f"# PWD: {cwd_str}\n")
            out.write(f"# TIME: {time_str}\n")
            if tmp_tsv.exists() and tmp_tsv.stat().st_size > 0:
                with tmp_tsv.open("r", encoding="utf-8", errors="replace") as inp:
                    shutil.copyfileobj(inp, out)

        # Best-effort cleanup
        if tmp_tsv.exists():
            tmp_tsv.unlink()

    def parse_general(self, tsv_path: Path) -> List[AlignmentHit]:
        hits: List[AlignmentHit] = []
        if not tsv_path.exists() or tsv_path.stat().st_size == 0:
            return hits
        with tsv_path.open("r", encoding="utf-8", errors="replace") as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 11:
                    continue
                try:
                    score = int(float(parts[0]))
                    name1 = parts[1]
                    strand1 = parts[2]
                    zstart1 = int(parts[3])
                    end1 = int(parts[4])
                    name2 = parts[5]
                    strand2 = parts[6]
                    zstart2p = int(parts[7])
                    end2p = int(parts[8])
                    idpct = self._parse_pct(parts[9])
                    covpct = self._parse_pct(parts[10])
                except ValueError:
                    continue
                hits.append(AlignmentHit(
                    score=score,
                    target_name=name1,
                    target_strand=strand1,
                    target_zstart=zstart1,
                    target_end=end1,
                    query_name=name2,
                    query_strand=strand2,
                    query_zstart_plus=zstart2p,
                    query_end_plus=end2p,
                    identity_pct=idpct,
                    coverage_pct=covpct,
                ))
        return hits


# -----------------------------
# samtools extractor
# -----------------------------
class FastaExtractor:
    def __init__(self, runner: CommandRunner, samtools: str = "samtools"):
        self.runner = runner
        self.samtools = samtools

    def has_samtools(self) -> bool:
        return shutil.which(self.samtools) is not None

    @staticmethod
    def revcomp(seq: str) -> str:
        comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
        return seq.translate(comp)[::-1]

    def ensure_fai(self, fasta_path: Path, log_path: Path) -> None:
        if fasta_path.with_suffix(fasta_path.suffix + ".fai").exists():
            return
        if not self.has_samtools():
            raise RuntimeError("samtools not found in PATH; cannot index FASTA.")
        self.runner.run([self.samtools, "faidx", str(fasta_path)], log_path=log_path)

    def fetch(self, fasta_path: Path, contig: str, start_1based: int, end_1based: int, log_path: Path) -> str:
        if not self.has_samtools():
            raise RuntimeError("samtools not found in PATH; cannot extract subsequences.")
        region = f"{contig}:{start_1based}-{end_1based}"
        cmd = [self.samtools, "faidx", str(fasta_path), region]
        if self.runner.dry_run:
            return "N" * max(0, end_1based - start_1based + 1)
        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            log_path.parent.mkdir(parents=True, exist_ok=True)
            with log_path.open("a", encoding="utf-8") as log:
                log.write(proc.stderr + "\n")
            raise RuntimeError(f"samtools faidx failed for region {region}")
        lines = [ln.strip() for ln in proc.stdout.splitlines() if ln and not ln.startswith(">")]
        return "".join(lines)


# -----------------------------
# Robust query-id parsing
# -----------------------------
class QueryIdParser:
    """
    Parse query IDs into genomic intervals.

    This pipeline supports both legacy underscore-delimited IDs and the new pipe-delimited IDs.

    Pipe-delimited formats (preferred; safe when contig/scaffold names contain underscores):
      - <build>|<contig>|<start>|<end>|<strand>
      - <build>|<contig>|<start>|<end>
      - <contig>|<start>|<end>|<strand>
      - <contig>|<start>|<end>

    Legacy underscore-delimited formats (backward compatible):
      - chr16_6018803_6019902
      - hg38_chr16_6018803_6019902
      - scaffold_foo_bar_100_500   (contig can include underscores; last two tokens must be integers)
    """
    _build_token = re.compile(
        r"^(hg\d+|mm\d+|rn\d+|dm\d+|ce\d+|bosTau\d+|canFam\d+|danRer\d+|galGal\d+|GRCh\d+|GRCm\d+)$",
        re.IGNORECASE,
    )

    @staticmethod
    def parse_extended(query_id: str) -> Tuple[Optional[str], str, int, int, str]:
        query_id = query_id.strip()
        if not query_id:
            raise RuntimeError("Empty query ID.")

        # Preferred: pipe-delimited
        if "|" in query_id:
            toks = [t for t in query_id.split("|") if t != ""]
            if len(toks) < 3:
                raise RuntimeError(f"Cannot parse pipe-delimited query ID (too few fields): {query_id}")

            strand = "+"
            if toks[-1] in {"+", "-"}:
                strand = toks[-1]
                toks = toks[:-1]

            if len(toks) < 3:
                raise RuntimeError(f"Cannot parse pipe-delimited query ID (missing contig/start/end): {query_id}")

            try:
                start = int(toks[-2])
                end = int(toks[-1])
            except ValueError:
                raise RuntimeError(f"Pipe-delimited query ID must end with |<start>|<end>[|<strand>]: {query_id}")

            contig = toks[-3]
            prefix = toks[:-3]
            build: Optional[str] = None
            if prefix:
                # If the first token looks like a build label, treat it as build.
                if QueryIdParser._build_token.match(prefix[0]):
                    build = prefix[0]
                else:
                    # Otherwise keep build undefined (but preserve backward compatibility).
                    build = None

            if start > end:
                start, end = end, start
            return build, contig, start, end, strand

        # Legacy: underscore-delimited (last two underscore tokens must be integers).
        toks = query_id.split("_")
        if len(toks) < 3:
            raise RuntimeError(f"Cannot parse query interval from query ID: {query_id}")
        if not (toks[-1].isdigit() and toks[-2].isdigit()):
            raise RuntimeError(f"Query ID must end with _<start>_<end>: {query_id}")

        start = int(toks[-2])
        end = int(toks[-1])
        if start > end:
            start, end = end, start

        contig_tokens = toks[:-2]
        if len(contig_tokens) >= 2 and QueryIdParser._build_token.match(contig_tokens[0]):
            build = contig_tokens[0]
            contig_tokens = contig_tokens[1:]
        else:
            build = None

        contig = "_".join(contig_tokens)
        return build, contig, start, end, "+"

    @staticmethod
    def parse(query_id: str) -> Tuple[str, int, int]:
        _build, contig, start, end, _strand = QueryIdParser.parse_extended(query_id)
        return contig, start, end

# -----------------------------
# RBH
# -----------------------------
class RecordParserRBH:
    @staticmethod
    def parse_species_from_assembly(assembly_id: str) -> str:
        m = re.match(r"^m\d+_(.+)$", assembly_id)
        return m.group(1) if m else assembly_id

    @staticmethod
    def needs_revcomp(target_strand: str, query_strand: str) -> bool:
        return (target_strand.strip() == "-") ^ (query_strand.strip() == "-")


class RBHPipeline:
    def __init__(
        self,
        lastz: LastzRunner,
        extractor: FastaExtractor,
        genome_dir: Path,
        work_dir: Path,
        query_fa: Path,
        ref_genome: Path,
        min_id_pct: float,
        min_cov_pct: float,
        min_score_ratio: float,
        min_ref_overlap_bp: int,
        min_ref_overlap_frac: float,
        skip_forward_lastz: bool,
        dry_run: bool,
    ):
        self.lastz = lastz
        self.extractor = extractor
        self.genome_dir = genome_dir
        self.work_dir = work_dir
        self.query_fa = query_fa
        self.ref_genome = ref_genome
        self.min_id_pct = min_id_pct
        self.min_cov_pct = min_cov_pct
        self.min_score_ratio = min_score_ratio
        self.min_ref_overlap_bp = min_ref_overlap_bp
        self.min_ref_overlap_frac = min_ref_overlap_frac
        self.skip_forward_lastz = skip_forward_lastz
        self.dry_run = dry_run
        self.query_records = {r.rid: r for r in FastaIO.read_fasta(self.query_fa)}

    @staticmethod
    def list_genomes(genome_dir: Path) -> List[Path]:
        return sorted([
            p for p in genome_dir.iterdir()
            if p.is_file() and p.suffix.lower() in {".fa", ".fasta", ".fna"} and not p.name.endswith(".fai")
        ])

    @staticmethod
    def pick_best2(hits: List[AlignmentHit]) -> Optional[Best2]:
        if not hits:
            return None
        srt = sorted(hits, key=lambda h: (h.score, h.coverage_pct, h.identity_pct), reverse=True)
        return Best2(best=srt[0], second=(srt[1] if len(srt) > 1 else None))

    @staticmethod
    def overlap_1based(a1: int, a2: int, b1: int, b2: int) -> int:
        lo = max(a1, b1)
        hi = min(a2, b2)
        return max(0, hi - lo + 1)

    def _reciprocal_pass(self, qid: str, extracted_fa: Path, out_dir: Path, log_dir: Path, qtag: str = "") -> bool:
        base = "reciprocal.lastz" if not qtag else f"reciprocal.{qtag}.lastz"
        reciprocal_tsv = out_dir / f"{base}.tsv"
        reciprocal_log = log_dir / f"{base}.log"
        self.lastz.run_lastz(target_fa=self.ref_genome, query_fa=extracted_fa, out_tsv=reciprocal_tsv, log_path=reciprocal_log)
        hits = self.lastz.parse_general(reciprocal_tsv)
        b2 = self.pick_best2(hits)
        if b2 is None:
            return False
        best = b2.best

        try:
            exp_contig, exp_start, exp_end = QueryIdParser.parse(qid)
        except RuntimeError:
            return True

        allowed = RefContigResolver.acceptable_ref_contig_names(exp_contig)
        if best.target_name not in allowed:
            return False

        ov = self.overlap_1based(best.target_start_1based, best.target_end_1based, exp_start, exp_end)
        exp_len = exp_end - exp_start + 1
        ov_frac = ov / float(exp_len) if exp_len > 0 else 0.0
        return bool(ov >= self.min_ref_overlap_bp and ov_frac >= self.min_ref_overlap_frac)

    def process_one_genome(self, genome_fa: Path) -> Tuple[str, int, int, int, List[Tuple[str, str]]]:
        assembly_id = genome_fa.stem
        species = RecordParserRBH.parse_species_from_assembly(assembly_id)

        out_dir = self.work_dir / assembly_id
        log_dir = out_dir / "logs"
        out_dir.mkdir(parents=True, exist_ok=True)
        log_dir.mkdir(parents=True, exist_ok=True)

        forward_tsv = out_dir / "forward.lastz.tsv"
        forward_log = log_dir / "forward.lastz.log"
        fai_log = log_dir / "samtools_faidx.log"
        fetch_log = log_dir / "samtools_fetch.log"

        if not self.skip_forward_lastz:
            self.lastz.run_lastz(target_fa=genome_fa, query_fa=self.query_fa, out_tsv=forward_tsv, log_path=forward_log)
        else:
            if (not forward_tsv.exists()) or (forward_tsv.stat().st_size == 0):
                return assembly_id, 0, 0, 0, []

        hits = self.lastz.parse_general(forward_tsv)
        if not hits:
            return assembly_id, 0, 0, 0, []

        by_qid: Dict[str, List[AlignmentHit]] = {}
        for h in hits:
            by_qid.setdefault(h.query_name, []).append(h)

        self.extractor.ensure_fai(genome_fa, log_path=fai_log)

        qids_total = len(by_qid)
        forward_pass = 0
        rbh_pass = 0
        accepted: List[Tuple[str, str]] = []
        multi_query = (qids_total > 1)

        for q_index, (qid, qhits) in enumerate(sorted(by_qid.items(), key=lambda x: x[0]), start=1):
            b2 = self.pick_best2(qhits)
            if b2 is None:
                continue
            best = b2.best
            ratio = b2.ratio()

            ok = True
            if best.identity_pct < self.min_id_pct:
                ok = False
            if best.coverage_pct < self.min_cov_pct:
                ok = False
            if ratio < self.min_score_ratio:
                ok = False
            if not ok:
                continue
            forward_pass += 1

            seq_plus = self.extractor.fetch(
                fasta_path=genome_fa,
                contig=best.target_name,
                start_1based=best.target_start_1based,
                end_1based=best.target_end_1based,
                log_path=fetch_log,
            )
            needs_rc = RecordParserRBH.needs_revcomp(best.target_strand, best.query_strand)
            seq_oriented = self.extractor.revcomp(seq_plus) if needs_rc else seq_plus
            qtag = f"q{q_index}" if multi_query else ""
            tmp_extracted = out_dir / (f"tmp_extracted.{qtag}.fasta" if qtag else "tmp_extracted.fasta")

            maf_strand = "-" if needs_rc else "+"

            tmp_id = (
                f"{species}|{assembly_id}|{best.target_name}|{best.target_start_1based}|{best.target_end_1based}|{maf_strand}|"
                f"tAlnStrand={best.target_strand}|qAlnStrand={best.query_strand}"
            )
            FastaIO.write_fasta([(tmp_id, seq_oriented)], tmp_extracted)

            rec_ok = False
            try:
                rec_ok = self._reciprocal_pass(qid=qid, extracted_fa=tmp_extracted, out_dir=out_dir, log_dir=log_dir, qtag=qtag)
            finally:
                tmp_extracted.unlink(missing_ok=True)

            if not rec_ok:
                continue

            rbh_pass += 1
            record_id = (
                f"{species}|{assembly_id}|{best.target_name}|{best.target_start_1based}|{best.target_end_1based}|{maf_strand}|"
                f"tAlnStrand={best.target_strand}|qAlnStrand={best.query_strand}|"
                f"score={best.score}|id%={best.identity_pct:.2f}|cov%={best.coverage_pct:.2f}"
            )
            accepted.append((record_id, seq_oriented))

        if accepted:
            accepted_sorted = sorted(accepted, key=lambda x: (OrthologSorter.id_cov_product(x[0]), x[0]), reverse=True)
            FastaIO.write_fasta(accepted_sorted, out_dir / "orthologs.fasta")

        return assembly_id, qids_total, forward_pass, rbh_pass, accepted

    def get_primary_query_id(self) -> str:
        # Prefer the first record, but keep backward compatibility for common patterns.
        for rid in self.query_records.keys():
            if "|" in rid:
                # pipe-delimited: <build>|<contig>|<start>|<end>[|<strand>]
                return rid
            if rid.startswith("hg") and "_chr" in rid:
                return rid
            if rid.startswith("chr"):
                return rid
        return sorted(self.query_records.keys())[0]


def _worker_rbh(args: dict, genome_path: str) -> Tuple[str, int, int, int, List[Tuple[str, str]]]:
    runner = CommandRunner(dry_run=args["dry_run"])
    lastz = LastzRunner(lastz_path=Path(args["lastz"]), runner=runner, extra_params=args["lastz_params"])
    extractor = FastaExtractor(runner=runner, samtools=args["samtools"])
    pipe = RBHPipeline(
        lastz=lastz,
        extractor=extractor,
        genome_dir=Path(args["genome_dir"]),
        work_dir=Path(args["work_dir"]),
        query_fa=Path(args["query_fa"]),
        ref_genome=Path(args["ref_genome"]),
        min_id_pct=args["min_id_pct"],
        min_cov_pct=args["min_cov_pct"],
        min_score_ratio=args["min_score_ratio"],
        min_ref_overlap_bp=args["min_ref_overlap_bp"],
        min_ref_overlap_frac=args["min_ref_overlap_frac"],
        skip_forward_lastz=args["skip_forward_lastz"],
        dry_run=args["dry_run"],
    )
    return pipe.process_one_genome(Path(genome_path))


# -----------------------------
# MUSCLE
# -----------------------------
class MuscleRunner:
    """Run MUSCLE with automatic CLI adaptation for MUSCLE v3 and v5.

    MUSCLE v3 uses:
        muscle -in <in.fa> -out <out.fa> [options]

    MUSCLE v5 uses:
        muscle -align <in.fa> -output <out.fa> [options]
    """

    _re_major = re.compile(r"(?:MUSCLE\s*(?:v|V)?\s*)?(\d+)(?:\.\d+)*")

    def __init__(self, muscle_bin: str = "muscle"):
        self.muscle_bin = muscle_bin
        if shutil.which(self.muscle_bin) is None:
            raise RuntimeError(f"MUSCLE not found in PATH: {self.muscle_bin}")

        self.version_text = self._detect_version_text()
        self.major_version = self._parse_major_version(self.version_text)

    @staticmethod
    def _parse_major_version(version_text: str) -> int:
        if not version_text:
            return 3
        m = MuscleRunner._re_major.search(version_text)
        if m:
            try:
                return int(m.group(1))
            except ValueError:
                return 3
        # Fallback: first digit in the output
        m2 = re.search(r"(\d+)", version_text)
        if m2:
            try:
                return int(m2.group(1))
            except ValueError:
                return 3
        return 3

    def _detect_version_text(self) -> str:
        """Best-effort MUSCLE version detection.

        MUSCLE v3 and v5 differ in supported flags. We attempt a few common version
        flags and return combined stdout/stderr.
        """
        candidates = [
            [self.muscle_bin, "-version"],
            [self.muscle_bin, "--version"],
            [self.muscle_bin, "-V"],
        ]
        for cmd in candidates:
            try:
                proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            except Exception:
                continue
            out = (proc.stdout or "") + (proc.stderr or "")
            if out.strip():
                return out.strip()

        # Final fallback: ask for help and try to infer major version from it
        for cmd in ([self.muscle_bin, "-h"], [self.muscle_bin, "--help"], [self.muscle_bin, "-help"]):
            try:
                proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            except Exception:
                continue
            out = (proc.stdout or "") + (proc.stderr or "")
            if out.strip():
                return out.strip()

        return ""

    @staticmethod
    def _strip_io_flags(tokens: List[str]) -> List[str]:
        """Remove input/output flags that would conflict with our own base command."""
        io_flags = {"-in", "-out", "-align", "-output"}
        out: List[str] = []
        i = 0
        while i < len(tokens):
            t = tokens[i]
            if t in io_flags:
                # Skip flag and its value if present
                i += 1
                if i < len(tokens) and not tokens[i].startswith("-"):
                    i += 1
                continue
            out.append(t)
            i += 1
        return out

    @staticmethod
    def _translate_v3_to_v5(tokens: List[str]) -> List[str]:
        """Translate a small subset of MUSCLE v3 flags to their v5 equivalents.

        This is intentionally conservative and only translates flags we know are
        common in ReAlignPro defaults. Other flags are passed through as-is.
        """
        out: List[str] = []
        i = 0
        while i < len(tokens):
            t = tokens[i]

            # v3: -maxiters N  -> v5: -refineiters N (closest analogue)
            if t == "-maxiters" and i + 1 < len(tokens):
                n = tokens[i + 1]
                out.extend(["-refineiters", n])
                i += 2
                continue

            # v3 optimization flag not supported in v5
            if t == "-diags":
                i += 1
                continue

            out.append(t)
            i += 1
        return out

    def run(self, in_fa: Path, out_fa: Path, extra_args: Optional[List[str]] = None) -> None:
        out_fa.parent.mkdir(parents=True, exist_ok=True)

        # Normalize extra args (remove in/out, then adapt if needed)
        extra: List[str] = []
        if extra_args:
            extra = self._strip_io_flags(list(extra_args))
            if self.major_version >= 5:
                extra = self._translate_v3_to_v5(extra)

        # Build MUSCLE command based on detected major version
        if self.major_version >= 5:
            # MUSCLE v5 CLI: -align <in.fa> -output <out.fa>
            cmd = [self.muscle_bin, "-align", str(in_fa), "-output", str(out_fa)]
        else:
            # MUSCLE v3 CLI: -in <in.fa> -out <out.fa>
            cmd = [self.muscle_bin, "-in", str(in_fa), "-out", str(out_fa), "-quiet"]

        if extra:
            cmd.extend(extra)

        proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        if proc.returncode != 0:
            raise RuntimeError(
                "MUSCLE failed.\n"
                f"Detected MUSCLE major version: {self.major_version}\n"
                f"CMD: {' '.join(cmd)}\n"
                f"STDERR:\n{proc.stderr}\n"
            )
def reorder_aligned_fasta_inplace(aligned_fa: Path, desired_order: List[str]) -> None:
    recs = FastaIO.read_fasta(aligned_fa)
    by_id = {r.rid: r.seq for r in recs}
    missing = [rid for rid in desired_order if rid not in by_id]
    if missing:
        raise RuntimeError(f"Aligned FASTA missing expected IDs: {missing[:5]}")
    FastaIO.write_fasta([(rid, by_id[rid]) for rid in desired_order], aligned_fa)


# -----------------------------
# Header tools: sanitize + header_map
# -----------------------------
class HeaderTools:
    # Legacy: contig:start-end
    _re_loc_legacy = re.compile(r"^([^:]+):(\d+)-(\d+)$")

    # New/legacy strand key/value forms
    _re_tstrand = re.compile(r"(?:tAlnStrand|tstrand)=([+-])")
    _re_qstrand = re.compile(r"(?:qAlnStrand|qstrand)=([+-])")

    _re_score = re.compile(r"score=(\d+)")
    _re_id = re.compile(r"id%=(\d+(?:\.\d+)?)")
    _re_cov = re.compile(r"cov%=(\d+(?:\.\d+)?)")

    @staticmethod
    def parse_query_interval(query_id: str) -> Tuple[str, int, int]:
        return QueryIdParser.parse(query_id)

    @staticmethod
    def parse_query_interval_extended(query_id: str) -> Tuple[Optional[str], str, int, int, str]:
        return QueryIdParser.parse_extended(query_id)

    @staticmethod
    def parse_ortholog_meta(original_id: str) -> Dict[str, Optional[str]]:
        """
        Parse ortholog headers.

        Supported formats:

        (A) New safe, pipe-delimited (no nested pipes in any single field):
            species|assembly|qcontig|qstart|qend|qstrand|tcontig|tstart|tend|maf_strand|tAlnStrand=+|qAlnStrand=-|score=...|id%=...|cov%=...

        (B) Legacy (pre-Feb 2026):
            species|assembly|qid|tcontig:tstart-tend|tstrand=...|qstrand=...|score=...|id%=...|cov%=...
        """
        parts = original_id.split("|")
        d: Dict[str, Optional[str]] = {
            "species": None, "assembly": None, "qid": None,
            "contig": None, "start_1based": None, "end_1based": None,
            "tstrand": None, "qstrand": None,
            "maf_strand": None,
            "score": None, "id_pct": None, "cov_pct": None,
        }

        # New format (A2): species|assembly|tcontig|tstart|tend|maf_strand|tAlnStrand=...|qAlnStrand=...|score=...|id%=...|cov%=...
        if (
            len(parts) >= 6
            and parts[3].isdigit() and parts[4].isdigit()
            and parts[5] in {"+", "-"}
        ):
            d["species"] = parts[0]
            d["assembly"] = parts[1]
            d["qid"] = ""
            d["contig"] = parts[2]
            d["start_1based"] = parts[3]
            d["end_1based"] = parts[4]
            d["maf_strand"] = parts[5]

            rest = "|".join(parts[6:])
            mt = HeaderTools._re_tstrand.search(rest)
            mq = HeaderTools._re_qstrand.search(rest)
            if mt:
                d["tstrand"] = mt.group(1)
            if mq:
                d["qstrand"] = mq.group(1)

            ms = HeaderTools._re_score.search(rest)
            if ms:
                d["score"] = ms.group(1)
            mid = HeaderTools._re_id.search(rest)
            if mid:
                d["id_pct"] = mid.group(1)
            mc = HeaderTools._re_cov.search(rest)
            if mc:
                d["cov_pct"] = mc.group(1)
            return d

        # New format (A)
        if (
            len(parts) >= 10
            and parts[3].isdigit() and parts[4].isdigit()
            and parts[7].isdigit() and parts[8].isdigit()
            and parts[9] in {"+", "-"}
        ):
            d["species"] = parts[0]
            d["assembly"] = parts[1]
            qcontig, qstart, qend, qstrand_id = parts[2], parts[3], parts[4], parts[5]
            d["qid"] = f"{qcontig}|{qstart}|{qend}|{qstrand_id}"
            d["contig"] = parts[6]
            d["start_1based"] = parts[7]
            d["end_1based"] = parts[8]
            d["maf_strand"] = parts[9]

            rest = "|".join(parts[10:])
            mt = HeaderTools._re_tstrand.search(rest)
            mq = HeaderTools._re_qstrand.search(rest)
            if mt:
                d["tstrand"] = mt.group(1)
            if mq:
                d["qstrand"] = mq.group(1)

            ms = HeaderTools._re_score.search(rest)
            if ms:
                d["score"] = ms.group(1)
            mid = HeaderTools._re_id.search(rest)
            if mid:
                d["id_pct"] = mid.group(1)
            mc = HeaderTools._re_cov.search(rest)
            if mc:
                d["cov_pct"] = mc.group(1)
            return d

        # Legacy format (B)
        if len(parts) < 4:
            return d
        d["species"] = parts[0]
        d["assembly"] = parts[1]
        d["qid"] = parts[2]

        loc = parts[3]
        mloc = HeaderTools._re_loc_legacy.match(loc)
        if mloc:
            d["contig"] = mloc.group(1)
            d["start_1based"] = mloc.group(2)
            d["end_1based"] = mloc.group(3)

        rest = "|".join(parts[4:])
        mt = HeaderTools._re_tstrand.search(rest)
        mq = HeaderTools._re_qstrand.search(rest)
        if mt:
            d["tstrand"] = mt.group(1)
        if mq:
            d["qstrand"] = mq.group(1)

        ms = HeaderTools._re_score.search(rest)
        if ms:
            d["score"] = ms.group(1)
        mid = HeaderTools._re_id.search(rest)
        if mid:
            d["id_pct"] = mid.group(1)
        mc = HeaderTools._re_cov.search(rest)
        if mc:
            d["cov_pct"] = mc.group(1)

        # Legacy IDs do not carry maf_strand explicitly.
        return d

    @staticmethod
    def sanitize_and_map(
        merged_fa: Path,
        sanitized_fa: Path,
        header_map_tsv: Path,
        query_id: str,
    ) -> Tuple[List[str], Dict[str, Dict[str, str]]]:
        recs = FastaIO.read_fasta(merged_fa)
        seen: Dict[str, int] = {}

        def unique(base: str) -> str:
            seen[base] = seen.get(base, 0) + 1
            return base if seen[base] == 1 else f"{base}|dup={seen[base]}"

        out_records: List[Tuple[str, str]] = []
        desired_order: List[str] = []
        meta_by_sanitized: Dict[str, Dict[str, str]] = {}

        header_map_tsv.parent.mkdir(parents=True, exist_ok=True)
        with header_map_tsv.open("w", encoding="utf-8") as hm:
            hm.write("\t".join([
                "sanitized_id", "original_id", "is_query",
                "assembly_id", "species", "qid",
                "contig", "start_1based", "end_1based",
                "tstrand", "qstrand", "maf_strand",
                "score", "id_pct", "cov_pct", "idcov",
            ]) + "\n")

            for r in recs:
                if r.rid == query_id:
                    # Keep the query header as-is (already pipe-delimited in the new convention).
                    out_records.append((query_id, r.seq))
                    desired_order.append(query_id)

                    _b, contig, s, e, qstrand_id = HeaderTools.parse_query_interval_extended(query_id)
                    meta = {
                        "sanitized_id": query_id,
                        "original_id": query_id,
                        "is_query": "1",
                        "assembly_id": "REF",
                        "species": "",
                        "qid": query_id,
                        "contig": contig,
                        "start_1based": str(s),
                        "end_1based": str(e),
                        "tstrand": qstrand_id,
                        "qstrand": qstrand_id,
                        "maf_strand": qstrand_id,
                        "score": "",
                        "id_pct": "",
                        "cov_pct": "",
                        "idcov": "",
                    }
                    meta_by_sanitized[query_id] = meta
                    hm.write("\t".join(meta[k] for k in [
                        "sanitized_id","original_id","is_query","assembly_id","species","qid",
                        "contig","start_1based","end_1based","tstrand","qstrand","maf_strand",
                        "score","id_pct","cov_pct","idcov"
                    ]) + "\n")
                    continue

                meta_raw = HeaderTools.parse_ortholog_meta(r.rid)
                asm = meta_raw.get("assembly") or "NA"
                contig = meta_raw.get("contig") or "NA"
                s1 = meta_raw.get("start_1based") or "0"
                e1 = meta_raw.get("end_1based") or "0"

                # Determine maf_strand (prefer explicit if present; otherwise infer from alignment strands).
                maf_strand = meta_raw.get("maf_strand") if meta_raw.get("maf_strand") in {"+", "-"} else None
                tstrand = meta_raw.get("tstrand") or "+"
                qstrand = meta_raw.get("qstrand") or "+"
                if maf_strand is None:
                    needs_rc = (tstrand == "-") ^ (qstrand == "-")
                    maf_strand = "-" if needs_rc else "+"

                # Sanitized IDs used for MUSCLE should remain unique and include genomic region info.
                sid_base = f"{asm}|{contig}|{s1}|{e1}|{maf_strand}"
                sid = unique(sid_base)

                out_records.append((sid, r.seq))
                desired_order.append(sid)

                id_pct = meta_raw.get("id_pct") or ""
                cov_pct = meta_raw.get("cov_pct") or ""
                idcov = ""
                try:
                    if id_pct and cov_pct:
                        idcov = f"{float(id_pct) * float(cov_pct):.6f}"
                except ValueError:
                    idcov = ""

                meta = {
                    "sanitized_id": sid,
                    "original_id": r.rid,
                    "is_query": "0",
                    "assembly_id": asm,
                    "species": meta_raw.get("species") or "",
                    "qid": meta_raw.get("qid") or "",
                    "contig": contig,
                    "start_1based": s1,
                    "end_1based": e1,
                    "tstrand": tstrand,
                    "qstrand": qstrand,
                    "maf_strand": maf_strand,
                    "score": meta_raw.get("score") or "",
                    "id_pct": id_pct,
                    "cov_pct": cov_pct,
                    "idcov": idcov,
                }
                meta_by_sanitized[sid] = meta
                hm.write("\t".join(meta[k] for k in [
                    "sanitized_id","original_id","is_query","assembly_id","species","qid",
                    "contig","start_1based","end_1based","tstrand","qstrand","maf_strand",
                    "score","id_pct","cov_pct","idcov"
                ]) + "\n")

        FastaIO.write_fasta(out_records, sanitized_fa)
        return desired_order, meta_by_sanitized

# -----------------------------
# TSS anchoring + alignment utilities
# -----------------------------
class TSSAnchorResolver:
    @staticmethod
    def parse_tss_genomic(tss_genomic: str) -> Tuple[str, int]:
        # Allow arbitrary scaffold/contig names.
        s = tss_genomic.strip()
        m2 = re.match(r"^(.+)[|:](\d+)$", s)
        if not m2:
            raise RuntimeError(
                f"Invalid --tss_genomic format: {tss_genomic} "
                "(expected like contig|pos or contig:pos)"
            )
        return m2.group(1), int(m2.group(2))

    @staticmethod
    def compute_tss_offset_0based(query_id: str, tss_genomic: str, query_seq_len: int) -> int:
        contig, s, _e = HeaderTools.parse_query_interval(query_id)
        tcontig, tpos = TSSAnchorResolver.parse_tss_genomic(tss_genomic)
        if contig != tcontig:
            raise RuntimeError(f"TSS contig mismatch: query={contig}, tss={tcontig}")
        off = tpos - s
        if off < 0 or off >= query_seq_len:
            raise RuntimeError(f"TSS offset out of range: offset={off}, query_len={query_seq_len}")
        return off


class AlnUtils:
    @staticmethod
    def count_ungapped_before(aln_seq: str, col: int) -> int:
        return sum(1 for c in aln_seq[:col] if c != "-")

    @staticmethod
    def build_ref_ungapped_to_col(ref_aln_seq: str) -> Dict[int, int]:
        d: Dict[int, int] = {}
        ref_i = -1
        for col, c in enumerate(ref_aln_seq):
            if c != "-":
                ref_i += 1
                d[ref_i] = col
        return d

    @staticmethod
    def slice_by_columns(aln_seq: str, cols: List[int]) -> str:
        return "".join(aln_seq[c] for c in cols)


class WindowExtractor:
    @staticmethod
    def extract_up_dn(ungapped_seq: str, tss_idx: int, window: int, pad_char: str) -> Tuple[str, str]:
        w = window
        up_start = max(0, tss_idx - w)
        up = ungapped_seq[up_start:tss_idx]
        if len(up) < w:
            up = (pad_char * (w - len(up))) + up
        dn_end = min(len(ungapped_seq), tss_idx + w)
        dn = ungapped_seq[tss_idx:dn_end]
        if len(dn) < w:
            dn = dn + (pad_char * (w - len(dn)))
        return up, dn

    @staticmethod
    def compute_window_coords(
        maf_strand: str,
        frag_start_1based: int,
        frag_end_1based: int,
        tss_idx: int,
        window: int,
        ungapped_len: int,
    ) -> Dict[str, int]:
        w = window
        up_avail = min(w, tss_idx)
        dn_avail = min(w, max(0, ungapped_len - tss_idx))
        up_pad = w - up_avail
        dn_pad = w - dn_avail

        if maf_strand == "+":
            tss_pos = frag_start_1based + tss_idx
        else:
            tss_pos = frag_end_1based - tss_idx

        out: Dict[str, int] = {
            "tss_pos_1based": tss_pos,
            "up_pad": up_pad,
            "dn_pad": dn_pad,
        }

        if up_avail > 0:
            i0 = tss_idx - up_avail
            i1 = tss_idx - 1
            if maf_strand == "+":
                up_start = frag_start_1based + i0
                up_end = frag_start_1based + i1
            else:
                up_start = frag_end_1based - i1
                up_end = frag_end_1based - i0
            out["up_start_1based"] = min(up_start, up_end)
            out["up_end_1based"] = max(up_start, up_end)
        else:
            out["up_start_1based"] = 0
            out["up_end_1based"] = 0

        if dn_avail > 0:
            j0 = tss_idx
            j1 = tss_idx + dn_avail - 1
            if maf_strand == "+":
                dn_start = frag_start_1based + j0
                dn_end = frag_start_1based + j1
            else:
                dn_start = frag_end_1based - j1
                dn_end = frag_end_1based - j0
            out["dn_start_1based"] = min(dn_start, dn_end)
            out["dn_end_1based"] = max(dn_start, dn_end)
        else:
            out["dn_start_1based"] = 0
            out["dn_end_1based"] = 0

        coords = []
        if out["up_start_1based"] and out["up_end_1based"]:
            coords.extend([out["up_start_1based"], out["up_end_1based"]])
        if out["dn_start_1based"] and out["dn_end_1based"]:
            coords.extend([out["dn_start_1based"], out["dn_end_1based"]])
        if coords:
            out["win_start_1based"] = min(coords)
            out["win_end_1based"] = max(coords)
        else:
            out["win_start_1based"] = 0
            out["win_end_1based"] = 0
        return out


# -----------------------------
# Conservation + dumps
# -----------------------------
class ConservationMatrix:
    def __init__(
        self,
        aligned_records: List[FastaRecord],
        query_id: str,
        tss_offset_0based: int,
        threshold: float,
        window: int,
        pad_char: str,
        exclude_query: bool,
    ):
        self.records = aligned_records
        self.query_id = query_id
        self.tss_offset = tss_offset_0based
        self.threshold = threshold
        self.window = window
        self.pad_char = pad_char
        self.exclude_query = exclude_query

        self._validate()
        self.query_rec = self._get_query_record()
        self.ref_ungapped_to_col = AlnUtils.build_ref_ungapped_to_col(self.query_rec.seq)
        self.tss_col = self.ref_ungapped_to_col[self.tss_offset]
        self.conserved_cols = self._compute_conserved_cols()

    def _validate(self) -> None:
        lens = {len(r.seq) for r in self.records}
        if len(lens) != 1:
            raise RuntimeError("Aligned FASTA has inconsistent record lengths.")
        if not (0.0 < self.threshold <= 1.0):
            raise ValueError("--ConservationThreshold must be in (0, 1].")

    def _get_query_record(self) -> FastaRecord:
        for r in self.records:
            if r.rid == self.query_id:
                return r
        raise RuntimeError("Query record not found in alignment.")

    def _compute_conserved_cols(self) -> List[int]:
        ref_seq = self.query_rec.seq
        aln_len = len(ref_seq)
        seqs = [r for r in self.records if (not self.exclude_query) or (r.rid != self.query_id)]
        denom = len(seqs)
        if denom == 0:
            raise RuntimeError("No assembly sequences available for conservation threshold calculation.")
        conserved: List[int] = []
        for col in range(aln_len):
            ref_base = ref_seq[col].upper()
            if ref_base not in DNA_BASES:
                continue
            same = 0
            for r in seqs:
                if r.seq[col].upper() == ref_base:
                    same += 1
            if (same / float(denom)) >= self.threshold:
                conserved.append(col)
        return conserved

    def assembly_tss_idx(self, aln_seq: str) -> int:
        return AlnUtils.count_ungapped_before(aln_seq, self.tss_col)

    def build_rows(self) -> List[Tuple[str, str, str, str, str]]:
        rows: List[Tuple[str, str, str, str, str]] = []
        for rec in self.records:
            ungapped = rec.seq.replace("-", "")
            tss_idx = self.assembly_tss_idx(rec.seq)
            up, dn = WindowExtractor.extract_up_dn(ungapped, tss_idx, self.window, self.pad_char)

            up_pos: List[int] = []
            dn_pos: List[int] = []
            ref_seq = self.query_rec.seq

            for col in self.conserved_cols:
                ref_base = ref_seq[col].upper()
                b = rec.seq[col].upper()
                if b != ref_base or b not in DNA_BASES:
                    continue
                asm_pos = AlnUtils.count_ungapped_before(rec.seq, col)
                rel = asm_pos - tss_idx
                if -self.window <= rel <= -1:
                    up_pos.append(rel + self.window)
                elif 0 <= rel <= (self.window - 1):
                    dn_pos.append(rel)

            rows.append((
                rec.rid,
                up,
                ",".join(map(str, sorted(set(up_pos)))) if up_pos else "",
                dn,
                ",".join(map(str, sorted(set(dn_pos)))) if dn_pos else "",
            ))
        return rows


class MatrixWriter:
    @staticmethod
    def write_conservation_matrix(rows: List[Tuple[str, str, str, str, str]], out_tsv: Path, window: int) -> None:
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        w = int(window)
        with out_tsv.open("w", encoding="utf-8") as out:
            out.write("\t".join([
                "sequence_id",
                f"upstream_{w}bp",
                "upstream_conserved_positions",
                f"downstream_{w}bp",
                "downstream_conserved_positions",
            ]) + "\n")
            for r in rows:
                out.write("\t".join(r) + "\n")

    @staticmethod
    def write_tss_window_coords(
        out_tsv: Path,
        ordered_ids: List[str],
        meta_by_id: Dict[str, Dict[str, str]],
        aligned_by_id: Dict[str, str],
        tss_col: int,
        window: int,
    ) -> None:
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        with out_tsv.open("w", encoding="utf-8") as out:
            out.write("\t".join([
                "sanitized_id",
                "original_id",
                "assembly_id",
                "contig",
                "maf_strand",
                "frag_start_1based",
                "frag_end_1based",
                "tss_col_in_alignment",
                "tss_idx_ungapped",
                "tss_pos_1based",
                "up_start_1based","up_end_1based","up_pad",
                "dn_start_1based","dn_end_1based","dn_pad",
                "win_start_1based","win_end_1based",
            ]) + "\n")

            for sid in ordered_ids:
                meta = meta_by_id[sid]
                aln_seq = aligned_by_id[sid]
                ungapped = aln_seq.replace("-", "")
                tss_idx = AlnUtils.count_ungapped_before(aln_seq, tss_col)

                fs = int(meta["start_1based"]) if meta["start_1based"].isdigit() else 0
                fe = int(meta["end_1based"]) if meta["end_1based"].isdigit() else 0

                coords = WindowExtractor.compute_window_coords(
                    maf_strand=meta["maf_strand"],
                    frag_start_1based=fs,
                    frag_end_1based=fe,
                    tss_idx=tss_idx,
                    window=window,
                    ungapped_len=len(ungapped),
                )

                out.write("\t".join([
                    sid,
                    meta["original_id"],
                    meta["assembly_id"],
                    meta["contig"],
                    meta["maf_strand"],
                    str(fs),
                    str(fe),
                    str(tss_col),
                    str(tss_idx),
                    str(coords["tss_pos_1based"]),
                    str(coords["up_start_1based"]), str(coords["up_end_1based"]), str(coords["up_pad"]),
                    str(coords["dn_start_1based"]), str(coords["dn_end_1based"]), str(coords["dn_pad"]),
                    str(coords["win_start_1based"]), str(coords["win_end_1based"]),
                ]) + "\n")


# -----------------------------
# MAF writer
# -----------------------------
class MafWriter:
    @staticmethod
    def maf_src_ref(ref_name: str, contig: str) -> str:
        return f"{ref_name}.{contig}"

    @staticmethod
    def maf_src_assembly(assembly_id: str, contig: str) -> str:
        return f"{assembly_id}.{contig}"

    @staticmethod
    def write_single_block_maf(out_path: Path, rows: List[Tuple[str, int, int, str, int, str]]) -> None:
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with out_path.open("w", encoding="utf-8") as out:
            out.write("##maf version=1\n")
            out.write("a score=0\n")
            for src, start0, size, strand, srcSize, text in rows:
                out.write(f"s {src} {start0} {size} {strand} {srcSize} {text}\n")
            out.write("\n")


# -----------------------------
# TSS slicing (no extra alignment)
# -----------------------------
class TSSWindowSlicer:
    @staticmethod
    def query_window_columns(query_aln_seq: str, tss_offset_0based: int, window: int) -> List[int]:
        ref_ungapped_to_col = AlnUtils.build_ref_ungapped_to_col(query_aln_seq)
        start_i = tss_offset_0based - window
        end_i = tss_offset_0based + window - 1
        if start_i < 0 or end_i > max(ref_ungapped_to_col.keys()):
            raise RuntimeError("Requested TSS window is out of range of the query fragment/alignment.")
        return [ref_ungapped_to_col[i] for i in range(start_i, end_i + 1)]

    @staticmethod
    def build_aligned_slice_fasta(aligned_records: List[FastaRecord], cols: List[int], out_fa: Path) -> None:
        recs = [(r.rid, AlnUtils.slice_by_columns(r.seq, cols)) for r in aligned_records]
        FastaIO.write_fasta(recs, out_fa)

    @staticmethod
    def build_ungapped_tss_fasta(
        aligned_records: List[FastaRecord],
        tss_col: int,
        window: int,
        pad_char: str,
        out_fa: Path,
        id_override_by_old_id: Optional[Dict[str, str]] = None,
    ) -> None:
        """Write a gap-free, fixed-length TSS-window FASTA.

        By default, record IDs are kept as-is. If `id_override_by_old_id` is provided,
        each output ID is replaced by the corresponding value in the mapping, allowing
        the FASTA headers to encode coordinate systems consistent with other outputs
        (for example, matching the gapless MAF start positions).
        """
        recs: List[Tuple[str, str]] = []
        for r in aligned_records:
            ungapped = r.seq.replace("-", "")
            tss_idx = AlnUtils.count_ungapped_before(r.seq, tss_col)
            up, dn = WindowExtractor.extract_up_dn(ungapped, tss_idx, window, pad_char)
            out_id = r.rid
            if id_override_by_old_id is not None and out_id in id_override_by_old_id:
                out_id = id_override_by_old_id[out_id]
            recs.append((out_id, up + dn))
        FastaIO.write_fasta(recs, out_fa)


# -----------------------------
# Genome FASTA discovery
# -----------------------------
class GenomeFastaIndex:
    @staticmethod
    def build(genome_dir: Path) -> Dict[str, Path]:
        d: Dict[str, Path] = {}
        for p in sorted(genome_dir.iterdir()):
            if p.is_file() and p.suffix.lower() in {".fa", ".fasta", ".fna"}:
                d[p.stem] = p
        return d


# -----------------------------
# Dataset processor
# -----------------------------
class DatasetProcessor:
    def __init__(
        self,
        work_dir: Path,
        runner: CommandRunner,
        muscle_bin: str,
        muscle_extra: Optional[List[str]],
        samtools: str,
        genome_dir: Path,
        ref_genome: Path,
        ref_name: str,
        query_contig: str,
        query_start_1based: int,
        query_end_1based: int,
        ref_contig_size: int,
        tss_pos_1based: int,
        tss_offset_0based: int,
        window: int,
        pad_char: str,
        conservation_threshold: float,
        include_query_in_threshold: bool,
    ):
        self.work_dir = work_dir
        self.runner = runner
        self.muscle = MuscleRunner(muscle_bin=muscle_bin)
        self.muscle_extra = muscle_extra
        self.samtools = samtools

        self.genome_index = GenomeFastaIndex.build(genome_dir)
        self.fai = FaiIndex()
        self.genome_dir = genome_dir
        self.ref_genome = ref_genome

        self.ref_name = ref_name
        self.query_contig = query_contig
        self.query_start_1based = query_start_1based
        self.query_end_1based = query_end_1based
        self.ref_contig_size = ref_contig_size

        self.tss_pos_1based = tss_pos_1based
        self.tss_offset_0based = tss_offset_0based
        self.window = window
        self.pad_char = pad_char
        self.conservation_threshold = conservation_threshold
        self.exclude_query = not include_query_in_threshold

    def process(
        self,
        tag: str,
        merged_like_fasta: Path,
        query_id: str,
        make_conservation_matrix: bool,
        outputs: Dict[str, Path],
    ) -> None:
        desired_order, meta_by_id = HeaderTools.sanitize_and_map(
            merged_fa=merged_like_fasta,
            sanitized_fa=outputs["sanitized_fa"],
            header_map_tsv=outputs["header_map_tsv"],
            query_id=query_id,
        )

        self.muscle.run(in_fa=outputs["sanitized_fa"], out_fa=outputs["aligned_fa"], extra_args=self.muscle_extra)
        reorder_aligned_fasta_inplace(outputs["aligned_fa"], desired_order)

        aligned_records = FastaIO.read_fasta(outputs["aligned_fa"])
        aligned_by_id = {r.rid: r.seq for r in aligned_records}

        cm = ConservationMatrix(
            aligned_records=aligned_records,
            query_id=query_id,
            tss_offset_0based=self.tss_offset_0based,
            threshold=self.conservation_threshold,
            window=self.window,
            pad_char=self.pad_char,
            exclude_query=self.exclude_query,
        )

        if make_conservation_matrix:
            MatrixWriter.write_conservation_matrix(cm.build_rows(), outputs["matrix_tsv"], window=self.window)

        MatrixWriter.write_tss_window_coords(
            out_tsv=outputs["tss_coords_tsv"],
            ordered_ids=desired_order,
            meta_by_id=meta_by_id,
            aligned_by_id=aligned_by_id,
            tss_col=cm.tss_col,
            window=self.window,
        )

        # TSS FASTA (gap-free, fixed length = 2*window with N padding)
        TSSWindowSlicer.build_ungapped_tss_fasta(
            aligned_records=aligned_records,
            tss_col=cm.tss_col,
            window=self.window,
            pad_char=self.pad_char,
            out_fa=outputs["tss_ungapped_fa"],
        )

        # Gap-free TSS-window pseudo-MAF (no '-' in the sequence text).
        if "maf_tss_gapless" in outputs:
            gapless_by_id: Dict[str, str] = {}
            for r in aligned_records:
                ungapped = r.seq.replace("-", "")
                tss_idx = AlnUtils.count_ungapped_before(r.seq, cm.tss_col)
                up, dn = WindowExtractor.extract_up_dn(ungapped, tss_idx, self.window, self.pad_char)
                gapless_by_id[r.rid] = (up + dn).upper()

            query_win_start1 = self.tss_pos_1based - self.window

            maf_rows_gapless: List[Tuple[str, int, int, str, int, str]] = []
            id_override_by_sid: Dict[str, str] = {}
            for sid in desired_order:
                meta = meta_by_id[sid]
                text = gapless_by_id[sid]
                size = len(text)

                if meta["is_query"] == "1":
                    src = MafWriter.maf_src_ref(self.ref_name, self.query_contig)
                    start0 = max(0, query_win_start1 - 1)
                    maf_rows_gapless.append((src, start0, size, "+", self.ref_contig_size, text))
                    # Sync FASTA headers to gapless MAF coordinates (1-based inclusive)
                    start1_fa = start0 + 1
                    end1_fa = start0 + size
                    if "." in src:
                        pfx, c = src.split(".", 1)
                        id_override_by_sid[sid] = f"{pfx}|{c}|{start1_fa}|{end1_fa}|+"
                    else:
                        id_override_by_sid[sid] = f"{src}|{start1_fa}|{end1_fa}|+"
                    continue

                asm = meta["assembly_id"]
                contig = meta["contig"]
                if asm not in self.genome_index:
                    raise RuntimeError(f"Genome FASTA not found for assembly_id={asm} in --genome_dir")

                gfa = self.genome_index[asm]
                gfa_log = self.work_dir / f"logs_faidx_genomes.{tag or 'default'}.log"
                src_size = self.fai.get_len(gfa, contig, self.samtools, self.runner, gfa_log)

                fs = int(meta["start_1based"])
                fe = int(meta["end_1based"])
                tss_idx = AlnUtils.count_ungapped_before(aligned_by_id[sid], cm.tss_col)

                if meta["maf_strand"] == "+":
                    start1 = fs + (tss_idx - self.window)
                else:
                    start1 = fe - (tss_idx + self.window - 1)

                if start1 < 1:
                    start1 = 1
                if start1 > src_size:
                    start1 = src_size

                start0 = start1 - 1
                src = MafWriter.maf_src_assembly(asm, contig)
                maf_rows_gapless.append((src, start0, size, meta["maf_strand"], src_size, text))
                # Sync FASTA headers to gapless MAF coordinates (1-based inclusive)
                start1_fa = start0 + 1
                end1_fa = start0 + size
                strand_fa = meta["maf_strand"]
                if "." in src:
                    pfx, c = src.split(".", 1)
                    id_override_by_sid[sid] = f"{pfx}|{c}|{start1_fa}|{end1_fa}|{strand_fa}"
                else:
                    id_override_by_sid[sid] = f"{src}|{start1_fa}|{end1_fa}|{strand_fa}"

            MafWriter.write_single_block_maf(outputs["maf_tss_gapless"], maf_rows_gapless)

            # Rewrite the gap-free TSS-window FASTA so its header coordinates match the gapless MAF.
            TSSWindowSlicer.build_ungapped_tss_fasta(
                aligned_records=aligned_records,
                tss_col=cm.tss_col,
                window=self.window,
                pad_char=self.pad_char,
                out_fa=outputs["tss_ungapped_fa"],
                id_override_by_old_id=id_override_by_sid,
            )

        # Build aligned-by-query slice FASTA
        cols = TSSWindowSlicer.query_window_columns(
            query_aln_seq=aligned_by_id[query_id],
            tss_offset_0based=self.tss_offset_0based,
            window=self.window,
        )

        TSSWindowSlicer.build_aligned_slice_fasta(
            aligned_records=aligned_records,
            cols=cols,
            out_fa=outputs["tss_aligned_slice_fa"],
        )

        # MAF full (entire MUSCLE alignment)
        maf_rows_full: List[Tuple[str, int, int, str, int, str]] = []
        for r in aligned_records:
            meta = meta_by_id[r.rid]
            aln_text = r.seq
            size = sum(1 for c in aln_text if c != "-")
            if meta["is_query"] == "1":
                src = MafWriter.maf_src_ref(self.ref_name, self.query_contig)
                start0 = self.query_start_1based - 1
                maf_rows_full.append((src, start0, size, "+", self.ref_contig_size, aln_text))
            else:
                asm = meta["assembly_id"]
                contig = meta["contig"]
                if asm not in self.genome_index:
                    raise RuntimeError(f"Genome FASTA not found for assembly_id={asm} in --genome_dir")
                gfa = self.genome_index[asm]
                gfa_log = self.work_dir / f"logs_faidx_genomes.{tag or 'default'}.log"
                src_size = self.fai.get_len(gfa, contig, self.samtools, self.runner, gfa_log)
                start_1 = int(meta["start_1based"])
                start0 = start_1 - 1
                src = MafWriter.maf_src_assembly(asm, contig)
                maf_rows_full.append((src, start0, size, meta["maf_strand"], src_size, aln_text))

        MafWriter.write_single_block_maf(outputs["maf_full"], maf_rows_full)

        # MAF TSS slice (aligned-by-query slice)
        slice_records = FastaIO.read_fasta(outputs["tss_aligned_slice_fa"])
        slice_by_id = {r.rid: r.seq for r in slice_records}

        query_win_start1 = self.tss_pos_1based - self.window
        maf_rows_tss: List[Tuple[str, int, int, str, int, str]] = []

        for sid in desired_order:
            meta = meta_by_id[sid]
            aln_text = slice_by_id[sid]
            size = sum(1 for c in aln_text if c != "-")

            if meta["is_query"] == "1":
                src = MafWriter.maf_src_ref(self.ref_name, self.query_contig)
                start0 = query_win_start1 - 1
                maf_rows_tss.append((src, start0, size, "+", self.ref_contig_size, aln_text))
                continue

            full_aln = aligned_by_id[sid]
            non_gap_cols = [c for c in cols if full_aln[c] != "-"]
            asm = meta["assembly_id"]
            contig = meta["contig"]
            if asm not in self.genome_index:
                raise RuntimeError(f"Genome FASTA not found for assembly_id={asm} in --genome_dir")
            gfa = self.genome_index[asm]
            gfa_log = self.work_dir / f"logs_faidx_genomes.{tag or 'default'}.log"
            src_size = self.fai.get_len(gfa, contig, self.samtools, self.runner, gfa_log)
            src = MafWriter.maf_src_assembly(asm, contig)

            if not non_gap_cols:
                maf_rows_tss.append((src, 0, 0, meta["maf_strand"], src_size, aln_text))
                continue

            first_col = non_gap_cols[0]
            last_col = non_gap_cols[-1]
            i_first = AlnUtils.count_ungapped_before(full_aln, first_col)
            i_last = AlnUtils.count_ungapped_before(full_aln, last_col)

            fs = int(meta["start_1based"])
            fe = int(meta["end_1based"])
            if meta["maf_strand"] == "+":
                start1 = fs + i_first
            else:
                start1 = fe - i_last
            start0 = start1 - 1

            maf_rows_tss.append((src, start0, size, meta["maf_strand"], src_size, aln_text))

        MafWriter.write_single_block_maf(outputs["maf_tss"], maf_rows_tss)


# -----------------------------
# Work directory layout helper
# -----------------------------
@dataclasses.dataclass
class WorkLayout:
    root: Path
    intermediates: Path
    output: Path

    @staticmethod
    def from_root(root: Path) -> "WorkLayout":
        root.mkdir(parents=True, exist_ok=True)
        inter = root / "intermediates"
        out = root / "output"
        inter.mkdir(parents=True, exist_ok=True)
        out.mkdir(parents=True, exist_ok=True)
        return WorkLayout(root=root, intermediates=inter, output=out)


class OutputPublisher:
    @staticmethod
    def copy_to_output(src: Path, dst: Path) -> None:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)


class RefNameInferer:
    @staticmethod
    def infer_from_query_id(query_id: str) -> str:
        # Pipe-delimited: <build>|<contig>|<start>|<end>[|<strand>]
        if "|" in query_id:
            first = query_id.split("|", 1)[0]
            if re.match(r"^([A-Za-z]+\d+)$", first):
                return first
            if re.match(r"^(GRCh\d+)$", first, re.IGNORECASE):
                return first

        # Legacy: hg38_chr..., hg19_chr..., mm10_chr...
        m = re.match(r"^([A-Za-z]+\d+)_", query_id)
        if m:
            return m.group(1)
        # Also accept IDs like "GRCh38_chr..." (rare in some environments)
        m2 = re.match(r"^(GRCh\d+)_", query_id, re.IGNORECASE)
        if m2:
            return m2.group(1)
        return "ref"


# -----------------------------
# Main
# -----------------------------
def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--genome_dir", required=True)
    ap.add_argument("--query_fa", required=True)
    ap.add_argument("--ref_genome", required=True)
    ap.add_argument("--work_dir", required=True)
    ap.add_argument("--lastz", default="lastz")
    ap.add_argument("--samtools", default="samtools")

    ap.add_argument("--threads", type=int, default=8)
    ap.add_argument("--lastz_params", default="T=1 --strand=both --ambiguous=iupac")
    ap.add_argument("--skip_rbh", action="store_true")
    ap.add_argument("--skip_forward_lastz", action="store_true")

    ap.add_argument("--min_id_pct", type=float, default=60.0)
    ap.add_argument("--min_cov_pct", type=float, default=60.0)
    ap.add_argument("--min_score_ratio", type=float, default=1.2)
    ap.add_argument("--min_ref_overlap_bp", type=int, default=100)
    ap.add_argument("--min_ref_overlap_frac", type=float, default=0.2)

    ap.add_argument("--no_progress", action="store_true")
    ap.add_argument("--dry_run", action="store_true")

    ap.add_argument("--tss_genomic", required=True)
    ap.add_argument("--ConservationThreshold", type=float, default=1.0)
    ap.add_argument("--include_query_in_threshold", action="store_true")
    ap.add_argument("--window", type=int, default=150)
    ap.add_argument("--pad_char", default="N")

    ap.add_argument("--muscle", default="muscle")
    ap.add_argument("--muscle_args", default="-maxiters 2 -diags")

    ap.add_argument("--tree_nwk", default="")
    ap.add_argument("--tree_ordered_fasta_name", default="merged_orthologs.tree_ordered.fasta")

    # New: allow explicit MAF reference prefix (defaults to inferred, preserving hg38 behavior)
    ap.add_argument("--ref_name", default="")

    args = ap.parse_args(argv)

    layout = WorkLayout.from_root(Path(args.work_dir))
    work_dir = layout.intermediates
    out_dir = layout.output

    runner = CommandRunner(dry_run=args.dry_run)
    merged_path = work_dir / "merged_orthologs.fasta"

    # RBH -> merged
    if not args.skip_rbh:
        lastz = LastzRunner(lastz_path=Path(args.lastz), runner=runner, extra_params=args.lastz_params)
        extractor = FastaExtractor(runner=runner, samtools=args.samtools)
        rbh = RBHPipeline(
            lastz=lastz,
            extractor=extractor,
            genome_dir=Path(args.genome_dir),
            work_dir=work_dir,
            query_fa=Path(args.query_fa),
            ref_genome=Path(args.ref_genome),
            min_id_pct=args.min_id_pct,
            min_cov_pct=args.min_cov_pct,
            min_score_ratio=args.min_score_ratio,
            min_ref_overlap_bp=args.min_ref_overlap_bp,
            min_ref_overlap_frac=args.min_ref_overlap_frac,
            skip_forward_lastz=args.skip_forward_lastz,
            dry_run=args.dry_run,
        )
        genomes = rbh.list_genomes(Path(args.genome_dir))
        prog = ProgressReporter(total=len(genomes), enabled=(not args.no_progress))
        prog.start(args.threads, "RBH", extra=f"lastz_params='{args.lastz_params}'")

        all_accepted: List[Tuple[str, str]] = []
        with concurrent.futures.ProcessPoolExecutor(max_workers=args.threads) as ex:
            futures = [
                ex.submit(
                    _worker_rbh,
                    {
                        "lastz": args.lastz,
                        "lastz_params": args.lastz_params,
                        "samtools": args.samtools,
                        "genome_dir": args.genome_dir,
                        "work_dir": str(work_dir),
                        "query_fa": args.query_fa,
                        "ref_genome": args.ref_genome,
                        "min_id_pct": args.min_id_pct,
                        "min_cov_pct": args.min_cov_pct,
                        "min_score_ratio": args.min_score_ratio,
                        "min_ref_overlap_bp": args.min_ref_overlap_bp,
                        "min_ref_overlap_frac": args.min_ref_overlap_frac,
                        "skip_forward_lastz": args.skip_forward_lastz,
                        "dry_run": args.dry_run,
                    },
                    str(g),
                )
                for g in genomes
            ]
            for fut in concurrent.futures.as_completed(futures):
                assembly_id, qids_total, fpass, rpass, accepted = fut.result()
                all_accepted.extend(accepted)
                prog.update(f"DONE {assembly_id} (qids={qids_total})", fpass, rpass, len(accepted))

        prog.finish("RBH")

        qrecs = FastaIO.read_fasta(Path(args.query_fa))
        query_id = qrecs[0].rid
        query_seq = qrecs[0].seq

        merged_sorted = OrthologSorter.sort_records(query_id=query_id, records=[(query_id, query_seq)] + all_accepted)
        FastaIO.write_fasta(merged_sorted, merged_path)
    else:
        if not merged_path.exists():
            print(f"ERROR: merged_orthologs.fasta not found: {merged_path}", file=sys.stderr)
            return 2
        merged_recs = FastaIO.read_fasta(merged_path)
        query_id = merged_recs[0].rid
        merged_sorted = OrthologSorter.sort_records(query_id=query_id, records=[(r.rid, r.seq) for r in merged_recs])
        FastaIO.write_fasta(merged_sorted, merged_path)

    merged_recs = FastaIO.read_fasta(merged_path)
    query_id = merged_recs[0].rid
    query_seq = merged_recs[0].seq
    if len(merged_recs) < 2:
        print("ERROR: merged_orthologs.fasta contains only query.", file=sys.stderr)
        return 3

    # shared coordinates from query id
    query_contig, query_start, query_end = HeaderTools.parse_query_interval(query_id)
    tss_contig, tss_pos = TSSAnchorResolver.parse_tss_genomic(args.tss_genomic)
    if tss_contig != query_contig:
        print("ERROR: --tss_genomic contig does not match contig parsed from query ID.", file=sys.stderr)
        return 2

    tss_offset = TSSAnchorResolver.compute_tss_offset_0based(query_id, args.tss_genomic, len(query_seq))

    # Determine reference name prefix in MAF src.
    ref_name = args.ref_name.strip() if args.ref_name.strip() else RefNameInferer.infer_from_query_id(query_id)

    # Resolve srcSize for reference contig from --ref_genome (no longer assumes hg38 naming).
    fai = FaiIndex()
    ref_fai_log = work_dir / "logs_ref_faidx.log"
    ref_contigs = fai.contigs(Path(args.ref_genome), args.samtools, runner, ref_fai_log)
    resolved_ref_contig = RefContigResolver.resolve(query_contig, ref_contigs)
    ref_contig_size = fai.get_len(Path(args.ref_genome), resolved_ref_contig, args.samtools, runner, ref_fai_log)

    # MUSCLE args
    muscle_extra = args.muscle_args.strip().split() if args.muscle_args.strip() else None

    processor = DatasetProcessor(
        work_dir=work_dir,
        runner=runner,
        muscle_bin=args.muscle,
        muscle_extra=muscle_extra,
        samtools=args.samtools,
        genome_dir=Path(args.genome_dir),
        ref_genome=Path(args.ref_genome),
        ref_name=ref_name,
        query_contig=query_contig,
        query_start_1based=query_start,
        query_end_1based=query_end,
        ref_contig_size=ref_contig_size,
        tss_pos_1based=tss_pos,
        tss_offset_0based=tss_offset,
        window=args.window,
        pad_char=args.pad_char,
        conservation_threshold=float(args.ConservationThreshold),
        include_query_in_threshold=bool(args.include_query_in_threshold),
    )

    tag_bp = tss_tag_2xwindow(args.window)

    default_outputs = {
        "sanitized_fa": work_dir / "merged_orthologs.assembly_only.fasta",
        "header_map_tsv": work_dir / "header_map.tsv",
        "aligned_fa": work_dir / "merged_orthologs_aligned.fasta",
        "maf_full": work_dir / "merged_orthologs_aligned.maf",
        "tss_ungapped_fa": work_dir / f"tss_anchored_5pTSS3p_{tag_bp}.fasta",
        "tss_aligned_slice_fa": work_dir / f"tss_anchored_5pTSS3p_{tag_bp}.aligned_by_query.fasta",
        "maf_tss": work_dir / f"tss_anchored_5pTSS3p_{tag_bp}.aligned_by_query.maf",
        "maf_tss_gapless": work_dir / f"tss_anchored_5pTSS3p_{tag_bp}.maf",
        "tss_coords_tsv": work_dir / "tss_window_coords.tsv",
        "matrix_tsv": work_dir / "conservation_matrix.tsv",
    }

    processor.process(
        tag="default",
        merged_like_fasta=merged_path,
        query_id=query_id,
        make_conservation_matrix=True,
        outputs=default_outputs,
    )

    tree_outputs: Optional[Dict[str, Path]] = None
    if args.tree_nwk.strip():
        tree_path = Path(args.tree_nwk)
        if not tree_path.exists():
            print(f"ERROR: tree file not found: {tree_path}", file=sys.stderr)
            return 2
        tips = NewickTipOrder.read_tip_order(tree_path)
        tree_fa = work_dir / args.tree_ordered_fasta_name
        NewickTipOrder.write_tree_ordered_fasta(
            merged_records=[(r.rid, r.seq) for r in merged_recs],
            query_id=query_id,
            tree_tips=tips,
            out_path=tree_fa,
        )

        tree_outputs = {
            "sanitized_fa": work_dir / "merged_orthologs.tree_ordered.assembly_only.fasta",
            "header_map_tsv": work_dir / "header_map.tree_ordered.tsv",
            "aligned_fa": work_dir / "merged_orthologs.tree_ordered_aligned.fasta",
            "maf_full": work_dir / "merged_orthologs.tree_ordered_aligned.maf",
            "tss_ungapped_fa": work_dir / f"tss_anchored_5pTSS3p_{tag_bp}.tree_ordered.fasta",
            "tss_aligned_slice_fa": work_dir / f"tss_anchored_5pTSS3p_{tag_bp}.tree_ordered.aligned_by_query.fasta",
            "maf_tss": work_dir / "merged_orthologs.tree_ordered_tss_slice.maf",
            "maf_tss_gapless": work_dir / f"tss_anchored_5pTSS3p_{tag_bp}.tree_ordered.maf",
            "tss_coords_tsv": work_dir / "tss_window_coords.tree_ordered.tsv",
            "matrix_tsv": work_dir / "conservation_matrix.tree_ordered.tsv",
        }

        processor.process(
            tag="tree_ordered",
            merged_like_fasta=tree_fa,
            query_id=query_id,
            make_conservation_matrix=True,
            outputs=tree_outputs,
        )

    # Publish final files into output/
    try:
        if args.tree_nwk.strip() and tree_outputs is not None:
            publish_items = [
                tree_outputs["matrix_tsv"],
                tree_outputs["aligned_fa"],
                tree_outputs["maf_full"],
                tree_outputs["tss_ungapped_fa"],
                tree_outputs["maf_tss_gapless"],
            ]
        else:
            publish_items = [
                default_outputs["matrix_tsv"],
                default_outputs["aligned_fa"],
                default_outputs["maf_full"],
                default_outputs["tss_ungapped_fa"],
                default_outputs["maf_tss_gapless"],
            ]

        for p in publish_items:
            OutputPublisher.copy_to_output(p, out_dir / p.name)
            print(f"[OUTPUT] {out_dir / p.name}")
    except FileNotFoundError as e:
        print(f"ERROR: failed to publish outputs: {e}", file=sys.stderr)
        return 2

    print(f"[INFO] ref_name={ref_name} query_contig={query_contig} resolved_ref_contig={resolved_ref_contig} ref_contig_size={ref_contig_size}")
    print(f"[DONE] intermediates dir: {work_dir}")
    print(f"[DONE] output dir:        {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
