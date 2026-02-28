#!/usr/bin/env python3
"""
ReAlignPro_tsv2fig.py

Purpose
- Create two PDF sequence-table plots (upstream/downstream from TSS) from a
  precomputed conservation matrix TSV, with optional marker highlighting and
  target motif/homopolymer highlighting.

Key features
1) Window-size independent TSV parsing
   - Auto-detect upstream/downstream sequence columns: upstream_XXXbp / downstream_XXXbp
   - If multiple matches exist, choose the largest XXX
   - Optional overrides: --upstream_seq_col / --downstream_seq_col

2) Labels (assembly/species IDs)
   - Shown to the left of sequences, right-aligned
   - Exactly 1 base-column gap before the first base
   - If sequence_id contains '|', only the part before the first '|' is used for:
       (a) display label
       (b) marker matching
     Example: "m11_Homo_sapiens|CP068262.2|6048623|6049722|+" -> "m11_Homo_sapiens"

3) Marker highlighting (optional)
   - Enable with: --marker marker_matrix.tsv
   - Marker matrix format (4 columns): assembly_id, template_seq, marker_5p, marker_3p
   - Rows with BOTH marker_5p and marker_3p empty are ignored
   - If one marker is empty, only the other is searched/highlighted
   - Markers are colored black with REGULAR weight (not bold)

4) Conservation shading
   - Default (no --marker):
       conserved positions: dimgray
       non-conserved: lightgrey
   - With --marker:
       conserved positions: grey
       non-conserved: lightgrey

5) Target highlight (motif / homopolymer)
   - Always drawn on top (overrides marker/conservation colors)
   - Base colors: A=red, T=fuchsia, C=blue, G=turquoise
   - Bold fontweight

Filtering behavior
- Default: NO assemblies are omitted, even when --marker is provided.
- Optional: --filter_unmarked plots ONLY rows that have a matched marker entry.
- Backward compatibility: --keep_unmarked disables filtering even if --filter_unmarked is set.

Examples
  python ReAlignPro_tsv2fig_v6.py conservation_matrix.tree_ordered.tsv --type variant --highlight CCG
  python ReAlignPro_tsv2fig_v6.py conservation_matrix.tree_ordered.tsv --type variant --highlight CCG --marker marker_matrix.tsv
  python ReAlignPro_tsv2fig_v6.py conservation_matrix.tree_ordered.tsv --type variant --highlight CCG --marker marker_matrix.tsv --filter_unmarked
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

import matplotlib.pyplot as plt


# -----------------------------
# Defaults and styling
# -----------------------------
DEFAULT_VARIANT_TYPE = "variant"
DEFAULT_HIGHLIGHT = "CCG"
DEFAULT_THRESHOLD = 10
DEFAULT_FONT_FAMILY = "monospace"
DEFAULT_FONTSIZE = 7
DEFAULT_DPI = 300

LIGHT_GREY = "lightgrey"
DARK_GREY = "dimgray"   # marker mode OFF
GREY = "grey"           # marker mode ON
MARKER_COLOR = "black"

# Base colors preserved from legacy ReAlignMAF_seqPlot.py
BASE_COLOR_MAP: Dict[str, str] = {
    "A": "red",
    "T": "fuchsia",
    "C": "blue",
    "G": "turquoise",
}

UPSTREAM_SEQ_COL_RE = re.compile(r"^upstream_(\d+)bp$", re.IGNORECASE)
DOWNSTREAM_SEQ_COL_RE = re.compile(r"^downstream_(\d+)bp$", re.IGNORECASE)


def eprint(msg: str) -> None:
    """Print to stderr."""
    print(msg, file=sys.stderr)


def normalize_sequence_id(seq_id: str) -> str:
    """If '|' exists, keep only the part before the first '|'."""
    return seq_id.split("|", 1)[0] if "|" in seq_id else seq_id


# -----------------------------
# Data models
# -----------------------------
@dataclass(frozen=True)
class SequenceRecord:
    sequence_id: str
    upstream_seq: str
    upstream_conserved: Set[int]
    downstream_seq: str
    downstream_conserved: Set[int]


@dataclass(frozen=True)
class MarkerEntry:
    assembly_id: str
    template_seq: str
    marker_5p: str
    marker_3p: str

    def has_any_marker(self) -> bool:
        return bool(self.marker_5p) or bool(self.marker_3p)


# -----------------------------
# Column detection
# -----------------------------
class ColumnDetector:
    """Detect upstream/downstream sequence columns such as upstream_100bp."""

    def detect(
        self,
        fieldnames: List[str],
        upstream_override: str = "",
        downstream_override: str = "",
    ) -> Tuple[str, str]:
        up = upstream_override.strip()
        dn = downstream_override.strip()

        if up:
            if up not in fieldnames:
                raise ValueError(f"--upstream_seq_col '{up}' not found in TSV header.")
        else:
            up = self._pick_largest(fieldnames, UPSTREAM_SEQ_COL_RE, label="upstream")

        if dn:
            if dn not in fieldnames:
                raise ValueError(f"--downstream_seq_col '{dn}' not found in TSV header.")
        else:
            dn = self._pick_largest(fieldnames, DOWNSTREAM_SEQ_COL_RE, label="downstream")

        return up, dn

    @staticmethod
    def _pick_largest(fieldnames: List[str], pattern: re.Pattern, label: str) -> str:
        candidates: List[Tuple[int, str]] = []
        for fn in fieldnames:
            m = pattern.match(fn)
            if m:
                candidates.append((int(m.group(1)), fn))

        if not candidates:
            raise ValueError(
                f"TSV header is missing a {label} sequence column matching '{pattern.pattern}'.\n"
                f"Found columns: {fieldnames}"
            )

        candidates.sort(key=lambda x: x[0], reverse=True)
        return candidates[0][1]


# -----------------------------
# TSV reading
# -----------------------------
class ConservationMatrixReader:
    FIXED_REQUIRED = ["sequence_id", "upstream_conserved_positions", "downstream_conserved_positions"]

    def __init__(self, upstream_seq_col: str = "", downstream_seq_col: str = "") -> None:
        self.detector = ColumnDetector()
        self.upstream_override = upstream_seq_col
        self.downstream_override = downstream_seq_col

    def read(self, tsv_path: Path) -> Tuple[List[SequenceRecord], str, str]:
        if not tsv_path.exists():
            raise FileNotFoundError(f"Input TSV not found: {tsv_path}")

        out: List[SequenceRecord] = []
        with tsv_path.open("r", newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            if reader.fieldnames is None:
                raise ValueError("TSV appears to have no header row.")

            missing = [c for c in self.FIXED_REQUIRED if c not in reader.fieldnames]
            if missing:
                raise ValueError(
                    f"TSV header is missing required columns: {', '.join(missing)}\nFound columns: {reader.fieldnames}"
                )

            up_col, dn_col = self.detector.detect(reader.fieldnames, self.upstream_override, self.downstream_override)
            eprint(f"[INFO] Using upstream sequence column: {up_col}")
            eprint(f"[INFO] Using downstream sequence column: {dn_col}")

            for line_no, row in enumerate(reader, start=2):
                sid = (row.get("sequence_id") or "").strip()
                up_seq = (row.get(up_col) or "").strip()
                dn_seq = (row.get(dn_col) or "").strip()
                up_pos = (row.get("upstream_conserved_positions") or "").strip()
                dn_pos = (row.get("downstream_conserved_positions") or "").strip()

                if not sid:
                    raise ValueError(f"Empty sequence_id at line {line_no}.")
                if not up_seq or not dn_seq:
                    raise ValueError(f"Empty upstream/downstream sequence at line {line_no} (sequence_id={sid}).")

                out.append(
                    SequenceRecord(
                        sequence_id=sid,
                        upstream_seq=up_seq,
                        upstream_conserved=self._parse_positions(up_pos),
                        downstream_seq=dn_seq,
                        downstream_conserved=self._parse_positions(dn_pos),
                    )
                )

        if not out:
            raise ValueError("No records parsed from TSV.")

        self._validate_lengths(out)
        self._validate_positions(out)
        return out, up_col, dn_col

    @staticmethod
    def _parse_positions(s: str) -> Set[int]:
        if s == "" or s.upper() == "NA" or s == ".":
            return set()
        out: Set[int] = set()
        for tok in s.split(","):
            tok = tok.strip()
            if tok:
                out.add(int(tok))
        return out

    @staticmethod
    def _validate_lengths(records: List[SequenceRecord]) -> None:
        up_len = len(records[0].upstream_seq)
        dn_len = len(records[0].downstream_seq)
        for r in records:
            if len(r.upstream_seq) != up_len:
                raise ValueError("Upstream sequences are not all the same length.")
            if len(r.downstream_seq) != dn_len:
                raise ValueError("Downstream sequences are not all the same length.")

    @staticmethod
    def _validate_positions(records: List[SequenceRecord]) -> None:
        up_len = len(records[0].upstream_seq)
        dn_len = len(records[0].downstream_seq)
        for r in records:
            if any(i < 0 or i >= up_len for i in r.upstream_conserved):
                raise ValueError(f"Upstream conserved_positions out of bounds for {r.sequence_id} (len={up_len}).")
            if any(i < 0 or i >= dn_len for i in r.downstream_conserved):
                raise ValueError(f"Downstream conserved_positions out of bounds for {r.sequence_id} (len={dn_len}).")


# -----------------------------
# Marker reading and matching
# -----------------------------
class MarkerMatrixReader:
    """
    Read marker matrix file with 4 columns:
      assembly_id, template_seq, marker_5p, marker_3p

    Supported delimiters:
      - tab
      - comma
      - whitespace (fallback)

    Exclusion rule:
      - Rows with BOTH marker columns empty are excluded
      - If one marker is empty, the other one is still used
    """

    def read(self, path: Path) -> List[MarkerEntry]:
        if not path.exists():
            raise FileNotFoundError(f"Marker matrix file not found: {path}")

        lines = [ln.rstrip("\n") for ln in path.open("r")]
        if not lines:
            raise ValueError(f"Marker matrix is empty: {path}")

        rows = self._parse_rows(lines)
        if rows and self._looks_like_header(rows[0]):
            rows = rows[1:]

        entries: List[MarkerEntry] = []
        for row in rows:
            if len(row) < 4:
                continue

            assembly_raw = row[0].strip()
            template = row[1].strip()
            m5 = row[2].strip()
            m3 = row[3].strip()

            if not assembly_raw or not template:
                continue

            assembly = normalize_sequence_id(assembly_raw)
            entry = MarkerEntry(assembly_id=assembly, template_seq=template, marker_5p=m5, marker_3p=m3)
            if entry.has_any_marker():
                entries.append(entry)

        if not entries:
            raise ValueError("No usable marker rows found (need at least one non-empty marker column).")
        return entries

    @staticmethod
    def _parse_rows(lines: List[str]) -> List[List[str]]:
        tab_rows = [ln.split("\t") for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
        if tab_rows and max(len(r) for r in tab_rows) >= 4:
            return tab_rows

        comma_rows = [ln.split(",") for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]
        if comma_rows and max(len(r) for r in comma_rows) >= 4:
            return comma_rows

        return [re.split(r"\s+", ln.strip()) for ln in lines if ln.strip() and not ln.lstrip().startswith("#")]

    @staticmethod
    def _looks_like_header(row: List[str]) -> bool:
        joined = " ".join(row).lower()
        return any(k in joined for k in ["assembly", "sequence_id", "template", "marker", "5marker", "3marker"])


class MarkerMatcher:
    """
    Match marker entries to (assembly_id, region_sequence).

    1) Exact match on (assembly_id AND template_seq == region_seq) (case-insensitive)
    2) Fallback: if there is exactly one entry for assembly_id, use it
    3) Otherwise: None
    """

    def __init__(self, entries: List[MarkerEntry]) -> None:
        self.by_assembly: Dict[str, List[MarkerEntry]] = {}
        for e in entries:
            self.by_assembly.setdefault(e.assembly_id, []).append(e)

    def find(self, assembly_id: str, region_seq: str) -> Optional[MarkerEntry]:
        cands = self.by_assembly.get(assembly_id)
        if not cands:
            return None

        region_u = region_seq.upper()
        for e in cands:
            if e.template_seq.upper() == region_u:
                return e

        if len(cands) == 1:
            return cands[0]
        return None


# -----------------------------
# Target highlighting (motif/homopolymer)
# -----------------------------
class Highlighter:
    """Target motif/homopolymer highlighting (drawn on top)."""

    def __init__(self, mode: str, motif: str, threshold: int) -> None:
        if mode not in ("variant", "homopolymer"):
            raise ValueError("mode must be 'variant' or 'homopolymer'.")
        self.mode = mode
        self.motif = (motif or "").upper()
        self.threshold = int(threshold)

        if self.mode == "variant" and not self.motif:
            raise ValueError("In variant mode, --highlight must be provided.")
        if self.mode == "homopolymer" and self.threshold <= 0:
            raise ValueError("In homopolymer mode, --threshold must be > 0.")

    def indices(self, seq: str) -> Set[int]:
        seq_u = seq.upper()
        if self.mode == "variant":
            return self._motif_indices(seq_u, self.motif)
        return self._homopolymer_indices(seq_u, self.threshold)

    @staticmethod
    def _motif_indices(seq: str, motif: str) -> Set[int]:
        idxs: Set[int] = set()
        m = len(motif)
        start = 0
        while True:
            start = seq.find(motif, start)
            if start == -1:
                break
            idxs.update(range(start, start + m))
            start += 1
        return idxs

    @staticmethod
    def _homopolymer_indices(seq: str, threshold: int) -> Set[int]:
        idxs: Set[int] = set()
        i = 0
        n = len(seq)
        while i < n:
            b = seq[i]
            if b not in "ATCG":
                i += 1
                continue
            j = i + 1
            while j < n and seq[j] == b:
                j += 1
            if (j - i) >= threshold:
                idxs.update(range(i, j))
            i = j
        return idxs


# -----------------------------
# Marker highlighting
# -----------------------------
class MarkerHighlighter:
    """Find indices for marker sequences in a template sequence (supports empty markers)."""

    @staticmethod
    def _indices_for_one(seq: str, marker: str) -> Set[int]:
        marker_u = marker.upper().strip()
        if not marker_u:
            return set()

        seq_u = seq.upper()
        idxs: Set[int] = set()
        m = len(marker_u)
        start = 0
        while True:
            start = seq_u.find(marker_u, start)
            if start == -1:
                break
            idxs.update(range(start, start + m))
            start += 1
        return idxs

    def indices(self, seq: str, entry: MarkerEntry) -> Set[int]:
        idxs = set()
        idxs |= self._indices_for_one(seq, entry.marker_5p)
        idxs |= self._indices_for_one(seq, entry.marker_3p)
        return idxs


# -----------------------------
# Plotting
# -----------------------------
class SequenceTablePlotter:
    def __init__(self, font: str, fontsize: int) -> None:
        self.font = font
        self.fontsize = fontsize
        self.marker_highlighter = MarkerHighlighter()

    @staticmethod
    def _figsize(n_rows: int, n_cols: int) -> Tuple[float, float]:
        width = max(2.0, n_cols / 12.0)
        height = max(2.0, n_rows / 8.0)
        return width, height

    @staticmethod
    def _compute_left_margin(display_ids: List[str]) -> int:
        max_len = max((len(s) for s in display_ids), default=0)
        return max_len + 1  # 1 base-column gap

    def plot(
        self,
        sequences: List[str],
        conserved_sets: List[Set[int]],
        display_ids: List[str],
        title: str,
        out_pdf: Path,
        target_highlighter: Highlighter,
        marker_matcher: Optional[MarkerMatcher],
        filter_unmarked: bool,
        dpi: int,
    ) -> None:
        plt.rcParams["font.family"] = self.font

        # Build marker entries list (matching uses normalized display_ids)
        marker_entries: List[Optional[MarkerEntry]] = []
        if marker_matcher is None:
            marker_entries = [None] * len(sequences)
        else:
            for sid, seq in zip(display_ids, sequences):
                marker_entries.append(marker_matcher.find(sid, seq))

        # Optional filtering (explicit only)
        if marker_matcher is not None and filter_unmarked:
            filtered = []
            for sid, seq, cons, mentry in zip(display_ids, sequences, conserved_sets, marker_entries):
                if mentry is None or not mentry.has_any_marker():
                    continue
                filtered.append((sid, seq, cons, mentry))
            if not filtered:
                raise ValueError("No sequences matched marker entries after --filter_unmarked.")
            display_ids = [t[0] for t in filtered]
            sequences = [t[1] for t in filtered]
            conserved_sets = [t[2] for t in filtered]
            marker_entries = [t[3] for t in filtered]

        conserved_color = GREY if marker_matcher is not None else DARK_GREY

        n_rows = len(sequences)
        n_cols_seq = len(sequences[0])

        left_margin = self._compute_left_margin(display_ids)
        label_x = left_margin - 1
        seq_start_x = left_margin
        n_cols_total = left_margin + n_cols_seq

        fig = plt.figure(figsize=self._figsize(n_rows, n_cols_total))

        for row_idx, (sid, seq, cons, mentry) in enumerate(zip(display_ids, sequences, conserved_sets, marker_entries)):
            # 1) conservation shading
            colors = [conserved_color if i in cons else LIGHT_GREY for i in range(len(seq))]
            weights = ["regular"] * len(seq)

            # 2) marker highlight: black, regular
            if mentry is not None and mentry.has_any_marker():
                for i in self.marker_highlighter.indices(seq, mentry):
                    colors[i] = MARKER_COLOR
                    weights[i] = "regular"

            # 3) target highlight on top: base colors + bold
            for i in target_highlighter.indices(seq):
                base = seq[i].upper()
                colors[i] = BASE_COLOR_MAP.get(base, "black")
                weights[i] = "bold"

            y = n_rows - row_idx

            # left label (right-aligned), 1-column gap before sequence
            plt.text(
                label_x,
                y,
                sid,
                fontsize=self.fontsize,
                ha="right",
                va="center",
                color="black",
                fontweight="regular",
            )

            # sequence (shifted)
            for x, base in enumerate(seq):
                plt.text(
                    seq_start_x + x,
                    y,
                    base,
                    fontsize=self.fontsize,
                    ha="center",
                    va="center",
                    color=colors[x],
                    fontweight=weights[x],
                )

        plt.xlim(-0.5, n_cols_total - 0.5)
        plt.ylim(0, n_rows + 1)
        plt.xticks([])
        plt.yticks([])
        plt.xlabel("DNA Sequence Table")
        plt.title(title)

        out_pdf.parent.mkdir(parents=True, exist_ok=True)
        plt.savefig(str(out_pdf), format="pdf", bbox_inches="tight", dpi=dpi)
        plt.close(fig)


# -----------------------------
# CLI app
# -----------------------------
class App:
    def __init__(self, args: argparse.Namespace) -> None:
        self.args = args
        self.reader = ConservationMatrixReader(args.upstream_seq_col, args.downstream_seq_col)
        self.plotter = SequenceTablePlotter(args.font, args.fontsize)
        self.target_highlighter = Highlighter(args.type, args.highlight, args.threshold)

        self.marker_matcher: Optional[MarkerMatcher] = None
        if args.marker:
            entries = MarkerMatrixReader().read(Path(args.marker))
            self.marker_matcher = MarkerMatcher(entries)
            eprint(f"[INFO] Loaded marker entries: {len(entries)} (rows with both markers empty were excluded).")

    def run(self) -> None:
        tsv = Path(self.args.input).resolve()
        outdir = Path(self.args.outdir).resolve() if self.args.outdir else tsv.parent.resolve()
        prefix = self.args.prefix or tsv.stem

        records, _, _ = self.reader.read(tsv)

        up_pdf = Path(self.args.upstream_pdf).resolve() if self.args.upstream_pdf else outdir / f"{prefix}.upstream.pdf"
        dn_pdf = Path(self.args.downstream_pdf).resolve() if self.args.downstream_pdf else outdir / f"{prefix}.downstream.pdf"

        # Normalized IDs for display and marker matching
        display_ids = [normalize_sequence_id(r.sequence_id) for r in records]

        # Backward compatible behavior:
        # - Default keeps all rows
        # - --filter_unmarked filters rows
        # - --keep_unmarked disables filtering even if --filter_unmarked is set
        effective_filter = bool(self.args.filter_unmarked) and (not bool(self.args.keep_unmarked))

        self.plotter.plot(
            sequences=[r.upstream_seq for r in records],
            conserved_sets=[r.upstream_conserved for r in records],
            display_ids=display_ids,
            title=self.args.title_upstream or "ReAlignPro plot (upstream from TSS)",
            out_pdf=up_pdf,
            target_highlighter=self.target_highlighter,
            marker_matcher=self.marker_matcher,
            filter_unmarked=effective_filter,
            dpi=self.args.dpi,
        )
        self.plotter.plot(
            sequences=[r.downstream_seq for r in records],
            conserved_sets=[r.downstream_conserved for r in records],
            display_ids=display_ids,
            title=self.args.title_downstream or "ReAlignPro plot (downstream from TSS)",
            out_pdf=dn_pdf,
            target_highlighter=self.target_highlighter,
            marker_matcher=self.marker_matcher,
            filter_unmarked=effective_filter,
            dpi=self.args.dpi,
        )

        eprint(f"[DONE] Wrote: {up_pdf}")
        eprint(f"[DONE] Wrote: {dn_pdf}")


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Plot ReAlignPro sequence tables from conservation TSV, with optional marker emphasis and simplified assembly labels."
    )
    p.add_argument("input", help="Input TSV (conservation_matrix.tree_ordered.tsv)")

    p.add_argument("--type", choices=["variant", "homopolymer"], default=DEFAULT_VARIANT_TYPE)
    p.add_argument("--highlight", default=DEFAULT_HIGHLIGHT)
    p.add_argument("--threshold", type=int, default=DEFAULT_THRESHOLD)

    p.add_argument("--outdir", default="")
    p.add_argument("--prefix", default="")
    p.add_argument("--upstream_pdf", default="")
    p.add_argument("--downstream_pdf", default="")

    p.add_argument("--title_upstream", default="")
    p.add_argument("--title_downstream", default="")

    p.add_argument("--upstream_seq_col", default="", help="Override upstream sequence column name.")
    p.add_argument("--downstream_seq_col", default="", help="Override downstream sequence column name.")

    p.add_argument("--marker", default="", help="Marker matrix file (4 columns: assembly_id, template_seq, 5p_marker, 3p_marker).")
    p.add_argument("--filter_unmarked", action="store_true",
                   help="When --marker is provided, plot ONLY sequences with matched marker info (default: off).")

    # Backward compatibility with older scripts
    p.add_argument("--keep_unmarked", action="store_true",
                   help="Backward-compatible flag: disables filtering even if --filter_unmarked is set.")

    p.add_argument("--dpi", type=int, default=DEFAULT_DPI)
    p.add_argument("--fontsize", type=int, default=DEFAULT_FONTSIZE)
    p.add_argument("--font", default=DEFAULT_FONT_FAMILY)
    return p


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    App(args).run()
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
