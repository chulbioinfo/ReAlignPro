from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

import pytest

from realignpro.maf2con_cov import call_major_base_cov, matrix2con_cov, scan_expected_depth

SRC = Path(__file__).resolve().parents[1] / "src"


# --------------------------------------------------------------------------- #
# unit: coverage-aware major-allele call
# --------------------------------------------------------------------------- #
def test_call_major_base_cov_denominator_is_ncov() -> None:
    # gap/N excluded from BOTH numerator and denominator: 99 A + 1 C (+ gaps) -> 99/100.
    base, ratio, n_cov = call_major_base_cov(["A"] * 99 + ["C"] + ["-"] * 900, min_major_similarity=0.99)
    assert n_cov == 100
    assert base is None  # 0.99 is not > 0.99
    assert ratio == pytest.approx(0.99)

    base, ratio, n_cov = call_major_base_cov(["A"] * 100 + ["C"] + ["N"] * 50, min_major_similarity=0.99)
    assert n_cov == 101
    assert base == "A"
    assert ratio == pytest.approx(100 / 101)


def test_call_major_base_cov_min_depth_floor() -> None:
    base, _ratio, n_cov = call_major_base_cov(["A", "A", "A"], min_major_similarity=0.99, min_depth=5)
    assert n_cov == 3
    assert base is None  # below coverage floor


def test_call_major_base_cov_fixed_only() -> None:
    col = ["A"] * 199 + ["C"]  # 0.995: passes default threshold, but NOT 100% fixed
    assert call_major_base_cov(col, 0.99, 1, fixed_only=False)[0] == "A"
    assert call_major_base_cov(col, 0.99, 1, fixed_only=True)[0] is None

    base, ratio, n_cov = call_major_base_cov(["A"] * 10, 0.99, 1, fixed_only=True)
    assert base == "A" and ratio == pytest.approx(1.0) and n_cov == 10

    base, _r, n_cov = call_major_base_cov(["A"] * 5 + ["-", "N"], 0.99, 1, fixed_only=True)
    assert base == "A" and n_cov == 5


# --------------------------------------------------------------------------- #
# unit: block-level, no whole-block gate
# --------------------------------------------------------------------------- #
def _block(seqs: dict) -> dict:
    return {
        sid: {"seq": s, "chr": "chr1", "start": 0, "size": len(s.replace("-", "")),
              "strand": "+", "srcSize": 1000}
        for sid, s in seqs.items()
    }


def test_matrix2con_cov_no_whole_block_gate() -> None:
    # Only 2 of the 3 declared targets are present; the block is still assessed.
    block = _block({"hg38": "ACGT", "hs1": "ACGT"})
    out = matrix2con_cov(block, "hg38", ["hg38", "hs1", "hs2"], [], min_depth=2)
    assert out == ["chr1\t0\t4\n"]


def test_matrix2con_cov_call_rate_floor_drops_low_coverage() -> None:
    # N_exp=3, call-rate 0.99 -> floor ceil(2.97)=3; a 2-of-3 column is dropped.
    block = _block({"hg38": "ACGT", "hs1": "ACGT"})
    out = matrix2con_cov(block, "hg38", ["hg38", "hs1", "hs2"], [], min_depth=1,
                         min_call_rate=0.99, expected_depth_global=3)
    assert out == []


# --------------------------------------------------------------------------- #
# scan: per-chromosome N_exp
# --------------------------------------------------------------------------- #
def test_scan_expected_depth_per_chrom(tmp_path: Path) -> None:
    maf = tmp_path / "s.maf"
    maf.write_text(
        "##maf version=1\n\n"
        "a\ns hg38.chr1 0 2 + 9 AC\ns hs1.chr1 0 2 + 9 AC\ns hs2.chr1 0 2 + 9 AC\n\n"
        "a\ns hg38.chr1 5 2 + 9 AC\ns hs1.chr1 5 2 + 9 AC\n\n"
        "a\ns hg38.chrY 0 2 + 9 AC\ns hs1.chrY 0 2 + 9 AC\n\n",
        encoding="utf-8",
    )
    n_exp = scan_expected_depth(str(maf), "hg38", {"hg38", "hs1", "hs2"}, set())
    assert n_exp == {"chr1": 3, "chrY": 2}


# --------------------------------------------------------------------------- #
# e2e (default engine is coverage-aware; call-rate 0.99 default)
# --------------------------------------------------------------------------- #
def _run(args, env_src: Path) -> None:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(env_src) + os.pathsep + env.get("PYTHONPATH", "")
    subprocess.run([sys.executable, "-m", "realignpro", "maf2con", *args],
                   check=True, timeout=60, env=env)


def _two_block_maf(path: Path) -> None:
    # block1: hg38+hs1+hs2 present (conserved); block2: hs2 MISSING (conserved).
    path.write_text(
        "##maf version=1\n\n"
        "a score=0\ns hg38.chr1 0 4 + 1000 ACGT\ns hs1.chr1 0 4 + 1000 ACGT\ns hs2.chr1 0 4 + 1000 ACGT\n\n"
        "a score=0\ns hg38.chr1 10 4 + 1000 ACGT\ns hs1.chr1 10 4 + 1000 ACGT\n\n",
        encoding="utf-8",
    )


def test_e2e_no_gate_and_call_rate_floor(tmp_path: Path) -> None:
    maf = tmp_path / "two.maf"
    _two_block_maf(maf)

    # --min-call-rate 0: no relative floor AND no whole-block gate -> block2 kept.
    bed0 = tmp_path / "cr0.bed"
    _run(["--input", str(maf), "--output", str(bed0), "-t", "3", "--ref-id", "hg38",
          "--target-ids", "all", "--min-call-rate", "0", "--min-depth", "2"], SRC)
    assert bed0.read_text(encoding="utf-8") == "chr1\t0\t4\nchr1\t10\t14\n"

    # default call-rate 0.99: N_exp=3 -> floor 3; block2 (N_cov=2) dropped.
    bed99 = tmp_path / "cr99.bed"
    _run(["--input", str(maf), "--output", str(bed99), "-t", "3", "--ref-id", "hg38",
          "--target-ids", "all"], SRC)
    assert bed99.read_text(encoding="utf-8") == "chr1\t0\t4\n"


def _write_1000_haplotype_maf(path: Path) -> None:
    def alt(b: str) -> str:
        return {"A": "C", "C": "G", "G": "T", "T": "A"}[b]

    ref_seq = "ACGT" * 25
    ids = ["hg38"] + [f"hs{i}" for i in range(1, 1000)]
    lines = ["##maf version=1\n\n", "a score=0\n"]
    for idx, sid in enumerate(ids):
        seq = list(ref_seq)
        if 1 <= idx <= 10:
            seq[20] = alt(seq[20])
        if 1 <= idx <= 9:
            seq[50] = alt(seq[50])
        lines.append(f"s {sid}.chr1 0 100 + 100 {''.join(seq)}\n")
    lines.append("\n")
    path.write_text("".join(lines), encoding="utf-8")


def test_e2e_all_present_conservation_only(tmp_path: Path) -> None:
    # --min-call-rate 0 isolates the conservation test; all present -> pos20 out, pos50 in.
    maf = tmp_path / "h1000.maf"
    _write_1000_haplotype_maf(maf)
    bed = tmp_path / "h1000.bed"
    _run(["--input", str(maf), "--output", str(bed), "-t", "3", "--ref-id", "hg38",
          "--target-ids", "all", "--min-call-rate", "0", "--min-depth", "1"], SRC)
    assert bed.read_text(encoding="utf-8") == "chr1\t0\t20\nchr1\t21\t100\n"


def test_e2e_emit_depth_columns(tmp_path: Path) -> None:
    maf = tmp_path / "two.maf"
    _two_block_maf(maf)
    bed = tmp_path / "emit.bed"
    _run(["--input", str(maf), "--output", str(bed), "-t", "3", "--ref-id", "hg38",
          "--target-ids", "all", "--min-call-rate", "0", "--min-depth", "2", "--emit-depth"], SRC)
    assert bed.read_text(encoding="utf-8") == "chr1\t0\t4\t3\t1.0000\nchr1\t10\t14\t2\t1.0000\n"


def test_e2e_fixed_only_excludes_non_monomorphic_columns(tmp_path: Path) -> None:
    # 4 assemblies; ref ACGT. col1 has one mismatch (C,C,C,A -> 3/4); cols 0,2,3 fixed.
    maf = tmp_path / "fixed.maf"
    maf.write_text(
        "##maf version=1\n\n"
        "a score=0\n"
        "s hg38.chr1 0 4 + 100 ACGT\n"
        "s hs1.chr1 0 4 + 100 ACGT\n"
        "s hs2.chr1 0 4 + 100 ACGT\n"
        "s hs3.chr1 0 4 + 100 AAGT\n\n",
        encoding="utf-8",
    )
    # loose threshold 0.5, no relative floor: col1 major C at 0.75 passes -> whole [0,4)
    bed_a = tmp_path / "loose.bed"
    _run(["--input", str(maf), "--output", str(bed_a), "-t", "3", "--ref-id", "hg38",
          "--target-ids", "all", "--min-call-rate", "0", "--min-depth", "2",
          "--min-major-similarity", "0.5"], SRC)
    assert bed_a.read_text(encoding="utf-8") == "chr1\t0\t4\n"

    # fixed-only: col1 not monomorphic -> split into [0,1) and [2,4)
    bed_b = tmp_path / "fixed.bed"
    _run(["--input", str(maf), "--output", str(bed_b), "-t", "3", "--ref-id", "hg38",
          "--target-ids", "all", "--min-call-rate", "0", "--min-depth", "2", "--fixed-only"], SRC)
    assert bed_b.read_text(encoding="utf-8") == "chr1\t0\t1\nchr1\t2\t4\n"


def test_e2e_auto_depth_call_rate_floor(tmp_path: Path) -> None:
    # block1 depth 10 (N_exp), block2 depth 3.
    lines = ["##maf version=1\n\n", "a score=0\n", "s hg38.chr1 0 2 + 100 AC\n"]
    lines += [f"s hs{i}.chr1 0 2 + 100 AC\n" for i in range(1, 10)]  # 9 + ref = 10
    lines.append("\n")
    lines += ["a score=0\ns hg38.chr1 20 2 + 100 AC\ns hs1.chr1 20 2 + 100 AC\ns hs2.chr1 20 2 + 100 AC\n\n"]
    maf = tmp_path / "depth.maf"
    maf.write_text("".join(lines), encoding="utf-8")

    # no relative floor: both blocks called
    bed_a = tmp_path / "a.bed"
    _run(["--input", str(maf), "--output", str(bed_a), "-t", "3", "--ref-id", "hg38",
          "--target-ids", "all", "--min-call-rate", "0", "--min-depth", "2"], SRC)
    assert bed_a.read_text(encoding="utf-8") == "chr1\t0\t2\nchr1\t20\t22\n"

    # call-rate 0.5 with auto N_exp=10 -> floor 5 -> block2 (depth 3) excluded
    bed_b = tmp_path / "b.bed"
    _run(["--input", str(maf), "--output", str(bed_b), "-t", "3", "--ref-id", "hg38",
          "--target-ids", "all", "--min-depth", "2", "--min-call-rate", "0.5"], SRC)
    assert bed_b.read_text(encoding="utf-8") == "chr1\t0\t2\n"
