from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

from realignpro.maf2con import call_major_base


def _alt_base(base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def _write_1000_haplotype_maf(path: Path) -> None:
    """
    Create a 100 bp, 1000-haplotype MAF:
      - hg38 + hs1..hs999
      - position 20 has 990/1000 major support and remains constrained with >=0.99
      - position 50 has 989/1000 major support and is not constrained
    """
    ref_seq = "ACGT" * 25
    ids = ["hg38"] + [f"hs{i}" for i in range(1, 1000)]

    lines = ["##maf version=1\n\n", "a score=0\n"]
    for idx, sid in enumerate(ids):
        seq = list(ref_seq)
        if 1 <= idx <= 10:
            seq[20] = _alt_base(seq[20])
        if 1 <= idx <= 11:
            seq[50] = _alt_base(seq[50])
        lines.append(f"s {sid}.chr1 0 100 + 100 {''.join(seq)}\n")
    lines.append("\n")

    path.write_text("".join(lines), encoding="utf-8")


def test_call_major_base_uses_inclusive_threshold() -> None:
    base, ratio = call_major_base(["A"] * 990 + ["C"] * 10, min_major_similarity=0.99)
    assert base == "A"
    assert ratio == pytest.approx(0.99)

    base, ratio = call_major_base(["A"] * 989 + ["C"] * 11, min_major_similarity=0.99)
    assert base is None
    assert ratio == pytest.approx(0.989)


def test_maf2con_all_targets_on_1000_haplotype_100bp_maf(tmp_path: Path) -> None:
    maf_path = tmp_path / "haplotypes_1000_100bp.maf"
    bed_path = tmp_path / "haplotypes_1000_100bp.con.bed"
    _write_1000_haplotype_maf(maf_path)

    cmd = [
        sys.executable,
        "-m",
        "realignpro",
        "maf2con",
        "--input",
        str(maf_path),
        "--output",
        str(bed_path),
        "--threads",
        "3",
        "--ref-id",
        "hg38",
        "--target-ids",
        "all",
    ]
    subprocess.run(cmd, check=True, timeout=30)

    assert bed_path.read_text(encoding="utf-8") == "chr1\t0\t50\nchr1\t51\t100\n"
