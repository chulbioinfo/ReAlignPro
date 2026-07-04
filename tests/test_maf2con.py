from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

SRC = Path(__file__).resolve().parents[1] / "src"


def _run(args) -> None:
    env = dict(os.environ)
    env["PYTHONPATH"] = str(SRC) + os.pathsep + env.get("PYTHONPATH", "")
    subprocess.run(
        [sys.executable, "-m", "realignpro", "maf2con", *args],
        check=True, timeout=60, env=env,
    )


def _alt_base(base: str) -> str:
    return {"A": "C", "C": "G", "G": "T", "T": "A"}[base]


def _write_1000_haplotype_maf(path: Path) -> None:
    """
    100 bp, 1000-haplotype single-block MAF (hg38 + hs1..hs999):
      - position 20: 990/1000 major -> ratio 0.99, NOT constrained (strict >).
      - position 50: 991/1000 major -> constrained.
    """
    ref_seq = "ACGT" * 25
    ids = ["hg38"] + [f"hs{i}" for i in range(1, 1000)]
    lines = ["##maf version=1\n\n", "a score=0\n"]
    for idx, sid in enumerate(ids):
        seq = list(ref_seq)
        if 1 <= idx <= 10:
            seq[20] = _alt_base(seq[20])
        if 1 <= idx <= 9:
            seq[50] = _alt_base(seq[50])
        lines.append(f"s {sid}.chr1 0 100 + 100 {''.join(seq)}\n")
    lines.append("\n")
    path.write_text("".join(lines), encoding="utf-8")


def test_maf2con_default_all_targets_1000_haplotypes(tmp_path: Path) -> None:
    # Default engine is coverage-aware (call-rate 0.99, similarity 0.99). With all 1000
    # present, N_exp=1000 -> floor 990; every column clears it, so the output is governed
    # by conservation: pos20 (990/1000=0.99) excluded, pos50 (991/1000) kept.
    maf_path = tmp_path / "haplotypes_1000_100bp.maf"
    bed_path = tmp_path / "haplotypes_1000_100bp.con.bed"
    _write_1000_haplotype_maf(maf_path)

    _run(["--input", str(maf_path), "--output", str(bed_path), "--threads", "3",
          "--ref-id", "hg38", "--target-ids", "all"])

    assert bed_path.read_text(encoding="utf-8") == "chr1\t0\t20\nchr1\t21\t100\n"
