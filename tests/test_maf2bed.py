"""Unit tests for the maf2bed calling logic.

Focus: a column must not be called when any target base is missing or ambiguous
(gap, N, or another IUPAC code). These run directly against ``matrix2var`` so no
MAF file or multiprocessing setup is needed.
"""

from __future__ import annotations

from realignpro.maf2bed import matrix2var


def _block(seqs: dict[str, str], strand: str = "+") -> dict[str, dict[str, object]]:
    return {
        sid: {
            "seq": seq,
            "chr": "chr1",
            "start": 100,
            "size": len(seq.replace("-", "")),
            "strand": strand,
            "srcSize": 1000,
        }
        for sid, seq in seqs.items()
    }


def _call(seqs: dict[str, str], targets: list[str], strand: str = "+") -> list[str]:
    block = _block(seqs, strand)
    return [line.strip() for line in matrix2var(block, "hg38", targets, [])]


def test_shared_allele_is_called() -> None:
    assert _call({"hg38": "AAA", "t2": "AAA", "o1": "CCC"}, ["hg38", "t2"]) == ["chr1\t100\t103"]


def test_n_run_in_all_targets_is_skipped() -> None:
    assert _call({"hg38": "NNN", "t2": "NNN", "o1": "ACG"}, ["hg38", "t2"]) == []


def test_gap_only_targets_are_skipped() -> None:
    # Reference is not a target here, so the reference-gap shortcut does not apply.
    assert _call({"hg38": "ACG", "t1": "---", "t2": "---", "o1": "TTT"}, ["t1", "t2"]) == []


def test_iupac_code_in_one_target_is_skipped() -> None:
    assert _call({"hg38": "ARA", "t2": "AAA", "o1": "CCC"}, ["hg38", "t2"]) == [
        "chr1\t100\t101",
        "chr1\t102\t103",
    ]


def test_skipped_column_breaks_the_merged_interval() -> None:
    # Without the flush, the N column would bridge positions 100 and 102 into one interval.
    assert _call({"hg38": "ANA", "t2": "ANA", "o1": "CCC"}, ["hg38", "t2"]) == [
        "chr1\t100\t101",
        "chr1\t102\t103",
    ]


def test_skipped_column_breaks_interval_on_minus_strand() -> None:
    assert _call({"hg38": "ANA", "t2": "ANA", "o1": "CCC"}, ["hg38", "t2"], strand="-") == [
        "chr1\t899\t900",
        "chr1\t897\t898",
    ]


def test_soft_masked_lowercase_bases_are_still_called() -> None:
    assert _call({"hg38": "aaa", "t2": "AAA", "o1": "CCC"}, ["hg38", "t2"]) == ["chr1\t100\t103"]
