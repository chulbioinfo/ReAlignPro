# ReAlignPro v0.2.1 — `maf2bed` skips missing/ambiguous target bases

**Release date:** 2026-08-06
**Scope:** `realignpro maf2bed` only. `fa2maf`, `maf2con`, and `tsv2fig` are unchanged.

## Highlights

`maf2bed` no longer reports positions where a target assembly has no unambiguous base.
A reference column is assessed only when **every** target carries an `A`, `C`, `G`, or `T`;
a gap, `N`, or any other IUPAC ambiguity code in any target skips the column.

This closes a false-positive class in which an **assembly-gap run shared by all targets was
emitted as a target-specific interval**. `maf2con` has always restricted calls to `A/C/G/T`
(`DNA_BASES` in `maf2con_cov.py`); `maf2bed` now uses the same rule.

## The bug

`matrix2var` upper-cased each column and used the characters as alleles directly, filtering
only the reference gap (`ref_nt == "-"`). A column therefore passed the
"targets share exactly one allele, and that allele is absent from the others" test whenever
the targets were uniformly `N` (or uniformly `-`, when the reference is not itself a target):

```
ref/target 1   ...ACGT NN ACGT...
target 2       ...ACGT NN ACGT...     tset == {"N"}, "N" not in oset  ->  called
outgroup       ...TTTT AA GGCC...
```

Because adjacent hits are merged, a single N-run also **fused otherwise separate intervals**
into one longer interval. The inflation scales with assembly-gap content, so it was largest
in exactly the regions where target data is weakest.

## The fix

`src/realignpro/maf2bed.py`:

- Added `DNA_BASES = {"A", "C", "G", "T"}` (same definition and name as `maf2con_cov.py`).
- In `matrix2var`, after the per-column allele lists are built:

```python
if any(nt not in DNA_BASES for nt in target_nts):
    flush_current()
    continue
```

`flush_current()` closes the open interval before skipping, so a masked column breaks the
merge run instead of bridging the hits on either side. Correct on both `+` and `-` strand
(the `-` strand walks high→low and merges on `base_end == cur_start`).

Outgroup / non-target species are unchanged: they are still evaluated as before, so
`N` in a non-target is still read as "this allele is absent here". Sequence is still
upper-cased before the check, so soft-masked lowercase bases are called normally.

## Output impact

Strictly fewer / shorter BED intervals than 0.2.0 on any alignment containing assembly gaps;
identical output on gap-free alignments. Every removed position is one where no unambiguous
target base existed, so none of them were supported calls.

Observed on a synthetic block (`ref = ACGTNNACGT`, both targets identical, outgroup
`TTTTAAGGCC`):

| | 0.2.0 | 0.2.1 |
|---|---|---|
| BED | `chr1 100 103`, `chr1 104 110` | `chr1 100 103`, `chr1 106 110` |

Per-column behavior:

| Column | 0.2.0 | 0.2.1 |
|---|---|---|
| All targets `N` | called | **skipped** |
| All targets `-` (reference not a target) | called | **skipped** |
| One target `N`, mid-interval | merged across | interval split at the `N` |
| One target `R` (IUPAC) | already skipped (2 alleles) | skipped |
| Target `A`, all others `N` | called | called (unchanged) |
| Soft-masked `aaa` vs `AAA` | called | called (unchanged) |

## Upgrade

```bash
pip install --upgrade realignpro
```

No CLI or option changes; no script changes required. If a downstream BED must be
byte-identical to a previous run, regenerate it — do not mix 0.2.0 and 0.2.1 outputs.
Pin `realignpro==0.2.0` to reproduce the old behavior.

## Internal changes

- `src/realignpro/maf2bed.py`: `DNA_BASES` constant, the target-column guard in `matrix2var`,
  and docstring updates to the "Variant definition" list.
- `tests/test_maf2bed.py` (new): unit tests for the `A/C/G/T`-only rule, gap-only targets,
  IUPAC codes, merge-breaking on both strands, and soft-mask handling.

See `CHANGELOG.md` for the concise entry.
