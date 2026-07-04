# ReAlignPro v0.2.0 — coverage-aware `maf2con`

**Release date:** 2026-07-03
**Scope:** `realignpro maf2con` only. `fa2maf`, `maf2bed`, and `tsv2fig` are unchanged.

## Highlights

`maf2con` now calls constrained regions with a **coverage-aware** engine, so it works
correctly when the set of aligned assemblies varies from block to block — most importantly
on **sex chromosomes** (chrX / chrY), where alignment depth changes along the chromosome.
The previous fixed-denominator / whole-block-gate engine (v0.1.2) has been removed and the
coverage-aware engine is now the default and only engine.

## What changed

- **Per-column, coverage-aware denominator.** For each reference base the major-allele
  fraction is `n_major / N_cov`, where `N_cov` is the number of *aligned* A/C/G/T at that
  column (gap / `N` / missing are excluded from both numerator and denominator). Previously
  the denominator was the fixed total target count.
- **No whole-block gate.** A block is assessed on whichever target assemblies are actually
  present, instead of being skipped unless *all* targets are present. This is what made the
  old engine emit **zero** intervals on chrX/chrY when `--target-ids all` spanned assemblies
  that only align to autosomes.
- **Per-column coverage floor.** A column is callable only when
  `N_cov >= --min-depth` and, when the relative floor is in effect,
  `N_cov >= ceil(--min-call-rate * N_exp)`, where `N_exp` is the reference chromosome's
  expected aligned depth.
- **Automatic per-chromosome `N_exp`.** When the relative floor is used and
  `--expected-depth` is not given, `N_exp` is derived per reference chromosome by a cheap
  pre-scan, so autosomes, chrX, and chrY are handled correctly in a **single run** with no
  manual numbers.

## New / changed options

| Option | Meaning |
|---|---|
| `--min-major-similarity` | Strict lower bound (`>`) for the major-allele fraction among aligned bases. Default **0.99**. |
| `--min-call-rate` | Relative coverage floor: require `N_cov >= ceil(rate * N_exp)`. Default **0.99**; set `0` to disable and use only `--min-depth`. |
| `--min-depth` | Absolute coverage floor (minimum aligned A/C/G/T). Default `2`. |
| `--expected-depth` | Pin `N_exp` to a single value and **skip the pre-scan** (saves one read pass). |
| `--fixed-only` | Strictest mode: report only 100 %-fixed (monomorphic) columns — a single distinct aligned allele, zero mismatch. Ignores `--min-major-similarity`. |
| `--emit-depth` | Append `N_cov` and major-allele-frequency columns to each BED interval (for downstream QC / filtering). |
| ~~`--coverage-aware`~~, ~~`--auto-depth`~~ | **Removed.** Both are now the default / implicit behavior. |

## Examples

```bash
# Default: coverage-aware, per-chromosome N_exp, call-rate 0.99, similarity 0.99.
realignpro maf2con -i merged.maf.gz --ref-id hg38 --target-ids all

# Strictest: only 100%-fixed columns, with depth columns for QC.
realignpro maf2con -i merged.maf.gz --ref-id hg38 --target-ids all --fixed-only --emit-depth

# Pin the expected depth to skip the pre-scan (single-chromosome file).
realignpro maf2con -i chrX.maf.gz --ref-id hg38 --target-ids all --expected-depth 363

# Reproduce a pure per-column call with no relative floor.
realignpro maf2con -i merged.maf.gz --ref-id hg38 --target-ids all --min-call-rate 0
```

## Breaking change / migration

- The output of `maf2con` and its CLI differ from v0.1.2. Since `maf2con` was introduced in
  0.1.2, this only affects that command. To reproduce v0.1.2 behavior, pin
  `realignpro==0.1.2`.
- Scripts passing `--coverage-aware` or `--auto-depth` must drop those flags (they are now
  the default). No other command is affected.

## Notes on sex chromosomes

Homologous cross-alignment between chrX and chrY (PAR / XTR / gametolog regions) is
**intentionally not separated**. Pooling homologous lineages only dilutes the major-allele
fraction, so the dominant effect is a conservative false negative (a missed constraint),
never a false positive at the core. Use `--emit-depth` to keep `N_cov` / major-allele
frequency for downstream region-aware filtering.

## Internal changes

- New module `src/realignpro/maf2con_cov.py` (coverage-aware engine: `call_major_base_cov`,
  `matrix2con_cov`, `scan_expected_depth`, and the reader/worker/writer orchestration).
- `src/realignpro/maf2con.py` is now a thin CLI + config layer that dispatches to the
  coverage-aware engine; the legacy `matrix2con` / `call_major_base` / multiprocessing code
  was removed.
- Tests: `tests/test_maf2con_cov.py` (unit + end-to-end); the legacy-only test in
  `tests/test_maf2con.py` was updated to the new default behavior.

See `CHANGELOG.md` for the concise entry and `docs/maf2con_plan.md` for the (historical)
v0.1.2 design.
