# Changelog

All notable changes to **ReAlignPro** will be documented in this file.

This project follows a lightweight semantic versioning scheme:
- **MAJOR**: incompatible CLI or output changes
- **MINOR**: new features that remain backward compatible
- **PATCH**: bug fixes and minor improvements

## [0.2.2] - 2026-08-20

### Added
- `realignpro maf2bed` now decompresses **BGZF** input through `bgzip -@` instead of
  Python's `gzip` module. BGZF is what `bgzip`, `cactus-hal2maf` and `taffy view -c`
  write, and although it is a valid gzip stream (so `gzip.open` reads it correctly),
  `gzip.open` decodes it as one serial stream while `bgzip` decodes its blocks in
  parallel. maf2bed has a single reader feeding every worker, so the reader's speed
  is the whole tool's ceiling.
  - Measured on `vgp-577way.mammals.maf.gz` (186 GB, BGZF): reader throughput
    **243 MB/s -> 843 MB/s (3.5x)**.
  - New `--reader-threads` (default 4) sets `bgzip -@`. It is *additional* to
    `--threads`, which still counts 1 reader + 1 writer + workers.
  - Automatic and safe: the BGZF magic (FEXTRA plus the "BC" subfield) is checked,
    and the code falls back to `gzip.open` for plain gzip, for non-BGZF input, or
    when `bgzip` is not on PATH. `fa2maf`, `maf2con` and `tsv2fig` are unchanged.

> **Output impact: none.** Verified byte-identical to 0.2.1 across BGZF, plain gzip
> and uncompressed forms of the same MAF, in both directions.

### Fixed
- A truncated or corrupt BGZF file now raises instead of being read as a short file:
  a non-zero `bgzip` exit becomes a `RuntimeError` (SIGPIPE from stopping the read
  early is excluded). Previously a decompression failure could look like a
  legitimate end of file and silently yield a partial BED.
  - The check waits up to `BGZIP_EXIT_GRACE_SEC` (60 s) for `bgzip` to exit on its
    own before terminating it, and never reports a failure for a process it killed
    itself. Reaching EOF on the pipe does not mean the child has been reaped: on a
    24 GB input that race made a completed read report `exit code -15` (SIGTERM
    sent by this very code path). It is also skipped while another exception is
    unwinding, so it cannot mask the real cause.
- The reader now pushes its worker stop sentinels from a `finally`. Previously a
  reader that died part-way left every worker blocked on an empty queue, so the
  run hung silently instead of failing -- one such hang burned 6 h of wall clock
  before it was noticed. Present in 0.2.1 as well (reproducible there with a
  truncated plain-gzip input); fixed here.

## [0.2.1] - 2026-08-06

### Fixed
- `realignpro maf2bed` no longer calls columns where a target base is missing or ambiguous.
  A column is assessed only when **every** target carries an unambiguous `A/C/G/T`; a gap, `N`,
  or any other IUPAC ambiguity code in any target skips the column. Previously `N` and `-` were
  treated as ordinary alleles, so an assembly-gap run shared by all targets (`tset == {"N"}`)
  satisfied the "targets share one allele, others lack it" rule and was emitted as a BED interval.
  This brings `maf2bed` in line with `maf2con`, which has always restricted calls to `A/C/G/T`.
  - Skipped columns also break interval merging, so a masked column no longer bridges two
    otherwise separate hits into one interval (verified on both `+` and `-` strand).
  - Soft-masked lowercase bases are unaffected — sequence is still upper-cased before the check.

> **Output impact:** `maf2bed` emits strictly fewer / shorter intervals than 0.2.0 on alignments
> containing assembly gaps. The positions removed are those with no unambiguous target base, which
> could not have been supported calls. `fa2maf`, `maf2con`, and `tsv2fig` are unchanged.

### Added
- `tests/test_maf2bed.py`: regression tests for the `A/C/G/T`-only target rule, interval-merge
  breaking on both strands, and soft-mask handling.

## [0.2.0] - 2026-07-03

### Changed (BREAKING)
> `maf2con` was introduced in 0.1.2, so this reworks a brand-new command; the other
> commands (`fa2maf`, `maf2bed`, `tsv2fig`) are unchanged. The behavior/CLI change is
> incompatible with 0.1.2's `maf2con` — pin `realignpro==0.1.2` to reproduce it.

- `realignpro maf2con` is now **coverage-aware by default and is the only engine**; the legacy fixed-denominator / whole-block-gate engine has been removed.
  - Per reference column the denominator is the aligned depth `N_cov` (gap/N excluded from numerator and denominator), the whole-block gate is gone, and a per-column coverage floor is applied so blocks with missing assemblies are still assessed.
  - New defaults: `--min-call-rate 0.99` (require `N_cov >= ceil(rate * N_exp)`, where `N_exp` is the reference chromosome's expected depth) and `--min-major-similarity 0.99`.
  - `N_exp` is auto-derived per reference chromosome by a cheap pre-scan when the relative floor is in effect and `--expected-depth` is not given, so autosomes / chrX / chrY are handled in one run.
  - Removed flags `--coverage-aware` and `--auto-depth` (both are now the default/implicit behavior). To reproduce the previous behavior, use the 0.1.2 release.

### Added
- `--min-depth` (absolute coverage floor), `--expected-depth` (pin `N_exp` and skip the pre-scan), `--min-call-rate` (relative floor; `0` disables it), `--fixed-only` (report only 100% fixed / monomorphic columns, ignoring `--min-major-similarity`), `--emit-depth` (append `N_cov` and major-allele frequency columns to each BED interval).
- New module `realignpro.maf2con_cov`; tests `tests/test_maf2con_cov.py`.

### Notes
- Homologous X/Y cross-alignment (PAR / XTR / gametolog) is intentionally not separated; pooling only dilutes the major-allele fraction, a conservative false negative. See `ConstraintVariantAnalysisPlan_v0.3.md`.

## [0.1.2] - 2026-05-08

### Added
- `realignpro maf2con`: conversion of MAF/MAF.GZ to BED3 constrained intervals using strict target-group major-allele similarity, with `--target-ids all` support.

## [0.1.1] - 2026-03-01

### Changed
- MUSCLE v5 is now supported without pinning by auto-detecting MUSCLE v3 versus v5 and adapting the CLI invocation accordingly.
- Synchronize TSS-anchored FASTA header coordinates with gapless MAF coordinates to ensure consistent genomic intervals across outputs.

### Fixed
- Removed tracked Python bytecode cache files (`__pycache__`) from the repository.


## [0.1.0] - 2026-02-27

### Added
- Initial public release of the `realignpro` command-line interface with three subcommands:
  - `realignpro fa2maf`: local alignment-based ortholog retrieval and multiple sequence alignment, exporting MAF and a conservation-matrix TSV.
  - `realignpro maf2bed`: conversion of MAF/MAF.GZ to BED3 intervals using required `--ref-id` and comma-separated `--target-ids`.
  - `realignpro tsv2fig`: generation of publication-ready upstream/downstream PDF figures from the conservation-matrix TSV, with optional motif highlighting.
- End-to-end demo dataset and runner script under `examples/test/` (`run_test.sh`, `README.md`).

### Changed
- Deprecated the prior CTL-file control mode for the MAF-to-BED step in favour of explicit command-line options.
- CLI imports are performed lazily to avoid importing heavy dependencies when printing help/version.

### Fixed
- Continuous integration end-to-end demo reliability:
  - Pinned MUSCLE to version 3.8.1551 for the CI workflow due to CLI incompatibilities with MUSCLE v5.
  - Added missing standard-library import required by the MAF-to-BED step.


[0.1.2]: https://github.com/chulbioinfo/ReAlignPro/releases/tag/v0.1.2
[0.1.1]: https://github.com/chulbioinfo/ReAlignPro/releases/tag/v0.1.1
[0.1.0]: https://github.com/chulbioinfo/ReAlignPro/releases/tag/v0.1.0
