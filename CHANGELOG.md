# Changelog

All notable changes to **ReAlignPro** will be documented in this file.

This project follows a lightweight semantic versioning scheme:
- **MAJOR**: incompatible CLI or output changes
- **MINOR**: new features that remain backward compatible
- **PATCH**: bug fixes and minor improvements

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


[0.1.1]: https://github.com/chulbioinfo/ReAlignPro/releases/tag/v0.1.1
[0.1.0]: https://github.com/chulbioinfo/ReAlignPro/releases/tag/v0.1.0
