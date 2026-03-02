# ReAlignPro Release Report (GitHub + Bioconda)

**Project:** ReAlignPro  
**Scope:** Public distribution of a Python CLI toolkit with external bioinformatics dependencies  
**Version covered:** v0.1.1  
**Prepared by:** Chul Lee (chul.bioinfo@gmail.com)  

---

## Executive summary

This report documents the end-to-end process used to package, test, and distribute **ReAlignPro** via:

1. **GitHub** (source-of-truth repository + tagged releases)
2. **Bioconda** (conda distribution for bioinformatics users, including external binary dependencies)
3. **PyPI (optional/parallel)** for `pip install realignpro`

The workflow prioritizes **reproducibility**, **robust CLI behavior**, and **user-facing reliability**, validated through a fast-running end-to-end demo (`examples/test/`).

Key technical milestones included in v0.1.1:
- **Automatic MUSCLE v3/v5 CLI adaptation** (improves portability across environments).
- **Coordinate consistency** between TSS-anchored FASTA headers and gapless MAF coordinates (improves interpretability and downstream validation).

---

## 1) Objectives and distribution targets

### 1.1 Goals
- Provide a **Python package** exposing a stable command-line interface:
  - `realignpro fa2maf`
  - `realignpro maf2bed`
  - `realignpro tsv2fig`
- Provide a **reproducible demo** under `examples/test/` (inputs + runner script + documentation).
- Publish to **GitHub** with clear tagging and release notes for manuscript review.
- Publish to **Bioconda** to maximize adoption among bioinformatics users and bundle external runtime dependencies.
- Maintain manuscript-friendly metadata:
  - `CHANGELOG.md`
  - `CITATION.cff`
  - Nature-style code availability statement

### 1.2 Definition of “successful distribution”
A release is considered successful when all of the following are true:
- `pip install realignpro` works (PyPI optional but recommended).
- `mamba install -c conda-forge -c bioconda realignpro` works (Bioconda).
- The demo workflow runs end-to-end and generates expected outputs (MAF, BED, TSV, PDF figures).

---

## 2) Repository structure and packaging layout

### 2.1 Source layout (Python “src layout”)
ReAlignPro uses the recommended packaging structure:

- `src/realignpro/`
  - `cli.py` (top-level CLI router)
  - `__main__.py` (enables `python -m realignpro`)
  - `__init__.py` (package metadata, `__version__`)
  - `fa2maf.py`
  - `maf2bed.py`
  - `tsv2fig.py`

This prevents accidental imports from the project root and improves packaging reliability.

### 2.2 Demo layout
A minimal, fast-running demo is organized as:

- `examples/test/`
  - `files/`
    - query FASTA(s)
    - target FASTA(s)
  - `run_test.sh` (or `run_demo.sh`)
  - `README.md` describing purpose, expected runtime/hardware, and outputs
  - outputs written to `examples/test/work/`

This structure ensures validation without private datasets and keeps reviewer reproduction straightforward.

### 2.3 Documentation layout
- Root `README.md`: installation, dependencies, quick start, and demo entry point.
- Optional: `docs/workflows/*.md` for command-specific workflows and “paper-ready” documentation.

---

## 3) CLI stabilization and testing strategy

### 3.1 CLI routing and ergonomics
- CLI standardized to `realignpro <subcommand> [options]`.
- Help text made consistent across subcommands.

### 3.2 Lazy imports
To avoid importing heavy dependencies when printing help/version, `cli.py` uses **lazy imports**. This improves:
- `realignpro --help`
- `realignpro <cmd> --help`

and reduces startup overhead and avoidable warnings.

### 3.3 CI strategy: PR tests vs E2E tests
Two complementary test tiers are recommended:

1. **PR-level tests (fast and deterministic)**
   - `pytest` unit/smoke tests (help output, import checks)

2. **E2E demo tests (integration-heavy)**
   - Run `examples/test/run_test.sh`
   - Often scheduled or manual to reduce CI time and avoid flaky dependencies

---

## 4) Core issue resolved in v0.1.1: MUSCLE version compatibility

### 4.1 Problem
MUSCLE v3 and v5 have incompatible CLIs:
- v3 typically uses `-in/-out` and supports flags such as `-diags`, `-maxiters`.
- v5 typically uses `-align/-output` and rejects v3-only flags.

Environments with MUSCLE v5 would fail with:
- `Invalid command line`
- `Unknown option`

### 4.2 Resolution (v0.1.1)
`fa2maf` now:
- Detects MUSCLE major version at runtime.
- Constructs a compatible MUSCLE command for v3 or v5.
- Translates or drops incompatible flags when necessary.

### 4.3 Distribution implication
Because MUSCLE v5 is supported, Bioconda can relax MUSCLE pinning to:
- `muscle >= 3.8.1551`

This prevents forcing all conda users onto MUSCLE v3.

---

## 5) Core issue resolved in v0.1.1: FASTA header vs MAF coordinate mismatch

### 5.1 Problem
TSS-anchored FASTA headers (e.g., `tss_anchored_*.fasta`) could contain coordinates differing from the corresponding gapless MAF (`tss_anchored_*.maf`). This complicates downstream validation.

### 5.2 Resolution (v0.1.1)
During gapless MAF construction, the `(start0, size, strand)` used for MAF `s` lines is used to rewrite FASTA headers consistently:
- MAF start is **0-based** (`start0`)
- FASTA headers are **1-based inclusive**:
  - `start1 = start0 + 1`
  - `end1 = start0 + size`

### 5.3 Verification (example)
If MAF contains:
- `s hg38.NC_000016.10 121 200 + 1100 ...`

FASTA header becomes:
- `>hg38|NC_000016.10|122|321|+`

---

## 6) Repository hygiene and build artifact control

### 6.1 `.gitignore`
A `.gitignore` was added/expanded to exclude:
- `__pycache__/`, `*.pyc`
- `dist/`, `build/`, `*.egg-info/`
- demo work directories: `work/`, `output/`, `intermediates/`

### 6.2 Removing previously tracked caches
Some caches were historically tracked (e.g., `tests/__pycache__`).
They were removed via:
- `git rm -r --cached <path>`

This stabilizes `git status`, `git pull --rebase`, and release tagging.

---

## 7) GitHub release procedure (v0.1.1)

### 7.1 Recommended commit sequence
1. Feature commit(s): MUSCLE v3/v5 support, FASTA/MAF coordinate sync
2. Non-functional cleanup: comment/header cleanup
3. Version bump:
   - `pyproject.toml`: `version = "0.1.1"`
   - `src/realignpro/__init__.py`: `__version__ = "0.1.1"`
4. `CHANGELOG.md` update:
   - Keep `0.1.0` section intact
   - Add a `0.1.1` section describing changes
5. Push to `main` and ensure CI is green
6. Tag + push tag:
   - `git tag -a v0.1.1 -m "ReAlignPro 0.1.1"`
   - `git push origin v0.1.1`
7. Create a GitHub Release:
   - Release notes derived from `CHANGELOG.md`

---

## 8) PyPI distribution (optional/parallel)

### 8.1 Build and upload
- Build: `python -m build`
- Validate: `twine check dist/*`
- Upload: `twine upload dist/*`

### 8.2 Post-upload validation
- Install in a fresh environment:
  - `pip install realignpro`
- Verify:
  - `realignpro --version`
  - help for each subcommand

---

## 9) Bioconda distribution procedure (v0.1.1)

### 9.1 Fork and clone
- Fork `bioconda/bioconda-recipes` to the maintainer account.
- Clone the fork locally (e.g., `~/bioconda_work/bioconda-recipes`).

### 9.2 Branch workflow
- Create a branch for the recipe update (e.g., `realignpro-0.1.1`).

### 9.3 Recipe authoring (`recipes/realignpro/meta.yaml`)
Key choices:
- Source tarball from GitHub tag:
  - `https://github.com/chulbioinfo/ReAlignPro/archive/refs/tags/v{{ version }}.tar.gz`
- Compute `sha256` from the release tarball.
- Runtime dependencies include:
  - `lastz`, `samtools`, `muscle >= 3.8.1551`
  - `matplotlib-base` (preferred over `matplotlib` for Bioconda linting)

### 9.4 sha256 verification pitfalls
A common failure mode is hashing a small HTML error page (e.g., 404). Always confirm:
- downloaded file size is non-trivial (MB-scale, not bytes)
- sha256 is computed on the real tarball
- `curl -L https://github.com/chulbioinfo/ReAlignPro/archive/refs/tags/v0.1.1.tar.gz -o /tmp/realignpro-0.1.1.tar.gz`
- `sha256sum /tmp/realignpro-0.1.1.tar.gz | awk "{print $1}"`

### 9.5 Tests in the recipe
Use robust, fast tests:
- `realignpro --help`
- `realignpro <subcmd> --help`
- `python -c "import realignpro; import realignpro.cli"`

### 9.6 PR labeling and review
- Request bot labeling:
  - `@BiocondaBot please add label`
- Wait for:
  - CI checks to pass
  - at least one approving review

### 9.7 CI troubleshooting
Some CI failures are transient network issues (e.g., `IncompleteRead`) and can be resolved by rerunning the job.

---

## 10) Post-release verification checklist

### 10.1 GitHub
- Tag exists and tarball downloads correctly (large file size).
- Release notes match `CHANGELOG.md`.

### 10.2 Bioconda / conda install
In a clean environment:
- `mamba install -c conda-forge -c bioconda realignpro`
- Confirm external binaries:
  - `which lastz`, `which samtools`, `which muscle`
- Verify CLI:
  - `realignpro --version`
  - `realignpro fa2maf --help` etc.

### 10.3 End-to-end demo (installed package + repo demo)
- Clone GitHub repo and run:
  - `bash examples/test/run_test.sh`
- Confirm outputs exist:
  - `merged_orthologs_aligned.maf`
  - `merged_orthologs_aligned.bed`
  - `conservation_matrix.tsv`
  - `conservation_matrix.upstream.pdf`
  - `conservation_matrix.downstream.pdf`

---

## 11) Final deliverables

### 11.1 GitHub
- Public repository
- Tagged release `v0.1.1`
- Demo dataset + scripts under `examples/test/`
- Documentation (README + optional workflow docs)
- Release notes + citation metadata

### 11.2 PyPI (optional/parallel)
- `pip install realignpro`

### 11.3 Bioconda
- Merged recipe in `bioconda/bioconda-recipes`
- Package visible on Anaconda.org under the `bioconda` channel

---

## Appendix: minimal command snippets

### A) GitHub tagging

```bash
git tag -a v0.1.1 -m "ReAlignPro 0.1.1"
git push origin v0.1.1
