# ReAlignPro end to end demo (fa2maf, maf2bed, tsv2fig)

This folder contains a small, self-contained dataset and a single runner script, `run_test.sh`, that executes the full ReAlignPro workflow from FASTA inputs to publication-ready PDFs.

Author: Chul Lee (chul.bioinfo@gmail.com)

## What this demo demonstrates

This demo is designed to be a realistic, minimal example of the core ReAlignPro pipeline:

1. **fa2maf**
   - Uses a short query fragment to find orthologous regions in multiple target genomes via local alignment (LASTZ).
   - Extracts orthologous sequences from each target.
   - Builds a multiple sequence alignment (MUSCLE) and writes:
     - a merged MAF alignment for downstream analysis
     - a conservation matrix TSV for visualization

2. **maf2bed**
   - Scans the merged MAF and generates BED3 intervals for sites that satisfy user-defined constraints, such as the presence of all target species in each processed MAF block.
   - Produces a BED file that can be used as a genome browser track or as input for downstream filtering.

3. **tsv2fig**
   - Converts the conservation matrix TSV into PDF figures for upstream and downstream windows, optionally highlighting a motif (for example, "CCG").

The dataset and defaults are tuned to be quick to run while still exercising the same code paths used for real analyses.

## Expected runtime and hardware

This demo is intended to complete quickly on common research hardware.

- **Typical runtime**: tens of seconds to a few minutes, depending on CPU, disk I/O, and external tool availability.
- **CPU**: 4 cores recommended (the default `THREADS=4` is used in `run_test.sh`).
- **Memory**: 2 to 4 GB RAM is typically sufficient for this demo dataset.
- **OS**: Linux and macOS are recommended for the standard bioinformatics toolchain.

If you increase the window size, number of target genomes, or the amount of sequence, runtime will increase.

## Reproducibility and version pinning

For reproducible results, record and pin both:
- the ReAlignPro version (or git commit) and
- the versions of external binaries (LASTZ, MUSCLE, samtools), plus Python dependencies.

Recommended practices:

1. **Use tagged releases**
   - Prefer running a tagged release of ReAlignPro for analyses that will be reported in a manuscript.
   - Record the tag name and commit hash in your notes.

2. **Pin your conda environment**
   - Create an environment file and pin versions explicitly.
   - Example (edit versions as needed):

   ```yaml
   name: realignpro-demo
   channels:
     - conda-forge
     - bioconda
   dependencies:
     - python=3.11
     - realignpro=0.1.0
     - lastz
     - muscle
     - samtools
     - matplotlib
   ```

3. **Capture an exact environment snapshot**
   - For an exact reproduction on the same platform:

   ```bash
   conda list --explicit > conda-spec.txt
   ```

   - For a shorter, more portable file (less strict, better across platforms):

   ```bash
   conda env export --from-history > environment.yml
   ```

4. **Record tool versions**
   - Save these outputs alongside your results:

   ```bash
   realignpro --version
   lastz --help | head
   muscle -version || true
   samtools --version | head
   ```

If you later publish a Zenodo DOI or another persistent archive for a specific release, cite that DOI in the manuscript and point users to the corresponding tag.

## Folder layout

```
examples/
  test/
    files/
      query_hg38genome.fa
      query_RBFOX1fragment.fa
      target/
        Gorilla_gorilla.fa
        Homo_sapiens.fa
        Pan_paniscus.fa
        Pan_troglodytes.fa
    run_test.sh
    README.md
    work/                 # created and used by the demo
```

Inputs:
- `files/query_hg38genome.fa`, reference genome FASTA used for coordinate system and indexing.
- `files/query_RBFOX1fragment.fa`, query fragment FASTA used as a seed for local alignment.
- `files/target/*.fa`, target genomes, one FASTA per species.

Outputs (generated):
- `work/output/`, the main output directory created by the workflow.

## Requirements

ReAlignPro relies on external command-line tools that must be available in `PATH`:

- `realignpro`
- `lastz`
- `muscle`
- `samtools`

### Recommended installation (conda)

```bash
mamba install -c conda-forge -c bioconda realignpro lastz muscle samtools
```

Verify installation:

```bash
realignpro --version
which lastz
which muscle
which samtools
```

## Quick start

From `examples/test/`:

```bash
chmod +x run_test.sh
./run_test.sh
```

The script runs the full pipeline:

1. `realignpro fa2maf`
2. `realignpro maf2bed`
3. `realignpro tsv2fig`

All outputs are written under:

- `work/output/`

## Pipeline stages and key options

### Stage 1, fa2maf

Purpose:
- Identify orthologous regions, build alignments, and write MAF plus conservation TSV.

Inputs used in this demo:
- `--genome_dir ./files/target/`
- `--query_fa ./files/query_RBFOX1fragment.fa`
- `--ref_genome ./files/query_hg38genome.fa`
- `--work_dir ./work`

Key options (defaults shown as used by `run_test.sh`):
- `--threads 4`
- `--lastz lastz`
- `--muscle muscle`
- `--samtools samtools`
- `--lastz_params "T=1 --strand=both --ambiguous=iupac"`
- `--tss_genomic "NC_000016.10:222"`
- `--window 100`
- `--ConservationThreshold 1`

Expected outputs under `work/output/`:
- `merged_orthologs_aligned.maf` (or `.maf.gz`)
- `conservation_matrix.tsv`

### Stage 2, maf2bed

Purpose:
- Convert the merged MAF into BED3 intervals based on target constraints.

Input:
- `--input work/output/merged_orthologs_aligned.maf` (or `.maf.gz`)

Required options:
- `--ref-id hg38`
- `--target-ids hg38,Homo_sapiens`

Important:
- `--target-ids` must be **comma-separated with no spaces**.
  - Valid: `hg38,Homo_sapiens`
  - Invalid: `hg38 Homo_sapiens`

Optional options:
- `--outgroup-ids <comma-separated>`
- `--threads 4` (must be >= 3)
- `--work-qsize 200`
- `--out-qsize 200`
- `--output <path>`, if omitted the tool derives a `.bed` path from the input MAF

Expected output under `work/output/`:
- `merged_orthologs_aligned.bed`

### Stage 3, tsv2fig

Purpose:
- Generate PDF figures from the conservation matrix TSV.

Input:
- `work/output/conservation_matrix.tsv`

Options used in this demo:
- `--type variant`
- `--highlight CCG`
- `--outdir work/output`
- `--prefix conservation_matrix`

Expected outputs under `work/output/`:
- `conservation_matrix.upstream.pdf`
- `conservation_matrix.downstream.pdf`

## Customization via environment variables

You can override parameters without editing `run_test.sh`:

```bash
# increase threads
THREADS=8 ./run_test.sh

# override maf2bed targets
REF_ID=hg38 TARGET_IDS=hg38,Homo_sapiens ./run_test.sh

# add an outgroup
OUTGROUP_IDS=Pan_troglodytes ./run_test.sh

# change window and highlight motif
WINDOW=200 HIGHLIGHT=CCG ./run_test.sh
```

Variables supported by `run_test.sh`:
- `THREADS` (default 4)
- `LASTZ_BIN` (default `lastz`)
- `MUSCLE_BIN` (default `muscle`)
- `SAMTOOLS_BIN` (default `samtools`)
- `LASTZ_PARAMS` (default `T=1 --strand=both --ambiguous=iupac`)
- `TSS_GENOMIC` (default `NC_000016.10:222`)
- `WINDOW` (default 100)
- `CONS_THRESHOLD` (default 1)
- `REF_ID` (default `hg38`)
- `TARGET_IDS` (default `hg38,Homo_sapiens`)
- `OUTGROUP_IDS` (default empty)
- `HIGHLIGHT` (default `CCG`)
- `VARIANT_TYPE` (default `variant`)

## Troubleshooting

### realignpro not found
Activate your environment, then confirm:

```bash
which realignpro
realignpro --help
```

### lastz, muscle, or samtools missing
Confirm binaries are in `PATH`:

```bash
which lastz
which muscle
which samtools
```

Install via conda if needed:

```bash
mamba install -c conda-forge -c bioconda lastz muscle samtools
```

### Missing outputs after Stage 1
If the MAF or TSV is not produced:
- confirm `files/target/` contains readable FASTA files
- confirm `files/query_hg38genome.fa` is valid and indexable (samtools may create a `.fai`)
- inspect logs under `work/`

## Repository hygiene

It is recommended to ignore generated outputs in git, such as:
- `work/`
- `*.fai`, `*.maf*`, `*.bed`, `*.tsv`, `*.pdf`, `*.log`
