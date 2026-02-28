#!/usr/bin/env bash
set -euo pipefail

# ReAlignPro end-to-end demo runner (fa2maf -> maf2bed -> tsv2fig)
# Intended to be executed from this directory:
#   examples/test/
#
# Requirements (recommended via conda):
#   - realignpro (this package)
#   - lastz, muscle, samtools available in PATH

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# ---- Inputs (relative to this demo directory) ----
TARGET_DIR="${HERE}/files/target"
QUERY_FRAG="${HERE}/files/query_RBFOX1fragment.fa"
REF_GENOME="${HERE}/files/query_hg38genome.fa"

# ---- Work/output directory ----
WORK_DIR="${HERE}/work"

# ---- Parameters (edit if needed) ----
THREADS="${THREADS:-4}"
LASTZ_BIN="${LASTZ_BIN:-lastz}"
MUSCLE_BIN="${MUSCLE_BIN:-muscle}"
SAMTOOLS_BIN="${SAMTOOLS_BIN:-samtools}"

LASTZ_PARAMS="${LASTZ_PARAMS:-T=1 --strand=both --ambiguous=iupac}"
TSS_GENOMIC="${TSS_GENOMIC:-NC_000016.10:222}"
WINDOW="${WINDOW:-100}"
CONS_THRESHOLD="${CONS_THRESHOLD:-1}"

# maf2bed requirements
REF_ID="${REF_ID:-hg38}"
TARGET_IDS="${TARGET_IDS:-hg38,Homo_sapiens}"   # comma-separated, no spaces
OUTGROUP_IDS="${OUTGROUP_IDS:-}"               # optional (comma-separated)

# tsv2fig
HIGHLIGHT="${HIGHLIGHT:-CCG}"
VARIANT_TYPE="${VARIANT_TYPE:-variant}"

echo "[INFO] Demo directory : ${HERE}"
echo "[INFO] Work directory : ${WORK_DIR}"
echo "[INFO] Using threads  : ${THREADS}"
echo

# ---- Basic checks ----
if ! command -v realignpro >/dev/null 2>&1; then
  echo "[ERROR] 'realignpro' not found in PATH. Did you install the conda package?" >&2
  exit 2
fi

echo "[INFO] realignpro version:"
realignpro --version || true
echo

# External binaries are typically installed via conda (bioconda/conda-forge).
for bin in "${LASTZ_BIN}" "${MUSCLE_BIN}" "${SAMTOOLS_BIN}"; do
  if ! command -v "${bin}" >/dev/null 2>&1; then
    echo "[WARN] Required binary not found in PATH: ${bin}"
  fi
done
echo

# ---- Step 1: fa2maf ----
echo "[STEP 1/3] realignpro fa2maf"
mkdir -p "${WORK_DIR}"

realignpro fa2maf   --genome_dir "${TARGET_DIR}"   --query_fa "${QUERY_FRAG}"   --ref_genome "${REF_GENOME}"   --work_dir "${WORK_DIR}"   --lastz "${LASTZ_BIN}"   --samtools "${SAMTOOLS_BIN}"   --muscle "${MUSCLE_BIN}"   --threads "${THREADS}"   --lastz_params "${LASTZ_PARAMS}"   --tss_genomic "${TSS_GENOMIC}"   --window "${WINDOW}"   --ConservationThreshold "${CONS_THRESHOLD}"

echo

# Determine expected outputs (fa2maf writes into WORK_DIR/output)
OUT_DIR="${WORK_DIR}/output"
MAF_IN="${OUT_DIR}/merged_orthologs_aligned.maf"
if [[ -f "${MAF_IN}.gz" ]]; then
  MAF_IN="${MAF_IN}.gz"
fi

TSV_IN="${OUT_DIR}/conservation_matrix.tsv"

if [[ ! -f "${MAF_IN}" ]]; then
  echo "[ERROR] Expected MAF not found: ${MAF_IN}" >&2
  echo "       Please check fa2maf logs under: ${WORK_DIR}" >&2
  exit 2
fi

if [[ ! -f "${TSV_IN}" ]]; then
  echo "[ERROR] Expected TSV not found: ${TSV_IN}" >&2
  echo "       Please check fa2maf logs under: ${WORK_DIR}" >&2
  exit 2
fi

# ---- Step 2: maf2bed ----
echo "[STEP 2/3] realignpro maf2bed"
BED_OUT="${OUT_DIR}/merged_orthologs_aligned.bed"

MAF2BED_ARGS=(
  --input "${MAF_IN}"
  --output "${BED_OUT}"
  --threads "${THREADS}"
  --ref-id "${REF_ID}"
  --target-ids "${TARGET_IDS}"
)

if [[ -n "${OUTGROUP_IDS}" ]]; then
  MAF2BED_ARGS+=( --outgroup-ids "${OUTGROUP_IDS}" )
fi

realignpro maf2bed "${MAF2BED_ARGS[@]}"
echo

if [[ ! -f "${BED_OUT}" ]]; then
  echo "[ERROR] Expected BED not found: ${BED_OUT}" >&2
  exit 2
fi

# ---- Step 3: tsv2fig ----
echo "[STEP 3/3] realignpro tsv2fig"
# Keep old-style output names: conservation_matrix.upstream.pdf / downstream.pdf
realignpro tsv2fig "${TSV_IN}"   --type "${VARIANT_TYPE}"   --highlight "${HIGHLIGHT}"   --outdir "${OUT_DIR}"   --prefix "conservation_matrix"

UP_PDF="${OUT_DIR}/conservation_matrix.upstream.pdf"
DN_PDF="${OUT_DIR}/conservation_matrix.downstream.pdf"

echo
echo "[DONE] Outputs:"
echo "  - MAF : ${MAF_IN}"
echo "  - BED : ${BED_OUT}"
echo "  - TSV : ${TSV_IN}"
echo "  - PDF : ${UP_PDF}"
echo "  - PDF : ${DN_PDF}"
