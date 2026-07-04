#!/usr/bin/env bash
set -euo pipefail

# ReAlignPro `maf2con` demo (coverage-aware) on a tiny chrY-style MAF.
# Demonstrates:
#   - per-chromosome expected depth (N_exp) + call-rate coverage floor,
#   - no whole-block gate (blocks with missing assemblies are still assessed),
#   - --fixed-only and --emit-depth.
# Run from anywhere; outputs go to ./work next to this script.

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
MAF="${HERE}/chrY_demo.maf"
WORK="${HERE}/work"
mkdir -p "${WORK}"

if ! command -v realignpro >/dev/null 2>&1; then
  echo "[ERROR] 'realignpro' not found in PATH. Install the package first (pip/conda)." >&2
  exit 2
fi

echo "[INFO] realignpro version: $(realignpro --version 2>&1 | head -1)"
echo "[INFO] MAF: ${MAF}"
echo

fail() { echo "[FAIL] $1"; echo "  expected: [$2]"; echo "  got:      [$3]"; exit 1; }
check() {  # name  file  expected(newline-stripped)
  local got; got="$(cat "$2")"
  if [ "${got}" = "$3" ]; then echo "[OK]   $1"; else fail "$1" "$3" "${got}"; fi
}

# ---- 1) list IDs ----
IDS="$(realignpro maf2con --ids "${MAF}")"
echo "[1] maf2con --ids -> ${IDS}"
[ "${IDS}" = "hg38, hs1, hs2, hs3, hs4" ] || fail "maf2con --ids" "hg38, hs1, hs2, hs3, hs4" "${IDS}"
echo "[OK]   maf2con --ids"
echo

# ---- 2) default: coverage-aware, auto N_exp (=5), call-rate 0.99 -> floor = 5 ----
# Block 2 (depth 4 < 5) is dropped by the coverage floor; only the full-depth block survives.
realignpro maf2con -i "${MAF}" -o "${WORK}/chrY.default.bed" --ref-id hg38 --target-ids all -t 3
echo "[2] default (call-rate 0.99):"; sed 's/^/    /' "${WORK}/chrY.default.bed"
check "default output" "${WORK}/chrY.default.bed" "$(printf 'chrY\t0\t8')"
echo

# ---- 3) --min-call-rate 0: no relative floor and no whole-block gate ----
# Block 2 is now assessed on its 4 present assemblies; the polymorphic column (chrY:24)
# splits the constrained interval.
realignpro maf2con -i "${MAF}" -o "${WORK}/chrY.allcov.bed" --ref-id hg38 --target-ids all -t 3 --min-call-rate 0
echo "[3] --min-call-rate 0 (no gate):"; sed 's/^/    /' "${WORK}/chrY.allcov.bed"
check "no-floor output" "${WORK}/chrY.allcov.bed" "$(printf 'chrY\t0\t8\nchrY\t20\t24\nchrY\t25\t28')"
echo

# ---- 4) --fixed-only (100% monomorphic) + --emit-depth (N_cov, major-allele freq) ----
realignpro maf2con -i "${MAF}" -o "${WORK}/chrY.fixed.bed" --ref-id hg38 --target-ids all -t 3 \
  --min-call-rate 0 --fixed-only --emit-depth
echo "[4] fixed-only + emit-depth:"; sed 's/^/    /' "${WORK}/chrY.fixed.bed"
awk -F'\t' 'NF!=5 || $5!="1.0000"{print "[FAIL] unexpected emit-depth line: "$0; exit 1}' "${WORK}/chrY.fixed.bed"
[ "$(cut -f1-3 "${WORK}/chrY.fixed.bed")" = "$(cat "${WORK}/chrY.allcov.bed")" ] \
  || { echo "[FAIL] fixed-only intervals differ from the no-floor intervals"; exit 1; }
echo "[OK]   fixed-only + emit-depth"
echo

echo "[DONE] maf2con chrY demo passed."
