#!/bin/bash
## ============================================================================
## 02_run_stage_a_by_matrix.sh
##
## ITEM 1: run the location model WITHIN each sample matrix.
##
## Deliberately reuses the UNCHANGED published Stage A script
## (Coding/Paper1_Submission2/Location_QuantileReg_Fixed.R) rather than
## reimplementing it. Only the input row set differs, so any difference in the
## result is attributable to matrix stratification and not to method drift.
##
## That script is written for the HPC (setwd to /home/<user>/..., flat
## ../Input/Harmonized_*.RData, output to ../Output/Paper1_Output_Submission2).
## Each stratum therefore gets a scratch tree of that exact shape, with
## Paper1_Submission2 symlinked so `source("Paper1_Submission2/...")` resolves.
##
## Writes ONLY under Output/Matrix_Stratified/. The published Pooled_Output is
## never touched.
## ============================================================================
set -u
export LANG=en_US.UTF-8
export LC_ALL=en_US.UTF-8

ROOT="${HHEAR_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)}"
OUT="$ROOT/Output/Matrix_Stratified"
STAGE_A="$ROOT/Coding/Paper1_Submission2/Location_QuantileReg_Fixed.R"

# "Other" (breast milk + saliva, 2 studies) is excluded: a location model needs
# more than two locations to mean anything.
for CLASS in Urine Blood; do
  IN="$OUT/Input_$CLASS"
  [ -d "$IN" ] || { echo "[skip] $CLASS - no prepped input at $IN"; continue; }

  WORK="$OUT/stageA_$CLASS"
  rm -rf "$WORK"
  mkdir -p "$WORK/Coding" "$WORK/Input" "$WORK/Output/Paper1_Output_Submission2"
  ln -s "$ROOT/Coding/Paper1_Submission2" "$WORK/Coding/Paper1_Submission2"
  cp "$IN/Harmonized_Targeted.RData" "$IN/Harmonized_SDOH.RData" "$WORK/Input/"

  # only change to the published script: where it thinks it is. Since
  # 2026-08-11 the script locates its root from Rscript's --file= argument,
  # in which R encodes spaces as "~+~" ("My Drive" here), so decode that.
  sed -e "s|setwd(\"/home/<user>/HHEAR/HHEAR_Quantreg/Coding\")|setwd(\"$WORK/Coding\")|" \
      -e 's|fa <- argv\[grep("^--file=", argv)\]|fa <- gsub("~+~", " ", argv[grep("^--file=", argv)], fixed = TRUE)|' \
      "$STAGE_A" > "$WORK/stage_a.R"
  grep -q 'gsub("~+~"' "$WORK/stage_a.R" || { echo "[error] root-path patch did not apply"; exit 1; }

  echo "[$(date +%T)] Stage A - $CLASS"
  # HHEAR_SPECS (optional) restricts the sweep, e.g. HHEAR_SPECS=Location for
  # the only model 03_compare_location_contribution.R reads. Unset = full sweep.
  ( cd "$WORK/Coding" && HHEAR_SPECS="${HHEAR_SPECS:-}" Rscript "$WORK/stage_a.R" > "$WORK/stage_a.log" 2>&1 )
  echo "[$(date +%T)]   exit=$?  results=$(find "$WORK/Output" -name '*.RData' | wc -l | tr -d ' ')"
done

echo "[$(date +%T)] done"
