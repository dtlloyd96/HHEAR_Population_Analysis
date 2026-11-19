#!/bin/bash
## ============================================================================
## run_location_unadjusted.sh
##
## Re-runs the published between-cohort Location model on the UNADJUSTED
## harmonized targeted data (no creatinine / SG correction), which puts every
## urine cohort on the same ng/mL scale. See HHEAR_Data_Clean_Targeted_V3.Rmd,
## DILUTION ADJUSTMENT SWITCH, for why.
##
## Reuses the unchanged Location_QuantileReg_Fixed.R in a scratch tree, exactly
## like Coding/Matrix_Stratified/02_run_stage_a_by_matrix.sh, so only the input
## rows differ from the published run. Writes ONLY under
## Output/Dilution_Unadjusted/.
##
## Prerequisite: Output/Dilution_Unadjusted/harmonized_unadjusted/ built by
##   (cd Coding && HHEAR_DILUTION_ADJUST=FALSE \
##     HHEAR_HARMONIZED_OUT=../Output/Dilution_Unadjusted/harmonized_unadjusted \
##     Rscript <purled HHEAR_Data_Clean_Targeted_V3.R>)
##
## Usage: bash run_location_unadjusted.sh [SPECS] [VARIANT]
##   SPECS   default Location
##   VARIANT unadjusted (default) or adjusted -- the adjusted rebuild gives an
##           exact like-for-like baseline from the same code
## ============================================================================
set -u
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8

ROOT="${HHEAR_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
OUT="$ROOT/Output/Dilution_Unadjusted"
STAGE_A="$ROOT/Coding/Paper1_Submission2/Location_QuantileReg_Fixed.R"
SPECS="${1:-Location}"
VARIANT="${2:-unadjusted}"

WORK="$OUT/stageA_$VARIANT"
rm -rf "$WORK"
mkdir -p "$WORK/Coding" "$WORK/Input" "$WORK/Output/Paper1_Output_Submission2"
ln -s "$ROOT/Coding/Paper1_Submission2" "$WORK/Coding/Paper1_Submission2"
cp "$OUT/harmonized_$VARIANT/Harmonized_Targeted.RData" "$WORK/Input/"
cp "$ROOT/Input/Harmonized_Datasets/Harmonized_SDOH.RData" "$WORK/Input/"

# R encodes spaces in --file= as "~+~" ("My Drive"); decode before the script
# walks up from its own path to find the project root.
sed -e 's|fa <- argv\[grep("^--file=", argv)\]|fa <- gsub("~+~", " ", argv[grep("^--file=", argv)], fixed = TRUE)|' \
    "$STAGE_A" > "$WORK/stage_a.R"
grep -q 'gsub("~+~"' "$WORK/stage_a.R" || { echo "[error] root-path patch did not apply"; exit 1; }

echo "[$(date +%T)] Stage A ($VARIANT), specs=$SPECS"
( cd "$WORK/Coding" && HHEAR_SPECS="$SPECS" Rscript "$WORK/stage_a.R" > "$WORK/stage_a.log" 2>&1 )
echo "[$(date +%T)] exit=$?  results=$(find "$WORK/Output" -name '*Quantile_Results.RData' | wc -l | tr -d ' ')"
