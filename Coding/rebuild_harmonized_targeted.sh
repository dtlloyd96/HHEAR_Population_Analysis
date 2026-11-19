#!/bin/bash
## ============================================================================
## rebuild_harmonized_targeted.sh
##
## Canonical rebuild of Input/Harmonized_Datasets/Harmonized_Targeted.RData
## (2026-11-03). Urine dilution basis:
##   creatinine-corrected cohorts -> ug/g creatinine  (x100 of ng/mL per mg/dL)
##   SG-corrected cohorts         -> ng/mL at SG 1.020
##   everything else              -> ng/mL as reported
## Units is left as reported (many scripts gate on "ng/mL"); the basis is
## recorded per row in the new Dilution_Adjustment column instead.
##
## Steps
##   1. purl HHEAR_Data_Clean_Targeted_V3.Rmd and run it in UGG mode -> Input/
##   2. run it again with no dilution adjustment -> scratch (same row order)
##   3. tag_dilution_adjustment.R compares the two row by row and adds
##      Dilution_Adjustment, so rows a cohort block silently left unadjusted
##      (missing creatinine / SG) are labelled "none" rather than assumed.
##
## Run from Coding/:  bash rebuild_harmonized_targeted.sh
## ============================================================================
set -euo pipefail
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
cd "$(dirname "$0")"

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
Rscript -e 'knitr::purl("HHEAR_Data_Clean_Targeted_V3.Rmd", output = commandArgs(TRUE)[1], quiet = TRUE)' "$TMP/v3.R" >/dev/null

echo "[$(date +%T)] 1/3 ug/g build -> ../Input/Harmonized_Datasets"
HHEAR_DILUTION_ADJUST=UGG HHEAR_HARMONIZED_OUT=../Input/Harmonized_Datasets \
  Rscript "$TMP/v3.R" > "$TMP/ugg.log" 2>&1 || { tail -30 "$TMP/ugg.log"; exit 1; }

echo "[$(date +%T)] 2/3 unadjusted reference build -> scratch"
HHEAR_DILUTION_ADJUST=FALSE HHEAR_HARMONIZED_OUT="$TMP/unadjusted" \
  Rscript "$TMP/v3.R" > "$TMP/unadj.log" 2>&1 || { tail -30 "$TMP/unadj.log"; exit 1; }

echo "[$(date +%T)] 3/3 tagging Dilution_Adjustment"
Rscript tag_dilution_adjustment.R ../Input/Harmonized_Datasets/Harmonized_Targeted.RData \
  "$TMP/unadjusted/Harmonized_Targeted.RData"
echo "[$(date +%T)] done"
