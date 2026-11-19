#!/bin/bash
## ============================================================================
## render_chain.sh -- re-render a named chain of analysis/figure scripts in
## order, each from its own directory (their paths are relative), logging each
## to Output/Dilution_Unadjusted/logs/. Continues past a failure and prints a
## summary, so one broken step does not hide the rest.
##
## Chains (2026-11-03, ug/g propagation):
##   fig3     ExWAS -> regime diagnostic -> Figure 3 panels -> Figure 3 numbers
##   matched  age/sex matched models -> matched comparison -> matched intervals
##   direct   figure inputs that read only the harmonized data
##   pooled   figure inputs that read Pooled_Output (run after the full sweep)
##
## Usage: bash render_chain.sh <chain> [<chain> ...]
## ============================================================================
set -u
export LANG=en_US.UTF-8 LC_ALL=en_US.UTF-8
ROOT="${HHEAR_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
P="$ROOT/Coding/Paper1_Submission2"
LOG="$ROOT/Output/Dilution_Unadjusted/logs"; mkdir -p "$LOG"

chain() {
  case "$1" in
    fig3)    echo "Targeted_ExWAS_Replication.Rmd cohort_regime_diagnostic.R Vis_Code/combined_realdata_figure.R Figure3_Results_Numbers.Rmd" ;;
    matched) echo "Age_Sex_Match.Rmd Vis_Code/Age_Sex_Match_Compare.Rmd Vis_Code/Matched_Intervals.Rmd" ;;
    direct)  echo "Vis_Code/Geography_PCA.Rmd Vis_Code/Cumulative_Conc_Plot_Map.Rmd Vis_Code/Figure2DemoPlot.Rmd Vis_Code/Mother_Child_Compare.Rmd Vis_Code/Within_Study_SDOH_Compare.Rmd" ;;
    pooled)  echo "Vis_Code/Location_Compare.Rmd Vis_Code/Intervals_ByCountry.Rmd Vis_Code/SDOH_Intervals.Rmd" ;;
    burden)  echo "Exposure_Burden_Creation_Perc_S2.Rmd ByCohort_BurdenExWAS_S2.Rmd Vis_Code/Burden_Plots_S2.Rmd Vis_Code/BurdenByStudy_CompareR2_S2.Rmd Vis_Code/Agg_Score_Patchplot_S2.Rmd" ;;
    corr)    echo "Correlations_EffectiveVars_S2.Rmd" ;;
    patch)   echo "Vis_Code/All_Patchplots.Rmd" ;;
    *) echo "unknown chain $1" >&2; exit 1 ;;
  esac
}

summary=()
# Targeted_ExWAS_Replication.Rmd resolves out_root by searching upward for an
# EXISTING directory and falls back to its own folder if none is found, which
# silently writes the tree to Coding/Paper1_Submission2/Health_Outcome_Analysis/.
# Make sure the real target exists (after the old tree has been moved aside).
mkdir -p "$ROOT/Health_Outcome_Analysis/Output_Targeted_ExWAS_R"
for c in "$@"; do
  for f in $(chain "$c"); do
    dir="$P/$(dirname "$f")"; base="$(basename "$f")"; log="$LOG/render_${base%.*}.log"
    echo "[$(date +%T)] $c: $f"
    if [[ "$base" == *.Rmd ]]; then
      ( cd "$dir" && Rscript -e "rmarkdown::render('$base', quiet = TRUE, envir = new.env())" > "$log" 2>&1 )
    else
      ( cd "$dir" && Rscript "$base" > "$log" 2>&1 )
    fi
    rc=$?
    summary+=("$(printf '%-4s %-45s %s' "$rc" "$f" "$(date +%T)")")
    [ $rc -ne 0 ] && echo "   FAILED (exit $rc) -- tail of $log:" && tail -5 "$log"
  done
done
echo; echo "exit  script                                         finished"; printf '%s\n' "${summary[@]}"
