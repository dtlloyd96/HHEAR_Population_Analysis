#!/bin/bash
## capture_all.sh -- run every script whose plot object Science_Supp_Figures.R
## re-exports, under capture_plots.R, so that
## Output/Paper1_Submission3/Science_Figures/objects/<key>.rds exists for figs S1-S22.
## Run from the repository root after the figure-input scripts (README, step 5):
##
##   bash Coding/Paper1_Submission2/Vis_Code/capture_all.sh            # all
##   bash Coding/Paper1_Submission2/Vis_Code/capture_all.sh S9 S10     # selected figures
##
## Note: Pooling_Simulation_Visualization.Rmd re-sources combined_realdata_figure.R,
## which rewrites the Figure 3 panels (figure3_panels.rds).
set -euo pipefail
ROOT="${HHEAR_ROOT:-$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)}"
cd "$ROOT"
export LANG=en_US.UTF-8
VIS=Coding/Paper1_Submission2/Vis_Code
RR=Coding/Paper1_Submission2/Reviewer_Round2

# figure -> script that saves its capture key (see the figure map in README.md)
declare -a MAP=(
  "S1  $VIS/All_Patchplots.Rmd"
  "S2  $VIS/LOD_File_Analysis.Rmd"
  "S3  $RR/B3_Analyte_Cohort_Overlap.R"
  "S4  $VIS/All_Patchplots.Rmd"
  "S5  $VIS/TEDDY_All_Figures.Rmd"
  "S6  $VIS/ZIP_All_Figures.Rmd"
  "S7  $VIS/Geography_PCA.Rmd"
  "S8  $RR/R4_5_supp_figures.R"
  "S9  $RR/R4_5_supp_figures.R"
  "S10 $RR/R4_5_supp_figures.R"
  "S11 $VIS/All_Patchplots.Rmd"
  "S12 $RR/C2_S8_assemble.R"
  "S13 $VIS/Between_Study_Power_Output_Visualization.Rmd"
  "S14 $VIS/Power_Output_Visualization.Rmd"
  "S15 $VIS/Pooling_Simulation_Visualization.Rmd"
  "S16 $VIS/Burden_Plots_S2.Rmd"
  "S17 $VIS/Burden_Plots_S2.Rmd"
  "S18 $VIS/Burden_Plots_S2.Rmd"
  "S19 $VIS/Agg_Score_Patchplot_S2.Rmd"
  "S20 Coding/Matrix_Stratified/04_supp_figure_matrix.R"
  "S21 $RR/R1_2_S17_assemble.R"
  "S22 $RR/Fig3_Estimator_Forest.R"
)

want=("$@")
done_scripts=" "
for row in "${MAP[@]}"; do
  fig="${row%% *}"; script="${row##* }"
  if [ ${#want[@]} -gt 0 ] && [[ ! " ${want[*]} " =~ " $fig " ]]; then continue; fi
  [[ "$done_scripts" == *" $script "* ]] && continue      # one run covers several figures
  echo "== $fig  $script"
  Rscript "$VIS/capture_plots.R" "$script"
  done_scripts+="$script "
done
