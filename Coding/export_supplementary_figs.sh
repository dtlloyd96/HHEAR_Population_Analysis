#!/bin/bash
## ============================================================================
## export_supplementary_figs.sh
##
## Writes Output/Supplementary/figs/S1.png ... S22.png -- the images the
## supplementary .docx builds (HHEAR_Manuscript1/Coding/Supplementary_Build)
## embed -- and the combined figure files in Output/Paper1_Submission3/.
##
## Since 2026-11-17 the figures are the Science-format exports
## (Output/Paper1_Submission3/Science_Figures/Supplementary/FigS#.pdf, made by
## Coding/Paper1_Submission2/Vis_Code/Science_Supp_Figures.R from objects captured
## by capture_plots.R). Each vector PDF is rasterized at 300 ppi for Word; the
## vector PDFs themselves are what to upload. The previous source-PDF mapping is
## kept in Output/Paper1_Submission3/_backup_2026-11-17_pre_S9split/.
##
## Run AFTER Science_Supp_Figures.R. Previous PNGs go to figs/_previous/.
## ============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
SCI="$ROOT/Output/Paper1_Submission3/Science_Figures/Supplementary"
FIGS="$ROOT/Output/Supplementary/figs"
N=22

mkdir -p "$FIGS/_previous"
cp -p "$FIGS"/S*.png "$FIGS/_previous/" 2>/dev/null || true
rm -f "$FIGS"/S*.png
missing=0
for i in $(seq 1 $N); do
  src="$SCI/FigS$i.pdf"
  if [ ! -f "$src" ]; then echo "[missing] S$i <- ${src#$ROOT/}"; missing=1; continue; fi
  pdftoppm -png -r 300 -singlefile "$src" "$FIGS/S$i"
  printf '%-8s <- %s  (%s)\n' "S$i.png" "${src#$ROOT/}" "$(date -r "$src" '+%Y-%m-%d %H:%M')"
done
[ $missing -eq 0 ] || { echo "not all figures exported"; exit 1; }
