#!/bin/bash
## ============================================================================
## rebuild_supplement.sh -- rebuild the supplementary Word documents from the
## current Science-format figure PDFs. No R, no re-analysis: it rasterizes
## Output/Paper1_Submission3/Science_Figures/Supplementary/FigS1..S22.pdf at
## 300 ppi (export_supplementary_figs.sh) for Word, runs the docx build,
## joins the vector PDFs into Output/Paper1_Submission3/Supplementary_Figures/
## All_Supplementary_Figures.pdf and copies the Word files to
## Output/Paper1_Submission3/Supplementary_Document/.
##
##   bash Coding/rebuild_supplement.sh
##
## To change a figure first, re-export it:
##   LANG=en_US.UTF-8 Rscript Coding/Paper1_Submission2/Vis_Code/Science_Supp_Figures.R S4
## To change Supplementary Text numbers, re-knit Coding/Hero_Numbers_Supp.Rmd
## (with LANG=en_US.UTF-8) before running this.
## ============================================================================
set -euo pipefail
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
BUILD="$ROOT/Coding/Supplementary_Build"

bash "$ROOT/Coding/export_supplementary_figs.sh"
for b in build_science.js; do
  echo "== $b"
  (cd "$BUILD" && HHEAR_ROOT="$ROOT" node "$b")
done
# Combined figure PDF: the 22 vector PDFs joined in order (no rasterizing)
SCI="$ROOT/Output/Paper1_Submission3/Science_Figures/Supplementary"
SUB="$ROOT/Output/Paper1_Submission3"
PDFS=(); for i in $(seq 1 22); do PDFS+=("$SCI/FigS$i.pdf"); done
pdfunite "${PDFS[@]}" "$SUB/Supplementary_Figures/All_Supplementary_Figures.pdf"
# Submission 3 copies of the Word files (previous copies kept in _previous/)
mkdir -p "$SUB/Supplementary_Document/_previous"
cp -p "$SUB/Supplementary_Document/"*.docx "$SUB/Supplementary_Document/_previous/" 2>/dev/null || true
cp -p "$ROOT/Output/Supplementary/"HHEAR_Supplementary_Science.docx "$SUB/Supplementary_Document/"
echo "Done:"
ls -lh "$SUB/Supplementary_Document/"*.docx "$SUB/Supplementary_Figures/All_Supplementary_Figures.pdf"
