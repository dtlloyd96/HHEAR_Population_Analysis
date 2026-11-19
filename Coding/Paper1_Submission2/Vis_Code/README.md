# Figures (steps 5–7)

Figures are made in three steps.

1. **Figure-input scripts** (`.Rmd` here, plus the assemblers in
   `Reviewer_Round2/` and `Matrix_Stratified/04_supp_figure_matrix.R`) build
   each panel and save the plot objects to `Output/Paper1_Output_Submission2/Plots/`.
   `All_Patchplots.Rmd` runs last.
2. **Capture.** `capture_plots.R` re-runs one figure script and keeps every
   ggplot/patchwork object it saves; `capture_all.sh` does this for every
   supplementary figure. Objects go to `Output/Paper1_Submission3/Science_Figures/objects/`.
3. **Export.** `Science_Figures.R` (Figs. 1–3) and `Science_Supp_Figures.R`
   (figs. S1–S22; pass figure names such as `S4 S9` to export a subset) apply
   the journal's specifications from `science_fig_style.R` (7.25 in wide,
   Helvetica, 6/7 pt text, 10 pt panel letters) and write PDF and SVG.

Helpers: `category_palette.R` (exposure-class colours), `secular_panel.R`
(the time panel in Fig. 2A and fig. S11C), `combined_realdata_figure.R`
(Fig. 3 panels).

Note that `Pooling_Simulation_Visualization.Rmd` re-sources
`combined_realdata_figure.R`, which rewrites the Fig. 3 panels. Fig. 1A is a
static map image and is not produced by code.

See the figure map in the [top-level README](../../../README.md#figure-and-table-map).
