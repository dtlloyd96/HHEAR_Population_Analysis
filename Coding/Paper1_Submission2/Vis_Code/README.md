# Figures (steps 5–7)

Figures are made in three steps.

1. **Build the panels.** The figure-input scripts (the `.Rmd` files here, the
   assemblers in `Reviewer_Round2/` and `Matrix_Stratified/04_supp_figure_matrix.R`)
   build each panel and save it, mostly to `Output/Paper1_Output_Submission2/Plots/`.
   Run `All_Patchplots.Rmd` last.
2. **Capture the plot objects.** `capture_plots.R` re-runs one figure script
   and saves every ggplot or patchwork object it produces;
   `bash Coding/Paper1_Submission2/Vis_Code/capture_all.sh` does this for every
   supplementary figure. The objects go to
   `Output/Paper1_Submission3/Science_Figures/objects/`.
3. **Export.** `Science_Figures.R` (Figs. 1–3) and `Science_Supp_Figures.R`
   (figs. S1–S22) apply the journal's figure format from `science_fig_style.R`
   (7.25 in wide, Helvetica, 6 and 7 pt text, 10 pt panel letters) and write
   PDF and SVG files. To export only some supplementary figures, name them:
   `Rscript Coding/Paper1_Submission2/Vis_Code/Science_Supp_Figures.R S4 S9`.

Helpers: `category_palette.R` (colours for exposure classes),
`secular_panel.R` (the time-trend panel in Fig. 2A and fig. S11C) and
`combined_realdata_figure.R` (Fig. 3 panels).

Two things to know:
- `Pooling_Simulation_Visualization.Rmd` re-runs `combined_realdata_figure.R`,
  which rewrites the Fig. 3 panels.
- Fig. 1A is a static map image and is not produced by code.

The [top-level README](../../../README.md#figure-and-table-map) maps every
figure to its scripts.
