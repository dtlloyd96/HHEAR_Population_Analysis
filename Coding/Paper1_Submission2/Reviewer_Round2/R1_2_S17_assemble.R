# Figure S17 -- Reviewer 1, comment 2: population structure changes inference.
#   A  who is "highly exposed": % above the POOLED HHEAR 90th percentile, serum
#      PFOA/PFOS, by cohort (C5_Reference_Changes_Inference.R)
#   B  PFOS vs child LDL cholesterol: within-cohort fits, cohort means, the
#      between-cohort trend and the naive pooled fit (R1_2_PFOS_LDL_Figure.R)
#   C  slope estimates: per cohort, cohort-adjusted, naive pooled
# The C5 kappa histogram is not included; its summary goes in the caption.
# Run both source scripts first. Writes Output/Reviewer_Round2/R1_2_combined/.
# Run from the project root.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
OUT <- "Output/Reviewer_Round2/R1_2_combined"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
ref <- readRDS("Output/Reviewer_Round2/C5/S17_panel_reference.rds")
pl  <- readRDS("Output/Reviewer_Round2/R1_2_PFOS_LDL/S17_panels_pfos_ldl.rds")
fig <- ref$pA / (pl$pA + pl$pB + plot_layout(widths = c(1.3, 1))) +
  plot_layout(heights = c(0.85, 1.15)) + plot_annotation(tag_levels = "A") &
  # one text scale for all panels (the source scripts used different base sizes)
  theme(plot.title = element_text(size = 12), axis.title = element_text(size = 10),
        axis.text = element_text(size = 9), strip.text = element_text(size = 9.5),
        legend.text = element_text(size = 9), legend.title = element_text(size = 9.5),
        plot.tag = element_text(size = 14, face = "bold"))
ggsave(file.path(OUT, "SuppFig_S17_Population_Structure_Inference.pdf"), fig, width = 11, height = 10.5)
ggsave(file.path(OUT, "SuppFig_S17_Population_Structure_Inference.png"), fig, width = 11, height = 10.5, dpi = 200)
cat("written\n")
