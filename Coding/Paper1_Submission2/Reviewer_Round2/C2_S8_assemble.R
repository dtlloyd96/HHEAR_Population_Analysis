# Figure S8 = single-chemical calendar-year trends (R4 comment 4).
#   A, B  within-cohort trends, legacy vs replacement; direction by class   (C2)
#   C     the same trends by country within TEDDY                           (C2b)
# Run C2_Chemical_Year_Trends.R and C2b_Cross_Country_Year_Trends.R first.
# Writes Output/Reviewer_Round2/C2/SuppFig_S8_Year_Trends_Combined.{pdf,png},
# which export_supplementary_figs.sh maps to S8. Run from the project root.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
ab <- readRDS("Output/Reviewer_Round2/C2/S8_panels_AB.rds")
cd <- readRDS("Output/Reviewer_Round2/C2b_cross_country/S8_panels_CD.rds")
## 2026-11-04 simplified: A curated (legacy + replacement in the same cohort),
## B unchanged, C eight chemicals by TEDDY country; cross-cohort panel dropped.
fig <- (ab$pA | (ab$pB / cd$pC)) + plot_layout(widths = c(1.35, 1)) +
  plot_annotation(tag_levels = "A") &
  theme(plot.tag = element_text(size = 14, face = "bold"))
h <- 9.5
ggsave("Output/Reviewer_Round2/C2/SuppFig_S8_Year_Trends_Combined.pdf", fig, width = 13, height = h, limitsize = FALSE)
ggsave("Output/Reviewer_Round2/C2/SuppFig_S8_Year_Trends_Combined.png", fig, width = 13, height = h, dpi = 150, limitsize = FALSE)
cat("written, height", h, "in\n")
