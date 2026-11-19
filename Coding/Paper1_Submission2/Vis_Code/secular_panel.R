# secular_panel.R -- the median-concentration-over-time panel (conc_year.RData,
# built in Within_Study_SDOH_Compare.Rmd) re-laid out for Science, for a chosen
# set of cohorts. Used by Figure 2A (six cohorts across continents) and
# fig. S11C (all cohorts spanning five or more years). IMHRD, a reference-only
# reconstruction, is always dropped. Needs science_fig_style.R and
# category_palette.R sourced first.
secular_panel <- function(cohorts = NULL, ncol = 3, title = NULL, pt_scale = 0.5) {
  e <- new.env(); load("Output/Paper1_Output_Submission2/Plots/conc_year.RData", envir = e)
  A <- suppressMessages(e$sum_yr +
    scale_fill_manual(values = CATEGORY_PALETTE, na.value = "grey80") +
    scale_colour_manual(values = CATEGORY_PALETTE, na.value = "grey80"))
  all_coh <- setdiff(sort(unique(A$data$Study)), grep("IMHRD", unique(A$data$Study), value = TRUE))
  if (is.null(cohorts)) cohorts <- all_coh
  keep <- function(d) if (is.data.frame(d) && "Study" %in% names(d)) d[d$Study %in% cohorts, ] else d
  A$data <- keep(A$data)
  for (i in seq_along(A$layers)) A$layers[[i]]$data <- keep(A$layers[[i]]$data)
  sci_prepare(A, pt_scale) +
    labs(title = title, subtitle = NULL,
         y = "Median concentration\n(ng/mL or µg/g creatinine; log scale)") +
    facet_wrap(vars(factor(Study, levels = cohorts)), ncol = ncol, scales = "free") +
    scale_x_continuous(breaks = scales::breaks_pretty(3)) +
    theme(legend.position = "bottom",
          panel.spacing.x = unit(8, "pt"), panel.spacing.y = unit(4, "pt")) +
    guides(size = "none", colour = guide_legend(nrow = 2), fill = guide_legend(nrow = 2))
}
