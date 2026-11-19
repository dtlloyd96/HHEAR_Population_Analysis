# Forest plot of the Figure 3 estimators (Fig3_Estimator_Comparison.R), Figure S18.
# S18 is the unadjusted version, aligned to Figure 3B: same 44 pairs, row order
# (read from figure3_panels.rds) and plotted within and naive values, asserted
# against cohort_regime_pairs.csv. As in 3B, within is scaled by within-cohort SDs
# and naive by total SDs; meta-analysis and Bb use the within-cohort scale.
# A: cohort fixed effects (= Bw) with 95% CI, naive pooled and within-cohort
#    meta-analysis as points on the same line (their CIs are in estimates_long.csv).
# B: the between-cohort term Bb on its own axis, points only: 37 of 44 pairs have
#    two cohorts, where Bb is a line through two means with no cohort-level
#    uncertainty.
# Points beyond an axis limit are drawn at the limit. Open and filled points match
# Figure 3B. Row labels give the number of cohorts only when it is more than two.
# The maternal-age-adjusted version is written as a sensitivity file only.
# Writes Output/Reviewer_Round2/Fig3_Estimators/. Run from the project root.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
OUT <- "Output/Reviewer_Round2/Fig3_Estimators"
res <- fread(file.path(OUT, "estimates_long.csv"))
LIM_A <- 0.3; LIM_B <- 1
OUTCOMES <- c(Gest_Age_Birth = "Gestational age", BMI = "BMI", Birthweight = "Birthweight")
PRETTY   <- c(Gest_Age_Birth = "Gestational Age", BMI = "BMI", Birthweight = "Birthweight")
DIAG <- "Output/Paper1_Output_Submission2/Cohort_Regime_Diagnostic"
fig3B <- readRDS(file.path(DIAG, "figure3_panels.rds"))$decomp$data
stopifnot(is.factor(fig3B$pair), nlevels(fig3B$pair) == 44)
regime <- as.data.table(fig3B)[, .(analyte, outcome, r_within, r_naive)]
LABS  <- c(fe = "Cohort fixed effects (= Bw)", naive = "Naive pooled", meta = "Meta-analysis")
PAL   <- c(fe = "#3B6FB6", naive = "#B0417A", meta = "#D17A22")
SHAPE <- c(fe = 19, naive = 1, meta = 17)
cap_to <- function(x, lim) pmax(pmin(x, lim), -lim)

draw <- function(d) {
  d <- copy(d)
  d[, lab := fifelse(K > 2, paste0(analyte, " (", K, ")"), analyte)]
  key <- function(x) x[, paste0(lab, "|", outcome)]
  d[, row := key(d)]
  d[, pair3 := paste0(analyte, " | ", PRETTY[outcome])]
  stopifnot(all(d$pair3 %chin% levels(fig3B$pair)))
  o <- d[estimator == "fe"][order(match(pair3, levels(fig3B$pair)))]
  d[, row := factor(row, unique(o$row))]
  d[, outcome_lab := factor(OUTCOMES[outcome], OUTCOMES)]
  # alternating row shading, numbered within each outcome
  bands <- unique(d[, .(row, outcome_lab)])[order(row)][, pos := seq_len(.N), by = outcome_lab][pos %% 2 == 0]

  a <- d[estimator %in% names(LABS)]
  a[, `:=`(estimator = factor(estimator, c("naive", "meta", "fe")), r_c = cap_to(r, LIM_A))]
  setorder(a, estimator)                       # fixed effects drawn last, on top
  b <- d[estimator == "bb"][, r_c := cap_to(r, LIM_B)]

  base <- function(p) p +
    geom_rect(data = bands, aes(ymin = pos - 0.5, ymax = pos + 0.5, xmin = -Inf, xmax = Inf),
              inherit.aes = FALSE, fill = "grey95") +
    geom_vline(xintercept = 0, colour = "grey45", linewidth = 0.45)
  thm <- theme_classic(base_size = 13) +
    theme(strip.background = element_blank(), panel.spacing.y = unit(10, "pt"),
          plot.title = element_text(size = 14), axis.ticks.y = element_blank())

  pA <- base(ggplot(a, aes(r_c, row, colour = estimator, shape = estimator))) +
    geom_errorbar(data = a[estimator == "fe"], aes(xmin = r_lo, xmax = r_hi), width = 0, linewidth = 0.7,
                  orientation = "y", show.legend = FALSE) +
    geom_point(size = 2.7, stroke = 0.9) +
    facet_grid(outcome_lab ~ ., scales = "free_y", space = "free_y") +
    scale_y_discrete(labels = function(x) sub("\\|.*$", "", x)) +
    scale_colour_manual(values = PAL, labels = LABS, breaks = names(LABS), name = NULL) +
    scale_shape_manual(values = SHAPE, labels = LABS, breaks = names(LABS), name = NULL) +
    scale_x_continuous(breaks = c(-0.3, -0.2, -0.1, 0, 0.1, 0.2, 0.3), labels = c("-0.3", "-0.2", "-0.1", "0", "0.1", "0.2", "0.3")) +
    coord_cartesian(xlim = c(-LIM_A, LIM_A)) +
    labs(title = "Pooled and within-cohort estimates", x = "Standardized association", y = NULL) +
    thm + theme(legend.position = "bottom", strip.text.y = element_blank(),
                legend.text = element_text(size = 12))

  pB <- base(ggplot(b, aes(r_c, row))) +
    geom_point(shape = 18, colour = "#2A9D8F", size = 3.2) +
    facet_grid(outcome_lab ~ ., scales = "free_y", space = "free_y") +
    scale_x_continuous(breaks = c(-1, 0, 1)) +
    coord_cartesian(xlim = c(-LIM_B * 1.08, LIM_B * 1.08)) +
    labs(title = "Between term (Bb)", x = "Standardized association", y = NULL) +
    thm + theme(axis.text.y = element_blank(), axis.line.y = element_blank(),
                strip.text.y = element_text(angle = -90, face = "bold", size = 13))

  pA + pB + plot_layout(widths = c(3.2, 1)) + plot_annotation(tag_levels = "A")
}

for (s in list(c("maternal age", "adj_maternal_age"), c("none", "unadjusted"))) {
  dat <- res[covariates == s[1]]
  if (s[1] == "none") {                        # plot exactly what Figure 3B plots
    chk <- merge(dat[estimator == "fe"], regime, by = c("analyte", "outcome"))
    stopifnot(nrow(chk) == 44, max(abs(chk$r - chk$r_within)) < 1e-8)
    dat <- merge(dat, regime[, .(analyte, outcome, r_naive)], by = c("analyte", "outcome"))
    dat[estimator == "naive", r := r_naive][, r_naive := NULL]
  }
  p <- draw(dat)
  ggsave(file.path(OUT, paste0("Fig3_Estimator_Forest_", s[2], ".pdf")), p, width = 10, height = 13.5)
  ggsave(file.path(OUT, paste0("Fig3_Estimator_Forest_", s[2], ".png")), p, width = 10, height = 13.5, dpi = 150)
  if (s[1] == "none") {   # Figure S18
    saveRDS(p, file.path(OUT, "panel_estimator_forest.rds"))
    ggsave(file.path(OUT, "SuppFig_S18_Estimator_Comparison.pdf"), p, width = 10, height = 13.5)
  }
}
