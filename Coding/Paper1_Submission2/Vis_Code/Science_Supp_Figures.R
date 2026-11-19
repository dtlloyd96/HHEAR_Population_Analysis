# Science_Supp_Figures.R -- supplementary figures re-exported to Science
# specifications (science_fig_style.R): vector PDF + SVG, Helvetica, 6/7 pt text with
# 10 pt bold panel letters, lines >= 0.28 pt, 7.25 in wide, at most 9 in tall.
# Plot objects come from capture_plots.R (Output/Paper1_Submission3/Science_Figures/objects);
# the height follows the aspect of the original figure unless set in TWEAK.
# Writes Output/Paper1_Submission3/Science_Figures/Supplementary/FigS#.pdf and a preview PNG.
# Run from the project root:  Rscript Coding/Paper1_Submission2/Vis_Code/Science_Supp_Figures.R [S1 S2 ...]
suppressPackageStartupMessages({ library(ggplot2); library(patchwork); library(forcats) })   # forcats: captured aes() calls fct_reorder
source("Coding/Paper1_Submission2/Vis_Code/science_fig_style.R")
OBJ <- "Output/Paper1_Submission3/Science_Figures/objects"
OUT <- "Output/Paper1_Submission3/Science_Figures/Supplementary"
SRC <- "Output/Paper1_Submission3/_backup_2026-11-17_supp_sources"   # originals, for aspect ratios

# figure number -> captured object (S9 and S10 are both built from the old S9 panels)
FIGS <- c(S1 = "Demo_Conc_Patch", S2 = "LOD_Patch", S3 = "SuppFig_Analyte_Cohort_Overlap",
          S4 = "Loc_Supp_Patch", S5 = "TEDDY_Figure", S6 = "ZIP_Figure", S7 = "Hier_Summary_Ridge",
          S8 = "Fig_R4_S1_Reported_Sources", S9 = "Fig_R4_S2_Air_Water_Linkages",
          S10 = "Fig_R4_S2_Air_Water_Linkages", S11 = "Time_Supp",
          S12 = "SuppFig_S8_Year_Trends_Combined", S13 = "Between_Study_Location_Power",
          S14 = "Wihtin_Cohort_Power_Patchplot", S15 = "Cohort_Pooling_Simulation_Patchplot",
          S16 = "Burden_RainCloud", S17 = "Burden_to_N", S18 = "Burden_qq", S19 = "FigureS14_All",
          S20 = "SuppFig_Matrix_Confound", S21 = "SuppFig_S17_Population_Structure_Inference",
          S22 = "SuppFig_S18_Estimator_Comparison")

# helpers for the per-figure fixes below
# drop figure-level titles (Science puts them in the legend) and move "A." / "B."
# typed into panel titles onto proper tags
retag <- function(p) {
  p <- sci_walk(p, function(q) {
    if (inherits(q, "ggplot") && is.character(q$labels$title))
      q$labels$title <- sub("^[A-Z]\\.\\s*", "", q$labels$title)
    q
  })
  p + plot_annotation(title = NULL, subtitle = NULL, tag_levels = "A")
}
# replace a string in every character/factor column of every plot and layer data
fix_text <- function(p, from, to) {
  f <- function(d) {
    if (!is.data.frame(d)) return(d)
    for (j in names(d)) {
      if (is.factor(d[[j]])) levels(d[[j]]) <- gsub(from, to, levels(d[[j]]), fixed = TRUE)
      else if (is.character(d[[j]])) d[[j]] <- gsub(from, to, d[[j]], fixed = TRUE)
    }
    d
  }
  sci_walk(p, function(q) {
    if (!inherits(q, "ggplot")) return(q)
    q$data <- f(q$data)
    for (i in seq_along(q$layers)) q$layers[[i]]$data <- f(q$layers[[i]]$data)
    q
  })
}
# set an axis title wherever the current one matches `pattern` (in labs() or in
# the scale's own name, which labs() cannot override)
relab <- function(p, pattern, new, aes = "x") sci_walk(p, function(q) {
  if (!inherits(q, "ggplot")) return(q)
  hit <- function(x) !is.null(x) && !inherits(x, "waiver") && grepl(pattern, paste(deparse(x), collapse = ""))
  if (hit(q$labels[[aes]])) q$labels[[aes]] <- new
  for (i in seq_along(q$scales$scales)) {
    sc <- q$scales$scales[[i]]
    if (aes %in% sc$aesthetics && hit(sc$name))
      q$scales$scales[[i]] <- local({ s <- sc; ggproto(NULL, s, name = new) })
  }
  q
})
# relabel legend entries (colour/fill scales) without touching the data, so
# manual palettes keyed on the original values keep working
relabel_legend <- function(p, fn) sci_walk(p, function(q) {
  if (!inherits(q, "ggplot")) return(q)
  for (i in seq_along(q$scales$scales)) {
    sc <- q$scales$scales[[i]]
    if (any(c("colour", "fill") %in% sc$aesthetics) && sc$is_discrete())
      q$scales$scales[[i]] <- local({ s <- sc; f <- fn; ggproto(NULL, s, labels = function(x) f(x)) })
  }
  q
})
DC_LAB <- "Distributional contribution (mean R²)"
# the TEDDY/ZIP figures: same fixes for both
site_fig <- function(p) {
  p <- relabel_legend(p, function(x) gsub("Georiga", "Georgia", x, fixed = TRUE))
  p <- relab(relab(p, "Contribution to prediction", DC_LAB), "Contribution to prediction", DC_LAB, "y")
  retag(p) & theme(legend.title = element_blank()) &
    guides(fill = guide_legend(nrow = 2))
}

# per-figure adjustments: height (in), point/line scale, and `build`, a function that
# re-assembles the figure from the captured script's own panel objects (attached as
# capture_env) when the original layout cannot be read at 7.25 in.
P <- function(x, k = 0.5) sci_prepare(x, k)
TWEAK <- list(
  # S1: the legend and text call this A (four demographic sub-panels, i-iv) and B
  # (cumulative curves); the demographic block gets most of the page so each
  # study row clears 6-pt text, and B's legend wraps to two rows
  S1 = list(height = 9, build = function() {
    # the panels exceed the capture size cap, so read them from the saved plot files
    pl <- "Output/Paper1_Output_Submission2/Plots/"
    load(paste0(pl, "Demo_Patch.RData")); load(paste0(pl, "cum_conc_plot.RData"))
    # long category names overran the bars at 6 pt; same categories, shorter words
    for (r in list(c("Four_Year_College_Degree", "College"), c("Some_Grad_School", "Some grad"),
                   c("Secondary_Education", "Secondary"), c("Primary_Education", "Primary"),
                   c("Primary_School", "Primary"), c("High_School", "High school"),
                   c("Some_College", "Some college"), c("No_Education", "None"),
                   c("Mixed_SouthAfrican", "Mixed"), c("Over_100k", ">$100k"), c("Over_49k", ">$49k"),
                   c("100_124k", "$100-124k"), c("25_49k", "$25-49k"), c("0_25k", "<$25k")))
      demo_patch <- fix_text(demo_patch, r[1], r[2])
    top <- sci_walk(P(demo_patch), function(q) if (inherits(q, "ggplot"))
      q + sci_dense() + theme(axis.title.y = element_blank()) else q)
    b <- P(cum_conc) + sci_dense() + scale_y_continuous(breaks = c(0, 100)) +
      guides(colour = guide_legend(nrow = 2)) +
      theme(legend.position = "bottom", axis.text.x = element_text(angle = 0, hjust = 0.5))
    wrap_plots(top + plot_layout(tag_level = "new"), b, ncol = 1, heights = c(2.1, 1.9)) +
      plot_annotation(tag_levels = c("A", "i"), tag_sep = "")
  }),
  # S4: two rows (the original reserved four); B and C colour legends repeat the
  # y axis and are dropped; one x-axis title for the same quantity in B to D; D's
  # legend uses readable comparison names
  S4 = list(height = 6.2, build = function() {
    x_lab <- labs(x = "Distributional contribution (mean R²)")
    no_leg <- theme(legend.position = "none")
    d <- P(r2_barplot_loc) + x_lab
    i <- which(vapply(d$scales$scales, function(s) "fill" %in% s$aesthetics, logical(1)))
    relabel <- function(x) c(United_States = "Within US", International = "International",
                             `United_States:International` = "US vs international")[x]
    if (length(i)) d$scales$scales[[i]] <- local({ s <- d$scales$scales[[i]]; ggproto(NULL, s, labels = relabel) })
    d <- d + labs(fill = NULL) + theme(legend.position = "bottom")
    ((P(average_map) + theme(legend.position = "bottom") | P(r2_barplot_pooled_Cohort) + x_lab + no_leg) /
       (P(r2_barplot_matched) + x_lab + no_leg | d)) +
      plot_layout(heights = c(1.1, 1)) + plot_annotation(tag_levels = "A")
  }),
  # S13: the legend has one entry (Location); S14: drop the figure-level title
  # (Science puts it in the legend) and the triplicated one-entry legend, and move
  # the "A." / "B." typed into panel titles onto proper tags
  S13 = list(fun = function(p) p & theme(legend.position = "none")),
  S14 = list(fun = function(p) retag(p) & theme(legend.position = "none")),
  S5 = list(fun = site_fig),
  S6 = list(fun = site_fig),
  S15 = list(height = 8.5, fun = retag),
  # S16/S17: the study legends repeat the y axis (S16); S17's must wrap to fit
  S16 = list(fun = function(p) retag(p) & theme(legend.position = "none")),
  S17 = list(height = 4.6, fun = function(p) retag(p) & guides(colour = guide_legend(ncol = 4)) &
               theme(legend.title = element_blank())),
  # S18 is one pooled QQ plot; a square panel is enough
  S18 = list(height = 4.2, fun = function(p) p + labs(title = "Q-Q plot of the cumulative exposure score")),
  # S19: D and E list 38 cohorts each, so they sit side by side at full height
  # below A to C; panels read from their saved files (larger than the capture cap)
  S19 = list(height = 9, build = function() {
    o <- "Output/Paper1_Output_Submission2/"
    for (f in c("Plots/Burden_Prop_Plot.RData", "R2_Barplot_Burden.RData", "cor_box.RData",
                "cor_bar.RData", "Plots/R2_Barplot_ByStudy.RData")) load(paste0(o, f))
    bot <- theme(legend.position = "bottom", legend.title = element_blank())
    two <- guides(fill = guide_legend(ncol = 2))
    D <- function(x) P(x) + sci_dense()
    d <- D(r2_barplot_agg) + labs(title = "Cumulative score, by cohort", x = DC_LAB) + bot + two
    e <- D(r2_barplot_bystudy) + labs(title = "Individual chemicals, by cohort", x = DC_LAB) + bot + two
    c_ <- relabel_legend(D(burden_prop), function(x) c(`0` = "20 or fewer", `1` = "More than 20")[x]) +
      labs(title = "Score vs share of high exposures", colour = "Exposures measured") +
      theme(legend.position = "bottom")
    top <- D(cor_box) | D(cor_bar) | c_
    wrap_plots(top, d | e, ncol = 1, heights = c(1, 2.1)) + plot_annotation(tag_levels = "A")
  }),
  # S3: 261 analyte rows cannot carry 6-pt names on one page; rows stay grouped by
  # class (the names are in table S6)
  S3 = list(height = 9, build = function() {
    no_y <- theme(axis.text.y = element_blank(), axis.text.y.left = element_blank(),
                  axis.ticks.y = element_blank(), axis.ticks.y.left = element_blank())
    two <- guides(fill = guide_legend(nrow = 2), alpha = guide_legend(nrow = 1))
    (P(pA) + no_y + two + theme(legend.box = "vertical") | (P(pB) / (P(pC) + two))) +
      plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A")
  }),
  # S8 and S11: right-hand legends stacked and long axis titles shortened so the
  # narrow right column is not clipped
  S8 = list(build = function() {
    V <- theme(legend.direction = "vertical")
    (P(pS1A) + labs(x = "Fold difference (95% CI)") + guides(colour = guide_legend(ncol = 2)) |
       ((P(pS1B) + V) / (P(pS1C) + labs(x = "Household water arsenic (µg/L)")))) +
      plot_layout(widths = c(1.25, 1)) + plot_annotation(tag_levels = "A")
  }),
  # S11: A collection years and B mother vs child means side by side, then C the
  # full secular-time grid (every cohort spanning five or more years; Fig. 2A
  # shows six of them). Panels read from their saved files.
  S11 = list(height = 9, build = function() {
    o <- "Output/Paper1_Output_Submission2/Plots/"
    for (f in c("year_hist.RData", "child_adult_mean_bar.RData")) load(paste0(o, f))
    V <- "Coding/Paper1_Submission2/Vis_Code"; source(file.path(V, "category_palette.R")); source(file.path(V, "secular_panel.R"))
    c_ <- secular_panel(ncol = 5, title = "Median concentration over time, all cohorts spanning five or more years")
    top <- (P(year_hist) + sci_dense()) | (P(child_adult_mean_bar) + sci_dense() + theme(legend.position = "bottom"))
    wrap_plots(top, c_, ncol = 1, heights = c(1, 1.25)) + plot_annotation(tag_levels = "A")
  }),
  S12 = list(build = function() {
    V <- theme(legend.direction = "vertical")
    (P(ab$pA) | ((P(ab$pB) + labs(title = "Trend direction by class") + V) / (P(cd$pC) + labs(title = "TEDDY, by country") + V))) +
      plot_layout(widths = c(1.45, 1)) + plot_annotation(tag_levels = "A")
  }),
  # S9: maps 2 x 2 and the four location-level tests 2 x 2 (old S9 A-H)
  S9 = list(height = 8.6, build = function() {
    D <- function(x) P(x) + sci_dense()
    M <- function(x) D(x) + theme(axis.title = element_blank(), legend.key.width = unit(13, "pt"))
    Y <- function(x, lab) D(x) + labs(y = lab)
    maps <- (M(mp$A) | M(mp$B)) / (M(mp$C) | M(mp$D))
    diet <- ((Y(pS2J, "Urinary 3-PBA, location median") | Y(pS2K, "Urinary TCPy, location median")) /
             (Y(pS2L, "Urinary arsenic, location median") | Y(pS2M, "Urinary lead, location median"))) +
            plot_layout(guides = "collect") & theme(legend.position = "bottom")
    wrap_plots(maps, diet, ncol = 1, heights = c(1.1, 1)) + plot_annotation(tag_levels = "A")
  }),
  # S10: household air and drinking water (old S9 I-M), A-B full width each half,
  # C-E below
  S10 = list(height = 6.2, build = function() {
    D <- function(x) P(x) + sci_dense()
    Y <- function(x, lab) D(x) + labs(y = lab)
    ((Y(pS2A, "1-hydroxypyrene, location median (ng/mL)") | Y(pS2B, "1-hydroxypyrene, location median (ng/mL)")) /
       (D(pS2C) | (Y(pS2D, "Urinary arsenic, location median") + labs(x = "P(groundwater As > 10 µg/L)")) |
        (Y(pS2E, "Blood PFAS, location median (ng/mL)") + labs(x = "Water systems with detection (%, state)")))) +
      plot_annotation(tag_levels = "A")
  })
)

orig_aspect <- function(key) {
  f <- file.path(SRC, paste0(key, ".pdf"))
  if (!file.exists(f)) return(NA)
  ps <- sub(".*Page size: *([0-9.]+) x ([0-9.]+).*", "\\1 \\2", grep("Page size", system2("pdfinfo", shQuote(f), stdout = TRUE), value = TRUE))
  v <- as.numeric(strsplit(ps, " ")[[1]]); v[2] / v[1]
}

which <- commandArgs(TRUE); if (!length(which)) which <- names(FIGS)
for (s in which) {
  key <- FIGS[[s]]; f <- file.path(OBJ, paste0(key, ".rds"))
  if (!file.exists(f)) { cat(sprintf("%-4s no captured object (%s)\n", s, key)); next }
  p <- readRDS(f)
  envf <- sub("\\.rds$", "__env.rds", f)                 # helpers the plot refers to
  # attached (not copied into the global env) so they cannot overwrite OUT, s, ...
  while ("capture_env" %in% search()) detach("capture_env")
  theme_set(theme_grey())
  if (file.exists(envf)) {
    ce <- readRDS(envf)
    if (!is.null(ce$.cap_theme)) theme_set(ce$.cap_theme)
    attach(ce, name = "capture_env", warn.conflicts = FALSE)
  }
  tw <- TWEAK[[s]] %||% list()
  # one failing figure is reported and skipped rather than ending the run
  tryCatch({
    p <- if (!is.null(tw$build)) tw$build() else sci_prepare(p, tw$scale %||% 0.5)
    if (!is.null(tw$fun)) p <- tw$fun(p)
    p <- p & sci_tags()
    h <- tw$height %||% min(9, round(SCI_W2 * orig_aspect(key), 2))
    if (is.na(h)) h <- 7
    sci_save(p, paste0("Fig", s), SCI_W2, h, OUT, svg = FALSE)   # Science takes PDF; SVG only for main figures
  }, error = function(e) cat(sprintf("%-4s FAILED: %s\n", s, conditionMessage(e))))
}
