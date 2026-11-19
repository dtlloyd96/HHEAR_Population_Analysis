# Science_Figures.R -- main-text figures re-laid out to Science specifications
# (see science_fig_style.R). Same panel objects as All_Patchplots.Rmd; only size,
# type and layout change. Writes Output/Paper1_Submission3/Science_Figures/.
# Run from the project root:  Rscript Coding/Paper1_Submission2/Vis_Code/Science_Figures.R [fig1 fig2 fig3]
suppressPackageStartupMessages({ library(tidyverse); library(patchwork); library(png) })
VIS <- "Coding/Paper1_Submission2/Vis_Code"
source(file.path(VIS, "science_fig_style.R"))
PL  <- "Output/Paper1_Output_Submission2/Plots"
OUT <- "Output/Paper1_Submission3/Science_Figures/Main"
which <- commandArgs(TRUE); if (!length(which)) which <- c("fig1", "fig2", "fig3")
ld <- function(f) { e <- new.env(); load(file.path(PL, f), envir = e); as.list(e) }

# ================================ Figure 3 =====================================
if ("fig3" %in% which) {
  suppressPackageStartupMessages(library(data.table))
  f3 <- readRDS("Output/Paper1_Output_Submission2/Cohort_Regime_Diagnostic/figure3_panels.rds")
  # 2026-11-18 (editor: "these outlines are small"): significance in A was a thin
  # black ring plus transparency, ~0.4 pt at print size. Now significant estimates
  # are larger, opaque and ringed with a 0.8-pt outline, and non-significant ones
  # small, faded and unringed, so size and opacity carry it even if the ring is
  # missed. B's points and rings are enlarged the same way. C labelled all 44
  # pairs (unreadable); it now labels pairs with I2 >= 50% or a >= 30-fold
  # difference. Full 9-in height (was 8).
  # 2026-11-19: A and B rebuilt from their data to be less busy. Both panels share
  # one row order (outcome groups, then the within-cohort estimate), so a pair reads
  # straight across and its label is printed once, on A, as the exposure alone
  # under an outcome heading. A is coloured by region (US vs international)
  # instead of 15 locations; significant cohort estimates are large and ringed.
  # B greys the concordant pairs so distortions and reversals stand out.
  OL <- c(Birthweight = "Birthweight", Gest_Age_Birth = "Gestational age", BMI = "Maternal BMI")
  b3 <- as.data.table(f3$decomp$data)[, `:=`(out = factor(OL[outcome], levels = OL), exp = analyte)]
  b3 <- b3[order(out, r_within)][, row := factor(paste(exp, out), levels = paste(exp, out))]
  a3 <- as.data.table(f3$spread$data)
  a3[, `:=`(out = factor(OL[outcome], levels = OL), exp = sub(" \\|.*$", "", as.character(pair)))]
  a3[, row := factor(paste(exp, out), levels = levels(b3$row))]
  stopifnot(!anyNA(a3$row), nrow(b3) == 44)
  intl <- c("South_Africa", "Mexico_City", "Thailand", "Tobago", "Suriname", "Central/SouthAmerica")
  a3[, region := factor(ifelse(location %in% intl, "International", "United States"), levels = c("United States", "International"))]
  a3[, sig := factor(ifelse(p_value < 0.05, "P < 0.05", "n.s."), levels = c("P < 0.05", "n.s."))]
  rng <- a3[, .(lo = min(z), hi = max(z)), by = .(row, out)]
  th3 <- theme_classic() + theme(panel.grid = element_blank(), strip.placement = "outside",
    strip.background = element_blank(), strip.text.y.left = element_text(angle = 90, face = "bold"),
    panel.spacing.y = unit(4, "pt"), axis.line.y = element_blank(), axis.ticks.y = element_blank())
  A3 <- ggplot(a3, aes(z, row)) +
    geom_vline(xintercept = 0, linetype = "dashed", colour = "grey60", linewidth = 0.3) +
    geom_segment(data = rng, aes(x = lo, xend = hi, y = row, yend = row), colour = "grey75", linewidth = 0.35, inherit.aes = FALSE) +
    geom_point(aes(fill = region, size = sig, alpha = sig, colour = sig), shape = 21, stroke = 0.6) +
    scale_fill_manual(values = c("United States" = "#3B6EA8", "International" = "#E07B39"), name = NULL) +
    scale_size_manual(values = c("P < 0.05" = 2.3, "n.s." = 1.2), name = NULL) +
    scale_alpha_manual(values = c("P < 0.05" = 1, "n.s." = 0.55), name = NULL) +
    scale_colour_manual(values = c("P < 0.05" = "black", "n.s." = "transparent"), name = NULL) +
    scale_y_discrete(labels = function(x) sub(" (Birthweight|Gestational age|Maternal BMI)$", "", x)) +
    facet_grid(out ~ ., scales = "free_y", space = "free_y", switch = "y") +
    guides(fill = guide_legend(override.aes = list(size = 2.3, colour = "transparent"), order = 1),
           size = guide_legend(override.aes = list(fill = "grey50", colour = c("black", "transparent"), alpha = c(1, 0.55)), order = 2),
           alpha = "none", colour = "none") +
    labs(title = "Cohort estimates disagree in direction", x = "Cohort effect estimate (z-statistic)", y = NULL) + th3
  B3 <- ggplot(b3, aes(y = row)) +
    geom_vline(xintercept = 0, colour = "grey60", linewidth = 0.3) +
    geom_segment(aes(x = r_within, xend = r_naive, yend = row, colour = class), linewidth = 0.8) +
    geom_point(aes(x = r_within), colour = "#14213D", size = 1.7) +
    geom_point(aes(x = r_naive), shape = 21, fill = "white", colour = "#14213D", size = 1.8, stroke = 0.6) +
    scale_colour_manual(values = c("Concordant" = "grey70", "Magnitude distortion" = "#D17A22",
                                   "Sign reversal (Simpson)" = "#C1121F"), name = NULL) +
    facet_grid(out ~ ., scales = "free_y", space = "free_y") +
    labs(title = "Pooling distorts or reverses within-cohort effects",
         x = "Standardized association (filled: within cohort; open: pooled)", y = NULL) +
    th3 + theme(axis.text.y = element_blank(), strip.text.y = element_blank()) + guides(colour = guide_legend(ncol = 1))
  A3 <- sci_prepare(A3, 1); B3 <- sci_prepare(B3, 1)
  f3$div <- sci_prepare(f3$div); f3$varshare <- sci_prepare(f3$varshare)
  lab <- length(f3$div$layers)
  f3$div$layers[[lab]]$data <- function(d) d[d$I2 >= 50 | d$med_fold >= 30, ]
  f3$div$layers[[lab]]$aes_params$size <- 2
  # point labels in C and D at the 6-pt minimum and in Helvetica (they were 4.9-5.7 pt Arial)
  for (nm in c("div", "varshare")) for (i in seq_along(f3[[nm]]$layers))
    if (inherits(f3[[nm]]$layers[[i]]$geom, c("GeomTextRepel", "GeomLabelRepel", "GeomText", "GeomLabel"))) {
      f3[[nm]]$layers[[i]]$aes_params$size <- SCI_SMALL / .pt
      f3[[nm]]$layers[[i]]$aes_params$family <- SCI_FONT }
  p <- A3 + B3 +
    free(f3$div + ggtitle("Heterogeneity and concentration") +
           labs(y = "I² (between-cohort effect heterogeneity, %)"), type = "label", side = "l") +
    free(f3$varshare + ggtitle("Effect distortion and exposure variance") +
           guides(fill = guide_legend(nrow = 2, override.aes = list(size = 2.5))), type = "label", side = "l") +
    plot_layout(design = "AB\nCD", heights = c(1.6, 1)) +
    plot_annotation(tag_levels = "A") &
    sci_tags() & theme(legend.position = "bottom", legend.box = "vertical")
  sci_save(p, "Figure3", SCI_W2, 9, OUT)
}

# ================================ Figure 1 =====================================
# A: composed map (raster, 3399 px across = 469 ppi at 7.25 in); B: exposure
# intervals by location; C: concentrations across matched strata.
if ("fig1" %in% which) {
  map_img <- readPNG(file.path(PL, "Figure1_PanelA.png"))
  # trim the white border, then draw at the image's own aspect ratio (it was
  # stretched to whatever box the layout gave it) with no plot margin
  ink <- map_img[, , 1] < 0.98 | map_img[, , 2] < 0.98 | map_img[, , 3] < 0.98
  rr <- range(which(rowSums(ink) > 0)); cc <- range(which(colSums(ink) > 0))
  map_img <- map_img[rr[1]:rr[2], cc[1]:cc[2], , drop = FALSE]
  map_asp <- nrow(map_img) / ncol(map_img)
  pA <- ggplot() + annotation_raster(map_img, 0, 1, 0, map_asp) +
    scale_x_continuous(limits = c(0, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(0, map_asp), expand = c(0, 0)) +
    coord_fixed() + theme_void() + theme(plot.margin = margin(0, 0, 4, 0)) + sci_tags()
  B <- ld("All_Country_Intervals.RData")$all_plot
  C <- ld("conc_ride_color.RData")$conc_ride_color
  # Panels B and C show the analytes we tried to link to a regional source
  # (Pb, PYR1, PFOA, TCP; Text S5, figs. S9-S10). MBZP had no such link and is
  # dropped (2026-11-17), from every layer of both panels.
  drop_analyte <- function(p, a) {
    keep_rows <- function(d) {
      keep <- rep(TRUE, nrow(d))
      for (col in intersect(c("curr_analyte", "Outcome", "Analyte_Code"), names(d))) keep <- keep & as.character(d[[col]]) != a
      if ("facet_label" %in% names(d)) keep <- keep & !startsWith(as.character(d$facet_label), paste0(a, "-"))
      keep
    }
    f <- function(d, keep) {
      d <- d[keep, , drop = FALSE]
      for (col in names(d)) if (is.factor(d[[col]])) d[[col]] <- droplevels(d[[col]])
      d
    }
    if (is.data.frame(p$data)) p$data <- f(p$data, keep_rows(p$data))
    for (i in seq_along(p$layers)) {
      l <- p$layers[[i]]
      if (!is.data.frame(l$data)) next
      k <- keep_rows(l$data); n <- nrow(l$data)
      l$data <- f(l$data, k)
      # per-row parameters (e.g. error-bar widths) were set as vectors over the old rows
      for (slot in c("aes_params", "geom_params")) {
        ap <- l[[slot]]
        for (nm in names(ap)) if (length(ap[[nm]]) == n && n > 1) ap[[nm]] <- ap[[nm]][k]
        l[[slot]] <- ap
      }
      p$layers[[i]] <- l
    }
    p
  }
  B <- drop_analyte(B, "MBZP"); C <- drop_analyte(C, "MBZP")
  # B relayout: short location names on each panel's axis (the colour legend only
  # repeated them, so it is dropped); fixed small points instead of the NHANES size
  # mapping; 1-2-3 ticks on the original scale (log1p axis) instead of six.
  short_loc <- function(x) {
    x <- dplyr::recode(x, "Northern Karst,Puerto Rico" = "Puerto Rico", "United_States" = "United States",
                       "Boston Massachusetts" = "Boston")
    trimws(gsub("_", " ", sub(",.*$", "", x)))
  }
  B$layers[[1]]$mapping$size <- NULL
  B$layers[[1]]$aes_params$size <- 0.35
  B$scales$scales <- Filter(function(sc) !"size" %in% sc$aesthetics, B$scales$scales)
  B <- B + scale_y_discrete(labels = short_loc) +
    scale_x_continuous(breaks = function(l) {
      v <- c(0, 1, 10, 100, 1000, 1e4, 1e5); b <- log1p(v)[log1p(v) >= l[1] & log1p(v) <= l[2]]
      if (length(b) > 3) b <- b[b != log1p(1)]          # 0 and 1 collide on a log1p axis
      if (length(b) > 3) b <- b[seq(1, length(b), by = 2)]        # every other decade
      b },
                       labels = function(b) format(expm1(b), drop0trailing = TRUE, scientific = FALSE, trim = TRUE)) +
    labs(subtitle = NULL)
  log1p_ticks <- scale_x_continuous(breaks = function(l) {
      v <- c(0, 1, 10, 100, 1000, 1e4, 1e5); b <- log1p(v)[log1p(v) >= l[1] & log1p(v) <= l[2]]
      if (length(b) > 3) b <- b[b != log1p(1)]
      if (length(b) > 3) b <- b[seq(1, length(b), by = 2)]
      b },
    labels = function(b) format(expm1(b), drop0trailing = TRUE, scientific = FALSE, trim = TRUE))
  B <- sci_prepare(B, 0.5) + theme(legend.position = "none", panel.spacing.x = unit(4, "pt"))
  # C relayout (2026-11-18): children on the top row and adults below (4 x 2, was
  # 5 columns with an empty half row); strata in age order (they sorted as text,
  # so 10-16 came before 5-10) and labelled "0-3 F"; ticks on the original scale
  # like B; and a location key, since the colours were otherwise unreadable.
  # Adult strata are mothers, so carry no sex.
  strata_lo <- function(s) as.numeric(sub("^\\((\\d+),.*$", "\\1", s))
  lv <- unique(as.character(C$data$Strata))
  lv <- lv[order(strata_lo(lv), lv)]
  C$data$Strata <- factor(as.character(C$data$Strata), levels = lv)
  C$data$facet_label <- factor(as.character(C$data$facet_label),
    levels = c(paste0(c("Pb", "PFOA", "PYR1", "TCP"), "-Child"), paste0(c("Pb", "PFOA", "PYR1", "TCP"), "-Adult")))
  # shorten the location names in the key only: the manual palette is keyed on the
  # original names, so the data must keep them
  fs <- C$scales$get_scales("fill")
  if (!is.null(fs)) { fs <- fs$clone(); fs$labels <- function(x) short_loc(x); C$scales$scales <-
    c(Filter(function(sc) !"fill" %in% sc$aesthetics, C$scales$scales), list(fs)) }
  C <- C + facet_wrap(vars(facet_label), ncol = 4, scales = "free") +
    scale_y_discrete(labels = function(x) sub(":", " ", gsub("^\\((\\d+),(\\d+)\\)", "\\1–\\2", x))) +
    log1p_ticks +
    labs(title = "Exposure concentrations across matched age and sex strata",
         x = "Concentration (ng/mL or µg/g creatinine; log scale)", y = "Age (years) and sex", fill = NULL)
  C <- sci_prepare(C, 0.45) +
    theme(legend.position = "bottom", legend.key.size = unit(5, "pt"),
          legend.text = element_text(size = SCI_SMALL), panel.spacing.x = unit(8, "pt")) +
    guides(fill = guide_legend(nrow = 3))
  map_h <- SCI_W2 * map_asp
  # free(): A is not aligned to B/C's axis-label margins, so the map uses the full width
  p <- free(pA) / B / C + plot_layout(heights = c(map_h, 2.2, 2.5)) +
    plot_annotation(tag_levels = "A") & sci_tags()
  # A's tag sits on the map's empty top-left corner so the map can use the full width
  p[[1]] <- p[[1]] + theme(plot.tag.location = "plot", plot.tag.position = c(0.005, 0.99))
  sci_save(p, "Figure1", SCI_W2, 9.0, OUT)
}

# ================================ Figure 2 =====================================
# A time series by cohort (full width); B year range vs year contribution and
# C child vs maternal concentration side by side; D child/mother log ratios (full
# width); E education and F race/ethnicity in CCREOH side by side.
if ("fig2" %in% which) {
  source(file.path(VIS, "category_palette.R"))
  recolor_cat <- function(p) suppressMessages(
    p + scale_fill_manual(values = CATEGORY_PALETTE, na.value = "grey80") +
        scale_colour_manual(values = CATEGORY_PALETTE, na.value = "grey80"))
  B  <- recolor_cat(ld("year_outcome_scatter.RData")$year_outcome_scatter)
  C  <- ld("mat_child_scatter.RData")$mat_child_scatter
  D  <- recolor_cat(ld("vi_plot.RData")$vi_plot)
  e  <- new.env(); load("Output/Paper1_Output_Submission2/madres_plots.RData", envir = e)
  E  <- e$madres_plots[["CCREOH__Mat_Edu"]] + labs(colour = "Education", title = "Maternal education, CCREOH (Suriname)", subtitle = NULL) +
        theme(axis.title.y = element_blank())
  F_ <- e$madres_plots[["CCREOH__Mat_Race"]] + labs(colour = "Race/ethnicity", title = "Race/ethnicity, CCREOH (Suriname)", subtitle = NULL) +
        theme(axis.title.y = element_blank())
  prep <- function(p, k = 0.5) sci_prepare(p, k) + theme(legend.position = "bottom")
  # A: six cohorts across four continents (US: SEARCH, OBSS, CHARGE; Europe and
  # US: TEDDY; South Africa: DCHS; Chile: GOC), chosen 2026-11-18 because the 15
  # cohorts spanning five or more years were unreadable at this size. All 15 are
  # in fig. S11C. IMHRD (reference only) is dropped inside secular_panel().
  source(file.path(VIS, "secular_panel.R"))
  A <- secular_panel(c("SEARCH", "OBSS", "CHARGE", "TEDDY", "DCHS", "GOC"), ncol = 3,
                     title = "Median concentration over time, selected cohorts across continents")
  # B (rebuilt 2026-11-18): every chemical-cohort sample-year contribution (small
  # points) under the class mean per cohort that the saved panel plotted (outlined),
  # on a square-root axis. Chemical values are the same tau-weighted pseudo-R2 as
  # Within_Study_SDOH_Compare.Rmd, read from the within-cohort year models. The
  # dashed line and the subtitle are the class-mean regression (the saved
  # subtitle had been typed in as 0.18 from the March data).
  # B and D use the same category colours as A, whose legend serves all three
  cls <- B$data
  pull_sy <- function(f) {
    e <- new.env(); load(f, envir = e); L <- get(ls(e)[1], envir = e)
    bind_rows(lapply(names(L), function(coh) bind_rows(lapply(names(L[[coh]]), function(an) {
      r <- L[[coh]][[an]][["Sample_Collection_Year"]]
      if (is.null(r) || !is.data.frame(r$Regression_Results)) return(NULL)
      r$Regression_Results %>% filter(Variable == "Sample_Collection_Year") %>%
        transmute(Study = coh, Outcome = an, Tau, pseudo_R)
    }))))
  }
  RES <- "Output/Paper1_Output_Submission2"
  cat_map <- read.csv("Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv") %>%
    filter(Category %in% unique(cls$Category)) %>% distinct(Outcome = Code, Category)
  chem <- bind_rows(pull_sy(file.path(RES, "By_Study_Single_Exp_Results_Child.RData")),
                    pull_sy(file.path(RES, "By_Study_Single_Exp_Results_Mat.RData"))) %>%
    filter(!Study %in% c("DISCOVER", "VAPS")) %>%
    inner_join(cat_map, by = "Outcome") %>%
    group_by(Study, Category, Outcome) %>%
    summarise(r2 = weighted.mean(pmax(pseudo_R, 0), exp(-((Tau - 0.5)^2) / (2 * 0.15^2)), na.rm = TRUE),
              .groups = "drop") %>%
    filter(!is.na(r2)) %>%
    inner_join(distinct(ungroup(cls), Study, year_range), by = "Study")
  B_lm  <- lm(mean_r2 ~ year_range, data = cls)
  B_line <- data.frame(year_range = seq(min(cls$year_range), max(cls$year_range), length.out = 100))
  B_line$r2 <- predict(B_lm, B_line)
  set.seed(1)
  B <- ggplot(chem, aes(year_range, r2)) +
    geom_jitter(aes(colour = Category), width = 0.25, height = 0, size = 0.9, alpha = 0.6, stroke = 0) +
    geom_point(data = cls, aes(y = mean_r2, fill = Category), shape = 21, size = 1.9,
               colour = "black", stroke = 0.25) +
    geom_line(data = B_line, colour = "black", linetype = "dashed", linewidth = 0.5) +
    scale_colour_manual(values = CATEGORY_PALETTE) + scale_fill_manual(values = CATEGORY_PALETTE) +
    scale_y_sqrt(breaks = c(0, 0.01, 0.05, 0.1, 0.2, 0.4, 0.6)) +
    labs(title = "Year range vs sample-year contribution",
         subtitle = sprintf("Adjusted R² = %.2f", summary(B_lm)$adj.r.squared),
         x = "Year range", y = "Sample-year contribution\n(square-root scale)") +
    theme_classic()
  B <- prep(B, 1) + theme(legend.position = "none")   # sizes above are already final
  C <- prep(C, 0.4) + labs(subtitle = NULL)
  # exposure names read vertically instead of at 45 degrees (less depth)
  D <- prep(D, 0.4) + theme(legend.position = "none",
                            axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
  few_x <- scale_x_continuous(breaks = scales::breaks_pretty(3))
  E <- prep(E, 0.4) + few_x + theme(legend.title = element_blank())
  F_ <- prep(F_, 0.4) + few_x + guides(colour = guide_legend(nrow = 2)) +
    theme(legend.title = element_blank(), plot.margin = margin(3, 3, 3, 18))   # room for the F tag
  p <- A / (B | C) / D / (E | F_) +
    plot_layout(heights = c(2.7, 1.8, 1.5, 1.8)) +
    plot_annotation(tag_levels = "A") & sci_tags()
  sci_save(p, "Figure2", SCI_W2, 9.0, OUT)
}
