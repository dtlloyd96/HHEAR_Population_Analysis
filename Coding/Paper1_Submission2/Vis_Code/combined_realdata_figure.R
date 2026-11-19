## ============================================================================
## combined_realdata_figure.R   (run from Coding/Paper1_Submission2/Vis_Code)
##
## Builds ONE cohesive real-data figure to follow the cohort-pooling simulation,
## combining the existing MS heterogeneity panels with the within/pooled/
## ecological decomposition, all on the SAME slice: maternal panel, ADJUSTED
## model, log1p scale.
##
##   A (spread)      : per-cohort effect (z-stat) by study LOCATION; shows that
##                     cohorts disagree (direction discordance).      [MS pipeline]
##   B (decomposition): within-cohort vs naive-pooled vs ecological standardized
##                     effect; shows that pooling distorts/reverses the
##                     contextual effect.                              [new]
##   C (divergence)  : I^2 vs between-cohort concentration fold-change; shows the
##                     heterogeneity is NOT explained by exposure spread. [MS pipeline]
##   D (varshare)    : between-cohort SHARE of exposure variance vs how far pooling
##                     moves the estimate -- the mechanism behind A-C, and unlike
##                     B's ecological slope it is estimable at 2 cohorts.
##                     Added 2026-08-10, replacing the budget-allocation
##                     simulation as panel D (see section 2c).          [new]
##
## Row order in A and B is BLOCKED BY OUTCOME, then sorted by each panel's own
## metric within the block (see the note above the lvlA/lvlB assignment).
##
## Saves:
##   Output/.../Plots/Final_Patchplots/MS_Fig_RealData_Combined.pdf
##   Output/.../Plots/Final_Patchplots/Supp_Pooling_Budget_Tradeoff.pdf  (retired D)
##   Output/.../Cohort_Regime_Diagnostic/figure3_panels.rds  (spread/decomp/div/
##                                                            varshare/cost)
##   Output/.../Cohort_Regime_Diagnostic/panelD_variance_share.csv
##   Output/.../Cohort_Regime_Diagnostic/combined_real_fig.rds   (for the Rmd)
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr); library(tidyr); library(forcats); library(ggplot2)
  library(ggrepel); library(patchwork); library(stringr); library(purrr)
  library(arrow); library(scales); library(metafor); library(readr)
})

# std_study(), to_ng_per_mL(). Added 2026-08-11; see CODE_AUDIT_2026-08-11.md.
local({
  cands <- c("../harmonization_helpers.R", "harmonization_helpers.R",
             "Coding/Paper1_Submission2/harmonization_helpers.R")
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) stop("Cannot locate harmonization_helpers.R")
  source(hit[[1]], local = FALSE)
})

exwas_root <- "../../../Health_Outcome_Analysis/Output_Targeted_ExWAS_R"
out_dir    <- "../../../Output/Paper1_Output_Submission2"
input_dir  <- "../../../Health_Outcome_Analysis/Input"
final_dir  <- file.path(out_dir, "Plots", "Final_Patchplots")
diag_dir   <- file.path(out_dir, "Cohort_Regime_Diagnostic")
MIN_STUDY_N <- 20

## ============================================================================
## 1. MS pipeline (reproduced from ExWAS_Location_Heterogeneity.Rmd)
## ============================================================================
result_paths <- list.files(exwas_root, pattern = "exposure_results_table\\.parquet$",
                           recursive = TRUE, full.names = TRUE)
result_paths <- result_paths[grepl("scale=", result_paths)]

results <- bind_rows(lapply(result_paths, read_parquet)) %>%
  # CCREOH EXCLUSION REMOVED 2026-08-10. The filter was `!study %in% c("CCREOH")`,
  # copied from ExWAS_Location_Heterogeneity.Rmd, with no comment at any of the
  # four sites that carried it and no rationale in any handoff. Checked before
  # removing: CCREOH is not an outlier here -- max |z| 1.87 (the mildest of the
  # large cohorts), 0% of its estimates significant, unremarkable SEs -- so the
  # obvious "it would dominate the axis" reason does not hold. Dropping the
  # largest maternal cohort (n = 656) from a cross-cohort heterogeneity figure
  # needs a stated reason, and there was none.
  # Effect of including it: 35 -> 44 pairs, 83 -> 101 estimates, 11 -> 13 cohorts
  # (CCREOH also gives THS a partner, so THS clears the 2-cohort minimum too).
  # Panel B is essentially unchanged (sign reversal 45.7% -> 45.5%); panel A's
  # direction discordance drops 68.6% -> 61.4%, i.e. the exclusion had been
  # inflating the headline number.
  # MATERNAL ONLY (reverted 2026-08-11 at the author's request). Opening the
  # figure to both panels was tried and rolled back: it took panel A from 44 to
  # 61 rows and weakened every headline (discordance excess over the permutation
  # null +5.4pp -> +3.4pp, panel C rho 0.11 -> 0.03) while mixing maternal and
  # child estimates of the same analyte, which are different
  # exposure-outcome relationships. The child-only cardiometabolic outcomes
  # (GOALS, NASH, SEARCH, CHARGE, MDUC) are therefore out of this figure again.
  # `panel` is still carried through the joins below -- with one panel it is
  # constant, but the keys stay correct if the filter is ever lifted.
  filter(status == "ok", model == "Adjusted", panel == "mat",
         !outcome %in% c("Asthma", "Severe_Asthma"),
         !is.na(beta)) %>%
  mutate(
    # Prefer the se the ExWAS now records directly (2026-08-11, audit 3.8). The
    # /3.92 reconstruction assumes a 1.96 critical value while the linear CIs use
    # qt(0.975, n-p), so it overestimates se by ~6% at df = 20 and biases the
    # rma() fits below. Kept as a fallback for older result trees.
    se = if ("se" %in% names(.)) {
      coalesce(se, case_when(
        outcome_type == "continuous" ~ (ci_high - ci_low) / 3.92,
        TRUE ~ (log(pmax(ci_high, 1e-9)) - log(pmax(ci_low, 1e-9))) / 3.92))
    } else {
      case_when(
        outcome_type == "continuous" ~ (ci_high - ci_low) / 3.92,
        TRUE ~ (log(pmax(ci_high, 1e-9)) - log(pmax(ci_low, 1e-9))) / 3.92)
    },
    outcome_clean = case_when(
      outcome == "Gest_Age_Birth" ~ "Gestational Age",
      outcome == "Preterm_Birth"  ~ "Preterm Birth",
      TRUE ~ outcome)
  ) %>%
  filter(is.finite(beta), is.finite(se), se > 0)

load(file.path("../../../Input/Harmonized_Datasets/Harmonized_SDOH.RData"))
load(file.path("../../../Input/Harmonized_Datasets/Harmonized_Targeted.RData"))

collapse_country <- function(loc, cty) {
  us <- c("United_States","NorthEast_US","New_Hampshire","Washington_State",
          "Midwestern_US","Denver","Baltimore","Sacramento","Southern_California",
          "Syracuse","Michigan","California","SanFrancisco","Massachusetts",
          "Boston","Los_Angeles","New_York",
          # Added 2026-08-10: ATLAA (Atlanta, Georgia). Its SDOH Country is NA, so
          # collapse_country falls through to this list; without the entry the
          # cohort was classified NON-US and every pair containing it was
          # mis-coloured in panel C ("US + international" instead of "All US
          # cohorts"). SWSD (Thailand) and GRAPHS (Accra) are correctly
          # international and need no entry.
          "Atlanta")
  if (!is.na(cty) && nzchar(cty)) return(cty)
  if (!is.na(loc) && loc %in% us) return("United_States")
  loc
}

study_loc <- SDOH_Data_All %>%
  # std_study() is the canonical cohort key (2026-08-11, audit section 3.5).
  mutate(study = std_study(Study),
         Location = as.character(Location), Country = as.character(Country),
         country_val = map2_chr(Location, Country, collapse_country)) %>%
  count(study, Location, country_val, name = "n") %>%
  group_by(study) %>% slice_max(n, n = 1, with_ties = FALSE) %>% ungroup() %>%
  transmute(study, location = Location, country = country_val)

# BOTH PANELS (2026-08-11). Was maternal-only. Keyed by panel as well as
# study/exposure, because the median concentration of an analyte differs between
# mothers and children -- collapsing them would blend two distributions.
conc_sum <- bind_rows(
    Targeted_Data_All %>% filter(!is.na(Mat_PID),   Mat_PID   != "") %>% mutate(panel = "mat"),
    Targeted_Data_All %>% filter(!is.na(Child_PID), Child_PID != "") %>% mutate(panel = "child")
  ) %>%
  transmute(panel, study = std_study(Study),
            exposure = as.character(Analyte_Code),
            conc = suppressWarnings(as.numeric(Concentration))) %>%
  filter(!is.na(conc), is.finite(conc), conc >= 0) %>%
  group_by(panel, study, exposure) %>%
  summarise(med = median(conc, na.rm = TRUE), .groups = "drop")

run_meta <- function(df) {
  fit <- tryCatch(rma(yi = df$beta, sei = df$se, method = "REML"),
                  error = function(e) tryCatch(rma(yi = df$beta, sei = df$se, method = "DL"),
                                               error = function(e) NULL))
  if (is.null(fit)) return(tibble())
  tibble(k = fit$k, I2 = fit$I2, tau2 = fit$tau2, QE = fit$QE, QEp = fit$QEp,
         pooled = as.numeric(fit$beta))
}

het_tbl <- results %>%
  filter(scale == "log1p") %>%
  # panel added 2026-08-11: a maternal and a child estimate for the same analyte
  # are different relationships and must not be meta-analysed together.
  group_by(panel, outcome_clean, exposure) %>%
  filter(n_distinct(study) >= 2) %>%
  group_modify(~ run_meta(.x)) %>% ungroup() %>%
  mutate(facet = paste0(panel, " | ", outcome_clean, " | ", exposure))

# TYPE SCALED UP 2026-08-10. Figure3_All.pdf is drawn at 20 x 18 inches, and the
# original 10pt base left every label small relative to the canvas. axis.text is
# held at 12 deliberately: panels A and B carry ~44 y-axis rows across roughly
# 10.5 inches, i.e. ~17pt of vertical space per row, so larger text collides.
ms_theme <- theme_classic(base_size = 14) +
  theme(axis.text = element_text(size = 12, color = "black"),
        axis.title = element_text(size = 15),
        strip.text = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 13),
        strip.background = element_blank(), legend.position = "bottom",
        panel.grid.major.y = element_line(color = "grey93", linewidth = 0.25),
        # Was element_blank(), which suppressed panel titles even when set.
        # Panels that pass no `title` to labs() still show nothing, so this only
        # affects the panels that ask for one (C and D).
        plot.title = element_text(size = 17, face = "bold", color = "#14213D",
                                  margin = margin(b = 8)))

loc_pal <- c("Boston"="#08306B","New_Hampshire"="#2171B5","New_York"="#6BAED6",
             "California"="#238B45","Los_Angeles"="#74C476","Mexico_City"="#D7301F",
             "Central/SouthAmerica"="#FD8D3C","South_Africa"="#8C510A",
             # Added 2026-08-10 for the three new cohorts. Keys are the RAW
             # SDOH Location values ("Accra", "Thailand"), NOT the display names
             # used in Intervals_ByCountry/Matched_Intervals ("Accra,Ghana",
             # "Chiang Mai,Thailand") -- this script joins Location directly.
             # Hex values reused from the Figure 1 city palette so a cohort keeps
             # one colour across figures.
             # Without these, scale_fill_manual sent them to na.value grey and
             # dropped them from the legend, so ATLAA -- which contributes 8 of
             # the 9 new pairs -- was an unlabelled grey dot.
             "Atlanta"="#4292C6","Accra"="#DFC27D","Thailand"="#762A83",
             # Michigan (CEHC) was the ONLY location left falling through to
             # na.value once the three new cohorts were added, so it is included
             # here to make panel A fully coloured. Green to match the other US
             # entries (the palette groups US = blues/greens, Mexico = red,
             # Central/South America = orange, Africa = brown/tan, Asia = purple).
             # These 12 keys now cover every location that actually reaches
             # panel A; other cohorts in Harmonized_SDOH are absent from this
             # panel entirely because they contribute no >=2-cohort pair.
             "Michigan"="#41AB5D",
             # Added 2026-08-10 with the CCREOH exclusion removal: CCREOH is
             # Suriname and THS is Tobago, neither of which had an entry. Teal and
             # magenta are the only hues not already spoken for -- the palette's
             # loose geographic grouping (US = blues/greens, Mexico = red,
             # Central/South America = orange, Africa = brown/tan, Asia = purple)
             # has no free slot in the Americas that would stay distinguishable
             # from Mexico_City and Central/SouthAmerica.
             "Suriname"="#35978F","Tobago"="#C51B7D")

## ---- Panel A data: cohort effect spread (z) by location ---------------------
## Shared label so A and B use identical row names; the row ORDER is derived
## below from panel B's within-cohort effect and applied to both panels.
pretty_out <- function(o) dplyr::recode(o,
  Gest_Age_Birth = "Gestational Age", Preterm_Birth = "Preterm Birth", .default = o)

strip_df <- results %>%
  filter(scale == "log1p") %>%
  mutate(z = beta / se,
         pair = paste0(exposure, " | ", pretty_out(outcome))) %>%
  left_join(study_loc %>% select(study, location), by = "study") %>%
  group_by(pair) %>% filter(n_distinct(study) >= 2) %>%
  mutate(spread = max(z) - min(z), discordant = any(z > 0) & any(z < 0)) %>%
  ungroup() %>%
  mutate(sig_lab = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"))
discord_pct <- strip_df %>% distinct(pair, discordant) %>% summarise(p = mean(discordant)) %>% pull(p)
# p_spread is built further down, once the shared row order is known.

## ---- Panel C: I^2 vs concentration divergence -------------------------------
div_tbl <- results %>% filter(scale == "log1p") %>%
  distinct(panel, outcome_clean, exposure, study) %>%
  left_join(conc_sum, by = c("panel", "study", "exposure")) %>%
  group_by(panel, outcome_clean, exposure) %>%
  summarise(med_fold = max(med, na.rm = TRUE) / pmax(min(med, na.rm = TRUE), 1e-9),
            .groups = "drop") %>% filter(is.finite(med_fold))
geo_cat_tbl <- results %>% filter(scale == "log1p") %>%
  distinct(panel, outcome_clean, exposure, study) %>%
  left_join(study_loc %>% select(study, country), by = "study") %>%
  group_by(panel, outcome_clean, exposure) %>%
  summarise(n_us = sum(country == "United_States", na.rm = TRUE), n_tot = dplyr::n(), .groups = "drop") %>%
  mutate(geo_cat = case_when(n_us == n_tot ~ "All US cohorts", n_us == 0 ~ "All international",
                             TRUE ~ "US + international"))
geo_pal <- c("All US cohorts"="#2171B5","US + international"="#7B5EA7","All international"="#D7301F")
# `panel` MUST be in both join keys. Without it each het_tbl row matches the
# maternal AND the child div_tbl row for the same analyte+outcome, doubling the
# panel C point cloud and corrupting the Spearman rho (2026-08-11).
het_div <- het_tbl %>%
  left_join(div_tbl, by = c("panel", "outcome_clean", "exposure")) %>%
  left_join(geo_cat_tbl %>% select(panel, outcome_clean, exposure, geo_cat),
            by = c("panel", "outcome_clean", "exposure")) %>%
  filter(is.finite(med_fold)) %>% mutate(geo_cat = factor(geo_cat, levels = names(geo_pal)))
stopifnot(nrow(het_div) <= nrow(het_tbl))   # guard against the join above re-multiplying
sp <- suppressWarnings(cor.test(log(het_div$med_fold), het_div$I2, method = "spearman"))

## ---------------------------------------------------------------------------
## Presentation fixes (2026-08):
##  1. The fitted line is now drawn ONLY when the correlation it depicts is
##     significant. The panel reports Spearman rho = -0.09, p = 0.66 -- a null
##     -- yet unconditionally drew an OLS fit with a confidence ribbon, which
##     reads as a trend the data do not support.
##  2. y is clamped to [0, 100]. I-squared is a percentage bounded at zero, but
##     the ribbon extended to roughly -12%, outside the statistic's range.
##  3. Repel tuned (labels collided in both bottom corners) and seeded so the
##     layout is reproducible across knits.
##  4. Subtitle now states the k = 2 caveat: I-squared is badly behaved with
##     two cohorts, and much of the mass on the zero line is k = 2 pairs.
## ---------------------------------------------------------------------------
div_trend_is_real <- is.finite(sp$p.value) && sp$p.value < 0.05
n_k2 <- sum(het_div$k == 2, na.rm = TRUE)

p_div <- ggplot(het_div, aes(x = med_fold, y = I2)) +
  (if (div_trend_is_real)
     geom_smooth(method = "lm", se = TRUE, color = "grey30", fill = "grey85",
                 linewidth = 0.6, alpha = 0.4)
   else NULL) +
  geom_point(aes(size = k, fill = geo_cat), shape = 21, color = "grey20", stroke = 0.7, alpha = 0.9) +
  geom_text_repel(aes(label = exposure), size = 3.6, color = "grey25",
                  max.overlaps = 20, box.padding = 0.45, point.padding = 0.2,
                  min.segment.length = 0.2, segment.color = "grey60", seed = 42) +
  scale_x_log10(labels = function(x) paste0(round(x), "x")) +
  scale_fill_manual(values = geo_pal, name = "Cohort geography",
                    guide = guide_legend(override.aes = list(size = 3.2))) +
  scale_size_continuous(range = c(2.5, 6), breaks = c(2, 3, 4), name = "N cohorts") +
  coord_cartesian(ylim = c(0, 100)) +
  labs(title = "Heterogeneity is not explained by concentration divergence",
       x = "Between-cohort median concentration fold-change",
       y = expression("I"^2*" (between-cohort effect heterogeneity, %)"),
       # CAVEAT, no longer stated on the figure (caption removed 2026-08-10 --
       # axis labels, ticks and titles only). It still applies and must travel
       # with the panel in the manuscript caption: the x axis is not fully
       # comparable across cohorts. 4 of the 13 (DCHS, SWSD, THS, ZIP) are
       # creatinine-corrected in HHEAR_Data_Clean_Targeted_V3.Rmd (concentration
       # divided by urinary creatinine, typically 50-200 mg/dL) and 9 are not, so
       # any analyte shared between the two groups carries a ~100x offset that is
       # a harmonisation difference rather than an exposure difference. Measured:
       # all 18 of THS's urinary trace elements sit below their peers by a median
       # factor of 179. No pair is dropped on this account.
       ) +
  ms_theme

## ============================================================================
## 2. Panel B: within vs naive-pooled decomposition (maternal, UNADJUSTED)
##    Read from the validated cohort_regime_diagnostic.R output. UNADJUSTED is
##    deliberate: it isolates the cohort-pooling mechanism (Study is the
##    confounder under study), and covariate completeness is itself
##    cohort-specific so adjustment would drop cohorts.
##    Pair set = exactly panel A's pairs (maternal, >=2 cohorts) so B matches A.
## ============================================================================
A_pairs <- results %>% filter(scale == "log1p") %>%
  group_by(panel, outcome, exposure) %>% filter(n_distinct(study) >= 2) %>%
  ungroup() %>% distinct(panel, outcome, exposure)

decomp <- readr::read_csv(file.path(diag_dir, "cohort_regime_pairs.csv"),
                          show_col_types = FALSE) %>%
  semi_join(A_pairs, by = c("panel", "outcome", "analyte" = "exposure")) %>%
  mutate(
    class = factor(recode(class,
                          concordant = "Concordant",
                          magnitude_distort = "Magnitude distortion",
                          sign_reversal = "Sign reversal (Simpson)"),
                   levels = c("Concordant", "Magnitude distortion", "Sign reversal (Simpson)")),
    pair = paste0(analyte, " | ", pretty_out(outcome))
  )

## ---- Per-panel row order: each panel sorted by its OWN metric ---------------
## A (spread): by cohort-effect disagreement (z-range).  B (decomposition): by
## within-cohort effect. The two panels are therefore in different row orders.
common_pairs <- intersect(unique(as.character(strip_df$pair)),
                          unique(as.character(decomp$pair)))
strip_df <- strip_df %>% filter(pair %in% common_pairs)
decomp   <- decomp   %>% filter(pair %in% common_pairs)

## ROW ORDER REWORKED 2026-08-10. Previously each panel was sorted purely by its
## own metric, which interleaved the three outcomes -- BMI, Birthweight and
## Gestational Age rows alternated all the way down both panels, so a reader could
## not see which outcome a block of rows belonged to without reading every label.
## Now rows are BLOCKED BY OUTCOME first and sorted by the panel's own metric
## WITHIN each block, which keeps the intended per-panel sort while making the
## outcome structure legible. Blocks are ordered by number of pairs (largest
## first, ties alphabetical) so the block order is data-driven, not hand-set.
pair_outcome <- function(p) sub("^.* \\| ", "", as.character(p))

blocks <- strip_df %>%
  distinct(pair) %>% mutate(o = pair_outcome(pair)) %>%
  count(o, name = "n_pairs") %>% arrange(desc(n_pairs), o) %>% pull(o)

# ggplot puts factor level 1 at the BOTTOM of a discrete y axis, so levels are
# built bottom-up: desc(blk) places blocks[1] at the top, and ascending metric
# within a block places the largest value at the top of that block.
order_pairs <- function(df, metric) {
  df %>% distinct(pair, {{ metric }}) %>%
    mutate(blk = match(pair_outcome(pair), blocks)) %>%
    arrange(desc(blk), {{ metric }}) %>% pull(pair)
}

lvlA <- order_pairs(strip_df, spread)
strip_df <- strip_df %>% mutate(pair = factor(pair, levels = lvlA))

lvlB <- order_pairs(decomp, r_within)
decomp <- decomp %>% mutate(pair = factor(pair, levels = lvlB))

# Boundaries between outcome blocks, for the separator rules drawn on A and B.
block_breaks <- function(lvls) {
  o <- pair_outcome(lvls)
  which(head(o, -1) != tail(o, -1)) + 0.5
}
breaksA <- block_breaks(lvlA)
breaksB <- block_breaks(lvlB)
stopifnot(length(breaksA) == length(blocks) - 1,
          length(breaksB) == length(blocks) - 1)

## ---- Panel A plot (polished: row banding, bigger points, clearer sig cue) ---
nA <- nlevels(strip_df$pair)
bandA <- data.frame(ypos = which(seq_len(nA) %% 2 == 0))

p_spread <- ggplot(strip_df, aes(x = z, y = pair)) +
  geom_rect(data = bandA, inherit.aes = FALSE,
            aes(ymin = ypos - 0.5, ymax = ypos + 0.5, xmin = -Inf, xmax = Inf),
            fill = "grey95") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey55", linewidth = 0.45) +
  # Rules separating the outcome blocks (see the row-order note above).
  geom_hline(yintercept = breaksA, color = "grey55", linewidth = 0.5) +
  geom_line(aes(group = pair), color = "grey75", linewidth = 0.5) +
  geom_point(aes(fill = location, color = sig_lab, alpha = sig_lab),
             shape = 21, size = 3.8, stroke = 0.9) +
  scale_fill_manual(values = loc_pal, name = "Cohort location",
                    guide = guide_legend(order = 1, nrow = 3,
                                         override.aes = list(color = "grey35", size = 3.4, alpha = 1))) +
  scale_color_manual(values = c("p < 0.05" = "grey5", "p >= 0.05" = "grey70"), name = NULL,
                     guide = guide_legend(order = 2, override.aes = list(fill = "grey85", size = 3.4, alpha = 1))) +
  scale_alpha_manual(values = c("p < 0.05" = 1, "p >= 0.05" = 0.45), guide = "none") +
  labs(title = "Cohort effects disagree in direction across locations",
       x = "Cohort effect estimate (z-statistic)", y = NULL) +
  ms_theme + theme(
                   panel.grid.major.y = element_blank(),
                   legend.box = "vertical")

## ---- Panel B plot (matching alternating row shading) -----------------------
class_pal <- c("Concordant"="#0F4C5C","Magnitude distortion"="#D17A22","Sign reversal (Simpson)"="#C1121F")
nB <- nlevels(decomp$pair)
bandB <- data.frame(ypos = which(seq_len(nB) %% 2 == 0))

p_decomp <- ggplot(decomp, aes(y = pair)) +
  geom_rect(data = bandB, inherit.aes = FALSE,
            aes(ymin = ypos - 0.5, ymax = ypos + 0.5, xmin = -Inf, xmax = Inf),
            fill = "grey95") +
  geom_vline(xintercept = 0, color = "grey55", linewidth = 0.45) +
  geom_hline(yintercept = breaksB, color = "grey55", linewidth = 0.5) +
  geom_segment(aes(x = r_within, xend = r_naive, yend = pair, color = class), linewidth = 1.0) +
  geom_point(aes(x = r_within), color = "#14213D", size = 2.8) +
  geom_point(aes(x = r_naive), shape = 21, fill = "white", color = "#14213D", size = 3.0, stroke = 0.9) +
  scale_color_manual(values = class_pal, name = NULL, drop = FALSE) +
  labs(title = "Pooling distorts or reverses the within-cohort effect",
       x = "Standardized association (partial-correlation scale)", y = NULL,
       # CAPTION REWORDED 2026-08-10. It previously said pooling pulls the estimate
       # toward the between-cohort (ecological) "relationship". That overstates
       # what the term can carry here: r_between is a correlation across K cohort
       # means, and at K = 2 -- which is 37 of the 44 pairs -- a correlation
       # between two points is +/-1 by construction (verified: all 37 are exactly
       # 1). So for most pairs it is the SIGN of the difference between two cohort
       # means, not a fitted gradient. Among the 5 pairs with K >= 4, where it is
       # estimable, |r_between| is at most 0.42. The reversal classification is
       # unaffected (it compares within vs naive pooled and never uses r_between),
       # so only the mechanistic wording needed fixing.
       # CAVEAT, no longer on the figure (caption removed 2026-08-10): pooling
       # pulls the estimate toward the between-cohort contrast, and where it
       # crosses zero (red) that is Simpson's paradox. For two-cohort pairs that
       # contrast is a direction only, not a gradient. Carry this in the
       # manuscript caption.
       ) +
  ms_theme + theme(panel.grid.major.y = element_blank())

## ============================================================================
## 2b. "The cost" panel: budget allocation (diverse units vs unplanned volume)
##     Read from the simulation output; styled with ms_theme to match A-C.
## ============================================================================
budget <- readr::read_csv(file.path(out_dir, "Pooling_Sim_Output",
                                    "budget_allocation_results.csv"),
                          show_col_types = FALSE)
budget_metric_levels <- c(
  "Population-average effect (RMSE)"  = "beta_within_rmse",
  "Cross-cohort heterogeneity (RMSE)" = "tau_rmse",
  "Single location (median SE)"       = "median_location_se")
budget_palette <- c(
  "Population-average effect (RMSE)"  = "#0F4C5C",
  "Cross-cohort heterogeneity (RMSE)" = "#9A4C95",
  "Single location (median SE)"       = "#D17A22")
budget_long <- budget %>%
  select(K, n_per, all_of(unname(budget_metric_levels))) %>%
  tidyr::pivot_longer(-c(K, n_per), names_to = "metric_key", values_to = "value") %>%
  mutate(metric = factor(names(budget_metric_levels)[match(metric_key, budget_metric_levels)],
                         levels = names(budget_metric_levels)))
## ---------------------------------------------------------------------------
## Presentation fixes (2026-08):
##  1. Break thinning. Every simulated K was labelled, so on a log x axis the
##     adjacent values overprinted -- "5" onto "6", "(1200)" onto "(1000)".
##     Keep a break only if it is far enough from the previous kept one in log
##     space; the endpoints are always kept so the axis still shows its range.
##  2. The crossover is the panel's actual finding -- below it, diverse units
##     win; above it, subjects-per-cohort wins -- so it is now marked rather
##     than left for the reader to locate.
## ---------------------------------------------------------------------------
.k_all  <- sort(unique(budget$K))
.keep   <- rep(FALSE, length(.k_all))
.last   <- -Inf
for (i in seq_along(.k_all)) {
  # 0.11 log10 units. K = 5 and 6 sit 0.079 apart, so a smaller gap still let
  # "5"/"6" and "(1200)"/"(1000)" overprint; 0.11 drops the 6 and leaves every
  # other break (closest remaining pair, 15 vs 20, is 0.125 apart).
  if (log10(.k_all[i]) - .last >= 0.11) { .keep[i] <- TRUE; .last <- log10(.k_all[i]) }
}
.keep[length(.k_all)] <- TRUE          # always show the largest K
k_breaks <- .k_all[.keep]
k_labels <- budget %>% distinct(K, n_per) %>% arrange(K) %>%
  dplyr::filter(K %in% k_breaks) %>%
  mutate(lab = paste0(K, "\n(", n_per, ")")) %>% pull(lab)

## Locate where "many cohorts" stops paying and "more subjects each" starts.
.cross_k <- local({
  w <- budget_long %>%
    dplyr::filter(metric %in% c("Population-average effect (RMSE)",
                                "Single location (median SE)")) %>%
    dplyr::select(K, metric, value) %>%
    tidyr::pivot_wider(names_from = metric, values_from = value) %>%
    dplyr::arrange(K)
  d <- w$`Population-average effect (RMSE)` - w$`Single location (median SE)`
  i <- which(diff(sign(d)) != 0)[1]
  if (is.na(i)) return(NA_real_)
  # linear interpolation in log10(K), where the curves are near-linear
  l1 <- log10(w$K[i]); l2 <- log10(w$K[i + 1])
  10^(l1 + (l2 - l1) * abs(d[i]) / (abs(d[i]) + abs(d[i + 1])))
})

p_cost <- ggplot(budget_long, aes(x = K, y = value, color = metric)) +
  (if (is.finite(.cross_k))
     list(
       geom_vline(xintercept = .cross_k, linetype = "dotted",
                  color = "grey45", linewidth = 0.5),
       annotate("text", x = .cross_k, y = Inf,
                label = sprintf("  crossover ~%d cohorts", round(.cross_k)),
                hjust = 0, vjust = 1.4, size = 3.8, color = "grey35"))
   else NULL) +
  geom_line(linewidth = 1.0) + geom_point(size = 2.0) +
  scale_x_log10(breaks = k_breaks, labels = k_labels) +
  scale_color_manual(values = budget_palette, name = NULL) +
  labs(x = "Number of cohorts  (subjects per cohort)",
       y = "Error / uncertainty (lower is better)") +
  ms_theme + theme(
                   axis.text.x = element_text(size = 10))

## ============================================================================
## 2c. Panel D (new, 2026-08-10): between-cohort share of exposure variance vs
##     how far pooling moves the estimate.
##
## This replaces the budget-allocation simulation as panel D (that panel moves to
## supplementary; it is still built above and still saved as `cost`).
##
## Rationale. Panels A-C show that cohorts disagree and that pooling reverses
## signs, but the mechanism offered for it -- r_between, the ecological slope --
## is a correlation across K cohort means, and at K = 2 (37 of 44 pairs) a
## correlation between two points is +/-1 by construction. It carries a sign and
## no magnitude. The ICC of log1p(concentration) across the contributing cohorts
## measures the same idea (how much of the exposure's variation is BETWEEN rather
## than WITHIN cohorts) but is well behaved with two groups, so it is available
## for 43 of the 44 pairs.
##
## The relationship is empirical and survives the creatinine confound:
##   all pairs            n = 43  rho = 0.59  p = 5.0e-05
##   same-status pairs    n = 34  rho = 0.61  p = 2.0e-04
##   all-uncorrected only n = 24  rho = 0.38  p = 0.069  (same direction, small n)
## Points are coloured BY CORRECTION STRATUM on purpose, so the harmonisation
## issue is visible in the figure rather than only in a caption.
##
## Deliberately plotted against the pooling SHIFT (continuous), not the reversal
## class: ICC predicts the shift robustly, but its association with the binary
## reversal label loses significance once mixing is controlled (Wilcoxon p = 0.12).
## ============================================================================
CRE_CORRECTED <- c("DCHS", "SWSD", "THS", "ZIP")

# One exposure value per participant per analyte, matching the pipeline's
# group_by(Study, PID, Analyte_Code) %>% mean() convention.
mat_exposure <- Targeted_Data_All %>%
  transmute(study = std_study(Study),
            Mat_PID = as.character(Mat_PID),
            exposure = as.character(Analyte_Code),
            conc = suppressWarnings(as.numeric(Concentration))) %>%
  filter(!is.na(Mat_PID), Mat_PID != "", !is.na(conc), is.finite(conc), conc >= 0) %>%
  mutate(x = log1p(conc)) %>%
  group_by(study, Mat_PID, exposure) %>%
  summarise(x = mean(x), .groups = "drop")

# Between-cohort share of variance = SSB / (SSB + SSW), on the cohorts that
# actually contribute to that pair.
icc_for <- function(an, sts) {
  d <- mat_exposure %>% filter(exposure == an, study %in% sts)
  if (dplyr::n_distinct(d$study) < 2 || nrow(d) < 10) return(NA_real_)
  g <- d %>% group_by(study) %>% summarise(n = dplyr::n(), m = mean(x), .groups = "drop")
  grand <- sum(g$n * g$m) / sum(g$n)
  ssb <- sum(g$n * (g$m - grand)^2)
  ssw <- d %>% group_by(study) %>% summarise(s = sum((x - mean(x))^2), .groups = "drop") %>%
    pull(s) %>% sum()
  if (!is.finite(ssb + ssw) || (ssb + ssw) == 0) return(NA_real_)
  ssb / (ssb + ssw)
}

## Colour key = BIOLOGICAL MATRIX CONCORDANCE (2026-08-10). Previously the points
## were coloured by creatinine-correction status. Matrix concordance is the more
## fundamental comparability question: whether the cohorts pooled for a pair even
## measured the analyte in the same medium. 16 of the 44 pairs do not -- 9 pool
## blood with urine and 7 pool blood with "Other" -- and pooling a blood
## concentration with a urine concentration of the same analyte is not a scale
## offset that correction can fix.
##
## DEPENDENCY: Matrix_Class was added to Harmonized_Targeted.RData on 2026-08-10 by
## the sample-matrix harmonisation work, derived from Sample_Matrix. This panel
## therefore depends on that column; if its definition changes, re-run this script.
## Only maternal rows with a non-missing Matrix_Class are used.
matrix_lookup <- Targeted_Data_All %>%
  transmute(study = std_study(Study),
            Mat_PID = as.character(Mat_PID),
            exposure = as.character(Analyte_Code),
            Matrix_Class = as.character(Matrix_Class)) %>%
  filter(!is.na(Mat_PID), Mat_PID != "", !is.na(Matrix_Class)) %>%
  distinct(study, exposure, Matrix_Class)

matrices_for <- function(an, sts) {
  v <- matrix_lookup %>% filter(exposure == an, study %in% sts) %>% pull(Matrix_Class) %>% unique()
  if (!length(v)) return(NA_character_)
  paste(sort(v), collapse = " + ")
}

pair_cohorts <- results %>% filter(scale == "log1p") %>%
  group_by(panel, outcome, exposure) %>%
  filter(dplyr::n_distinct(study) >= 2) %>%
  summarise(k = dplyr::n_distinct(study), studies = list(unique(study)),
            n_corr = dplyr::n_distinct(study[study %in% CRE_CORRECTED]),
            n_raw  = dplyr::n_distinct(study[!study %in% CRE_CORRECTED]),
            .groups = "drop") %>%
  mutate(cre_stratum = case_when(n_corr > 0 & n_raw == 0 ~ "All creatinine-corrected",
                                 n_raw  > 0 & n_corr == 0 ~ "All uncorrected",
                                 TRUE ~ "Mixed correction status"),
         icc = mapply(icc_for, exposure, studies),
         matrices = mapply(matrices_for, exposure, studies),
         stratum = case_when(is.na(matrices) ~ NA_character_,
                             grepl("[+]", matrices) ~ "Different matrices",
                             TRUE ~ "Same matrix"))

var_tbl <- decomp %>%
  transmute(outcome, exposure = analyte, class, r_within, r_naive,
            shift = abs(r_naive - r_within)) %>%
  inner_join(pair_cohorts %>% select(outcome, exposure, k, icc, stratum,
                                     matrices, cre_stratum),
             by = c("outcome", "exposure")) %>%
  filter(is.finite(icc), is.finite(shift), !is.na(stratum))

strat_pal <- c("Same matrix" = "#2171B5", "Different matrices" = "#D7301F")
var_tbl <- var_tbl %>% mutate(stratum = factor(stratum, levels = names(strat_pal)))

sp_var      <- suppressWarnings(cor.test(var_tbl$icc, var_tbl$shift, method = "spearman"))
sp_var_same <- suppressWarnings(with(dplyr::filter(var_tbl, stratum == "Same matrix"),
                                     cor.test(icc, shift, method = "spearman")))
var_trend_is_real <- is.finite(sp_var$p.value) && sp_var$p.value < 0.05

## LABEL THINNING (2026-08-10): labelling all 43 points made the low-share/
## low-distortion corner unreadable. Only the high-distortion pairs carry a label
## now -- those are the ones the panel is about -- and the cut is data-driven.
SHIFT_LABEL_CUT <- 0.10
var_lab <- var_tbl %>% filter(shift > SHIFT_LABEL_CUT)

p_varshare <- ggplot(var_tbl, aes(x = icc, y = shift)) +
  (if (var_trend_is_real)
     geom_smooth(method = "lm", se = TRUE, color = "grey30", fill = "grey85",
                 linewidth = 0.6, alpha = 0.4)
   else NULL) +
  geom_point(aes(size = k, fill = stratum), shape = 21, color = "grey20",
             stroke = 0.7, alpha = 0.9) +
  ggrepel::geom_text_repel(data = var_lab, aes(label = exposure), size = 3.9,
                           color = "grey25", max.overlaps = Inf, box.padding = 0.5,
                           point.padding = 0.25, min.segment.length = 0.2,
                           segment.color = "grey60", seed = 42) +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1)) +
  scale_fill_manual(values = strat_pal, name = "Biological matrix",
                    guide = guide_legend(override.aes = list(size = 3.2))) +
  scale_size_continuous(range = c(2.5, 6), breaks = c(2, 3, 4), name = "N cohorts") +
  # Kept to ~55 characters: this panel is half the figure width and at 17pt a
  # longer title clips at the right edge rather than wrapping.
  labs(title = "Distortion scales with between-cohort exposure variance",
       x = "Between-cohort share of exposure variance",
       y = "Distortion from pooling  |naive - within|",
       # CAVEAT, no longer on the figure (caption removed 2026-08-10). Between-
       # cohort share is the ICC of log1p(concentration) across the cohorts
       # contributing to each pair -- unlike the ecological slope in B it is
       # estimable at two cohorts. Colour marks whether those cohorts measured the
       # analyte in the same biological medium; the pairs that do not (see
       # panelD_variance_share.csv, stratum column) carry ~4x the median
       # distortion of same-matrix pairs. Carry this in the manuscript caption --
       # it is the reason the overall correlation is stronger than the
       # same-matrix one.
       ) +
  ms_theme

## ============================================================================
## 3. Assemble + save
## ============================================================================
# Save the individual panels for the manuscript assembler (All_Patchplots.Rmd).
# `varshare` is the new panel D; `cost` is the budget simulation, retained for
# the supplementary figure.
saveRDS(list(spread = p_spread, decomp = p_decomp, div = p_div,
             varshare = p_varshare, cost = p_cost),
        file.path(diag_dir, "figure3_panels.rds"))

# Standalone supplementary PDF for the retired budget panel.
ggsave(file.path(final_dir, "Supp_Pooling_Budget_Tradeoff.pdf"), p_cost,
       height = 6, width = 9)

readr::write_csv(var_tbl %>% select(outcome, exposure, k, stratum, matrices, cre_stratum, icc, r_within, r_naive, shift, class),
                 file.path(diag_dir, "panelD_variance_share.csv"))

combined <- (p_spread / p_decomp / p_div) +
  plot_layout(heights = c(1.35, 1.35, 1.0)) +
  plot_annotation(
    tag_levels = "A",
    title = "Cross-cohort heterogeneity and the consequences of pooling in the targeted ExWAS (maternal)",
    theme = theme(plot.title = element_text(face = "bold", size = 18, color = "#14213D"))) &
  theme(plot.tag = element_text(size = 16, face = "bold"))

dir.create(final_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(diag_dir, recursive = TRUE, showWarnings = FALSE)
ggsave(file.path(final_dir, "MS_Fig_RealData_Combined.pdf"), combined,
       width = 9, height = 17, dpi = 300)
saveRDS(combined, file.path(diag_dir, "combined_real_fig.rds"))

cat("decomposition (maternal pairs matching panel A):\n")
print(as.data.frame(decomp %>% mutate(across(where(is.numeric), ~round(., 3)))))
cat("\nSaved combined figure to:\n  ", normalizePath(final_dir), "\n", sep = "")
