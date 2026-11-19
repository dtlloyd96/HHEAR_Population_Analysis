# R4_5_supp_figures.R
# Reviewer 4 -- two paneled supplementary figures. Reads saved R4_5 outputs only
# (plus the raw MEBCO table for the arsenic panel), so run it after the other
# R4_5 scripts.
#   Fig_R4_S1_Reported_Sources  sources reported by participants, within cohorts:
#     A key within-cohort results (one row per source and main biomarker)
#     B same source pooled across cohorts (random effects)
#     C MEBCO household drinking-water arsenic vs urinary arsenic
#   Fig_R4_S2_Air_Water_Linkages  locations linked to regional air and water data:
#     A-D maps of the four layers with location biomarker medians (panels from R4_5_map.R)
#     I-K air: 1-hydroxypyrene vs ambient PM2.5 and vs household polluting cooking fuel;
#             PM2.5 effect before/after adding fuel (individual-level models)
#     E-H pesticides, diet and leaded petrol: 3-PBA and TCPy vs insecticide use, As vs fish supply, Pb vs years since ban
#     L-M water: urinary As vs groundwater As; blood PFOS/PFOA vs US drinking-water
#             PFAS detections (UCMR3)
# Writes PNG + PDF (cairo) to Output/Reviewer_Round2/R4_5_determinants/supp_figures/.
# Run from project root.

suppressMessages({library(dplyr); library(tidyr)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
det <- "Output/Reviewer_Round2/R4_5_determinants"
out <- file.path(det, "supp_figures"); dir.create(out, showWarnings = FALSE)
rd  <- function(...) read.csv(file.path(det, ...), check.names = FALSE)
lab_loc <- function(x) gsub("_", " ", x)

source("Coding/Paper1_Submission2/Reviewer_Round2/R4_5_fig_style.R")
save_fig <- function(g, name, w, h) {
  ggsave(file.path(out, paste0(name, ".png")), g, width = w, height = h, dpi = 300)
  ggsave(file.path(out, paste0(name, ".pdf")), g, width = w, height = h, device = cairo_pdf)
}
FDR <- rd("register", "location_level_tests_fdr.csv")
rho_lab <- function(a, layer, st = "all") { z <- FDR[FDR$analyte == a & FDR$layer == layer & FDR$set == st, ]
  sprintf("ρ = %.2f, p = %.2g, q = %.2g (n = %d)", z$rho, z$perm_p, z$q_all, as.integer(z$n_locations)) }
pal_verdict <- c("FDR < 0.05" = "#c0392b", "p < 0.05" = "#e08e0b", "No association" = "grey55", "Opposite direction" = "#3b6fb6")

# ======================= Figure S1: reported sources ==========================
# A. One row per source and its main biomarker. Cotinine (the positive control,
# 7-97x) is left out so the axis stays readable; it is in the text and Appendix.
K <- rd("within_cohort_sources", "within_cohort_key_results.csv")
keep <- tribble(
  ~cohort, ~panel,   ~source_label,                          ~analyte, ~group,
  "DCHS",  "Mother", "Smoking in pregnancy",                 "NAP1",   "Tobacco",
  "DCHS",  "Mother", "Smoking in pregnancy",                 "PYR1",   "Tobacco",
  "SICAS", "Child",  "Smoker(s) at home",                    "NNAL",   "Tobacco",
  "DCHS",  "Mother", "Fossil-fuel cooking stove",            "NAP1",   "Combustion",
  "DCHS",  "Mother", "Fossil-fuel heating",                  "NAP1",   "Combustion",
  "DCHS",  "Mother", "Home PM10 (modelled)",                 "PHEN1",  "Combustion",
  "MADRES","Mother", "Any seafood",                          "As",     "Diet",
  "MADRES","Mother", "Any seafood",                          "Hg",     "Diet",
  "NFGS",  "Child",  "No effort to use BPA-free containers", "MBZP",   "Packaging",
  "NFGS",  "Child",  "Plastic-packaged foods (24 h)",        "MEHHP",  "Packaging",
  "SICAS", "Child",  "Hair spray",                           "MEP",    "Personal care",
  "SICAS", "Child",  "Hair oil",                             "PRPB",   "Personal care",
  "ZIP",   "Mother", "Mosquito repellent use",               "PBA",    "Pests",
  "EEICC", "Child",  "Cockroaches seen",                     "PBA",    "Pests",
  "EEICC", "Child",  "Hardwood floors (vs carpet/vinyl)",    "MBZP",   "Housing")
an_lab <- c(NAP1 = "1-naphthol", PYR1 = "1-OH-pyrene", NNAL = "NNAL", PHEN1 = "1-OH-phenanthrene", As = "arsenic", Hg = "mercury",
            MBZP = "MBzP", MEHHP = "MEHHP", MEP = "MEP", PRPB = "propylparaben", PBA = "3-PBA")
s1a <- keep %>% inner_join(K, by = c("cohort", "panel", "source_label", "analyte")) %>%
  mutate(verdict = case_when(!direction_ok & p < 0.05 ~ "Opposite direction", q < 0.05 ~ "FDR < 0.05", p < 0.05 ~ "p < 0.05", TRUE ~ "No association"),
         verdict = factor(verdict, levels = names(pal_verdict)),
         adj = ifelse(model == "adjusted for smoking", ", smoking-adj.", ""),
         label = sprintf("%s → %s (%s%s)", sub(" \\(vs carpet/vinyl\\)", "", source_label), an_lab[analyte], cohort, adj),
         group = factor(group, levels = unique(keep$group)))
stopifnot(nrow(s1a) == nrow(keep))
s1a <- s1a %>% mutate(label = factor(label, levels = rev(label)))
pS1A <- ggplot(s1a, aes(fold, label, colour = verdict)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_linerange(aes(xmin = fold_lo, xmax = fold_hi), linewidth = 0.6) + geom_point(size = 1.8) +
  facet_grid(group ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_x_log10(breaks = c(0.5, 1, 2, 4), labels = function(x) format(x, drop0trailing = TRUE)) + scale_colour_manual(values = pal_verdict, drop = FALSE) +
  labs(title = "Within-cohort sources: key results", subtitle = "One row per source and main biomarker; cotinine positive control not shown",
       x = "Fold difference, exposed vs not or per doubling (95% CI)", y = NULL, colour = NULL) +
  theme(legend.position = "bottom", strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1, size = 8.5)) + guides(colour = guide_legend(nrow = 2))

# B. Pooled across cohorts (cotinine excluded for scale).
MA <- rd("source_meta_analysis", "pooled_effects.csv")
s1b <- MA %>% filter(outcome != "Cotinine") %>%
  mutate(label = sprintf("%s\n%d cohorts, I² = %.0f%%", sub("NAP1", "1-naphthol", sub("NAP2", "2-naphthol", sub("PYR1", "1-OH-pyrene", pair))), as.integer(k), I2),
         sig = ifelse(lo > 1, "95% CI excludes 1", "95% CI includes 1")) %>% arrange(pooled_fold) %>% mutate(label = factor(label, levels = label))
pS1B <- ggplot(s1b, aes(pooled_fold, label, colour = sig)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.6) + geom_point(shape = 18, size = 3) + scale_x_log10() +
  scale_colour_manual(values = c("95% CI excludes 1" = "black", "95% CI includes 1" = "grey55")) +
  labs(title = "Same source, pooled across cohorts", subtitle = sprintf("Random effects; smoking → cotinine (not shown) %.0f×", MA$pooled_fold[MA$outcome == "Cotinine"]),
       x = "Pooled fold difference (95% CI)", y = NULL, colour = NULL) + theme(legend.position = "bottom", axis.text.y = element_text(size = 7))

# C. MEBCO household water arsenic.
m <- new.env(); load("Input/MEBCO_34_Data_All.RData", envir = m); ML <- m[[ls(m)[1]]]
mb <- ML$MEBCO_Targeted %>% transmute(pid = as.character(PID), uAs = As) %>%
  inner_join(ML$MEBCO_Epi %>% transmute(pid = as.character(PID), wAs = wateras_c1) %>% distinct(pid, .keep_all = TRUE), by = "pid") %>%
  filter(!is.na(wAs), !is.na(uAs), uAs > 0) %>% distinct(pid, .keep_all = TRUE) %>% mutate(wAs_plot = pmax(wAs, 0.5))
W1a <- rd("water", "1a_MEBCO_water_vs_urine_As.csv")[1, ]
pS1C <- ggplot(mb, aes(wAs_plot, uAs)) + annotate("rect", xmin = 10, xmax = Inf, ymin = 0, ymax = Inf, fill = pal$red, alpha = 0.04) +
  geom_vline(xintercept = c(10, 50), linetype = 3, colour = pal$grey) +
  geom_point(alpha = 0.6, size = 1.6, colour = pal$blue) + geom_smooth(method = "lm", formula = y ~ x, colour = pal$dark, fill = pal$light, linewidth = 0.7) +
  annotate("text", x = c(10, 50), y = min(mb$uAs), label = c("WHO 10 µg/L", "Bangladesh 50 µg/L"), angle = 90, hjust = 0, vjust = -0.5, size = 2.8, colour = pal$mid) +
  scale_x_log10() + scale_y_log10() +
  labs(title = "Household drinking water and urinary arsenic",
       subtitle = sprintf("MEBCO, Bangladesh: %.1f-fold (95%% CI %.1f–%.1f) higher urinary\narsenic per 10-fold more water arsenic; ρ = %.2f, n = %d",
                          W1a$fold_per_10x_water, W1a$fold_lo, W1a$fold_hi, W1a$spearman, as.integer(W1a$n)),
       x = "Household water arsenic (µg/L; < 0.5 plotted at 0.5)", y = "Urinary arsenic, unadjusted (µg/L)")

g1 <- (pS1A | (pS1B / pS1C)) + plot_layout(widths = c(1.15, 1)) + plot_annotation(tag_levels = "A")
save_fig(g1, "Fig_R4_S1_Reported_Sources", 12, 8.5)

# ======================= Figure S2: air and water linkages =====================
M <- rd("geospatial", "biomarker_location_medians.csv"); L <- rd("geospatial", "location_layers_time_matched.csv")
pyr <- M %>% filter(Analyte_Code == "PYR1") %>% inner_join(L %>% select(loc_id, PM25_ugm3, polluting_fuel_pct), by = "loc_id") %>%
  mutate(set = ifelse(set == "all", "All locations", "Dilution-adjusted urine"))
pyr_pts <- pyr %>% group_by(loc_id) %>% slice_max(set == "Dilution-adjusted urine", n = 1, with_ties = FALSE) %>% ungroup() %>%
  mutate(urine = ifelse(set == "Dilution-adjusted urine", "Dilution-adjusted", "Unadjusted"))
pyr_scatter <- function(xvar, layer, title, xlab, col) {
  ggplot(pyr_pts %>% filter(!is.na(.data[[xvar]])), aes(.data[[xvar]], median)) +
    geom_smooth(data = pyr %>% filter(set == "Dilution-adjusted urine", !is.na(.data[[xvar]])), method = "lm", formula = y ~ x,
                se = FALSE, colour = col, linewidth = 0.5) +
    geom_point(aes(shape = urine), size = 2.3, colour = col, fill = "white", stroke = 0.8) +
    geom_text_repel(aes(label = lab_loc(loc_id)), size = 2.3, max.overlaps = 30, segment.size = 0.2) +
    scale_shape_manual(values = c(`Dilution-adjusted` = 16, Unadjusted = 21)) + scale_y_log10() +
    annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.3, size = 2.3,
             label = paste0("Adjusted urine: ", rho_lab("PYR1", layer, "dilution_adjusted_only"), "\nAll locations: ", rho_lab("PYR1", layer, "all"))) +
    labs(title = title, x = xlab, y = "1-hydroxypyrene, location median (ng/mL)", shape = "Urine") + theme(legend.position = "bottom")
}
pS2A <- pyr_scatter("PM25_ugm3", "PM25_ugm3", "Ambient air: PM2.5", "PM2.5, population-weighted, sampling years (µg/m³)", "grey35")
pS2B <- pyr_scatter("polluting_fuel_pct", "polluting_fuel_pct", "Household air: polluting cooking fuel", "% of population cooking with polluting fuels (WHO)", "#c0392b")

DM <- rd("geospatial", "pah_pm25", "D_cleanfuel_individual_models.csv")
fd <- DM %>% filter(term == "lpm_loc", model %in% c("M1_PM25", "M3_both")) %>%
  mutate(model = factor(recode(model, M1_PM25 = "PM2.5 alone", M3_both = "PM2.5 + cooking fuel"), levels = c("PM2.5 alone", "PM2.5 + cooking fuel")),
         analyte = recode(analyte, NAP1 = "1-naphthol", NAP2 = "2-naphthol", PYR1 = "1-OH-pyrene"),
         subset = recode(subset, `All locations` = "All locations", `Dilution-adjusted urine only` = "Adjusted urine"))
pS2C <- ggplot(fd, aes(effect, analyte, colour = model)) + geom_vline(xintercept = 1, linetype = 2, colour = "grey50") +
  geom_linerange(aes(xmin = effect_lo, xmax = effect_hi), linewidth = 0.6, position = position_dodge(0.5)) +
  geom_point(size = 1.8, position = position_dodge(0.5)) + facet_wrap(~subset, ncol = 1) + scale_x_log10() +
  scale_colour_manual(values = c("PM2.5 alone" = "grey45", "PM2.5 + cooking fuel" = "#c0392b")) +
  labs(title = "Cooking fuel absorbs the PM2.5 association", subtitle = "Individual-level mixed models",
       x = "Fold change per doubling of PM2.5 (95% CI)", y = NULL, colour = NULL) + theme(legend.position = "bottom")

LM <- rd("diet_pesticide_lead", "location_medians_with_layers.csv")
as_loc <- LM %>% filter(Analyte_Code == "As", set == "all", !is.na(P_As))
pS2D <- ggplot(as_loc, aes(P_As, median)) + geom_point(size = 2.1, colour = "#3b6fb6") +
  geom_text_repel(aes(label = lab_loc(loc_id)), size = 2.3, segment.size = 0.2) + scale_y_log10(expand = expansion(mult = c(0.05, 0.2))) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.4, size = 2.3, label = paste0(rho_lab("As", "GW_As_prob_gt10"), "\nDriven by Bangladesh")) +
  labs(title = "Groundwater arsenic", x = "Probability groundwater As > 10 µg/L (Podgorski & Berg 2020)", y = "Urinary arsenic, location median")

pf <- rd("water", "2_PFAS_location_medians_UCMR3.csv") %>% filter(is.finite(year_max), year_max >= 2000) %>%
  mutate(pct_own = ifelse(Analyte_Code == "PFOS", pct_PFOS, pct_PFOA))
pS2E <- ggplot(pf, aes(pct_own, median, colour = Analyte_Code)) + geom_point(size = 2.1) +
  geom_text_repel(aes(label = lab_loc(loc_id)), size = 2.2, segment.size = 0.2, show.legend = FALSE) +
  scale_y_log10(expand = expansion(mult = c(0.05, 0.25))) + scale_colour_manual(values = c(PFOS = "#3b6fb6", PFOA = "#e08e0b")) +
  annotate("text", x = -Inf, y = Inf, hjust = -0.03, vjust = 1.3, size = 2.3,
           label = paste0("PFOS: ", rho_lab("PFOS", "UCMR3 detections (own analyte)"), "\nPFOA: ", rho_lab("PFOA", "UCMR3 detections (own analyte)"))) +
  labs(title = "US drinking-water PFAS (EPA UCMR3, 2013–15)", x = "% of public water systems with the analyte detected (state)",
       y = "Blood PFAS, location median (ng/mL)", colour = NULL) + theme(legend.position = "bottom")

# Pre-specified diet, pesticide and leaded-petrol linkages (location level).
us_locs <- c("Atlanta", "Sacramento", "Denver", "Baltimore", "SanFrancisco", "Boston", "Los_Angeles", "Syracuse", "New_York", "Michigan",
             "California", "Washington_State", "New_Hampshire", "Massachusetts", "NorthEast_US", "Midwestern_US", "United_States", "TEDDY_United_States")
rho_txt <- function(a, layer) { z <- FDR[FDR$analyte == a & FDR$layer == layer & FDR$set == "all", ]
  if (!nrow(z)) return(NULL)
  if (is.na(z$q_prespec)) return(sprintf("ρ = %.2f, p = %.2g, q (all tests) = %.2g, n = %d; exploratory", z$rho, z$perm_p, z$q_all, as.integer(z$n_locations)))
  sprintf("ρ = %.2f, p = %.2g, q (pre-specified) = %.2g, n = %d", z$rho, z$perm_p, z$q_prespec, as.integer(z$n_locations)) }
LMa <- LM %>% filter(set == "all")
dpl <- function(a, xv, title, xlab, ylab, logx = FALSE) {
  d <- LMa %>% filter(Analyte_Code == a, !is.na(.data[[xv]])) %>% mutate(reg = ifelse(loc_id %in% us_locs, "US", "Non-US"))
  p <- ggplot(d, aes(.data[[xv]], median)) + geom_smooth(method = "lm", formula = y ~ x, se = FALSE, colour = pal$grey, linewidth = 0.5) +
    geom_point(aes(colour = reg), size = 2.3) + geom_text_repel(aes(label = lab_loc(loc_id)), size = 2.4, segment.size = 0.2, max.overlaps = 20, seed = 3) +
    scale_y_log10() + scale_colour_manual(values = pal_region, name = NULL) +
    labs(title = title, subtitle = rho_txt(a, xv), x = xlab, y = ylab)
  if (logx) p + scale_x_log10() else p
}
pS2J <- dpl("PBA", "insecticide_kg_ha", "Pyrethroids: 3-PBA vs insecticide use", "Insecticide use (kg a.i./ha cropland)", "Urinary 3-PBA, location median", TRUE)
pS2K <- dpl("TCP", "insecticide_kg_ha", "Chlorpyrifos: TCPy vs insecticide use", "Insecticide use (kg a.i./ha cropland)", "Urinary TCPy, location median", TRUE)
pS2L <- dpl("As", "fish_kg_cap", "Diet: arsenic vs fish and seafood supply", "Fish and seafood supply (kg/person/year, FAOSTAT)", "Urinary arsenic, location median")
pS2M <- dpl("Pb", "years_since_lead_ban", "Leaded petrol: lead vs years since ban", "Years since leaded-petrol ban (OWID)", "Urinary lead, location median")
diet <- (pS2J | pS2K | pS2L | pS2M) + plot_layout(guides = "collect") & theme(legend.position = "bottom")

# Maps (A-D, built by R4_5_map.R) on top; location-level scatter and models (E-I) below.
mp <- readRDS(file.path(det, "map", "map_panels.rds"))
maps <- (mp$A | mp$B) / (mp$C | mp$D)
links <- (pS2A | pS2B | pS2C) / (pS2D | pS2E)
g2 <- wrap_plots(maps, diet, links, ncol = 1, heights = c(1.05, 0.42, 1)) +
  plot_annotation(tag_levels = "A",
                  caption = paste("A-D: HHEAR locations (unit interior point) over each layer; filled points = location median of the paired biomarker (log scale), grey = not measured.",
                                  "US multi-site cohorts sit at the US interior point.\nLayers: ACAG V6.GL.02.04 PM2.5; WHO/World Bank SDG 7.1.2; Podgorski & Berg 2020; FAOSTAT; EPA UCMR3; OWID leaded-petrol ban.\nE-H: location-level tests (G, arsenic vs fish, is exploratory); grey line = least-squares fit; q (pre-specified) = Benjamini-Hochberg across the pre-specified location-level tests; G and I-M show q across all location-level tests."),
                  theme = theme(plot.caption = element_text(size = 9, hjust = 0, colour = pal$mid))) &
  theme(plot.tag = element_text(face = "bold", size = 14))
save_fig(g2, "Fig_R4_S2_Air_Water_Linkages", 15, 25)
cat("written:", list.files(out), sep = "\n")
