# R4_5_PAH_PM25.R
# Reviewer 4, determinants of geography. Proof of concept: does ambient PM2.5
# explain between-location differences in urinary PAH metabolites?
# Builds on R4_5_geospatial_extract.R / _scatter.R (run those first). Run from
# the project root. Writes Output/Reviewer_Round2/R4_5_determinants/geospatial/pah_pm25/.
#
# Analyses
#   A. Exposure contrast: PM2.5 at the PAH locations, US vs non-US (fold range).
#   B. Location level (ecological): log10 location median ~ log10 PM2.5, OLS with
#      a location bootstrap CI; Spearman; leave-one-location-out; US-only and
#      non-US-only correlations.
#   C. Individual level, between/within decomposition (Mundlak):
#        log10(PAH) ~ PM_between + PM_within + dilution basis + (1|location) + (1|participant)
#      PM is each sample's location-year PM2.5 (log10). PM_between = location mean;
#      PM_within = deviation of the sampling year from it. The within slope uses
#      only year-to-year change inside a location (same cohort, lab, protocol),
#      so it is free of between-cohort confounding -- but PM2.5 varies little
#      over a cohort's sampling window, so it is expected to be imprecise.
#      Effects are reported as fold change in PAH per doubling of PM2.5 (2^beta).
#
# Data rules match R4_5_geospatial_scatter.R: urine; to_ng_per_mL(); <LOD -> LOD/sqrt(2);
# VDAART "PYR" remapped to PYR1 for this analysis only (canonical data unchanged).
# Cohorts with no sampling year (ELEMENT; ZIP/SWSD/ATLAA/MEBCO have no PAH) get
# the location mean PM2.5, i.e. contribute to the between term only.

suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(lme4); library(patchwork)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
source("Coding/Paper1_Submission2/harmonization_helpers.R")
geo <- "Output/Reviewer_Round2/R4_5_determinants/geospatial"
out <- file.path(geo, "pah_pm25"); dir.create(out, showWarnings = FALSE)
set.seed(20261104)

PAH <- c("NAP1", "NAP2", "PYR1")
us_locs <- c("Atlanta","Sacramento","Denver","Baltimore","SanFrancisco","Boston","Los_Angeles","Syracuse","New_York",
             "Michigan","California","Washington_State","New_Hampshire","Massachusetts","NorthEast_US","Midwestern_US",
             "United_States","TEDDY_United_States")
L  <- read.csv(file.path(geo, "location_layers_time_matched.csv"))
PM <- read.csv(file.path(geo, "layer_PM25_location_year.csv")) %>% select(loc_id, year, pm = value)
M  <- read.csv(file.path(geo, "biomarker_location_medians.csv")) %>% filter(Analyte_Code %in% PAH)

# ---- sample-level data -------------------------------------------------------
e <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = e)
s <- new.env(); load("Input/Harmonized_Datasets/Harmonized_SDOH.RData",   envir = s)
loc_map <- s$SDOH_Data_All %>% distinct(Study, Child_PID, Mat_PID, Location, Country) %>%
  mutate(Study = sub("_SDOH$", "_Targeted", Study),
         loc_id = case_when(Study == "TEDDY_Targeted" ~ paste0("TEDDY_", Location),
                            Study == "ZIP_Targeted" & !is.na(Country) ~ paste0("ZIP_", Country), TRUE ~ Location))
ok_units <- c("ng/mL", "ng/ml", "ng/L", "pg/mL", "pg/ml", "ug/dL", "mg/L")
d <- e$Targeted_Data_All %>%
  mutate(Analyte_Code = ifelse(Study == "VDAART_Targeted" & Analyte_Code == "PYR", "PYR1", Analyte_Code)) %>%
  filter(Analyte_Code %in% PAH, Matrix_Class == "Urine", is.na(Units) | Units %in% ok_units) %>%
  mutate(conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
         det = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
         val = ifelse(det, conc, lod / sqrt(2)), pid = paste(Study, coalesce(Mat_PID, Child_PID))) %>%
  filter(!is.na(val), val > 0)
d <- bind_rows(
  d %>% filter(Study != "TEDDY_Targeted") %>% left_join(loc_map %>% distinct(Study, loc_id), by = "Study"),
  d %>% filter(Study == "TEDDY_Targeted") %>% left_join(loc_map %>% filter(Study == "TEDDY_Targeted") %>% distinct(Child_PID, loc_id), by = "Child_PID")
) %>% filter(!is.na(loc_id)) %>%
  left_join(PM, by = c("loc_id", "Sample_Collection_Year" = "year")) %>%
  left_join(L %>% select(loc_id, PM25_loc = PM25_ugm3), by = "loc_id") %>%
  mutate(pm = coalesce(pm, PM25_loc), lpm = log10(pm), lval = log10(val),
         region = ifelse(loc_id %in% us_locs, "US", "Non-US"),
         adjusted = Dilution_Adjustment != "none")

# Same rule as R4_5_geospatial_scatter.R: drop location-analytes with < 20% detection
# (only PARENTs/Los Angeles NAP1, 0% detected).
d <- d %>% group_by(loc_id, Analyte_Code) %>% filter(mean(det) >= 0.2) %>% ungroup()
saveRDS(d, file.path(out, "pah_sample_level.rds"))  # reused by R4_5_PAH_cleanfuel.R

# ---- A. exposure contrast ----------------------------------------------------
pah_locs <- M %>% filter(set == "all") %>% distinct(loc_id) %>% pull()
A <- L %>% filter(loc_id %in% pah_locs) %>% transmute(loc_id, PM25 = PM25_ugm3, NO2 = NO2_ppb,
       region = ifelse(loc_id %in% us_locs, "US", "Non-US"))
us_med <- median(A$PM25[A$region == "US"])
A <- A %>% mutate(fold_vs_US_median = PM25 / us_med,
                  group = case_when(region == "US" ~ "US", grepl("^TEDDY_(Finland|Sweden|Germany)", loc_id) ~ "Europe",
                                    TRUE ~ "Latin America & Africa"))
write.csv(A, file.path(out, "A_PM25_at_PAH_locations.csv"), row.names = FALSE)

# ---- B. location level -------------------------------------------------------
loc_fit <- function(z) {
  z <- z %>% filter(is.finite(log10(median)), is.finite(log10(PM25)))
  if (nrow(z) < 4) return(tibble(n = nrow(z), slope = NA, lo = NA, hi = NA, rho = NA))
  f <- function(zz) coef(lm(log10(median) ~ log10(PM25), zz))[2]
  bs <- replicate(5000, { i <- sample(nrow(z), replace = TRUE); if (length(unique(z$loc_id[i])) < 3) NA else f(z[i, ]) })
  tibble(n = nrow(z), slope = f(z), lo = quantile(bs, .025, na.rm = TRUE), hi = quantile(bs, .975, na.rm = TRUE),
         rho = cor(z$median, z$PM25, method = "spearman"))
}
Bd <- M %>% inner_join(A %>% select(loc_id, PM25, region, group), by = "loc_id")
B <- bind_rows(lapply(PAH, function(a) bind_rows(
  loc_fit(Bd %>% filter(Analyte_Code == a, set == "all"))                         %>% mutate(subset = "All locations"),
  loc_fit(Bd %>% filter(Analyte_Code == a, set == "dilution_adjusted_only"))      %>% mutate(subset = "Dilution-adjusted urine only"),
  loc_fit(Bd %>% filter(Analyte_Code == a, set == "all", region == "US"))         %>% mutate(subset = "US only"),
  loc_fit(Bd %>% filter(Analyte_Code == a, set == "all", region == "Non-US"))     %>% mutate(subset = "Non-US only")
) %>% mutate(analyte = a))) %>%
  mutate(fold_per_doubling = 2^slope, fold_lo = 2^lo, fold_hi = 2^hi)
write.csv(B, file.path(out, "B_location_level_slopes.csv"), row.names = FALSE)

LOO <- bind_rows(lapply(PAH, function(a) {
  z <- Bd %>% filter(Analyte_Code == a, set == "all")
  bind_rows(lapply(z$loc_id, function(l) { zz <- z %>% filter(loc_id != l)
    tibble(analyte = a, dropped = l, rho = cor(zz$median, zz$PM25, method = "spearman"),
           slope = coef(lm(log10(median) ~ log10(PM25), zz))[2]) }))
}))
write.csv(LOO, file.path(out, "B_leave_one_location_out.csv"), row.names = FALSE)

# ---- C. individual-level between/within --------------------------------------
fit_bw <- function(dd, label) {
  if (n_distinct(dd$loc_id) < 4) return(NULL)  # same minimum as the location-level fits
  dd <- dd %>% group_by(loc_id) %>% mutate(pm_between = mean(lpm), pm_within = lpm - pm_between) %>% ungroup()
  has_within <- sd(dd$pm_within) > 0
  fml <- if (length(unique(dd$Dilution_Adjustment)) > 1) "lval ~ pm_between + pm_within + Dilution_Adjustment + (1|loc_id) + (1|pid)"
         else "lval ~ pm_between + pm_within + (1|loc_id) + (1|pid)"
  if (!has_within) fml <- sub(" \\+ pm_within", "", fml)
  m <- lmer(as.formula(fml), data = dd, REML = TRUE, control = lmerControl(calc.derivs = FALSE))
  cf <- summary(m)$coefficients; ci <- suppressMessages(confint(m, method = "Wald"))
  terms <- intersect(c("pm_between", "pm_within"), rownames(cf))
  tibble(analyte = unique(dd$Analyte_Code), model = label, term = terms, beta = cf[terms, 1], se = cf[terms, 2],
         lo = ci[terms, 1], hi = ci[terms, 2], n_samples = nrow(dd), n_participants = n_distinct(dd$pid),
         n_locations = n_distinct(dd$loc_id),
         sd_within_log10pm = sd(dd$pm_within), locations_with_within_variation = sum(tapply(dd$pm_within, dd$loc_id, sd) > 0, na.rm = TRUE))
}
C <- bind_rows(lapply(PAH, function(a) {
  da <- d %>% filter(Analyte_Code == a)
  bind_rows(fit_bw(da, "All locations"),
            fit_bw(da %>% filter(adjusted), "Dilution-adjusted urine only"),
            fit_bw(da %>% filter(region == "US"), "US only"),
            fit_bw(da %>% filter(region == "Non-US"), "Non-US only"))
})) %>% mutate(fold_per_doubling = 2^(beta * log10(2)), fold_lo = 2^(lo * log10(2)), fold_hi = 2^(hi * log10(2)))
# beta is per log10(PM2.5); per doubling = beta*log10(2) on the log10 PAH scale -> fold = 10^(beta*log10(2))
C <- C %>% mutate(fold_per_doubling = 10^(beta * log10(2)), fold_lo = 10^(lo * log10(2)), fold_hi = 10^(hi * log10(2)))
write.csv(C, file.path(out, "C_individual_between_within.csv"), row.names = FALSE)
# make B consistent: location slope is log10 PAH per log10 PM -> fold per doubling = 10^(slope*log10 2)
B <- B %>% mutate(fold_per_doubling = 10^(slope * log10(2)), fold_lo = 10^(lo * log10(2)), fold_hi = 10^(hi * log10(2)))
write.csv(B, file.path(out, "B_location_level_slopes.csv"), row.names = FALSE)

print(A %>% arrange(desc(PM25)) %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())
print(B %>% select(analyte, subset, n, rho, slope, fold_per_doubling, fold_lo, fold_hi) %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())
print(C %>% select(analyte, model, term, beta, lo, hi, fold_per_doubling, fold_lo, fold_hi, n_locations, locations_with_within_variation) %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())
print(LOO %>% group_by(analyte) %>% summarise(rho_min = min(rho), rho_max = max(rho), most_influential = dropped[which.max(abs(rho - median(rho)))]))

# ---- Figure -----------------------------------------------------------------
cols <- c(US = "#3b6fb6", Europe = "#2e8b57", `Latin America & Africa` = "#c0392b")
pA <- ggplot(A, aes(PM25, reorder(gsub("_", " ", loc_id), PM25), colour = group)) +
  geom_vline(xintercept = us_med, linetype = 2, colour = "grey50") + geom_point(size = 2.6) +
  geom_text(aes(label = sprintf("%.1f×", fold_vs_US_median)), hjust = -0.35, size = 2.6, show.legend = FALSE) +
  scale_colour_manual(values = cols, name = NULL) + scale_x_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  labs(title = "A. PM2.5 at the urinary-PAH locations", x = "PM2.5, pop-weighted, sampling years (µg/m³)", y = NULL,
       subtitle = "Labels: fold vs US median") + theme_bw(base_size = 9) + theme(legend.position = "bottom")
lab_B <- B %>% filter(subset == "All locations") %>% transmute(Analyte_Code = analyte,
         lab = sprintf("ρ = %.2f\n×%.2f per doubling\n(95%% CI %.2f–%.2f)", rho, fold_per_doubling, fold_lo, fold_hi))
pB <- ggplot(Bd %>% filter(set == "all"), aes(PM25, median)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", fill = "grey85", linewidth = 0.5) +
  geom_point(aes(colour = group, shape = df < 0.5), size = 2.3) +
  ggrepel::geom_text_repel(aes(label = gsub("_", " ", loc_id)), size = 2.1, colour = "grey30", seed = 1, max.overlaps = 30) +
  geom_text(data = lab_B, aes(label = lab), x = -Inf, y = Inf, hjust = -0.03, vjust = 1.15, size = 2.3, lineheight = 0.9, inherit.aes = FALSE) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("≥50% detected", "<50% detected"), name = NULL) +
  facet_wrap(~Analyte_Code, scales = "free_y") +
  labs(title = "B. Location medians vs PM2.5 (log–log, one point per location)", x = "PM2.5 (µg/m³)", y = "Location median (ng/mL or µg/g creatinine)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
fd <- bind_rows(B %>% transmute(analyte, model = subset, level = "Location level (between, ecological)", fold_per_doubling, fold_lo, fold_hi),
                C %>% filter(term == "pm_between") %>% transmute(analyte, model, level = "Individual level: between locations", fold_per_doubling, fold_lo, fold_hi),
                C %>% filter(term == "pm_within")  %>% transmute(analyte, model, level = "Individual level: within location (year-to-year)", fold_per_doubling, fold_lo, fold_hi)) %>%
  filter(is.finite(fold_per_doubling))
fd$model <- factor(fd$model, levels = rev(c("All locations", "Dilution-adjusted urine only", "US only", "Non-US only")))
pC <- ggplot(fd, aes(fold_per_doubling, model, colour = level)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = pmax(fold_lo, 0.05), xmax = pmin(fold_hi, 20)), position = position_dodge(width = 0.6), size = 0.25) +
  scale_x_log10() + facet_wrap(~analyte) + scale_colour_manual(values = c("#555555", "#c0392b", "#e08e0b"), name = NULL) +
  labs(title = "C. Fold change in urinary PAH per doubling of PM2.5", x = "Fold change (log scale; CIs truncated at 0.05–20)", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom") + guides(colour = guide_legend(nrow = 2))
g <- (pA | pB) / pC + plot_layout(widths = c(1, 2.2), heights = c(1.15, 1))
ggsave(file.path(out, "SuppFig_PAH_PM25_ProofOfConcept.png"), g, width = 13, height = 9.5, dpi = 200)
ggsave(file.path(out, "SuppFig_PAH_PM25_ProofOfConcept.pdf"), g, width = 13, height = 9.5, device = cairo_pdf)
