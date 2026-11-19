# R4_5_MEBCO_creatinine_sensitivity.R
# Sensitivity: how much does creatinine division of MEBCO (Bangladesh) urine
# drive the arsenic / lead / trace-element geography?
# MEBCO creatinine is very low (median 29.6 mg/dL, min 0.3) and rises with water
# As, so ug/g = ug/L * 100 / CRE inflates and re-ranks values (R4_5_water.R).
# Canonical MEBCO values are verified to equal raw / CRE * 100 exactly.
#
# Versions of the MEBCO urine values (all other cohorts unchanged):
#   S0 canonical   ug/g creatinine (raw * 100 / CRE)
#   S1 unadjusted  ug/L (raw)
#   S2 CAS         covariate-adjusted standardisation (O'Brien et al. 2016,
#                  Environ Health Perspect): b = slope of log(raw) on log(CRE)
#                  within MEBCO (for As, adjusted for log water As);
#                  value = raw * (100 / CRE)^b, i.e. standardised to CRE = 100 mg/dL
#                  so it sits on the same numeric scale as ug/g but without
#                  assuming b = 1.
# Re-run for each version:
#   A. location medians (urine) for As, Cd, Mn, Mo, Pb
#   B. location-level correlations: As ~ P(groundwater As>10), As ~ P(drinking
#      groundwater As>10), Pb ~ years since leaded-petrol ban
#   C. B4 location contribution (urine, child panel -- where MEBCO sits -- and
#      the analyte mean over panels as in B4), same engine/aggregation as B4;
#      canonical refit checked against B4's published values.
# Writes Output/Reviewer_Round2/R4_5_determinants/MEBCO_creatinine/. Canonical data untouched.
# Run from project root.

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(quantreg); library(survival); library(ggplot2)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(file.path(root, "Coding"))
source("Paper1_Submission2/censored_quantile_helpers.R")
source("Paper1_Submission2/harmonization_helpers.R")
out <- file.path(root, "Output/Reviewer_Round2/R4_5_determinants/MEBCO_creatinine"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
MET <- c("As", "Cd", "Mn", "Mo", "Pb")

# ---- MEBCO raw values and the three versions --------------------------------
m <- new.env(); load(file.path(root, "Input/MEBCO_34_Data_All.RData"), envir = m); ML <- m[[ls(m)[1]]]
raw <- as.data.frame(ML$MEBCO_Targeted) %>% select(PID, all_of(MET), CRE) %>%
  left_join(as.data.frame(ML$MEBCO_Epi) %>% distinct(PID, .keep_all = TRUE) %>% select(PID, wAs = wateras_c1), by = "PID")
bfit <- bind_rows(lapply(MET, function(a) {
  z <- raw %>% filter(.data[[a]] > 0, CRE > 0)
  f <- if (a == "As") lm(log(z[[a]]) ~ log(z$CRE) + log(z$wAs)) else lm(log(z[[a]]) ~ log(z$CRE))
  ci <- confint(f)[2, ]
  tibble(analyte = a, b = coef(f)[2], b_lo = ci[1], b_hi = ci[2], n = nobs(f), adjusted_for = ifelse(a == "As", "log water As", "none"))
}))
write.csv(bfit, file.path(out, "CAS_creatinine_exponents.csv"), row.names = FALSE)
long <- raw %>% pivot_longer(all_of(MET), names_to = "Analyte_Code", values_to = "raw") %>% left_join(bfit %>% select(Analyte_Code = analyte, b), by = "Analyte_Code") %>%
  mutate(S0 = raw * 100 / CRE, S1 = raw, S2 = raw * (100 / CRE)^b, PID = as.character(PID))

e <- new.env(); load(file.path(root, "Output/Matrix_Stratified/Input_Urine/Harmonized_Targeted.RData"), envir = e)
s <- new.env(); load(file.path(root, "Output/Matrix_Stratified/Input_Urine/Harmonized_SDOH.RData"), envir = s)
TD0 <- e$Targeted_Data_All; SD <- s$SDOH_Data_All %>% mutate(Study = sub("_.*$", "", Study))
make_version <- function(v) {
  rep <- long %>% select(Child_PID = PID, Analyte_Code, new = all_of(v))
  TD0 %>% left_join(rep, by = c("Child_PID", "Analyte_Code")) %>%
    mutate(Concentration = ifelse(Study == "MEBCO_Targeted" & Analyte_Code %in% MET & !is.na(new), new, Concentration),
           Dilution_Adjustment = ifelse(Study == "MEBCO_Targeted" & Analyte_Code %in% MET, paste0("MEBCO_", v), Dilution_Adjustment)) %>% select(-new)
}
V <- list(S0_canonical_ugg = make_version("S0"), S1_unadjusted_ugL = make_version("S1"), S2_CAS = make_version("S2"))
chk <- V$S0_canonical_ugg %>% filter(Study == "MEBCO_Targeted", Analyte_Code %in% MET) %>% select(Child_PID, Analyte_Code, c0 = Concentration) %>%
  left_join(TD0 %>% filter(Study == "MEBCO_Targeted", Analyte_Code %in% MET) %>% select(Child_PID, Analyte_Code, c1 = Concentration), by = c("Child_PID", "Analyte_Code"))
stopifnot(isTRUE(all.equal(chk$c0, chk$c1, tolerance = 1e-8)))

# ---- A. location medians ------------------------------------------------------
loc_map <- SD %>% distinct(Study, Child_PID, Mat_PID, Location)
loc_medians <- function(TD, label) {
  d <- TD %>% filter(Analyte_Code %in% MET, Matrix_Class == "Urine") %>%
    mutate(conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
           det = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
           val = ifelse(det, conc, lod / sqrt(2)), pid = paste(Study, coalesce(Mat_PID, Child_PID)), Study = sub("_.*$", "", Study)) %>%
    filter(!is.na(val), val > 0, !grepl("^TEDDY$|^ZIP$", Study)) %>% left_join(loc_map %>% distinct(Study, Location), by = "Study")
  # TEDDY by child (ZIP has no metals)
  dt <- TD %>% filter(Analyte_Code %in% MET, Matrix_Class == "Urine", Study == "TEDDY_Targeted") %>%
    mutate(conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
           det = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
           val = ifelse(det, conc, lod / sqrt(2)), pid = paste(Study, Child_PID), Study = "TEDDY") %>% filter(!is.na(val), val > 0) %>%
    left_join(loc_map %>% filter(Study == "TEDDY") %>% distinct(Child_PID, Location), by = "Child_PID") %>% mutate(Location = paste0("TEDDY_", Location))
  bind_rows(d, dt) %>% filter(!is.na(Location)) %>% group_by(Analyte_Code, Location, pid) %>% summarise(v = median(val), det = mean(det), .groups = "drop") %>%
    group_by(Analyte_Code, loc_id = Location) %>% summarise(median = median(v), df = mean(det), n = n(), .groups = "drop") %>%
    filter(n >= 20, df >= 0.2) %>% group_by(Analyte_Code) %>% mutate(rank = rank(-median)) %>% ungroup() %>% mutate(version = label)
}
A <- bind_rows(lapply(names(V), function(v) loc_medians(V[[v]], v)))
write.csv(A, file.path(out, "A_location_medians_by_version.csv"), row.names = FALSE)

# ---- B. location-level correlations ------------------------------------------
gw  <- read.csv(file.path(root, "Output/Reviewer_Round2/R4_5_determinants/water/1b_groundwater_As_exposure_by_location.csv"))
lay <- read.csv(file.path(root, "Output/Reviewer_Round2/R4_5_determinants/diet_pesticide_lead/layers_by_location.csv"))
perm_rho <- function(x, y, B = 20000) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  r0 <- cor(x, y, method = "spearman"); set.seed(1); rp <- replicate(B, cor(x, sample(y), method = "spearman"))
  c(rho = r0, p = (1 + sum(abs(rp) >= abs(r0) - 1e-12)) / (B + 1), n = length(x)) }
AB <- A %>% left_join(gw %>% select(loc_id, P_As_popw, P_drink_As_gt10), by = "loc_id") %>% left_join(lay %>% select(loc_id, years_since_lead_ban), by = "loc_id")
B <- bind_rows(lapply(names(V), function(v) bind_rows(
  { z <- AB %>% filter(version == v, Analyte_Code == "As"); r <- perm_rho(z$P_As_popw, z$median); tibble(pair = "As ~ P(groundwater As>10)", rho = r[1], p = r[2], n = r[3]) },
  { z <- AB %>% filter(version == v, Analyte_Code == "As"); r <- perm_rho(z$P_drink_As_gt10, z$median); tibble(pair = "As ~ P(drinking groundwater As>10)", rho = r[1], p = r[2], n = r[3]) },
  { z <- AB %>% filter(version == v, Analyte_Code == "Pb"); r <- perm_rho(z$years_since_lead_ban, z$median); tibble(pair = "Pb ~ years since leaded-petrol ban", rho = r[1], p = r[2], n = r[3]) }
) %>% mutate(version = v)))
write.csv(B, file.path(out, "B_location_correlations_by_version.csv"), row.names = FALSE)

# ---- C. B4 location contribution for the MEBCO metals -------------------------
r2_weighted <- function(res) {
  if (is.null(res) || inherits(res, "try-error")) return(NA_real_)
  t <- res$Regression_Results %>% distinct(Tau, pseudo_R); if (!nrow(t)) return(NA_real_)
  w <- exp(-((t$Tau - 0.5)^2) / (2 * 0.15^2)); 100 * weighted.mean(pmax(t$pseudo_R, 0, na.rm = TRUE), w, na.rm = TRUE)
}
fitC <- function(TD, label) bind_rows(lapply(c(Child = "Child_PID", Mat = "Mat_PID"), function(id) {
  panel <- build_censored_targeted_panel(targeted_df = TD, id_var = id, unit_converter = to_ng_per_mL)
  an <- intersect(MET, unique(panel$Analyte_Code))
  lst <- build_censored_targeted_list(panel_df = panel, id_var = id, analytes = an)
  ads <- SD %>% filter(!is.na(.data[[id]])) %>% distinct(.data[[id]], Study, Location)
  bind_rows(lapply(an, function(a) {
    dat <- lst[[a]] %>% mutate(Study = sub("_.*$", "", Study)) %>% select_outcome_data(join_var = id, outcome_col = a) %>%
      filter(!is.na(.data[[a]])) %>% left_join(ads, by = c(id, "Study")) %>% filter(!is.na(Location)) %>%
      mutate(!!a := as.numeric(.data[[a]]), Location = factor(Location)) %>% distinct()
    nl <- n_distinct(dat$Location)
    r <- if (nl >= 2) try(suppressWarnings(suppressMessages(run_censored_quantile_reg(NULL, dat, a, "Location", id, summary_R = 2))), silent = TRUE) else NULL
    tibble(version = label, panel = sub("_PID", "", id), analyte = a, n = nrow(dat), locations = nl, has_MEBCO = "MEBCO" %in% dat$Study, R2 = r2_weighted(r))
  }))
}))
C <- bind_rows(lapply(names(V), function(v) fitC(V[[v]], v)))
write.csv(C, file.path(out, "C_location_R2_by_panel.csv"), row.names = FALSE)
b4 <- read.csv(file.path(root, "Output/Reviewer_Round2/B4/location_by_analyte_pooled_vs_matrix.csv")) %>% filter(Analyte_Code %in% MET) %>% select(analyte = Analyte_Code, B4_published_urine = Urine)
Cw <- C %>% filter(!is.na(R2)) %>% group_by(version, analyte) %>% summarise(R2 = mean(R2), .groups = "drop") %>%
  pivot_wider(names_from = version, values_from = R2) %>% left_join(b4, by = "analyte")
write.csv(Cw, file.path(out, "C_location_R2_comparison.csv"), row.names = FALSE)

print(bfit, width = Inf)
print(A %>% filter(loc_id == "Bangladesh") %>% select(Analyte_Code, version, median, rank) %>% pivot_wider(names_from = version, values_from = c(median, rank)) %>% as.data.frame(), digits = 3)
print(B %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())
print(Cw %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())
print(C %>% filter(panel == "Child") %>% select(version, analyte, R2) %>% pivot_wider(names_from = version, values_from = R2) %>% as.data.frame(), digits = 3)

# ---- figure: location medians by version (As, Pb) ----------------------------
fd <- A %>% filter(Analyte_Code %in% c("As", "Pb", "Cd", "Mn")) %>%
  mutate(version = factor(version, names(V), c("µg/g creatinine (canonical)", "Unadjusted µg/L", "Covariate-adjusted (CAS)")),
         bd = loc_id == "Bangladesh")
g <- ggplot(fd, aes(median, reorder(gsub("_", " ", loc_id), median), colour = version, shape = bd)) +
  geom_point(size = 2.2, position = position_dodge(width = 0.5)) + scale_x_log10() + facet_wrap(~Analyte_Code, scales = "free", nrow = 1) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17), guide = "none") + scale_colour_manual(values = c("grey40", "#e08e0b", "#7d3c98"), name = "MEBCO version") +
  labs(title = "Urinary metal location medians with MEBCO (Bangladesh, triangles) under three dilution treatments",
       subtitle = "Other cohorts unchanged", x = "Location median (log)", y = NULL) + theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out, "SuppFig_MEBCO_Creatinine_Sensitivity.png"), g, width = 12, height = 4.8, dpi = 200)
ggsave(file.path(out, "SuppFig_MEBCO_Creatinine_Sensitivity.pdf"), g, width = 12, height = 4.8, device = cairo_pdf)
