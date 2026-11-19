# R4_5_analysis_register.R
# Reviewer 4 -- labels every R4_5 analysis as pre-specified or exploratory, and
# applies Benjamini-Hochberg FDR across all location-level correlations.
#
# "Pre-specified" = named, with its expected direction, in
# Output/Paper1_Submission3/Text/PLAN_R4_Geographic_Determinants.md (sections 0,
# 3 and 5.3) before any result was seen. Substitute layers that replaced a
# planned one without looking at results (NO2 for PAH emissions; FAOSTAT
# insecticide use for PEST-CHEMGRIDS) keep the pre-specified label.
# Everything added after a result was seen is "Exploratory".
#
# FDR families, computed separately for each urine set (all locations = primary;
# dilution-adjusted urine only = sensitivity):
#   q_all     BH across every location-level test in that set
#   q_prespec BH across the pre-specified tests only in that set
# Writes Output/Reviewer_Round2/R4_5_determinants/register/. Run from project root.

suppressMessages({library(dplyr); library(tidyr); library(readr)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
det <- "Output/Reviewer_Round2/R4_5_determinants"
out <- file.path(det, "register"); dir.create(out, showWarnings = FALSE)
rd  <- function(...) read_csv(file.path(det, ...), show_col_types = FALSE)
# "0.50 (p = 0.12, n = 11)" -> rho, p, n
parse_res <- function(s) tibble(rho = as.numeric(sub("^\\s*(-?[0-9.]+).*", "\\1", s)),
                                perm_p = as.numeric(sub(".*p = ([0-9.eE-]+).*", "\\1", s)),
                                n_locations = as.integer(sub(".*n = ([0-9]+).*", "\\1", s)))

pre <- function(analyte, layer) case_when(
  layer %in% c("PM25_ugm3", "NO2_ppb") & analyte %in% c("NAP1", "NAP2", "PYR1", "Pb", "Cd") ~ "Pre-specified",
  layer == "GW_As_prob_gt10" & analyte == "As"                                        ~ "Pre-specified",
  layer == "insecticide_kg_ha" & analyte %in% c("TCP", "PBA")                         ~ "Pre-specified",
  layer == "rice_kg_cap" & analyte == "As"                                            ~ "Pre-specified",
  layer == "fish_kg_cap" & analyte == "Hg"                                            ~ "Pre-specified",
  layer == "years_since_lead_ban" & analyte == "Pb"                                   ~ "Pre-specified",
  TRUE                                                                                ~ "Exploratory")

lab <- c(PM25_ugm3 = "Ambient PM2.5", NO2_ppb = "Traffic NO2", pop_density_km2 = "Population density",
         GW_As_prob_gt10 = "Groundwater As", insecticide_kg_ha = "Insecticide use", pesticide_total_kg_ha = "Total pesticide use",
         fish_kg_cap = "Fish supply", rice_kg_cap = "Rice supply", years_since_lead_ban = "Years since leaded-petrol ban",
         polluting_fuel_pct = "Polluting cooking fuel")

tests_both <- bind_rows(
  rd("geospatial", "prespecified_pairs_spearman.csv") %>% 
    transmute(analysis = "Linked layers (geospatial)", set, analyte, layer, n_locations, rho, perm_p),
  rd("diet_pesticide_lead", "location_correlations.csv") %>% filter(block != "robustness", !grepl("excluded", layer)) %>%
    transmute(analysis = "Pesticide use, diet, lead", set, analyte, layer, n_locations, rho, perm_p),
  rd("geospatial", "pah_pm25", "D_cleanfuel_location_level.csv") %>% 
    transmute(analysis = "Cooking fuel", set, analyte, layer = "polluting_fuel_pct", n_locations = n, rho = rho_fuel, perm_p = p_fuel),
  rd("water", "1b_As_location_correlations.csv") %>% filter(grepl("^P\\(drinking", layer)) %>%
    bind_cols(parse_res(.$res)) %>% transmute(analysis = "Groundwater reliance", set = ifelse(set == "all", "all", "dilution_adjusted_only"), analyte = "As", layer = "GW As × reliance", n_locations, rho, perm_p),
  rd("water", "2_PFAS_location_correlations.csv") %>% filter(metric == "% PWS with own-analyte detected") %>%
    bind_cols(parse_res(.$res)) %>% transmute(analysis = "UCMR3 PFAS", set = "all", analyte, layer = "UCMR3 detections (own analyte)", n_locations, rho, perm_p)
) %>% filter(!is.na(perm_p), n_locations >= 4) %>%
  mutate(status = pre(analyte, layer), layer_label = coalesce(lab[layer], layer)) %>%
  group_by(set) %>% mutate(q_all = p.adjust(perm_p, "BH")) %>%
  group_by(set, status) %>% mutate(q_prespec = if_else(status == "Pre-specified", p.adjust(perm_p, "BH"), NA_real_)) %>% ungroup() %>%
  arrange(set, status != "Pre-specified", perm_p)
write_csv(tests_both, file.path(out, "location_level_tests_fdr.csv"))
tests <- tests_both %>% filter(set == "all")

# ---- analysis-level register -------------------------------------------------
WC <- rd("within_cohort_sources", "within_cohort_source_results.csv")
n_wc <- nrow(distinct(WC %>% filter(model == "unadjusted"), cohort, panel, source, analyte))
MA <- rd("source_meta_analysis", "pooled_effects.csv")
reg <- tribble(
  ~Part, ~Analysis, ~Status, ~Tests, ~`Multiple testing`,
  "A", "Composition vs context (age, sex, education, year)", "Pre-specified", "1 per class", "Descriptive (% of location R² retained)",
  "A", "ZIP household pesticide use, pests, jobs, water source", "Pre-specified", "Country effects before vs after", "Descriptive (attenuation)",
  "A", "MEBCO household water As → urinary As", "Pre-specified", "1", "Single test",
  "A", "Within-cohort reported sources (6 cohorts)", "Exploratory", as.character(n_wc), "BH FDR across all tests",
  "A", "Source-question inventory (23 cohorts)", "Exploratory", "None (descriptive)", "–",
  "A", "Same source across cohorts (meta-analysis)", "Exploratory", as.character(nrow(MA)), "Pooled; heterogeneity tests reported",
  "A", "Biomarker fingerprints", "Exploratory", "Descriptive + 1 validation per ratio", "Descriptive",
  "B", "PAH, Pb, Cd vs PM2.5; PAH vs NO2; As vs groundwater As", "Pre-specified", as.character(sum(tests$analysis == "Linked layers (geospatial)" & tests$status == "Pre-specified")), "BH FDR (q_prespec, q_all)",
  "B", "Pyrethroid/chlorpyrifos vs insecticide use; As vs rice; Hg vs fish; Pb vs lead ban", "Pre-specified", as.character(sum(tests$analysis == "Pesticide use, diet, lead" & tests$status == "Pre-specified")), "BH FDR (q_prespec, q_all)",
  "B", "Cooking fuel vs PAH", "Exploratory", as.character(sum(tests$analysis == "Cooking fuel")), "BH FDR (q_all)",
  "B", "Other OP metabolites, total pesticide use, As vs fish, population density", "Exploratory", as.character(sum(tests$analysis %in% c("Pesticide use, diet, lead", "Linked layers (geospatial)") & tests$status == "Exploratory")), "BH FDR (q_all)",
  "B", "Groundwater As × reliance; UCMR3 PFAS", "Exploratory", as.character(sum(tests$analysis %in% c("Groundwater reliance", "UCMR3 PFAS"))), "BH FDR (q_all)",
  "S", "Time matching, MEBCO creatinine, VDAART PYR, leave-one-out", "Sensitivity", "–", "Not tested for significance"
)
write_csv(reg, file.path(out, "analysis_register.csv"))

cat("location-level tests:", nrow(tests), "| pre-specified:", sum(tests$status == "Pre-specified"), "\n")
print(as.data.frame(tests_both %>% filter(q_all < 0.25 | perm_p < 0.05) %>% select(set, status, analyte, layer_label, n_locations, rho, perm_p, q_all, q_prespec)), digits = 2)
