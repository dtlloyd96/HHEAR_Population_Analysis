# R4_5_B4_PYR_sensitivity.R
# Sensitivity for B4 (location contribution within urine): VDAART reports
# 1-hydroxypyrene as "PYR"; every other cohort uses "PYR1". In the canonical
# data they are separate analytes, so VDAART (576 children, US national) is
# missing from the PYR1 location model. This script refits the urine PAH
# location models with and without the PYR -> PYR1 remap, using the same engine
# and aggregation as B4 / R4_5_composition_context.R, and first checks that the
# un-remapped fit reproduces B4's published values.
# Canonical data and B4 outputs are NOT modified.
# Writes Output/Reviewer_Round2/R4_5_determinants/B4_PYR_sensitivity/. Run from project root.

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(quantreg); library(survival)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(file.path(root, "Coding"))
source("Paper1_Submission2/censored_quantile_helpers.R")
source("Paper1_Submission2/harmonization_helpers.R")
out <- file.path(root, "Output/Reviewer_Round2/R4_5_determinants/B4_PYR_sensitivity"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
CW <- read.csv(file.path(root, "Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv"))
PAH <- CW$Code[CW$Category == "PAH"]

r2_weighted <- function(res) {
  if (is.null(res) || inherits(res, "try-error")) return(NA_real_)
  t <- res$Regression_Results %>% distinct(Tau, pseudo_R); if (!nrow(t)) return(NA_real_)
  w <- exp(-((t$Tau - 0.5)^2) / (2 * 0.15^2)); 100 * weighted.mean(pmax(t$pseudo_R, 0, na.rm = TRUE), w, na.rm = TRUE)
}
e <- new.env(); load(file.path(root, "Output/Matrix_Stratified/Input_Urine/Harmonized_Targeted.RData"), envir = e)
s <- new.env(); load(file.path(root, "Output/Matrix_Stratified/Input_Urine/Harmonized_SDOH.RData"), envir = s)
SD <- s$SDOH_Data_All %>% mutate(Study = sub("_.*$", "", Study))

fit_version <- function(TD, label) {
  bind_rows(lapply(c(Child = "Child_PID", Mat = "Mat_PID"), function(id) {
    panel <- build_censored_targeted_panel(targeted_df = TD, id_var = id, unit_converter = to_ng_per_mL)
    an <- intersect(sort(unique(panel$Analyte_Code)), PAH)
    lst <- build_censored_targeted_list(panel_df = panel, id_var = id, analytes = an)
    ads <- SD %>% filter(!is.na(.data[[id]])) %>% distinct(.data[[id]], Study, Location)
    bind_rows(lapply(an, function(a) {
      dat <- lst[[a]] %>% mutate(Study = sub("_.*$", "", Study)) %>% select_outcome_data(join_var = id, outcome_col = a) %>%
        filter(!is.na(.data[[a]])) %>% left_join(ads, by = c(id, "Study")) %>% filter(!is.na(Location)) %>%
        mutate(!!a := as.numeric(.data[[a]]), Location = factor(Location)) %>% distinct()
      nl <- n_distinct(dat$Location)
      r <- if (nl >= 2) try(suppressWarnings(suppressMessages(run_censored_quantile_reg(NULL, dat, a, "Location", id, summary_R = 2))), silent = TRUE) else NULL
      tibble(version = label, panel = sub("_PID", "", id), analyte = a, n = nrow(dat), locations = nl,
             cohorts = paste(sort(unique(dat$Study)), collapse = "+"), R2 = r2_weighted(r))
    }))
  }))
}
TD0 <- e$Targeted_Data_All
TD1 <- TD0 %>% mutate(Analyte_Code = ifelse(Study == "VDAART_Targeted" & Analyte_Code == "PYR", "PYR1", Analyte_Code))
R <- bind_rows(fit_version(TD0, "canonical (PYR separate)"), fit_version(TD1, "PYR remapped to PYR1"))
write.csv(R, file.path(out, "PAH_urine_location_R2_by_panel.csv"), row.names = FALSE)

# B4 aggregation: analyte = mean over panels; class = mean over analytes within panel, then over panels
an_tab <- R %>% filter(!is.na(R2)) %>% group_by(version, analyte) %>% summarise(R2 = mean(R2), panels = n(), .groups = "drop")
cls <- R %>% filter(!is.na(R2)) %>% group_by(version, panel) %>% summarise(R2 = mean(R2), .groups = "drop") %>%
  group_by(version) %>% summarise(PAH_class_R2 = round(mean(R2), 1), .groups = "drop")
b4 <- read.csv(file.path(root, "Output/Reviewer_Round2/B4/location_by_analyte_pooled_vs_matrix.csv")) %>% filter(Category == "PAH") %>%
  transmute(analyte = Analyte_Code, B4_published = Urine)
cmp <- an_tab %>% select(version, analyte, R2) %>% pivot_wider(names_from = version, values_from = R2) %>% left_join(b4, by = "analyte")
write.csv(cmp, file.path(out, "PAH_urine_location_R2_comparison.csv"), row.names = FALSE)
write.csv(cls, file.path(out, "PAH_urine_class_R2.csv"), row.names = FALSE)
print(R %>% filter(analyte %in% c("PYR", "PYR1")) %>% as.data.frame(), digits = 3)
print(cmp %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())
print(cls); cat("B4 published PAH urine class:", read.csv(file.path(root, "Output/Reviewer_Round2/B4/location_by_class_pooled_vs_matrix.csv")) %>% filter(Category == "PAH") %>% pull(Urine), "\n")
