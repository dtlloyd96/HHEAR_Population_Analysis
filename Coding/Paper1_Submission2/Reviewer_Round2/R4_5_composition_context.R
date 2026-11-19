# R4_5_composition_context.R
# Reviewer 4, determinants of geography -- Part 1: composition vs context.
# Does location still explain concentrations once we account for WHO was
# sampled (age, sex, sampling year, education)?
#
# Same engine, data and pseudo-R2 aggregation as the published location models
# and B4: run_censored_quantile_reg() (censored quantile regression, taus chosen
# by choose_censored_taus), within a single matrix, on the matrix-stratified
# inputs (Output/Matrix_Stratified/Input_{Urine,Blood}; ug/g creatinine build of
# 2026-11-03); per-analyte R2 = Gaussian-tau-weighted (sd 0.15 around 0.5)
# mean of max(pseudo_R, 0), x100, exactly as B4_Matrix_Location_ugg.R.
#
# For each analyte x panel (Child / Mat) x matrix, on ONE complete-case sample:
#   L   ~ Location
#   C   ~ covariates
#   LC  ~ Location + covariates
# Location total (on this sample)  = R2(L)
# Location unique (beyond composition) = R2(LC) - R2(C)
# % of location signal explained by composition = 100 * (1 - unique / total)
# Covariate tiers (sampling year kept separate: it is largely a COHORT attribute
# -- era -- and is heavily aliased with location, so it is not "composition"):
#   T1 age_sex     : Child: age + sex;  Mat: age
#   T2 + education : T1 + education (3 levels: <= primary / secondary / > secondary;
#                    child panel uses household max education, else maternal)
#   T3 + year      : T2 + sampling year (era)
# Race/ethnicity is NOT used: it is a US-specific construct, harmonized only for
# US cohorts, and would drop nearly every non-US location from complete cases.
# Adult sex is omitted: the adult panel is mostly pregnant women and Adult_Sex is
# mostly unrecorded. Sampling year = participant median Sample_Collection_Year,
# as in Location_QuantileReg_Fixed.R (the SDOH `Year` field is a cohort constant).
# Inference is not needed (pseudo-R2 only), so summary_R = 2 skips the bootstrap.
#
# Requires >= 3 locations and >= 30 participants in the complete-case sample.
# Writes Output/Reviewer_Round2/R4_5_determinants/composition/. Run from project root.
# Env: COMP_CORES (default 8), COMP_TEST=1 to time a single analyte.

suppressPackageStartupMessages({library(dplyr); library(tidyr); library(quantreg); library(survival); library(parallel)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(file.path(root, "Coding"))
source("Paper1_Submission2/censored_quantile_helpers.R")
source("Paper1_Submission2/harmonization_helpers.R")
out <- file.path(root, "Output/Reviewer_Round2/R4_5_determinants/composition"); dir.create(out, recursive = TRUE, showWarnings = FALSE)
CORES <- as.integer(Sys.getenv("COMP_CORES", "8")); TEST <- Sys.getenv("COMP_TEST", "") == "1"

CW <- read.csv(file.path(root, "Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv"))
CLASSES_11 <- c("FlameRetardants","PAH","Parabens","PCB","Pesticides","PFAS","Phenols","Phthalates","Smoking","TraceElements","VOCs")
class_of <- setNames(CW$Category, CW$Code)

edu3 <- function(x) dplyr::case_when(
  x %in% c("No_Education", "Primary_Education", "Primary_School") ~ "1_Primary_or_less",
  x %in% c("High_School", "Secondary_Education") ~ "2_Secondary",
  x %in% c("Some_College", "Four_Year_College_Degree", "Some_Grad_School", "Graduate_School", "Some_Graduate_School") ~ "3_Post_secondary",
  TRUE ~ NA_character_)

r2_weighted <- function(res) {
  if (is.null(res) || inherits(res, "try-error")) return(NA_real_)
  t <- res$Regression_Results %>% distinct(Tau, pseudo_R)
  if (!nrow(t)) return(NA_real_)
  w <- exp(-((t$Tau - 0.5)^2) / (2 * 0.15^2))
  100 * weighted.mean(pmax(t$pseudo_R, 0, na.rm = TRUE), w, na.rm = TRUE)
}
fit_r2 <- function(dd, outcome, covs, id) {
  res <- try(suppressWarnings(suppressMessages(run_censored_quantile_reg(
    exp_var = NULL, dataset = dd, my_outcome = outcome, my_covariates = covs, join_var = id, summary_R = 2))), silent = TRUE)
  r2_weighted(res)
}

run_matrix <- function(mx) {
  e <- new.env(); load(file.path(root, "Output/Matrix_Stratified", paste0("Input_", mx), "Harmonized_Targeted.RData"), envir = e)
  s <- new.env(); load(file.path(root, "Output/Matrix_Stratified", paste0("Input_", mx), "Harmonized_SDOH.RData"), envir = s)
  TD <- e$Targeted_Data_All; SD <- s$SDOH_Data_All
  SD <- SD %>% mutate(Edu3 = edu3(dplyr::coalesce(Max_Edu, Mat_Edu)), Mat_Age = suppressWarnings(as.numeric(Mat_Age)),
                      Child_Age = suppressWarnings(as.numeric(Child_Age)), Study = sub("_SDOH$", "", Study))
  yr <- function(id) TD %>% filter(!is.na(.data[[id]])) %>% group_by(.data[[id]]) %>%
    summarise(Samp_Year = { y <- suppressWarnings(as.numeric(Sample_Collection_Year)); if (all(is.na(y))) NA_real_ else median(y, na.rm = TRUE) }, .groups = "drop")
  panels <- list(
    Child = list(id = "Child_PID", covT1 = c("Child_Age", "Child_Sex"),
                 ads = SD %>% filter(!is.na(Child_PID)) %>% distinct(Child_PID, Study, .keep_all = TRUE) %>%
                   select(Child_PID, Study, Location, Child_Age, Child_Sex, Edu3) %>% left_join(yr("Child_PID"), by = "Child_PID")),
    Mat   = list(id = "Mat_PID", covT1 = c("Mat_Age"),
                 ads = SD %>% filter(!is.na(Mat_PID)) %>% distinct(Mat_PID, Study, .keep_all = TRUE) %>%
                   select(Mat_PID, Study, Location, Mat_Age, Edu3) %>% left_join(yr("Mat_PID"), by = "Mat_PID")))
  bind_rows(lapply(names(panels), function(pn) {
    P <- panels[[pn]]; id <- P$id
    panel <- build_censored_targeted_panel(targeted_df = TD, id_var = id, unit_converter = to_ng_per_mL)
    an <- sort(unique(panel$Analyte_Code)); an <- an[class_of[an] %in% CLASSES_11]
    if (TEST) an <- intersect(c("Pb", "PFOS", "PYR1"), an)[1]
    lst <- build_censored_targeted_list(panel_df = panel, id_var = id, analytes = an)
    ads <- P$ads %>% mutate(Study = sub("_.*$", "", Study))
    one <- function(a) {
      dat <- lst[[a]] %>% mutate(Study = sub("_.*$", "", Study)) %>% select_outcome_data(join_var = id, outcome_col = a) %>%
        filter(!is.na(.data[[a]])) %>% left_join(ads, by = c(id, "Study")) %>% mutate(!!a := as.numeric(.data[[a]])) %>% distinct()
      dat$Location <- factor(dat$Location)
      full_L <- dat %>% filter(!is.na(Location))
      r_full <- if (n_distinct(full_L$Location) >= 3) fit_r2(full_L, a, "Location", id) else NA_real_
      tiers <- list(T1_age_sex = P$covT1, T2_plus_education = c(P$covT1, "Edu3"), T3_plus_year = c(P$covT1, "Edu3", "Samp_Year"))
      bind_rows(lapply(names(tiers), function(tn) {
        cv <- tiers[[tn]]
        cc <- dat %>% filter(!is.na(Location), if_all(all_of(cv), ~ !is.na(.))) %>% droplevels()
        base <- tibble(analyte = a, class = class_of[[a]], matrix = mx, panel = pn, tier = tn,
                       n_full = nrow(full_L), locations_full = n_distinct(full_L$Location), R2_L_fullsample = r_full,
                       n_cc = nrow(cc), locations_cc = n_distinct(cc$Location),
                       locations_dropped = paste(setdiff(levels(droplevels(full_L$Location)), unique(as.character(cc$Location))), collapse = ";"))
        if (nrow(cc) < 30 || n_distinct(cc$Location) < 3) return(base)
        base %>% mutate(R2_L = fit_r2(cc, a, "Location", id), R2_C = fit_r2(cc, a, cv, id), R2_LC = fit_r2(cc, a, c("Location", cv), id))
      }))
    }
    t0 <- Sys.time()
    res <- if (CORES > 1 && !TEST) mclapply(an, function(a) tryCatch(one(a), error = function(err) tibble(analyte = a, error = conditionMessage(err))), mc.cores = CORES)
           else lapply(an, one)
    message(sprintf("[%s %s] %d analytes in %.1f min", mx, pn, length(an), as.numeric(difftime(Sys.time(), t0, units = "mins"))))
    bind_rows(res)
  }))
}

R <- bind_rows(lapply(c("Urine", "Blood"), run_matrix)) %>%
  mutate(R2_L_unique = R2_LC - R2_C, pct_location_explained_by_composition = 100 * (1 - R2_L_unique / R2_L),
         # censored pseudo-R2 is not strictly nested (LC can fall below C), so a clipped
         # version bounded to [0, 100] is used for summaries; raw values are kept.
         pct_explained_clipped = pmin(pmax(pct_location_explained_by_composition, 0), 100))
if (TEST) { print(as.data.frame(R)); quit(save = "no") }
write.csv(R, file.path(out, "composition_by_analyte.csv"), row.names = FALSE)

# class-level summary (mean over analytes and panels), as in B4
S <- R %>% filter(!is.na(R2_L), !is.na(R2_LC), !is.na(R2_C)) %>% group_by(class, matrix, tier) %>%
  summarise(n_analyte_panels = n(), R2_L_fullsample = mean(R2_L_fullsample, na.rm = TRUE), R2_L = mean(R2_L), R2_C = mean(R2_C),
            R2_LC = mean(R2_LC), R2_L_unique = mean(R2_L_unique),
            median_pct_explained = median(pct_location_explained_by_composition, na.rm = TRUE),
            median_pct_explained_clipped = median(pct_explained_clipped, na.rm = TRUE),
            median_locations_cc = median(locations_cc), .groups = "drop") %>%
  mutate(pct_explained_of_mean = 100 * (1 - R2_L_unique / R2_L))
write.csv(S, file.path(out, "composition_by_class.csv"), row.names = FALSE)
print(S %>% mutate(across(where(is.numeric), ~ round(.x, 1))) %>% as.data.frame())

# ---- figure: location R2 before vs after each covariate tier (class means) ----
suppressPackageStartupMessages(library(ggplot2))
S2 <- S %>% mutate(pct_retained = 100 * R2_L_unique / R2_L, label = paste0(class, " (", tolower(matrix), ")"),
                   tier = factor(tier, c("T1_age_sex", "T2_plus_education", "T3_plus_year"),
                                 c("Age + sex", "+ education", "+ sampling year")))
write.csv(S2 %>% select(class, matrix, tier, n_analyte_panels, R2_L, R2_L_unique, pct_retained, median_pct_explained_clipped, median_locations_cc),
          file.path(out, "composition_pct_location_retained.csv"), row.names = FALSE)
pd <- S2 %>% select(label, tier, n_analyte_panels, `Location alone` = R2_L, `Location beyond covariates` = R2_L_unique) %>%
  pivot_longer(c(`Location alone`, `Location beyond covariates`), names_to = "what", values_to = "R2")
g <- ggplot(pd, aes(pmax(R2, 0), reorder(label, R2), fill = what)) + geom_col(position = position_dodge(width = 0.75), width = 0.7) +
  geom_text(data = S2, aes(x = pmax(R2_L, R2_L_unique) + 0.6, y = label, label = sprintf("%.0f%% kept (n=%d)", pmax(pct_retained, 0), n_analyte_panels)),
            inherit.aes = FALSE, hjust = 0, size = 2.4) +
  facet_wrap(~tier, nrow = 1) + scale_fill_manual(values = c("grey70", "#c0392b"), name = NULL) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.45))) +
  labs(title = "Location contribution before and after adjusting for who was sampled",
       subtitle = "Class means of tau-weighted pseudo-R² (%), within matrix, each tier on its own complete-case sample",
       x = "Pseudo-R² (%)", y = NULL,
       caption = "'Kept' = location R² beyond covariates as % of location R² alone. Education and sampling year are largely cohort-level\nattributes with 3–6 locations per class, so the drops in the right-hand panels are partly aliasing with location, not composition.") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
ggsave(file.path(out, "SuppFig_Composition_vs_Context.png"), g, width = 12, height = 5.2, dpi = 200)
ggsave(file.path(out, "SuppFig_Composition_vs_Context.pdf"), g, width = 12, height = 5.2, device = cairo_pdf)
