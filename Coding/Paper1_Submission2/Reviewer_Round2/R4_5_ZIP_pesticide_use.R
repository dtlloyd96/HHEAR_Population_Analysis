# R4_5_ZIP_pesticide_use.R
# Reviewer 4, determinants of geography -- within-study test in ZIP (Brazil,
# Guatemala, Nicaragua, Puerto Rico; one lab, one protocol, n ~1,360 women).
# Question: do household pesticide use, pest pressure, recent use and
# agricultural/gardening work explain the between-country differences in
# urinary pesticide metabolites (DAPs, TCPy, 3-PBA, 4F-3-PBA, PNP)?
#
# Determinants (participant level, decoded via ZIP_Codebook / harmonized table):
#   indoor_use   "Uses Pesticide Inside the Home" (code 1088) at baseline OR follow-up 1
#   pest_nuis    "Insects Are A Nuisance In The Home" (1086) at baseline OR follow-up 1
#   use_48h      "Uses Pesticide" in the 48 h before baseline or follow-up (1021)
#   ag_house     mother's or partner's job title matches agriculture/gardening
#                (agricult*, campesin*, lavrador, finca, fazend*, ganad*, cultiv*,
#                productor, jardiner*, fumig*); "jornalero" (day labourer) NOT counted.
#   water        primary water source (utility / private well / bottled), from R4_5_water.R
#   Outdoor (lawn/garden) use is too rare (11/755) to model.
# These questions were answered by ~half of ZIP; analyses are complete-case on
# the determinant set, with the "country only" model refitted on the SAME rows.
#
# Models per metabolite (detection >= 50%), sample level, log10, <LOD -> LOD/sqrt(2):
#   M0: lval ~ country + (1|Mat_PID)
#   M1: lval ~ country + indoor_use + pest_nuis + use_48h + ag_house + water + (1|Mat_PID)
# Reported: determinant fold changes, LRT M1 vs M0, % attenuation of each
# country effect, and the share of between-country variance in fitted country
# means removed by the determinants.
# Writes Output/Reviewer_Round2/R4_5_determinants/zip_pesticide_use/. Run from project root.

suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(lme4); library(patchwork)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
source("Coding/Paper1_Submission2/harmonization_helpers.R")
out <- "Output/Reviewer_Round2/R4_5_determinants/zip_pesticide_use"; dir.create(out, showWarnings = FALSE, recursive = TRUE)

z <- new.env(); load("Input/ZIP_1839_Data_All.RData", envir = z); ZL <- z[[ls(z)[1]]]
E1 <- ZL$ZIP_EPI_1; H <- ZL$ZIP_Harmonized
hv <- function(pat) grep(pat, names(H), value = TRUE)
anyval <- function(cols, code) { m <- sapply(cols, function(cn) as.character(H[[cn]])); m <- matrix(m, nrow = nrow(H))
  ifelse(rowSums(!is.na(m)) == 0, NA, rowSums(m == as.character(code), na.rm = TRUE) > 0) }
indoor_cols <- hv("Indoor-Pesticide-Usage-Exposure-Pregnancy-(Baseline|Followup-1)$")
nuis_cols   <- hv("Pest-Infestation-Pregnancy-(Baseline|Followup-1)$")
use48_cols  <- hv("Pesticide-Usage-Exposure-Forty-Eight-Hours-Prior")
stopifnot(length(indoor_cols) == 2, length(nuis_cols) == 2, length(use48_cols) >= 1)
# Harmonized rows are aligned with ZIP_EPI_1 rows (same 1,464/1,468 order is NOT assumed):
# link through the water-source column, which exists in both, is not reliable either, so
# use the Mother/Child PID columns if present; otherwise fall back to EPI_1 raw variables.
id_col <- intersect(c("Mother-Human-Identifier", "Human-ID"), names(H))
cat("Harmonized ID columns:", id_col, "\n")
H_id <- if ("Mother-Human-Identifier" %in% names(H)) as.character(H[["Mother-Human-Identifier"]]) else as.character(H[["Human-ID"]])
det_h <- tibble(Mat_PID = H_id, indoor_use = anyval(indoor_cols, 1088), pest_nuis = anyval(nuis_cols, 1086),
                use_48h = anyval(use48_cols, 1021))
cat("Harmonized IDs matching EPI_1 MOTHER_PID:", sum(det_h$Mat_PID %in% as.character(E1$MOTHER_PID)), "of", nrow(det_h), "\n")

cl <- function(x) tolower(iconv(x, "latin1", "ASCII//TRANSLIT", sub = ""))
ag_pat <- "agricult|campesin|lavrador|finca|fazend|ganad|cultiv|productor|jardiner|fumig"
det_e <- E1 %>% transmute(Mat_PID = as.character(MOTHER_PID), country = tri_country,
  ag_house = grepl(ag_pat, cl(tri_husb_job_title)) | grepl(ag_pat, cl(tri_prim_job_title)),
  water = case_when(tri_prim_water_src == 1 ~ "Water utility", tri_prim_water_src == 2 ~ "Private well",
                    tri_prim_water_src == 7 ~ "Bottled/delivered", TRUE ~ NA_character_)) %>%
  filter(!is.na(Mat_PID)) %>% distinct(Mat_PID, .keep_all = TRUE)
D <- det_e %>% left_join(det_h %>% distinct(Mat_PID, .keep_all = TRUE), by = "Mat_PID")

prev <- D %>% filter(!is.na(country)) %>% group_by(country) %>%
  summarise(n = n(), answered_use_questions = mean(!is.na(indoor_use)), indoor_use = mean(indoor_use, na.rm = TRUE), pest_nuis = mean(pest_nuis, na.rm = TRUE),
            use_48h = mean(use_48h, na.rm = TRUE), ag_house = mean(ag_house), private_well = mean(water == "Private well", na.rm = TRUE),
            bottled = mean(water == "Bottled/delivered", na.rm = TRUE), .groups = "drop")
write.csv(prev, file.path(out, "determinant_prevalence_by_country.csv"), row.names = FALSE)

e <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = e)
ok_units <- c("ng/mL", "ng/ml", "ng/L", "pg/mL", "pg/ml", "ug/dL", "mg/L")
U <- e$Targeted_Data_All %>% filter(Study == "ZIP_Targeted", is.na(Units) | Units %in% ok_units) %>%
  mutate(conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
         det = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
         val = ifelse(det, conc, lod / sqrt(2))) %>% filter(!is.na(val), val > 0) %>%
  inner_join(D, by = "Mat_PID") %>% filter(!is.na(country)) %>%
  mutate(lval = log10(val), country = relevel(factor(country), "Puerto Rico"))
an <- U %>% group_by(Analyte_Code) %>% summarise(df = mean(det), .groups = "drop") %>% filter(df >= 0.5) %>% pull(Analyte_Code)

fitpair <- function(a) {
  dd <- U %>% filter(Analyte_Code == a, !is.na(indoor_use), !is.na(pest_nuis), !is.na(use_48h), !is.na(water)) %>%
    mutate(water = relevel(factor(water), "Water utility"))
  if (n_distinct(dd$country) < 2) return(NULL)
  ctl <- lmerControl(calc.derivs = FALSE)
  f0 <- lmer(lval ~ country + (1 | Mat_PID), dd, REML = FALSE, control = ctl)
  f1 <- lmer(lval ~ country + indoor_use + pest_nuis + use_48h + ag_house + water + (1 | Mat_PID), dd, REML = FALSE, control = ctl)
  b0 <- fixef(f0); b1 <- fixef(f1); ci <- suppressMessages(confint(f1, method = "Wald"))
  ct <- grep("^country", names(b0), value = TRUE); dt <- setdiff(names(b1), c("(Intercept)", ct))
  # spread of fitted country means (Puerto Rico = 0), before vs after
  v0 <- var(c(0, b0[ct])); v1 <- var(c(0, b1[ct]))
  list(
    det = tibble(analyte = a, term = dt, fold = 10^b1[dt], lo = 10^ci[dt, 1], hi = 10^ci[dt, 2]),
    cty = tibble(analyte = a, country = sub("country", "", ct), fold_country_only = 10^b0[ct], fold_adjusted = 10^b1[ct],
                 pct_attenuation = 100 * (1 - abs(b1[ct]) / abs(b0[ct]))),
    summ = tibble(analyte = a, n_participants = n_distinct(dd$Mat_PID), n_samples = nrow(dd), detection = mean(dd$det),
                  LRT_p = anova(f0, f1)$`Pr(>Chisq)`[2], pct_between_country_variance_removed = 100 * (1 - v1 / v0)))
}
R <- lapply(an, fitpair); R <- R[!sapply(R, is.null)]
DET <- bind_rows(lapply(R, `[[`, "det")); CTY <- bind_rows(lapply(R, `[[`, "cty")); SUM <- bind_rows(lapply(R, `[[`, "summ"))
write.csv(DET, file.path(out, "determinant_effects.csv"), row.names = FALSE)
write.csv(CTY, file.path(out, "country_effects_before_after.csv"), row.names = FALSE)
write.csv(SUM, file.path(out, "model_summary.csv"), row.names = FALSE)
print(prev %>% mutate(across(where(is.numeric) & !n, ~ round(.x, 3))) %>% as.data.frame())
print(SUM %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())
print(DET %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())
print(CTY %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())

# ---- figure --------------------------------------------------------------
lab <- c(indoor_useTRUE = "Indoor pesticide use", pest_nuisTRUE = "Insects a nuisance at home", use_48hTRUE = "Pesticide use in prior 48 h",
         ag_houseTRUE = "Agricultural / gardening job", `waterBottled/delivered` = "Bottled water (vs utility)", `waterPrivate well` = "Private well (vs utility)")
pA <- ggplot(DET %>% mutate(term = lab[term]), aes(fold, analyte)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi), size = 0.2, colour = "#7d3c98") + scale_x_log10() + facet_wrap(~term, nrow = 2) +
  labs(title = "A. Household determinants (mutually adjusted, within country)", x = "Fold difference in urinary metabolite (log)", y = NULL) +
  theme_bw(base_size = 9)
cd <- CTY %>% pivot_longer(c(fold_country_only, fold_adjusted), names_to = "model", values_to = "fold") %>%
  mutate(model = factor(model, c("fold_country_only", "fold_adjusted"), c("Country only", "+ household determinants")))
pB <- ggplot(cd, aes(fold, analyte, colour = model)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_point(position = position_dodge(width = 0.5), size = 1.8) + scale_x_log10() + facet_wrap(~country) +
  scale_colour_manual(values = c("grey40", "#e08e0b"), name = NULL) +
  labs(title = "B. Country differences vs Puerto Rico, before and after the determinants", x = "Fold vs Puerto Rico (log)", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
g <- pA / pB + plot_layout(heights = c(1.2, 1))
ggsave(file.path(out, "SuppFig_ZIP_Pesticide_Determinants.png"), g, width = 11, height = 9, dpi = 200)
ggsave(file.path(out, "SuppFig_ZIP_Pesticide_Determinants.pdf"), g, width = 11, height = 9, device = cairo_pdf)
