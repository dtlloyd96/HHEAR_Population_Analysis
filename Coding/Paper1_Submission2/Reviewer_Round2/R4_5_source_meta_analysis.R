# R4_5_source_meta_analysis.R
# Reviewer 4 -- does the same source raise the same chemical in every cohort?
# Recurring source -> chemical pairs (from R4_5_source_inventory.R) are tested
# within each cohort with the same participant-level model, then pooled with a
# random-effects meta-analysis (metafor::rma, REML). Heterogeneity (I2, Q) and a
# US vs non-US moderator test whether the effect differs by region.
#
# Smoking = any reported smoking (code >= 1) vs none (0). Definitions differ by
# cohort (current, during pregnancy, ever, caregiver, household), taken from:
#   harmonized SDOH Mat_Smoking / Mat_Smoker (CEHC, CHARGE, ELEMENT, MADRES,
#   NHBCS_Cot, PARENTs, SICAS, VIVA); DCHS raw smoke_during_pregnancy;
#   VAPS raw smoke_exposure; TSIC: child's tobacco-smoke exposure, cigarettes
#   per day at baseline > 0 vs 0 (TSIC caregivers are almost all smokers, so
#   caregiver smoking status has no contrast). MDUC excluded (3 smokers).
# Panel: the mother panel is used where a cohort has both (smoking is maternal).
# Cohort-analyte pairs with < 20% detection are dropped (as elsewhere).
# Pairs: smoking -> cotinine (one measure per cohort: urine COTT, else COT);
#        smoking -> NNAL; smoking -> 1-naphthol, 2-naphthol, 1-hydroxypyrene;
#        smoking -> cadmium (urine or blood, as measured);
#        pests -> 3-PBA (ZIP "insects a nuisance"; EEICC cockroaches seen).
# Model per cohort: lm(median log10 biomarker per participant ~ source); effect =
# fold difference (10^beta). Same value rules as the other R4_5 scripts.
# Writes Output/Reviewer_Round2/R4_5_determinants/source_meta_analysis/. Run from project root.

suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(metafor)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
source("Coding/Paper1_Submission2/harmonization_helpers.R")
out <- "Output/Reviewer_Round2/R4_5_determinants/source_meta_analysis"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
ld <- function(f) { e <- new.env(); load(file.path("Input", paste0(f, "_Data_All.RData")), envir = e); e[[ls(e)[1]]] }
t <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = t); TD <- t$Targeted_Data_All
s <- new.env(); load("Input/Harmonized_Datasets/Harmonized_SDOH.RData", envir = s); SD <- s$SDOH_Data_All
region <- c(DCHS = "Non-US", ELEMENT = "Non-US", ZIP = "Non-US")

vals <- function(study, analytes) TD %>% filter(Study == study, Analyte_Code %in% analytes) %>%
  mutate(id = coalesce(Mat_PID, Child_PID), panel = ifelse(!is.na(Mat_PID), "Mat", "Child"),
         conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
         det = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
         val = ifelse(det, conc, lod / sqrt(2))) %>% filter(!is.na(val), val > 0) %>%
  group_by(panel, id, Analyte_Code, Matrix_Class) %>% summarise(y = median(log10(val)), det = mean(det), .groups = "drop")

# smoking per participant, by both IDs, from harmonized SDOH
smk_sdoh <- SD %>% mutate(cohort = sub("_SDOH$", "", Study), sm = suppressWarnings(as.numeric(coalesce(as.character(Mat_Smoking), as.character(Mat_Smoker))))) %>%
  filter(!is.na(sm)) %>% mutate(x = as.integer(sm >= 1)) %>% select(cohort, Child_PID, Mat_PID, x)
smk_for <- function(cohort) {
  z <- smk_sdoh %>% filter(cohort == !!cohort)
  bind_rows(z %>% filter(!is.na(Mat_PID)) %>% transmute(panel = "Mat", id = as.character(Mat_PID), x),
            z %>% filter(!is.na(Child_PID)) %>% transmute(panel = "Child", id = as.character(Child_PID), x)) %>% distinct(panel, id, .keep_all = TRUE)
}
C <- ld("DCHS_00496")$DCHS_EPI
smk_raw <- list(
  DCHS = bind_rows(C %>% transmute(panel = "Mat", id = as.character(mother_pid), x = smoke_during_pregnancy),
                   C %>% transmute(panel = "Child", id = as.character(child_pid), x = smoke_during_pregnancy)),
  VAPS = ld("VAPS_2120")$VAPS_Epi %>% transmute(panel = "Child", id = as.character(PID), x = smoke_exposure),
  TSIC = { h <- ld("TSIC_1762")$TSIC_Harmonized
           tibble(panel = "Child", id = as.character(h[["Human-ID"]]),
                  x = as.integer(h[["Child-Human-Environmental-Tobacco-Smoke-Exposure-Measurement-Cigarettes-Per-Day-Childhood-Baseline"]] > 0)) %>%
             filter(!is.na(x)) %>% distinct(id, .keep_all = TRUE) })

fit_cohort <- function(cohort, study, analytes, expo, pair, outcome) {
  v <- vals(study, analytes); if (!nrow(v)) return(NULL)
  d <- v %>% inner_join(expo %>% filter(!is.na(x)), by = c("panel", "id"))
  if (!nrow(d)) return(NULL)
  # one analyte per cohort for the pair (first available in preference order); panel with most data
  a <- analytes[analytes %in% d$Analyte_Code][1]; d <- d %>% filter(Analyte_Code == a)
  pn <- if ("Mat" %in% d$panel) "Mat" else "Child"; d <- d %>% filter(panel == pn)
  if (mean(d$det) < 0.2) return(NULL)
  if (nrow(d) < 20 || min(table(d$x)) < 5) return(NULL)
  f <- lm(y ~ x, d); cf <- summary(f)$coefficients
  tibble(pair = pair, outcome = outcome, cohort = cohort, region = ifelse(cohort %in% names(region), "Non-US", "US"),
         analyte = a, matrix = first(d$Matrix_Class), panel = pn, n = nrow(d), n_exposed = sum(d$x == 1),
         detection = mean(d$det), yi = cf["x", 1], sei = cf["x", 2],
         exposure_type = ifelse(grepl("^Smoking", pair), ifelse(pn == "Mat", "Own smoking (mother sampled)", "Second-hand (child sampled)"), "Pests"))
}
smk <- function(co) if (co %in% names(smk_raw)) smk_raw[[co]] else smk_for(co)
st <- function(co) unique(TD$Study[sub("_Targeted$", "", TD$Study) == co])[1]

E <- list()
add <- function(cohorts, analytes, pair, outcome, expo_fun = smk) for (co in cohorts) {
  r <- tryCatch(fit_cohort(co, st(co), analytes, expo_fun(co), pair, outcome), error = function(e) NULL); if (!is.null(r)) E[[length(E) + 1]] <<- r }
add(c("DCHS", "SICAS", "NHBCS_Cot", "TSIC", "VAPS"), c("COTT", "COT"), "Smoking → cotinine", "Cotinine")
add(c("SICAS", "NHBCS_Cot", "VAPS"), "NNAL", "Smoking → NNAL", "NNAL")
for (a in c("NAP1", "NAP2", "PYR1")) add(c("DCHS", "SICAS", "ELEMENT", "PARENTs"), a, paste("Smoking →", a), a)
add(c("CEHC", "CHARGE", "MADRES", "VIVA"), "Cd", "Smoking → cadmium", "Cadmium")
# pests -> 3-PBA
Z <- ld("ZIP_1839"); H <- Z$ZIP_Harmonized
nuis <- grep("Pest-Infestation-Pregnancy-(Baseline|Followup-1)$", names(H), value = TRUE)
zip_pest <- tibble(panel = "Mat", id = as.character(H[["Human-ID"]]),
                   x = as.integer(apply(sapply(nuis, function(cn) as.character(H[[cn]]) == "1086"), 1, function(r) if (all(is.na(r))) NA else any(r, na.rm = TRUE))))
ee <- ld("EEICC_2512")$EEICC_Epi %>% transmute(panel = "Child", id = as.character(PID), x = roach_any)
add("ZIP", "PBA", "Pests at home → 3-PBA", "3-PBA", function(co) zip_pest)
add("EEICC", "PBA", "Pests at home → 3-PBA", "3-PBA", function(co) ee)

ES <- bind_rows(E) %>% mutate(fold = 10^yi, lo = 10^(yi - 1.96 * sei), hi = 10^(yi + 1.96 * sei))
write.csv(ES, file.path(out, "cohort_effects.csv"), row.names = FALSE)
print(ES %>% select(pair, cohort, analyte, matrix, panel, n, n_exposed, detection, yi, sei) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())

MA <- bind_rows(lapply(split(ES, ES$pair), function(z) {
  if (nrow(z) < 2) return(NULL)
  m <- tryCatch(rma(yi, sei = sei, data = z, method = "REML"), error = function(e) rma(yi, sei = sei, data = z, method = "DL"))
  mod <- if (n_distinct(z$region) == 2 && nrow(z) >= 3) tryCatch(rma(yi, sei = sei, mods = ~ region, data = z, method = "REML"),
                                                                 error = function(e) tryCatch(rma(yi, sei = sei, mods = ~ region, data = z, method = "DL"), error = function(e) NULL)) else NULL
  tibble(pair = z$pair[1], outcome = z$outcome[1], k = nrow(z), n_total = sum(z$n), cohorts = paste(z$cohort, collapse = ", "),
         pooled_fold = 10^m$b[1], lo = 10^m$ci.lb, hi = 10^m$ci.ub, p = m$pval, I2 = m$I2, tau2 = m$tau2, Q_p = m$QEp,
         region_p = if (!is.null(mod)) mod$QMp else NA_real_,
         exposure_type_p = { if (n_distinct(z$exposure_type) == 2 && nrow(z) >= 3) { me <- tryCatch(rma(yi, sei = sei, mods = ~ exposure_type, data = z, method = "REML"), error = function(e) NULL)
                               if (is.null(me)) NA_real_ else me$QMp } else NA_real_ },
         own_vs_secondhand_ratio = { if (n_distinct(z$exposure_type) == 2 && nrow(z) >= 3) { me <- tryCatch(rma(yi, sei = sei, mods = ~ exposure_type, data = z, method = "REML"), error = function(e) NULL)
                                       if (is.null(me)) NA_real_ else 10^(-me$b[2]) } else NA_real_ },
         region_ratio_nonUS_vs_US = if (!is.null(mod)) 10^(-mod$b[2]) else NA_real_)   # refs alphabetically: Non-US is reference
}))
write.csv(MA, file.path(out, "pooled_effects.csv"), row.names = FALSE)
print(ES %>% select(pair, cohort, region, analyte, matrix, panel, n, n_exposed, fold, lo, hi) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())
print(MA %>% select(-cohorts) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())

# ---- forest plot ----------------------------------------------------------------
lvl <- c("Smoking → cotinine", "Smoking → NNAL", "Smoking → NAP1", "Smoking → NAP2", "Smoking → PYR1", "Smoking → cadmium", "Pests at home → 3-PBA")
nice <- c("Smoking → cotinine" = "Smoking → cotinine", "Smoking → NNAL" = "Smoking → NNAL (tobacco nitrosamine)",
          "Smoking → NAP1" = "Smoking → 1-naphthol", "Smoking → NAP2" = "Smoking → 2-naphthol", "Smoking → PYR1" = "Smoking → 1-hydroxypyrene",
          "Smoking → cadmium" = "Smoking → cadmium", "Pests at home → 3-PBA" = "Pests at home → 3-PBA (pyrethroids)")
fp <- bind_rows(
  ES %>% transmute(pair, row = paste0(cohort, " (", region, ifelse(exposure_type == "Pests", "", ifelse(grepl("^Own", exposure_type), ", own", ", second-hand")), ", n = ", n, ")"), fold, lo, hi, type = "Cohort"),
  MA %>% transmute(pair, row = sprintf("Pooled (I² = %.0f%%%s%s)", I2, ifelse(is.na(region_p), "", sprintf("; region p = %.2g", region_p)),
                                       ifelse(is.na(exposure_type_p), "", sprintf("; own vs 2nd-hand p = %.2g", exposure_type_p))),
                   fold = pooled_fold, lo, hi, type = "Pooled")) %>%
  mutate(pair = factor(nice[pair], nice[lvl]))
fp$row <- factor(fp$row, levels = rev(c(unique(fp$row[fp$type == "Cohort"]), unique(fp$row[fp$type == "Pooled"]))))
g <- ggplot(fp, aes(fold, row)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = lo, xmax = hi, colour = type, shape = type), size = 0.3) +
  scale_shape_manual(values = c(Cohort = 16, Pooled = 18), name = NULL) + scale_colour_manual(values = c(Cohort = "grey35", Pooled = "#c0392b"), name = NULL) +
  scale_x_log10() + facet_wrap(~pair, scales = "free", ncol = 2) +
  labs(title = "Same source, several cohorts: random-effects meta-analysis", x = "Fold difference, exposed vs not (95% CI, log scale)", y = NULL) +
  theme_bw(base_size = 8.5) + theme(legend.position = "bottom")
ggsave(file.path(out, "Fig_Source_Meta_Analysis.png"), g, width = 10, height = 8.5, dpi = 200)
ggsave(file.path(out, "Fig_Source_Meta_Analysis.pdf"), g, width = 10, height = 8.5, device = cairo_pdf)
