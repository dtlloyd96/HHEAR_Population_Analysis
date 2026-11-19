# R4_5_within_cohort_sources.R
# Reviewer 4 -- within-cohort source attribution. Uses source questions each
# cohort collected from its own participants (same cohort, same lab), so lab /
# protocol / era differences cannot drive the associations.
#
# Cohort      Source variables (raw study file)                          Biomarkers (urine)
# DCHS        fossilstove, fossilheating (0/1); modelled home PM10, NO2,   PAHs (NAP1 NAP2 PYR1 FLUO2 PHEN1);
#             benzene, toluene (log2); smoke_during_pregnancy (0/1)        cotinine (COTT) for smoking
# MADRES      2nd-trimester seafood items (0/1; -99 = missing): fish        As, Hg
#             sticks, oily, fresh, canned, fried, shellfish; count of types
# NFGS        canned food, plastic-packaged food (FOODPACK), plastic-       phthalate metabolites
#             packaged beverages, all counts in prior 24 h; BPA-free
#             container use (code 594 vs 787425, ZIP/NFGS codebook)
# SICAS       hair_oil, hair_spray, lotion (0/1); live_home_cig            phthalates, parabens, phenols;
#             (smokers at home, 0 vs >=1)                                  COTT, NNAL, NAP1, NAP2, PYR1
# EEICC       hardwood_any (0/1; hardwood vs carpet/vinyl)                 phthalates
#             roach_any, roachbait_k (0/1)                                 pyrethroid / OP metabolites
# SICAS (+)   deodorant (0/1)                                              parabens, phenols, phthalates
# ZIP         tri_mosq_rep_ind_1 (uses mosquito repellent, baseline; 1 vs  pyrethroid metabolites (PBA, FPBA)
#             0; code 2 dropped)
# (+) items added after the systematic inventory (R4_5_source_inventory.R).
# DISCOVER is NOT used: its diet / product frequency items are coded 1-4 with
# no labels in the codebook, so direction cannot be interpreted. MADRES rice is
# not used (237/241 ate rice). EEICC has no lead, so peeling paint is not testable.
#
# Model: one value per participant (median of log10 values across samples;
# <LOD -> LOD/sqrt(2); to_ng_per_mL; canonical dilution basis of the cohort),
# lm(log10 biomarker ~ source). DCHS PAH models also adjusted for smoking in
# pregnancy (M2). Binary sources: fold difference exposed vs not; counts: per
# item; ambient: per doubling. Benjamini-Hochberg FDR over all tests. Each test
# has a pre-specified expected direction (`expect`).
# Writes Output/Reviewer_Round2/R4_5_determinants/within_cohort_sources/. Run from project root.

suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
source("Coding/Paper1_Submission2/harmonization_helpers.R")
out <- "Output/Reviewer_Round2/R4_5_determinants/within_cohort_sources"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
ld <- function(f) { e <- new.env(); load(file.path("Input", paste0(f, "_Data_All.RData")), envir = e); e[[ls(e)[1]]] }
t <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = t); TD <- t$Targeted_Data_All

part_vals <- function(study, id) TD %>% filter(Study == study, !is.na(.data[[id]])) %>%
  mutate(conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
         det = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
         val = ifelse(det, conc, lod / sqrt(2))) %>% filter(!is.na(val), val > 0) %>%
  group_by(pid = as.character(.data[[id]]), Analyte_Code) %>% summarise(y = median(log10(val)), det = mean(det), .groups = "drop")

fit_one <- function(dd, x, covs = character(0)) {
  dd <- dd %>% filter(!is.na(.data[[x]]), if_all(all_of(covs), ~ !is.na(.)))
  if (nrow(dd) < 20 || n_distinct(dd[[x]]) < 2) return(NULL)
  if (is.numeric(dd[[x]]) && n_distinct(dd[[x]]) == 2) if (min(table(dd[[x]])) < 5) return(NULL)
  f <- lm(as.formula(paste("y ~", paste(c(x, covs), collapse = " + "))), dd)
  cf <- summary(f)$coefficients; tn <- rownames(cf)[2]; ci <- confint(f)[tn, ]
  tibble(beta = cf[tn, 1], lo = ci[1], hi = ci[2], p = cf[tn, 4], n = nobs(f),
         n_exposed = if (n_distinct(dd[[x]]) == 2) sum(dd[[x]] == max(dd[[x]])) else NA_integer_)
}
run <- function(cohort, panel, vals, expo, tests, covs = character(0), adj_label = "unadjusted") {
  dat <- vals %>% inner_join(expo, by = "pid")
  bind_rows(lapply(seq_len(nrow(tests)), function(i) { tt <- tests[i, ]
    r <- fit_one(dat %>% filter(Analyte_Code == tt$analyte), tt$source, covs)
    if (is.null(r)) return(NULL)
    bind_cols(tibble(cohort = cohort, panel = panel, category = tt$category, source = tt$source, source_label = tt$label,
                     unit = tt$unit, analyte = tt$analyte, expect = tt$expect, model = adj_label,
                     detection = mean(dat$det[dat$Analyte_Code == tt$analyte])), r)
  }))
}
tests_for <- function(src, labels, units, analytes, category, expect) {
  expand_grid(i = seq_along(src), analyte = analytes) %>%
    mutate(source = src[i], label = labels[i], unit = units[i], category = category, expect = expect[i]) %>% select(-i)
}
an_present <- function(study, an) intersect(an, unique(TD$Analyte_Code[TD$Study == study]))

res <- list()
# ---- DCHS: household combustion, ambient air, smoking -> PAHs ------------------
C <- ld("DCHS_00496")$DCHS_EPI
PAH <- an_present("DCHS_Targeted", c("NAP1", "NAP2", "PYR1", "FLUO2", "PHEN1"))
for (pn in c("Mother" = "Mat_PID", "Child" = "Child_PID")) {
  idcol <- ifelse(pn == "Mat_PID", "mother_pid", "child_pid"); panel <- names(which(c(Mother = "Mat_PID", Child = "Child_PID") == pn))
  lg2 <- function(x) ifelse(!is.na(x) & x > 0, log2(x), NA_real_)  # non-positive modelled values treated as missing
  ex <- C %>% transmute(pid = as.character(.data[[idcol]]), fossilstove, fossilheating, smoke_preg = smoke_during_pregnancy,
                        log2_PM10 = lg2(pm10concgm31), log2_NO2 = lg2(ambientconcno2gm3_ant),
                        log2_benzene = lg2(ambientconcbenzenegm3), log2_toluene = lg2(ambientconctoluenegm3)) %>%
    filter(!is.na(pid)) %>% distinct(pid, .keep_all = TRUE)
  v <- part_vals("DCHS_Targeted", pn)
  comb <- tests_for(c("fossilstove", "fossilheating", "log2_PM10", "log2_NO2", "log2_benzene", "log2_toluene"),
                    c("Fossil-fuel cooking stove", "Fossil-fuel heating", "Home PM10 (modelled)", "Home NO2 (modelled)", "Home benzene (modelled)", "Home toluene (modelled)"),
                    c("yes vs no", "yes vs no", "per doubling", "per doubling", "per doubling", "per doubling"), PAH, "Household combustion & air", rep("higher", 6))
  smk <- tests_for("smoke_preg", "Smoking in pregnancy", "yes vs no", c(PAH, "COTT"), "Tobacco smoke", "higher")
  res[[length(res) + 1]] <- run("DCHS", panel, v, ex, bind_rows(comb, smk))
  res[[length(res) + 1]] <- run("DCHS", panel, v, ex, comb, covs = "smoke_preg", adj_label = "adjusted for smoking")
}
# ---- MADRES: seafood -> As, Hg ---------------------------------------------
Mx <- ld("MADRES_1945")$MADRES_Epi
sf <- c("t2_seafood_fsticks", "t2_seafood_oily", "t2_seafood_fresh", "t2_seafood_can", "t2_seafood_fried", "t2_seafood_shell")
ex <- Mx %>% transmute(pid = as.character(child_pid), across(all_of(sf), ~ ifelse(.x %in% c(0, 1), .x, NA))) %>%
  mutate(seafood_types = rowSums(across(all_of(sf))), any_seafood = as.integer(seafood_types > 0))
res[[length(res) + 1]] <- run("MADRES", "Mother", part_vals("MADRES_Targeted", "Mat_PID"), ex,
  tests_for(c("any_seafood", "seafood_types", "t2_seafood_shell", "t2_seafood_oily", "t2_seafood_can", "t2_seafood_fried", "t2_seafood_fsticks"),
            c("Any seafood", "Number of seafood types", "Shellfish", "Oily fish", "Canned fish", "Fried fish", "Fish sticks"),
            c("yes vs no", "per type", "yes vs no", "yes vs no", "yes vs no", "yes vs no", "yes vs no"),
            an_present("MADRES_Targeted", c("As", "Hg")), "Diet: seafood", rep("higher", 7)))
# ---- NFGS: packaging -> phthalates ---------------------------------------------
N <- ld("NFGS_2537"); h <- N$NFGS_Harmonized
ex <- N$NFGS_Epi %>% transmute(pid = as.character(PID), plastic_food = FOODPACK) %>%
  left_join(h %>% transmute(pid = as.character(`Human-ID`),
    canned_food = `Child-Human-Count-Of-Canned-Food-Consumption-Twenty-Four-Hours-Prior-To-Baseline`,
    plastic_bev = `Child-Human-Count-Of-Plastic-Packaged-Beverage-Consumption-Twenty-Four-Hours-Prior-To-Baseline`,
    no_bpa_free = case_when(`Child-Human-Uses-Plastic-Containers-For-Food-Or-Beverage-Storage-Twenty-Four-Hours-Prior-To-Baseline` == 787425 ~ 1L,
                            `Child-Human-Uses-Plastic-Containers-For-Food-Or-Beverage-Storage-Twenty-Four-Hours-Prior-To-Baseline` == 594 ~ 0L)), by = "pid") %>%
  mutate(across(c(plastic_food, canned_food, plastic_bev), ~ suppressWarnings(as.numeric(.x))))
PHT_N <- an_present("NFGS_Targeted", c("MEP", "MNBP", "MIBP", "MBZP", "MEHP", "MEHHP", "MEOHP", "MECPP", "MCPP", "MCOCH", "MINP", "MECPTP", "MEHHTP"))
res[[length(res) + 1]] <- run("NFGS", "Child", part_vals("NFGS_Targeted", "Child_PID"), ex,
  tests_for(c("plastic_food", "plastic_bev", "canned_food", "no_bpa_free"),
            c("Plastic-packaged foods (24 h)", "Plastic-packaged drinks (24 h)", "Canned foods (24 h)", "No effort to use BPA-free containers"),
            c("per item", "per item", "per item", "yes vs no"), PHT_N, "Diet: packaging", rep("higher", 4)))
# ---- SICAS: personal care and home smoking -------------------------------------
S <- ld("SICAS_1407")$SICAS_Epi
ex <- S %>% transmute(pid = as.character(PID), hair_oil, hair_spray, lotion, deodorant, smokers_home = as.integer(live_home_cig >= 1))
PC <- an_present("SICAS_Targeted", c("MEP", "MNBP", "MIBP", "MBZP", "MEHP", "MEHHP", "MEPB", "PRPB", "BUPB", "ETPB", "BP3", "TCS", "BPA"))
res[[length(res) + 1]] <- run("SICAS", "Child", part_vals("SICAS_Targeted", "Child_PID"), ex,
  bind_rows(tests_for(c("hair_oil", "hair_spray", "lotion", "deodorant"), c("Hair oil", "Hair spray", "Lotion", "Deodorant"), rep("yes vs no", 4), PC, "Personal care", rep("higher", 4)),
            tests_for("smokers_home", "Smoker(s) at home", "yes vs no", an_present("SICAS_Targeted", c("COTT", "NNAL", "NAP1", "NAP2", "PYR1")), "Tobacco smoke", "higher")))
# ---- EEICC: flooring -> phthalates ----------------------------------------------
E <- ld("EEICC_2512")$EEICC_Epi
ex <- E %>% transmute(pid = as.character(PID), hardwood_any, roach_any, roach_bait = as.integer(roachbait_k == "Yes"))
vE <- part_vals("EEICC_Targeted", "Child_PID")
res[[length(res) + 1]] <- run("EEICC", "Child", vE, ex,
  tests_for("hardwood_any", "Hardwood floors (vs carpet/vinyl)", "yes vs no",
            an_present("EEICC_Targeted", c("MBZP", "MEHP", "MEHHP", "MNBP", "MIBP", "MEP", "MINP")), "Housing", "lower"))
res[[length(res) + 1]] <- run("EEICC", "Child", vE, ex,
  tests_for(c("roach_any", "roach_bait"), c("Cockroaches seen", "Roach bait used"), c("yes vs no", "yes vs no"),
            an_present("EEICC_Targeted", c("PBA", "FPBA", "TCP", "PNP")), "Pests & pesticides", c("higher", "higher")))
# ---- ZIP: mosquito repellent -> pyrethroids ---------------------------------------
Z3 <- ld("ZIP_1839")$ZIP_EPI_3
ex <- Z3 %>% transmute(pid = as.character(MOTHER_PID), mosq_repellent = ifelse(tri_mosq_rep_ind_1 %in% c(0, 1), tri_mosq_rep_ind_1, NA)) %>%
  filter(!is.na(pid)) %>% distinct(pid, .keep_all = TRUE)
res[[length(res) + 1]] <- run("ZIP", "Mother", part_vals("ZIP_Targeted", "Mat_PID"), ex,
  tests_for("mosq_repellent", "Mosquito repellent use", "yes vs no", an_present("ZIP_Targeted", c("PBA", "FPBA")), "Pests & pesticides", "higher"))

R <- bind_rows(res) %>% mutate(fold = 10^beta, fold_lo = 10^lo, fold_hi = 10^hi)
R <- R %>% group_by(model) %>% mutate(q = p.adjust(p, "BH")) %>% ungroup() %>%
  mutate(direction_ok = (expect == "higher" & beta > 0) | (expect == "lower" & beta < 0),
         verdict = case_when(q < 0.05 & direction_ok ~ "Supported (FDR < 0.05)", p < 0.05 & direction_ok ~ "Suggestive (p < 0.05)",
                             p < 0.05 & !direction_ok ~ "Opposite direction", TRUE ~ "No association"))
write.csv(R, file.path(out, "within_cohort_source_results.csv"), row.names = FALSE)

S1 <- R %>% filter(model == "unadjusted") %>% group_by(category, cohort, source_label) %>%
  summarise(tests = n(), supported = sum(verdict == "Supported (FDR < 0.05)"), suggestive = sum(verdict == "Suggestive (p < 0.05)"),
            opposite = sum(verdict == "Opposite direction"),
            best = analyte[which.min(p)], best_fold = fold[which.min(p)], best_p = p[which.min(p)], .groups = "drop")
write.csv(S1, file.path(out, "within_cohort_source_summary.csv"), row.names = FALSE)
print(R %>% filter(p < 0.1) %>% arrange(category, p) %>% select(cohort, panel, model, source_label, analyte, fold, fold_lo, fold_hi, p, q, n, verdict) %>%
        mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())
print(S1 %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame())
cat("tests:", nrow(R %>% filter(model == "unadjusted")), "\n")

# ---- figure: forest of all unadjusted tests, grouped by source ------------------
pd <- R %>% filter(model == "unadjusted") %>%
  mutate(row = paste0(source_label, " — ", analyte, ifelse(cohort == "DCHS", paste0(" (", tolower(panel), ")"), "")),
         facet = paste0(category, "\n", cohort), verdict = factor(verdict, c("Supported (FDR < 0.05)", "Suggestive (p < 0.05)", "No association", "Opposite direction")))
g <- ggplot(pd, aes(fold, reorder(row, fold), colour = verdict)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = fold_lo, xmax = fold_hi), size = 0.18) + scale_x_log10() +
  facet_grid(facet ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("#c0392b", "#e08e0b", "grey55", "#3b6fb6"), name = NULL, drop = FALSE) +
  labs(title = "Within-cohort source attribution", subtitle = "Fold difference in urinary biomarker (95% CI); one model per source × biomarker",
       x = "Fold difference (log scale)", y = NULL) +
  theme_bw(base_size = 7.5) + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, size = 7), panel.spacing = unit(2, "pt"))
ggsave(file.path(out, "SuppFig_Within_Cohort_Sources.png"), g, width = 9, height = 0.13 * nrow(pd) + 2, dpi = 200, limitsize = FALSE)
ggsave(file.path(out, "SuppFig_Within_Cohort_Sources.pdf"), g, width = 9, height = 0.13 * nrow(pd) + 2, device = cairo_pdf, limitsize = FALSE)

# ---- compact key-results figure (pre-specified primary biomarker(s) per source) ----
key <- tribble(
  ~cohort, ~panel, ~source, ~analyte, ~model,
  "DCHS", "Mother", "smoke_preg", "COTT", "unadjusted",
  "DCHS", "Mother", "smoke_preg", "PYR1", "unadjusted",
  "DCHS", "Mother", "smoke_preg", "NAP1", "unadjusted",
  "SICAS", "Child", "smokers_home", "COTT", "unadjusted",
  "SICAS", "Child", "smokers_home", "NNAL", "unadjusted",
  "SICAS", "Child", "smokers_home", "PYR1", "unadjusted",
  "DCHS", "Mother", "fossilstove", "PYR1", "adjusted for smoking",
  "DCHS", "Mother", "fossilstove", "NAP1", "adjusted for smoking",
  "DCHS", "Child", "fossilstove", "PYR1", "adjusted for smoking",
  "DCHS", "Mother", "fossilheating", "NAP1", "adjusted for smoking",
  "DCHS", "Mother", "fossilheating", "PYR1", "adjusted for smoking",
  "DCHS", "Mother", "log2_PM10", "PHEN1", "unadjusted",
  "DCHS", "Mother", "log2_PM10", "PYR1", "unadjusted",
  "DCHS", "Mother", "log2_benzene", "PYR1", "unadjusted",
  "MADRES", "Mother", "any_seafood", "As", "unadjusted",
  "MADRES", "Mother", "any_seafood", "Hg", "unadjusted",
  "MADRES", "Mother", "t2_seafood_shell", "As", "unadjusted",
  "MADRES", "Mother", "t2_seafood_can", "Hg", "unadjusted",
  "MADRES", "Mother", "t2_seafood_fried", "As", "unadjusted",
  "NFGS", "Child", "no_bpa_free", "MEP", "unadjusted",
  "NFGS", "Child", "no_bpa_free", "MBZP", "unadjusted",
  "NFGS", "Child", "no_bpa_free", "MEHHP", "unadjusted",
  "NFGS", "Child", "canned_food", "MEP", "unadjusted",
  "NFGS", "Child", "plastic_food", "MEHHP", "unadjusted",
  "SICAS", "Child", "hair_spray", "MEP", "unadjusted",
  "SICAS", "Child", "hair_oil", "PRPB", "unadjusted",
  "SICAS", "Child", "lotion", "MEPB", "unadjusted",
  "SICAS", "Child", "deodorant", "MEPB", "unadjusted",
  "EEICC", "Child", "roach_any", "PBA", "unadjusted",
  "EEICC", "Child", "roach_bait", "PBA", "unadjusted",
  "ZIP", "Mother", "mosq_repellent", "PBA", "unadjusted",
  "EEICC", "Child", "hardwood_any", "MBZP", "unadjusted")
alab <- c(COTT = "cotinine", NNAL = "NNAL (tobacco)", PYR1 = "1-OH-pyrene", NAP1 = "1-naphthol", PHEN1 = "1-OH-phenanthrene",
          As = "arsenic", Hg = "mercury", PBA = "3-PBA (pyrethroids)", MEP = "MEP (DEP)", MBZP = "MBzP (BBzP)", MEHHP = "MEHHP (DEHP)", PRPB = "propylparaben", MEPB = "methylparaben")
kd <- key %>% inner_join(R, by = c("cohort", "panel", "source", "analyte", "model")) %>%
  mutate(row = paste0(source_label, " → ", alab[analyte], " (", cohort, if_else(cohort == "DCHS", paste0(", ", tolower(panel)), ""),
                      if_else(model == "adjusted for smoking", ", smoking-adj.", ""), ")"),
         category = factor(category, c("Tobacco smoke", "Household combustion & air", "Diet: seafood", "Diet: packaging", "Personal care", "Pests & pesticides", "Housing")),
         verdict = factor(verdict, c("Supported (FDR < 0.05)", "Suggestive (p < 0.05)", "No association", "Opposite direction")))
kd$row <- factor(kd$row, levels = rev(kd$row))
write.csv(kd, file.path(out, "within_cohort_key_results.csv"), row.names = FALSE)
gk <- ggplot(kd, aes(fold, row, colour = verdict)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = fold_lo, xmax = fold_hi), size = 0.25) + scale_x_log10(breaks = c(0.5, 1, 2, 5, 10, 30, 100)) +
  facet_grid(category ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = c("#c0392b", "#e08e0b", "grey55", "#3b6fb6"), name = NULL, drop = FALSE) +
  labs(title = "Within-cohort source attribution: key results", x = "Fold difference in urinary biomarker (95% CI, log scale)", y = NULL,
       caption = "Binary sources: exposed vs not; counts: per item; modelled air: per doubling. FDR (Benjamini-Hochberg) within each model set.") +
  theme_bw(base_size = 8.5) + theme(legend.position = "bottom", strip.text.y = element_text(angle = 0, size = 7.5))
ggsave(file.path(out, "Fig_Within_Cohort_Sources_Key.png"), gk, width = 9, height = 8.6, dpi = 200)
ggsave(file.path(out, "Fig_Within_Cohort_Sources_Key.pdf"), gk, width = 9, height = 8.6, device = cairo_pdf)
