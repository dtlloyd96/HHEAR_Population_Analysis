# C2 -- Calendar-year trends at the level of single chemicals
# (Round 2: Reviewer 4, comment 4 -- "Aggregating at the class level is likely
#  to mask temporal trends that run in opposite directions for individual
#  chemicals ... would be more informative at the level of single chemicals ...
#  linking the observed trends to known regulatory actions")
#
# No new models. Reads the existing within-cohort sampling-year models,
# By_Study_Single_Exp_Results_{Child,Mat}.RData (Within_Cohort_SDOH_QuantileReg.Rmd:
# censored quantile regression of concentration on Sample_Collection_Year, fit
# separately in every cohort). These are within-cohort, so the 2026-11-03
# ug/g creatinine change does not affect them (verified: coefficients identical
# up to the constant cohort scale).
#
# Slopes are on the concentration scale, so for comparison across chemicals each
# is expressed as % of that cohort's median concentration per year (scale-free,
# hence also invariant to the ug/g x100). Only cohorts spanning >= MIN_SPAN
# distinct sampling years are used; a slope over a 1-2 year window says nothing
# about secular change.
#
# Writes Output/Reviewer_Round2/C2/. Run from the project root.

suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("Coding/Paper1_Submission2/harmonization_helpers.R")

OUT <- "Output/Reviewer_Round2/C2"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
RES <- "Output/Paper1_Output_Submission2"
MIN_SPAN <- 4          # distinct sampling years in the cohort-panel
TAU      <- 0.5
ALPHA    <- 0.05

CW <- fread("Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv")
CLASSES_11 <- c("FlameRetardants","PAH","Parabens","PCB","Pesticides","PFAS",
                "Phenols","Phthalates","Smoking","TraceElements","VOCs")

# ---- 1. year slopes from the existing within-cohort models ---------------------
pull <- function(file, panel) {
  e <- new.env(); load(file, envir = e); L <- get(ls(e)[1], envir = e)
  rbindlist(lapply(names(L), function(coh) rbindlist(lapply(names(L[[coh]]), function(an) {
    r <- L[[coh]][[an]][["Sample_Collection_Year"]]
    if (is.null(r) || !is.data.frame(r$Regression_Results)) return(NULL)
    d <- as.data.table(r$Regression_Results)[Variable == "Sample_Collection_Year"]
    if (!nrow(d)) return(NULL)
    d[, .(Cohort = coh, Analyte_Code = an, Panel = panel, Tau, Beta, Lower, Upper, P_Value,
          pseudo_R, Pseudo_R_Overall, N = Model_N)]
  }), fill = TRUE)), fill = TRUE)
}
S <- rbind(pull(file.path(RES, "By_Study_Single_Exp_Results_Child.RData"), "Child"),
           pull(file.path(RES, "By_Study_Single_Exp_Results_Mat.RData"),   "Adult"), fill = TRUE)
S <- merge(S, unique(CW[Category %in% CLASSES_11, .(Analyte_Code = Code, Category)]), by = "Analyte_Code")

# ---- 2. cohort medians and year spans ------------------------------------------
# The medians MUST be on the scale the models were fit on. The year models were
# fit 2026-08-11 on the pre-ug/g data (creatinine cohorts divided by mg/dL), so
# read that build; dividing an old-scale slope by a ug/g median would understate
# % change 100-fold for every creatinine-adjusted cohort. Re-point this at the
# canonical file once By_Study_Single_Exp_Results_* are re-run on ug/g data.
# 2026-11-18: the year models were re-run on the ug/g build (after the fix that
# had dropped 11 cohorts), so the medians now come from the canonical file.
MODEL_DATA <- "Input/Harmonized_Datasets/Harmonized_Targeted.RData"
stopifnot(file.mtime(MODEL_DATA) < file.mtime(file.path(RES, "By_Study_Single_Exp_Results_Mat.RData")))
load(MODEL_DATA)
H <- as.data.table(Targeted_Data_All)
H[, `:=`(Cohort = std_study(Study), Panel = fifelse(!is.na(Child_PID), "Child", "Adult"),
         pid = fcoalesce(as.character(Child_PID), as.character(Mat_PID)),
         conc = to_ng_per_mL(Concentration, fifelse(is.na(Units), "ng/mL", Units)))]
med  <- H[!is.na(conc), .(conc = mean(conc)), by = .(Cohort, Panel, Analyte_Code, pid)][
          , .(median_conc = median(conc)), by = .(Cohort, Panel, Analyte_Code)]
span <- H[!is.na(Sample_Collection_Year), .(n_years = uniqueN(Sample_Collection_Year),
          first_year = min(Sample_Collection_Year), last_year = max(Sample_Collection_Year)),
          by = .(Cohort, Panel)]
S <- merge(S, med,  by = c("Cohort", "Panel", "Analyte_Code"), all.x = TRUE)
S <- merge(S, span, by = c("Cohort", "Panel"), all.x = TRUE)
S <- S[n_years >= MIN_SPAN & median_conc > 0]
S[, `:=`(pct_per_year = 100 * Beta / median_conc,
         pct_lo = 100 * Lower / median_conc, pct_hi = 100 * Upper / median_conc)]
fwrite(S[order(Category, Analyte_Code, Cohort, Panel, Tau)], file.path(OUT, "year_slopes_all_taus.csv"))

M <- S[Tau == TAU]
# a fit with no inference (P_Value NA) counts as no significant trend
M[, direction := fifelse(P_Value %in% P_Value[P_Value < ALPHA] & Beta < 0, "Decreasing",
                 fifelse(P_Value %in% P_Value[P_Value < ALPHA] & Beta > 0, "Increasing", "No significant trend"))]
M[is.na(P_Value) | P_Value >= ALPHA, direction := "No significant trend"]
fwrite(M[order(Category, Analyte_Code, Cohort)], file.path(OUT, "year_slopes_median.csv"))

# ---- 3. does the class average hide opposite-direction chemicals? ---------------
by_class <- M[, .(chemical_cohort_pairs = .N, chemicals = uniqueN(Analyte_Code),
                  cohorts = uniqueN(Cohort),
                  decreasing = sum(direction == "Decreasing"),
                  increasing = sum(direction == "Increasing"),
                  not_significant = sum(direction == "No significant trend")), by = Category]
by_class[, both_directions := decreasing > 0 & increasing > 0]
# Within a single cohort-panel-class: chemicals moving in opposite directions at once
opp <- M[direction != "No significant trend",
         .(dec = sum(direction == "Decreasing"), inc = sum(direction == "Increasing"),
           dec_chems = paste(sort(Analyte_Code[direction == "Decreasing"]), collapse = ","),
           inc_chems = paste(sort(Analyte_Code[direction == "Increasing"]), collapse = ",")),
         by = .(Cohort, Panel, Category)][dec > 0 & inc > 0]
fwrite(by_class[order(-chemical_cohort_pairs)], file.path(OUT, "direction_by_class.csv"))
fwrite(opp, file.path(OUT, "opposite_direction_within_cohort_class.csv"))

# ---- 4. legacy vs replacement chemicals with regulatory context ----------------
# Dates are the commonly cited milestones; confirm each against the primary
# source before quoting in the manuscript.
FAMILIES <- rbindlist(list(
  data.table(family = "PFAS", role = "Legacy",
             Analyte_Code = c("PFOS", "PFOA"),
             milestone = "PFOS: 3M phase-out 2000-02; PFOA: EPA Stewardship 2006-15"),
  data.table(family = "PFAS", role = "Other / replacement",
             Analyte_Code = c("PFHXS", "PFNA", "PFDA", "PFUNDA", "PFHPA", "PFBS", "PFHXA"),
             milestone = ""),
  data.table(family = "Phthalates", role = "Legacy (DEHP, DBP, BBP)",
             Analyte_Code = c("MEHP", "MEHHP", "MEOHP", "MECPP", "MNBP", "MIBP", "MBZP"),
             milestone = "US CPSIA 2008 (toys); EU REACH authorisation 2015"),
  data.table(family = "Phthalates", role = "Replacement (DINP, DEHTP, DINCH)",
             Analyte_Code = c("MCOP", "MINP", "MCINP", "MCIOP", "MECPTP", "MEHHTP", "MHNCH", "MCOCH"),
             milestone = ""),
  data.table(family = "Bisphenols", role = "Legacy (BPA)", Analyte_Code = "BPA",
             milestone = "EU 2011 / FDA 2012 baby-bottle restrictions"),
  data.table(family = "Bisphenols", role = "Replacement (BPS, BPF)", Analyte_Code = c("BPS", "BPF", "BPAF"),
             milestone = ""),
  data.table(family = "Organophosphates", role = "Chlorpyrifos (TCP)", Analyte_Code = "TCP",
             milestone = "US residential phase-out 2000-01; food tolerances revoked 2021"),
  data.table(family = "Organophosphates", role = "Dialkyl phosphates",
             Analyte_Code = c("DEP", "DETP", "DMP", "DMTP", "DEDP", "DMDP"), milestone = ""),
  data.table(family = "Pyrethroids", role = "Pyrethroids (3-PBA)", Analyte_Code = c("PBA", "DCCA", "TDCCA"),
             milestone = "rising use as OP replacements")))
F <- merge(M, FAMILIES, by = "Analyte_Code")
fwrite(F[order(family, role, Analyte_Code, Cohort)], file.path(OUT, "year_slopes_legacy_vs_replacement.csv"))

# ---- 5. figures ------------------------------------------------------------------
dir_pal <- c("Decreasing" = "#2166AC", "Increasing" = "#B2182B", "No significant trend" = "grey60")
cap <- function(x, lim = 60) pmax(pmin(x, lim), -lim)   # keep a few huge CIs from flattening the panel

# fits that returned no slope (no year variation after filtering) have nothing to draw
F <- F[is.finite(pct_per_year)]
F[, lab := paste0(Analyte_Code, " | ", Cohort, " ", tolower(Panel), " (", first_year, "-", last_year, ")")]
F[, role := factor(role, unique(FAMILIES$role))]
# capped values stored as columns so the saved panel needs no helper functions
F[, `:=`(x_c = cap(pct_per_year), lo_c = cap(pct_lo), hi_c = cap(pct_hi))]
pA <- ggplot(F, aes(x_c, reorder(lab, pct_per_year), colour = direction)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey40") +
  geom_errorbarh(aes(xmin = lo_c, xmax = hi_c), height = 0, alpha = 0.6) +
  geom_point(size = 1.8) +
  facet_grid(family + role ~ ., scales = "free_y", space = "free_y", switch = "y",
             labeller = labeller(.multi_line = FALSE)) +
  scale_colour_manual(values = dir_pal) +
  labs(title = "Annual change, legacy and replacement chemicals",
       x = "% change per year", y = NULL, colour = NULL) +
  theme_minimal(base_size = 8) +
  theme(strip.placement = "outside", strip.text.y.left = element_text(angle = 0, hjust = 1),
        legend.position = "bottom", panel.spacing.y = unit(2, "pt"))

bc <- melt(by_class[, .(Category, decreasing, increasing, not_significant)], id.vars = "Category",
           variable.name = "direction", value.name = "n")
bc[, direction := factor(direction, c("decreasing", "not_significant", "increasing"),
                         c("Decreasing", "No significant trend", "Increasing"))]
pB <- ggplot(bc, aes(n, reorder(Category, n, sum), fill = direction)) +
  geom_col() +
  scale_fill_manual(values = dir_pal) +
  labs(title = "Trend direction by exposure class", x = "Chemical x cohort pairs", y = NULL, fill = NULL) +
  theme_minimal(base_size = 8) + theme(legend.position = "bottom")


# Panel titles/subtitles and in-plot labels removed 2026-11-03 (author: easier to read);
# panels carry A/B/C tags only.
fig <- pA + pB + plot_layout(widths = c(1.6, 1)) + plot_annotation(tag_levels = "A")
h <- max(7, 0.12 * nrow(F) + 2)
# ---- curated panel A for Figure S8 (2026-11-04) ---------------------------------
# The full panel A (~110 rows) is kept in the standalone figure and the tables.
# S8 shows only cohorts where a legacy chemical AND its replacement were both
# measured, and shows EVERY such row whatever its result (including BPS falling
# in DCHS adults, PFNA falling in SEARCH, no trend for MCOP/BPF), so the panel
# is not selected on outcome.
CUR <- rbindlist(list(
  data.table(cmp = "Phthalates: DEHP vs DEHTP/DINP", Cohort = "nuMoM2b",
             Analyte_Code = c("MEHP","MEHHP","MEOHP","MECPP","MECPTP","MEHHTP","MCOP"),
             role = c(rep("Legacy", 4), rep("Replacement", 3))),
  data.table(cmp = "Bisphenols: BPA vs BPS/BPF", Cohort = rep(c("TEDDY","DCHS"), each = 3),
             Analyte_Code = rep(c("BPA","BPS","BPF"), 2), role = rep(c("Legacy","Replacement","Replacement"), 2)),
  data.table(cmp = "PFAS: PFOS/PFOA vs PFNA/PFHxS", Cohort = rep(c("CTS","SEARCH"), each = 4),
             Analyte_Code = rep(c("PFOS","PFOA","PFNA","PFHXS"), 2), role = rep(c("Legacy","Legacy","Replacement","Replacement"), 2)),
  data.table(cmp = "Organophosphate vs pyrethroid", Cohort = "DCHS",
             Analyte_Code = c("TCP","PBA"), role = c("Legacy","Replacement"))))
Fc <- merge(M[is.finite(pct_per_year)], CUR, by = c("Cohort", "Analyte_Code"))
Fc[, `:=`(x_c = cap(pct_per_year), lo_c = cap(pct_lo), hi_c = cap(pct_hi))]
Fc[, lab := paste0(Analyte_Code, "  ", Cohort, " ", tolower(Panel), " ", first_year, "-", last_year)]
Fc[, cmp := factor(cmp, unique(CUR$cmp))]
setorder(Fc, cmp, Cohort, Panel, role, Analyte_Code)
Fc[, lab := factor(lab, rev(unique(lab)))]
fwrite(Fc[, .(cmp, Cohort, Panel, Analyte_Code, role, first_year, last_year,
              pct_per_year = round(pct_per_year, 1), pct_lo = round(pct_lo, 1), pct_hi = round(pct_hi, 1),
              p = signif(P_Value, 2), direction)], file.path(OUT, "S8A_curated_rows.csv"))
pA_cur <- ggplot(Fc, aes(x_c, lab, colour = direction, shape = role)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo_c, xmax = hi_c), width = 0, alpha = 0.6, orientation = "y") +
  geom_point(size = 2.6) +
  facet_grid(cmp ~ ., scales = "free_y", space = "free_y") +
  scale_colour_manual(values = dir_pal) +
  scale_shape_manual(values = c(Legacy = 16, Replacement = 17)) +
  labs(title = "Legacy and replacement chemicals measured in the same cohort",
       x = "% change per year", y = NULL, colour = NULL, shape = NULL) +
  theme_minimal(base_size = 10) +
  theme(strip.text.y = element_text(angle = 0, hjust = 0), legend.position = "bottom",
        legend.box = "vertical", panel.spacing.y = unit(6, "pt"))
saveRDS(list(pA = pA_cur, pB = pB), file.path(OUT, "S8_panels_AB.rds"))
ggsave(file.path(OUT, "SuppFig_Chemical_Year_Trends.pdf"), fig, width = 13, height = h, limitsize = FALSE)
ggsave(file.path(OUT, "SuppFig_Chemical_Year_Trends.png"), fig, width = 13, height = h, dpi = 200, limitsize = FALSE)

# ---- 6. console ------------------------------------------------------------------
cat("Cohort-panels with >=", MIN_SPAN, "sampling years:\n")
print(unique(M[, .(Cohort, Panel, first_year, last_year, n_years)])[order(Cohort)])
cat("\nDirection by class (tau = 0.5, p <", ALPHA, "):\n"); print(by_class[order(-chemical_cohort_pairs)])
cat("\nWithin one cohort and class, chemicals moving in opposite directions:\n"); print(opp)
cat("\nLegacy vs replacement:\n")
print(F[order(family, role, pct_per_year), .(family, role, Analyte_Code, Cohort, Panel, years = paste0(first_year, "-", last_year),
        pct_per_year = round(pct_per_year, 1), p = signif(P_Value, 2), direction)])
