# C5 -- Does the choice of reference population change who is "highly exposed"?
# (Round 2: Reviewer 1, option 1 -- "show a concrete example where
#  population-aware percentiles change a clinically or epidemiologically
#  meaningful inference")
#
# The inference examined is the standard biomonitoring one: flagging the
# participants above the 90th percentile of a reference distribution as the
# highly exposed group (for follow-up, intervention, or as a binary exposure in
# an association analysis). Two references are compared, per analyte, within a
# single sample matrix and population panel (child / adult):
#   pooled  - the 90th percentile of all HHEAR cohorts combined
#   own     - the 90th percentile of the participant's own cohort
# By construction "own" flags ~10% of every cohort. The question is what the
# pooled reference does instead, and how far the two flagged groups agree.
#
# Comparability rules (see Coding/Paper1_Submission2/harmonization_helpers.R):
#   * one Sample_Matrix per comparison; units converted with to_ng_per_mL();
#     rows in pass-through units (ng/g, pmol/mL, nM, unknown) are dropped,
#     because those are not on a common scale across cohorts
#   * PRIMARY set = serum/plasma PFAS: same matrix, same units, no urinary
#     dilution correction needed -- the cleanest cross-cohort comparison HHEAR has
#   * SECONDARY set = urinary analytes; reported with the caveat that HHEAR
#     urine values here are not dilution-corrected
#   * non-detects enter at the LOD (Exposure_Obs); a P90 cutoff sits above the
#     LOD for every analyte kept (checked below)
#
# Read-only on Input/. Writes only to Output/Reviewer_Round2/C5/.
# Run from the project root:  Rscript Coding/Paper1_Submission2/Reviewer_Round2/C5_Reference_Changes_Inference.R

suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(ggplot2); library(patchwork)
})
source("Coding/Paper1_Submission2/harmonization_helpers.R")
source("Coding/Paper1_Submission2/censored_quantile_helpers.R")

HARM <- "Input/Harmonized_Datasets"
OUT  <- "Output/Reviewer_Round2/C5"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

load(file.path(HARM, "Harmonized_Targeted.RData"))          # Targeted_Data_All
CW <- fread(file.path(HARM, "Targeted_Analyte_Category_Mapping.csv"))
e  <- new.env(); load(file.path(HARM, "Harmonized_SDOH.RData"), envir = e)
S  <- as.data.table(e$SDOH_Data_All)

CLASSES_11 <- c("FlameRetardants","PAH","Parabens","PCB","Pesticides","PFAS",
                "Phenols","Phthalates","Smoking","TraceElements","VOCs")
CONVERTIBLE <- c("ng/mL", "ug/dL", "ng/L", "pg/mL", "pg/ml", "mg/L")
P_CUT   <- 0.90
MIN_COH <- 3      # cohorts per analyte x matrix x panel
MIN_N   <- 30     # participants per cohort

T <- as.data.table(Targeted_Data_All)
T <- T[Analyte_Code %in% CW[Category %in% CLASSES_11, Code]]
T[is.na(Units), Units := "ng/mL"]
n_dropped_units <- T[!Units %in% CONVERTIBLE, .N]
T <- T[Units %in% CONVERTIBLE]

# Build the censored panel separately per sample matrix so a participant with
# both serum and DBS values (ATLAA) is never averaged across media.
panel_for <- function(id_var) {
  rbindlist(lapply(split(T[!is.na(get(id_var))], by = "Sample_Matrix"), function(d) {
    p <- as.data.table(build_censored_targeted_panel(as.data.frame(d), id_var, to_ng_per_mL))
    p[, Sample_Matrix := d$Sample_Matrix[1]]
    setnames(p, id_var, "pid")
  }))
}
P <- rbind(panel_for("Child_PID")[, Panel := "Child"],
           panel_for("Mat_PID")[,   Panel := "Adult"])
P[, Cohort := std_study(Study)]
P <- P[!is.na(Exposure_Obs)]
P <- merge(P, CW[, .(Analyte_Code = Code, Category)], by = "Analyte_Code")

# Sampling year (cohort-level median) -- used to describe what the pooled flag
# is really picking up.
S[, Cohort := std_study(Study)]
S[, pid := fcoalesce(as.character(Child_PID), as.character(Mat_PID))]
yr <- unique(T[!is.na(Sample_Collection_Year),
               .(Cohort = std_study(Study),
                 pid = fcoalesce(as.character(Child_PID), as.character(Mat_PID)),
                 Sample_Collection_Year)])
coh_year <- yr[, .(median_year = median(Sample_Collection_Year)), by = Cohort]

# ---- eligible comparisons -------------------------------------------------------
P[, n_coh := .N, by = .(Analyte_Code, Sample_Matrix, Panel, Cohort)]
P <- P[n_coh >= MIN_N]
P[, k := uniqueN(Cohort), by = .(Analyte_Code, Sample_Matrix, Panel)]
P <- P[k >= MIN_COH]

P[, cut_pooled := quantile(Exposure_Obs, P_CUT, type = 7), by = .(Analyte_Code, Sample_Matrix, Panel)]
P[, cut_own    := quantile(Exposure_Obs, P_CUT, type = 7), by = .(Analyte_Code, Sample_Matrix, Panel, Cohort)]
# ">" not ">=": ties at a cutoff (mostly non-detects sitting at the LOD) are not flagged
P[, high_pooled := Exposure_Obs > cut_pooled]
P[, high_own    := Exposure_Obs > cut_own]

# a pooled P90 that sits at or below the LOD would make "high" mean "detected"
lod_check <- P[, .(cut_pooled = cut_pooled[1], max_LOD = suppressWarnings(max(LOD, na.rm = TRUE)),
                   pct_nd = mean(!Detected %in% TRUE)), by = .(Analyte_Code, Sample_Matrix, Panel)]
bad_lod <- lod_check[is.finite(max_LOD) & cut_pooled <= max_LOD]
P <- P[!bad_lod, on = .(Analyte_Code, Sample_Matrix, Panel)]

kappa2 <- function(a, b) {
  po <- mean(a == b); pa <- mean(a); pb <- mean(b)
  pe <- pa * pb + (1 - pa) * (1 - pb)
  if (pe == 1) NA_real_ else (po - pe) / (1 - pe)
}

# ---- per comparison --------------------------------------------------------------
per_cmp <- P[, {
  by_coh <- .SD[, .(pct_flag_pooled = 100 * mean(high_pooled), n = .N), by = Cohort]
  top    <- .SD[high_pooled == TRUE, .N, by = Cohort][order(-N)]
  list(
    Category          = Category[1],
    n_cohorts         = uniqueN(Cohort),
    n_participants    = .N,
    kappa             = kappa2(high_pooled, high_own),
    pct_own_high_missed_by_pooled = 100 * mean(!high_pooled[high_own]),
    top_cohort        = top$Cohort[1],
    top_cohort_share_of_pooled_high = 100 * top$N[1] / sum(top$N),
    cohorts_le1pct_flagged = sum(by_coh$pct_flag_pooled <= 1),
    cohorts_ge30pct_flagged = sum(by_coh$pct_flag_pooled >= 30),
    max_cohort_pct_flagged = max(by_coh$pct_flag_pooled)
  )
}, by = .(Analyte_Code, Sample_Matrix, Panel)]
per_cmp[, Set := fifelse(Category == "PFAS" & Sample_Matrix %in% c("Serum", "Plasma"),
                         "Primary: serum/plasma PFAS",
                         fifelse(Sample_Matrix == "Urine", "Secondary: urine", "Other matrix"))]
fwrite(per_cmp[order(Set, kappa)], file.path(OUT, "reference_comparison_by_analyte.csv"))

per_cohort <- P[, .(n = .N, pct_flag_pooled = round(100 * mean(high_pooled), 1),
                    pct_flag_own = round(100 * mean(high_own), 1),
                    median_conc = signif(median(Exposure_Obs), 3)),
                by = .(Analyte_Code, Sample_Matrix, Panel, Cohort)]
per_cohort <- merge(per_cohort, coh_year, by = "Cohort", all.x = TRUE)
fwrite(per_cohort[order(Analyte_Code, Sample_Matrix, Panel, -pct_flag_pooled)],
       file.path(OUT, "reference_comparison_by_cohort.csv"))

summ <- per_cmp[Set != "Other matrix", .(
  comparisons = .N,
  median_kappa = round(median(kappa, na.rm = TRUE), 2),
  median_pct_own_high_missed = round(median(pct_own_high_missed_by_pooled), 1),
  median_top_cohort_share = round(median(top_cohort_share_of_pooled_high), 1),
  pct_cmp_with_a_cohort_le1pct = round(100 * mean(cohorts_le1pct_flagged >= 1), 1),
  pct_cmp_with_a_cohort_ge30pct = round(100 * mean(cohorts_ge30pct_flagged >= 1), 1)),
  by = Set]
fwrite(summ, file.path(OUT, "reference_comparison_summary.csv"))

# ---- figure: the PFOA example plus the distribution over all comparisons ---------
ex <- per_cohort[Analyte_Code %in% c("PFOA", "PFOS") & Sample_Matrix %in% c("Serum", "Plasma")]
ex[, lab := paste0(Cohort, fifelse(is.na(median_year), " (year n/r)", paste0(" (", round(median_year), ")")))]
ex[, facet := paste0(Analyte_Code, ", ", tolower(Panel))]
pA <- ggplot(ex, aes(pct_flag_pooled, reorder(paste(facet, lab, sep = "___"), pct_flag_pooled))) +
  geom_vline(xintercept = 100 * (1 - P_CUT), linetype = 2, colour = "grey40") +
  geom_col(fill = "#4C72B0") +
  facet_wrap(~ facet, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(x) sub("^.*___", "", x)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05))) +
  labs(title = "Serum PFAS above the pooled 90th percentile",
       x = "% flagged (pooled P90)", y = NULL) +
  theme_minimal(base_size = 8)

pB <- ggplot(per_cmp[Set != "Other matrix"], aes(kappa, fill = Set)) +
  geom_histogram(binwidth = 0.05, boundary = 0, colour = "white") +
  facet_wrap(~ Set, ncol = 1, scales = "free_y") +
  scale_fill_manual(values = c("#4C72B0", "#DD8452"), guide = "none") +
  labs(title = "Pooled vs cohort-specific classification",
       x = "Cohen's kappa", y = "Comparisons") +
  theme_minimal(base_size = 8)


# Panel titles/subtitles and in-plot labels removed 2026-11-03 (author: easier to read);
# panels carry A/B/C tags only.
# panel A is reused in the combined Figure S17 (R1_2_S17_assemble.R)
saveRDS(list(pA = pA), file.path(OUT, "S17_panel_reference.rds"))
fig <- pA + pB + plot_layout(widths = c(1.2, 1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT, "Fig_Reference_Changes_HighExposure.pdf"), fig, width = 10, height = 4.5)
ggsave(file.path(OUT, "Fig_Reference_Changes_HighExposure.png"), fig, width = 10, height = 4.5, dpi = 300)

# ---- console ---------------------------------------------------------------------
cat("Rows dropped for non-convertible units:", n_dropped_units, "\n")
cat("Comparisons dropped because pooled P90 <= LOD:", nrow(bad_lod), "\n\n")
print(summ)
cat("\nPrimary set, per analyte:\n")
print(per_cmp[Set == "Primary: serum/plasma PFAS",
              .(Analyte_Code, Sample_Matrix, Panel, n_cohorts, n_participants, kappa = round(kappa, 2),
                pct_own_high_missed_by_pooled = round(pct_own_high_missed_by_pooled, 1),
                top_cohort, top_share = round(top_cohort_share_of_pooled_high, 1),
                cohorts_le1pct_flagged, max_cohort_pct_flagged = round(max_cohort_pct_flagged, 1))])
cat("\nPFOA example:\n"); print(ex[order(-pct_flag_pooled)])
