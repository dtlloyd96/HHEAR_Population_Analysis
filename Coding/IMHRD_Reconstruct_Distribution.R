## =====================================================================
## IMHRD PFAS — reconstruct a SYNTHETIC participant-level distribution
## from the published year x analyte summary statistics.
##
## IMHRD ships ONLY aggregated tables (median / IQR / LOQ / n per
## sample_collection_year x analyte); there are no participant rows. This
## script reconstructs a *plausible* per-participant distribution so the
## cohort can be explored/visualized alongside the others. It is SYNTHETIC
## reference data (Study = "IMHRD_reconstructed") and must NOT be merged into
## Harmonized_Targeted or treated as real measurements.
##
## Assumption: each year x analyte is lognormal (the usual model for
## environmental biomarker concentrations). A lognormal is uniquely pinned by
## its median and IQR:
##   median = exp(mu)                          => mu = log(median)
##   Q1 = median * exp(-z*sigma),  Q3 = median * exp(+z*sigma),  z = qnorm(0.75)
##   IQR = Q3 - Q1 = 2*median*sinh(z*sigma)    => sigma = asinh(IQR/(2*median))/z
## We then place n values on a deterministic lognormal quantile lattice (so the
## target median/IQR are reproduced even for small n), and (to mirror the
## LOD/sqrt(2) convention used by the other cohorts) floor anything below LOQ at
## LOQ/sqrt(2). These are REPRESENTATIVE points, not random draws / real people.
##
## Charleston, SC cohort. Output: Input/Reference/IMHRD_PFAS_Reconstructed.RData
## + a CSV, both under Input/Reference/ (outside the *All.RData auto-load glob).
## =====================================================================
suppressPackageStartupMessages({library(tidyverse)})
set.seed(20260803)  # reproducible

setwd_try <- function(p) tryCatch({setwd(p); TRUE}, error=function(e) FALSE)
if (!file.exists("../Input/Reference/IMHRD_PFAS_YearSummary.RData"))
  setwd_try(file.path(Sys.getenv("HHEAR_ROOT", ".."), "Coding"))

load("../Input/Reference/IMHRD_PFAS_YearSummary.RData")  # -> IMHRD_PFAS_YearSummary (list)

year_col <- names(IMHRD_PFAS_YearSummary$median)[1]   # "sample_collection_year"
to_long <- function(df, value_name) df %>%
  pivot_longer(-all_of(year_col), names_to = "Analyte_Code", values_to = value_name) %>%
  mutate(!!value_name := suppressWarnings(as.numeric(.data[[value_name]])))

summ <- to_long(IMHRD_PFAS_YearSummary$median, "median") %>%
  left_join(to_long(IMHRD_PFAS_YearSummary$iqr,    "iqr"),    by = c(year_col,"Analyte_Code")) %>%
  left_join(to_long(IMHRD_PFAS_YearSummary$loq,    "loq"),    by = c(year_col,"Analyte_Code")) %>%
  left_join(to_long(IMHRD_PFAS_YearSummary$n,      "n"),      by = c(year_col,"Analyte_Code")) %>%
  rename(Sample_Collection_Year = !!year_col) %>%
  filter(!is.na(n), n > 0, !is.na(median), median > 0)

z <- qnorm(0.75)  # 0.6744898

simulate_cell <- function(med, iqr, loq, n) {
  n <- as.integer(round(n))
  if (is.na(iqr) || iqr <= 0) {                      # no spread -> replicate median
    vals <- rep(med, n)
  } else {
    sigma <- asinh(iqr / (2 * med)) / z
    mu    <- log(med)
    # Deterministic quantile lattice: place the n points at the theoretical
    # lognormal quantiles (i-0.5)/n. This reproduces the target median/IQR
    # closely for ANY n (random draws are too noisy for the small-n cells).
    probs <- (seq_len(n) - 0.5) / n
    vals  <- qlnorm(probs, meanlog = mu, sdlog = sigma)
  }
  if (!is.na(loq)) vals <- ifelse(vals < loq, loq / sqrt(2), vals)  # LOD/sqrt(2)-style floor
  vals
}

## IMPORTANT — MARGINALS ONLY, NO JOINT STRUCTURE.
## The source provides only per-(year, analyte) summary statistics, so nothing is
## known about how analytes co-vary within a person. An earlier version emitted an
## ascending lattice bound to a SHARED sequential participant id, which made every
## analyte perfectly rank-correlated (Spearman = 1.000) — entirely an artifact.
## We therefore (a) shuffle the lattice within each cell and (b) give each
## (year, analyte) cell its OWN row ids rather than a shared person axis, so no
## cross-analyte linkage is implied. Use this for marginal/distributional
## comparisons only; do NOT use it for correlation, PCA, network or co-exposure
## analyses.
IMHRD_Reconstructed <- summ %>%
  mutate(.vals = pmap(list(median, iqr, loq, n), simulate_cell)) %>%
  transmute(
    Study          = "IMHRD_reconstructed",
    Chemical_Group = "SPFAS",
    Analyte_Code,
    Row_ID         = pmap(list(Sample_Collection_Year, Analyte_Code, n),
                          ~ paste0("IMHRD_", ..1, "_", ..2, "_", seq_len(round(..3)))),
    Concentration  = map(.vals, ~ sample(.x)),   # shuffle within cell (set.seed above)
    Units          = "ng/mL",     # IMHRD PFAS reported in ng/mL (serum)
    LOD            = loq,         # LOQ used as the reporting limit
    Sample_Collection_Year,
    Reconstructed  = TRUE
  ) %>%
  unnest(c(Row_ID, Concentration)) %>%
  mutate(Detected = Concentration >= LOD)

## ---- QC: do reconstructed summaries match the targets? ----
qc <- IMHRD_Reconstructed %>%
  group_by(Sample_Collection_Year, Analyte_Code) %>%
  summarise(median_sim = median(Concentration),
            iqr_sim = IQR(Concentration), n_sim = n(), .groups = "drop") %>%
  left_join(summ, by = c("Sample_Collection_Year","Analyte_Code")) %>%
  mutate(median_pct_err = 100*(median_sim-median)/median,
         iqr_pct_err = ifelse(iqr>0, 100*(iqr_sim-iqr)/iqr, NA))

dir.create("../Input/Reference", recursive = TRUE, showWarnings = FALSE)
save(IMHRD_Reconstructed, summ, qc,
     file = "../Input/Reference/IMHRD_PFAS_Reconstructed.RData")
readr::write_csv(IMHRD_Reconstructed, "../Input/Reference/IMHRD_PFAS_Reconstructed.csv")

cat("Reconstructed rows:", nrow(IMHRD_Reconstructed),
    " | year x analyte cells:", nrow(summ),
    " | analytes:", n_distinct(IMHRD_Reconstructed$Analyte_Code),
    " | years:", n_distinct(IMHRD_Reconstructed$Sample_Collection_Year), "\n")
cat("QC — median reconstruction error: median|%|=",
    round(median(abs(qc$median_pct_err), na.rm=TRUE),2),
    " max|%|=", round(max(abs(qc$median_pct_err), na.rm=TRUE),2), "\n")
cat("QC — IQR reconstruction error:    median|%|=",
    round(median(abs(qc$iqr_pct_err), na.rm=TRUE),2),
    " max|%|=", round(max(abs(qc$iqr_pct_err), na.rm=TRUE),2), "\n")
## The IQR is systematically UNDER-reproduced because the LOQ/sqrt(2) floor
## truncates the lower tail after the lattice is built (most severe in cells where
## a large share of values sit below the LOQ). Warn rather than fail silently.
if (isTRUE(median(abs(qc$iqr_pct_err), na.rm = TRUE) > 5)) {
  warning("IMHRD reconstruction: median |IQR error| > 5% (LOQ flooring compresses ",
          "the lower tail). Reconstructed spread is conservative — treat IQR-based ",
          "comparisons with care.")
}
cat("NOTE: marginal distributions only — cross-analyte correlations are NOT meaningful.\n")
cat("Saved: Input/Reference/IMHRD_PFAS_Reconstructed.RData (+ .csv)\n")
