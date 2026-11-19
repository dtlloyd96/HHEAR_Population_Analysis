# Figure 3 -- sensitivity of each exposure-outcome association to how cohort is
# handled. For the 44 maternal pairs of Figure 3B, fits on the same participants:
#   naive    Y ~ X + C                                  (cohort ignored)
#   fe       Y ~ X + C + cohort                         (cohort fixed effects)
#   meta     Y ~ X + C within each cohort, combined by random-effects meta-analysis
#   bw, bb   Y ~ Bw (X - Xbar_j) + Bb Xbar_j + C        (within-between decomposition)
# Without covariates bw is algebraically fe and bb is the n-weighted slope through
# the cohort means, so both are asserted against cohort_regime_pairs.csv.
#
# Specs: C = none (as Figure 3B) or maternal age (the only covariate every cohort
# has). CTS has no BMI since 2026-11-09 (it was coded 0/1/2/9, not kg/m2).
# Same prep as cohort_regime_diagnostic.R: ng/mL, log1p, mean per study x person.
# Effects share one per-pair scale (within-cohort SDs of X and Y), so differences
# between estimators are differences in slope alone.
# Writes Output/Reviewer_Round2/Fig3_Estimators/. Run from the project root.
suppressPackageStartupMessages({ library(data.table); library(metafor); library(arrow) })
source("Coding/Paper1_Submission2/harmonization_helpers.R")
OUT <- "Output/Reviewer_Round2/Fig3_Estimators"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
IN  <- "Health_Outcome_Analysis/Input"
MIN_STUDY_N <- 20

# ---- the 44 pairs: Figure 3 panel A slice joined to the decomposition -----------
f <- list.files("Health_Outcome_Analysis/Output_Targeted_ExWAS_R", "exposure_results_table\\.parquet$",
                recursive = TRUE, full.names = TRUE)
f <- f[grepl("scale=log1p", f) & grepl("panel=mat", f)]
ex <- rbindlist(lapply(f, read_parquet), fill = TRUE)
ex <- ex[status == "ok" & model == "Adjusted" & panel == "mat" & !outcome %in% c("Asthma", "Severe_Asthma") &
           is.finite(beta) & is.finite(se) & se > 0]
A <- unique(ex[, .(k = uniqueN(study)), by = .(outcome, exposure)][k >= 2, .(outcome, analyte = exposure)])
regime <- fread("Output/Paper1_Output_Submission2/Cohort_Regime_Diagnostic/cohort_regime_pairs.csv")[panel == "mat"]
regime <- regime[A, on = .(outcome, analyte), nomatch = NULL]
stopifnot(nrow(regime) == 44)

# ---- inputs -----------------------------------------------------------------------
tg <- fread(file.path(IN, "Harmonized_Targeted.csv"), encoding = "Latin-1", colClasses = list(character = "Mat_PID"))
tg <- tg[Analyte_Code %in% regime$analyte & !is.na(Mat_PID) & Mat_PID != ""]
tg[, `:=`(Study = std_study(Study), conc = to_ng_per_mL(Concentration, Units))]
X <- tg[!is.na(conc) & conc >= 0, .(x = mean(log1p(conc))), by = .(Study, PID = Mat_PID, analyte = Analyte_Code)]

oc <- fread(file.path(IN, "Harmonized_Outcomes.csv"), colClasses = "character",
            select = c("Study", "Mat_PID", unique(regime$outcome)))
oc <- oc[!is.na(Mat_PID) & Mat_PID != ""]
oc <- melt(oc[, Study := std_study(Study)], id.vars = c("Study", "Mat_PID"), variable.name = "outcome",
           value.name = "y", variable.factor = FALSE)
oc[, y := suppressWarnings(as.numeric(y))]
oc <- oc[!is.na(y)]
setnames(oc, "Mat_PID", "PID")

e <- new.env(); load("Input/Harmonized_Datasets/Harmonized_SDOH.RData", envir = e)
S <- as.data.table(e$SDOH_Data_All)[!is.na(Mat_PID) & Mat_PID != "",
       .(Study = std_study(Study), PID = as.character(Mat_PID), age = suppressWarnings(as.numeric(Mat_Age)))]
S <- S[!is.na(age) & age >= 12 & age <= 60, .(age = median(age)), by = .(Study, PID)]

D <- merge(X, oc, by = c("Study", "PID"), allow.cartesian = TRUE)
D <- merge(D, S, by = c("Study", "PID"), all.x = TRUE)
D <- D[paste(analyte, outcome) %chin% regime[, paste(analyte, outcome)]]

# ---- estimators for one pair ------------------------------------------------------
coef_row <- function(fit, term, name) {
  s <- summary(fit)$coefficients
  if (!term %in% rownames(s)) return(NULL)
  data.table(estimator = name, est = s[term, 1], se = s[term, 2])
}

fit_pair <- function(d, adjust) {
  if (adjust) d <- d[!is.na(age)]
  d <- d[, if (.N >= MIN_STUDY_N && var(x) > 0) .SD, by = Study]
  K <- uniqueN(d$Study)
  if (K < 2) return(NULL)
  cv <- if (adjust) " + age" else ""
  d[, `:=`(xbar = mean(x)), by = Study]
  d[, xdev := x - xbar]
  d[, Study := factor(Study)]
  # one scale for every estimator in this pair
  sx <- sqrt(sum(d$xdev^2) / (nrow(d) - K))
  sy <- d[, y - mean(y), by = Study][, sqrt(sum(V1^2) / (nrow(d) - K))]

  fn <- lm(as.formula(paste0("y ~ x", cv)), d)
  ff <- lm(as.formula(paste0("y ~ x", cv, " + Study")), d)
  fw <- lm(as.formula(paste0("y ~ xdev + xbar", cv)), d)
  per <- d[, { f <- lm(as.formula(paste0("y ~ x", cv)), .SD); s <- summary(f)$coefficients
               .(b = s["x", 1], s = s["x", 2], n = .N) }, by = Study]
  ma <- tryCatch(rma(yi = per$b, sei = per$s, method = "REML", control = list(stepadj = 0.5, maxiter = 1000)),
                 error = function(e) rma(yi = per$b, sei = per$s, method = "DL"))
  r <- rbind(coef_row(fn, "x", "naive"), coef_row(ff, "x", "fe"),
             data.table(estimator = "meta", est = as.numeric(ma$beta), se = ma$se),
             coef_row(fw, "xdev", "bw"), coef_row(fw, "xbar", "bb"))
  r[, `:=`(K = K, N = nrow(d), std = sx / sy, I2 = ma$I2, QEp = ma$QEp,
           n_pos = sum(per$b > 0), n_neg = sum(per$b < 0))]
  r
}

stopifnot(D[Study == "CTS" & outcome == "BMI", .N] == 0)
res <- rbindlist(lapply(c(FALSE, TRUE), function(adj) {
  r <- D[, fit_pair(copy(.SD), adj), by = .(analyte, outcome)]
  r[, covariates := if (adj) "maternal age" else "none"]
}))
res[, `:=`(lo = est - 1.96 * se, hi = est + 1.96 * se, r = est * std, r_lo = (est - 1.96 * se) * std,
           r_hi = (est + 1.96 * se) * std, sig = abs(est / se) > 1.96)]

# ---- reproduce the existing decomposition (no covariates) -------------------------
chk <- dcast(res[covariates == "none"], analyte + outcome + K + N ~ estimator, value.var = "est")
chk <- merge(chk, regime[, .(analyte, outcome, K0 = K, N0 = N, b_within, b_naive, b_between)], by = c("analyte", "outcome"))
stopifnot(nrow(chk) == 44, all(chk$K == chk$K0), all(chk$N == chk$N0),
          max(abs(chk$fe - chk$b_within)) < 1e-6, max(abs(chk$bw - chk$b_within)) < 1e-6,
          max(abs(chk$naive - chk$b_naive)) < 1e-6, max(abs(chk$bb - chk$b_between), na.rm = TRUE) < 1e-6)
cat("Reproduced cohort_regime_pairs.csv: fe = bw = b_within, naive, bb = b_between (44 pairs)\n")

fwrite(res, file.path(OUT, "estimates_long.csv"))

# ---- does the conclusion depend on the estimator? ---------------------------------
summ <- function(x) {
  w <- dcast(x, analyte + outcome + K ~ estimator, value.var = c("est", "sig", "lo", "hi"))
  rev_fn <- sign(w$est_fe) != sign(w$est_naive)
  rev_mn <- sign(w$est_meta) != sign(w$est_naive)
  ratio  <- w$est_naive / w$est_fe
  data.table(pairs = nrow(w), K2 = sum(w$K == 2),
    fe_vs_naive_reversed = sum(rev_fn), meta_vs_naive_reversed = sum(rev_mn),
    fe_vs_meta_reversed = sum(sign(w$est_fe) != sign(w$est_meta)),
    magnitude_only = sum(!rev_fn & (ratio < 0.5 | ratio > 2)),
    concordant = sum(!rev_fn & ratio >= 0.5 & ratio <= 2),
    naive_sig = sum(w$sig_naive), fe_sig = sum(w$sig_fe), meta_sig = sum(w$sig_meta),
    sig_differs_naive_fe = sum(w$sig_naive != w$sig_fe),
    sig_differs_naive_meta = sum(w$sig_naive != w$sig_meta),
    reversed_and_either_sig = sum(rev_fn & (w$sig_naive | w$sig_fe)),
    reversed_and_both_sig = sum(rev_fn & w$sig_naive & w$sig_fe),
    bb_opposes_bw = sum(sign(w$est_bb) != sign(w$est_bw), na.rm = TRUE),
    reversed_K_ge3 = sum(rev_fn & w$K >= 3), pairs_K_ge3 = sum(w$K >= 3))
}
S_tab <- res[, summ(.SD), by = covariates]
fwrite(S_tab, file.path(OUT, "summary_by_spec.csv"))
options(width = 220)
print(t(S_tab))

w <- dcast(res[covariates == "none"], analyte + outcome + K + N ~ estimator, value.var = c("r", "sig"))
fwrite(w, file.path(OUT, "wide_unadjusted.csv"))
cat("\nPairs whose sign or significance changes by estimator (no covariates):\n")
print(w[sign(r_naive) != sign(r_fe) | sig_naive != sig_fe | sig_fe != sig_meta,
        .(analyte, outcome, K, N, naive = round(r_naive, 3), fe = round(r_fe, 3), meta = round(r_meta, 3),
          bb = round(r_bb, 3), sig_naive, sig_fe, sig_meta)][order(outcome, analyte)])
