## ============================================================================
## cohort_regime_diagnostic.R
##
## Purpose
## -------
## Apply the cohort-pooling simulation's diagnostic to the REAL targeted ExWAS
## data: for every continuous outcome x targeted-analyte pair that appears in
## >= MIN_STUDIES studies, decompose the unadjusted exposure->outcome
## association into the three estimands from the simulation and classify which
## "regime" the pair is in:
##
##   within-pooled  : common within-study slope (cohort-adjusted / fixed-effect)
##                    -- the contextual effect, the truth-analog
##   naive-pooled   : slope ignoring the Study label -- what pooling reports
##   between-study  : ecological slope of study-mean(y) on study-mean(x)
##
## Classification per pair (on log1p exposure scale, unadjusted):
##   concordant        : naive ~ within (same sign, ratio in [0.5, 2])
##   magnitude_distort : same sign, naive materially larger/smaller
##   sign_reversal     : naive and within have OPPOSITE signs (Simpson-like)
## plus a heterogeneity measure (I^2) and the within-vs-between sign agreement.
##
## Everything is computed from per-study sufficient statistics (n, Sx, Sy, Sxx,
## Sxy, Syy), so no per-pair lm() is fit -- the full grid is feasible.
##
## Mirrors the data prep in Targeted_ExWAS_Replication.Rmd (std_study, IL marker
## drop, unit -> ng/mL, log1p, mean per study-PID-analyte).
##
## Output: Output/Paper1_Output_Submission2/Cohort_Regime_Diagnostic/
##   cohort_regime_pairs.csv     - one row per (panel, outcome, analyte) pair
##   cohort_regime_summary.csv   - class counts / fractions overall & by panel
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(stringr)
  library(tibble)
})

## ---- Config -----------------------------------------------------------------
## NOTE: empirically almost every (analyte, outcome) pair is observed in only
## ONE cohort at n>=20, so we decompose at >=2 cohorts to capture the usable set
## and separately flag the >=3-cohort subset. The coverage table records how
## many cohorts each pair is seen in (the sparsity itself is a key result).
MIN_STUDIES   <- 2     # pair must appear in at least this many studies to decompose
MIN_STUDY_N   <- 20    # a study contributes only if it has >= this many cc rows
RATIO_LO      <- 0.5   # concordance band on naive/within slope ratio
RATIO_HI      <- 2.0
SIG_T         <- 1.96  # |t| threshold for "significant" (strong reversal flag)
MAX_ANALYTES  <- Inf   # set small (e.g. 8) for a quick smoke test

## ---- Permutation null (added 2026-08-11) ------------------------------------
## CODE_AUDIT_2026-08-11.md section 1.6.
##
## `class` compares sign(b_within) against sign(b_naive) with no precision
## condition. Under a true null those signs are close to independent, so a large
## "sign_reversal" rate is expected FROM NOISE ALONE -- especially at K = 2,
## which is most of the grid. Reporting the raw rate without a null benchmark
## makes it uninterpretable.
##
## Null: permute the cohort label within each pair. This destroys any real
## cohort-level structure (so between-cohort and within-cohort slopes agree in
## expectation) while preserving N, K, the cohort size profile, and the marginal
## x-y association. The fraction of permutations landing in each class is the
## rate noise produces for that pair's exact geometry.
##
## Set N_PERM <- 0 to skip (the sweep is roughly N_PERM x slower with it on).
## N_PERM = 50 for routine reruns (2026-08-11). At 50 the null rate for a single
## pair has a standard error of ~7pp, which is too coarse to quote per-pair, but
## the summary tables average the null across hundreds of pairs, so the reported
## null_pct_* and excess_pct_* are stable to well under a point. Raise to 200+
## if you ever want to quote a per-pair null.
N_PERM      <- 50
PERM_SEED   <- 20260811L

## ---- Constants copied from the ExWAS pipeline -------------------------------
IL_MARKERS <- c(
  "IL10","IL17A","IL1B","IL5","IL6","IL8","IL12P70","IL13","IL15","IL17EIL25",
  "IL17F","IL2","IL21","IL22","IL23","IL27","IL28A","IL31","IL33","IL4","IL9",
  "MIP3A","MMP1","MMP10","MMP2","MMP7","MMP9","TIMP1","TIMP2","TCHL","TNON",
  "CHOL","TG","HAPTOGLOBIN","FLT3L","GMCSF","nan","TNFB"
)

# std_study() and to_ng_per_mL() now come from harmonization_helpers.R
# (2026-08-11). The local to_ng_per_mL passed mg/dL, ug/mL, ng/g, pmol/mL, nM
# and kg/m3 through unconverted as if they were ng/mL -- 5.9% of rows, mis-scaled
# by 10^3-10^6. Since this script produces the sign_reversal classification behind
# Figure 3 Panel B, that error fed the headline result directly.
local({
  cands <- c("harmonization_helpers.R",
             "Coding/Paper1_Submission2/harmonization_helpers.R",
             "../harmonization_helpers.R")
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) stop("Cannot locate harmonization_helpers.R")
  source(hit[[1]], local = FALSE)
})

## Outcome whitelist ----------------------------------------------------------
## Previously ANY column with >2 distinct numeric values was treated as a
## continuous outcome. Against Harmonized_Outcomes.csv that swept in 193 of 345
## columns, 64 of them with <=5 distinct values -- questionnaire response codes
## such as chapstick_month, mouthwash_last, deodorant_last, polish_remover_last,
## cologne_last and birth_due_date (values 1-3). All were OLS-regressed and
## counted in the concordant / magnitude_distort / sign_reversal percentages.
##
## Set OUTCOME_WHITELIST to NULL to restore the old cardinality-only behaviour
## (not recommended); MIN_DISTINCT_CONTINUOUS is a floor applied in either case.
##
## !! THIS LIST IS A SCIENTIFIC SCOPING DECISION AND NEEDS AUTHOR SIGN-OFF.
## !! It is the set of harmonized, cross-cohort, genuinely continuous clinical
## !! outcomes present in Harmonized_Outcomes.csv. Add or remove deliberately;
## !! do not go back to inferring outcome status from cardinality.
OUTCOME_WHITELIST <- c(
  # birth
  "Gest_Age_Birth", "Birthweight", "Birth_Weight_Z", "APGAR_1", "APGAR_5",
  # anthropometry
  "BMI", "BMI_Z", "Waist_Circumference",
  # cardiometabolic
  "Systolic_BP", "Diastolic_BP", "Cholesterol_HDL", "Cholesterol_LDL",
  "Cholesterol_Total", "Triglycerides", "Glucose", "Insulin",
  "Insulin_Resistance", "HBA1C", "HSCRP", "ALT",
  # neurodevelopment
  "EarlyLearning_Composite", "Language_Tscore", "Fine_Motor_Tscore",
  "Visual_Tscore", "Cognitive_Composite", "Motor_Composite", "Social_Composite"
)
MIN_DISTINCT_CONTINUOUS <- 10   # a "continuous" outcome needs real spread

## ---- Paths ------------------------------------------------------------------
resolve_existing_dir <- function(candidates) {
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) == 0) stop("Could not find a candidate dir:\n",
                             paste(candidates, collapse = "\n"))
  hit[[1]]
}

input_dir <- resolve_existing_dir(c(
  "../../Health_Outcome_Analysis/Input",
  "../../../Health_Outcome_Analysis/Input",
  "Health_Outcome_Analysis/Input"
))
submission_out <- resolve_existing_dir(c(
  "../../Output/Paper1_Output_Submission2",
  "../../../Output/Paper1_Output_Submission2",
  "Output/Paper1_Output_Submission2"
))
out_dir <- file.path(submission_out, "Cohort_Regime_Diagnostic")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Load & prep ------------------------------------------------------------
message("Reading inputs ...")
targeted <- readr::read_csv(file.path(input_dir, "Harmonized_Targeted.csv"),
                            show_col_types = FALSE,
                            locale = readr::locale(encoding = "latin1"))
outcomes <- readr::read_csv(file.path(input_dir, "Harmonized_Outcomes.csv"),
                            show_col_types = FALSE, guess_max = 100000)

targeted$Study <- std_study(targeted$Study)
outcomes$Study <- std_study(outcomes$Study)
targeted <- targeted %>% filter(!(Analyte_Code %in% IL_MARKERS))

unit_col <- c("Concentration_Units", "Concentration_Unit", "Units", "Unit") %>%
  purrr::keep(~ .x %in% names(targeted)) %>% dplyr::first()

# Exposure in ng/mL -> log1p, then mean per (study, PID, analyte), long format.
build_exposure_long <- function(pid_col) {
  d <- targeted %>%
    filter(!is.na(.data[[pid_col]]), !is.na(Analyte_Code)) %>%
    transmute(
      Study,
      PID = .data[[pid_col]],
      analyte = Analyte_Code,
      conc = to_ng_per_mL(Concentration, .data[[unit_col]])
    ) %>%
    filter(!is.na(conc), conc >= 0) %>%
    mutate(x = log1p(conc))
  d %>%
    group_by(Study, PID, analyte) %>%
    summarise(x = mean(x, na.rm = TRUE), .groups = "drop")
}

# Continuous outcomes (n_distinct > 2 within panel), long format.
build_outcome_long <- function(pid_col) {
  oc <- outcomes %>% filter(!is.na(.data[[pid_col]]))
  meta <- c("Study", "Child_PID", "Mat_PID")
  cand <- setdiff(names(oc), meta)
  long <- oc %>%
    transmute(Study, PID = .data[[pid_col]],
              across(all_of(cand), as.character)) %>%
    pivot_longer(all_of(cand), names_to = "outcome", values_to = "y_raw") %>%
    mutate(y = suppressWarnings(as.numeric(y_raw))) %>%
    filter(!is.na(y)) %>%
    select(Study, PID, outcome, y)
  # Keep only genuinely continuous outcomes. Cardinality alone is not enough
  # (see OUTCOME_WHITELIST above): it admitted 64 questionnaire response codes.
  nu_tbl <- long %>% group_by(outcome) %>%
    summarise(nu = dplyr::n_distinct(y), .groups = "drop")

  cont <- nu_tbl %>% filter(nu >= MIN_DISTINCT_CONTINUOUS) %>% pull(outcome)
  if (!is.null(OUTCOME_WHITELIST)) {
    cont <- intersect(cont, OUTCOME_WHITELIST)
  }

  dropped <- setdiff(nu_tbl$outcome, cont)
  message(sprintf("[%s] outcomes: %d candidate columns -> %d analysed, %d dropped (not whitelisted or <%d distinct values)",
                  pid_col, nrow(nu_tbl), length(cont), length(dropped),
                  MIN_DISTINCT_CONTINUOUS))
  if (length(cont) == 0) {
    stop("No outcomes survived the whitelist. Check OUTCOME_WHITELIST spellings against names(outcomes).")
  }

  long %>% filter(outcome %in% cont)
}

## ---- Per-pair decomposition from per-study sufficient statistics ------------
# stats: per-study tibble with n, Sxx, Sxy, Syy, xbar, ybar (qualified studies)
decompose_pair <- function(st) {
  K <- nrow(st)
  N <- sum(st$n)

  # within-pooled (common slope, study fixed effects)
  sumSxx <- sum(st$Sxx); sumSxy <- sum(st$Sxy); sumSyy <- sum(st$Syy)
  b_within <- sumSxy / sumSxx
  rss_w <- sumSyy - b_within * sumSxy
  df_w  <- N - K - 1
  se_within <- if (df_w > 0 && sumSxx > 0) sqrt((rss_w / df_w) / sumSxx) else NA_real_

  # global means and naive pooled (within + between mixed)
  Xbar <- sum(st$n * st$xbar) / N
  Ybar <- sum(st$n * st$ybar) / N
  tot_Sxx <- sumSxx + sum(st$n * (st$xbar - Xbar)^2)
  tot_Sxy <- sumSxy + sum(st$n * (st$xbar - Xbar) * (st$ybar - Ybar))
  tot_Syy <- sumSyy + sum(st$n * (st$ybar - Ybar)^2)
  b_naive <- tot_Sxy / tot_Sxx
  rss_n <- tot_Syy - b_naive * tot_Sxy
  df_n  <- N - 2
  se_naive <- if (df_n > 0 && tot_Sxx > 0) sqrt((rss_n / df_n) / tot_Sxx) else NA_real_

  # between-study ecological slope (n-weighted regression of ybar on xbar)
  bxx <- sum(st$n * (st$xbar - Xbar)^2)
  bxy <- sum(st$n * (st$xbar - Xbar) * (st$ybar - Ybar))
  b_between <- if (bxx > 0) bxy / bxx else NA_real_

  # per-study slopes + SEs for heterogeneity (DerSimonian-Laird I^2)
  b_s  <- st$Sxy / st$Sxx
  rss_s <- st$Syy - b_s * st$Sxy
  s2_s  <- rss_s / (st$n - 2)
  se_s  <- sqrt(s2_s / st$Sxx)
  ok <- is.finite(b_s) & is.finite(se_s) & se_s > 0
  if (sum(ok) >= 2) {
    w <- 1 / se_s[ok]^2
    bbar <- sum(w * b_s[ok]) / sum(w)
    Q <- sum(w * (b_s[ok] - bbar)^2)
    dfQ <- sum(ok) - 1
    I2 <- max(0, (Q - dfQ) / Q) * 100
    n_pos <- sum(b_s[ok] > 0); n_neg <- sum(b_s[ok] < 0)
  } else {
    I2 <- NA_real_; n_pos <- NA_integer_; n_neg <- NA_integer_
  }

  # Standardized (partial-correlation-scale) effects so outcomes on different
  # units are comparable in one figure. within uses within-study SDs (cohort
  # partialled out); naive uses total SDs (= marginal correlation); between
  # uses the cohort-mean SDs.
  sd_xw <- sqrt(sumSxx / (N - K)); sd_yw <- sqrt(sumSyy / (N - K))
  sd_xt <- sqrt(tot_Sxx / (N - 1)); sd_yt <- sqrt(tot_Syy / (N - 1))
  byy <- sum(st$n * (st$ybar - Ybar)^2)
  r_within  <- if (sd_yw > 0) b_within  * sd_xw / sd_yw else NA_real_
  r_naive   <- if (sd_yt > 0) b_naive   * sd_xt / sd_yt else NA_real_
  r_between <- if (byy  > 0) b_between * sqrt(bxx / byy) else NA_real_

  tibble(
    K = K, N = N,
    b_within = b_within, se_within = se_within,
    b_naive = b_naive, se_naive = se_naive,
    b_between = b_between,
    r_within = r_within, r_naive = r_naive, r_between = r_between,
    I2 = I2, n_study_pos = n_pos, n_study_neg = n_neg
  )
}

## ---- Fast sufficient statistics + the permutation null ----------------------

# Per-group (n, Sxx, Sxy, Syy, xbar, ybar) via rowsum -- fast enough to run
# inside a permutation loop, unlike a dplyr group_by.
suff_stats <- function(g, x, y) {
  g <- as.character(g)
  n    <- as.vector(rowsum(rep(1, length(x)), g))
  sx   <- as.vector(rowsum(x, g));     sy  <- as.vector(rowsum(y, g))
  sxx  <- as.vector(rowsum(x * x, g)); syy <- as.vector(rowsum(y * y, g))
  sxy  <- as.vector(rowsum(x * y, g))
  xbar <- sx / n; ybar <- sy / n
  tibble(
    n = n,
    Sxx = sxx - n * xbar^2,
    Sxy = sxy - n * xbar * ybar,
    Syy = syy - n * ybar^2,
    xbar = xbar, ybar = ybar
  )
}

# Classify one pair from its sufficient statistics, returning just the label.
classify_from_stats <- function(st) {
  st <- st[st$n >= MIN_STUDY_N & st$Sxx > 0, , drop = FALSE]
  if (nrow(st) < MIN_STUDIES) return(NA_character_)
  d <- decompose_pair(st)
  if (!is.finite(d$b_within) || !is.finite(d$b_naive)) return(NA_character_)
  ratio <- d$b_naive / d$b_within
  if (sign(d$b_within) != sign(d$b_naive)) return("sign_reversal")
  if (ratio < RATIO_LO || ratio > RATIO_HI) return("magnitude_distort")
  "concordant"
}

# Null class rates for one pair, by permuting the cohort label.
permute_pair_null <- function(study, x, y, n_perm = N_PERM) {
  if (n_perm <= 0) {
    return(c(concordant = NA_real_, magnitude_distort = NA_real_,
             sign_reversal = NA_real_, n_valid = NA_real_))
  }
  labs <- replicate(n_perm, sample(study), simplify = FALSE)
  cls  <- vapply(labs, function(g) classify_from_stats(suff_stats(g, x, y)),
                 character(1))
  ok <- !is.na(cls)
  n_ok <- sum(ok)
  if (n_ok == 0) {
    return(c(concordant = NA_real_, magnitude_distort = NA_real_,
             sign_reversal = NA_real_, n_valid = 0))
  }
  c(concordant        = mean(cls[ok] == "concordant"),
    magnitude_distort = mean(cls[ok] == "magnitude_distort"),
    sign_reversal     = mean(cls[ok] == "sign_reversal"),
    n_valid           = n_ok)
}

## ---- Sweep the grid, one analyte at a time (memory-safe) ---------------------
run_panel <- function(pid_col, panel_name) {
  message(sprintf("[%s] building tables ...", panel_name))
  expo <- build_exposure_long(pid_col)
  outc <- build_outcome_long(pid_col)

  analytes <- sort(unique(expo$analyte))
  if (is.finite(MAX_ANALYTES)) analytes <- head(analytes, MAX_ANALYTES)
  message(sprintf("[%s] %d analytes x %d continuous outcomes",
                  panel_name, length(analytes), dplyr::n_distinct(outc$outcome)))

  res <- vector("list", length(analytes))
  cov <- vector("list", length(analytes))
  for (i in seq_along(analytes)) {
    a <- analytes[i]
    ex_a <- expo %>% filter(analyte == a) %>% select(Study, PID, x)
    joined <- inner_join(ex_a, outc, by = c("Study", "PID"))
    if (nrow(joined) == 0) next

    st <- joined %>%
      group_by(outcome, Study) %>%
      summarise(
        n = dplyr::n(),
        Sxx = sum((x - mean(x))^2),
        Sxy = sum((x - mean(x)) * (y - mean(y))),
        Syy = sum((y - mean(y))^2),
        xbar = mean(x), ybar = mean(y),
        .groups = "drop"
      ) %>%
      filter(n >= MIN_STUDY_N, Sxx > 0)

    if (nrow(st) == 0) next

    # coverage: how many cohorts is each pair observed in (the sparsity result)
    cov[[i]] <- st %>% count(outcome, name = "n_studies") %>%
      mutate(panel = panel_name, analyte = a)

    st_ok <- st %>% group_by(outcome) %>% filter(dplyr::n() >= MIN_STUDIES) %>%
      ungroup()
    if (nrow(st_ok) == 0) next

    dec <- st_ok %>% group_by(outcome) %>%
      group_modify(~ decompose_pair(.x)) %>% ungroup() %>%
      mutate(panel = panel_name, analyte = a, .before = 1)

    # Permutation null for each decomposed pair (2026-08-11). Uses the SAME rows
    # the observed statistics came from, restricted to the qualifying cohorts.
    if (N_PERM > 0 && nrow(dec) > 0) {
      qualifying <- st_ok %>% distinct(outcome, Study)
      null_rows <- purrr::map_dfr(dec$outcome, function(oc) {
        keep_studies <- qualifying$Study[qualifying$outcome == oc]
        sub <- joined[joined$outcome == oc & joined$Study %in% keep_studies, ]
        nl <- permute_pair_null(sub$Study, sub$x, sub$y)
        tibble(outcome = oc,
               null_pct_concordant    = 100 * nl[["concordant"]],
               null_pct_magnitude     = 100 * nl[["magnitude_distort"]],
               null_pct_sign_reversal = 100 * nl[["sign_reversal"]],
               n_perm_valid           = nl[["n_valid"]])
      })
      dec <- dec %>% left_join(null_rows, by = "outcome")
    }

    res[[i]] <- dec
    if (i %% 50 == 0) message(sprintf("[%s] %d/%d analytes", panel_name, i, length(analytes)))
  }
  list(pairs = bind_rows(res), coverage = bind_rows(cov))
}

child_run <- run_panel("Child_PID", "child")
mat_run   <- run_panel("Mat_PID",   "mat")
pairs_raw <- bind_rows(child_run$pairs, mat_run$pairs)
coverage  <- bind_rows(child_run$coverage, mat_run$coverage)

## ---- Classify ---------------------------------------------------------------
pairs <- pairs_raw %>%
  mutate(
    t_within = b_within / se_within,
    t_naive  = b_naive  / se_naive,
    ratio    = b_naive / b_within,
    sign_within  = sign(b_within),
    sign_naive   = sign(b_naive),
    sign_between = sign(b_between),
    class = case_when(
      sign_within != sign_naive ~ "sign_reversal",
      ratio < RATIO_LO | ratio > RATIO_HI ~ "magnitude_distort",
      TRUE ~ "concordant"
    ),
    # strong reversal: opposite signs AND both estimates significant
    strong_reversal = sign_within != sign_naive &
      is.finite(t_within) & is.finite(t_naive) &
      abs(t_within) > SIG_T & abs(t_naive) > SIG_T,
    between_opposes_within = is.finite(sign_between) & sign_between != sign_within
  )

## ---- Coverage / sparsity (the headline structural result) -------------------
coverage_summary <- coverage %>%
  group_by(panel) %>%
  summarise(
    pairs_total      = dplyr::n(),
    pairs_1_cohort   = sum(n_studies == 1),
    pairs_ge2_cohort = sum(n_studies >= 2),
    pairs_ge3_cohort = sum(n_studies >= 3),
    pairs_ge5_cohort = sum(n_studies >= 5),
    max_cohorts      = max(n_studies),
    pct_single_cohort = round(100 * mean(n_studies == 1), 1),
    .groups = "drop"
  )

## ---- Regime classification summary (decomposed pairs, >=2 and >=3 cohorts) ---
summarise_class <- function(df, label) {
  if (nrow(df) == 0) {
    return(tibble(subset = label, n_pairs = 0))
  }
  df %>% summarise(
    subset = label,
    n_pairs = dplyr::n(),
    pct_concordant        = round(100 * mean(class == "concordant"), 1),
    pct_magnitude_distort = round(100 * mean(class == "magnitude_distort"), 1),
    pct_sign_reversal     = round(100 * mean(class == "sign_reversal"), 1),
    # Null benchmark: the rate cohort-label permutation produces for these same
    # pairs. The OBSERVED rate is only evidence of real cohort structure to the
    # extent it exceeds this. Added 2026-08-11.
    null_pct_concordant    = if ("null_pct_concordant" %in% names(df))
      round(mean(null_pct_concordant, na.rm = TRUE), 1) else NA_real_,
    null_pct_magnitude     = if ("null_pct_magnitude" %in% names(df))
      round(mean(null_pct_magnitude, na.rm = TRUE), 1) else NA_real_,
    null_pct_sign_reversal = if ("null_pct_sign_reversal" %in% names(df))
      round(mean(null_pct_sign_reversal, na.rm = TRUE), 1) else NA_real_,
    excess_pct_sign_reversal = round(
      100 * mean(class == "sign_reversal") -
        (if ("null_pct_sign_reversal" %in% names(df))
          mean(null_pct_sign_reversal, na.rm = TRUE) else NA_real_), 1),
    pct_strong_reversal   = round(100 * mean(strong_reversal), 1),
    pct_between_opposes   = round(100 * mean(between_opposes_within, na.rm = TRUE), 1),
    median_I2 = round(median(I2, na.rm = TRUE), 1),
    pct_I2_over_50 = round(100 * mean(I2 > 50, na.rm = TRUE), 1)
  )
}

summary_tbl <- bind_rows(
  summarise_class(pairs, "all panels, >=2 cohorts"),
  summarise_class(filter(pairs, panel == "child"), "child, >=2 cohorts"),
  summarise_class(filter(pairs, panel == "mat"),   "mat, >=2 cohorts"),
  summarise_class(filter(pairs, K >= 3), "all panels, >=3 cohorts")
)

readr::write_csv(pairs,            file.path(out_dir, "cohort_regime_pairs.csv"))
readr::write_csv(summary_tbl,      file.path(out_dir, "cohort_regime_summary.csv"))
readr::write_csv(coverage,         file.path(out_dir, "cohort_regime_coverage.csv"))
readr::write_csv(coverage_summary, file.path(out_dir, "cohort_regime_coverage_summary.csv"))

cat("\nWrote diagnostic to:\n  ", normalizePath(out_dir), "\n\n", sep = "")
cat("=== COVERAGE (how many cohorts is each pair seen in, n>=20) ===\n")
print(as.data.frame(coverage_summary))
cat("\n=== REGIME CLASSIFICATION (decomposed pairs) ===\n")
print(as.data.frame(summary_tbl))

if (N_PERM > 0) {
  cat("\n=== READ THIS BEFORE QUOTING A REVERSAL RATE ===\n")
  cat("null_pct_* is the rate produced by permuting cohort labels within each\n")
  cat("pair (", N_PERM, " permutations/pair), holding N, K and the cohort size\n", sep = "")
  cat("profile fixed. Sign reversal is common under the null -- at K = 2 the two\n")
  cat("slopes are nearly independent in sign -- so quote excess_pct_sign_reversal\n")
  cat("(observed minus null), not pct_sign_reversal alone.\n")
}
