## ============================================================================
## budget_allocation_simulation.R
##
## Purpose
## -------
## Design argument behind the cohort-pooling simulation: under a FIXED budget,
## a deliberately DIVERSE set of locations beats unplanned accumulation of
## subjects from a few convenient sites. This script makes the tradeoff
## quantitative.
##
## Total sample size is held constant (N_total). We sweep the allocation from
## "few large cohorts" to "many small cohorts" -- i.e. K cohorts of n_per each
## with K * n_per = N_total -- and, at each split, measure three things you
## might want from the study:
##
##   (1) Average contextual effect  -- RMSE of the cohort-adjusted (within)
##       slope vs the true mean effect. A fixed N buys this at almost any split.
##   (2) Heterogeneity recovery     -- RMSE of the estimated cross-cohort effect
##       SD (DerSimonian-Laird tau) vs the truth, plus the probability of even
##       DETECTING heterogeneity. Needs MANY units; few large cohorts are blind
##       to it.
##   (3) Single-location precision  -- median standard error of a single
##       cohort's slope. Needs ENOUGH subjects per unit; many tiny cohorts make
##       every location's estimate noise.
##
## The point: (1) is flat, but (2) and (3) pull in opposite directions, so the
## defensible design is a diverse set of units each with enough subjects -- not
## "smaller cohorts," and not "pile everyone into a few convenient sites."
##
## This script ONLY writes results as CSVs; the figure is built in
##   Coding/Paper1_Submission2/Vis_Code/Pooling_Simulation_Visualization.Rmd
##
## Outputs (Output/Paper1_Output_Submission2/Pooling_Sim_Output/):
##   budget_allocation_results.csv  - one row per allocation (K, n_per) summary
##   budget_allocation_params.csv   - generator parameters (provenance)
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

set.seed(20260616)

resolve_existing_dir <- function(candidates) {
  hit <- candidates[dir.exists(candidates)]
  if (length(hit) == 0) {
    stop("Could not find Output/Paper1_Output_Submission2. Checked:\n",
         paste(candidates, collapse = "\n"))
  }
  hit[[1]]
}

submission_out <- resolve_existing_dir(c(
  "../../Output/Paper1_Output_Submission2",
  "../../../Output/Paper1_Output_Submission2",
  "Output/Paper1_Output_Submission2"
))
out_dir <- file.path(submission_out, "Pooling_Sim_Output")
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

## ---- Parameters -------------------------------------------------------------
## K_grid are divisors of N_total so every split uses the full budget exactly.
PARAMS <- list(
  N_total  = 6000,   # fixed budget: total subjects, constant across splits
  beta_bar = 0.30,   # super-population mean exposure effect
  tau_true = 0.30,   # super-population SD of cohort effects (heterogeneity)
  sigma_x  = 1,      # within-cohort exposure SD
  sigma_y  = 1,      # residual outcome SD
  mu_min   = 0,      # cohort exposure means spread (diverse profiles)
  mu_max   = 12,
  n_rep    = 150     # Monte-Carlo replicates per allocation
)
K_GRID <- c(2, 3, 5, 6, 10, 15, 20, 30, 50, 75, 100, 150, 300)

## ---- One allocation, one replicate ------------------------------------------
## Vectorized generation + grouped per-cohort regression summaries (slope b_k
## and its SE) computed in closed form -- no per-cohort lm() calls.
one_alloc <- function(K, n_per, p) {
  cohort  <- rep(seq_len(K), each = n_per)
  mu_k    <- runif(K, p$mu_min, p$mu_max)
  beta_k  <- rnorm(K, p$beta_bar, p$tau_true)
  alpha_k <- rnorm(K, 0, 2)

  x <- rnorm(K * n_per, mean = mu_k[cohort], sd = p$sigma_x)
  y <- alpha_k[cohort] + beta_k[cohort] * x + rnorm(K * n_per, sd = p$sigma_y)

  stats <- tibble(cohort, x, y) %>%
    group_by(cohort) %>%
    summarise(
      n   = dplyr::n(),
      Sxx = sum((x - mean(x))^2),
      Sxy = sum((x - mean(x)) * (y - mean(y))),
      Syy = sum((y - mean(y))^2),
      .groups = "drop"
    ) %>%
    mutate(
      b      = Sxy / Sxx,
      sigma2 = (Syy - b * Sxy) / (n - 2),
      se     = sqrt(sigma2 / Sxx)
    )

  # (1) average contextual effect = common within-cohort slope estimator
  beta_within <- sum(stats$Sxy) / sum(stats$Sxx)

  # (2) heterogeneity: DerSimonian-Laird tau + Cochran's Q detection test
  w     <- 1 / stats$se^2
  bbar  <- sum(w * stats$b) / sum(w)
  Q     <- sum(w * (stats$b - bbar)^2)
  Cterm <- sum(w) - sum(w^2) / sum(w)
  tau2  <- max(0, (Q - (K - 1)) / Cterm)
  tau_hat   <- sqrt(tau2)
  detect_het <- as.numeric(Q > stats::qchisq(0.95, df = K - 1))

  # (3) single-location precision
  med_se <- median(stats$se)

  tibble(beta_within = beta_within, tau_hat = tau_hat,
         detect_het = detect_het, med_se = med_se)
}

## ---- Sweep the allocation grid ----------------------------------------------
results <- purrr::map_dfr(K_GRID, function(K) {
  n_per <- PARAMS$N_total / K
  reps  <- purrr::map_dfr(seq_len(PARAMS$n_rep),
                          function(i) one_alloc(K, n_per, PARAMS))
  tibble(
    K = K,
    n_per = n_per,
    # (1) average contextual effect: RMSE vs truth (want LOW; ~flat in K)
    beta_within_rmse = sqrt(mean((reps$beta_within - PARAMS$beta_bar)^2)),
    # (2) heterogeneity recovery: RMSE of tau and detection probability
    tau_rmse   = sqrt(mean((reps$tau_hat - PARAMS$tau_true)^2)),
    tau_hat_mean = mean(reps$tau_hat),
    detect_rate = mean(reps$detect_het),
    # (3) single-location precision: median SE of one cohort's slope (want LOW)
    median_location_se = mean(reps$med_se)
  )
})

params_tbl <- tibble(parameter = names(PARAMS),
                     value = unlist(PARAMS, use.names = FALSE))

readr::write_csv(results,    file.path(out_dir, "budget_allocation_results.csv"))
readr::write_csv(params_tbl, file.path(out_dir, "budget_allocation_params.csv"))

cat("Wrote budget-allocation results to:\n  ", normalizePath(out_dir), "\n", sep = "")
print(results, n = nrow(results), digits = 3)
