## ============================================================================
## cohort_pooling_simulation.R
##
## Purpose
## -------
## Reviewer critique: "Your cohorts differ so much by construction that the
## sources of variation are contrived." Our rebuttal: an exposure-outcome
## association only has meaning relative to the UNIT (cohort/location) it
## emerged from, and pooling across units can DESTROY or even REVERSE the
## very signal that is real within each unit.
##
## A colleague suggested an analogy to their longitudinal 'n=1 -> population'
## work. That framework is the within-person vs. between-person decomposition.
## The cross-sectional analog here is identical in structure, with COHORT
## playing the role that SUBJECT plays in longitudinal data:
##
##     longitudinal: within-person  vs. between-person  vs. naive-pooled
##     this study  : within-cohort  vs. between-cohort   vs. naive-pooled
##
## So we are NOT importing their repeated-measures method; we are borrowing the
## one idea that does transfer: associations are non-collapsible across the
## clustering unit (Simpson's paradox / ecological-aggregation bias).
##
## This script ONLY generates and writes the simulation results as CSVs.
## All tables and figures are built downstream in
##   Coding/Paper1_Submission2/Vis_Code/Pooling_Simulation_Visualization.Rmd
##
## Outputs (written to Output/Paper1_Output_Submission2/Pooling_Sim_Output/):
##   pooling_sim_individual.csv   - individual-level simulated data (both regimes)
##   pooling_sim_per_cohort.csv   - per-cohort estimated & true slopes
##   pooling_sim_estimands.csv    - within / naive-pooled / between summaries
##   pooling_sim_curve.csv        - "toward the limit" pooling curve
##   pooling_sim_params.csv       - generator parameters used (provenance)
## ============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
})

set.seed(20260616)

## ---- Resolve output directory (mirrors the other Submission2 scripts) -------
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

## ---- Generator parameters (single source of truth) --------------------------
PARAMS <- list(
  K           = 50,    # number of cohorts / locations
  n_per       = 250,   # subjects per cohort
  sigma_x     = 1,     # within-cohort exposure SD
  sigma_y     = 1,     # residual outcome SD
  beta_bar    = 0.30,  # averaging regime: population-mean exposure effect
  tau         = 0.45,  # averaging regime: SD of cohort effects (heterogeneity)
  beta_within = 0.50,  # simpson regime: common within-cohort effect (>0)
  gamma       = 1.20,  # simpson regime: strength of location confounding
  a0          = 0,     # simpson regime: baseline intercept anchor
  mu_min      = 0,     # cohort exposure means spread (different profiles)
  mu_max      = 12,
  n_draw      = 200    # random subsets per pooling-curve point
)

## ----------------------------------------------------------------------------
## Data generator (see header for the model). regime in {averaging, simpson}.
## ----------------------------------------------------------------------------
simulate_cohorts <- function(regime, p = PARAMS) {
  mu_k <- runif(p$K, min = p$mu_min, max = p$mu_max)  # per-location profiles

  if (regime == "averaging") {
    beta_k  <- rnorm(p$K, mean = p$beta_bar, sd = p$tau)  # heterogeneous, some<0
    alpha_k <- rnorm(p$K, mean = 0, sd = 2)               # baseline ⟂ exposure
    # Orthogonalize cohort effects against the cohort exposure level so there is
    # ZERO cross-cohort confounding by construction. The naive pooled slope picks
    # up not just cov(alpha_k, mu_k) but also the beta_k x mu_k^2 interaction, so
    # we residualize on BOTH mu_k and mu_k^2 (then restore the intended mean).
    # This makes the between-cohort trend equal the within-cohort effect for any
    # K/seed, so the regime demonstrates honest blurring rather than incidental
    # bias: pooled ≈ within, only the cohort heterogeneity is averaged away.
    beta_k  <- as.numeric(residuals(lm(beta_k ~ mu_k + I(mu_k^2))))  + p$beta_bar
    alpha_k <- as.numeric(residuals(lm(alpha_k ~ mu_k + I(mu_k^2))))
  } else if (regime == "simpson") {
    beta_k  <- rep(p$beta_within, p$K)                    # positive everywhere
    alpha_k <- p$a0 - p$gamma * mu_k                      # confounding by site
  } else {
    stop("Unknown regime: ", regime)
  }

  purrr::map_dfr(seq_len(p$K), function(k) {
    x <- rnorm(p$n_per, mean = mu_k[k], sd = p$sigma_x)
    y <- alpha_k[k] + beta_k[k] * x + rnorm(p$n_per, sd = p$sigma_y)
    tibble(regime = regime, cohort = k, x = x, y = y,
           mu_k = mu_k[k], beta_k = beta_k[k], alpha_k = alpha_k[k])
  })
}

## ---- Three estimands of "the association" -----------------------------------
##   naive   : lm(y ~ x)            -> mixes within + between
##   within  : lm(y ~ x + cohort)   -> contextual (cohort-adjusted)
##   between : lm(ybar ~ xbar)       -> ecological / location-level
estimands <- function(df) {
  d <- df %>% mutate(cohort = factor(cohort))
  agg <- d %>% group_by(cohort) %>%
    summarise(xbar = mean(x), ybar = mean(y), .groups = "drop")
  tibble(
    regime          = df$regime[[1]],
    within_true     = mean(unique(df$beta_k)),
    within_adjusted = coef(lm(y ~ x + cohort, data = d))[["x"]],
    naive_pooled    = coef(lm(y ~ x, data = d))[["x"]],
    between_cohort  = coef(lm(ybar ~ xbar, data = agg))[["xbar"]]
  )
}

per_cohort_slopes <- function(df) {
  df %>% group_by(regime, cohort) %>%
    summarise(beta_hat  = coef(lm(y ~ x))[["x"]],
              beta_true = first(beta_k),
              mu_k      = first(mu_k),
              n         = dplyr::n(),
              .groups = "drop")
}

## ---- "Sampling unit approaches the limit" curve -----------------------------
## For m = 1..K pooled cohorts, draw n_draw random subsets, fit naive lm on each.
pooling_curve <- function(df, p = PARAMS) {
  cohorts <- unique(df$cohort)
  purrr::map_dfr(seq_len(p$K), function(m) {
    slopes <- purrr::map_dbl(seq_len(p$n_draw), function(i) {
      pick <- sample(cohorts, m)
      sub  <- df[df$cohort %in% pick, ]
      tryCatch(coef(lm(y ~ x, data = sub))[["x"]], error = function(e) NA_real_)
    })
    tibble(regime = df$regime[[1]], m = m,
           mean_slope = mean(slopes, na.rm = TRUE),
           lo = stats::quantile(slopes, 0.05, na.rm = TRUE),
           hi = stats::quantile(slopes, 0.95, na.rm = TRUE))
  })
}

## ---- Run both regimes and assemble result tables ----------------------------
regimes <- c("averaging", "simpson")
sim_list <- purrr::map(regimes, simulate_cohorts)

individual <- bind_rows(sim_list)
per_cohort <- purrr::map_dfr(sim_list, per_cohort_slopes)
estimand_tbl <- purrr::map_dfr(sim_list, estimands)
curve_tbl <- purrr::map_dfr(sim_list, pooling_curve)
params_tbl <- tibble(parameter = names(PARAMS),
                     value = unlist(PARAMS, use.names = FALSE))

## ---- Write -----------------------------------------------------------------
readr::write_csv(individual,   file.path(out_dir, "pooling_sim_individual.csv"))
readr::write_csv(per_cohort,   file.path(out_dir, "pooling_sim_per_cohort.csv"))
readr::write_csv(estimand_tbl, file.path(out_dir, "pooling_sim_estimands.csv"))
readr::write_csv(curve_tbl,    file.path(out_dir, "pooling_sim_curve.csv"))
readr::write_csv(params_tbl,   file.path(out_dir, "pooling_sim_params.csv"))

cat("Wrote simulation results to:\n  ", normalizePath(out_dir), "\n", sep = "")
print(estimand_tbl, digits = 3)
