median_or_na <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  stats::median(x)
}

outcome_support_cols <- function(outcome_col) {
  c(
    outcome_col,
    paste0(outcome_col, c("_LOD", "_Detected", "_ExposureObs"))
  )
}

select_outcome_data <- function(df, join_var, outcome_col) {
  df %>%
    dplyr::select(all_of(c(join_var, "Study", outcome_support_cols(outcome_col))))
}

build_censored_targeted_panel <- function(targeted_df, id_var, unit_converter) {
  targeted_df %>%
    dplyr::filter(!is.na(.data[[id_var]]), !is.na(Analyte_Code)) %>%
    dplyr::mutate(
      Units = dplyr::if_else(is.na(Units), "ng/mL", Units),
      Exposure = unit_converter(Concentration, Units),
      LOD = unit_converter(LOD, Units),
      Detected = dplyr::case_when(
        !is.na(Detected) ~ as.logical(Detected),
        !is.na(Exposure) & !is.na(LOD) ~ Exposure >= LOD,
        !is.na(Exposure) & is.na(LOD) ~ TRUE,
        TRUE ~ NA
      )
    ) %>%
    dplyr::group_by(Study, .data[[id_var]], Analyte_Code) %>%
    dplyr::summarise(
      Exposure = mean(Exposure, na.rm = TRUE),
      LOD = median_or_na(LOD),
      # CENSORING PRESERVED (fixed 2026-08-11). The lab-reported Detected flag was
      # computed in the mutate above and then SILENTLY DROPPED here -- summarise()
      # did not carry it -- so the mutate below re-derived it as Exposure >= LOD.
      #
      # That comparison destroys the censoring, because the harmonized data stores
      # non-detects AT the LOD, not below it:
      #     PFOA: conc < LOD = 0 rows, conc == LOD = 192 rows (the 191 lab non-detects)
      # Across all analytes, `Exposure >= LOD` reclassified 154,413 of 154,939
      # lab-flagged non-detects (99.7%) as DETECTED. Every analyte then reached
      # crq() with 0% censoring, which makes the Portnoy fit degenerate and
      # summary.crq() throw -- silently caught, yielding NA inference throughout.
      #
      # Majority rule across replicates; for the usual one-measurement-per-person
      # case it is just the lab flag itself.
      Detected_lab = if (all(is.na(Detected))) NA else mean(as.logical(Detected), na.rm = TRUE) >= 0.5,
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      Exposure = dplyr::if_else(is.nan(Exposure) | is.infinite(Exposure), NA_real_, Exposure),
      LOD = dplyr::if_else(is.nan(LOD) | is.infinite(LOD), NA_real_, LOD),
      Detected = dplyr::case_when(
        # authoritative: what the lab actually reported
        !is.na(Detected_lab) ~ Detected_lab,
        # fallback for the ~6% of rows with no flag. STRICT `>`: a value sitting
        # exactly at the LOD is a non-detect by the substitution convention above.
        !is.na(Exposure) & !is.na(LOD) ~ Exposure > LOD,
        !is.na(Exposure) & is.na(LOD) ~ TRUE,
        TRUE ~ NA
      ),
      Exposure_Obs = dplyr::case_when(
        !is.na(Exposure) & (Detected %in% TRUE | is.na(LOD)) ~ Exposure,
        !is.na(LOD) ~ LOD,
        TRUE ~ NA_real_
      )
    ) %>%
    dplyr::select(-Detected_lab)
}

build_censored_targeted_list <- function(panel_df, id_var, analytes = sort(unique(panel_df$Analyte_Code))) {
  id_sym <- rlang::sym(id_var)

  purrr::map(
    analytes,
    function(analyte) {
      panel_df %>%
        dplyr::filter(Analyte_Code == analyte) %>%
        dplyr::mutate(Study = std_study(Study)) %>%
        dplyr::transmute(
          Study,
          !!id_sym := .data[[id_var]],
          !!analyte := Exposure,
          !!paste0(analyte, "_LOD") := LOD,
          !!paste0(analyte, "_Detected") := Detected,
          !!paste0(analyte, "_ExposureObs") := Exposure_Obs
        ) %>%
        dplyr::distinct()
    }
  ) %>%
    rlang::set_names(analytes)
}

rho_tau <- function(u, tau) {
  u * (tau - (u < 0))
}

censored_quantile_loss <- function(exposure_obs, detected, fitted_quantile, tau) {
  detected_idx <- detected %in% TRUE
  censored_idx <- detected %in% FALSE

  detected_loss <- sum(
    rho_tau(exposure_obs[detected_idx] - fitted_quantile[detected_idx], tau),
    na.rm = TRUE
  )

  censored_loss <- sum(
    dplyr::if_else(
      fitted_quantile[censored_idx] > exposure_obs[censored_idx],
      rho_tau(exposure_obs[censored_idx] - fitted_quantile[censored_idx], tau),
      0
    ),
    na.rm = TRUE
  )

  detected_loss + censored_loss
}

build_df_cens <- function(df, outcome_col) {
  lod_col <- paste0(outcome_col, "_LOD")
  detected_col <- paste0(outcome_col, "_Detected")
  observed_col <- paste0(outcome_col, "_ExposureObs")

  df %>%
    dplyr::mutate(
      exposure = .data[[outcome_col]],
      LOD = .data[[lod_col]],
      detected = dplyr::case_when(
        !is.na(.data[[detected_col]]) ~ as.logical(.data[[detected_col]]),
        !is.na(exposure) & !is.na(LOD) ~ exposure >= LOD,
        !is.na(exposure) & is.na(LOD) ~ TRUE,
        TRUE ~ NA
      ),
      exposure_obs = dplyr::coalesce(
        .data[[observed_col]],
        dplyr::case_when(
          !is.na(exposure) & (detected %in% TRUE | is.na(LOD)) ~ exposure,
          !is.na(LOD) ~ LOD,
          TRUE ~ NA_real_
        )
      )
    ) %>%
    dplyr::filter(!is.na(exposure_obs), !is.na(detected))
}

choose_censored_taus <- function(df_cens) {
  if (mean(!df_cens$detected, na.rm = TRUE) > 0.50) {
    c(0.50, 0.75, 0.95)
  } else {
    c(0.25, 0.50, 0.75, 0.95)
  }
}

keep_varying_terms <- function(data, rhs_formula) {
  rhs_terms <- attr(stats::terms(rhs_formula), "term.labels")
  if (length(rhs_terms) == 0) {
    return(rhs_formula)
  }

  varying_terms <- rhs_terms[
    vapply(
      rhs_terms,
      function(term) dplyr::n_distinct(data[[term]][!is.na(data[[term]])]) > 1,
      logical(1)
    )
  ]

  if (length(varying_terms) == 0) {
    stats::as.formula("~ 1")
  } else {
    stats::reformulate(varying_terms)
  }
}

make_crq_predictions <- function(fit, rhs_formula, data, taus_use) {
  predicted_vals <- tryCatch(
    predict(fit, newdata = data),
    error = function(e) NULL
  )

  if (!is.null(predicted_vals)) {
    pred_mat <- as.matrix(predicted_vals)
    if (nrow(pred_mat) == nrow(data)) {
      if (ncol(pred_mat) == 1 && length(taus_use) > 1) {
        pred_mat <- matrix(pred_mat, nrow = nrow(data), ncol = length(taus_use))
      }
      if (ncol(pred_mat) == length(taus_use)) {
        colnames(pred_mat) <- as.character(taus_use)
        return(pred_mat)
      }
    }
  }

  design_mat <- stats::model.matrix(rhs_formula, data = data)
  coef_raw <- stats::coef(fit, taus = taus_use)
  if (is.null(dim(coef_raw))) {
    coef_mat <- matrix(
      coef_raw,
      nrow = ncol(design_mat),
      dimnames = list(colnames(design_mat), as.character(taus_use))
    )
  } else {
    coef_mat <- as.matrix(coef_raw)
  }

  if (is.null(rownames(coef_mat)) && nrow(coef_mat) == ncol(design_mat)) {
    rownames(coef_mat) <- colnames(design_mat)
  }

  coef_terms <- rownames(coef_mat)
  if (is.null(coef_terms)) {
    stop("Coefficient matrix does not align with design matrix")
  }

  coef_lookup <- stats::setNames(seq_along(coef_terms), coef_terms)
  aligned_coef <- matrix(
    0,
    nrow = ncol(design_mat),
    ncol = ncol(coef_mat),
    dimnames = list(colnames(design_mat), colnames(coef_mat))
  )

  matched_terms <- intersect(colnames(design_mat), names(coef_lookup))
  if (length(matched_terms) == 0) {
    stop("Coefficient matrix does not align with design matrix")
  }

  aligned_coef[matched_terms, ] <- coef_mat[coef_lookup[matched_terms], , drop = FALSE]
  predictions <- design_mat %*% aligned_coef
  colnames(predictions) <- as.character(taus_use)
  predictions
}

# IDENTIFICATION GUARD ON pseudo_R (added 2026-08-11).
#
# censored_quantile_loss() charges a censored observation NOTHING when the fitted
# quantile sits at or below its detection limit -- correct, since any value under
# the bound is consistent with the data, but it means the loss is NOT a proper
# scoring rule: it is minimised by predicting arbitrarily low. Whenever tau falls
# at or below the censored fraction, the fitted quantile lies inside the censored
# region for essentially every observation, full_loss collapses to 0, and
#
#     pseudo_R = 1 - 0 / null_loss = 1
#
# mechanically -- not because the model explains anything. Demonstrated directly:
# at 83% censoring, full_loss = 0.000 and pseudo_R = 1.0000 at tau = 0.25/0.50/0.75.
#
# This was invisible until now only because the censoring bug (see
# build_censored_targeted_panel) meant nothing was ever marked censored. With
# censoring restored, 14.7% of model rows came back at exactly 1.0, rising to
# 63.6% of rows for analytes above 80% censoring.
#
# A quantile at tau <= the censored fraction is not identified from the data, so
# report NA there rather than a number that is guaranteed to be 1. Above that
# threshold the full-vs-null comparison is meaningful: both models face the same
# censored observations.
compute_censored_pseudo_r <- function(df_cens, fitted_full, fitted_null, taus_use) {
  purrr::map2_dfr(
    seq_along(taus_use),
    taus_use,
    function(idx, tau) {
      ff <- fitted_full[, idx]
      fn <- fitted_null[, idx]

      # ESTIMABILITY (2026-08-11). crq() returns NA fitted values at quantiles it
      # cannot identify -- typically at or below the censored fraction, where the
      # tau-th quantile lies inside the region we only know as "below the LOD".
      # That is the estimator being honest.
      #
      # censored_quantile_loss() sums with na.rm = TRUE, so an all-NA fit used to
      # score full_loss = 0 -- i.e. a PERFECT fit. When the simpler null model
      # still converged, pseudo_R came out as exactly 1.0: 99.6% of the 1.0 rows
      # had no coefficients at all. Subset to the rows both models actually fitted
      # and mark the tau unestimable when there are none, so a failed fit reads as
      # NA rather than as a perfect score.
      ok <- !is.na(ff) & !is.na(fn)
      if (!any(ok)) {
        return(tibble::tibble(
          Tau = tau, pseudo_R = NA_real_, Full_Loss = NA_real_, Null_Loss = NA_real_,
          Tau_Estimable = FALSE, N_Used = 0L
        ))
      }

      full_loss <- censored_quantile_loss(
        exposure_obs = df_cens$exposure_obs[ok],
        detected = df_cens$detected[ok],
        fitted_quantile = ff[ok],
        tau = tau
      )
      null_loss <- censored_quantile_loss(
        exposure_obs = df_cens$exposure_obs[ok],
        detected = df_cens$detected[ok],
        fitted_quantile = fn[ok],
        tau = tau
      )

      estimable <- is.finite(full_loss) && is.finite(null_loss) && null_loss > 0

      tibble::tibble(
        Tau = tau,
        pseudo_R = if (estimable) 1 - full_loss / null_loss else NA_real_,
        Full_Loss = if (estimable) full_loss else NA_real_,
        Null_Loss = if (estimable) null_loss else NA_real_,
        Tau_Estimable = estimable,
        N_Used = sum(ok)
      )
    }
  )
}

#' One goodness-of-fit number for a model, across the quantiles it could fit.
#'
#' Weighted by each quantile's null loss, which makes the result identical to the
#' pooled R1 computed from summed losses:
#'
#'   sum_t( L_null_t * R_t ) / sum_t( L_null_t )
#'     = sum_t( L_null_t - L_full_t ) / sum_t( L_null_t )
#'     = 1 - sum_t(L_full_t) / sum_t(L_null_t)
#'
#' so quantiles with more variation to explain carry proportionally more weight,
#' and the aggregate is on the same 0-1 scale as the per-tau values. Quantiles
#' that could not be estimated are excluded entirely rather than treated as zero.
#'
#' Also returns the unweighted mean for comparison: if the two diverge sharply,
#' the fit is being carried by one quantile and the aggregate should be read with
#' the per-tau values beside it.
#' @param require_taus quantiles that MUST all be estimable, and the only ones
#'   aggregated over. Defaults to the top three, c(0.5, 0.75, 0.95).
#'
#'   COMPARABILITY (2026-08-11). Without this, the aggregate averages over
#'   whatever each model happened to fit, so a heavily censored analyte resting on
#'   tau = 0.95 alone was being compared against one averaging four quantiles --
#'   different functionals wearing the same column name. Measured on TEDDY: 93% of
#'   models had all four taus, 12 had {0.5,0.75,0.95}, 12 had {0.95} only.
#'
#'   Fixing the grid at the top three keeps 96% of models (336/350), drops the
#'   quantile most often lost to censoring, and makes every reported value the
#'   same quantity. Models missing any required tau return NA rather than a
#'   number computed on a different basis -- which also excludes the k = 1 fits
#'   that produced nonsense (down to -8.4).
#'
#'   Pass require_taus = NULL for the old "whatever was estimable" behaviour.
#'
#'   NOTE this does NOT make the value comparable across model SPECIFICATIONS.
#'   Coefficient count drives it hard (Spearman 0.748 on TEDDY: median R2 0.046
#'   at 9 coefficients vs 0.008 at 1-3). Compare analytes within a fixed
#'   covariate set; use a nested-model test, not R2, to compare specifications.
aggregate_pseudo_r <- function(r1_df, require_taus = c(0.5, 0.75, 0.95)) {
  none <- list(overall = NA_real_, equal_weight = NA_real_,
               n_taus = 0L, taus = NA_character_, basis = NA_character_)

  ok <- r1_df$Tau_Estimable %in% TRUE

  if (!is.null(require_taus)) {
    # every required tau must have been fitted AND estimable
    have <- r1_df$Tau[ok]
    if (!all(require_taus %in% have)) return(none)
    ok <- ok & (r1_df$Tau %in% require_taus)
  }

  if (!any(ok)) return(none)

  list(
    overall      = 1 - sum(r1_df$Full_Loss[ok]) / sum(r1_df$Null_Loss[ok]),
    equal_weight = mean(r1_df$pseudo_R[ok]),
    n_taus       = sum(ok),
    taus         = paste(r1_df$Tau[ok], collapse = ","),
    basis        = if (is.null(require_taus)) "estimable" else "fixed_grid"
  )
}

extract_crq_table <- function(summary_obj, outcome_col, taus_use, exposure_col = NULL) {
  purrr::map_dfr(
    seq_along(summary_obj),
    function(idx) {
      curr_summary <- summary_obj[[idx]]
      coef_df <- as.data.frame(curr_summary$coefficients) %>%
        tibble::rownames_to_column("Variable") %>%
        dplyr::mutate(
          Tau = curr_summary$tau,
          Outcome = outcome_col
        )

      coef_df <- coef_df %>%
        dplyr::rename(
          Beta = Value,
          Lower = `Lower Bd`,
          Upper = `Upper Bd`,
          Standard_Error = `Std Error`,
          T_Value = `T Value`,
          P_Value = `Pr(>|t|)`
        )

      if (is.null(exposure_col)) {
        coef_df %>%
          dplyr::select(
            Variable,
            Beta,
            Standard_Error,
            T_Value,
            P_Value,
            Tau,
            Outcome,
            Lower,
            Upper
          )
      } else {
        coef_df %>%
          dplyr::mutate(Exposure = exposure_col) %>%
          dplyr::select(
            Variable,
            Beta,
            Standard_Error,
            T_Value,
            P_Value,
            Tau,
            Outcome,
            Exposure,
            Lower,
            Upper
          )
      }
    }
  )
}

fallback_crq_summary <- function(fit, rhs_formula, data, taus_use) {
  design_mat <- stats::model.matrix(rhs_formula, data = data)
  coef_raw <- stats::coef(fit, taus = taus_use)
  if (is.null(dim(coef_raw))) {
    coef_mat <- matrix(
      coef_raw,
      nrow = ncol(design_mat),
      dimnames = list(colnames(design_mat), as.character(taus_use))
    )
  } else {
    coef_mat <- as.matrix(coef_raw)
  }
  if (is.null(rownames(coef_mat)) && nrow(coef_mat) == ncol(design_mat)) {
    rownames(coef_mat) <- colnames(design_mat)
  }

  purrr::map(
    seq_along(taus_use),
    function(idx) {
      coeff_mat <- cbind(
        Value = coef_mat[, idx],
        `Lower Bd` = NA_real_,
        `Upper Bd` = NA_real_,
        `Std Error` = NA_real_,
        `T Value` = NA_real_,
        `Pr(>|t|)` = NA_real_
      )
      list(
        tau = taus_use[[idx]],
        coefficients = coeff_mat
      )
    }
  )
}

# Plain quantile regression for analytes with no censoring. Returns the same
# shape as fit_censored_quantile_models() so everything downstream is unchanged.
# summary.rq(se = "boot") reports Value / Std. Error / t value / Pr(>|t|); the
# Lower Bd / Upper Bd columns extract_crq_table() expects are derived from the
# normal approximation, matching how the crq path reports its bootstrap interval.
fit_uncensored_quantile_models <- function(df_cens, rhs_formula, taus_use, summary_R) {
  model_formula <- stats::update(rhs_formula, exposure_obs ~ .)
  full_fit <- suppressWarnings(tryCatch(
    do.call(quantreg::rq, list(formula = model_formula, data = df_cens, tau = taus_use)),
    error = function(e) NULL))
  null_fit <- suppressWarnings(tryCatch(
    do.call(quantreg::rq, list(formula = exposure_obs ~ 1, data = df_cens, tau = taus_use)),
    error = function(e) NULL))
  if (is.null(full_fit) || is.null(null_fit)) return(NULL)

  full_summary <- NULL
  if (summary_R >= 25) {
    full_summary <- suppressWarnings(tryCatch({
      s <- summary(full_fit, se = "boot", R = summary_R)
      if (length(taus_use) == 1) s <- list(s)
      lapply(s, function(si) {
        cf <- si$coefficients
        est <- cf[, "Value"]; se <- cf[, "Std. Error"]
        list(tau = si$tau, coefficients = cbind(
          Value = est,
          `Lower Bd` = est - stats::qnorm(0.975) * se,
          `Upper Bd` = est + stats::qnorm(0.975) * se,
          `Std Error` = se,
          `T Value` = cf[, 3],
          `Pr(>|t|)` = cf[, 4]))
      })
    }, error = function(e) NULL))
  }
  if (is.null(full_summary)) {
    full_summary <- fallback_crq_summary(full_fit, rhs_formula, df_cens, taus_use)
  }

  fitted_full <- tryCatch(as.matrix(stats::predict(full_fit, newdata = df_cens)), error = function(e) NULL)
  fitted_null <- tryCatch(as.matrix(stats::predict(null_fit, newdata = df_cens)), error = function(e) NULL)
  if (is.null(fitted_full) || is.null(fitted_null)) return(NULL)
  if (ncol(fitted_full) != length(taus_use)) fitted_full <- matrix(fitted_full, nrow = nrow(df_cens), ncol = length(taus_use))
  if (ncol(fitted_null) != length(taus_use)) fitted_null <- matrix(fitted_null, nrow = nrow(df_cens), ncol = length(taus_use))

  list(full_fit = full_fit, null_fit = null_fit, full_summary = full_summary,
       fitted_full = fitted_full, fitted_null = fitted_null)
}

fit_censored_quantile_models <- function(df_cens, rhs_formula, taus_use, summary_R = 2) {
  # Minimum bootstrap replicates required before we report standard errors / CIs /
  # p-values. A bootstrap with only a couple of replicates is not meaningful, so
  # below this threshold we return point estimates with NA inference (via
  # fallback_crq_summary) rather than unstable numbers. This also preserves the
  # historical behavior at the default summary_R = 2 (NA inference); raise
  # summary_R (e.g. 200+) to obtain valid bootstrap CIs / p-values.
  min_boot_R <- 25

  # UNCENSORED ANALYTES USE rq(), NOT crq() (added 2026-08-11).
  #
  # crq(method = "Portnoy") is degenerate when nothing is censored: summary.crq()
  # then throws "argument is of length zero", which the tryCatch below swallows,
  # and the whole model silently returns NA inference. Some analytes legitimately
  # have no non-detects, so route those to plain quantile regression, which is the
  # correct estimator for that case and supports a bootstrap summary.
  n_censored <- sum(!df_cens$detected, na.rm = TRUE)
  if (n_censored == 0) {
    return(fit_uncensored_quantile_models(df_cens, rhs_formula, taus_use, summary_R))
  }

  surv_formula <- stats::update(rhs_formula, survival::Surv(exposure_obs, detected, type = "left") ~ .)
  null_formula <- survival::Surv(exposure_obs, detected, type = "left") ~ 1

  # Fits are built with do.call(), splicing the actual formula and data objects
  # into the stored model call. crq() records its call by symbol name, and its
  # bootstrap SE routine (summary(..., se = "boot")) re-evaluates that call; when
  # the formula/data are referenced by local variable names, the bootstrap cannot
  # resolve them from its own evaluation frame and silently returns NA standard
  # errors regardless of summary_R. Embedding the objects makes the call
  # self-contained so bootstrap inference works when requested.
  full_fit <- suppressWarnings(
    tryCatch(
      do.call(quantreg::crq,
              list(formula = surv_formula, data = df_cens, taus = taus_use, method = "Portnoy")),
      error = function(e) NULL
    )
  )
  null_fit <- suppressWarnings(
    tryCatch(
      do.call(quantreg::crq,
              list(formula = null_formula, data = df_cens, taus = taus_use, method = "Portnoy")),
      error = function(e) NULL
    )
  )

  if (is.null(full_fit) || is.null(null_fit)) {
    return(NULL)
  }

  full_summary <- if (summary_R >= min_boot_R) {
    suppressWarnings(
      tryCatch(
        summary(full_fit, taus = taus_use, se = "boot", R = summary_R),
        error = function(e) NULL
      )
    )
  } else {
    NULL
  }

  if (is.null(full_summary)) {
    full_summary <- fallback_crq_summary(full_fit, rhs_formula, df_cens, taus_use)
  }

  fitted_full <- tryCatch(make_crq_predictions(full_fit, rhs_formula, df_cens, taus_use), error = function(e) NULL)
  fitted_null <- tryCatch(make_crq_predictions(null_fit, ~ 1, df_cens, taus_use), error = function(e) NULL)

  if (is.null(fitted_full) || is.null(fitted_null)) {
    return(NULL)
  }

  list(
    full_fit = full_fit,
    null_fit = null_fit,
    full_summary = full_summary,
    fitted_full = fitted_full,
    fitted_null = fitted_null
  )
}

run_censored_quantile_reg <- function(exp_var = NULL,
                                      dataset,
                                      my_outcome,
                                      my_covariates,
                                      join_var,
                                      summary_R = 2) {
  outcome_cols <- outcome_support_cols(my_outcome)
  base_required_cols <- unique(c(my_outcome, my_covariates))

  if (is.null(exp_var)) {
    mydf_ua <- dataset %>%
      dplyr::select(all_of(c(join_var, my_outcome, outcome_cols[-1], my_covariates))) %>%
      dplyr::ungroup() %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(base_required_cols), ~ !is.na(.)))

    if ("Study" %in% colnames(mydf_ua)) {
      table_n <- table(mydf_ua$Study)
      max_category <- names(table_n)[which.max(table_n)]

      mydf_ua <- mydf_ua %>%
        dplyr::mutate(Study = factor(Study)) %>%
        dplyr::mutate(Study = stats::relevel(as.factor(Study), ref = max_category))
    }

    if (nrow(mydf_ua) == 0) {
      return(NULL)
    }

    rhs_formula <- stats::as.formula(paste("~", paste0(my_covariates, collapse = "+")))

    if ("Study" %in% colnames(mydf_ua)) {
      tocheck <- my_covariates[!my_covariates %in% "Study"]
      clean_df <- mydf_ua %>%
        dplyr::group_by(Study) %>%
        dplyr::filter(dplyr::if_all(dplyr::all_of(tocheck), ~ dplyr::n_distinct(.) > 1)) %>%
        dplyr::ungroup()
    } else {
      clean_df <- mydf_ua
    }
  } else {
    mydf_ua <- dataset %>%
      dplyr::select(all_of(c(join_var, my_outcome, outcome_cols[-1], exp_var, my_covariates))) %>%
      dplyr::filter(dplyr::if_all(dplyr::all_of(unique(c(base_required_cols, exp_var))), ~ !is.na(.)))

    if (nrow(mydf_ua) == 0) {
      return(NULL)
    }

    table_n <- table(mydf_ua$Study)
    max_category <- names(table_n)[which.max(table_n)]

    mydf_ua <- mydf_ua %>%
      dplyr::mutate(Study = factor(Study)) %>%
      dplyr::mutate(Study = stats::relevel(as.factor(Study), ref = max_category))

    if (length(unique(mydf_ua[[exp_var]])) > 10) {
      mydf_ua <- mydf_ua %>%
        dplyr::group_by(!!rlang::sym(exp_var)) %>%
        dplyr::ungroup()
    } else {
      mydf_ua <- mydf_ua %>%
        dplyr::group_by(!!rlang::sym(exp_var)) %>%
        dplyr::filter(dplyr::n() >= 5) %>%
        dplyr::ungroup()
    }

    if (nrow(mydf_ua) == 0) {
      return(NULL)
    }

    clean_df <- mydf_ua %>%
      dplyr::group_by(Study) %>%
      dplyr::filter(dplyr::n_distinct(!!rlang::sym(exp_var)) > 1) %>%
      dplyr::ungroup() %>%
      droplevels()

    n_study <- length(unique(clean_df$Study))
    if (n_study == 1) {
      rhs_formula <- stats::as.formula(paste("~", exp_var))
    } else {
      rhs_formula <- stats::as.formula(paste("~", exp_var, "+", paste0(my_covariates, collapse = "+")))
    }
  }

  if (nrow(clean_df) == 0) {
    return(NULL)
  }

  df_cens <- build_df_cens(clean_df, my_outcome)
  if (nrow(df_cens) == 0) {
    return(NULL)
  }

  rhs_formula <- keep_varying_terms(df_cens, rhs_formula)
  if (length(attr(stats::terms(rhs_formula), "term.labels")) == 0) {
    return(NULL)
  }

  taus_use <- choose_censored_taus(df_cens)
  censoring_summary <- df_cens %>%
    dplyr::summarise(
      n = dplyr::n(),
      pct_detected = mean(detected, na.rm = TRUE),
      pct_censored = mean(!detected, na.rm = TRUE)
    )

  fit_obj <- fit_censored_quantile_models(
    df_cens = df_cens,
    rhs_formula = rhs_formula,
    taus_use = taus_use,
    summary_R = summary_R
  )

  if (is.null(fit_obj)) {
    return(NULL)
  }

  r1_all_df <- compute_censored_pseudo_r(
    df_cens = df_cens,
    fitted_full = fit_obj$fitted_full,
    fitted_null = fit_obj$fitted_null,
    taus_use = taus_use
  )
  r1_overall <- aggregate_pseudo_r(r1_all_df)

  ctable_final <- extract_crq_table(
    summary_obj = fit_obj$full_summary,
    outcome_col = my_outcome,
    taus_use = taus_use,
    exposure_col = exp_var
  ) %>%
    dplyr::left_join(r1_all_df, by = "Tau") %>%
    dplyr::mutate(
      Model_N = censoring_summary$n[[1]],
      Pct_Detected = censoring_summary$pct_detected[[1]],
      Pct_Censored = censoring_summary$pct_censored[[1]],
      Tau_Set_Used = paste(taus_use, collapse = ","),
      # Single goodness-of-fit number per model, pooled over the quantiles that
      # were actually estimable (2026-08-11). Null-loss weighted -- see
      # aggregate_pseudo_r(). Constant down the rows of a model, like Model_N.
      Pseudo_R_Overall      = r1_overall$overall,
      Pseudo_R_EqualWeight  = r1_overall$equal_weight,
      N_Taus_Estimable      = r1_overall$n_taus,
      Taus_Estimable        = r1_overall$taus,
      Pseudo_R_Basis        = r1_overall$basis
    )

  quantile_results <- list(
    df_cens,
    as.data.frame(fit_obj$fitted_full)
  )
  attr(quantile_results, "censoring_summary") <- censoring_summary %>%
    dplyr::mutate(Tau_Set_Used = paste(taus_use, collapse = ","))

  final_results <- list(ctable_final, quantile_results)
  names(final_results) <- c("Regression_Results", "Quantile_Results")
  final_results
}
