#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
  library(quantreg)
  library(survival)
})

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      next
    }
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    key <- parts[[1]]
    value <- paste(parts[-1], collapse = "=")
    out[[key]] <- value
  }
  out
}

`%||%` <- function(x, y) if (is.null(x) || length(x) == 0 || is.na(x) || x == "") y else x

find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir, winslash = "/", mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(cur, "Coding")) && dir.exists(file.path(cur, "Input"))) {
      return(cur)
    }
    parent <- dirname(cur)
    if (identical(parent, cur)) {
      break
    }
    cur <- parent
  }
  normalizePath(start_dir, winslash = "/", mustWork = TRUE)
}

to_bool <- function(x, default = FALSE) {
  if (is.null(x) || is.na(x) || x == "") {
    return(default)
  }
  tolower(as.character(x)) %in% c("true", "t", "1", "yes", "y")
}

to_num <- function(x, default = NA_real_) {
  if (is.null(x) || is.na(x) || x == "") {
    return(default)
  }
  suppressWarnings(as.numeric(x))
}

to_int <- function(x, default = NA_integer_) {
  val <- to_num(x, default = default)
  if (is.na(val)) default else as.integer(val)
}

parse_numeric_vector <- function(x, default) {
  if (is.null(x) || is.na(x) || x == "") {
    return(default)
  }
  parts <- trimws(unlist(strsplit(as.character(x), ",")))
  parts <- parts[nzchar(parts)]
  vals <- suppressWarnings(as.numeric(parts))
  vals[!is.na(vals)]
}

# std_study() and to_ng_per_mL() now come from harmonization_helpers.R
# (2026-08-11). The local to_ng_per_mL left mg/dL, ug/mL, ng/g, pmol/mL, nM and
# kg/m3 unconverted; the local std_study used a third regex variant.
# See CODE_AUDIT_2026-08-11.md sections 1.1 and 3.5.
source_harmonization_helpers <- function() {
  argv <- commandArgs(trailingOnly = FALSE)
  file_arg <- argv[grep("^--file=", argv)]
  here <- if (length(file_arg) > 0) {
    dirname(normalizePath(sub("^--file=", "", file_arg[[1]]), mustWork = FALSE))
  } else {
    getwd()
  }
  cands <- c(
    file.path(here, "harmonization_helpers.R"),
    "harmonization_helpers.R",
    "Coding/Paper1_Submission2/harmonization_helpers.R"
  )
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) {
    stop("Cannot locate harmonization_helpers.R. Looked in:\n  ",
         paste(cands, collapse = "\n  "))
  }
  source(hit[[1]])
}
source_harmonization_helpers()

coerce_predictors <- function(df, covariates) {
  out <- df
  for (nm in covariates) {
    if (is.character(out[[nm]])) {
      out[[nm]] <- as.factor(out[[nm]])
    }
  }
  droplevels(out)
}

varying_covariates <- function(df, covariates) {
  covariates[
    covariates %in% names(df) &
      vapply(
        covariates[covariates %in% names(df)],
        function(nm) dplyr::n_distinct(df[[nm]][!is.na(df[[nm]])]) > 1,
        logical(1)
      )
  ]
}

extract_rq_min_pvalue <- function(summary_obj) {
  if (is.null(summary_obj) || is.null(summary_obj$coefficients)) {
    return(NA_real_)
  }
  coeff_df <- as.data.frame(summary_obj$coefficients)
  if (!"Pr(>|t|)" %in% colnames(coeff_df)) {
    return(NA_real_)
  }
  coeff_df$Variable <- rownames(coeff_df)
  pvals <- coeff_df %>%
    dplyr::filter(.data$Variable != "(Intercept)") %>%
    dplyr::pull(`Pr(>|t|)`)
  if (length(pvals) == 0 || all(is.na(pvals))) {
    NA_real_
  } else {
    min(pvals, na.rm = TRUE)
  }
}

make_crq_predictions <- function(fit, rhs_formula, data, tau) {
  predicted_vals <- tryCatch(
    predict(fit, newdata = data),
    error = function(e) NULL
  )

  if (!is.null(predicted_vals)) {
    pred_mat <- as.matrix(predicted_vals)
    if (nrow(pred_mat) == nrow(data)) {
      if (ncol(pred_mat) == 1) {
        pred_mat <- matrix(pred_mat[, 1], ncol = 1)
      }
      colnames(pred_mat) <- as.character(tau)
      return(pred_mat)
    }
  }

  design_mat <- stats::model.matrix(rhs_formula, data = data)
  coef_raw <- stats::coef(fit, taus = tau)
  coef_mat <- matrix(
    coef_raw,
    nrow = ncol(design_mat),
    dimnames = list(colnames(design_mat), as.character(tau))
  )
  if (is.null(rownames(coef_mat))) {
    rownames(coef_mat) <- colnames(design_mat)
  }
  coef_lookup <- stats::setNames(seq_len(nrow(coef_mat)), rownames(coef_mat))
  aligned_coef <- matrix(
    0,
    nrow = ncol(design_mat),
    ncol = 1,
    dimnames = list(colnames(design_mat), as.character(tau))
  )
  matched_terms <- intersect(colnames(design_mat), names(coef_lookup))
  if (length(matched_terms) == 0) {
    stop("Coefficient matrix does not align with design matrix")
  }
  aligned_coef[matched_terms, 1] <- coef_mat[coef_lookup[matched_terms], 1]
  predictions <- design_mat %*% aligned_coef
  colnames(predictions) <- as.character(tau)
  predictions
}

extract_crq_min_pvalue <- function(summary_obj) {
  if (is.null(summary_obj)) {
    return(NA_real_)
  }
  summary_list <- if (inherits(summary_obj, "summary.crqs")) summary_obj else list(summary_obj)
  pvals <- unlist(
    lapply(
      summary_list,
      function(curr_summary) {
        coeff <- curr_summary$coefficients
        if (is.null(coeff)) {
          return(NA_real_)
        }
        coeff_df <- as.data.frame(coeff)
        if (!"Pr(>|t|)" %in% colnames(coeff_df)) {
          return(NA_real_)
        }
        coeff_df$Variable <- rownames(coeff_df)
        coeff_df %>%
          dplyr::filter(.data$Variable != "(Intercept)") %>%
          dplyr::pull(`Pr(>|t|)`)
      }
    ),
    use.names = FALSE
  )
  if (length(pvals) == 0 || all(is.na(pvals))) {
    NA_real_
  } else {
    min(pvals, na.rm = TRUE)
  }
}

build_df_cens <- function(df, outcome_col) {
  lod_col <- paste0(outcome_col, "_LOD")
  detected_col <- paste0(outcome_col, "_Detected")
  observed_col <- paste0(outcome_col, "_Obs")
  df %>%
    dplyr::mutate(
      exposure = .data[[outcome_col]],
      LOD = .data[[lod_col]],
      detected = dplyr::case_when(
        !is.na(.data[[detected_col]]) ~ as.logical(.data[[detected_col]]),
        !is.na(exposure) & !is.na(LOD) ~ exposure > LOD,
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

fit_censored_quantile_comparison <- function(df, outcome_col, covariates, tau, summary_R = 20L) {
  df_cens <- build_df_cens(df, outcome_col)
  if (nrow(df_cens) == 0) {
    return(list(full_fit = NULL, summary_obj = NULL, fitted_full = NULL, min_p = NA_real_, df_cens = df_cens))
  }

  covariates_use <- varying_covariates(df_cens, covariates)
  if (length(covariates_use) == 0) {
    return(list(full_fit = NULL, summary_obj = NULL, fitted_full = NULL, min_p = NA_real_, df_cens = df_cens))
  }

  df_cens <- coerce_predictors(df_cens, covariates_use)
  rhs_formula <- stats::reformulate(covariates_use)
  surv_formula <- stats::update(rhs_formula, survival::Surv(exposure_obs, detected, type = "left") ~ .)

  full_fit <- suppressWarnings(
    tryCatch(quantreg::crq(surv_formula, data = df_cens, taus = tau, method = "Portnoy"), error = function(e) NULL)
  )
  if (is.null(full_fit)) {
    return(list(full_fit = NULL, summary_obj = NULL, fitted_full = NULL, min_p = NA_real_, df_cens = df_cens))
  }

  summary_obj <- suppressWarnings(
    tryCatch(summary(full_fit, taus = tau, se = "boot", R = summary_R), error = function(e) NULL)
  )
  if (is.null(summary_obj)) {
    summary_obj <- suppressWarnings(
      tryCatch(summary(full_fit, taus = tau), error = function(e) NULL)
    )
  }

  fitted_full <- tryCatch(make_crq_predictions(full_fit, rhs_formula, df_cens, tau), error = function(e) NULL)
  min_p_use <- extract_crq_min_pvalue(summary_obj)

  if (is.na(min_p_use) || is.null(fitted_full)) {
    fallback_formula <- stats::reformulate(covariates_use, response = "exposure_obs")
    fallback_fit <- suppressWarnings(
      tryCatch(quantreg::rq(fallback_formula, tau = tau, data = df_cens), error = function(e) NULL)
    )
    if (!is.null(fallback_fit)) {
      fallback_summary <- suppressWarnings(
        tryCatch(summary(fallback_fit, se = "nid"), error = function(e) NULL)
      )
      fallback_pred <- suppressWarnings(
        tryCatch(as.matrix(stats::predict(fallback_fit, newdata = df_cens)), error = function(e) NULL)
      )
      if (!is.null(fallback_pred)) {
        if (is.null(dim(fallback_pred))) {
          fallback_pred <- matrix(fallback_pred, ncol = 1)
        }
        colnames(fallback_pred) <- as.character(tau)
        fitted_full <- fallback_pred
      }
      if (is.null(summary_obj)) {
        summary_obj <- fallback_summary
      }
      min_p_fallback <- extract_rq_min_pvalue(fallback_summary)
      if (!is.na(min_p_fallback)) {
        min_p_use <- min_p_fallback
      }
    }
  }

  list(
    full_fit = full_fit,
    summary_obj = summary_obj,
    fitted_full = fitted_full,
    min_p = min_p_use,
    df_cens = df_cens,
    covariates_use = covariates_use
  )
}

estimate_censored_power_one_tau <- function(df, outcome_col, covariates, tau, nsim, alpha, summary_R = 20L) {
  base_fit <- fit_censored_quantile_comparison(df, outcome_col, covariates, tau, summary_R = summary_R)

  if (is.null(base_fit$full_fit) || is.na(base_fit$min_p) || is.null(base_fit$fitted_full)) {
    return(tibble(
      tau = tau,
      observed_p = NA_real_,
      power = NA_real_,
      mcse = NA_real_,
      successful_sims = 0L
    ))
  }

  df_cens <- base_fit$df_cens
  fitted_vals <- as.numeric(base_fit$fitted_full[, 1])
  residuals_tau <- df_cens$exposure_obs - fitted_vals

  sim_pvalues <- rep(NA_real_, nsim)
  for (b in seq_len(nsim)) {
    latent_sim <- fitted_vals + sample(residuals_tau, size = nrow(df_cens), replace = TRUE)
    latent_sim <- pmax(latent_sim, 0)
    detected_sim <- dplyr::case_when(
      is.na(df_cens$LOD) ~ TRUE,
      TRUE ~ latent_sim > df_cens$LOD
    )
    exposure_obs_sim <- dplyr::if_else(detected_sim, latent_sim, df_cens$LOD)

    sim_df <- df
    sim_df[[outcome_col]] <- exposure_obs_sim
    sim_df[[paste0(outcome_col, "_Obs")]] <- exposure_obs_sim
    sim_df[[paste0(outcome_col, "_Detected")]] <- detected_sim

    sim_fit <- fit_censored_quantile_comparison(sim_df, outcome_col, covariates, tau, summary_R = summary_R)
    sim_pvalues[[b]] <- sim_fit$min_p
  }

  successful_sims <- sum(!is.na(sim_pvalues))

  # Same upward-bias fix as Targeted_Quantile_Power_Array.R (2026-08-11,
  # CODE_AUDIT_2026-08-11.md section 1.7): failed replicates were dropped from
  # the denominator, and failures concentrate where power is low.
  n_reject <- sum(sim_pvalues < alpha, na.rm = TRUE)
  power_est <- if (nsim > 0) n_reject / nsim else NA_real_
  power_completed <- if (successful_sims > 0) n_reject / successful_sims else NA_real_
  failure_rate <- if (nsim > 0) (nsim - successful_sims) / nsim else NA_real_
  mcse_est <- if (nsim > 0) sqrt(power_est * (1 - power_est) / nsim) else NA_real_

  tibble(
    tau = tau,
    observed_p = base_fit$min_p,
    power = power_est,
    power_completed = power_completed,
    failure_rate = failure_rate,
    mcse = mcse_est,
    successful_sims = successful_sims
  )
}

prepare_targeted_lists <- function(targeted_df) {
  child_data <- build_censored_targeted_panel(
    targeted_df = targeted_df,
    id_var = "Child_PID",
    unit_converter = to_ng_per_mL
  )
  child_analytes <- sort(unique(child_data$Analyte_Code))
  child_list <- build_censored_targeted_list(
    panel_df = child_data,
    id_var = "Child_PID",
    analytes = child_analytes
  )

  mat_data <- build_censored_targeted_panel(
    targeted_df = targeted_df,
    id_var = "Mat_PID",
    unit_converter = to_ng_per_mL
  )
  mat_analytes <- sort(unique(mat_data$Analyte_Code))
  mat_list <- build_censored_targeted_list(
    panel_df = mat_data,
    id_var = "Mat_PID",
    analytes = mat_analytes
  )

  list(
    child_data = child_data,
    child_analytes = child_analytes,
    child_list = child_list,
    mat_data = mat_data,
    mat_analytes = mat_analytes,
    mat_list = mat_list
  )
}

build_sdoh_panels <- function(sdoh_df, targeted_df) {
  exwas_sdoh <- sdoh_df %>%
    mutate(
      Edu_Any = dplyr::coalesce(Max_Edu, Mat_Edu),
      Mat_Child_Comb_Race = dplyr::coalesce(Child_Race, Mat_Race),
      Mat_Age = suppressWarnings(as.numeric(Mat_Age)),
      Location = stats::relevel(as.factor(Location), ref = "United_States"),
      Location_Country = dplyr::case_when(
        Location %in% c(
          "United_States", "NorthEast_US", "New_Hampshire", "Washington_State",
          "Midwestern_US", "Denver", "Baltimore", "Sacramento",
          "Southern_California", "Syracuse", "Michigan", "California",
          "SanFrancisco", "Massachusetts"
        ) ~ "United_States",
        TRUE ~ "Worldwide"
      ),
      Mat_Edu_Binary = dplyr::case_when(
        Mat_Edu %in% c("No_Education", "Primary_Education", "Primary_School", "High_School") ~ "No_College",
        Mat_Edu %in% c("Some_College", "Four_Year_College_Degree", "Some_Grad_School") ~ "Some_College",
        Mat_Edu == "Not_Reported" ~ NA_character_,
        TRUE ~ NA_character_
      ),
      Mat_Income_Binary = dplyr::case_when(
        Mat_Income %in% c("0_25k", "25_49k") ~ "Under_49k",
        Mat_Income %in% c("50_74k", "75_99k", "100_124k", "Over_100k", "Over_125k", "Over_49k") ~ "Over_49k",
        TRUE ~ NA_character_
      )
    )

  child_sample_yr <- targeted_df %>%
    filter(!is.na(Child_PID)) %>%
    group_by(Child_PID) %>%
    summarise(Sample_Collection_Year = median_or_na(Sample_Collection_Year), .groups = "drop")

  mat_sample_yr <- targeted_df %>%
    filter(!is.na(Mat_PID)) %>%
    group_by(Mat_PID) %>%
    summarise(Sample_Collection_Year = median_or_na(Sample_Collection_Year), .groups = "drop")

  child_sdoh <- exwas_sdoh %>%
    filter(!is.na(Child_PID)) %>%
    select(Study, Child_PID, Mat_PID, everything()) %>%
    left_join(child_sample_yr, by = "Child_PID") %>%
    mutate(Sample_Collection_Year = as.numeric(Sample_Collection_Year))

  mat_sdoh <- exwas_sdoh %>%
    filter(!is.na(Mat_PID)) %>%
    select(Study, Child_PID, Mat_PID, everything()) %>%
    left_join(mat_sample_yr, by = "Mat_PID") %>%
    mutate(Sample_Collection_Year = as.numeric(Sample_Collection_Year))

  child_ads <- child_sdoh %>%
    select(Child_PID, Study, Location, Location_Region, Year, Center, Child_Age, Child_Sex)

  mat_ads <- mat_sdoh %>%
    select(Mat_PID, Study, Location, Location_Region, Year, Center, Mat_Age)

  list(
    child_ads = child_ads,
    mat_ads = mat_ads
  )
}

child_model_specs <- list(
  Study = "Study",
  Year = "Year",
  Location = "Location",
  Location_Region = "Location_Region",
  Center = "Center",
  Child_Age = "Child_Age",
  Child_Sex = "Child_Sex",
  Child_Age_Child_Sex_Study = c("Child_Age", "Child_Sex", "Study"),
  Location_Center = c("Location", "Center"),
  Location_Year = c("Location", "Year"),
  Location_Study = c("Location", "Study")
)

mat_model_specs <- list(
  Study = "Study",
  Year = "Year",
  Location = "Location",
  Location_Region = "Location_Region",
  Center = "Center",
  Mat_Age = "Mat_Age",
  Location_Center = c("Location", "Center"),
  Location_Year = c("Location", "Year"),
  Location_Study = c("Location", "Study")
)

build_between_study_dataset <- function(exposure_data,
                                        outcome_data,
                                        join_var,
                                        analyte,
                                        covariates,
                                        min_n) {
  curr_exposure <- exposure_data %>%
    mutate(Study = str_remove(Study, "[_]\\w+|:")) %>%
    select(all_of(c(join_var, "Study", covariates))) %>%
    distinct()

  curr_outcome <- outcome_data[[analyte]] %>%
    mutate(Study = str_remove(Study, "[_]\\w+|:")) %>%
    select_outcome_data(join_var = join_var, outcome_col = analyte) %>%
    filter(!is.na(.data[[analyte]]))

  analysis_df <- left_join(curr_outcome, curr_exposure, by = c(join_var, "Study")) %>%
    mutate(!!analyte := as.numeric(.data[[analyte]])) %>%
    distinct(.keep_all = TRUE)

  required_cols <- unique(c(join_var, "Study", analyte, covariates, paste0(analyte, c("_LOD", "_Detected", "_ExposureObs"))))
  keep_cols <- required_cols[required_cols %in% names(analysis_df)]
  analysis_df <- analysis_df %>%
    select(all_of(keep_cols)) %>%
    drop_na(all_of(c(join_var, "Study", analyte)))

  available_covariates <- varying_covariates(analysis_df, covariates)
  if (length(available_covariates) == 0) {
    return(list(status = "no_varying_covariates", n = nrow(analysis_df), covariates = character(0), n_studies = n_distinct(analysis_df$Study), data = NULL))
  }

  analysis_df <- analysis_df %>%
    drop_na(all_of(available_covariates)) %>%
    coerce_predictors(available_covariates)

  available_covariates <- varying_covariates(analysis_df, available_covariates)
  if (length(available_covariates) == 0) {
    return(list(status = "no_varying_covariates", n = nrow(analysis_df), covariates = character(0), n_studies = n_distinct(analysis_df$Study), data = NULL))
  }

  if (nrow(analysis_df) < min_n) {
    return(list(status = "below_min_n", n = nrow(analysis_df), covariates = available_covariates, n_studies = n_distinct(analysis_df$Study), data = NULL))
  }

  if (n_distinct(analysis_df$Study) < 2) {
    return(list(status = "fewer_than_two_studies", n = nrow(analysis_df), covariates = available_covariates, n_studies = n_distinct(analysis_df$Study), data = NULL))
  }

  list(
    status = "ok",
    n = nrow(analysis_df),
    covariates = available_covariates,
    n_studies = n_distinct(analysis_df$Study),
    data = analysis_df
  )
}

run_task_power <- function(task_row,
                           child_ads,
                           mat_ads,
                           child_list,
                           mat_list,
                           taus,
                           nsim,
                           alpha,
                           min_n,
                           summary_R) {
  if (task_row$panel == "child") {
    exposure_data <- child_ads
    outcome_data <- child_list
    join_var <- "Child_PID"
    model_specs <- child_model_specs
  } else {
    exposure_data <- mat_ads
    outcome_data <- mat_list
    join_var <- "Mat_PID"
    model_specs <- mat_model_specs
  }

  covariates <- model_specs[[task_row$model_name]]
  prep <- build_between_study_dataset(
    exposure_data = exposure_data,
    outcome_data = outcome_data,
    join_var = join_var,
    analyte = task_row$analyte,
    covariates = covariates,
    min_n = min_n
  )

  screening <- tibble(
    task_id = task_row$task_id,
    panel = task_row$panel,
    analysis = "between_study",
    model_name = task_row$model_name,
    analyte = task_row$analyte,
    status = prep$status,
    n = prep$n,
    n_studies = prep$n_studies,
    covariates_used = paste(prep$covariates, collapse = ";")
  )

  if (!identical(prep$status, "ok")) {
    return(list(screening = screening, power = tibble()))
  }

  power <- purrr::map_dfr(
    taus,
    function(curr_tau) {
      estimate_censored_power_one_tau(
        df = prep$data,
        outcome_col = task_row$analyte,
        covariates = prep$covariates,
        tau = curr_tau,
        nsim = nsim,
        alpha = alpha,
        summary_R = summary_R
      ) %>%
        mutate(
          task_id = task_row$task_id,
          panel = task_row$panel,
          analysis = "between_study",
          model_name = task_row$model_name,
          analyte = task_row$analyte,
          n = prep$n,
          n_studies = prep$n_studies,
          covariates_used = paste(prep$covariates, collapse = ";")
        ) %>%
        select(task_id, panel, analysis, model_name, analyte, tau, n, n_studies, covariates_used, observed_p, power, mcse, successful_sims)
    }
  )

  list(screening = screening, power = power)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- script_arg[grep("^--file=", script_arg)]
script_file <- if (length(script_file) > 0) sub("^--file=", "", script_file[[1]]) else "."
# Rscript encodes spaces in the --file= argument as "~+~". The cluster path has
# no spaces so this never mattered there, but the local checkout lives under
# "My Drive" and the encoded path fails to resolve.
script_file <- gsub("~+~", " ", script_file, fixed = TRUE)
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
source(file.path(script_dir, "censored_quantile_helpers.R"))
project_root_default <- if (dir.exists(Sys.getenv("HHEAR_ROOT", path.expand("~/HHEAR/HHEAR_Quantreg")))) {
  Sys.getenv("HHEAR_ROOT", path.expand("~/HHEAR/HHEAR_Quantreg"))
} else {
  file.path(script_dir, "..", "..")
}
project_root <- find_project_root(args$project_root %||% project_root_default)

mode <- args$mode %||% "run_task"
targeted_rdata <- args$targeted_rdata %||% file.path(project_root, "Input", "Harmonized_Targeted.RData")
sdoh_rdata <- args$sdoh_rdata %||% file.path(project_root, "Input", "Harmonized_SDOH.RData")
out_dir <- args$out_dir %||% file.path(project_root, "Output", "Paper1_Output_Submission2", "Between_Study_Power_Output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "task_results"), recursive = TRUE, showWarnings = FALSE)

taus <- parse_numeric_vector(args$taus, c(0.05, 0.25, 0.5, 0.75, 0.95))
nsim <- to_int(args$nsim, 200L)
alpha <- to_num(args$alpha, 0.05)
min_n <- to_int(args$min_n, 60L)
seed <- to_int(args$seed, 20260511L)
summary_R <- to_int(args$summary_R, 20L)
run_child <- to_bool(args$run_child, TRUE)
run_maternal <- to_bool(args$run_maternal, TRUE)
max_analytes <- to_int(args$max_analytes, NA_integer_)
task_id <- to_int(args$task_id, to_int(Sys.getenv("SLURM_ARRAY_TASK_ID"), NA_integer_))

# Per-task seed (2026-08-11). See Targeted_Quantile_Power_Array.R for the full
# note: a bare set.seed(seed) gave every SLURM array task the identical random
# stream, so Monte Carlo errors were perfectly correlated across the array.
task_seed <- if (!is.na(task_id)) seed + as.integer(task_id) else seed
set.seed(task_seed)
message(sprintf("[seed] base=%d task_id=%s -> task_seed=%d",
                seed, ifelse(is.na(task_id), "NA", task_id), task_seed))

load(targeted_rdata)
load(sdoh_rdata)
stopifnot(exists("Targeted_Data_All"), exists("SDOH_Data_All"))

targeted_objs <- prepare_targeted_lists(Targeted_Data_All)
sdoh_panels <- build_sdoh_panels(SDOH_Data_All, Targeted_Data_All)

make_task_manifest <- function(panel_name, panel_df, analytes, model_specs) {
  analyte_study_counts <- panel_df %>%
    mutate(Study = std_study(Study)) %>%
    distinct(Study, Analyte_Code) %>%
    count(Analyte_Code, name = "n_studies") %>%
    filter(n_studies >= 2) %>%
    arrange(desc(n_studies), Analyte_Code)

  analytes_keep <- analyte_study_counts$Analyte_Code
  if (!is.na(max_analytes)) {
    analytes_keep <- head(analytes_keep, max_analytes)
  }

  tidyr::expand_grid(
    panel = panel_name,
    model_name = names(model_specs),
    analyte = analytes_keep
  )
}

manifest <- bind_rows(
  if (run_child) make_task_manifest("child", targeted_objs$child_data, targeted_objs$child_analytes, child_model_specs) else tibble(),
  if (run_maternal) make_task_manifest("maternal", targeted_objs$mat_data, targeted_objs$mat_analytes, mat_model_specs) else tibble()
) %>%
  mutate(task_id = row_number()) %>%
  select(task_id, panel, model_name, analyte)

manifest_file <- file.path(out_dir, "task_manifest.csv")
readr::write_csv(manifest, manifest_file)

if (mode == "build_manifest") {
  cat("Manifest written to:", manifest_file, "\n")
  cat("Task count:", nrow(manifest), "\n")
  quit(save = "no", status = 0)
}

if (mode == "run_task") {
  if (is.na(task_id) || task_id < 1L || task_id > nrow(manifest)) {
    stop("A valid task_id or SLURM_ARRAY_TASK_ID is required for mode=run_task.")
  }

  task_row <- manifest[task_id, , drop = FALSE]
  result <- run_task_power(
    task_row = task_row,
    child_ads = sdoh_panels$child_ads,
    mat_ads = sdoh_panels$mat_ads,
    child_list = targeted_objs$child_list,
    mat_list = targeted_objs$mat_list,
    taus = taus,
    nsim = nsim,
    alpha = alpha,
    min_n = min_n,
    summary_R = summary_R
  )

  saveRDS(
    list(
      task = task_row,
      screening = result$screening,
      power = result$power,
      config = list(
        taus = taus,
        nsim = nsim,
        alpha = alpha,
        min_n = min_n,
        seed = seed,
        task_seed = task_seed,   # recorded 2026-08-11 so a task is reproducible alone
        summary_R = summary_R
      )
    ),
    file.path(out_dir, "task_results", sprintf("task_%06d.rds", task_id))
  )

  cat("Completed task", task_id, "of", nrow(manifest), "\n")
  quit(save = "no", status = 0)
}

if (mode == "combine_results") {
  task_files <- list.files(file.path(out_dir, "task_results"), pattern = "^task_[0-9]{6}\\.rds$", full.names = TRUE)
  if (length(task_files) == 0) {
    stop("No task result files found to combine.")
  }

  task_results <- lapply(task_files, readRDS)
  screening_results <- bind_rows(lapply(task_results, `[[`, "screening"))
  power_results <- bind_rows(lapply(task_results, `[[`, "power"))
  task_table <- bind_rows(lapply(task_results, `[[`, "task"))

  readr::write_csv(screening_results, file.path(out_dir, "between_study_quantile_power_screening.csv"))
  readr::write_csv(power_results, file.path(out_dir, "between_study_quantile_power_results.csv"))
  readr::write_csv(task_table, file.path(out_dir, "between_study_quantile_power_completed_tasks.csv"))
  saveRDS(
    list(
      manifest = manifest,
      screening_results = screening_results,
      power_results = power_results,
      completed_tasks = task_table
    ),
    file.path(out_dir, "between_study_quantile_power_results.rds")
  )

  cat("Combined", length(task_files), "task files.\n")
  quit(save = "no", status = 0)
}

stop("Unknown mode: ", mode)
