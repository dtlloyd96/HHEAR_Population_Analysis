#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(stringr)
  library(readr)
  library(tibble)
  library(quantreg)
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
  if (is.na(val)) {
    default
  } else {
    as.integer(val)
  }
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
# kg/m3 unconverted; the local std_study used a third regex variant. See
# CODE_AUDIT_2026-08-11.md sections 1.1 and 3.5.
source_harmonization_helpers <- function() {
  # Look next to this script first (works under Rscript and SLURM), then in the
  # usual relative locations for an interactive session.
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

median_or_na <- function(x) {
  x <- suppressWarnings(as.numeric(x))
  x <- x[is.finite(x)]
  if (length(x) == 0) {
    return(NA_real_)
  }
  stats::median(x)
}

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

# A factor level carrying only one or two observations leaves the rank-deficient
# covariance in summary.rq singular, and anova() then aborts with "singular
# matrix in 'backsolve'". The whole model returns NO p-value rather than a weak
# one, so the analyte silently contributes nothing to the power estimate. Drop
# rows sitting in levels below min_level_n; in practice this costs a handful of
# rows out of several hundred and makes the model estimable again.
MIN_LEVEL_N <- 5L

drop_rare_levels <- function(df, covariates, min_level_n = MIN_LEVEL_N) {
  for (v in covariates) {
    if (!v %in% names(df) || !is.factor(df[[v]])) next
    tt <- table(df[[v]])
    keep <- names(tt)[tt >= min_level_n]
    df <- df[as.character(df[[v]]) %in% keep, , drop = FALSE]
    df[[v]] <- droplevels(factor(df[[v]]))
  }
  df
}

extract_global_pvalue <- function(anova_obj) {
  if (is.null(anova_obj) || is.null(anova_obj$table) || nrow(anova_obj$table) == 0) {
    return(NA_real_)
  }
  as.numeric(anova_obj$table$pvalue[[1]])
}

fit_quantile_comparison <- function(df, outcome_col, covariates, tau) {
  full_formula <- stats::reformulate(covariates, response = outcome_col)
  null_formula <- stats::as.formula(paste(outcome_col, "~ 1"))

  fit_full <- suppressWarnings(
    tryCatch(quantreg::rq(full_formula, tau = tau, data = df), error = function(e) NULL)
  )
  fit_null <- suppressWarnings(
    tryCatch(quantreg::rq(null_formula, tau = tau, data = df), error = function(e) NULL)
  )

  if (is.null(fit_full) || is.null(fit_null)) {
    return(list(fit_full = fit_full, fit_null = fit_null, global_p = NA_real_))
  }

  anova_obj <- tryCatch(anova(fit_full, fit_null), error = function(e) NULL)

  list(
    fit_full = fit_full,
    fit_null = fit_null,
    global_p = extract_global_pvalue(anova_obj)
  )
}

estimate_power_one_tau <- function(df, outcome_col, covariates, tau, nsim, alpha, nonnegative_outcome = TRUE) {
  base_fit <- fit_quantile_comparison(df, outcome_col, covariates, tau)

  if (is.null(base_fit$fit_full) || is.na(base_fit$global_p)) {
    # Say WHY there is no estimate. The dominant cause is a censoring mass point:
    # when most of the outcome sits on one value (everything below the LOD), the
    # quantiles at or below that mass are not estimable by plain rq, and anova()
    # fails on a singular density estimate. That is a real limit of the method,
    # not a data error -- but it must not look like a silent zero.
    y <- df[[outcome_col]]
    modal_share <- if (length(y)) max(table(y)) / length(y) else NA_real_
    reason <- if (is.null(base_fit$fit_full)) "model_did_not_converge"
              else if (!is.na(modal_share) && modal_share >= tau) "quantile_below_censoring_mass"
              else "anova_singular"
    return(tibble(
      tau = tau,
      observed_p = NA_real_,
      power = NA_real_,
      power_completed = NA_real_,
      failure_rate = NA_real_,
      mcse = NA_real_,
      successful_sims = 0L,
      skip_reason = reason
    ))
  }

  fitted_vals <- as.numeric(stats::predict(base_fit$fit_full, newdata = df))
  residuals_tau <- df[[outcome_col]] - fitted_vals

  sim_pvalues <- rep(NA_real_, nsim)
  for (b in seq_len(nsim)) {
    sim_df <- df
    sim_df[[outcome_col]] <- fitted_vals + sample(residuals_tau, size = nrow(df), replace = TRUE)
    if (nonnegative_outcome) {
      sim_df[[outcome_col]] <- pmax(sim_df[[outcome_col]], 0)
    }
    sim_fit <- fit_quantile_comparison(sim_df, outcome_col, covariates, tau)
    sim_pvalues[[b]] <- sim_fit$global_p
  }

  successful_sims <- sum(!is.na(sim_pvalues))

  # POWER ESTIMATE BIAS (fixed 2026-08-11, CODE_AUDIT_2026-08-11.md section 1.7).
  #
  # This was `mean(sim_pvalues < alpha, na.rm = TRUE)`, which drops failed
  # replicates from the DENOMINATOR. A replicate fails when the model does not
  # converge or anova() hits a singular density -- i.e. disproportionately when
  # the signal is weak, which is exactly where power is low. Dropping them
  # therefore biases the estimate UPWARD.
  #
  # Reported both ways now:
  #   power           - failures counted as non-rejections (conservative, and the
  #                     honest answer to "will this analysis detect the effect?",
  #                     since a failed fit does not yield a detection)
  #   power_completed - the old among-successes rate, kept for comparability
  # If failure_rate is non-trivial the two diverge and the run needs inspecting.
  n_reject <- sum(sim_pvalues < alpha, na.rm = TRUE)
  power_est <- if (nsim > 0) n_reject / nsim else NA_real_
  power_completed <- if (successful_sims > 0) n_reject / successful_sims else NA_real_
  failure_rate <- if (nsim > 0) (nsim - successful_sims) / nsim else NA_real_

  mcse_est <- if (nsim > 0) sqrt(power_est * (1 - power_est) / nsim) else NA_real_

  tibble(
    tau = tau,
    observed_p = base_fit$global_p,
    power = power_est,
    power_completed = power_completed,
    failure_rate = failure_rate,
    mcse = mcse_est,
    successful_sims = successful_sims,
    skip_reason = NA_character_
  )
}

build_sdoh_panel <- function(sdoh_df, targeted_df, id_var, panel) {
  year_df <- targeted_df %>%
    filter(!is.na(.data[[id_var]])) %>%
    mutate(
      Study = std_study(Study),
      Sample_Collection_Year = suppressWarnings(as.numeric(Sample_Collection_Year))
    ) %>%
    group_by(Study, .data[[id_var]]) %>%
    summarise(Sample_Collection_Year = median_or_na(Sample_Collection_Year), .groups = "drop")

  base_sdoh <- sdoh_df %>%
    mutate(
      Study = std_study(Study),
      Edu_Any = dplyr::coalesce(Max_Edu, Mat_Edu),
      Mat_Child_Comb_Race = dplyr::coalesce(Child_Race, Mat_Race),
      Mat_Age = suppressWarnings(as.numeric(Mat_Age)),
      Child_Age = suppressWarnings(as.numeric(Child_Age))
    )

  if (panel == "child") {
    base_sdoh <- base_sdoh %>%
      filter(!is.na(Child_PID)) %>%
      distinct(Study, Child_PID, .keep_all = TRUE) %>%
      left_join(year_df, by = c("Study", "Child_PID"))
  } else {
    base_sdoh <- base_sdoh %>%
      filter(!is.na(Mat_PID)) %>%
      distinct(Study, Mat_PID, .keep_all = TRUE) %>%
      left_join(year_df, by = c("Study", "Mat_PID"))
  }

  base_sdoh
}

build_targeted_panel <- function(targeted_df, id_var, use_log1p = TRUE) {
  targeted_df %>%
    filter(!is.na(.data[[id_var]]), !is.na(Analyte_Code)) %>%
    mutate(
      Study = std_study(Study),
      Concentration_ng_mL = to_ng_per_mL(Concentration, Units)
    ) %>%
    group_by(Study, .data[[id_var]], Analyte_Code) %>%
    summarise(Concentration_ng_mL = mean(Concentration_ng_mL, na.rm = TRUE), .groups = "drop") %>%
    mutate(
      Exposure = if (use_log1p) log1p(pmax(Concentration_ng_mL, 0)) else Concentration_ng_mL
    )
}

build_within_study_dataset <- function(panel_sdoh, panel_targeted, study, analyte, id_var, preferred_covariates, min_n) {
  exp_df <- panel_targeted %>%
    filter(Study == study, Analyte_Code == analyte) %>%
    select(Study, all_of(id_var), Analyte_Code, Exposure, Concentration_ng_mL)

  if (nrow(exp_df) == 0) {
    return(list(status = "no_targeted_rows", n = 0L, n_studies = 1L, covariates = character(0), data = NULL))
  }

  analysis_df <- panel_sdoh %>%
    filter(Study == study) %>%
    select(Study, all_of(id_var), any_of(preferred_covariates)) %>%
    inner_join(exp_df, by = c("Study", id_var)) %>%
    distinct()

  available_covariates <- varying_covariates(analysis_df, preferred_covariates)
  if (length(available_covariates) == 0) {
    return(list(status = "no_varying_covariates", n = nrow(analysis_df), n_studies = 1L, covariates = character(0), data = NULL))
  }

  analysis_df <- analysis_df %>%
    select(Study, all_of(id_var), all_of(available_covariates), Exposure, Concentration_ng_mL) %>%
    tidyr::drop_na()

  available_covariates <- varying_covariates(analysis_df, available_covariates)
  analysis_df <- coerce_predictors(analysis_df, available_covariates)

  # Prune near-empty factor levels, then re-check that each covariate still
  # varies -- pruning can collapse a two-level factor down to one.
  analysis_df <- drop_rare_levels(analysis_df, available_covariates)
  available_covariates <- varying_covariates(analysis_df, available_covariates)
  if (length(available_covariates) == 0) {
    return(list(status = "no_varying_covariates_after_pruning", n = nrow(analysis_df),
                n_studies = 1L, covariates = character(0), data = NULL))
  }

  if (nrow(analysis_df) < min_n) {
    return(list(status = "below_min_n", n = nrow(analysis_df), n_studies = 1L, covariates = available_covariates, data = NULL))
  }
  if (dplyr::n_distinct(analysis_df$Exposure) < 5) {
    return(list(status = "too_little_outcome_variation", n = nrow(analysis_df), n_studies = 1L, covariates = available_covariates, data = NULL))
  }

  list(status = "ok", n = nrow(analysis_df), n_studies = 1L, covariates = available_covariates, data = analysis_df)
}

build_pooled_study_dataset <- function(panel_targeted, analyte, id_var, min_n) {
  analysis_df <- panel_targeted %>%
    filter(Analyte_Code == analyte) %>%
    select(Study, all_of(id_var), Analyte_Code, Exposure, Concentration_ng_mL) %>%
    distinct() %>%
    tidyr::drop_na(Study, Exposure)

  if (nrow(analysis_df) == 0) {
    return(list(status = "no_targeted_rows", n = 0L, n_studies = 0L, covariates = character(0), data = NULL))
  }

  n_studies <- dplyr::n_distinct(analysis_df$Study)
  if (n_studies < 2) {
    return(list(status = "fewer_than_two_studies", n = nrow(analysis_df), n_studies = n_studies, covariates = "Study", data = NULL))
  }
  if (nrow(analysis_df) < min_n) {
    return(list(status = "below_min_n", n = nrow(analysis_df), n_studies = n_studies, covariates = "Study", data = NULL))
  }
  if (dplyr::n_distinct(analysis_df$Exposure) < 5) {
    return(list(status = "too_little_outcome_variation", n = nrow(analysis_df), n_studies = n_studies, covariates = "Study", data = NULL))
  }

  analysis_df <- analysis_df %>% mutate(Study = droplevels(as.factor(Study)))
  list(status = "ok", n = nrow(analysis_df), n_studies = n_studies, covariates = "Study", data = analysis_df)
}

run_task_power <- function(task_row, child_sdoh, maternal_sdoh, child_targeted, maternal_targeted, child_covariates, maternal_covariates, taus, nsim, alpha, min_n) {
  if (task_row$panel == "child") {
    panel_sdoh <- child_sdoh
    panel_targeted <- child_targeted
    id_var <- "Child_PID"
    preferred_covariates <- child_covariates
  } else {
    panel_sdoh <- maternal_sdoh
    panel_targeted <- maternal_targeted
    id_var <- "Mat_PID"
    preferred_covariates <- maternal_covariates
  }

  if (task_row$analysis == "within_study_demographics") {
    prep <- build_within_study_dataset(
      panel_sdoh = panel_sdoh,
      panel_targeted = panel_targeted,
      study = task_row$study,
      analyte = task_row$analyte,
      id_var = id_var,
      preferred_covariates = preferred_covariates,
      min_n = min_n
    )
  } else {
    prep <- build_pooled_study_dataset(
      panel_targeted = panel_targeted %>% filter(Study %in% task_row$study_set),
      analyte = task_row$analyte,
      id_var = id_var,
      min_n = min_n
    )
  }

  screening <- tibble(
    task_id = task_row$task_id,
    panel = task_row$panel,
    analysis = task_row$analysis,
    study = task_row$study,
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
      estimate_power_one_tau(
        df = prep$data,
        outcome_col = "Exposure",
        covariates = prep$covariates,
        tau = curr_tau,
        nsim = nsim,
        alpha = alpha,
        nonnegative_outcome = TRUE
      ) %>%
        mutate(
          task_id = task_row$task_id,
          panel = task_row$panel,
          analysis = task_row$analysis,
          study = task_row$study,
          analyte = task_row$analyte,
          n = prep$n,
          n_studies = prep$n_studies,
          covariates_used = paste(prep$covariates, collapse = ";")
        ) %>%
        select(task_id, panel, analysis, study, analyte, tau, n, n_studies, covariates_used,
               observed_p, power, power_completed, failure_rate, mcse, successful_sims, skip_reason)
    }
  )

  list(screening = screening, power = power)
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

script_arg <- commandArgs(trailingOnly = FALSE)
script_file <- script_arg[grep("^--file=", script_arg)]
script_file <- if (length(script_file) > 0) sub("^--file=", "", script_file[[1]]) else "."
script_dir <- dirname(normalizePath(script_file, winslash = "/", mustWork = FALSE))
project_root_default <- if (dir.exists(Sys.getenv("HHEAR_ROOT", path.expand("~/HHEAR/HHEAR_Quantreg")))) {
  Sys.getenv("HHEAR_ROOT", path.expand("~/HHEAR/HHEAR_Quantreg"))
} else {
  file.path(script_dir, "..", "..")
}

project_root <- find_project_root(args$project_root %||% project_root_default)

mode <- args$mode %||% "run_task"
targeted_rdata <- args$targeted_rdata %||% file.path(project_root, "Input", "Harmonized_Targeted.RData")
sdoh_rdata <- args$sdoh_rdata %||% file.path(project_root, "Input", "Harmonized_SDOH.RData")
out_dir <- args$out_dir %||% file.path(project_root, "Output", "Power_Output")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(out_dir, "task_results"), recursive = TRUE, showWarnings = FALSE)

taus <- parse_numeric_vector(args$taus, c(0.05, 0.25, 0.5, 0.75, 0.95))
nsim <- to_int(args$nsim, 200L)
alpha <- to_num(args$alpha, 0.05)
min_n <- to_int(args$min_n, 60L)
seed <- to_int(args$seed, 20260504L)
use_log1p <- to_bool(args$use_log1p, TRUE)
run_child <- to_bool(args$run_child, TRUE)
run_maternal <- to_bool(args$run_maternal, TRUE)
run_pooled_study_analysis <- to_bool(args$run_pooled_study_analysis, TRUE)
max_studies <- to_int(args$max_studies, NA_integer_)
max_analytes <- to_int(args$max_analytes, NA_integer_)
task_id <- to_int(args$task_id, to_int(Sys.getenv("SLURM_ARRAY_TASK_ID"), NA_integer_))

# SEED INDEPENDENCE (fixed 2026-08-11, CODE_AUDIT_2026-08-11.md section 1.7).
#
# This was a bare set.seed(seed) with a constant default. Each SLURM array task
# is a separate R process that sets the same seed and then runs exactly one task,
# so EVERY task drew the identical random stream: the residual-resampling indices
# in estimate_power_one_tau() were the same sequence for every analyte and every
# cohort. Monte Carlo errors were perfectly correlated across the array rather
# than independent, so averaging across tasks did not reduce them.
#
# Offset by task_id so each task gets its own stream while the whole run stays
# reproducible from `seed`. Note the manifest is built BEFORE this point and is
# deterministic, so it is unaffected.
task_seed <- if (!is.na(task_id)) seed + as.integer(task_id) else seed
set.seed(task_seed)
message(sprintf("[seed] base=%d task_id=%s -> task_seed=%d",
                seed, ifelse(is.na(task_id), "NA", task_id), task_seed))

load(targeted_rdata)
load(sdoh_rdata)
stopifnot(exists("Targeted_Data_All"), exists("SDOH_Data_All"))

child_covariates <- c("Child_Age", "Child_Sex", "Mat_Child_Comb_Race", "Mat_Edu", "Sample_Collection_Year")
maternal_covariates <- c("Mat_Age", "Mat_Race", "Mat_Edu", "Sample_Collection_Year")

child_sdoh <- if (run_child) build_sdoh_panel(SDOH_Data_All, Targeted_Data_All, "Child_PID", "child") else NULL
maternal_sdoh <- if (run_maternal) build_sdoh_panel(SDOH_Data_All, Targeted_Data_All, "Mat_PID", "maternal") else NULL
child_targeted <- if (run_child) build_targeted_panel(Targeted_Data_All, "Child_PID", use_log1p) else NULL
maternal_targeted <- if (run_maternal) build_targeted_panel(Targeted_Data_All, "Mat_PID", use_log1p) else NULL

make_task_manifest <- function(panel_name, panel_targeted, id_var) {
  study_order <- panel_targeted %>%
    distinct(Study, .data[[id_var]]) %>%
    count(Study, sort = TRUE)
  keep_studies <- if (is.na(max_studies)) study_order$Study else head(study_order$Study, max_studies)

  analyte_order <- panel_targeted %>%
    filter(Study %in% keep_studies) %>%
    count(Analyte_Code, sort = TRUE)
  keep_analytes <- if (is.na(max_analytes)) analyte_order$Analyte_Code else head(analyte_order$Analyte_Code, max_analytes)

  within_tasks <- panel_targeted %>%
    filter(Study %in% keep_studies, Analyte_Code %in% keep_analytes) %>%
    distinct(Study, Analyte_Code) %>%
    arrange(Study, Analyte_Code) %>%
    mutate(
      panel = panel_name,
      analysis = "within_study_demographics",
      study_set = map(Study, ~ .x)
    ) %>%
    rename(study = Study, analyte = Analyte_Code)

  pooled_tasks <- if (run_pooled_study_analysis) {
    tibble(
      study = "Pooled",
      analyte = keep_analytes,
      panel = panel_name,
      analysis = "pooled_study",
      study_set = rep(list(keep_studies), length(keep_analytes))
    )
  } else {
    tibble()
  }

  bind_rows(within_tasks, pooled_tasks)
}

manifest <- bind_rows(
  if (run_child) make_task_manifest("child", child_targeted, "Child_PID") else tibble(),
  if (run_maternal) make_task_manifest("maternal", maternal_targeted, "Mat_PID") else tibble()
) %>%
  mutate(task_id = dplyr::row_number()) %>%
  select(task_id, panel, analysis, study, analyte, study_set)

manifest_write <- manifest %>%
  mutate(study_set = vapply(study_set, paste, character(1), collapse = ";")) %>%
  select(task_id, panel, analysis, study, analyte, study_set)

manifest_file <- file.path(out_dir, "task_manifest.csv")
readr::write_csv(manifest_write, manifest_file)

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
    child_sdoh = child_sdoh,
    maternal_sdoh = maternal_sdoh,
    child_targeted = child_targeted,
    maternal_targeted = maternal_targeted,
    child_covariates = child_covariates,
    maternal_covariates = maternal_covariates,
    taus = taus,
    nsim = nsim,
    alpha = alpha,
    min_n = min_n
  )

  saveRDS(
    list(
      task = task_row %>% mutate(study_set = vapply(study_set, paste, character(1), collapse = ";")),
      screening = result$screening,
      power = result$power,
      config = list(
        taus = taus,
        nsim = nsim,
        alpha = alpha,
        min_n = min_n,
        seed = seed,
        task_seed = task_seed,   # recorded 2026-08-11 so a task is reproducible alone
        use_log1p = use_log1p
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

  readr::write_csv(screening_results, file.path(out_dir, "targeted_quantile_power_screening.csv"))
  readr::write_csv(power_results, file.path(out_dir, "targeted_quantile_power_results.csv"))
  saveRDS(
    list(
      manifest = manifest_write,
      screening_results = screening_results,
      power_results = power_results
    ),
    file.path(out_dir, "targeted_quantile_power_results.rds")
  )

  cat("Combined", length(task_files), "task files.\n")
  quit(save = "no", status = 0)
}

stop("Unknown mode: ", mode)
