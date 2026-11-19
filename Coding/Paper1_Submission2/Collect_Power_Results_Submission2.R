#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
  library(purrr)
  library(tibble)
})

parse_args <- function(args) {
  out <- list()
  for (arg in args) {
    if (!grepl("=", arg, fixed = TRUE)) {
      next
    }
    parts <- strsplit(arg, "=", fixed = TRUE)[[1]]
    out[[parts[[1]]]] <- paste(parts[-1], collapse = "=")
  }
  out
}

`%||%` <- function(x, y) {
  if (is.null(x) || length(x) == 0 || is.na(x) || x == "") y else x
}

args <- parse_args(commandArgs(trailingOnly = TRUE))

default_out_dir <- "Output/Paper1_Output_Submission2/Power_Output"
out_dir <- normalizePath(args$out_dir %||% default_out_dir, winslash = "/", mustWork = FALSE)
task_dir <- file.path(out_dir, "task_results")
manifest_file <- file.path(out_dir, "task_manifest.csv")

if (!dir.exists(task_dir)) {
  stop("Task results directory not found: ", task_dir)
}

task_files <- list.files(task_dir, pattern = "^task_[0-9]{6}\\.rds$", full.names = TRUE)
if (length(task_files) == 0) {
  stop("No task result files found in: ", task_dir)
}

message("Reading ", length(task_files), " task result files from ", task_dir)

task_results <- lapply(task_files, readRDS)

screening_results <- bind_rows(lapply(task_results, `[[`, "screening"))
power_results <- bind_rows(lapply(task_results, `[[`, "power"))
task_table <- bind_rows(lapply(task_results, `[[`, "task"))

manifest_tbl <- if (file.exists(manifest_file)) {
  readr::read_csv(manifest_file, show_col_types = FALSE)
} else {
  tibble()
}

completed_task_ids <- sort(unique(task_table$task_id))
missing_task_tbl <- if (nrow(manifest_tbl) > 0) {
  manifest_tbl %>%
    filter(!task_id %in% completed_task_ids)
} else {
  tibble()
}

summary_tbl <- tibble(
  manifest_tasks = if (nrow(manifest_tbl) > 0) nrow(manifest_tbl) else NA_integer_,
  completed_tasks = length(completed_task_ids),
  missing_tasks = nrow(missing_task_tbl),
  screening_rows = nrow(screening_results),
  power_rows = nrow(power_results)
)

screening_file <- file.path(out_dir, "targeted_quantile_power_screening.csv")
power_file <- file.path(out_dir, "targeted_quantile_power_results.csv")
task_file <- file.path(out_dir, "targeted_quantile_power_completed_tasks.csv")
missing_file <- file.path(out_dir, "targeted_quantile_power_missing_tasks.csv")
summary_file <- file.path(out_dir, "targeted_quantile_power_collection_summary.csv")
rds_file <- file.path(out_dir, "targeted_quantile_power_collected_results.rds")

readr::write_csv(screening_results, screening_file)
readr::write_csv(power_results, power_file)
readr::write_csv(task_table, task_file)
readr::write_csv(summary_tbl, summary_file)

if (nrow(missing_task_tbl) > 0) {
  readr::write_csv(missing_task_tbl, missing_file)
}

saveRDS(
  list(
    summary = summary_tbl,
    manifest = manifest_tbl,
    completed_tasks = task_table,
    missing_tasks = missing_task_tbl,
    screening_results = screening_results,
    power_results = power_results
  ),
  rds_file
)

message("Wrote:")
message("  ", screening_file)
message("  ", power_file)
message("  ", task_file)
message("  ", summary_file)
if (nrow(missing_task_tbl) > 0) {
  message("  ", missing_file)
}
message("  ", rds_file)

print(summary_tbl)
