#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(dplyr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)
root_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = TRUE) else normalizePath('.', mustWork = TRUE)
in_rdata <- if (length(args) >= 2) args[[2]] else file.path(root_dir, 'Input', 'Health_Outcomes_From_DataDicts', 'Health_Outcomes_Stacked_From_DataDicts.RData')
out_dir <- if (length(args) >= 3) args[[3]] else file.path(root_dir, 'Input', 'Health_Outcomes_From_DataDicts')

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

load(in_rdata)
if (!exists('Health_Outcomes_Stacked')) stop('Expected object `Health_Outcomes_Stacked` in: ', in_rdata)

x <- Health_Outcomes_Stacked
names(x) <- make.unique(names(x))

id_cols <- intersect(names(x), c('Study', 'Study_Number', 'Source_Object', 'Child_PID', 'Mat_PID', 'ID'))

na_tokens <- c('', 'na', 'n/a', 'nan', 'null', 'missing', 'unk', 'unknown')
yes_tokens <- c('1','yes','y','true','t','case','present','positive')
no_tokens  <- c('0','no','n','false','f','control','absent','negative')

clean_log <- list()
outlier_log <- list()

safe_chr <- function(v) {
  out <- as.character(v)
  out <- iconv(out, from = "", to = "UTF-8", sub = "")
  out
}

coerce_binary_text <- function(v) {
  y <- tolower(trimws(safe_chr(v)))
  y[y %in% na_tokens] <- NA_character_
  if (!all(is.na(y) | y %in% c(yes_tokens, no_tokens))) return(NULL)
  out <- ifelse(is.na(y), NA_integer_, ifelse(y %in% yes_tokens, 1L, 0L))
  out
}

coerce_binary_two_level <- function(v) {
  y <- trimws(safe_chr(v))
  y[tolower(y) %in% na_tokens] <- NA_character_
  vals <- sort(unique(y[!is.na(y) & y != '']))
  if (length(vals) != 2) return(NULL)
  if (max(nchar(vals)) > 40) return(NULL)
  out <- ifelse(is.na(y), NA_integer_, ifelse(y == vals[2], 1L, 0L))
  list(vec = out, map = paste0(vals[1], '->0; ', vals[2], '->1'))
}

coerce_numeric_like <- function(v, threshold = 0.9) {
  y <- trimws(safe_chr(v))
  y[tolower(y) %in% na_tokens] <- NA_character_
  suppressWarnings(num <- readr::parse_number(y, na = c('', 'NA', 'N/A', 'NaN', 'NULL')))
  nonmiss <- !is.na(y)
  if (sum(nonmiss) == 0) return(NULL)
  parse_rate <- mean(!is.na(num[nonmiss]))
  if (parse_rate < threshold) return(NULL)
  num
}

for (nm in names(x)) {
  if (nm %in% id_cols) next

  v <- x[[nm]]
  v_chr <- as.character(v)
  conv <- 'unchanged'
  detail <- ''

  btxt <- coerce_binary_text(v_chr)
  if (!is.null(btxt)) {
    x[[nm]] <- btxt
    conv <- 'binary_0_1_text'
  } else {
    btwo <- coerce_binary_two_level(v_chr)
    if (!is.null(btwo)) {
      x[[nm]] <- btwo$vec
      conv <- 'binary_0_1_two_level'
      detail <- btwo$map
    } else {
      num <- coerce_numeric_like(v_chr)
      if (!is.null(num)) {
        x[[nm]] <- num
        conv <- 'numeric'
      }
    }
  }

  clean_log[[length(clean_log) + 1]] <- data.frame(
    column = nm,
    conversion = conv,
    detail = detail,
    non_missing = sum(!is.na(x[[nm]])),
    stringsAsFactors = FALSE
  )
}

# Outlier summary using 3*IQR fences on numeric non-binary columns
num_cols <- names(x)[vapply(x, is.numeric, logical(1))]
for (nm in setdiff(num_cols, id_cols)) {
  v <- x[[nm]]
  v_nonmiss <- v[!is.na(v)]
  if (length(v_nonmiss) < 20) next
  if (length(unique(v_nonmiss)) <= 2) next

  q1 <- as.numeric(stats::quantile(v_nonmiss, 0.25, na.rm = TRUE, type = 7))
  q3 <- as.numeric(stats::quantile(v_nonmiss, 0.75, na.rm = TRUE, type = 7))
  iqr <- q3 - q1
  lower <- q1 - 3 * iqr
  upper <- q3 + 3 * iqr

  out_idx <- which(v < lower | v > upper)
  n_out <- length(out_idx)

  if (n_out > 0) {
    flag_col <- paste0(nm, '_Outlier_Flag')
    x[[flag_col]] <- ifelse(is.na(v), NA_integer_, ifelse(v < lower | v > upper, 1L, 0L))
  }

  outlier_log[[length(outlier_log) + 1]] <- data.frame(
    column = nm,
    n = length(v_nonmiss),
    q1 = q1,
    q3 = q3,
    iqr = iqr,
    lower = lower,
    upper = upper,
    n_outliers = n_out,
    pct_outliers = round(100 * n_out / length(v_nonmiss), 3),
    stringsAsFactors = FALSE
  )
}

cleaning_log_df <- bind_rows(clean_log)
outlier_log_df <- bind_rows(outlier_log)

Health_Outcomes_Stacked_Cleaned <- x

save(Health_Outcomes_Stacked_Cleaned, file = file.path(out_dir, 'Health_Outcomes_Stacked_Cleaned.RData'))
write_csv(Health_Outcomes_Stacked_Cleaned, file.path(out_dir, 'Health_Outcomes_Stacked_Cleaned.csv'))
write_csv(cleaning_log_df, file.path(out_dir, 'health_outcomes_column_cleaning_log.csv'))
write_csv(outlier_log_df, file.path(out_dir, 'health_outcomes_outlier_summary.csv'))

message('Saved cleaned outcomes: ', file.path(out_dir, 'Health_Outcomes_Stacked_Cleaned.RData'))
message('Saved cleaned outcomes CSV: ', file.path(out_dir, 'Health_Outcomes_Stacked_Cleaned.csv'))
message('Saved cleaning log: ', file.path(out_dir, 'health_outcomes_column_cleaning_log.csv'))
message('Saved outlier summary: ', file.path(out_dir, 'health_outcomes_outlier_summary.csv'))
