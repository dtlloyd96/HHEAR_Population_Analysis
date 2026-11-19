#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
})

root <- normalizePath('.', mustWork = TRUE)
old_path <- file.path(root, 'Input', 'Harmonized_Datasets', 'Harmonized_Outcomes.csv')
new_path <- file.path(root, 'Input', 'Health_Outcomes_From_DataDicts', 'Health_Outcomes_Stacked_Cleaned.csv')
out_dir <- file.path(root, 'Health_Outcome_Analysis', 'Input')
out_path <- file.path(out_dir, 'Harmonized_Outcomes.csv')
out_rdata <- file.path(out_dir, 'Harmonized_Outcomes.RData')
report_path <- file.path(out_dir, 'Harmonized_Outcomes_Compatibility_Report.csv')

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

old <- read_csv(old_path, show_col_types = FALSE)
new <- read_csv(new_path, show_col_types = FALSE)

# canonical name for loose matching
canon <- function(x) tolower(gsub('[^a-z0-9]', '', x))
new_keys <- canon(names(new))

# Ensure Study values use *_Outcomes convention
if ('Study' %in% names(new)) {
  new$Study <- gsub('^Study_', '', new$Study)
  new$Study <- ifelse(grepl('_Outcomes$', new$Study), new$Study, paste0(new$Study, '_Outcomes'))
}

# Build compatibility frame with old columns first
compat <- new[0, , drop = FALSE]
compat <- data.frame(matrix(nrow = nrow(new), ncol = 0))

map_rows <- list()
for (oc in names(old)) {
  if (oc %in% names(new)) {
    compat[[oc]] <- new[[oc]]
    map_rows[[length(map_rows)+1]] <- data.frame(old_col = oc, source_col = oc, status = 'exact', stringsAsFactors = FALSE)
  } else {
    key <- canon(oc)
    idx <- match(key, new_keys)
    if (!is.na(idx)) {
      src <- names(new)[idx]
      compat[[oc]] <- new[[src]]
      map_rows[[length(map_rows)+1]] <- data.frame(old_col = oc, source_col = src, status = 'canonical_match', stringsAsFactors = FALSE)
    } else {
      compat[[oc]] <- NA
      map_rows[[length(map_rows)+1]] <- data.frame(old_col = oc, source_col = NA_character_, status = 'missing_filled_na', stringsAsFactors = FALSE)
    }
  }
}

# Add new-only columns at end (preserve new content for future use)
extra_cols <- setdiff(names(new), names(compat))
if (length(extra_cols)) {
  compat <- bind_cols(compat, new[, extra_cols, drop = FALSE])
}

# Keep Study first to match downstream expectations
if ('Study' %in% names(compat)) {
  ordered <- c('Study', setdiff(names(compat), 'Study'))
  compat <- compat[, ordered, drop = FALSE]
}

# Coerce PID cols to character
for (idc in intersect(c('Child_PID','Mat_PID','ID','Study'), names(compat))) {
  compat[[idc]] <- as.character(compat[[idc]])
}

Harmonized_Outcomes <- compat
write_csv(Harmonized_Outcomes, out_path)
save(Harmonized_Outcomes, file = out_rdata)
write_csv(bind_rows(map_rows), report_path)

cat('Wrote:', out_path, '\n')
cat('Wrote:', out_rdata, '\n')
cat('Wrote:', report_path, '\n')
cat('Rows:', nrow(Harmonized_Outcomes), ' Cols:', ncol(Harmonized_Outcomes), '\n')
