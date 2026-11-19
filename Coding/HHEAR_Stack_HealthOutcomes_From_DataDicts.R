#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(readxl)
  library(dplyr)
  library(stringr)
  library(readr)
})

args <- commandArgs(trailingOnly = TRUE)

find_project_root <- function(start_dir) {
  cur <- normalizePath(start_dir, mustWork = TRUE)
  repeat {
    has_datasets <- dir.exists(file.path(cur, 'Datasets'))
    has_input <- dir.exists(file.path(cur, 'Input')) || dir.exists(file.path(cur, 'input'))
    if (has_datasets && has_input) return(cur)
    parent <- dirname(cur)
    if (identical(parent, cur)) break
    cur <- parent
  }
  normalizePath(start_dir, mustWork = TRUE)
}

clean_hdr <- function(x) {
  x <- tolower(trimws(as.character(x)))
  x <- gsub('[^a-z0-9]+', '_', x)
  gsub('(^_+|_+$)', '', x)
}

to_readable_name <- function(x) {
  y <- trimws(as.character(x))
  y <- gsub('[^A-Za-z0-9]+', '_', y)
  y <- gsub('(^_+|_+$)', '', y)
  y <- gsub('_+', '_', y)
  ifelse(y == '' | is.na(y), NA_character_, y)
}

standardize_pid_columns <- function(df) {
  nms <- names(df)
  nms_low <- tolower(nms)

  child_idx <- which(nms_low %in% c('child_pid', 'childpid', 'child_id', 'childid', 'child_chear_pid'))
  mat_idx <- which(nms_low %in% c('mat_pid', 'mother_pid', 'motherpid', 'maternal_pid', 'mom_pid', 'pid_mother'))
  pid_idx <- which(nms_low %in% c('pid', 'chear_pid', 'subject_id', 'participant_id'))

  if (length(child_idx)) names(df)[child_idx[1]] <- 'Child_PID'
  if (length(mat_idx)) names(df)[mat_idx[1]] <- 'Mat_PID'
  if (!'Child_PID' %in% names(df) && !'Mat_PID' %in% names(df) && length(pid_idx)) names(df)[pid_idx[1]] <- 'Child_PID'

  if (!'Child_PID' %in% names(df)) df$Child_PID <- NA_character_
  if (!'Mat_PID' %in% names(df)) df$Mat_PID <- NA_character_
  df$ID <- dplyr::coalesce(as.character(df$Child_PID), as.character(df$Mat_PID))
  df
}

header_row_index <- function(df, required_patterns) {
  idx <- which(apply(df, 1, function(r) {
    row_vals <- tolower(as.character(r))
    all(vapply(required_patterns, function(p) any(grepl(p, row_vals, perl = TRUE), na.rm = TRUE), logical(1)))
  }))
  if (length(idx) == 0) NA_integer_ else idx[1]
}

find_dictionary_files <- function(datasets_dir) {
  files <- list.files(datasets_dir, pattern = '\\.xls[x]?$', full.names = TRUE, recursive = TRUE)
  files[grepl('ddcb|data[_ -]?dictionary', basename(files), ignore.case = TRUE)]
}

extract_study_num_from_rdata <- function(path) {
  m <- regexec('^[^_]+_([0-9]+)_Data_All\\.RData$', basename(path), perl = TRUE)
  r <- regmatches(basename(path), m)[[1]]
  if (length(r) >= 2) r[2] else NA_character_
}

extract_component <- function(dict_path, study_num = NA_character_) {
  nm <- toupper(tools::file_path_sans_ext(basename(dict_path)))
  nm <- gsub('DATA[_ -]?DICTIONARY', '', nm, perl = TRUE)
  nm <- gsub('DDCB', '', nm, perl = TRUE)
  nm <- gsub('[- ]+', '_', nm, perl = TRUE)

  if (!is.na(study_num) && nzchar(study_num)) {
    nm <- gsub(paste0('^.*', study_num, '_?'), '', nm, perl = TRUE)
  } else {
    nm <- gsub('^[0-9_]+', '', nm, perl = TRUE)
  }

  nm <- gsub('^_+|_+$', '', nm, perl = TRUE)
  nm <- gsub('_+', '_', nm, perl = TRUE)

  if (grepl('EPI_[0-9]+', nm, perl = TRUE)) return(str_extract(nm, 'EPI_[0-9]+'))
  if (grepl('(^|_)EPI($|_)', nm, perl = TRUE)) return('EPI')
  nm
}

match_object_name <- function(data_list_names, study, component, dlist = NULL, dict_vars = character(0)) {
  upper_names <- toupper(data_list_names)
  study_prefix <- paste0(toupper(study), '_')
  candidates <- data_list_names[startsWith(upper_names, study_prefix)]
  if (!length(candidates)) candidates <- data_list_names

  comp <- toupper(component)
  if (comp == '' || comp == 'EPI') {
    epi <- candidates[grepl('_EPI($|_[0-9]+$)', toupper(candidates))]
    if (length(epi)) return(epi[1])
    return(candidates[1])
  }

  target <- candidates[grepl(paste0('_', comp, '($|_)'), toupper(candidates))]
  if (length(target)) return(target[1])

  if (grepl('EPI', comp)) {
    epi <- candidates[grepl('_EPI($|_[0-9]+$)', toupper(candidates))]
    if (length(epi)) return(epi[1])
  }

  if (!is.null(dlist) && length(dict_vars)) {
    df_candidates <- candidates[vapply(candidates, function(nm) is.data.frame(dlist[[nm]]), logical(1))]
    if (length(df_candidates)) {
      dv <- tolower(dict_vars)
      score <- vapply(df_candidates, function(nm) sum(tolower(names(dlist[[nm]])) %in% dv), numeric(1))
      if (max(score, na.rm = TRUE) > 0) return(df_candidates[which.max(score)])
    }
  }

  if (length(candidates)) candidates[1] else NA_character_
}

# Heuristic outcome detector from dictionary metadata
is_health_outcome <- function(varname, vardesc, question) {
  txt <- tolower(paste(varname, vardesc, question, sep = ' '))

  include <- c(
    'bmi', 'body mass', 'asthma', 'wheez', 'eczema', 'allerg', 'diabet', 'glucose', 'insulin',
    'hba1c', 'cholesterol', 'triglycer', 'blood pressure', 'systolic', 'diastolic', 'birthweight',
    'low birthweight', 'preterm', 'gest age', 'gestational age', 'apgar', 'autis', 'ados',
    'depress', 'anxiety', 'cognitive', 'language composite', 'motor composite', 'socialemotional',
    'motor', 'fine motor', 'visual', 'vigilance', 'processing', 'tscore', 't score', 'score',
    'cancer', 'disease', 'syndrome', 'fibrosis', 'steatosis', 'nash', 'nafld', 'homa', 'homa_ir',
    'crp', 'hscrp', 'fvc', 'fev1', 'bronchodilator', 'lung function', 'severity', 'symptom', 'outcome',
    'sbp', 'dbp', 'ldl', 'hdl', 'trig', 'dm type', 'diagnos'
  )

  exclude <- c(
    'pid', 'chear pid', 'subject id', 'matrix', 'analyte',
    'concentration', 'unit', 'lod', 'lq', 'assay', 'method', 'lab', 'specimen', 'barcode',
    'race', 'ethnicity', 'income', 'education', 'zip code', 'address', 'site location'
  )

  has_inc <- any(vapply(include, function(k) grepl(k, txt, fixed = TRUE), logical(1)))
  has_exc <- any(vapply(exclude, function(k) grepl(k, txt, fixed = TRUE), logical(1)))

  has_inc && !has_exc
}

read_dict_outcome_vars <- function(dict_path) {
  sheets <- excel_sheets(dict_path)
  dd_sheet <- sheets[grepl('data dictionary|dictionary', sheets, ignore.case = TRUE)][1]
  if (is.na(dd_sheet)) dd_sheet <- sheets[1]

  raw <- suppressMessages(read_excel(dict_path, sheet = dd_sheet, col_names = FALSE))
  hdr_i <- header_row_index(raw, c('varname'))
  if (is.na(hdr_i)) return(data.frame(varname = character(), stringsAsFactors = FALSE))

  hdr <- clean_hdr(unlist(raw[hdr_i, ], use.names = FALSE))
  hdr[is.na(hdr) | hdr == ''] <- paste0('x', seq_along(hdr))[is.na(hdr) | hdr == '']
  names(raw) <- make.unique(hdr)

  dd <- raw[(hdr_i + 1):nrow(raw), , drop = FALSE] %>% mutate(across(everything(), as.character))

  vn_col <- names(dd)[grepl('varname', names(dd))][1]
  vd_col <- names(dd)[grepl('vardesc|description', names(dd))][1]
  q_col <- names(dd)[grepl('question', names(dd))][1]
  # Fallback for odd header parsing in some dictionaries
  if (is.na(vn_col) && ncol(dd) >= 1) vn_col <- names(dd)[1]
  if (is.na(vd_col) && ncol(dd) >= 4) vd_col <- names(dd)[4]
  if (is.na(q_col) && ncol(dd) >= 5) q_col <- names(dd)[5]
  if (is.na(vn_col)) return(data.frame(varname = character(), readable = character(), stringsAsFactors = FALSE))

  dd$varname <- trimws(dd[[vn_col]])
  dd$vardesc <- if (!is.na(vd_col)) dd[[vd_col]] else ''
  dd$question <- if (!is.na(q_col)) dd[[q_col]] else ''

  dd <- dd %>% filter(!is.na(varname), varname != '')
  dd <- dd %>%
    rowwise() %>%
    mutate(is_outcome = is_health_outcome(varname, vardesc, question)) %>%
    ungroup() %>%
    filter(is_outcome) %>%
    mutate(readable = to_readable_name(vardesc)) %>%
    distinct(varname, readable)

  out <- dd %>% select(varname, readable)
  out$readable <- as.character(out$readable)

  # Fallback for sparse/poorly described dictionaries (e.g., some studies where vardesc is thin):
  # keep most non-admin variables from the Data Dictionary so outcomes are not dropped.
  if (nrow(out) < 5) {
    dd_all <- raw[(hdr_i + 1):nrow(raw), , drop = FALSE] %>% mutate(across(everything(), as.character))
    vn2 <- names(dd_all)[grepl('varname', names(dd_all))][1]
    vd2 <- names(dd_all)[grepl('vardesc|description', names(dd_all))][1]
    q2 <- names(dd_all)[grepl('question', names(dd_all))][1]
    if (is.na(vn2) && ncol(dd_all) >= 1) vn2 <- names(dd_all)[1]
    if (is.na(vd2) && ncol(dd_all) >= 4) vd2 <- names(dd_all)[4]
    if (is.na(q2) && ncol(dd_all) >= 5) q2 <- names(dd_all)[5]

    if (!is.na(vn2)) {
      base <- dd_all %>%
        transmute(
          varname = trimws(.data[[vn2]]),
          vardesc = if (!is.na(vd2)) .data[[vd2]] else '',
          question = if (!is.na(q2)) .data[[q2]] else ''
        ) %>%
        filter(!is.na(varname), varname != '')

      is_admin <- function(varname, vardesc, question) {
        txt <- tolower(paste(varname, vardesc, question, sep = ' '))
        ex <- c(
          'pid', 'chear pid', 'subject id', 'sample', 'collection', 'site', 'location',
          'race', 'ethnicity', 'income', 'education', 'address', 'zip', 'season',
          'visit date', 'timepoint', 'specimen', 'barcode', 'assay', 'method'
        )
        any(vapply(ex, function(k) grepl(k, txt, fixed = TRUE), logical(1)))
      }

      base <- base %>%
        rowwise() %>%
        mutate(admin = is_admin(varname, vardesc, question)) %>%
        ungroup() %>%
        filter(!admin) %>%
        mutate(readable = to_readable_name(vardesc)) %>%
        select(varname, readable) %>%
        distinct()
      base$readable <- as.character(base$readable)

      out <- bind_rows(out, base) %>% distinct()
    }
  }

  out
}

root_arg <- if (length(args) >= 1) args[[1]] else '.'
project_dir <- find_project_root(root_arg)
input_dir <- if (dir.exists(file.path(project_dir, 'Input'))) file.path(project_dir, 'Input') else file.path(project_dir, 'input')
datasets_dir <- file.path(project_dir, 'Datasets')
out_dir <- if (length(args) >= 2) args[[2]] else file.path(input_dir, 'Health_Outcomes_From_DataDicts')

dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

rdata_files <- list.files(input_dir, pattern = '_Data_All\\.RData$', full.names = TRUE)
if (!length(rdata_files)) {
  cleaned_input_dir <- file.path(input_dir, 'Cleaned_DDCB')
  if (dir.exists(cleaned_input_dir)) {
    rdata_files <- list.files(cleaned_input_dir, pattern = '_Data_All_Cleaned\\.RData$', full.names = TRUE)
  }
}
dict_files <- find_dictionary_files(datasets_dir)

if (!length(rdata_files)) stop('No *_Data_All.RData files found in: ', input_dir)
if (!length(dict_files)) stop('No dictionary files found under: ', datasets_dir)

stack_list <- list()
sel_log <- list()

for (rf in sort(rdata_files)) {
  study_raw <- basename(rf)
  study_raw <- sub('_Data_All_Cleaned\\.RData$', '', study_raw)
  study_raw <- sub('_Data_All\\.RData$', '', study_raw)
  study <- sub('_[0-9]+$', '', study_raw)
  study_num <- extract_study_num_from_rdata(rf)
  if (is.na(study_num)) {
    study_num <- str_match(basename(rf), '_([0-9]+)_Data_All')[,2]
  }

  env <- new.env(parent = emptyenv())
  load(rf, envir = env)

  dl_name <- ls(env, pattern = '_data_list$')[1]
  if (is.na(dl_name)) next
  dlist <- env[[dl_name]]
  if (!is.list(dlist)) next

  study_dicts <- dict_files[grepl(paste0('/', study_num, '/'), dict_files)]
  if (!length(study_dicts)) next

  for (dfi in sort(study_dicts)) {
    dmeta <- read_dict_outcome_vars(dfi)
    dvars <- dmeta$varname
    dread <- dmeta$readable
    dkeep <- unique(c(dvars, dread[!is.na(dread) & dread != '']))
    if (!length(dkeep)) {
      sel_log[[length(sel_log)+1]] <- data.frame(
        study = study, study_num = study_num, dict_file = basename(dfi), object_name = NA_character_,
        selected_vars = 0L, rows = 0L, status = 'no_outcomes_detected', stringsAsFactors = FALSE
      )
      next
    }

    comp <- extract_component(dfi, study_num)
    obj <- match_object_name(names(dlist), study, comp, dlist = dlist, dict_vars = dvars)

    if (is.na(obj) || !obj %in% names(dlist) || !is.data.frame(dlist[[obj]])) {
      sel_log[[length(sel_log)+1]] <- data.frame(
        study = study, study_num = study_num, dict_file = basename(dfi), object_name = NA_character_,
        selected_vars = 0L, rows = 0L, status = 'no_matching_dataframe_object', stringsAsFactors = FALSE
      )
      next
    }

    dat <- standardize_pid_columns(dlist[[obj]])
    keep <- names(dat)[tolower(names(dat)) %in% tolower(dkeep)]

    # Keep subject IDs if present for linkage
    id_keep <- intersect(names(dat), c('Child_PID','Mat_PID'))
    keep <- unique(c(id_keep, keep))
    outcome_keep <- setdiff(keep, id_keep)

    if (!length(outcome_keep)) {
      sel_log[[length(sel_log)+1]] <- data.frame(
        study = study, study_num = study_num, dict_file = basename(dfi), object_name = obj,
        selected_vars = 0L, rows = nrow(dat), status = 'no_selected_outcome_columns_present', stringsAsFactors = FALSE
      )
      next
    }

    out <- dat[, keep, drop = FALSE] %>% mutate(across(everything(), as.character))
    out <- out[, setdiff(names(out), 'ID'), drop = FALSE]
    out$Study <- paste0(study, '_Outcomes')
    out$Study_Number <- study_num
    out$Source_Object <- obj

    stack_list[[length(stack_list)+1]] <- out

    sel_log[[length(sel_log)+1]] <- data.frame(
      study = study, study_num = study_num, dict_file = basename(dfi), object_name = obj,
      selected_vars = length(keep), rows = nrow(out), status = 'stacked', stringsAsFactors = FALSE
    )
  }
}

selection_log <- bind_rows(sel_log)
write_csv(selection_log, file.path(out_dir, 'health_outcome_selection_log.csv'))

if (!length(stack_list)) stop('No stacked health-outcome data was produced. Check selection log.')

Health_Outcomes_Stacked <- bind_rows(stack_list)
save(Health_Outcomes_Stacked, file = file.path(out_dir, 'Health_Outcomes_Stacked_From_DataDicts.RData'))
write_csv(Health_Outcomes_Stacked, file.path(out_dir, 'Health_Outcomes_Stacked_From_DataDicts.csv'))

message('Saved: ', file.path(out_dir, 'Health_Outcomes_Stacked_From_DataDicts.RData'))
message('Saved: ', file.path(out_dir, 'Health_Outcomes_Stacked_From_DataDicts.csv'))
message('Saved: ', file.path(out_dir, 'health_outcome_selection_log.csv'))
