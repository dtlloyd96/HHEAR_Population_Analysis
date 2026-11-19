#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

script_path <- local({
  file_arg <- "--file="
  args <- commandArgs(trailingOnly = FALSE)
  hit <- grep(file_arg, args, value = TRUE)
  if (length(hit) > 0) {
    raw_path <- sub(file_arg, "", hit[[1]])
    raw_path <- gsub("~\\+~", " ", raw_path, perl = TRUE)
    normalizePath(raw_path)
  } else {
    normalizePath(".")
  }
})

project_root <- normalizePath(file.path(dirname(script_path), ".."))
source(file.path(project_root, "Coding", "targeted_analyte_mapping.R"))

args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(args) >= 1) normalizePath(args[[1]], mustWork = FALSE) else {
  file.path(project_root, "HHEAR_Man1_App")
}

dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

message("Writing HHEAR app bundle to: ", output_dir)

std_study <- function(x) {
  stringr::str_remove(x, "[_]\\w+|:")
}

to_ng_per_mL <- function(value, units) {
  dplyr::case_when(
    is.na(units) ~ value,
    units == "ng/mL" ~ value,
    units == "ug/dL" ~ value * 10,
    units %in% c("ng/L") ~ value / 1000,
    units %in% c("pg/mL", "pg/ml") ~ value / 1000,
    units == "mg/L" ~ value * 1000,
    # 2026-11-19: ug/mL (MDUC haptoglobin) was passed through unconverted.
    units == "ug/mL" ~ value * 1000,
    # ng/ml, pmol/mL, pmol/ml, nM (= pmol/mL), ng/g, mg/dL, kg/m3 are left as
    # they are; unit_label() below says what scale they end up on.
    TRUE ~ value
  )
}

# Display unit for a value AFTER to_ng_per_mL(). The harmonized data keeps each
# cohort's lab unit in `Units` and records urine dilution handling in
# `Dilution_Adjustment`; creatinine-corrected values were x100 so that
# (ng/mL)/(mg/dL) reads as ug/g creatinine (HHEAR_Data_Clean_Targeted_V3.Rmd,
# CRE_SCALE). The same arithmetic turns pmol/mL into nmol/g creatinine.
# Added 2026-11-19 so the app can label axes, hovers and downloads.
unit_label <- function(units, matrix_class, dilution) {
  base <- dplyr::case_when(
    is.na(units) ~ "ng/mL",
    units %in% c("ng/mL", "ng/ml", "ng/L", "pg/mL", "pg/ml", "ug/dL", "mg/L", "ug/mL") ~ "ng/mL",
    units %in% c("pmol/mL", "pmol/ml", "nM") ~ "pmol/mL",
    TRUE ~ units
  )
  dilution <- dplyr::coalesce(dilution, "none")
  is_urine <- !is.na(matrix_class) & matrix_class == "Urine"
  dplyr::case_when(
    is_urine & dilution == "creatinine_ug_per_g" & base == "ng/mL" ~ "µg/g creatinine",
    is_urine & dilution == "creatinine_ug_per_g" & base == "pmol/mL" ~ "nmol/g creatinine",
    is_urine & grepl("^SG", dilution) ~ paste0(base, ", SG-adjusted"),
    is_urine ~ paste0(base, ", unadjusted"),
    TRUE ~ base
  )
}

# NHANES extract (NHANES_Subset_Intervals.RData) stores these three urinary PAH
# metabolites in ng/L; every other analyte there is ng/mL. See project memory
# nhanes-extract-pah-units-are-ng-per-l. NHANES urine is NOT creatinine-adjusted
# (the extract carries no creatinine), so it is labelled "unadjusted".
NHANES_NG_PER_L <- c("NAP1", "NAP2", "PYR1")

safe_interval_width <- function(upper, lower) {
  ifelse(!is.na(upper) & !is.na(lower), upper - lower, NA_real_)
}

cap_with_na <- function(x, cap) {
  ifelse(is.na(x), NA_real_, pmin(x, cap))
}

expected_taus <- c("Tau05", "Tau25", "Tau50", "Tau75", "Tau95")
tau_name_map <- c(
  "0.05" = "Tau05",
  "0.25" = "Tau25",
  "0.5" = "Tau50",
  "0.50" = "Tau50",
  "0.75" = "Tau75",
  "0.95" = "Tau95",
  "Tau05" = "Tau05",
  "Tau25" = "Tau25",
  "Tau50" = "Tau50",
  "Tau75" = "Tau75",
  "Tau95" = "Tau95"
)

get_quant_object <- function(x) {
  if (is.null(x) || length(x) == 0 || is.character(x)) {
    return(NULL)
  }
  if ("Quantile_Results" %in% names(x)) {
    return(x)
  }
  if (length(x[[1]]) == 0 || is.character(x[[1]])) {
    return(NULL)
  }
  x[[1]]
}

get_quant_table <- function(x) {
  q_obj <- get_quant_object(x)
  if (is.null(q_obj) || !"Quantile_Results" %in% names(q_obj) || length(q_obj[["Quantile_Results"]]) < 2) {
    return(NULL)
  }
  q_obj[["Quantile_Results"]][[2]]
}

calc_weighted_taus <- function(quant_table) {
  out <- c(Tau05 = NA_real_, Tau25 = NA_real_, Tau50 = NA_real_, Tau75 = NA_real_, Tau95 = NA_real_)
  if (is.null(quant_table)) {
    return(out)
  }

  tau_data <- as.data.frame(quant_table)
  if (ncol(tau_data) == 0) {
    return(out)
  }

  if (is.null(names(tau_data))) {
    if (ncol(tau_data) == 4) {
      names(tau_data) <- c("Tau25", "Tau50", "Tau75", "Tau95")
    } else if (ncol(tau_data) == 5) {
      names(tau_data) <- expected_taus
    } else {
      names(tau_data) <- paste0("V", seq_len(ncol(tau_data)))
    }
  }

  mapped_names <- unname(tau_name_map[names(tau_data)])
  valid_cols <- !is.na(mapped_names)
  tau_data <- tau_data[, valid_cols, drop = FALSE]
  if (ncol(tau_data) == 0) {
    return(out)
  }

  names(tau_data) <- mapped_names[valid_cols]
  tau_data <- tau_data[, !duplicated(names(tau_data)), drop = FALSE]

  for (tau in intersect(names(out), names(tau_data))) {
    curr_values <- suppressWarnings(as.numeric(tau_data[[tau]]))
    curr_values <- curr_values[!is.na(curr_values)]
    if (length(curr_values) > 0) {
      out[tau] <- mean(curr_values)
    }
  }

  out
}

adjust_interval_df <- function(df) {
  df %>%
    mutate(across(where(is.numeric), ~ round(.x, 3))) %>%
    mutate(
      IQR = safe_interval_width(Tau75, Tau25),
      Max_Min = safe_interval_width(Tau95, Tau05),
      Has_IQR = !is.na(Tau25) & !is.na(Tau75),
      Has_Outer = !is.na(Tau05) & !is.na(Tau95),
      Has_Median = !is.na(Tau50)
    ) %>%
    mutate(
      IQR_unadj = IQR,
      IQR = cap_with_na(IQR, 20),
      Max_Min_unadj = Max_Min,
      Max_Min = cap_with_na(Max_Min, 50),
      Flag = as.factor(ifelse(!is.na(IQR_unadj) & IQR_unadj >= 20, 1, 0))
    ) %>%
    mutate(
      Max_Min_Log = log10(Max_Min + 1),
      IQR_Log = log10(IQR + 1)
    ) %>%
    select(-any_of(c("Has_IQR", "Has_Outer", "Has_Median")))
}

build_targeted_lists <- function(targeted_data) {
  child_data <- targeted_data %>%
    filter(!is.na(Child_PID), !is.na(Analyte_Code)) %>%
    mutate(
      Units = if_else(is.na(Units), "ng/mL", Units),
      Concentration = to_ng_per_mL(Concentration, Units)
    ) %>%
    select(-Mat_PID) %>%
    group_by(Study, Child_PID, Analyte_Code) %>%
    summarise(Concentration = mean(Concentration, na.rm = TRUE), .groups = "drop")

  mat_data <- targeted_data %>%
    filter(!is.na(Mat_PID), !is.na(Analyte_Code)) %>%
    mutate(
      Units = if_else(is.na(Units), "ng/mL", Units),
      Concentration = to_ng_per_mL(Concentration, Units)
    ) %>%
    select(-Child_PID) %>%
    group_by(Study, Mat_PID, Analyte_Code) %>%
    summarise(Concentration = mean(Concentration, na.rm = TRUE), .groups = "drop")

  child_list <- purrr::map(
    sort(unique(child_data$Analyte_Code)),
    ~ child_data %>%
      filter(Analyte_Code == .x) %>%
      mutate(Study = std_study(Study)) %>%
      select(Study, Child_PID, Analyte_Code, Concentration) %>%
      tidyr::pivot_wider(names_from = Analyte_Code, values_from = Concentration) %>%
      select(Child_PID, Study, all_of(.x))
  ) %>% rlang::set_names(sort(unique(child_data$Analyte_Code)))

  mat_list <- purrr::map(
    sort(unique(mat_data$Analyte_Code)),
    ~ mat_data %>%
      filter(Analyte_Code == .x) %>%
      mutate(Study = std_study(Study)) %>%
      select(Study, Mat_PID, Analyte_Code, Concentration) %>%
      tidyr::pivot_wider(names_from = Analyte_Code, values_from = Concentration) %>%
      select(Mat_PID, Study, all_of(.x))
  ) %>% rlang::set_names(sort(unique(mat_data$Analyte_Code)))

  list(
    child_data = child_data,
    mat_data = mat_data,
    child_list = child_list,
    mat_list = mat_list
  )
}

build_base_sdoh <- function(sdoh_data) {
  sdoh_data %>%
    mutate(
      Edu_Any = case_when(!is.na(Max_Edu) ~ Max_Edu, is.na(Max_Edu) ~ Mat_Edu, !is.na(Mat_Edu) ~ Mat_Edu),
      Mat_Child_Comb_Race = case_when(!is.na(Child_Race) ~ Child_Race, is.na(Child_Race) ~ Mat_Race, !is.na(Mat_Race) ~ Mat_Race),
      Mat_Age = as.numeric(Mat_Age),
      Study = std_study(Study)
    )
}

build_location_sdoh <- function(sdoh_data, input_dir) {
  load(file.path(input_dir, "TEDDY_SDOH_Clean.RData"))
  load(file.path(input_dir, "ZIP_SDOH_Clean.RData"))

  zip_data <- ZIP_SDOH %>%
    mutate(
      Country = ifelse(
        is.na(Country),
        Country[match(Site_Location, Site_Location[!is.na(Country)])],
        Country
      ),
      Study = "ZIP"
    ) %>%
    select(Mat_PID, Study, Location = Country)

  teddy_data <- TEDDY_SDOH %>%
    mutate(Study = "TEDDY") %>%
    select(Child_PID, Study, Location)

  exwas_sdoh <- build_base_sdoh(sdoh_data) %>%
    filter(Study != "ZIP", Study != "TEDDY") %>%
    select(Child_PID, Mat_PID, Study, Location, Child_Age, Child_Sex, Mat_Age)

  child_sdoh <- bind_rows(
    exwas_sdoh %>% filter(!is.na(Child_PID)) %>% select(-Mat_PID),
    teddy_data
  ) %>% distinct()

  mat_sdoh <- bind_rows(
    exwas_sdoh %>% filter(!is.na(Mat_PID)) %>% select(-Child_PID),
    zip_data
  ) %>% distinct()

  list(child_sdoh = child_sdoh, mat_sdoh = mat_sdoh)
}

build_app_sdoh <- function(sdoh_data) {
  ref_vars <- c("Mat_Age", "Mat_Edu", "Mat_Income", "Max_Edu", "Child_Age", "Mat_Race", "Child_Race", "Child_Sex")
  exwas_sdoh <- build_base_sdoh(sdoh_data)
  list(
    child = exwas_sdoh %>% filter(!is.na(Child_PID)) %>% select(Study, Child_PID, Mat_PID, all_of(ref_vars)) %>% distinct(),
    mat = exwas_sdoh %>% filter(!is.na(Mat_PID)) %>% select(Study, Child_PID, Mat_PID, all_of(ref_vars)) %>% distinct()
  )
}

build_reference_intervals <- function(output_dir, pooled_dir) {
  load(file.path(pooled_dir, "Child_Child_Age_Quantile_Results.RData"))
  load(file.path(pooled_dir, "Child_Child_Sex_Quantile_Results.RData"))
  load(file.path(pooled_dir, "Mat_Mat_Age_Quantile_Results.RData"))

  shared_outcomes <- intersect(names(Child_Child_Age_Quantile_Results), names(Child_Child_Sex_Quantile_Results))

  child_intervals <- purrr::map(shared_outcomes, function(outcome_name) {
    curr_quant <- list(
      Child_Age = get_quant_table(Child_Child_Age_Quantile_Results[[outcome_name]]),
      Child_Sex = get_quant_table(Child_Child_Sex_Quantile_Results[[outcome_name]])
    )

    weighted <- purrr::map(curr_quant, calc_weighted_taus)
    weighted_df <- bind_rows(weighted)
    weighted_df %>%
      summarise(across(all_of(expected_taus), ~ mean(.x, na.rm = TRUE))) %>%
      mutate(across(all_of(expected_taus), ~ ifelse(is.nan(.x), NA_real_, .x)), Outcome = outcome_name) %>%
      select(Outcome, everything())
  }) %>% rlang::set_names(shared_outcomes)

  mat_intervals <- purrr::map(names(Mat_Mat_Age_Quantile_Results), function(outcome_name) {
    weighted <- calc_weighted_taus(get_quant_table(Mat_Mat_Age_Quantile_Results[[outcome_name]]))
    as.data.frame(t(weighted)) %>%
      mutate(Outcome = outcome_name) %>%
      select(Outcome, all_of(expected_taus))
  }) %>% rlang::set_names(names(Mat_Mat_Age_Quantile_Results))

  Final_Intervals_Df_Adjusted <- adjust_interval_df(bind_rows(child_intervals))
  Final_Intervals_Df_Adjusted_Mat <- adjust_interval_df(bind_rows(mat_intervals))

  save(Final_Intervals_Df_Adjusted, file = file.path(output_dir, "Child_Reference_Intervals.RData"))
  save(Final_Intervals_Df_Adjusted_Mat, file = file.path(output_dir, "Mat_Reference_Intervals.RData"))

  list(child = Final_Intervals_Df_Adjusted, mat = Final_Intervals_Df_Adjusted_Mat)
}

# Must stay in step with Location_Mapping in
# Coding/Paper1_Submission2/Vis_Code/Matched_Intervals.Rmd, which writes some of
# the same bundle files. The two drifting apart is how the app ended up labelling
# the same cohort differently in different tabs.
location_mapping <- tibble::tibble(
  Location = c(
    "NorthEast_US", "Washington_State", "Bangladesh", "Chile", "Finland", "Germany", "Sweden",
    "South_Africa", "Uganda", "Tobago", "Ecuador", "Puerto Rico", "Brazil", "Nicaragua", "Guatemala",
    "Accra", "Thailand"
  ),
  # NorthEast_US is SICAS, which the paper (Table S11) calls "Northeastern US";
  # it was labelled "Boston", colliding with OBSS. Changed 2026-11-19.
  New_Location = c(
    "Northeastern US", "Yakima Valley,Washington", "Pabna,Bangladesh", "Santiago,Chile", "Turku,Finland", "Munich,Germany", "Malmö,Sweden",
    "Cape Town,South Africa", "Kampala,Uganda", "Tobago,Trinidad and Tobago", "Pedro Moncayo, Ecuador", "Northern Karst,Puerto Rico",
    "Bahia, Brazil", "Managua,Nicaragua", "Suchitepéquez,Guatemala",
    "Accra,Ghana", "Chiang Mai,Thailand"
  )
)

build_location_bundle <- function(output_dir, interval_refs, targeted_lists, location_sdoh, nhanes_path) {
  load(nhanes_path)
  lookup <- all_plot_data %>%
    group_by(curr_analyte, HHEAR_Name) %>%
    summarise(n = n(), .groups = "drop")
  rename_vec <- setNames(lookup$curr_analyte, lookup$HHEAR_Name)
  set.seed(123)
  all_plot_data_NHANES_subset <- all_plot_data %>%
    mutate(Study = "NHANES") %>%
    rename(any_of(rename_vec)) %>%
    select(-curr_analyte, -ID, -PID, -SDDSRVYR, -curr_desc) %>%
    rename(curr_analyte = HHEAR_Name) %>%
    group_by(curr_analyte) %>%
    group_modify(~ dplyr::slice_sample(.x, n = min(500, nrow(.x)))) %>%
    ungroup() %>%
    # ng/L -> ng/mL for NAP1/NAP2/PYR1 (2026-11-19); without this NHANES sat
    # ~1000x above every HHEAR cohort for those analytes.
    mutate(
      plot_value = ifelse(curr_analyte %in% NHANES_NG_PER_L, plot_value / 1000, plot_value),
      across(any_of(NHANES_NG_PER_L), ~ .x / 1000)
    )
  pah_check <- all_plot_data_NHANES_subset %>%
    filter(curr_analyte %in% NHANES_NG_PER_L) %>%
    group_by(curr_analyte) %>%
    summarise(median = median(plot_value, na.rm = TRUE), .groups = "drop")
  message("NHANES PAH medians after ng/L -> ng/mL: ",
          paste0(pah_check$curr_analyte, "=", signif(pah_check$median, 3), collapse = ", "))

  build_location_side <- function(targeted_list, interval_df, join_df, id_col, population) {
    join_df <- join_df %>% distinct(.data[[id_col]], Study, .keep_all = TRUE)
    eligible_analytes <- purrr::keep(
      purrr::map(interval_df$Outcome, function(curr_analyte) {
        curr_interval <- interval_df %>% filter(Outcome == curr_analyte)
        curr_tau75 <- curr_interval$Tau75[[1]]
        curr_data <- targeted_list[[curr_analyte]]
        if (is.null(curr_data)) {
          return(NULL)
        }
        curr_data <- curr_data %>%
          mutate(
            Above_75 = ifelse(!is.na(curr_tau75) & .data[[curr_analyte]] >= curr_tau75 * 2, 1, 0)
          ) %>%
          rownames_to_column("ID") %>%
          mutate(
            plot_value = .data[[curr_analyte]],
            ID = as.numeric(ID),
            curr_analyte = curr_analyte
          )
        curr_data
      }) %>% rlang::set_names(interval_df$Outcome),
      ~ !is.null(.x)
    )

    overlap_names <- names(eligible_analytes)

    plot_data <- bind_rows(eligible_analytes[overlap_names])
    interval_data <- bind_rows(purrr::map(overlap_names, ~ interval_df %>% filter(Outcome == .x) %>% mutate(curr_analyte = .x)))

    error_df <- plot_data %>%
      group_by(curr_analyte) %>%
      summarise(error_y = mean(ID, na.rm = TRUE), bar_height = n() + 500, .groups = "drop")

    interval_data <- interval_data %>% left_join(error_df, by = "curr_analyte")

    max_ids <- plot_data %>%
      group_by(curr_analyte) %>%
      summarise(max_id = suppressWarnings(max(as.numeric(ID), na.rm = TRUE)), .groups = "drop") %>%
      mutate(max_id = ifelse(is.finite(max_id), max_id, 0))

    nhanes_aug <- all_plot_data_NHANES_subset %>%
      left_join(max_ids, by = "curr_analyte") %>%
      group_by(curr_analyte) %>%
      mutate(ID = max_id + row_number()) %>%
      select(-max_id) %>%
      ungroup()

    combined_plot <- bind_rows(plot_data, nhanes_aug) %>%
      filter(!is.na(ID))

    joined <- left_join(combined_plot, join_df, by = c(id_col, "Study")) %>%
      mutate(Location = ifelse(Study == "NHANES", "NHANES-USA", Location)) %>%
      left_join(location_mapping, by = "Location") %>%
      mutate(Location = coalesce(New_Location, Location))

    # These vectors set the plotting ORDER only. They must never decide
    # membership: a location absent from the list used to become NA and get
    # filtered out, which silently dropped every cohort sited somewhere new --
    # ATLAA, GRAPHS, SWSD, NYUWHS and THS were all missing from the maternal
    # panels for this reason. Anything unlisted is appended instead.
    preferred_levels <- if (population == "Child") {
      c(
        "NHANES-USA", "Baltimore", "Denver", "Sacramento", "SanFrancisco", "Syracuse", "Northeastern US", "New_York",
        "Midwestern_US",
        "United_States", "Yakima Valley,Washington", "Pabna,Bangladesh", "Santiago,Chile", "Pedro Moncayo, Ecuador",
        "Turku,Finland", "Munich,Germany", "Mexico_City", "Cape Town,South Africa", "Malmö,Sweden",
        "Kampala,Uganda", "Tobago,Trinidad and Tobago"
      )
    } else {
      c(
        "NHANES-USA", "Boston", "California", "Los_Angeles", "Massachusetts", "Michigan", "New_Hampshire",
        "Northern Karst,Puerto Rico", "Sacramento", "United_States", "Bahia, Brazil", "Suchitepéquez,Guatemala",
        "Mexico_City", "Managua,Nicaragua", "Cape Town,South Africa", "Suriname"
      )
    }

    observed <- unique(as.character(joined$Location))
    observed <- observed[!is.na(observed)]
    unlisted <- sort(setdiff(observed, preferred_levels))
    if (length(unlisted)) {
      message(population, ": locations not in the preferred order, appended: ",
              paste(unlisted, collapse = ", "))
    }

    n_before <- nrow(joined)
    joined <- joined %>%
      mutate(Location = factor(Location, levels = c(intersect(preferred_levels, observed), unlisted))) %>%
      filter(!is.na(Location))
    if (nrow(joined) != n_before) {
      message(population, ": dropped ", n_before - nrow(joined), " rows with a missing location.")
    }

    interval_joined <- joined %>%
      group_by(curr_analyte) %>%
      summarise(n_unique_locations = n_distinct(Location) / 2, .groups = "drop") %>%
      mutate(level_height = n_unique_locations * 2) %>%
      left_join(interval_data, by = c("curr_analyte" = "curr_analyte"))

    list(plot = joined, interval = interval_joined)
  }

  build_lod_side <- function(id_col, join_df, population) {
    join_df <- join_df %>% distinct(.data[[id_col]], Study, .keep_all = TRUE)

    lod_df <- Targeted_Data_All %>%
      filter(!is.na(.data[[id_col]]), !is.na(Analyte_Code), !is.na(LOD)) %>%
      mutate(
        Study = std_study(Study),
        Units = if_else(is.na(Units), "ng/mL", Units),
        LOD = to_ng_per_mL(LOD, Units),
        curr_analyte = Analyte_Code
      ) %>%
      left_join(join_df, by = c(id_col, "Study")) %>%
      mutate(Location = ifelse(Study == "NHANES", "NHANES-USA", Location)) %>%
      left_join(location_mapping, by = "Location") %>%
      mutate(
        Location = coalesce(New_Location, Location),
        Population = population
      ) %>%
      filter(!is.na(Location)) %>%
      group_by(Population, curr_analyte, Study, Location) %>%
      summarise(LOD_value = median(LOD, na.rm = TRUE), .groups = "drop")

    lod_df
  }

  child_side <- build_location_side(
    targeted_lists$child_list,
    interval_refs$child,
    location_sdoh$child_sdoh,
    "Child_PID",
    "Child"
  )
  mat_side <- build_location_side(
    targeted_lists$mat_list,
    interval_refs$mat,
    location_sdoh$mat_sdoh,
    "Mat_PID",
    "Adult"
  )

  Location_Plotdata_Child <- child_side$plot
  Location_interval_data_Child <- child_side$interval
  Location_Plotdata_Mat <- mat_side$plot
  Location_interval_data_Mat <- mat_side$interval
  Location_LOD_Data_Child <- build_lod_side("Child_PID", location_sdoh$child_sdoh, "Child") %>%
    mutate(Location = factor(Location, levels = levels(Location_Plotdata_Child$Location)))
  Location_LOD_Data_Mat <- build_lod_side("Mat_PID", location_sdoh$mat_sdoh, "Adult") %>%
    mutate(Location = factor(Location, levels = levels(Location_Plotdata_Mat$Location)))

  save(Location_Plotdata_Child, file = file.path(output_dir, "Child_Plot_Data_Shiny.RData"))
  save(Location_interval_data_Child, Location_LOD_Data_Child, file = file.path(output_dir, "Child_Interval_Data_Shiny.RData"))
  save(Location_Plotdata_Mat, file = file.path(output_dir, "MatPlot_Data_Shiny.RData"))
  save(Location_Plotdata_Mat, file = file.path(output_dir, "Mat_Plot_Data_Shiny.RData"))
  save(Location_interval_data_Mat, Location_LOD_Data_Mat, file = file.path(output_dir, "Mat_Interval_Data_Shiny.RData"))
}

build_sdoh_plot_dfs <- function(output_dir, interval_refs, targeted_lists, app_sdoh) {
  mystudies_mat <- unique(app_sdoh$mat$Study)
  myanalytes_mat <- unique(interval_refs$mat$Outcome)
  Plot_Dfs_Mat <- list()

  for (study_name in mystudies_mat) {
    curr_sdoh_data <- app_sdoh$mat %>%
      select(-Child_PID) %>%
      filter(Study == study_name) %>%
      select(where(~ !all(is.na(.)))) %>%
      select(where(~ dplyr::n_distinct(.) > 1)) %>%
      distinct()

    if (!"Mat_PID" %in% names(curr_sdoh_data)) {
      next
    }

    curr_demo_features <- setdiff(names(curr_sdoh_data), "Mat_PID")
    vals_list <- lapply(myanalytes_mat, function(an) {
      df <- targeted_lists$mat_list[[an]] %>% select(Mat_PID, !!sym(an))
      if (nrow(df) == 0) return(NULL)
      names(df)[names(df) == an] <- "Value"
      df$Outcome <- an
      df
    })

    ValuesLong <- bind_rows(vals_list) %>% filter(!is.na(Value))
    if (nrow(ValuesLong) == 0) next

    valid_outcomes <- unique(ValuesLong$Outcome)
    Intervals <- interval_refs$mat %>%
      filter(Outcome %in% valid_outcomes) %>%
      select(Outcome, Tau05, Tau25, Tau50, Tau75, Tau95)

    DF <- curr_sdoh_data %>%
      left_join(ValuesLong, by = "Mat_PID", relationship = "many-to-many") %>%
      left_join(Intervals, by = "Outcome") %>%
      filter(!is.na(Value))

    if (nrow(DF) == 0 || n_distinct(DF$Outcome) == 0) next

    for (group_col in curr_demo_features) {
      DF_plot <- DF
      if (group_col %in% c("Mat_Age", "Child_Age")) {
        DF_plot <- DF_plot %>% mutate(quantile_group = as.factor(dplyr::ntile(.data[[group_col]], 4)))
        group_col <- "quantile_group"
      }
      DF_plot <- DF_plot %>% filter(!is.na(.data[[group_col]]))
      if (nrow(DF_plot) == 0) next
      key <- paste0(study_name, "__", group_col)
      Plot_Dfs_Mat[[key]] <- DF_plot
    }
  }

  mystudies_child <- unique(targeted_lists$child_data %>% mutate(Study = std_study(Study)) %>% pull(Study))
  myanalytes_child <- unique(interval_refs$child$Outcome)
  Plot_Dfs_Child <- list()

  for (study_name in mystudies_child) {
    curr_sdoh_data <- app_sdoh$child %>%
      select(-Mat_PID) %>%
      filter(Study == study_name) %>%
      select(where(~ !all(is.na(.)))) %>%
      select(where(~ dplyr::n_distinct(.) > 1)) %>%
      distinct()

    if (!"Child_PID" %in% names(curr_sdoh_data)) {
      next
    }

    curr_demo_features <- setdiff(names(curr_sdoh_data), "Child_PID")
    vals_list <- lapply(myanalytes_child, function(an) {
      df <- targeted_lists$child_list[[an]] %>% select(Child_PID, !!sym(an))
      if (nrow(df) == 0) return(NULL)
      names(df)[names(df) == an] <- "Value"
      df$Outcome <- an
      df
    })

    ValuesLong <- bind_rows(vals_list) %>% filter(!is.na(Value))
    if (nrow(ValuesLong) == 0) next

    valid_outcomes <- unique(ValuesLong$Outcome)
    Intervals <- interval_refs$child %>%
      filter(Outcome %in% valid_outcomes) %>%
      select(Outcome, Tau05, Tau25, Tau50, Tau75, Tau95)

    DF <- curr_sdoh_data %>%
      left_join(ValuesLong, by = "Child_PID", relationship = "many-to-many") %>%
      left_join(Intervals, by = "Outcome") %>%
      filter(!is.na(Value))

    if (nrow(DF) == 0 || n_distinct(DF$Outcome) == 0) next

    for (group_col in curr_demo_features) {
      DF_plot <- DF
      if (group_col %in% c("Mat_Age", "Child_Age")) {
        DF_plot <- DF_plot %>% mutate(quantile_group = as.factor(dplyr::ntile(.data[[group_col]], 4)))
        group_col <- "quantile_group"
      }
      DF_plot <- DF_plot %>% filter(!is.na(.data[[group_col]]))
      if (nrow(DF_plot) == 0) next
      key <- paste0(study_name, "__", group_col)
      Plot_Dfs_Child[[key]] <- DF_plot
    }
  }

  save(Plot_Dfs_Child, file = file.path(output_dir, "Child_SDOH_Plot_Dfs.RData"))
  save(Plot_Dfs_Mat, file = file.path(output_dir, "Mat_SDOH_Plot_Dfs.RData"))
}

build_match_frame <- function(targeted_list, exposure_df, id_col, age_col, sex_col = NULL, adult = FALSE) {
  purrr::imap(targeted_list, function(df, analyte) {
    pheno_df <- df %>%
      select(all_of(id_col), Study, !!sym(analyte)) %>%
      filter(!is.na(.data[[analyte]]))

    curr <- left_join(pheno_df, exposure_df %>% distinct(.data[[id_col]], Study, .keep_all = TRUE), by = c(id_col, "Study")) %>%
      mutate(value = as.numeric(.data[[analyte]])) %>%
      distinct(.keep_all = TRUE) %>%
      filter(!is.na(Location), !is.na(.data[[age_col]]))

    if (!is.null(sex_col)) {
      curr <- curr %>% filter(!is.na(.data[[sex_col]]))
    }

    if (nrow(curr) == 0 || length(unique(curr$Location)) <= 1) {
      return(NULL)
    }

    band_yrs <- if (adult) {
      (max(curr[[age_col]], na.rm = TRUE) - min(curr[[age_col]], na.rm = TRUE)) / 10
    } else {
      max(curr[[age_col]], na.rm = TRUE) / 4
    }
    if (!is.finite(band_yrs) || band_yrs == 0) band_yrs <- 1

    curr <- curr %>%
      mutate(age_band = cut(.data[[age_col]], breaks = seq(0, 100, by = band_yrs), include.lowest = TRUE))

    if (is.null(sex_col)) {
      curr <- curr %>%
        group_by(age_band) %>%
        filter(n_distinct(Location) > 1) %>%
        mutate(stratum_id = cur_group_id()) %>%
        ungroup()
    } else {
      curr <- curr %>%
        group_by(age_band, .data[[sex_col]]) %>%
        filter(n_distinct(Location) > 1) %>%
        mutate(stratum_id = cur_group_id()) %>%
        ungroup()
    }

    if (nrow(curr) == 0) {
      return(NULL)
    }

    curr %>% mutate(plot_value = value, curr_analyte = analyte)
  }) %>% purrr::keep(~ !is.null(.x)) %>% bind_rows()
}

build_ridge_plot_data <- function(output_dir, targeted_lists, location_sdoh) {
  Child_Match_Plot_Data <- build_match_frame(
    targeted_lists$child_list,
    location_sdoh$child_sdoh,
    "Child_PID",
    "Child_Age",
    "Child_Sex",
    adult = FALSE
  )

  Mat_Match_Plot_Data <- build_match_frame(
    targeted_lists$mat_list,
    location_sdoh$mat_sdoh,
    "Mat_PID",
    "Mat_Age",
    adult = TRUE
  )

  Ridge_Plot_Data <- bind_rows(
    Child_Match_Plot_Data %>%
      mutate(Strata_clean = age_band) %>%
      mutate(age_band_num = str_remove_all(as.character(age_band), "[\\(\\]\\[]")) %>%
      separate(age_band_num, into = c("low", "high"), sep = ",", convert = TRUE) %>%
      mutate(low = round(low), high = round(high), Strata_rounded = paste0("(", low, ",", high, ")")) %>%
      mutate(Strata = paste0(Strata_rounded, ":", Child_Sex), Population = "Child", facet_label = paste0(curr_analyte, "-Child")),
    Mat_Match_Plot_Data %>%
      mutate(Strata_clean = age_band) %>%
      mutate(age_band_num = str_remove_all(as.character(age_band), "[\\(\\]\\[]")) %>%
      separate(age_band_num, into = c("low", "high"), sep = ",", convert = TRUE) %>%
      mutate(low = round(low), high = round(high), Strata_rounded = paste0("(", low, ",", high, ")")) %>%
      mutate(Strata = Strata_rounded, Population = "Adult", facet_label = paste0(curr_analyte, "-Adult"))
  ) %>%
    left_join(location_mapping, by = "Location") %>%
    mutate(Location = coalesce(New_Location, Location))

  save(Ridge_Plot_Data, file = file.path(output_dir, "Ridge_Plot_Data_Shiny.RData"))
}

build_mat_child_data <- function(output_dir, targeted_lists, input_dir, mapping_path) {
  load(file.path(input_dir, "DCHS_SDOH_Clean.RData"))
  load(file.path(input_dir, "ELEMENT_SDOH_Clean.RData"))
  load(file.path(input_dir, "Harmonized_Targeted.RData"))

  targeted_mapping <- read_targeted_analyte_mapping(mapping_path)
  plot_categories <- c("TraceElements", "Phenols", "Phthalates", "PAH", "PCB", "Parabens", "Smoking", "Pesticides", "VOCs", "FlameRetardants", "PFAS")
  targeted_category_lookup <- targeted_mapping %>%
    filter(Category %in% plot_categories) %>%
    select(Code, Category) %>%
    distinct()

  add_plot_category <- function(data, code_col = "Analyte_Code") {
    left_join(data, targeted_category_lookup, by = stats::setNames("Code", code_col))
  }

  studies <- "DCHS|ELEMENT"
  Mapping <- bind_rows(
    DCHS_SDOH %>% select(Mat_PID, Child_PID) %>% mutate(Study = "DCHS_Targeted"),
    ELEMENT_SDOH %>% select(Mat_PID, Child_PID) %>% mutate(Study = "ELEMENT_Targeted")
  )

  Child_Data <- targeted_lists$child_data %>%
    filter(str_detect(Study, studies)) %>%
    rename(Concentration_Child = Concentration)

  Mat_Data <- targeted_lists$mat_data %>%
    filter(str_detect(Study, studies)) %>%
    rename(Concentration_Mat = Concentration)

  Mat_Sample_Yr <- Targeted_Data_All %>%
    filter(!is.na(Mat_PID), str_detect(Study, studies)) %>%
    group_by(Mat_PID, Sample_Collection_Year) %>%
    summarise(n = n(), .groups = "drop")

  All_Data <- Mat_Data %>%
    left_join(Mapping, by = c("Mat_PID", "Study")) %>%
    left_join(Child_Data, by = c("Child_PID", "Study", "Analyte_Code")) %>%
    mutate(Conc_Diff = abs(Concentration_Mat - Concentration_Child)) %>%
    left_join(Mat_Sample_Yr, by = "Mat_PID")

  params <- All_Data %>%
    summarise(
      global_mean = mean(c(Concentration_Mat, Concentration_Child), na.rm = TRUE),
      global_sd = sd(c(Concentration_Mat, Concentration_Child), na.rm = TRUE)
    )

  plot_data <- All_Data %>%
    bind_cols(params) %>%
    mutate(log_ratio = log1p(Concentration_Child) - log1p(Concentration_Mat)) %>%
    mutate(
      highlight = case_when(
        log_ratio >= quantile(log_ratio, 0.90, na.rm = TRUE) ~ "High",
        log_ratio <= quantile(log_ratio, 0.10, na.rm = TRUE) ~ "Low",
        TRUE ~ "Normal"
      )
    ) %>%
    add_plot_category(code_col = "Analyte_Code") %>%
    filter(!is.na(Category))

  mean_tbl <- plot_data %>%
    group_by(Analyte_Code) %>%
    summarise(mean_log_ratio = mean(log_ratio, na.rm = TRUE), .groups = "drop")

  spread_tbl <- plot_data %>%
    group_by(Analyte_Code) %>%
    summarise(spread_IQR = IQR(log_ratio, na.rm = TRUE), .groups = "drop")

  plot_data_vi2 <- plot_data %>%
    left_join(mean_tbl, by = "Analyte_Code") %>%
    mutate(Analyte_Code = reorder(Analyte_Code, mean_log_ratio)) %>%
    left_join(spread_tbl, by = "Analyte_Code") %>%
    mutate(Analyte_Code = fct_reorder(Analyte_Code, spread_IQR, .desc = FALSE))

  # The app shows Study in the dyad legend; strip the "_Targeted" suffix (2026-11-19).
  plot_data_vi2 <- plot_data_vi2 %>% mutate(Study = std_study(Study))
  save(plot_data_vi2, file = file.path(output_dir, "Mat_Child_Data.RData"))
}

write_support_files <- function(output_dir, project_root, mapping_path) {
  file.copy(file.path(project_root, "Input", "Study_Lat_Longs.csv"), file.path(output_dir, "Study_Lat_Longs.csv"), overwrite = TRUE)
  file.copy(file.path(project_root, "Input", "Harmonized_Datasets", "Harmonized_SDOH.RData"), file.path(output_dir, "Harmonized_SDOH.RData"), overwrite = TRUE)
  file.copy(file.path(project_root, "Input", "Harmonized_Datasets", "Harmonized_Targeted.RData"), file.path(output_dir, "Harmonized_Targeted.RData"), overwrite = TRUE)

  chem_name_mapping <- readr::read_csv(mapping_path, show_col_types = FALSE) %>%
    select(Code, Full_Name, CAS) %>%
    distinct()
  readr::write_csv(chem_name_mapping, file.path(output_dir, "Chem_Name_Mapping.csv"))
}

input_dir <- file.path(project_root, "Input", "Harmonized_Datasets")
pooled_dir <- file.path(project_root, "Output", "Paper1_Output_Submission2", "Pooled_Output")
mapping_path <- file.path(project_root, "Input", "Harmonized_Datasets", "Targeted_Analyte_Category_Mapping.csv")
nhanes_path <- file.path(project_root, "Coding", "Final_Quantile_Reg_Pipeline", "NHANES_Replication", "Output", "NHANES_Subset_Intervals.RData")

if (!file.exists(mapping_path)) {
  stop("Missing analyte mapping: ", mapping_path)
}

load(file.path(input_dir, "Harmonized_Targeted.RData"))
targeted_lists <- build_targeted_lists(Targeted_Data_All)
child_list_Targeted <- targeted_lists$child_list
mat_list_Targeted <- targeted_lists$mat_list
save(child_list_Targeted, file = file.path(output_dir, "child_list_Targeted.RData"))
save(mat_list_Targeted, file = file.path(output_dir, "mat_list_Targeted.RData"))

# Specimen matrix per study x analyte. The app previously inferred blood vs urine
# from the location string, which was wrong for every cohort not named in that
# pattern. Matrix_Class is unique within study x analyte in the harmonized data
# (checked at build time below), so a lookup is exact.
Matrix_Lookup <- Targeted_Data_All %>%
  filter(!is.na(Analyte_Code), !is.na(Matrix_Class)) %>%
  mutate(Study = std_study(Study)) %>%
  distinct(Study, Analyte_Code, Matrix_Class)
dup_matrix <- Matrix_Lookup %>% count(Study, Analyte_Code) %>% filter(n > 1)
if (nrow(dup_matrix) > 0) {
  stop("Matrix_Class is not unique within study x analyte for ", nrow(dup_matrix),
       " combinations; the app's tissue labels would be ambiguous.")
}

# Display unit per study x analyte (2026-11-19). Within a creatinine-corrected
# cohort, samples with no creatinine stay unadjusted, so a study x analyte can
# mix labels; the majority label wins and Unit_Share records how much of the
# data it covers, so the app can flag the mixed ones.
# Rows with no Matrix_Class (NHBCS PFAS, 2,990 rows) are kept: they still need a
# unit, and the app falls back to the analyte's commonest matrix for them.
Unit_Lookup <- Targeted_Data_All %>%
  filter(!is.na(Analyte_Code)) %>%
  mutate(Study = std_study(Study),
         Unit = unit_label(Units, Matrix_Class, Dilution_Adjustment)) %>%
  count(Study, Analyte_Code, Unit) %>%
  group_by(Study, Analyte_Code) %>%
  mutate(Unit_Share = n / sum(n)) %>%
  arrange(desc(n), .by_group = TRUE) %>%
  # Where >10% of samples sit on another scale (DISCOVER and MARBLES are ~half
  # SG-adjusted, half unadjusted because SG is missing), say so in the label
  # rather than present the whole cohort as one scale.
  mutate(Minor_Unit = dplyr::lead(Unit), Minor_Share = dplyr::lead(Unit_Share)) %>%
  slice_head(n = 1) %>%
  ungroup() %>%
  mutate(Unit = dplyr::case_when(
    Unit_Share >= 0.9 | is.na(Minor_Unit) ~ Unit,
    grepl(", unadjusted$", Minor_Unit) ~
      paste0(Unit, " (", round(100 * (1 - Unit_Share)), "% of samples unadjusted)"),
    sub(",.*", "", Minor_Unit) == sub(",.*", "", Unit) ~
      paste0(Unit, " (", round(100 * (1 - Unit_Share)), "% of samples ", sub("^[^,]+, ", "", Minor_Unit), ")"),
    TRUE ~ paste0(Unit, " (", round(100 * (1 - Unit_Share)), "% of samples in ", Minor_Unit, ")")
  )) %>%
  select(Study, Analyte_Code, Unit, Unit_Share)
mixed_units <- Unit_Lookup %>% filter(Unit_Share < 0.95)
if (nrow(mixed_units)) {
  message("Study x analyte with <95% of rows on the majority unit (", nrow(mixed_units), "): ",
          paste0(mixed_units$Study, "/", mixed_units$Analyte_Code, " ",
                 round(100 * mixed_units$Unit_Share), "%", collapse = ", "))
}
message("Display units in use: ", paste(sort(unique(Unit_Lookup$Unit)), collapse = " | "))
Matrix_Lookup <- Matrix_Lookup %>% full_join(Unit_Lookup, by = c("Study", "Analyte_Code"))
save(Matrix_Lookup, file = file.path(output_dir, "Matrix_Lookup.RData"))

load(file.path(input_dir, "Harmonized_SDOH.RData"))
location_sdoh <- build_location_sdoh(SDOH_Data_All, input_dir)
app_sdoh <- build_app_sdoh(SDOH_Data_All)

interval_refs <- build_reference_intervals(output_dir, pooled_dir)
build_location_bundle(output_dir, interval_refs, targeted_lists, location_sdoh, nhanes_path)
build_sdoh_plot_dfs(output_dir, interval_refs, targeted_lists, app_sdoh)
build_ridge_plot_data(output_dir, targeted_lists, location_sdoh)
build_mat_child_data(output_dir, targeted_lists, input_dir, mapping_path)
write_support_files(output_dir, project_root, mapping_path)

message("Finished writing app bundle files.")
