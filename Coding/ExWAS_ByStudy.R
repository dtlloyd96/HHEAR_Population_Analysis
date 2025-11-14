# Load packages quietly ---------------------------------------------------
suppressPackageStartupMessages({
  library(tidyverse)   # dplyr, tidyr, purrr, stringr, tibble, etc.
  library(quantreg)    # quantile regression
  library(broom)       # tidy() for model outputs / CIs
  library(stringr)
  library(tibble)
})
setwd("/Users/dillon/Library/CloudStorage/GoogleDrive-2015.dillon@gmail.com/My Drive/Chirag_Patel/HHEAR_Analysis/HHEAR_Manuscript1")
# Load harmonized inputs --------------------------------------------------
load("../Input/Harmonized_Datasets/Harmonized_SDOH.RData")
load("../Input/Harmonized_Datasets/Harmonized_Targeted.RData")


# Helper to standardize Study labels by stripping suffixes like "_XYZ" or ":" ----
std_study <- function(x) {
  stringr::str_remove(x, '[_]\\w+|:')  # CHECK: ensure this collapse is intended
}

# Safely set a single name at position i (collapses multi-length labels)
safe_name_i <- function(lst, i, label) {
  nm <- names(lst)
  if (length(label) != 1L) {
    label <- paste(unique(as.character(label)), collapse = "|")
  }
  if (is.null(nm)) nm <- character(length(lst))
  nm[i] <- label
  names(lst) <- nm
  lst
}

# Safely fix names after subsetting a list (lengths must match)
fix_names_after_subset <- function(lst, desired_names = NULL, prefix = "Study") {
  if (is.null(names(lst))) names(lst) <- rep("", length(lst))
  if (!is.null(desired_names)) {
    desired_names <- as.character(desired_names)
    desired_names <- desired_names[seq_along(lst)]  # trim if too long
    if (length(desired_names) < length(lst)) {
      desired_names <- c(
        desired_names,
        paste0(prefix, "_", seq_len(length(lst) - length(desired_names)))
      )
    }
    names(lst) <- desired_names
  } else {
    # make sure names vector exists & matches length
    names(lst) <- ifelse(
      nzchar(names(lst)),
      names(lst),
      paste0(prefix, "_", seq_along(lst))
    )
  }
  # make unique to avoid duplicates after std_study collapses
  names(lst) <- make.unique(names(lst), sep = "_")
  lst
}


# Build SDOH exposure dataframe with helper recodes -----------------------
ExWAS_SDOH_Exposures <- SDOH_Data_All %>%
  mutate(
    # Prefer Max_Edu but fall back to Mat_Edu
    Edu_Any = dplyr::coalesce(Max_Edu, Mat_Edu),
    # Prefer Child_Race but fall back to Mat_Race
    Mat_Child_Comb_Race = dplyr::coalesce(Child_Race, Mat_Race),
    Mat_Age = suppressWarnings(as.numeric(Mat_Age)),  # NOTE: coercion; NAs may be introduced
    # Make Location a factor with reference "United_States"
    Location = stats::relevel(as.factor(Location), ref = "United_States"),
    # Collapse various US-specific locations to a single "United_States"
    Location_Country = dplyr::case_when(
      Location %in% c(
        "United_States","NorthEast_US","New_Hampshire","Washington_State",
        "Midwestern_US","Denver","Baltimore","Sacramento",
        "Southern_California","Syracuse","Michigan","California",
        "SanFrancisco","Massachusetts"
      ) ~ "United_States",
      TRUE ~ "Worldwide"
    ),
    # Binary maternal education bucket
    Mat_Edu_Binary = dplyr::case_when(
      Mat_Edu %in% c("No_Education","Primary_Education","Primary_School","High_School") ~ "No_College",
      Mat_Edu %in% c("Some_College","Four_Year_College_Degree","Some_Grad_School") ~ "Some_College",
      Mat_Edu == "Not_Reported" ~ NA_character_,
      TRUE ~ NA_character_  # CHECK: any unexpected levels?
    ),
    # Binary maternal income bucket
    Mat_Income_Binary = dplyr::case_when(
      Mat_Income %in% c("0_25k","25_49k") ~ "Under_49k",
      Mat_Income %in% c("50_74k","75_99k","100_124k","Over_100k","Over_125k","Over_49k") ~ "Over_49k",
      TRUE ~ NA_character_
    )
  )

# Unit conversion helper: convert to ng/mL --------------------------------
to_ng_per_mL <- function(value, units) {
  dplyr::case_when(
    is.na(units)          ~ value,          # assume already ng/mL if missing
    units == "ng/mL"      ~ value,
    units == "ug/dL"      ~ value * 10,     # 1 ug/dL = 10 ng/mL
    units %in% c("ng/L")  ~ value / 1000,   # 1 ng/L = 0.001 ng/mL
    units %in% c("pg/mL","pg/ml") ~ value / 1000,  # 1000 pg/mL = 1 ng/mL
    units == "mg/L"       ~ value * 1000,   # 1 mg/L = 1000 ng/mL
    TRUE ~ value          # CHECK: any other units present?
  )
}

# CHILD targeted data: long -> wide per-analyte list ----------------------
Child_Data_Targeted <- Targeted_Data_All %>%
  filter(!is.na(Child_PID)) %>%                           # keep child samples
  mutate(
    Units = if_else(is.na(Units), "ng/mL", Units),        # default units
    Concentration = to_ng_per_mL(Concentration, Units)    # normalize units
  ) %>%
  select(-Mat_PID) %>%                                    # not needed for child
  filter(!is.na(Analyte_Code)) %>%                        # require analyte code
  group_by(Study, Child_PID, Analyte_Code) %>%            # average replicates
  summarise(Concentration = mean(Concentration, na.rm = TRUE), .groups = "drop")

child_analytes <- sort(unique(Child_Data_Targeted$Analyte_Code))  # vector of analytes

# Build a list of per-analyte wide data.frames keyed by Child_PID ---------
child_list_Targeted <- purrr::map(
  child_analytes,
  ~ Child_Data_Targeted %>%
    filter(Analyte_Code == .x) %>%                         # one analyte
    mutate(Study = std_study(Study)) %>%                   # standardize label
    select(Study, Child_PID, Analyte_Code, Concentration) %>%
    tidyr::pivot_wider(names_from = Analyte_Code, values_from = Concentration) %>%
    select(Child_PID, Study, all_of(.x))                   # only current analyte col
) %>%
  rlang::set_names(child_analytes)

# MATERNAL targeted data: long -> wide per-analyte list -------------------
Mat_Data_Targeted <- Targeted_Data_All %>%
  filter(!is.na(Mat_PID)) %>%                              # keep maternal samples
  select(-Child_PID) %>%
  mutate(
    Units = if_else(is.na(Units), "ng/mL", Units),
    Concentration = to_ng_per_mL(Concentration, Units)
  ) %>%
  filter(!is.na(Analyte_Code)) %>%
  group_by(Study, Mat_PID, Analyte_Code) %>%
  summarise(Concentration = mean(Concentration, na.rm = TRUE), .groups = "drop") %>%
  mutate(across(everything(), ~ replace(., is.infinite(.), NA_real_)))  

mat_analytes <- sort(unique(Mat_Data_Targeted$Analyte_Code))

mat_list_Targeted <- purrr::map(
  mat_analytes,
  ~ Mat_Data_Targeted %>%
    filter(Analyte_Code == .x) %>%
    mutate(Study = std_study(Study)) %>%
    select(Study, Mat_PID, Analyte_Code, Concentration) %>%
    tidyr::pivot_wider(names_from = Analyte_Code, values_from = Concentration) %>%
    select(Mat_PID, Study, all_of(.x))
) %>%
  rlang::set_names(mat_analytes)

# Split SDOH into child and maternal panels -------------------------------
Child_SDOH <- ExWAS_SDOH_Exposures %>%
  filter(!is.na(Child_PID)) %>%
  select(Study, Child_PID, Mat_PID, dplyr::everything())

Mat_SDOH <- ExWAS_SDOH_Exposures %>%
  mutate(Adult_Sex = if_else(is.na(Adult_Sex), "F", Adult_Sex)) %>%  # default to F if missing
  filter(Adult_Sex != "M", !is.na(Mat_PID)) %>%                      # keep likely maternal rows
  select(Study, Child_PID, Mat_PID, dplyr::everything())

# Core quantile regression wrapper ---------------------------------------
Quantile_Reg <- function(exp_var = NULL,
                         dataset,
                         my_outcome,
                         Reg_Type,         # kept for API compatibility
                         my_covariates,
                         join_var) {
  
  # Build analysis df (with/without explicit exposure predictor)
  if (is.null(exp_var)) {
    mydf <- dataset %>%
      select(all_of(c(join_var, my_outcome, my_covariates))) %>%
      ungroup() %>%
      filter(complete.cases(.))                                # complete-case analysis
  } else {
    mydf <- dataset %>%
      select(all_of(c(join_var, my_outcome, exp_var, my_covariates))) %>%
      filter(!is.na(.data[[exp_var]])) %>%
      filter(complete.cases(.))
  }
  
  if (!nrow(mydf)) return(NULL)                                # early exit if empty
  
  # If Study present, set most frequent level as reference
  if ("Study" %in% names(mydf)) {
    max_category <- names(which.max(table(mydf$Study)))
    mydf <- mydf %>%
      mutate(Study = relevel(as.factor(Study), ref = max_category))
  }
  
  # Drop groups with no variation (either covariates or exposure), by Study
  if ("Study" %in% names(mydf)) {
    if (is.null(exp_var)) {
      tocheck <- setdiff(my_covariates, "Study")
      clean_df <- mydf %>%
        group_by(Study) %>%
        filter(if_all(all_of(tocheck), ~ dplyr::n_distinct(.) > 1)) %>%
        ungroup()
    } else {
      clean_df <- mydf %>%
        group_by(Study) %>%
        filter(dplyr::n_distinct(.data[[exp_var]]) > 1) %>%
        ungroup() %>%
        droplevels()
    }
  } else {
    clean_df <- mydf
  }
  
  if (!nrow(clean_df)) return(NULL)                            # nothing left to model
  
  # Quantiles to estimate
  quantiles <- c(0.05, 0.25, 0.5, 0.75, 0.95)
  
  # Build model formula: include Study if multiple studies remain
  if (is.null(exp_var)) {
    myform <- stats::as.formula(paste(my_outcome, "~", paste(my_covariates, collapse = "+")))
  } else {
    n_study <- if ("Study" %in% names(clean_df)) dplyr::n_distinct(clean_df$Study) else 0
    if (n_study <= 1) {
      myform <- stats::as.formula(paste(my_outcome, "~", exp_var))
    } else {
      myform <- stats::as.formula(paste(my_outcome, "~", paste(c(exp_var, my_covariates), collapse = "+")))
    }
  }
  
  # Fit full and null quantile regressions
  rq_model <- quantreg::rq(myform, data = clean_df, tau = quantiles)
  fit_null <- quantreg::rq(stats::as.formula(paste(my_outcome, "~ 1")), data = clean_df, tau = quantiles)
  
  # Predicted values (aligned with clean_df rows)
  yhat_full <- predict(rq_model)
  yhat_null <- predict(fit_null)
  y <- clean_df[[my_outcome]]  
  
  # Koenker's rho loss for pseudo-R1
  rho_tau <- function(u, tau) u * (tau - (u < 0))
  
  # Compute pseudo-R1 per quantile tau
  r1_all_df <- purrr::map2_dfr(
    seq_along(quantiles), quantiles,
    ~{
      tau <- .y
      u_full <- y - yhat_full[, .x]
      u_null <- y - yhat_null[, .x]
      num <- sum(rho_tau(u_full, tau), na.rm = TRUE)
      den <- sum(rho_tau(u_null, tau), na.rm = TRUE)
      tibble::tibble(Tau = tau, pseudo_R = if (den == 0) NA_real_ else 1 - num/den)
    }
  )
  
  # Bootstrap SEs for coefficients (R small to keep runtime comparable)
  ctab <- try(summary(rq_model, se = "boot", R = 2, conf.int = TRUE), silent = TRUE)  # CHECK: raise R for inference
  if (inherits(ctab, "try-error")) return(NULL)
  
  # Stack coefficient tables across taus, add Outcome label
  ctable_final <- purrr::map_dfr(ctab, function(x) {
    as.data.frame(coef(x)) %>%
      rownames_to_column("Variable") %>%
      mutate(Tau = x$tau, Outcome = my_outcome)
  })
  names(ctable_final) <- c("Variable","Beta","Standard_Error","T_Value","P_Value","Tau","Outcome")
  
  # Bind CI columns from tidy()
  conf_int_table <- broom::tidy(rq_model, conf.int = TRUE) %>%
    select(Lower = conf.low, Upper = conf.high)
  ctable_final <- dplyr::bind_cols(ctable_final, conf_int_table)
  
  # Add Exposure label when provided
  if (!is.null(exp_var)) {
    ctable_final <- ctable_final %>% mutate(Exposure = exp_var, .before = 1)
  }
  
  # Attach pseudo-R1 by tau
  ctable_final <- dplyr::left_join(ctable_final, r1_all_df, by = "Tau")
  
  # Return results + the analysis dataframe used
  list(
    Regression_Results = ctable_final,
    Quantile_Results   = list(data = clean_df)
  )
}

# Utility: silence console/message/warnings for wrapped expression --------
suppress_all_output <- function(expr) {
  tf <- tempfile()
  con <- file(tf, open = "wt")
  sink(con); sink(con, type = "message")  # divert stdout and messages
  ow <- getOption("warn"); options(warn = -1)
  on.exit({
    options(warn = ow)
    sink(NULL); sink(NULL, type = "message")
    close(con)
  }, add = TRUE)
  suppressWarnings(suppressMessages(tryCatch(eval(expr), error = function(e) NULL)))
}

# Run ExWAS for a set of outcomes within a single study -------------------
ExWAS_Function_Metab <- function(Exposure_Data, myexposures, myoutcomes,
                                 Outcome_Data, mycovars, join_var, Reg_Type, curr_study) {
  out_list <- vector("list", length(myoutcomes))
  names(out_list) <- myoutcomes
  
  for (j in seq_along(myoutcomes)) {
    pheno <- myoutcomes[j]
    
    # Outcome table for this outcome/study
    curr_pheno <- Outcome_Data[[j]] %>%
      mutate(Study = std_study(Study)) %>%
      select(all_of(c(join_var, "Study", pheno))) %>%
      filter(Study %in% curr_study) %>%
      filter(!is.na(.data[[pheno]])) %>%
      mutate(!!pheno := as.numeric(.data[[pheno]]))  # numeric outcome
    
    if (!nrow(curr_pheno)) {
      out_list[[j]] <- list()  # nothing to analyze for this outcome
      next
    }
    
    if (is.null(myexposures)) {
      # Covariates-only model (no explicit exposure predictor)
      qe <- Quantile_Reg(
        exp_var      = NULL,
        dataset      = curr_pheno %>% left_join(
          Exposure_Data %>% mutate(Study = std_study(Study)) %>% distinct(),
          by = c(join_var, "Study")
        ),
        my_outcome   = pheno,
        Reg_Type     = Reg_Type,
        my_covariates= mycovars,
        join_var     = join_var
      )
      out_list[[j]] <- list(CovOnly = qe)
      next
    }
    
    # Otherwise: loop exposures as predictors
    res_j <- vector("list", length(myexposures))
    names(res_j) <- myexposures
    
    for (i in seq_along(myexposures)) {
      curr_exp <- myexposures[i]
      
      # Merge covariates + current exposure
      curr_exposure <- Exposure_Data %>%
        mutate(Study = std_study(Study)) %>%
        select(all_of(c(join_var, "Study", curr_exp, mycovars))) %>%
        distinct()
      
      dat <- curr_pheno %>%
        left_join(curr_exposure, by = c(join_var, "Study")) %>%
        distinct()
      
      if (!nrow(dat)) {
        res_j[[i]] <- NULL
        next
      }
      
      # Run quantile regression for this exposure-outcome pair
      qe <- Quantile_Reg(
        exp_var       = curr_exp,
        dataset       = dat,
        my_outcome    = pheno,
        Reg_Type      = Reg_Type,
        my_covariates = mycovars,
        join_var      = join_var
      )
      res_j[[i]] <- qe
    }
    
    # Keep successful fits only
    res_j <- res_j[!purrr::map_lgl(res_j, is.null)]
    out_list[[j]] <- res_j
  }
  
  out_list
}

# Study-level runner: covariate-only models per study ---------------------
Run_ExWAS_By_Study <- function(to_keep, sdoh_df, mystudies, join_var, analytes_touse, targeted_list) {
  all_res <- vector("list", length(mystudies))
  
  for (i in seq_along(mystudies)) {
    # Keep study rows and drop all-NA or no-variation cols
    curr_sdoh <- sdoh_df %>%
      filter(Study == mystudies[i]) %>%
      select(where(~ !all(is.na(.)))) %>%
      select(where(~ dplyr::n_distinct(.) > 1), Study) %>%
      mutate(Study = std_study(Study))
    
    # Keep covariates requested and present in this study
    curr_social_exps <- intersect(setdiff(names(curr_sdoh), "Study"), to_keep)
    curr_study <- unique(curr_sdoh$Study)
    
    # Run ExWAS with covariates-only (myexposures=NULL)
    res <- try(ExWAS_Function_Metab(
      Exposure_Data = curr_sdoh,
      Outcome_Data  = targeted_list,
      myexposures   = NULL,
      myoutcomes    = analytes_touse,
      mycovars      = curr_social_exps,
      join_var      = join_var,
      Reg_Type      = "Quantile",
      curr_study    = curr_study
    ), silent = TRUE)
    
    all_res[[i]] <- if (inherits(res, "try-error")) list() else res
    names(all_res)[i] <- curr_study
  }
  
  # Keep non-empty studies and drop unexpected character elements
  out <- all_res[purrr::map_int(all_res, length) > 0]
  out <- out[!vapply(out, is.character, logical(1))]
  out
}

# Save helper: ensure directory exists -----------------------------------
save_rdata <- function(obj, path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  save(obj, file = path)
}

# Wrapper to run CHILD covariate sets across studies ----------------------
run_child_models <- function(to_keep, label) {
  mystudies <- sort(unique(Child_SDOH$Study))
  res <- suppress_all_output(quote(
    Run_ExWAS_By_Study(
      to_keep        = to_keep,
      sdoh_df        = Child_SDOH,
      mystudies      = mystudies,
      join_var       = "Child_PID",
      analytes_touse = child_analytes,
      targeted_list  = child_list_Targeted
    )
  ))
  res <- res[purrr::map_int(res, length) > 0]
  save_rdata(res, sprintf("Output/%s.RData", label))
  invisible(res)
}

# Wrapper to run MATERNAL sets across studies -----------------------------
run_mat_models <- function(to_keep, label, covonly = FALSE) {
  mystudies <- sort(unique(Mat_SDOH$Study))
  
  if (covonly) {
    # Covariate-only (no explicit exposure predictor)
    res <- suppress_all_output(quote(
      Run_ExWAS_By_Study(
        to_keep        = to_keep,
        sdoh_df        = Mat_SDOH,
        mystudies      = mystudies,
        join_var       = "Mat_PID",
        analytes_touse = mat_analytes,
        targeted_list  = mat_list_Targeted
      )
    ))
  } else {
    # Full SDOH sets per study (myexposures=NULL + mycovars = present SDOH)
    all_res <- vector("list", length(mystudies))
    for (i in seq_along(mystudies)) {
      curr_sdoh <- Mat_SDOH %>%
        filter(Study == mystudies[i]) %>%
        select(where(~ !all(is.na(.)))) %>%
        select(where(~ dplyr::n_distinct(.) > 1), Study) %>%
        mutate(Study = std_study(Study))
      
      curr_social_exps <- intersect(
        setdiff(names(curr_sdoh), "Study"),
        c("Mat_Edu","Mat_Income","Mat_Race","Mat_Age","Location","Site_Location")
      )
      curr_study <- unique(curr_sdoh$Study)
      
      res_i <- suppress_all_output(quote(ExWAS_Function_Metab(
        Exposure_Data = curr_sdoh,
        Outcome_Data  = mat_list_Targeted,
        myexposures   = NULL,
        myoutcomes    = mat_analytes,
        mycovars      = curr_social_exps,
        join_var      = "Mat_PID",
        Reg_Type      = "Quantile",
        curr_study    = curr_study
      )))
      
      all_res[[i]] <- res_i[purrr::map_int(res_i, length) > 0]
      names(all_res)[i] <- curr_study
    }
    res <- all_res[purrr::map_int(all_res, length) > 0]
  }
  
  res <- res[purrr::map_int(res, length) > 0]
  save_rdata(res, sprintf("Output/%s.RData", label))
  invisible(res)
}

# CHILD: per-study single-exposure models (exposures as predictors, Study as covariate)
mystudies <- sort(unique(Child_SDOH$Study))
All_Study_Single_Exp_Results <- vector("list", length(mystudies))
for (i in seq_along(mystudies)) {
  curr_sdoh <- Child_SDOH %>%
    filter(Study == mystudies[i]) %>%
    select(where(~ !all(is.na(.)))) %>%
    select(where(~ dplyr::n_distinct(.) > 1), Study) %>%
    mutate(Study = std_study(Study))
  
  tokeep <- c("Child_Age","Child_Race","Child_Sex","Mat_Edu","Mat_Income","Mat_Race","Mat_Age","Location")
  curr_social_exps <- names(curr_sdoh)[-1]
  curr_social_exps <- curr_social_exps[curr_social_exps %in% tokeep]  # only desired vars present
  curr_study <- unique(curr_sdoh$Study)
  
  Curr_Results <- suppress_all_output(quote(ExWAS_Function_Metab(
    Exposure_Data = curr_sdoh,
    Outcome_Data  = child_list_Targeted,
    myexposures   = curr_social_exps,  # exposures loop
    myoutcomes    = child_analytes,
    mycovars      = "Study",           # adjust for Study when multiple
    join_var      = "Child_PID",
    Reg_Type      = "Quantile",
    curr_study    = curr_study
  )))
  
  Final_Res <- Curr_Results[purrr::map_int(Curr_Results, length) > 0]
  All_Study_Single_Exp_Results[[i]] <- Final_Res
  names(All_Study_Single_Exp_Results)[i] <- curr_study
}
By_Study_Single_Exp_Results_Child <- All_Study_Single_Exp_Results[purrr::map_int(All_Study_Single_Exp_Results, length) > 0]
save(By_Study_Single_Exp_Results_Child, file = "Output/By_Study_Single_Exp_Results_Child.RData")

# CHILD: covariate-set runs per study ------------------------------------
run_child_models(to_keep = c("Child_Age","Child_Sex"), label = "By_Study_AgeSex_Child")
run_child_models(to_keep =c("Child_Age","Child_Sex","Mat_Child_Comb_Race"),label = "By_Study_AgeSexRace_Child")
run_child_models(to_keep =c("Child_Age","Child_Sex","Mat_Edu"), label ="By_Study_AgeSexEdu_Child")
run_child_models(to_keep =c("Child_Age","Child_Sex","Mat_Child_Comb_Race","Mat_Edu"),label = "By_Study_AgeSexRaceEdu_Child")

# CHILD: full SDOH within each study -------------------------------------
{
  mystudies <- sort(unique(Child_SDOH$Study))
  all_res <- list()
  
  for (i in seq_along(mystudies)) {
    print(i)
    curr_sdoh <- Child_SDOH %>%
      filter(Study == mystudies[i]) %>%
      select(where(~ !all(is.na(.)))) %>%
      select(where(~ dplyr::n_distinct(.) > 1), Study) %>%
      mutate(Study = std_study(Study))
    
    tokeep <- c("Child_Age","Child_Race","Child_Sex","Mat_Edu","Mat_Income","Mat_Race","Mat_Age","Location")
    curr_social_exps <- intersect(setdiff(names(curr_sdoh), "Study"), tokeep)  # keep present vars only
    curr_study <- unique(curr_sdoh$Study)
    
    res_i <- suppress_all_output(quote(ExWAS_Function_Metab(
      Exposure_Data = curr_sdoh,
      Outcome_Data  = child_list_Targeted,
      myexposures   = NULL,              # covariates-only model
      myoutcomes    = child_analytes,
      mycovars      = curr_social_exps,
      join_var      = "Child_PID",
      Reg_Type      = "Quantile",
      curr_study    = curr_study
    )))
    
    all_res[[i]] <- res_i[purrr::map_int(res_i, length) > 0]
    names(all_res)[i] <- curr_study
  }
  
  By_Study_AllSDOH_Child <- all_res[purrr::map_int(all_res, length) > 0]
  save(By_Study_AllSDOH_Child, "Output/By_Study_AllSDOH_Child.RData")
}

# MATERNAL: per-study single-exposure models ------------------------------
{
  mystudies <- sort(unique(Mat_SDOH$Study))
  all_res <- vector("list", length(mystudies))
  
  for (i in seq_along(mystudies)) {
    curr_sdoh <- Mat_SDOH %>%
      filter(Study == mystudies[i]) %>%
      select(where(~ !all(is.na(.)))) %>%
      select(where(~ dplyr::n_distinct(.) > 1), Study) %>%
      mutate(Study = std_study(Study))
    
    tokeep <- c("Mat_Edu","Mat_Income","Mat_Race","Mat_Age","Location","Site_Location")
    curr_social_exps <- intersect(setdiff(names(curr_sdoh), "Study"), tokeep)
    curr_study <- unique(curr_sdoh$Study)
    
    res_i <- suppress_all_output(quote(ExWAS_Function_Metab(
      Exposure_Data = curr_sdoh,
      Outcome_Data  = mat_list_Targeted,
      myexposures   = curr_social_exps,  # exposures loop
      myoutcomes    = mat_analytes,
      mycovars      = "Study",
      join_var      = "Mat_PID",
      Reg_Type      = "Quantile",
      curr_study    = curr_study
    )))
    
    all_res[[i]] <- res_i[purrr::map_int(res_i, length) > 0]
    names(all_res)[i] <- curr_study
  }
  
  By_Study_Single_Exp_Results_Mat <- all_res[purrr::map_int(all_res, length) > 0]
  save(By_Study_Single_Exp_Results_Mat, file = "Output/By_Study_Single_Exp_Results_Mat.RData")
}

# MATERNAL: full SDOH and covariate sets ---------------------------------
run_mat_models(to_keep =NULL, label = "By_Study_AllSDOH_Mat", covonly = FALSE)
run_mat_models(to_keep =c("Mat_Age","Mat_Race"), label =  "By_Study_AgeRace_Mat", covonly = TRUE)
run_mat_models(to_keep =c("Mat_Age","Mat_Edu"),  label =  "By_Study_AgeEdu_Mat",  covonly = TRUE)
run_mat_models(to_keep =c("Mat_Age","Mat_Race","Mat_Edu"), label =  "By_Study_AgeRaceEdu_Mat", covonly = TRUE)
