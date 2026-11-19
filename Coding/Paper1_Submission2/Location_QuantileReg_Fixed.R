library(tidyverse)
library(quantreg)
library(survival)
# PORTABLE WORKING DIRECTORY (2026-08-11, CODE_AUDIT_2026-08-11.md section 3.2).
# This was a hardcoded setwd("/home/<user>/HHEAR/HHEAR_Quantreg/Coding"), which
# exists only on the cluster account -- the script could not run anywhere else.
# Resolve the project root by walking up to the directory that contains both
# Coding/ and Input/, then chdir to its Coding/ subdirectory, which is the
# working directory the relative paths below assume.
local({
  argv <- commandArgs(trailingOnly = FALSE)
  fa <- argv[grep("^--file=", argv)]
  start <- if (length(fa)) dirname(normalizePath(sub("^--file=", "", fa[[1]]), mustWork = FALSE)) else getwd()
  cur <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (dir.exists(file.path(cur, "Coding")) && dir.exists(file.path(cur, "Input"))) break
    parent <- dirname(cur)
    if (identical(parent, cur)) stop("Could not locate the project root (a dir containing Coding/ and Input/) from ", start)
    cur <- parent
  }
  setwd(file.path(cur, "Coding"))
})

# The harmonized .RData live directly in Input/ on the cluster but in
# Input/Harmonized_Datasets/ locally. Resolve either layout.
harmonized <- function(f) {
  cands <- c(file.path("../Input", f), file.path("../Input/Harmonized_Datasets", f))
  hit <- cands[file.exists(cands)]
  if (!length(hit)) stop("Cannot find ", f, " in ../Input or ../Input/Harmonized_Datasets")
  hit[[1]]
}

source("Paper1_Submission2/censored_quantile_helpers.R")

#load("../../Input/Harmonized_Datasets/Harmonized_SDOH.RData")
#load("../../Input/Harmonized_Datasets/Harmonized_Targeted.RData")


load(harmonized("Harmonized_SDOH.RData"))
load(harmonized("Harmonized_Targeted.RData"))

std_study <- function(x) {
  stringr::str_remove(x, '[_]\\w+|:')  # CHECK: ensure this collapse is intended
}

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

# UNIT CONVERSION (2026-08-11, CODE_AUDIT_2026-08-11.md section 1.1).
# The local to_ng_per_mL() removed here ended in `TRUE ~ value`, so mg/dL,
# ug/mL, ng/g, pmol/mL, nM and kg/m3 fell through UNCONVERTED and were treated
# as ng/mL -- 39,475 rows (5.9%), wrong by factors of 10^3 to 10^6.
local({
  cands <- c("harmonization_helpers.R",
             "Paper1_Submission2/harmonization_helpers.R",
             "Coding/Paper1_Submission2/harmonization_helpers.R")
  hit <- cands[file.exists(cands)]
  if (length(hit) == 0) stop("Cannot locate harmonization_helpers.R")
  source(hit[[1]], local = FALSE)
})
# OPTIONAL SPEC RESTRICTION (added 2026-08-11).
# The full sweep fits 20 covariate specifications x ~217 analytes x 2 panels,
# which is the bulk of the runtime. Set the HHEAR_SPECS environment variable to a
# comma-separated list of covariate-set names to fit only those, e.g.
#
#     HHEAR_SPECS=Location Rscript Paper1_Submission2/Location_QuantileReg_Fixed.R
#
# Unset (the default) runs everything, so this changes nothing unless asked for.
SPECS_TO_RUN <- {
  s <- Sys.getenv("HHEAR_SPECS", "")
  if (nzchar(s)) trimws(strsplit(s, ",")[[1]]) else NULL
}
if (!is.null(SPECS_TO_RUN)) {
  message("[spec filter] restricted to: ", paste(SPECS_TO_RUN, collapse = ", "))
}

Child_Data_Targeted <- build_censored_targeted_panel(
  targeted_df = Targeted_Data_All,
  id_var = "Child_PID",
  unit_converter = to_ng_per_mL
)

child_analytes <- sort(unique(Child_Data_Targeted$Analyte_Code))  # vector of analytes

# Build a list of per-analyte wide data.frames keyed by Child_PID ---------
child_list_Targeted <- build_censored_targeted_list(
  panel_df = Child_Data_Targeted,
  id_var = "Child_PID",
  analytes = child_analytes
)

# MATERNAL targeted data: long -> wide per-analyte list -------------------
Mat_Data_Targeted <- build_censored_targeted_panel(
  targeted_df = Targeted_Data_All,
  id_var = "Mat_PID",
  unit_converter = to_ng_per_mL
)

mat_analytes <- sort(unique(Mat_Data_Targeted$Analyte_Code))

mat_list_Targeted <- build_censored_targeted_list(
  panel_df = Mat_Data_Targeted,
  id_var = "Mat_PID",
  analytes = mat_analytes
)

# SAMPLE-YEAR JOIN (fixed 2026-08-11, CODE_AUDIT_2026-08-11.md section 1.3).
# This was group_by(PID, Sample_Collection_Year) %>% summarise(n = n()), i.e. ONE
# ROW PER (PID, YEAR), then left_join()ed on PID alone. 2,355 Child_PIDs and
# 1,569 Mat_PIDs carry more than one collection year, so the join was
# many-to-many and multiplied every affected participant downstream. Collapsed to
# one row per participant using the median year, matching the ExWAS and the power
# arrays. The `n` column is gone, so the paired dplyr::select(-n) was dropped too.
Child_Sample_Yr <- Targeted_Data_All %>%
  filter(!is.na(Child_PID)) %>%
  group_by(Child_PID) %>%
  summarise(
    Sample_Collection_Year = {
      yr <- suppressWarnings(as.numeric(Sample_Collection_Year))
      if (all(is.na(yr))) NA_real_ else stats::median(yr, na.rm = TRUE)
    },
    .groups = "drop")

# SAMPLE-YEAR JOIN (fixed 2026-08-11, CODE_AUDIT_2026-08-11.md section 1.3).
# This was group_by(PID, Sample_Collection_Year) %>% summarise(n = n()), i.e. ONE
# ROW PER (PID, YEAR), then left_join()ed on PID alone. 2,355 Child_PIDs and
# 1,569 Mat_PIDs carry more than one collection year, so the join was
# many-to-many and multiplied every affected participant downstream. Collapsed to
# one row per participant using the median year, matching the ExWAS and the power
# arrays. The `n` column is gone, so the paired dplyr::select(-n) was dropped too.
Mat_Sample_Yr <- Targeted_Data_All %>%
  filter(!is.na(Mat_PID)) %>%
  group_by(Mat_PID) %>%
  summarise(
    Sample_Collection_Year = {
      yr <- suppressWarnings(as.numeric(Sample_Collection_Year))
      if (all(is.na(yr))) NA_real_ else stats::median(yr, na.rm = TRUE)
    },
    .groups = "drop")
# Split SDOH into child and maternal panels -------------------------------
Child_SDOH <- ExWAS_SDOH_Exposures %>%
  filter(!is.na(Child_PID)) %>%
  select(Study, Child_PID, Mat_PID, dplyr::everything()) %>% 
  left_join(.,Child_Sample_Yr, by = c("Child_PID")) %>% 

  mutate(Sample_Collection_Year = as.numeric(Sample_Collection_Year))

Mat_SDOH <- ExWAS_SDOH_Exposures %>%
  filter(!is.na(Mat_PID)) %>%                      # keep likely maternal rows
  select(Study, Child_PID, Mat_PID, dplyr::everything()) %>% 
  left_join(.,Mat_Sample_Yr, by = c("Mat_PID")) %>% 

  mutate(Sample_Collection_Year = as.numeric(Sample_Collection_Year))


# BOOTSTRAP REPLICATES FOR INFERENCE (pinned 2026-08-11).
# CODE_AUDIT_2026-08-11.md section 1.5.
#
# This wrapper did not pass summary_R, so run_censored_quantile_reg() used its
# default of 2. That is below fit_censored_quantile_models()'s min_boot_R = 25,
# so inference fell back to fallback_crq_summary() and EVERY Standard_Error,
# T_Value, P_Value, Lower and Upper came back NA. Downstream that failed
# silently rather than loudly: p.adjust(NA) is NA, and
# ifelse(NA <= 0.10, 1, 0) is NA, so significance-filtered figures emptied
# instead of erroring.
#
# RUNTIME: 200 replicates is materially slower than the old no-op default,
# because the bootstrap now actually runs. Lower it for exploratory sweeps, but
# anything below 25 silently disables inference entirely -- do not treat values
# in that range as "a bit less precise".
# FAST-REGENERATION SETTING (2026-08-11): 50 bootstrap replicates, not 200.
# Enough for point estimates and an indicative CI, NOT enough for a publication
# standard error -- raise to 200+ for the final run. Below 25 inference silently
# becomes NA (see the note above), so do not lower this further.
SUMMARY_R <- 50

Quantile_Reg <- function(exp_var = NULL,
                         dataset,
                         my_outcome,
                         Reg_Type,
                         my_covariates,
                         join_var,
                         summary_R = SUMMARY_R) {
  stopifnot(summary_R >= 25)
  run_censored_quantile_reg(
    exp_var = exp_var,
    dataset = dataset,
    my_outcome = my_outcome,
    my_covariates = my_covariates,
    join_var = join_var,
    summary_R = summary_R
  )
}



ExWAS_Function_Metab <- function(Exposure_Data, myexposures,myoutcomes, Outcome_Data,mycovars,join_var,Reg_Type){
  
  all_pheno_unadjusted <- list()
  all_pheno_adjusted <- list()
  
  #Run Regression for all outcomes
  for(j in seq_along(myoutcomes)){
    
    
    pheno <- myoutcomes[j]
    adjusted_results <- list()
    unadjusted_results <- list()
    
    if(is.null(myexposures)){
      loop_n <- "CovOnly"
    }else(
      loop_n = myexposures
    )
    
    for(i in seq_along(loop_n)){
      print(i)
      my_covariates <- mycovars
      
      
      #curr_metab <- child_list_Targeted[[i]]
      # curr_pheno <- Child_Outcomes %>% 
      #   dplyr::select(Child_PID,all_of(pheno))
      
      
      curr_exposure <- Exposure_Data %>% 
        mutate(Study =  str_remove(Study, '[_]\\w+|:')) %>% 
        dplyr::select(join_var,Study,myexposures[i],all_of(my_covariates)) %>% 
        #dplyr::filter(Study == "CEHC") %>% 
        distinct()
      
      curr_pheno <- Outcome_Data[[j]] %>% 
        mutate(Study =  str_remove(Study, '[_]\\w+|:')) %>% 
        select_outcome_data(join_var = join_var, outcome_col = pheno) %>% 
        #dplyr::filter(Study == "CEHC") %>% 
        dplyr::filter(!is.na(.data[[pheno]]))
      
      curr_outcome_pheno <- left_join(curr_pheno,curr_exposure,by = c(join_var,'Study')) %>%
        mutate(!! pheno := as.numeric(!! rlang::sym(pheno))) %>%
        distinct(.keep_all = T)
      
      # curr_outcome_pheno <- left_join(curr_pheno,curr_exposure,by = c(join_var,'Study')) %>% 
      #   mutate(!! pheno := scale(as.numeric(!! rlang::sym(pheno)))) %>% 
      #   distinct(.keep_all = T)
      
      
      
      #  curr_outcome_pheno <- curr_outcome_pheno
      
      
      curr_exp <- myexposures[i]
      
      
      
      #             in_silence(
      # {
      
      Curr_Unadjusted_Results <- try(Quantile_Reg(exp_var = curr_exp,
                                                  dataset = curr_outcome_pheno,
                                                  my_outcome = pheno,
                                                  Reg_Type = Reg_Type,
                                                  my_covariates = my_covariates,
                                                  join_var = join_var))
      
      unadjusted_results[[i]] <- Curr_Unadjusted_Results
      
      names(unadjusted_results)[i] <- curr_exp
      
      #})
    }
    
    
    #filtered_unadjusted_list <- lapply(unadjusted_results, function(df) if (nrow(df) >= 2) df else NULL)
    
    
    #unadjusted_df <- bind_rows(filtered_unadjusted_list)
    
    
    
    
    all_pheno_unadjusted[[j]] <- unadjusted_results
    
    names(all_pheno_unadjusted)[j] <- myoutcomes[j]
    
  }
  
  
  return(all_pheno_unadjusted)
  
  
}


suppress_all_output <- function(expr) {
  temp_file <- tempfile()
  con <- file(temp_file, open = "wt")
  
  # Redirect stdout and messages
  sink(con)
  sink(con, type = "message")
  
  # Suppress warnings and messages globally
  old_warn <- getOption("warn")
  options(warn = -1)
  
  # Evaluate expression quietly
  result <- suppressWarnings(
    suppressMessages(
      tryCatch(
        eval(expr),
        error = function(e) NULL  # suppress error output
      )
    )
  )
  
  # Restore options and sinks
  options(warn = old_warn)
  sink(NULL)
  sink(NULL, type = "message")
  close(con)
  
  invisible(result)
}


Child_Ads <- Child_SDOH %>% 
  dplyr::select(Child_PID,Study,Location,Location_Region,Year,Center,Child_Age,Child_Sex)

Mat_Ads <- Mat_SDOH %>% 
  dplyr::select(Mat_PID,Study,Location,Location_Region,Year,Center,Mat_Age)

output_dir <- "../Output/Paper1_Output_Submission2"
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Study_Only <- c("Location")
# 
# Child_Location_Results <- suppress_all_output(quote(ExWAS_Function_Metab(Exposure_Data = Child_Ads,
#                                                                                     Outcome_Data = child_list_Targeted,
#                                                                                     myexposures = NULL,
#                                                                                     myoutcomes = child_analytes,mycovars = Study_Only,join_var = "Child_PID",
#                                                                                     Reg_Type = "Quantile")))
# 
# save(Child_Location_Results, file = file.path(output_dir, "Child_Study_Quantile_Results.RData"))

Child_Age_Sex <- suppress_all_output(quote(ExWAS_Function_Metab(Exposure_Data = Child_Ads,
                                                                                    Outcome_Data = child_list_Targeted,
                                                                                    myexposures = c("Child_Age","Child_Sex"),
                                                                                    myoutcomes = child_analytes,mycovars = "Study",join_var = "Child_PID",
                                                                                    Reg_Type = "Quantile")))

save(Child_Age_Sex, file = file.path(output_dir, "Child_Age_Sex.RData"))


Mat_Age <- suppress_all_output(quote(ExWAS_Function_Metab(Exposure_Data = Mat_Ads,
                                                                Outcome_Data = mat_list_Targeted,
                                                                myexposures = c("Mat_Age"),
                                                                myoutcomes = mat_analytes,mycovars = "Study",join_var = "Mat_PID",
                                                                Reg_Type = "Quantile")))

save(Mat_Age, file = file.path(output_dir, "Mat_Age.RData"))

###Year
my_covars_list <- list(
  "Study",
  "Year",
  "Location",
  "Location_Region",
  "Center",
  "Child_Age",
  "Child_Sex",
  c("Child_Age","Child_Sex","Study"),
  c("Location","Center"),
  c("Location","Year"),
  c("Location","Study")
)

# Give each covariate list a name (required for object naming)
# LABEL MISALIGNMENT (fixed 2026-08-11).
#
# The child list has ELEVEN entries but only TEN names were supplied, so every
# label from the 8th entry onward was attached to the wrong covariate set and the
# last one became NA:
#
#   entry  covariates                     label it received
#     8    Child_Age + Child_Sex + Study  "Location_Center"   <- wrong
#     9    Location + Center              "Location_Year"     <- wrong
#    10    Location + Year                "Location_Study"    <- wrong
#    11    Location + Study               NA  -> Child_NA_Quantile_Results.RData
#
# Those names become the output filenames, and Location_Compare.Rmd reads
# Child_Location_Center_Quantile_Results believing it holds the Location+Center
# model -- it actually held the demographic model. The child R2 comparison across
# location specifications was therefore built on mislabelled series.
#
# The maternal block below has 9 entries and 9 names and was never affected.
#
# Named inline so a future edit to the list cannot silently desynchronise them,
# plus a hard check.
names(my_covars_list) <- c("Study", "Year", "Location", "Location_Region",
                           "Center", "Child_Age", "Child_Sex",
                           "Child_Age_Sex_Study",
                           "Location_Center", "Location_Year", "Location_Study")
stopifnot(length(my_covars_list) == length(names(my_covars_list)),
          !any(is.na(names(my_covars_list))))

#-----------------------------------------------------------
# Loop through each covariate specification
#-----------------------------------------------------------
for (nm in names(my_covars_list)) {
  if (!is.null(SPECS_TO_RUN) && !(nm %in% SPECS_TO_RUN)) next

  
  covars <- my_covars_list[[nm]]
  
  message("Running model for: ", nm)
  
  # Run your function with output suppressed
  result_obj <- suppress_all_output(
    quote(
      ExWAS_Function_Metab(
        Exposure_Data = Child_Ads,
        Outcome_Data = child_list_Targeted,
        myexposures  = NULL,
        myoutcomes   = child_analytes,
        mycovars     = covars,
        join_var     = "Child_PID",
        Reg_Type     = "Quantile"
      )
    )
  )
  
  #---------------------------------------------------------
  # Assign result to an object named using the covariate label
  # e.g., Child_StudyOnly_Results
  #---------------------------------------------------------
  obj_name <- paste0("Child_", nm, "_Quantile_Results")
  assign(obj_name, result_obj)
  
  #---------------------------------------------------------
  # Save each result as its own .RData file
  #---------------------------------------------------------
  save_path <- file.path(output_dir,
                         paste0(obj_name, ".RData"))
  
  save(list = obj_name, file = save_path)
  
  message("Saved: ", save_path)
}

#Adult Results


my_covars_list <- list(
  "Study",
  "Year",
  "Location",
  "Location_Region",
  "Center",
  "Mat_Age",
  c("Location","Center"),
  c("Location","Year"),
  c("Location","Study")
  
)

# Give each covariate list a name (required for object naming)
names(my_covars_list) <- c("Study", "Year", "Location", "Location_Region",
                           "Center", "Mat_Age", "Location_Center",
                           "Location_Year", "Location_Study")
# Same guard as the child block above (2026-08-11). This one was already aligned.
stopifnot(length(my_covars_list) == length(names(my_covars_list)),
          !any(is.na(names(my_covars_list))))

for (nm in names(my_covars_list)) {
  if (!is.null(SPECS_TO_RUN) && !(nm %in% SPECS_TO_RUN)) next

  
  covars <- my_covars_list[[nm]]
  
  message("Running model for: ", nm)
  
  # Run your function with output suppressed
  result_obj <- suppress_all_output(
    quote(
      ExWAS_Function_Metab(
        Exposure_Data = Mat_Ads,
        Outcome_Data = mat_list_Targeted,
        myexposures  = NULL,
        myoutcomes   = mat_analytes,
        mycovars     = covars,
        join_var     = "Mat_PID",
        Reg_Type     = "Quantile"
      )
    )
  )
  
  #---------------------------------------------------------
  # Assign result to an object named using the covariate label
  # e.g., Child_StudyOnly_Results
  #---------------------------------------------------------
  obj_name <- paste0("Mat_", nm, "_Quantile_Results")
  assign(obj_name, result_obj)
  
  #---------------------------------------------------------
  # Save each result as its own .RData file
  #---------------------------------------------------------
  save_path <- file.path(output_dir,
                         paste0(obj_name, ".RData"))
  
  save(list = obj_name, file = save_path)
  
  message("Saved: ", save_path)
}

#-----------------------------------------------------------
# SYNC THE FIGURE CACHE (added 2026-08-11)
#-----------------------------------------------------------
# This script writes its results to output_dir, but the figure scripts that
# consume them -- Location_Compare.Rmd, SDOH_Intervals.Rmd,
# Intervals_ByCountry.Rmd, Supp_ExWAS_Intervals.Rmd -- all read from
# output_dir/Pooled_Output/. Nothing in the codebase ever wrote to that
# directory: it was a hand-populated cache, frozen at 2026-08-10 10:08.
#
# The consequence was that re-running this script had NO effect on Figure 1 or
# Figure 2. The figures kept reading the old snapshot while the fresh results sat
# unread one directory up.
#
# Copy the results into the cache so a re-run actually propagates. Stale files
# from earlier runs are reported rather than silently left in place -- notably
# Child_NA_Quantile_Results.RData, which was the Location+Study model orphaned
# under a nonsense name by the label bug fixed above.
pooled_dir <- file.path(output_dir, "Pooled_Output")
dir.create(pooled_dir, recursive = TRUE, showWarnings = FALSE)

fresh <- list.files(output_dir, pattern = "_Quantile_Results\\.RData$", full.names = TRUE)
ok <- file.copy(fresh, pooled_dir, overwrite = TRUE)
message(sprintf("[cache sync] copied %d/%d result files into %s",
                sum(ok), length(fresh), pooled_dir))

cached <- list.files(pooled_dir, pattern = "_Quantile_Results\\.RData$")
orphans <- setdiff(cached, basename(fresh))
if (length(orphans)) {
  if (is.null(SPECS_TO_RUN))
    warning(sprintf(
      "[cache sync] %d file(s) in Pooled_Output/ were NOT produced by this run and are stale: %s. Delete them before building figures.",
      length(orphans), paste(orphans, collapse = ", ")), call. = FALSE)
}
