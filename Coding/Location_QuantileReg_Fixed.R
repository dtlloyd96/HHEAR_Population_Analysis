library(tidyverse)
library(quantreg)
#setwd("/Users/dillon/Library/CloudStorage/GoogleDrive-2015.dillon@gmail.com/My Drive/Chirag_Patel/HHEAR_Analysis/Coding/Paper1_AnalysisV2")
setwd("/home/dil402/HHEAR/HHEAR_Quantreg/Coding")

#load("../../Input/Harmonized_Datasets/Harmonized_SDOH.RData")
#load("../../Input/Harmonized_Datasets/Harmonized_Targeted.RData")


load("../Input/Harmonized_SDOH.RData")
load("../Input/Harmonized_Targeted.RData")

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

Child_Sample_Yr <- Targeted_Data_All %>%
  filter(!is.na(Child_PID)) %>% 
  select(-Mat_PID) %>%                                   
  group_by(Child_PID,Sample_Collection_Year) %>% 
  summarise(n = n())

Mat_Sample_Yr <- Targeted_Data_All %>% 
  filter(!is.na(Mat_PID)) %>% 
  select(-Child_PID) %>%                                   
  group_by(Mat_PID,Sample_Collection_Year) %>% 
  summarise(n = n())
# Split SDOH into child and maternal panels -------------------------------
Child_SDOH <- ExWAS_SDOH_Exposures %>%
  filter(!is.na(Child_PID)) %>%
  select(Study, Child_PID, Mat_PID, dplyr::everything()) %>% 
  left_join(.,Child_Sample_Yr, by = c("Child_PID")) %>% 
  dplyr::select(-n) %>% 
  mutate(Sample_Collection_Year = as.numeric(Sample_Collection_Year))

Mat_SDOH <- ExWAS_SDOH_Exposures %>%
  filter(!is.na(Mat_PID)) %>%                      # keep likely maternal rows
  select(Study, Child_PID, Mat_PID, dplyr::everything()) %>% 
  left_join(.,Mat_Sample_Yr, by = c("Mat_PID")) %>% 
  dplyr::select(-n)%>% 
  mutate(Sample_Collection_Year = as.numeric(Sample_Collection_Year))


Quantile_Reg <- function(exp_var = NULL,
                         dataset,
                         my_outcome,
                         Reg_Type,
                         my_covariates,
                         join_var) {
  
  
  
  if(is.null(exp_var)){
    
    
    mydf_ua <- dataset %>%
      dplyr::select(all_of(c(join_var,my_outcome,my_covariates)))
    
    mydf_ua <- mydf_ua %>%
      ungroup() %>%
      dplyr::filter(complete.cases(.))
    
    if ("Study" %in% colnames(mydf_ua)) {
      
      
      
      table_n <- table(mydf_ua$Study)
      max_category <- names(table_n)[which.max(table_n)]
      
      mydf_ua <- mydf_ua %>%
        dplyr::mutate(Study = factor(Study)) %>%
        dplyr::mutate(Study = relevel(as.factor(Study), ref = max_category))
    }
    
    
    #End Loop if no data
    if(nrow(mydf_ua) == 0){
      next
    }
    
    myform <- formula(paste(my_outcome, "~",  paste0(my_covariates,collapse = "+")))
    
    #Set quantiles
    quantiles <- c(0.05,0.25, 0.5, 0.75,0.95)
    
    if ("Study" %in% colnames(mydf_ua)) {
      
      tocheck <- my_covariates[!my_covariates %in% "Study"]
      clean_df <- mydf_ua %>%
        group_by(Study) %>%
        filter(if_all(all_of(tocheck), ~ n_distinct(.) > 1)) %>%
        ungroup()
    }else{
      clean_df <- mydf_ua
    }
    
    
    #Run Quantreg
    rq_model <- rq(myform, data = clean_df, tau = quantiles)
    
    
    myform_null <- formula(paste(my_outcome, "~",  1))
    
    fit_null <- rq(myform_null, tau = quantiles, data = clean_df)
    yhat_null <- predict(fit_null)
    
    # Predict quantile values for each observation
    predicted_quantiles <- predict(rq_model, newdata = clean_df)
    predicted_quantiles <- predicted_quantiles[complete.cases(predicted_quantiles), ]
    
    yhat_full  <- predict(rq_model)        # fitted at each tau
    yhat_null  <- predict(fit_null)        # intercept-only at each tau
    
    y <- clean_df[[my_outcome]]
    
    
    rho_tau <- function(u, tau) u * (tau - (u < 0))
    # r1_all <- list()
    # for(z in seq_along(quantiles)){
    #
    #   numerator <- sum(rho_tau(y - predicted_quantiles, tau = quantiles[z]))
    #   denominator <- sum(rho_tau(y - yhat_null, tau = quantiles[z]))
    #
    #     # Pseudo R1
    #   r1 <- 1 - numerator / denominator
    #   r1_df <- data.frame("Tau" = quantiles[z],pseudo_R = r1)
    #   r1_all[[z]] <- r1_df
    #
    # }
    #
    # r1_all_df <- bind_rows(r1_all)
    
    r1_all_df <- purrr::map2_dfr(
      seq_along(quantiles), quantiles,
      ~{
        tau <- .y
        u_full <- y - yhat_full[, .x]
        u_null <- y - yhat_null[, .x]
        
        num <- sum(rho_tau(u_full, tau), na.rm = TRUE)
        den <- sum(rho_tau(u_null, tau), na.rm = TRUE)
        
        data.frame(Tau = tau, pseudo_R = if (den == 0) NA_real_ else 1 - num/den)
      }
    )
    
    Quantile_Results <- list(mydf_ua,predicted_quantiles)
    
    ctable_all <- try(summary(rq_model, se = "boot", R = 2,conf.int = TRUE))  # bootstrapped SEs
    
    
    #Create Regression Results
    ctable_final <- data.frame()
    for(z in seq_along(ctable_all)){
      
      curr_ctable <- ctable_all[[z]]
      ctable_df <- as.data.frame(coef(curr_ctable)) %>%
        rownames_to_column("Variable") %>%
        mutate(Tau = curr_ctable$tau,
               Outcome = my_outcome)
      
      names(ctable_df) <- c("Variable","Beta","Standard_Error","T_Value","P_Value","Tau","Outcome")
      ctable_final <- bind_rows(ctable_final,ctable_df)
    }
    conf_int_table <- broom::tidy(rq_model,conf.int = T) %>%
      dplyr::select(Lower = conf.low,Upper = conf.high)
    
    ctable_final <- bind_cols(ctable_final,conf_int_table)
    
    ctable_final <- left_join(ctable_final,r1_all_df, by ="Tau")
    
    
    Final_Results <- list(ctable_final,Quantile_Results)
    
    
    names(Final_Results) <- c("Regression_Results","Quantile_Results")
    
  }else{
    #Select Variables for model
    mydf_ua <- dataset %>%
      dplyr::select(all_of(c(join_var,my_outcome, exp_var,my_covariates))) %>%
      dplyr::filter(!is.na(!!sym(exp_var)))
    
    #End Loop if no data
    if(nrow(mydf_ua) == 0){
      next
    }
    
    
    table_n <- table(mydf_ua$Study)
    max_category <- names(table_n)[which.max(table_n)]
    
    mydf_ua <- mydf_ua %>%
      dplyr::mutate(Study = factor(Study)) %>%
      dplyr::mutate(Study = relevel(as.factor(Study),ref = max_category))
    
    #Remove low count
    if(length(unique(mydf_ua[[exp_var]])) > 10){
      
      mydf_ua <- mydf_ua[complete.cases(mydf_ua), ] %>%
        group_by(!!sym(exp_var)) %>%
        #filter(n() >= 5) %>%
        ungroup()
    }else{
      mydf_ua <- mydf_ua[complete.cases(mydf_ua), ] %>%
        group_by(!!sym(exp_var)) %>%
        filter(n() >= 5) %>%
        ungroup()
    }
    
    #End Loop if no data
    if(nrow(mydf_ua) == 0){
      next
    }
    
    
    tocheck <- exp_var[!exp_var %in% "Study"]
    clean_df <- mydf_ua %>%
      group_by(Study) %>%
      filter(if_all(all_of(tocheck), ~ n_distinct(.) > 1)) %>%
      ungroup()
    
    clean_df <- mydf_ua %>%
      group_by(Study) %>%
      filter(n_distinct(!!sym(exp_var)) > 1) %>%
      ungroup() %>%
      droplevels()
    
    
    
    n_study <- length(unique(clean_df$Study))
    
    #Adjust for study if more than one study has the information
    if(n_study == 1){
      
      myform <- formula(paste(my_outcome, "~",  exp_var))
    }else{
      
      myform <- formula(paste(my_outcome, "~",  exp_var, "+",my_covariates))
      
    }
    
    #Set quantiles
    quantiles <- c(0.05,0.25, 0.5, 0.75,0.95)
    
    
    
    
    #Run Quantreg
    rq_model <- rq(myform, data = clean_df, tau = quantiles)
    
    
    myform_null <- formula(paste(my_outcome, "~",  1))
    
    fit_null <- rq(myform_null, tau = quantiles, data = clean_df)
    yhat_null <- predict(fit_null)
    
    # Predict quantile values for each observation
    predicted_quantiles <- predict(rq_model, newdata = clean_df)
    predicted_quantiles <- predicted_quantiles[complete.cases(predicted_quantiles), ]
    
    yhat_full  <- predict(rq_model)        # fitted at each tau
    yhat_null  <- predict(fit_null)        # intercept-only at each tau
    
    y <- mydf_ua[[my_outcome]]
    
    
    rho_tau <- function(u, tau) u * (tau - (u < 0))
    # r1_all <- list()
    # for(z in seq_along(quantiles)){
    #
    #   numerator <- sum(rho_tau(y - predicted_quantiles, tau = quantiles[z]))
    #   denominator <- sum(rho_tau(y - yhat_null, tau = quantiles[z]))
    #
    #     # Pseudo R1
    #   r1 <- 1 - numerator / denominator
    #   r1_df <- data.frame("Tau" = quantiles[z],pseudo_R = r1)
    #   r1_all[[z]] <- r1_df
    #
    # }
    #
    # r1_all_df <- bind_rows(r1_all)
    
    r1_all_df <- purrr::map2_dfr(
      seq_along(quantiles), quantiles,
      ~{
        tau <- .y
        u_full <- y - yhat_full[, .x]
        u_null <- y - yhat_null[, .x]
        
        num <- sum(rho_tau(u_full, tau), na.rm = TRUE)
        den <- sum(rho_tau(u_null, tau), na.rm = TRUE)
        
        data.frame(Tau = tau, pseudo_R = if (den == 0) NA_real_ else 1 - num/den)
      }
    )
    
    
    Quantile_Results <- list(mydf_ua,predicted_quantiles)
    
    ctable_all <- try(summary(rq_model, se = "boot", R = 2,conf.int = TRUE))  # bootstrapped SEs
    
    
    
    #Create Regression Results
    ctable_final <- data.frame()
    for(z in seq_along(ctable_all)){
      
      curr_ctable <- ctable_all[[z]]
      ctable_df <- as.data.frame(coef(curr_ctable)) %>%
        rownames_to_column("Variable") %>%
        mutate(Tau = curr_ctable$tau,
               Outcome = my_outcome,
               Exposure = exp_var)
      
      names(ctable_df) <- c("Variable","Beta","Standard_Error","T_Value","P_Value","Tau","Outcome","Exposure")
      ctable_final <- bind_rows(ctable_final,ctable_df)
    }
    conf_int_table <- broom::tidy(rq_model,conf.int = T) %>%
      dplyr::select(Lower = conf.low,Upper = conf.high)
    
    ctable_final <- bind_cols(ctable_final,conf_int_table)
    
    ctable_final <- left_join(ctable_final,r1_all_df, by ="Tau")
    
    
    Final_Results <- list(ctable_final,Quantile_Results)
    
    
    names(Final_Results) <- c("Regression_Results","Quantile_Results")
    
  }
  return(Final_Results)
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
        dplyr::select(join_var,Study,all_of(pheno)) %>% 
        #dplyr::filter(Study == "CEHC") %>% 
        dplyr::filter(!is.na(pheno))
      
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

# Study_Only <- c("Location")
# 
# Child_Location_Results <- suppress_all_output(quote(ExWAS_Function_Metab(Exposure_Data = Child_Ads,
#                                                                                     Outcome_Data = child_list_Targeted,
#                                                                                     myexposures = NULL,
#                                                                                     myoutcomes = child_analytes,mycovars = Study_Only,join_var = "Child_PID",
#                                                                                     Reg_Type = "Quantile")))
# 
# save(Child_Location_Results, file = "../Output/Paper1_OutputV2/Child_Study_Quantile_Results.RData")

Child_Age_Sex <- suppress_all_output(quote(ExWAS_Function_Metab(Exposure_Data = Child_Ads,
                                                                                    Outcome_Data = child_list_Targeted,
                                                                                    myexposures = c("Child_Age","Child_Sex"),
                                                                                    myoutcomes = child_analytes,mycovars = "Study",join_var = "Child_PID",
                                                                                    Reg_Type = "Quantile")))

save(Child_Age_Sex, file = "../Output/Paper1_OutputV2/Child_Age_Sex.RData")


Mat_Age <- suppress_all_output(quote(ExWAS_Function_Metab(Exposure_Data = Mat_Ads,
                                                                Outcome_Data = mat_list_Targeted,
                                                                myexposures = c("Mat_Age"),
                                                                myoutcomes = mat_analytes,mycovars = "Study",join_var = "Mat_PID",
                                                                Reg_Type = "Quantile")))

save(Mat_Age, file = "../Output/Paper1_OutputV2/Mat_Age.RData")

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
names(my_covars_list) <- c("Study", "Year", "Location", "Location_Region",
                           "Center","Child_Age","Child_Sex","Location_Center",
                           "Location_Year","Location_Study")

#-----------------------------------------------------------
# Loop through each covariate specification
#-----------------------------------------------------------
for (nm in names(my_covars_list)) {
  
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
  save_path <- file.path("../Output/Paper1_OutputV2",
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
                           "Center","Mat_Age","Location_Center",
                           "Location_Year","Location_Study")

for (nm in names(my_covars_list)) {
  
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
  save_path <- file.path("../Output/Paper1_OutputV2",
                         paste0(obj_name, ".RData"))
  
  save(list = obj_name, file = save_path)
  
  message("Saved: ", save_path)
}
