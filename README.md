HHEAR Population Analysis Repository for the Manuscript : Cross-population comparisons of the chemical exposome in participants of the Human Health Exposure Analysis Resource: Recommendations for assembling future exposomic cohorts 

HHEAR Exposomics Analysis Scripts

This repository contains the full set of R Markdown (Rmd) and R scripts used to clean, harmonize, and analyze targeted exposomics data across multiple HHEAR cohorts. Analyses include data cleaning, PCA, exposome-wide association studies (ExWAS), quantile regression, burden scores, demographic comparisons, and study-specific figure generation.

Below is a description of each file and how it fits into the analysis pipeline.

Data Cleaning & Preparation
HHEAR_Data_Clean_V2.Rmd

Cleans, harmonizes, and merges all core HHEAR datasets (demographics, exposures, SDOH) 

HHEAR_Data_Clean_Targeted_V2.Rmd

Version focused specifically on the targeted chemical exposome. Handles  filtering, LOQ handling, and dataset creation.

Exploratory & Descriptive Analysis
Sociodemo_PCA.Rmd

Runs PCA on sociodemographic  variables across cohorts. 

Geography_PCA.Rmd

Conducts PCA on geographic variables


Correlatons_EffectiveVars.Rmd

Computes correlation matrices for exposures and evaluates effective number of variables (Meff) using eigenvalue variance.

Study-Level Comparison and ExWAS
ExWAS_ByStudy.R

Runs study-specific exposome-wide association models using exposures as predictors and selected outcomes, producing standardized result tables.

ByStudy_CompareR2.Rmd

Compares R² and pseudo-R² across studies for key models, summarizing heterogeneity and variability in explained variance.

Quantile Regression
ByCohort_QuantileReg.Rmd

Runs quantile regression by cohort for selected exposures/outcomes, generating τ-specific effect estimates and plots.

Pooled_QuantileReg.R

Runs pooled quantile regression analyses across all studies, with cohort fixed effects and harmonized covariates.

Aggregate Exposure Metrics
ByCohort_AggExWAS.Rmd

Computes burden scores (e.g., IDW, aggregate percentile scores, mixture metrics) within each cohort and performs burden-based ExWAS.

Intervals_ByCountry.Rmd

Exposome Intervals By Country

SDOH_Intervals.Rmd

Summarizes demographic exposure Intervals 

TEDDY and ZIP Study Analyses

TEDDY_Analysis.Rmd

Full TEDDY cohort analysis pipeline 

TEDDY_All_Figures.Rmd

Generates all study-specific figures for the TEDDY cohort

ZIP_Analysis.Rmd

Complete ZIP cohort analysis 

ZIP_All_Figures.Rmd

Produces figures for the ZIP cohort



Cross-Study PCA & Visualization
ByStudy_CompareR2.Rmd

Visualizes and compares model fit 


Additional Files
Correlatons_EffectiveVars.Rmd

Calculates correlation matrices and effective number of variables 
  
