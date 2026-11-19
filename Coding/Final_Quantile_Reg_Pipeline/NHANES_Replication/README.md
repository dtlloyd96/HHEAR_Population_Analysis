# NHANES reference intervals

Nationally representative US intervals, shown alongside the HHEAR cohorts in
Fig. 1B and in the Shiny app. Run the scripts from `NHANES_Replication/Coding/`;
paths below are relative to `NHANES_Replication/`.

1. `download_nhanesa.R` and `download_nhanesa_exposure.R` download the
   demographic and chemical tables with the `nhanesA` package into `Input/`.
2. `NHANES_Quantile_Reg_Survey_Design.Rmd` fits survey-weighted quantile models.
3. `NHANES_Ref_Intervals_Plot.Rmd` writes `Output/NHANES_Subset_Intervals.RData`.

Units: in the merged NHANES extract (`Input/NHANES_ExWAS_Data.Rdata`), the
urinary PAH metabolites NAP1, NAP2 (naphthols) and PYR1 (1-hydroxypyrene) are
in ng/L; every other analyte is in ng/mL.
