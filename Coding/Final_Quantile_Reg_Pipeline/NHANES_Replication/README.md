# NHANES reference intervals

Population-representative US intervals overlaid in Fig. 1B and the Shiny app.
Run from `Coding/`:

1. `download_nhanesa.R` and `download_nhanesa_exposure.R`: demographics and
   chemical tables via the `nhanesA` package, saved to `../Input/`.
2. `NHANES_Quantile_Reg_Survey_Design.Rmd`: survey-weighted quantile models.
3. `NHANES_Ref_Intervals_Plot.Rmd`: writes `Output/NHANES_Subset_Intervals.RData`.

The merged NHANES extract (`Input/NHANES_ExWAS_Data.Rdata`) holds urinary
2-naphthol, 1-naphthol and 1-hydroxypyrene (NAP1, NAP2, PYR1) in ng/L; every
other analyte is in ng/mL.
