# Matrix-stratified sensitivity analysis (fig. S20)

Biomarkers of the same chemical can be measured in blood in some cohorts and
urine in others, and no location contributes more than one matrix for a given
exposure class, so matrix and location are confounded in the pooled models.
These scripts refit the location model within each matrix and compare the
location contribution with the pooled result. Run from the repository root:

1. `01_prep_matrix_aware_input.R`: splits the harmonized data by matrix into
   `Output/Matrix_Stratified/Input_<matrix>/`.
2. `02_run_stage_a_by_matrix.sh`: runs `Paper1_Submission2/Location_QuantileReg_Fixed.R`
   once per matrix.
3. `03_compare_location_contribution.R`: writes `location_contribution_by_matrix.csv`.
4. `04_supp_figure_matrix.R`: fig. S20.

Everything is written under `Output/Matrix_Stratified/`; nothing else is overwritten.
