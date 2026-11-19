# Models, simulations and power (step 3)

Render `.Rmd` files and run `.R` scripts from this folder unless the header
says otherwise. Inputs are the harmonized datasets in
`Input/Harmonized_Datasets/`; outputs go to `Output/Paper1_Output_Submission2/`.

| Script | What it does |
|---|---|
| `harmonization_helpers.R` | shared helpers (study names, unit conversion) |
| `censored_quantile_helpers.R` | censored quantile regression for values below the limit of detection |
| `Location_QuantileReg_Fixed.R` | pooled between-cohort models (location, age, sex, …); select specifications with `HHEAR_SPECS` |
| `Age_Sex_Match.Rmd` | age- and sex-matched models |
| `Within_Cohort_SDOH_QuantileReg.Rmd` | within-cohort sociodemographic and sampling-year models |
| `TEDDY_Analysis.Rmd`, `ZIP_Analysis.Rmd` | multi-site cohorts: site, country and state models |
| `Exposure_Burden_Creation_Perc_S2.Rmd` → `ByCohort_BurdenExWAS_S2.Rmd` | cumulative exposure score and its within-cohort models |
| `Correlations_EffectiveVars_S2.Rmd` | within-class correlation and effective number of exposures |
| `Targeted_ExWAS_Replication.Rmd` | exposure–outcome associations within and across cohorts (Fig. 3) |
| `cohort_regime_diagnostic.R` | within- vs between-cohort decomposition of each association |
| `cohort_pooling_simulation.R`, `budget_allocation_simulation.R` | pooling and sampling-budget simulations (Fig. 3, fig. S15) |
| `Figure3_Results_Numbers.Rmd` | every number quoted for Fig. 3 |
| `Targeted_Quantile_Power_Array.R`, `submit_targeted_quantile_power_array.sh`, `Collect_Power_Results_Submission2.R` | within-cohort power (SLURM array; fig. S14) |
| `BetweenStudy_Quantile_Power_Array.R`, `submit_betweenstudy_quantile_power_array.sh`, `Collect_Between_Study_Power_Results_Submission2.R` | between-cohort power (SLURM array; fig. S13) |

The SLURM submitters read the project location from `HHEAR_ROOT` (default
`$HOME/HHEAR/HHEAR_Quantreg`) and load site-specific modules; adjust the
`#SBATCH` and `module load` lines for your cluster.

Subfolders: [`Reviewer_Round2/`](Reviewer_Round2/) (revision analyses) and
[`Vis_Code/`](Vis_Code/) (figures).
