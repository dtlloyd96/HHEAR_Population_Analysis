# Models, simulations and power (step 3)

Run these after harmonization (step 2). Render `.Rmd` files and run `.R`
scripts from this folder unless a script's header says otherwise. Inputs are
the harmonized datasets in `Input/Harmonized_Datasets/`; results go to
`Output/Paper1_Output_Submission2/`.

| Script | What it does |
|---|---|
| `harmonization_helpers.R` | Shared helpers: study names, unit conversion. Sourced by most scripts. |
| `censored_quantile_helpers.R` | Censored quantile regression, which handles values below the limit of detection. |
| `Location_QuantileReg_Fixed.R` | Pooled models across cohorts (location, age, sex, …). `HHEAR_SPECS` selects which models to fit. |
| `Age_Sex_Match.Rmd` | The same models on age- and sex-matched subsets. |
| `Within_Cohort_SDOH_QuantileReg.Rmd` | Models within each cohort: sociodemographic factors and sampling year. |
| `TEDDY_Analysis.Rmd`, `ZIP_Analysis.Rmd` | Site, country and state models for the two multi-site cohorts. |
| `Exposure_Burden_Creation_Perc_S2.Rmd` → `ByCohort_BurdenExWAS_S2.Rmd` | Cumulative exposure score, then models of the score within each cohort. |
| `Correlations_EffectiveVars_S2.Rmd` | Correlation among chemicals in the same class and the effective number of independent exposures. |
| `Targeted_ExWAS_Replication.Rmd` | Exposure–outcome associations within and across cohorts (Fig. 3). |
| `cohort_regime_diagnostic.R` | Splits each association into within- and between-cohort parts. |
| `cohort_pooling_simulation.R`, `budget_allocation_simulation.R` | Pooling and sampling-budget simulations (Fig. 3, fig. S15). |
| `Figure3_Results_Numbers.Rmd` | Every number quoted for Fig. 3. |

## Power simulations (figs. S13 and S14)

Both run as SLURM array jobs. The submitters expect a project folder on the
cluster (`HHEAR_ROOT`, default `$HOME/HHEAR/HHEAR_Quantreg`) containing
`Coding/` and, directly in `Input/`, `Harmonized_Targeted.RData` and
`Harmonized_SDOH.RData`. Edit the `#SBATCH` and `module load` lines for your
cluster.

| Power | Simulation | Submitter | Collector |
|---|---|---|---|
| Within cohort (fig. S14) | `Targeted_Quantile_Power_Array.R` | `submit_targeted_quantile_power_array.sh` | `Collect_Power_Results_Submission2.R` |
| Between cohorts (fig. S13) | `BetweenStudy_Quantile_Power_Array.R` | `submit_betweenstudy_quantile_power_array.sh` | `Collect_Between_Study_Power_Results_Submission2.R` |

Submit with `sbatch <submitter>`, then run the collector once every task has
finished.

Subfolders: [`Reviewer_Round2/`](Reviewer_Round2/) (revision analyses) and
[`Vis_Code/`](Vis_Code/) (figures).
