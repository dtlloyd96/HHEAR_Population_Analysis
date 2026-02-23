# HHEAR Manuscript Analysis

Repository for analyses supporting the manuscript:
**Population Structure of the Chemical Exposome Across Global Cohorts**.

## License
This repository is released under the **MIT License**.
See `LICENSE`.

## How to Cite
If you use this code or derived outputs, please cite both:
1. The manuscript (when published).
2. This repository.

### Manuscript citation 
Lloyd D.T, et al. *Population Structure of the Chemical Exposome Across Global Cohorts*. Year;Journal:Pages. DOI.

### Repository/software citation 
Lloyd D.T., Patel C.J *Population Structure of the Chemical Exposome Across Global Cohorts* (Version 1) [R Software]. 2026.  
`https://github.com/dtlloyd96/HHEAR_Population_Analysis`

## Folder layout
- `Coding/`: data cleaning, modeling, interval generation, and study-specific analyses.
- `Coding/Vis_Code`: Data visualizations
- `Input/`: Underlying data must be downloaded from [HHEAR Data Website](https://hheardatacenter.mssm.edu/), and fed through cleaning scripts.
- `Output/`: analysis outputs (`.RData`, `.pdf`, `.csv`).

## Suggested run order
1. Run harmonization scripts.
2. Run pooled and by-study quantile regression scripts.
3. Build interval objects and matched/comparison datasets.
4. Run scripts in `Coding/Vis_Code/` to assemble final figures.

## Script annotations (`Coding/`)

### `Coding/HHEAR_Data_Clean_V2.Rmd`
- **Purpose:** Harmonize cohort-level SDOH/demographic fields into one common schema.
- **Main inputs:** Raw cohort `*All.RData` files in `Input/`.
- **Main outputs:** `Input/Harmonized_Datasets/Harmonized_SDOH.RData`, `TEDDY_SDOH_Clean.RData`, `ZIP_SDOH_Clean.RData`.

### `Coding/HHEAR_Data_Clean_Targeted_V2.Rmd`
- **Purpose:** Harmonize targeted analyte data across cohorts; apply LOD imputation and SG/creatinine corrections.
- **Main inputs:** Raw cohort `*All.RData` files in `Input/`.
- **Main outputs:** `Input/Harmonized_Datasets/Harmonized_Targeted.RData`.

### `Coding/Pooled_QuantileReg.R`
- **Purpose:** Primary pooled quantile-regression pipeline for child and maternal analyses.
- **Main inputs:** `Harmonized_SDOH.RData`, `Harmonized_Targeted.RData`.
- **Main outputs:** Multiple `Output/Final_Paper1_Output/*.RData` model result objects.

### `Coding/ExWAS_ByStudy.R`
- **Purpose:** Run study-specific ExWAS/quantile models and export standardized per-study result sets.
- **Main inputs:** `Harmonized_SDOH.RData`, `Harmonized_Targeted.RData`.
- **Main outputs:** `Output/By_Study_*.RData` result objects.

### `Coding/Within_Cohort_SDOH_QuantileReg.Rmd`
- **Purpose:** Within-cohort quantile regressions across multiple child/maternal covariate specifications.
- **Main inputs:** `Harmonized_SDOH.RData`, `Harmonized_Targeted.RData`.
- **Main outputs:** `Output/Paper1_OutputV2/By_Study_*.RData` objects.

### `Coding/Location_QuantileReg_Fixed.R`
- **Purpose:** Location-focused quantile-regression variant with demographic adjustment sets.
- **Main inputs:** Harmonized SDOH/targeted datasets.
- **Main outputs:** Location and age/sex result objects in `Output/Paper1_OutputV2/`.

### `Coding/Location_Sample_Year_Reg.R`
- **Purpose:** Sensitivity workflow adding sample-year structure to location-oriented models.
- **Main inputs:** Harmonized SDOH/targeted datasets.
- **Main outputs:** Child/maternal result objects in `Output/Paper1_OutputV2/`.

### `Coding/Age_Sex_Match.Rmd`
- **Purpose:** Create age/sex-matched child and maternal datasets and rerun quantile models on matched subsets.
- **Main inputs:** Harmonized SDOH/targeted datasets.
- **Main outputs:** `Age_Sex_Child_Match.RData`, `Age_Sex_Mat_Match.RData`.

### `Coding/Correlatons_EffectiveVars.Rmd`
- **Purpose:** Compute correlation structure among exposures and effective-variable summaries; create correlation plots.
- **Main inputs:** Harmonized SDOH/targeted datasets.
- **Main outputs:** Correlation figures and derived correlation tables in `Output/Final_Paper1_Output/Plots/`.

### `Coding/Intervals_ByCountry.Rmd`
- **Purpose:** Build country/location-weighted reference interval objects from pooled quantile model outputs.
- **Main inputs:** Quantile results + harmonized datasets.
- **Main outputs:** Country-specific interval objects and plotting-ready analyte lists.

### `Coding/SDOH_Intervals.Rmd`
- **Purpose:** Build SDOH-stratified interval summaries and interval-linked plotting data.
- **Main inputs:** Quantile results + harmonized datasets.
- **Main outputs:** SDOH interval objects used by downstream visualization scripts.

### `Coding/TEDDY_Analysis.Rmd`
- **Purpose:** TEDDY-only quantile-regression workflow (single-exposure and multivariable models).
- **Main inputs:** `TEDDY_SDOH_Clean.RData`, `Harmonized_Targeted.RData`.
- **Main outputs:** `Output/Final_Paper1_Output/TEDDY_*.RData`.

### `Coding/ZIP_Analysis.Rmd`
- **Purpose:** ZIP-only quantile-regression workflow with race/education/country/state/site model variants.
- **Main inputs:** `ZIP_SDOH_Clean.RData`, `Harmonized_Targeted.RData`.
- **Main outputs:** `Output/Final_Paper1_Output/ZIP_*.RData`.

## Script annotations (`Coding/Vis_Code/`)

### Figure component builders
- `Coding/Vis_Code/Cumulative_Conc_Plot_Map.Rmd`: cumulative concentration summaries + map assets.
- `Coding/Vis_Code/Sociodemo_PCA.Rmd`: within-study sociodemographic PCA distributions (`summ_socio_ridge.RData`).
- `Coding/Vis_Code/Geography_PCA.Rmd`: geography PCA/ranking summaries (`Avg_World_Map.RData` + ridge/rank plots).
- `Coding/Vis_Code/Within_Study_SDOH_Compare.Rmd`: by-study R2 comparisons and year-linked components.
- `Coding/Vis_Code/Location_Compare.Rmd`: pooled location/demographic explanatory-power barplots.
- `Coding/Vis_Code/Age_Sex_Match_Compare.Rmd`: matched-model explanatory-power comparisons.
- `Coding/Vis_Code/Mother_Child_Compare.Rmd`: maternal vs child comparison visuals.
- `Coding/Vis_Code/Figure2DemoPlot.Rmd`: demographic patch components.

### Interval visualization
- `Coding/Vis_Code/Intervals_ByCountry.Rmd`: visualization-oriented country interval workflow.
- `Coding/Vis_Code/Matched_Intervals.Rmd`: matched-cohort interval visuals.
- `Coding/Vis_Code/SDOH_Intervals.Rmd`: SDOH interval visual workflow.

### Final figure assembly
- `Coding/Vis_Code/Time_Patchplot.Rmd`: time-effect patch plot.
- `Coding/Vis_Code/Age_Match_Patch.Rmd`: age-match patch figure.
- `Coding/Vis_Code/Geography_Patchplot.Rmd`: geography multi-panel figure assembly.
- `Coding/Vis_Code/SocioDemo_Patchplot.Rmd`: sociodemographic multi-panel assembly.
- `Coding/Vis_Code/All_Patchplots.Rmd`: master final manuscript patchwork assembly.

## Notes
- Most scripts assume relative execution paths from their own directory.
- Many scripts write intermediate `.RData` files consumed by downstream scripts.

