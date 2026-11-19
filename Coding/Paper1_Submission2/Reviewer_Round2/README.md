# Revision analyses (Round 2)

Analyses added during peer review. Run them after harmonization (step 2) and
the models (step 3), from the repository root:
`Rscript Coding/Paper1_Submission2/Reviewer_Round2/<script>`. Results go to
`Output/Reviewer_Round2/`.

| Script(s) | Analysis | Used in |
|---|---|---|
| `B3_Analyte_Cohort_Overlap.R` | analyte × cohort overlap; class × continent counts | fig. S3, table S9 |
| `B4_Matrix_Location_ugg.R` | pooled vs within-matrix location R² (creatinine-adjusted urine) | text |
| `C1_prep_common_analyte_inputs.R` → `C1_Common_Analyte_Sensitivity.R` | restricting to analytes measured in many cohorts | text |
| `C2_Chemical_Year_Trends.R`, `C2b_Cross_Country_Year_Trends.R` → `C2_S8_assemble.R` | within-cohort and cross-country time trends | fig. S12 |
| `C5_Reference_Changes_Inference.R`, `R1_2_PFOS_LDL_Figure.R` → `R1_2_S17_assemble.R` | choice of reference population; PFOS–LDL example | fig. S21 |
| `Fig3_Estimator_Comparison.R` → `Fig3_Estimator_Forest.R` | naive vs fixed-effects vs meta-analytic estimators for the Fig. 3 pairs | fig. S22 |
| `R4_3_prep_within_lab.R` → `R4_3_Lab_Batch_Robustness.R` | single-laboratory subsets | text |
| `location_range_by_analyte.R` | spread of location contributions by analyte | text |
| `run_location_unadjusted.sh` → `compare_location_unadjusted.R` | location models on unadjusted, adjusted and urine-only builds | text; the urine run feeds the Pb overlay in Fig. 1B |
| `R4_5_*.R` | geographic determinants (below) | figs. S8–S10, tables S10–S12, data S1 |
| `R4_5_supp_tables.py` | writes tables S10–S12 and data S1 | |
| `render_chain.sh` | re-renders named chains of model and figure scripts (`fig3`, `matched`, `direct`, `pooled`, `burden`, `corr`, `patch`) | step 5 |

## Geographic determinants (`R4_5_*`)

Run in this order: `geospatial_extract` → `geospatial_no2` →
`geospatial_scatter` → `PAH_PM25` → `PAH_cleanfuel` → `water` →
`diet_pesticide_lead` → `ZIP_pesticide_use` → `source_inventory` →
`within_cohort_sources` → `source_meta_analysis` → `fingerprints` →
`composition_context` → `MEBCO_creatinine_sensitivity` → `B4_PYR_sensitivity`
→ `analysis_register` → `map` → `supp_figures`.

These scripts read public data layers from a folder outside the repository,
`$HHEAR_GEO_CACHE` (default `~/HHEAR_geo_cache`). Download each layer from its
provider into the subfolder shown; full references are in the paper.

| Folder | Layer | Source |
|---|---|---|
| `pm25/V6GL0204_PM25_0p10_<year>.nc` | annual PM2.5, 0.1°, ACAG V6.GL.02.04 | Atmospheric Composition Analysis Group, Washington University (van Donkelaar *et al.* 2021; Shen *et al.* 2024) |
| `no2/*_final_1km.nc` | annual surface NO2, 1 km | Anenberg *et al.* 2022 |
| `ghsl/` | GHS-POP R2023A, 30 arcsec, 1990–2020 epochs (population weights) | European Commission JRC, doi:10.2905/2FF68A52-5B5B-4A22-8F40-C41DA8332CFE |
| `as_gw/` | probability of groundwater As > 10 µg/L, and `data/data/groundwater_utilization_by_country.csv` | Podgorski & Berg 2020 |
| `cleanfuel/EG.CFT.ACCS.ZS_*.json` | access to clean cooking fuels and technologies (SDG 7.1.2) | WHO / World Bank indicator `EG.CFT.ACCS.ZS` |
| `ucmr3/UCMR3_All.txt` | PFAS in US public water systems | US EPA, Third Unregulated Contaminant Monitoring Rule |
| `faostat/` | Pesticides Use, Land Use, Food Balances (current and historic) bulk CSVs | FAO, FAOSTAT |
| `lead/lead-petrol-ban.csv` | year leaded petrol was banned, by country | Our World in Data |
