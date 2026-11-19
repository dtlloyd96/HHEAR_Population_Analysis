## Refreshes Health_Outcome_Analysis/Input/Harmonized_Targeted.csv (the ExWAS /
## Figure 3 input) from the canonical Input/Harmonized_Datasets/
## Harmonized_Targeted.RData, WITHOUT changing its schema.
##
## The CSV is not a plain dump of the RData: it keeps the raw Sample_Matrix
## strings ("maternal serum") that Figure3_Results_Numbers.Rmd relies on, and has
## no Matrix_Class. Verified 2026-11-03 that its rows are the RData's rows in the
## same order with identical concentrations, so this script swaps in
## Concentration, LOD and Detected by position (after asserting the order still
## matches) and appends Dilution_Adjustment. Everything else is left untouched.
##
## Run from the project root:  Rscript Coding/sync_health_outcome_targeted_csv.R
suppressPackageStartupMessages(library(data.table))

CSV <- "Health_Outcome_Analysis/Input/Harmonized_Targeted.csv"
e <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = e)
R <- as.data.table(e$Targeted_Data_All)
x <- fread(CSV, colClasses = list(character = c("Child_PID", "Mat_PID", "Timepoint", "SID", "trimester")))

stopifnot(nrow(x) == nrow(R),
          identical(x$Study, R$Study),
          identical(x$Analyte_Code, R$Analyte_Code),
          identical(fcoalesce(x$Child_PID, x$Mat_PID), fcoalesce(as.character(R$Child_PID), as.character(R$Mat_PID))))

bak <- sub("\\.csv$", sprintf("_backup_%s.csv", format(Sys.Date())), CSV)
if (!file.exists(bak)) file.copy(CSV, bak)

x[, `:=`(Concentration = R$Concentration, LOD = R$LOD, Detected = R$Detected,
         Dilution_Adjustment = R$Dilution_Adjustment)]
fwrite(x, CSV)
cat("Wrote", CSV, "(backup:", bak, ")\n")
print(x[Analyte_Code == "Mo", .(med = signif(median(Concentration, na.rm = TRUE), 3)), by = Study][order(med)])
