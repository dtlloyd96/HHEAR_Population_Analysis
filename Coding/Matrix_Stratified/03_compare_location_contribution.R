## ============================================================================
## 03_compare_location_contribution.R
##
## The question this pipeline exists to answer:
##   the published analysis reports location explaining ~49% of trace-element
##   variation. Matrix is perfectly confounded with location (zero locations
##   contribute more than one matrix) and the same analyte differs by a median
##   56x between blood and urine. How much of that 49% survives once location
##   is compared only WITHIN a matrix?
##
## Reads the published Pooled_Output and the two stratified runs. Writes only
## under Output/Matrix_Stratified/.
## ============================================================================
suppressPackageStartupMessages({library(tidyverse); library(data.table)})

ROOT <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
OUT  <- file.path(ROOT, "Output", "Matrix_Stratified")
CW   <- fread(file.path(ROOT, "Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv"))

## Same weighted-pseudo-R2 aggregation the published figures use
## (Age_Sex_Match_Compare.Rmd / Location_Compare.Rmd): weight each tau by a
## gaussian centred on the median, average within analyte, then average within
## category.
collect <- function(z, acc = list()) {
  if (is.data.frame(z)) {
    if (all(c("pseudo_R","Tau") %in% names(z))) acc[[length(acc)+1]] <- as.data.table(z)
    return(acc)
  }
  if (is.list(z)) for (el in z) acc <- collect(el, acc)
  acc
}

category_r2 <- function(file) {
  if (!file.exists(file)) return(NULL)
  e <- new.env(); load(file, envir = e); L <- get(ls(e)[1], envir = e)
  per <- rbindlist(lapply(names(L), function(nm) {
    fr <- collect(L[[nm]]); if (!length(fr)) return(NULL)
    d <- rbindlist(fr, fill = TRUE)
    d[, .(Analyte_Code = nm, pseudo_R, Tau)]
  }), fill = TRUE)
  if (!nrow(per)) return(NULL)
  per[, plot_R2 := pmax(pseudo_R, 0, na.rm = TRUE)]
  per[, w := exp(-((Tau - 0.5)^2) / (2 * 0.15^2))]
  a <- per[, .(R2 = weighted.mean(plot_R2, w, na.rm = TRUE)), by = Analyte_Code]
  a <- merge(a, unique(CW[Category != "", .(Analyte_Code = Code, Category)]), by = "Analyte_Code")
  a[, .(R2 = mean(R2, na.rm = TRUE), n_analytes = .N), by = Category]
}

sources <- list(
  published = file.path(ROOT, "Output/Paper1_Output_Submission2/Pooled_Output"),
  Urine     = file.path(OUT, "stageA_Urine/Output/Paper1_Output_Submission2"),
  Blood     = file.path(OUT, "stageA_Blood/Output/Paper1_Output_Submission2")
)

res <- rbindlist(lapply(names(sources), function(s) {
  out <- rbindlist(lapply(c("Mat_Location_Quantile_Results.RData",
                            "Child_Location_Quantile_Results.RData"), function(f) {
    r <- category_r2(file.path(sources[[s]], f)); if (is.null(r)) return(NULL)
    r[, panel := sub("_.*$", "", f)][]
  }), fill = TRUE)
  if (!nrow(out)) return(NULL)
  out[, .(R2 = mean(R2), n_analytes = sum(n_analytes)), by = Category][, source := s][]
}), fill = TRUE)

wide <- dcast(res, Category ~ source, value.var = "R2")
for (j in setdiff(names(wide), "Category"))
  set(wide, j = j, value = round(100 * wide[[j]], 1))
if (all(c("published","Urine") %in% names(wide)))
  wide[, urine_minus_published := round(Urine - published, 1)]
setorderv(wide, "published", -1, na.last = TRUE)

cat("\n=== Location contribution (%): published pooled vs within-matrix ===\n")
cat("Published pools blood and urine cohorts in one model. Urine/Blood columns\n")
cat("compare locations only within a matrix, so the medium cannot masquerade\n")
cat("as geography.\n\n")
print(wide)

fwrite(wide, file.path(OUT, "location_contribution_by_matrix.csv"))
fwrite(res,  file.path(OUT, "location_contribution_long.csv"))
cat("\nwritten: Output/Matrix_Stratified/location_contribution_by_matrix.csv\n")
