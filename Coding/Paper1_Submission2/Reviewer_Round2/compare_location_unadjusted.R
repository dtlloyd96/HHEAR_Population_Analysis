## ============================================================================
## compare_location_unadjusted.R
##
## Location contribution per exposure class: published run (urine cohorts mixed
## across creatinine-divided, SG-normalized and unadjusted scales) vs the same
## model on unadjusted concentrations (one ng/mL scale for every urine cohort).
##
## Aggregation is identical to Coding/Matrix_Stratified/03_compare_location_
## contribution.R and the published figures: gaussian tau weights (sd 0.15,
## centred on the median), mean within analyte, mean within class, mean of the
## child and adult panels.
##
## Classes measured only in blood/serum are untouched by dilution adjustment, so
## they should not move: that is the internal control.
##
## Writes only under Output/Dilution_Unadjusted/.
## Run from the project root.
## ============================================================================
suppressPackageStartupMessages({library(data.table)})

OUT <- "Output/Dilution_Unadjusted"
CW  <- fread("Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv")

collect <- function(z, acc = list()) {
  if (is.data.frame(z)) {
    if (all(c("pseudo_R", "Tau") %in% names(z))) acc[[length(acc) + 1]] <- as.data.table(z)
    return(acc)
  }
  if (is.list(z)) for (el in z) acc <- collect(el, acc)
  acc
}

analyte_r2 <- function(file) {
  if (!file.exists(file)) return(NULL)
  e <- new.env(); load(file, envir = e); L <- get(ls(e)[1], envir = e)
  per <- rbindlist(lapply(names(L), function(nm) {
    fr <- collect(L[[nm]]); if (!length(fr)) return(NULL)
    rbindlist(fr, fill = TRUE)[, .(Analyte_Code = nm, pseudo_R, Tau)]
  }), fill = TRUE)
  if (!nrow(per)) return(NULL)
  per[, plot_R2 := pmax(pseudo_R, 0, na.rm = TRUE)]
  per[, w := exp(-((Tau - 0.5)^2) / (2 * 0.15^2))]
  a <- per[, .(R2 = weighted.mean(plot_R2, w, na.rm = TRUE)), by = Analyte_Code]
  merge(a, unique(CW[Category != "", .(Analyte_Code = Code, Category)]), by = "Analyte_Code")
}

# adjusted_rebuild = same code, adjusted rebuild of the data, run through the
# same scratch runner; it should reproduce `published` and confirms that the
# published -> unadjusted difference is the dilution adjustment alone.
sources <- c(published        = "Output/Paper1_Output_Submission2/Pooled_Output",
             adjusted_rebuild = file.path(OUT, "stageA_adjusted/Output/Paper1_Output_Submission2"),
             unadjusted       = file.path(OUT, "stageA_unadjusted/Output/Paper1_Output_Submission2"),
             # sensitivity: creatinine-corrected cohorts x100 (ug/g creatinine),
             # SG-normalized and raw cohorts unchanged
             ugg              = file.path(OUT, "stageA_ugg/Output/Paper1_Output_Submission2"))

per_analyte <- rbindlist(lapply(names(sources), function(s) {
  rbindlist(lapply(c(Adult = "Mat_Location_Quantile_Results.RData",
                     Child = "Child_Location_Quantile_Results.RData"), function(f) {
    r <- analyte_r2(file.path(sources[[s]], f)); if (is.null(r)) return(NULL)
    r[, panel := if (startsWith(f, "Mat")) "Adult" else "Child"][]
  }))[, source := s][]
}), fill = TRUE)
fwrite(per_analyte, file.path(OUT, "location_contribution_by_analyte.csv"))

# Which classes are touched at all? (any urine cohort that was CRE/SG adjusted)
e <- new.env(); load(file.path(OUT, "harmonized_unadjusted/Harmonized_Targeted.RData"), envir = e)
U <- as.data.table(e$Targeted_Data_All)
# The adjusted rebuild comes from the same script run, so rows line up by
# position (verified: identical to Input/ in every concentration).
e2 <- new.env(); load(file.path(OUT, "harmonized_adjusted/Harmonized_Targeted.RData"), envir = e2)
A <- as.data.table(e2$Targeted_Data_All)
stopifnot(nrow(A) == nrow(U), identical(A$Analyte_Code, U$Analyte_Code), identical(A$Study, U$Study))
differs <- !((is.na(A$Concentration) & is.na(U$Concentration)) |
             abs(A$Concentration - U$Concentration) <= 1e-9 * pmax(1, abs(U$Concentration)))
touched <- unique(U$Analyte_Code[differs %in% TRUE])
per_analyte[, dilution_affected := Analyte_Code %in% touched]

cat_tab <- per_analyte[, .(R2 = mean(R2, na.rm = TRUE), n_analytes = .N,
                           n_affected = sum(dilution_affected)), by = .(Category, panel, source)][
  , .(R2 = mean(R2), n_analytes = sum(n_analytes), n_affected = sum(n_affected)), by = .(Category, source)]
wide <- dcast(cat_tab, Category ~ source, value.var = c("R2", "n_analytes", "n_affected"))
wide[, `:=`(published_pct = round(100 * R2_published, 1),
            rebuild_pct = round(100 * R2_adjusted_rebuild, 1),
            unadjusted_pct = round(100 * R2_unadjusted, 1),
            ugg_pct = round(100 * R2_ugg, 1))]
wide[, change_pts := unadjusted_pct - published_pct]
out <- wide[, .(Category, published_pct, rebuild_pct, unadjusted_pct, ugg_pct, change_pts,
                analytes = n_analytes_published, dilution_affected_analytes = n_affected_published)][
  order(-published_pct)]
fwrite(out, file.path(OUT, "location_contribution_published_vs_unadjusted.csv"))

cat("\n=== Location contribution (%): published vs unadjusted urine ===\n")
print(out)
cat("\nInternal control -- analytes with no dilution adjustment anywhere should not move:\n")
ctrl <- dcast(per_analyte[dilution_affected == FALSE], Analyte_Code + panel ~ source, value.var = "R2")
print(ctrl[, .(n = .N, max_abs_diff = max(abs(unadjusted - published), na.rm = TRUE))])
