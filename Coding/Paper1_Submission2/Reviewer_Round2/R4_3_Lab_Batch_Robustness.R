# R4 comment 3 -- are the location estimates robust to analytical lab / batch?
# "If the authors have access to some of this info, it would be good that they
#  show that the geographic and temporal variance components are robust to
#  these factors."
#
# HHEAR records the analysis CENTER (laboratory) per cohort; sample-level batch
# and storage time are not available.
#
# WHY NOT ADJUST FOR CENTER: the Table S8 sensitivity models (Location + Center,
# Location + Year) were checked 2026-11-03 and cannot be used -- center is
# almost perfectly aliased with location, so those fits largely fail (e.g. Mn
# returns no Location + Center result; class means came out BELOW location
# alone, which is impossible for a nested model that fitted).
#
# INSTEAD, STRATIFY: location contribution among cohorts analysed in the SAME
# HHEAR lab, within ONE matrix (R4_3_prep_within_lab.R builds the strata; the
# runs are run_location_unadjusted.sh Location lab_<Center>_<Matrix>). Lab and
# matrix are then constant within each comparison. Compared against the
# all-labs within-matrix estimate (B4) on the SAME analytes.
#
# Writes Output/Reviewer_Round2/R4_3_batch/. Run from the project root.
suppressPackageStartupMessages(library(data.table))
OUT <- "Output/Reviewer_Round2/R4_3_batch"
CW  <- fread("Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv")
CLASSES_11 <- c("FlameRetardants","PAH","Parabens","PCB","Pesticides","PFAS",
                "Phenols","Phthalates","Smoking","TraceElements","VOCs")

collect <- function(z, acc = list()) {
  if (is.data.frame(z)) { if (all(c("pseudo_R","Tau") %in% names(z))) acc[[length(acc)+1]] <- as.data.table(z); return(acc) }
  if (is.list(z)) for (el in z) acc <- collect(el, acc); acc
}
analyte_r2 <- function(dir) rbindlist(lapply(c("Child", "Mat"), function(p) {
  f <- file.path(dir, paste0(p, "_Location_Quantile_Results.RData")); if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e); L <- get(ls(e)[1], envir = e)
  d <- rbindlist(lapply(names(L), function(nm) { fr <- collect(L[[nm]]); if (!length(fr)) return(NULL)
    rbindlist(fr, fill = TRUE)[, .(Analyte_Code = nm, pseudo_R, Tau)] }), fill = TRUE)
  if (!nrow(d)) return(NULL)
  d[, w := exp(-((Tau - 0.5)^2) / (2 * 0.15^2))]
  d[, .(R2 = 100 * weighted.mean(pmax(pseudo_R, 0, na.rm = TRUE), w, na.rm = TRUE)), by = Analyte_Code][, panel := p][]
}))

variants <- readLines(file.path(OUT, "within_lab_variants.txt"))
lab <- rbindlist(lapply(variants, function(v) {
  r <- analyte_r2(file.path("Output/Dilution_Unadjusted", paste0("stageA_", v), "Output/Paper1_Output_Submission2"))
  if (!nrow(r)) return(NULL)
  r[, `:=`(Center = sub("^lab_(.*)_(Urine|Blood)$", "\\1", v), Matrix = sub("^.*_(Urine|Blood)$", "\\1", v))][]
}))
ref <- rbindlist(list(
  analyte_r2("Output/Matrix_Stratified/stageA_Urine/Output/Paper1_Output_Submission2")[, Matrix := "Urine"],
  analyte_r2("Output/Matrix_Stratified/stageA_Blood/Output/Paper1_Output_Submission2")[, Matrix := "Blood"]))
setnames(ref, "R2", "R2_all_labs")
m <- merge(lab, ref, by = c("Analyte_Code", "panel", "Matrix"))
m <- merge(m, unique(CW[Category %in% CLASSES_11, .(Analyte_Code = Code, Category)]), by = "Analyte_Code")
fwrite(m, file.path(OUT, "within_lab_by_analyte.csv"))

cls <- m[, .(analytes = uniqueN(Analyte_Code),
             within_lab = round(mean(R2), 1), all_labs_same_analytes = round(mean(R2_all_labs), 1)),
         by = .(Matrix, Center, Category)][order(Matrix, Center, -all_labs_same_analytes)]
fwrite(cls, file.path(OUT, "within_lab_by_class.csv"))
overall <- m[, .(pairs = .N, within_lab_median = round(median(R2), 1), all_labs_median = round(median(R2_all_labs), 1),
                 spearman = round(cor(R2, R2_all_labs, method = "spearman"), 2)), by = Matrix]
fwrite(overall, file.path(OUT, "within_lab_overall.csv"))
options(width = 200); print(cls); cat("\nOverall (analyte x panel pairs):\n"); print(overall)
