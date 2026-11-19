## B4 -- location contribution pooled across matrices vs within a single matrix,
## on the ug/g-creatinine canonical data (Round 2: Reviewer 4, comment 3).
## Pooled = Output/Dilution_Unadjusted/stageA_ugg (same data as the canonical
## build); within-matrix = Output/Matrix_Stratified/stageA_{Urine,Blood}.
## Same tau-weighted pseudo-R2 aggregation as the published figures.
## Writes Output/Reviewer_Round2/B4/. Run from the project root.
suppressPackageStartupMessages(library(data.table))
OUT <- "Output/Reviewer_Round2/B4"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CW  <- fread("Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv")
CLASSES_11 <- c("FlameRetardants","PAH","Parabens","PCB","Pesticides","PFAS",
                "Phenols","Phthalates","Smoking","TraceElements","VOCs")
collect <- function(z, acc = list()) {
  if (is.data.frame(z)) { if (all(c("pseudo_R","Tau") %in% names(z))) acc[[length(acc)+1]] <- as.data.table(z); return(acc) }
  if (is.list(z)) for (el in z) acc <- collect(el, acc); acc
}
analyte_r2 <- function(file) {
  if (!file.exists(file)) return(NULL)
  e <- new.env(); load(file, envir = e); L <- get(ls(e)[1], envir = e)
  per <- rbindlist(lapply(names(L), function(nm) {
    fr <- collect(L[[nm]]); if (!length(fr)) return(NULL)
    rbindlist(fr, fill = TRUE)[, .(Analyte_Code = nm, pseudo_R, Tau)] }), fill = TRUE)
  if (!nrow(per)) return(NULL)
  per[, w := exp(-((Tau - 0.5)^2) / (2 * 0.15^2))]
  per[, .(R2 = 100 * weighted.mean(pmax(pseudo_R, 0, na.rm = TRUE), w, na.rm = TRUE)), by = Analyte_Code]
}
src <- c(pooled = "Output/Dilution_Unadjusted/stageA_ugg/Output/Paper1_Output_Submission2",
         Urine  = "Output/Matrix_Stratified/stageA_Urine/Output/Paper1_Output_Submission2",
         Blood  = "Output/Matrix_Stratified/stageA_Blood/Output/Paper1_Output_Submission2")
pa <- rbindlist(lapply(names(src), function(s) rbindlist(lapply(c("Mat","Child"), function(p) {
  r <- analyte_r2(file.path(src[[s]], paste0(p, "_Location_Quantile_Results.RData")))
  if (!is.null(r)) r[, `:=`(panel = p, source = s)][] }))))
pa <- merge(pa, unique(CW[Category %in% CLASSES_11, .(Analyte_Code = Code, Category)]), by = "Analyte_Code")
an <- pa[, .(R2 = mean(R2)), by = .(Category, Analyte_Code, source)]
aw <- dcast(an, Category + Analyte_Code ~ source, value.var = "R2")
fwrite(aw[order(Category, -pooled)], file.path(OUT, "location_by_analyte_pooled_vs_matrix.csv"))
cls <- pa[, .(R2 = mean(R2), n = .N), by = .(Category, panel, source)][, .(R2 = round(mean(R2), 1)), by = .(Category, source)]
cw <- dcast(cls, Category ~ source, value.var = "R2")[order(-pooled)]
fwrite(cw, file.path(OUT, "location_by_class_pooled_vs_matrix.csv"))
print(cw)
cat("\nTrace elements, per analyte (pooled vs within urine / blood):\n")
print(aw[Category == "TraceElements"][order(-pooled), .(Analyte_Code, pooled = round(pooled,1), Urine = round(Urine,1), Blood = round(Blood,1))])
