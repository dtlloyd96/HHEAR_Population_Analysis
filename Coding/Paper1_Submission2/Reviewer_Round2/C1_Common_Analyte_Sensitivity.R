# C1 (part 2) -- location contribution under the Reviewer 1 sensitivity rules,
# on the canonical ug/g data. Same tau-weighted pseudo-R2 aggregation as the
# published figures (analyte = gaussian-weighted mean over taus; class = mean
# over analytes, then mean of child and adult panels).
#
# Columns
#   all_pooled          every analyte, every cohort, matrices pooled (= main analysis)
#   common_pooled       common analyte set only (>=3 cohorts, >=2 continents)
#   all_urine           every analyte, within urine           (B4)
#   common_urine        common analytes, within urine
#   hiov_pooled         low-overlap cohorts removed, refit, matrices pooled
#   hiov_common_pooled  ... and common analytes only
#   hiov_common_urine   ... within urine, common analytes only  (strictest)
#
# Prerequisites: C1_prep_common_analyte_inputs.R, then
#   run_location_unadjusted.sh Location c1_highoverlap
#   run_location_unadjusted.sh Location c1_highoverlap_urine
# Writes Output/Reviewer_Round2/C1/. Run from the project root.
suppressPackageStartupMessages(library(data.table))
OUT <- "Output/Reviewer_Round2/C1"
CW  <- fread("Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv")
CLASSES_11 <- c("FlameRetardants","PAH","Parabens","PCB","Pesticides","PFAS",
                "Phenols","Phthalates","Smoking","TraceElements","VOCs")
common <- fread(file.path(OUT, "common_analyte_set.csv"))$Analyte_Code

collect <- function(z, acc = list()) {
  if (is.data.frame(z)) { if (all(c("pseudo_R","Tau") %in% names(z))) acc[[length(acc)+1]] <- as.data.table(z); return(acc) }
  if (is.list(z)) for (el in z) acc <- collect(el, acc); acc
}
analyte_r2 <- function(dir) rbindlist(lapply(c(Adult = "Mat", Child = "Child"), function(p) {
  f <- file.path(dir, paste0(p, "_Location_Quantile_Results.RData")); if (!file.exists(f)) return(NULL)
  e <- new.env(); load(f, envir = e); L <- get(ls(e)[1], envir = e)
  d <- rbindlist(lapply(names(L), function(nm) { fr <- collect(L[[nm]]); if (!length(fr)) return(NULL)
    rbindlist(fr, fill = TRUE)[, .(Analyte_Code = nm, pseudo_R, Tau)] }), fill = TRUE)
  d[, w := exp(-((Tau - 0.5)^2) / (2 * 0.15^2))]
  d[, .(R2 = 100 * weighted.mean(pmax(pseudo_R, 0, na.rm = TRUE), w, na.rm = TRUE)), by = Analyte_Code][, panel := p][]
}))

src <- list(
  all      = "Output/Paper1_Output_Submission2/Pooled_Output",
  urine    = "Output/Matrix_Stratified/stageA_Urine/Output/Paper1_Output_Submission2",
  hiov     = "Output/Dilution_Unadjusted/stageA_c1_highoverlap/Output/Paper1_Output_Submission2",
  hiov_ur  = "Output/Dilution_Unadjusted/stageA_c1_highoverlap_urine/Output/Paper1_Output_Submission2")
A <- rbindlist(lapply(names(src), function(s) analyte_r2(src[[s]])[, source := s]))
A <- merge(A, unique(CW[Category %in% CLASSES_11, .(Analyte_Code = Code, Category)]), by = "Analyte_Code")
A[, in_common := Analyte_Code %in% common]
fwrite(A, file.path(OUT, "location_by_analyte_all_rules.csv"))

cls <- function(d) d[, .(R2 = mean(R2)), by = .(Category, panel)][, .(R2 = round(mean(R2), 1)), by = Category]
cols <- list(
  all_pooled         = A[source == "all"],
  common_pooled      = A[source == "all" & in_common],
  all_urine          = A[source == "urine"],
  common_urine       = A[source == "urine" & in_common],
  hiov_pooled        = A[source == "hiov"],
  hiov_common_pooled = A[source == "hiov" & in_common],
  hiov_common_urine  = A[source == "hiov_ur" & in_common])
tab <- Reduce(function(x, y) merge(x, y, by = "Category", all = TRUE),
              lapply(names(cols), function(n) setnames(cls(cols[[n]]), "R2", n)))
n_an <- A[source == "all", .(analytes = uniqueN(Analyte_Code), common_analytes = uniqueN(Analyte_Code[in_common])), by = Category]
tab <- merge(n_an, tab, by = "Category")[order(-all_pooled)]
fwrite(tab, file.path(OUT, "location_by_class_sensitivity.csv"))

options(width = 200)
cat("Location contribution (%) by class under each rule:\n"); print(tab)
cat("\nRank agreement with the main analysis (Spearman, classes present in both):\n")
for (n in setdiff(names(cols), "all_pooled")) {
  x <- tab[!is.na(get(n)) & !is.na(all_pooled)]
  cat(sprintf("  %-20s rho = %5.2f  (%d classes)\n", n, cor(x$all_pooled, x[[n]], method = "spearman"), nrow(x)))
}
