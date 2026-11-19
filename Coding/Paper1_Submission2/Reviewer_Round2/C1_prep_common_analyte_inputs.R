# C1 (part 1) -- define the common analyte set and the high-overlap cohort set,
# and write filtered harmonized inputs for the refits.
# (Round 2: Reviewer 1 -- "repeat key analyses using only the largest common
#  analyte set. Sensitivity analyses excluding low-overlap cohorts ...")
#
# Common analyte set: measured in >= MIN_COH cohorts spanning >= MIN_CONT
#   continents (from B3's analyte_breadth.csv, i.e. the same exposure rule and
#   site-level geography as the paper).
# Low-overlap cohort: contributes <= MAX_LOW of the common analytes.
#
# Writes (canonical ug/g data, filtered):
#   Output/Dilution_Unadjusted/harmonized_c1_highoverlap/Harmonized_Targeted.RData
#   Output/Dilution_Unadjusted/harmonized_c1_highoverlap_urine/Harmonized_Targeted.RData
# so run_location_unadjusted.sh can fit them (VARIANT = c1_highoverlap[_urine]),
# plus the set definitions under Output/Reviewer_Round2/C1/.
# Run from the project root.
suppressPackageStartupMessages(library(data.table))
source("Coding/Paper1_Submission2/harmonization_helpers.R")

MIN_COH <- 3; MIN_CONT <- 2; MAX_LOW <- 3
OUT <- "Output/Reviewer_Round2/C1"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

br <- fread("Output/Reviewer_Round2/B3/analyte_breadth.csv")
common <- br[n_cohorts >= MIN_COH & n_continents >= MIN_CONT]
fwrite(common[order(Category, -n_cohorts)], file.path(OUT, "common_analyte_set.csv"))

load("Input/Harmonized_Datasets/Harmonized_Targeted.RData")
H <- as.data.table(Targeted_Data_All)
H[, Cohort := std_study(Study)]
coh <- H[Analyte_Code %in% common$Analyte_Code, .(common_analytes = uniqueN(Analyte_Code)), by = Cohort]
all_coh <- data.table(Cohort = unique(H$Cohort))
coh <- merge(all_coh, coh, by = "Cohort", all.x = TRUE)[is.na(common_analytes), common_analytes := 0L]
coh[, keep := common_analytes > MAX_LOW]
fwrite(coh[order(-common_analytes)], file.path(OUT, "cohort_overlap.csv"))

write_variant <- function(D, name) {
  d <- file.path("Output/Dilution_Unadjusted", paste0("harmonized_", name))
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  Targeted_Data_All <- D
  save(Targeted_Data_All, file = file.path(d, "Harmonized_Targeted.RData"))
  cat(sprintf("%-28s %8d rows  %2d cohorts  %3d analytes\n", name, nrow(D), uniqueN(std_study(D$Study)), uniqueN(D$Analyte_Code)))
}
keep <- H[Cohort %in% coh[keep == TRUE, Cohort]][, Cohort := NULL]
write_variant(keep, "c1_highoverlap")
write_variant(keep[Matrix_Class == "Urine" | Analyte_Code %in% c("CRE", "SG")], "c1_highoverlap_urine")

cat("\nCommon analyte set:", nrow(common), "analytes; by class:\n"); print(common[, .N, by = Category][order(-N)])
cat("\nCohorts dropped as low-overlap (<=", MAX_LOW, "common analytes):\n"); print(coh[keep == FALSE][order(-common_analytes)])
