# R4 comment 3 -- location contribution WITHIN a single HHEAR laboratory.
#
# Adjusting for analysis center (Location + Center, Table S8) does not work:
# center is almost perfectly aliased with location, and those fits mostly fail.
# Stratifying instead removes between-laboratory differences by design: compare
# locations only among cohorts analysed by the same HHEAR lab, and within one
# sample matrix (so the matrix confound is also excluded).
#
# Cohorts analysed by more than one center ("Emory:Minnesota", ...) and cohorts
# with no recorded center are excluded, since their lab cannot be assigned.
# Writes Output/Dilution_Unadjusted/harmonized_lab_<Center>_<Matrix>/ for every
# lab x matrix stratum with >= 2 locations, and the list of variants to
# Output/Reviewer_Round2/R4_3_batch/within_lab_variants.txt.
# Run from the project root.
suppressPackageStartupMessages(library(data.table))
source("Coding/Paper1_Submission2/harmonization_helpers.R")
OUT <- "Output/Reviewer_Round2/R4_3_batch"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

meta <- fread("Output/Hero_Numbers/SuppTable_StudyMeta_FINAL.tsv", sep = "\t")
setnames(meta, gsub(" ", "_", names(meta)))
meta[, Center := gsub("\"", "", HHEAR_Center)][, Center := gsub(" ", "_", Center)]
meta <- meta[Center != "" & !grepl(":", Center), .(Cohort = Study, Center)]

load("Input/Harmonized_Datasets/Harmonized_Targeted.RData")
H <- as.data.table(Targeted_Data_All)
H[, Cohort := std_study(Study)]
H <- merge(H, meta, by = "Cohort")

e <- new.env(); load("Input/Harmonized_Datasets/Harmonized_SDOH.RData", envir = e)
S <- as.data.table(e$SDOH_Data_All)[, .(Cohort = std_study(Study), Location = as.character(Location))]
locs <- unique(S[!is.na(Location)])

variants <- c(); summ <- list()
for (ctr in sort(unique(H$Center))) for (mx in c("Urine", "Blood")) {
  d <- H[Center == ctr & (Matrix_Class == mx | Analyte_Code %in% c("CRE", "SG"))]
  coh <- unique(d[Matrix_Class == mx, Cohort])
  n_loc <- uniqueN(locs[Cohort %in% coh, Location])
  summ[[length(summ) + 1]] <- data.table(Center = ctr, Matrix = mx, cohorts = paste(sort(coh), collapse = ", "),
                                         n_cohorts = length(coh), n_locations = n_loc)
  if (n_loc < 2) next
  v <- paste0("lab_", ctr, "_", mx)
  dir.create(file.path("Output/Dilution_Unadjusted", paste0("harmonized_", v)), recursive = TRUE, showWarnings = FALSE)
  Targeted_Data_All <- d[, !c("Cohort", "Center"), with = FALSE]
  save(Targeted_Data_All, file = file.path("Output/Dilution_Unadjusted", paste0("harmonized_", v), "Harmonized_Targeted.RData"))
  variants <- c(variants, v)
}
summ <- rbindlist(summ); fwrite(summ, file.path(OUT, "within_lab_strata.csv"))
writeLines(variants, file.path(OUT, "within_lab_variants.txt"))
print(summ)
