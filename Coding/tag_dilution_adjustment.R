## Adds Dilution_Adjustment to the canonical harmonized targeted data by
## comparing it row by row with an unadjusted build of the same script.
## Called by rebuild_harmonized_targeted.sh; see that file for context.
##
## Usage: Rscript tag_dilution_adjustment.R <canonical.RData> <unadjusted.RData>
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 2)

ld <- function(f) { e <- new.env(); load(f, envir = e); as.data.table(e$Targeted_Data_All) }
D <- ld(args[1]); U <- ld(args[2])

# Same script, same inputs: rows must line up exactly.
stopifnot(nrow(D) == nrow(U),
          identical(D$Study, U$Study),
          identical(D$Analyte_Code, U$Analyte_Code),
          identical(fcoalesce(D$Child_PID, D$Mat_PID), fcoalesce(U$Child_PID, U$Mat_PID)))

differs <- !is.na(D$Concentration) & !is.na(U$Concentration) &
  abs(D$Concentration - U$Concentration) > 1e-9 * pmax(1, abs(U$Concentration))

# Which kind of adjustment a cohort uses is fixed by its block in V3. Identify
# it from the data rather than a hand list: SG correction rescales by
# 0.020/(SG-1), creatinine by 100/CRE; the ratio range separates them cleanly
# (SG ratios sit near 1, creatinine ratios spread over ~0.3-20).
ratio <- D$Concentration / U$Concentration
kind <- data.table(Study = D$Study, ratio, differs)[differs == TRUE,
  .(q05 = quantile(ratio, .05, na.rm = TRUE), q95 = quantile(ratio, .95, na.rm = TRUE)), by = Study]
SG_STUDIES <- c("DISCOVER", "MARBLES", "NFGS", "nuMoM2b", "VDAART")   # SG blocks in V3
kind[, method := fifelse(sub("_Targeted$", "", Study) %in% SG_STUDIES,
                         "SG_normalized_1.020", "creatinine_ug_per_g")]

# "none" inside an adjusted cohort means the value was left as reported: either
# no marker for that participant, or a marker exactly at the reference
# (creatinine 100 mg/dL, SG 1.020), where adjustment is the identity and so is
# indistinguishable here. Both cases are on the cohort's own scale.
D[, Dilution_Adjustment := "none"]
D[kind, on = "Study", method := i.method]
D[differs %in% TRUE, Dilution_Adjustment := method]
D[, method := NULL]

cat("Adjustment ratio range by cohort (adjusted / unadjusted):\n")
print(kind[order(method, Study)])
cat("\nRows by Dilution_Adjustment:\n")
print(D[, .N, by = Dilution_Adjustment])
cat("\nAdjusted cohorts with rows left unadjusted (missing marker or non-urine):\n")
print(D[Study %in% kind$Study, .(adjusted = sum(Dilution_Adjustment != "none"),
                                left_as_reported = sum(Dilution_Adjustment == "none")),
        by = Study][left_as_reported > 0][order(-left_as_reported)])

Targeted_Data_All <- D
save(Targeted_Data_All, file = args[1])
cat("\nWrote", args[1], "\n")
