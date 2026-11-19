## ============================================================================
## 01_prep_matrix_aware_input.R
##
## Matrix-aware analysis pipeline -- SEPARATE from the published Paper1 pipeline.
## This script is READ-ONLY with respect to everything that already exists:
## it reads Input/Harmonized_Datasets/*.RData and writes only under
## Output/Matrix_Stratified/. Nothing in Input/ or Output/Paper1_* is modified.
##
## Motivation. Sample_Matrix / Matrix_Class were added to Harmonized_Targeted in
## Aug 2026. They expose a confound the published analysis could not see:
## matrix is PERFECTLY confounded with location. Across every exposure class
## measured in more than one matrix (TraceElements, PFAS, Pesticides,
## FlameRetardants, Smoking, Inflammation) exactly ZERO locations contribute
## more than one matrix. Meanwhile the same analyte differs by a median 56x
## between blood and urine, and up to 965x (Mn), 781x (Pb), 398x (Se).
##
## So a model asked "does location explain concentration?" is partly answering
## "was this cohort sampled in blood or urine?". This pipeline separates the two
## by analysing WITHIN matrix.
##
## Implements items 2-5 of the agreed plan; item 1 (the stratified location
## model) is 02_run_stage_a_by_matrix.sh, which consumes what this writes.
## ============================================================================

suppressPackageStartupMessages({
  library(tidyverse); library(data.table)
})

ROOT <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
HARM <- file.path(ROOT, "Input", "Harmonized_Datasets")
OUT  <- file.path(ROOT, "Output", "Matrix_Stratified")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

load(file.path(HARM, "Harmonized_Targeted.RData"))   # Targeted_Data_All
e <- new.env(); load(file.path(HARM, "Harmonized_SDOH.RData"), envir = e)
SDOH_Data_All <- get(ls(e)[1], envir = e)

H <- as.data.table(Targeted_Data_All)
stopifnot(all(c("Sample_Matrix", "Matrix_Class") %in% names(H)))
log_lines <- c()
say <- function(...) { m <- sprintf(...); message(m); log_lines <<- c(log_lines, m) }

say("input: %s rows, %s studies", format(nrow(H), big.mark=","), uniqueN(H$Study))

## ---------------------------------------------------------------------------
## ITEM 3 -- retire the _BR / _DBS suffix hack.
##
## Those suffixes exist only because there was no matrix column: PFOA measured
## in breast milk had to be renamed PFOA_BR to stop it averaging into serum
## PFOA. With Matrix_Class the same separation is expressed by the data model
## instead of by mangling the identifier, which also means these analytes stop
## being invisible to the crosswalk and to every category-based analysis.
## ---------------------------------------------------------------------------
n_suffixed <- H[grepl("_(BR|DBS)$", Analyte_Code), .N]
H[, Analyte_Code_Original := Analyte_Code]
H[, Analyte_Code := sub("_(BR|DBS)$", "", Analyte_Code)]
say("ITEM 3: restored %s rows from _BR/_DBS suffixes to canonical codes (%s codes)",
    format(n_suffixed, big.mark=","),
    uniqueN(H[grepl("_(BR|DBS)$", Analyte_Code_Original), Analyte_Code_Original]))

## ---------------------------------------------------------------------------
## ITEM 4 -- drop non-biological media.
## PROGRESS ships house dust. It is a legitimate exposure measurement but it is
## not a biomarker of internal dose, and pooling it into concentration
## distributions alongside urine and serum is a category error.
## ---------------------------------------------------------------------------
n_dust <- H[Matrix_Class == "Non-biological", .N]
dust_studies <- H[Matrix_Class == "Non-biological", paste(unique(Study), collapse=", ")]
H <- H[is.na(Matrix_Class) | Matrix_Class != "Non-biological"]
say("ITEM 4: dropped %s non-biological rows (%s)", format(n_dust, big.mark=","),
    if (nzchar(dust_studies)) dust_studies else "none")

## ---------------------------------------------------------------------------
## ITEM 5 -- dilution-correction guard.
##
## Creatinine and specific-gravity correction are only interpretable for urine.
## The harmonized table is already corrected, so this cannot be undone here --
## it is reported as a CHECK. Any cohort that carries a dilution marker while
## being a non-urine matrix is flagged for review rather than silently used.
## ---------------------------------------------------------------------------
dil <- H[Analyte_Code %in% c("CRE", "SG"),
         .(dilution_rows = .N), by = .(Study, Matrix_Class)]
bad <- dil[!is.na(Matrix_Class) & Matrix_Class != "Urine"]
if (nrow(bad)) {
  say("ITEM 5: %d cohort/matrix combinations carry a dilution marker outside urine:", nrow(bad))
  for (i in seq_len(nrow(bad)))
    say("         %s (%s): %s rows", bad$Study[i], bad$Matrix_Class[i], bad$dilution_rows[i])
} else {
  say("ITEM 5: PASS -- every dilution marker sits on a urine row")
}
fwrite(dil, file.path(OUT, "dilution_marker_check.csv"))

## ---------------------------------------------------------------------------
## ITEM 2 -- one matrix per analyte.
##
## Reference intervals pool concentrations across cohorts. An analyte measured
## in two matrices produces a bimodal "interval" that describes no real
## population. For each analyte the dominant matrix (most rows) is kept and the
## minority matrix is split into its own stratum rather than discarded.
## ---------------------------------------------------------------------------
per_analyte <- H[!is.na(Matrix_Class) & !Analyte_Code %in% c("CRE","SG"),
                 .N, by = .(Analyte_Code, Matrix_Class)]
multi <- per_analyte[, .(n_matrix = uniqueN(Matrix_Class)), by = Analyte_Code][n_matrix > 1]
say("ITEM 2: %d analytes span more than one matrix and would produce bimodal intervals if pooled",
    nrow(multi))
fwrite(per_analyte[Analyte_Code %in% multi$Analyte_Code][order(Analyte_Code, -N)],
       file.path(OUT, "analytes_spanning_matrices.csv"))

## ---------------------------------------------------------------------------
## ITEM 1 input -- write one harmonized set per matrix class.
## Each is a drop-in replacement for Harmonized_Targeted.RData so the UNCHANGED
## Stage A script can be pointed at it. That keeps the method provably identical
## to the published pipeline; only the row set differs.
## ---------------------------------------------------------------------------
classes <- sort(unique(na.omit(H$Matrix_Class)))
say("writing per-matrix inputs for: %s", paste(classes, collapse = ", "))

for (cl in classes) {
  sub <- H[Matrix_Class == cl | Analyte_Code %in% c("CRE","SG")]
  # keep dilution markers only for the cohorts present in this stratum
  sub <- sub[Study %in% unique(H[Matrix_Class == cl, Study])]
  d <- file.path(OUT, paste0("Input_", gsub("[^A-Za-z]", "", cl)))
  dir.create(d, showWarnings = FALSE, recursive = TRUE)
  Targeted_Data_All <- as.data.frame(sub)
  save(Targeted_Data_All, file = file.path(d, "Harmonized_Targeted.RData"))
  save(SDOH_Data_All,     file = file.path(d, "Harmonized_SDOH.RData"))
  say("  %-15s %8s rows  %2d studies  %3d analytes",
      cl, format(nrow(sub), big.mark=","), uniqueN(sub$Study),
      uniqueN(sub[!Analyte_Code %in% c("CRE","SG"), Analyte_Code]))
}

## full matrix-aware table, for the report
Targeted_Data_All <- as.data.frame(H)
save(Targeted_Data_All, file = file.path(OUT, "Harmonized_Targeted_MatrixAware.RData"))
writeLines(log_lines, file.path(OUT, "01_prep_log.txt"))
say("done. nothing outside Output/Matrix_Stratified was written.")
