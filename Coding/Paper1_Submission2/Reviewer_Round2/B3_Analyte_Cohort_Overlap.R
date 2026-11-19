# B3 -- Analyte x cohort overlap (Round 2: Reviewer 1; Reviewer 4 comment 2)
#
# R1: "provide a clear analyte-by-cohort overlap matrix, quantify how many
#      chemicals are measured in sufficient numbers across regions"
# R4: "show, for each chemical class, how many cohorts and participants come
#      from each continent"
#
# Uses the same exposure rule as Hero_Numbers.Rmd (A4): drop CRE/SG, keep only
# codes the crosswalk places in the 11 exposure classes, fold NHBCS_Cot into
# NHBCS. TEDDY and ZIP participants are assigned to their own country, not the
# study's, so the continent counts follow the paper's site-level geography.
#
# Read-only on Input/. Writes only to Output/Reviewer_Round2/B3/.
# Run from the project root:  Rscript Coding/Paper1_Submission2/Reviewer_Round2/B3_Analyte_Cohort_Overlap.R

suppressPackageStartupMessages({
  library(data.table); library(ggplot2); library(patchwork)
})

HARM <- "Input/Harmonized_Datasets"
OUT  <- "Output/Reviewer_Round2/B3"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

load(file.path(HARM, "Harmonized_Targeted.RData"))          # Targeted_Data_All
e <- new.env(); load(file.path(HARM, "Harmonized_SDOH.RData"), envir = e)
S  <- as.data.table(e$SDOH_Data_All)
CW <- fread(file.path(HARM, "Targeted_Analyte_Category_Mapping.csv"))

CLASSES_11 <- c("FlameRetardants","PAH","Parabens","PCB","Pesticides","PFAS",
                "Phenols","Phthalates","Smoking","TraceElements","VOCs")
fold_cot <- function(x) sub("_Cot$", "", x)

H <- as.data.table(Targeted_Data_All)
H[, Std := fold_cot(sub("_Targeted$", "", Study))]
H[, Population := fifelse(!is.na(Child_PID), "Child", "Adult")]
H[, pid := fcoalesce(as.character(Child_PID), as.character(Mat_PID))]
H <- H[!Analyte_Code %in% c("CRE", "SG")]
exposure_codes <- intersect(unique(H$Analyte_Code), CW[Category %in% CLASSES_11, Code])
H <- H[Analyte_Code %in% exposure_codes & !is.na(Concentration)]
H <- merge(H, CW[, .(Analyte_Code = Code, Category)], by = "Analyte_Code")

# ---- geography ---------------------------------------------------------------
# Study-level location, then overridden per participant for the two multi-
# country cohorts. ZIP Country is missing for ~100 participants; Site_Location
# codes recover it (2 Brazil, 6 Guatemala, 7 Nicaragua, 8 Puerto Rico).
CONTINENT <- c(
  Accra = "Africa", South_Africa = "Africa", Uganda = "Africa",
  Bangladesh = "Asia", Thailand = "Asia",
  Finland = "Europe", Germany = "Europe", Sweden = "Europe",
  Mexico_City = "Latin America", Tobago = "Latin America", Suriname = "Latin America",
  Chile = "Latin America", Ecuador = "Latin America", Brazil = "Latin America",
  Guatemala = "Latin America", Nicaragua = "Latin America",
  Puerto_Rico = "Latin America")            # US territory; grouped geographically
continent_of <- function(loc) fifelse(loc %in% names(CONTINENT), CONTINENT[loc], "North America (US)")

S[, Std := fold_cot(sub("_SDOH$", "", Study))]
S[, pid := fcoalesce(as.character(Child_PID), as.character(Mat_PID))]
zip_site <- c(`2` = "Brazil", `6` = "Guatemala", `7` = "Nicaragua", `8` = "Puerto_Rico")
S[, Place := as.character(Location)]
S[Std == "ZIP", Place := fcoalesce(gsub(" ", "_", Country), zip_site[as.character(Site_Location)])]
person_place <- unique(S[Std %in% c("TEDDY", "ZIP") & !is.na(Place), .(Std, pid, Place)])
study_place  <- unique(S[!Std %in% c("TEDDY", "ZIP"), .(Std, StudyPlace = as.character(Location))])
stopifnot(!anyDuplicated(study_place$Std))

H <- merge(H, person_place, by = c("Std", "pid"), all.x = TRUE)
H <- merge(H, study_place,  by = "Std", all.x = TRUE)
H[is.na(Place), Place := StudyPlace]
n_unplaced <- uniqueN(H[is.na(Place), .(Std, pid)])
H[is.na(Place), Place := fifelse(Std == "TEDDY", "TEDDY_unknown_site", "ZIP_unknown_country")]
H[, Continent := fifelse(grepl("unknown", Place), "Unassigned", continent_of(Place))]

# Cohort label used in the matrix: TEDDY/ZIP split by country, like the paper.
H[, Unit := fifelse(Std %in% c("TEDDY", "ZIP"), paste0(Std, "-", Place), Std)]

# ---- 1. long overlap table ---------------------------------------------------
cell <- H[, .(N = uniqueN(paste(Population, pid))),
          by = .(Analyte_Code, Category, Std, Unit, Place, Continent)]
fwrite(cell, file.path(OUT, "overlap_long_analyte_by_cohort.csv"))

# Wide presence matrix (N participants per cell), one row per analyte.
wide <- dcast(cell[, .(N = sum(N)), by = .(Category, Analyte_Code, Std)],
              Category + Analyte_Code ~ Std, value.var = "N", fill = 0)
fwrite(wide, file.path(OUT, "overlap_matrix_analyte_by_cohort.csv"))

# ---- 2. per-analyte breadth --------------------------------------------------
breadth <- cell[, .(n_cohorts   = uniqueN(Std),
                    n_places    = uniqueN(Place),
                    n_continents = uniqueN(Continent[Continent != "Unassigned"]),
                    has_nonUS   = any(!Continent %in% c("North America (US)", "Unassigned")),
                    has_US      = any(Continent == "North America (US)"),
                    N_total     = sum(N),
                    N_nonUS     = sum(N[!Continent %in% c("North America (US)", "Unassigned")])),
                by = .(Analyte_Code, Category)]
fwrite(breadth[order(-n_cohorts)], file.path(OUT, "analyte_breadth.csv"))

thresholds <- breadth[, .(
  analytes                    = .N,
  in_ge2_cohorts              = sum(n_cohorts >= 2),
  in_ge3_cohorts              = sum(n_cohorts >= 3),
  in_ge5_cohorts              = sum(n_cohorts >= 5),
  in_ge10_cohorts             = sum(n_cohorts >= 10),
  in_ge2_continents           = sum(n_continents >= 2),
  in_ge3_continents           = sum(n_continents >= 3),
  US_and_nonUS                = sum(has_US & has_nonUS),
  US_only                     = sum(has_US & !has_nonUS),
  nonUS_only                  = sum(!has_US & has_nonUS))]
fwrite(thresholds, file.path(OUT, "overlap_thresholds_overall.csv"))

# ---- 3. per-class table by continent (R4 comment 2) --------------------------
class_cont <- H[, .(participants = uniqueN(paste(Std, Population, pid)),
                    cohorts      = uniqueN(Std)),
                by = .(Category, Continent)]
class_cont_w <- dcast(class_cont, Category ~ Continent,
                      value.var = c("cohorts", "participants"), fill = 0)
class_sum <- breadth[, .(analytes = .N,
                         analytes_ge2_cohorts    = sum(n_cohorts >= 2),
                         analytes_ge2_continents = sum(n_continents >= 2),
                         analytes_US_and_nonUS   = sum(has_US & has_nonUS)),
                     by = Category]
class_tot <- H[, .(cohorts_total = uniqueN(Std),
                   participants_total = uniqueN(paste(Std, Population, pid)),
                   pct_participants_nonUS = round(100 * uniqueN(paste(Std, Population, pid)[
                     !Continent %in% c("North America (US)", "Unassigned")]) /
                     uniqueN(paste(Std, Population, pid)), 1)),
               by = Category]
class_table <- Reduce(function(a, b) merge(a, b, by = "Category"),
                      list(class_sum, class_tot, class_cont_w))[order(-participants_total)]
fwrite(class_table, file.path(OUT, "class_by_continent_table.csv"))

# Cohort-level: which cohorts contribute to more than one class? (R4: "only two
# cohorts outside the US contribute multiple chemicals")
cohort_breadth <- H[, .(n_analytes = uniqueN(Analyte_Code), n_classes = uniqueN(Category),
                        classes = paste(sort(unique(Category)), collapse = ", "),
                        participants = uniqueN(paste(Population, pid)),
                        continents = paste(sort(unique(Continent)), collapse = "; ")),
                    by = Std][order(-n_analytes)]
fwrite(cohort_breadth, file.path(OUT, "cohort_breadth.csv"))

# ---- 4. figure ---------------------------------------------------------------
cont_levels <- c("North America (US)", "Latin America", "Europe", "Africa", "Asia")
cont_pal <- c("North America (US)" = "#4C72B0", "Latin America" = "#DD8452",
              "Europe" = "#55A868", "Africa" = "#C44E52", "Asia" = "#8172B3")

# Each cohort's primary continent (TEDDY is the only cohort that spans two; the
# matrix shows it once and the per-class table carries its split).
cohort_cont <- H[Continent != "Unassigned", .N, by = .(Std, Continent)][
  order(Std, -N)][, .SD[1], by = Std]
cohort_cont[Std == "TEDDY", Continent := "Europe"]
cohort_order <- cohort_breadth[, .(Std, n_analytes)][cohort_cont, on = "Std"][
  , Continent := factor(Continent, cont_levels)][order(Continent, -n_analytes), Std]

an_order <- breadth[order(Category, -n_cohorts, Analyte_Code), Analyte_Code]
mat <- cell[, .(N = sum(N)), by = .(Analyte_Code, Category, Std)]
mat <- merge(mat, cohort_cont[, .(Std, Continent)], by = "Std")
mat[, `:=`(Std = factor(Std, cohort_order), Analyte_Code = factor(Analyte_Code, rev(an_order)),
           Continent = factor(Continent, cont_levels))]

pA <- ggplot(mat, aes(Std, Analyte_Code, fill = Continent, alpha = log10(N))) +
  geom_tile(colour = "white", linewidth = 0.1) +
  facet_grid(Category ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_manual(values = cont_pal, drop = FALSE) +
  scale_alpha_continuous(range = c(0.35, 1), name = "log10 N") +
  labs(title = "Analytes measured by cohort", x = NULL, y = NULL) +
  theme_minimal(base_size = 7) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 3.2), panel.grid = element_blank(),
        strip.placement = "outside", strip.text.y.left = element_text(angle = 0, size = 6),
        panel.spacing.y = unit(1, "pt"), legend.position = "bottom",
        legend.box = "vertical", legend.box.just = "left")   # side by side, Asia was clipped

br <- copy(breadth)
br[, scope := fifelse(has_US & has_nonUS, "US and non-US",
                      fifelse(has_US, "US only", "Non-US only"))]
pB <- ggplot(br, aes(pmin(n_cohorts, 10), fill = scope)) +
  geom_bar() +
  scale_x_continuous(breaks = 1:10, labels = c(1:9, "10+")) +
  scale_fill_manual(values = c("US only" = "#4C72B0", "US and non-US" = "#937860",
                               "Non-US only" = "#DD8452")) +
  labs(title = "Number of cohorts measuring each analyte", x = "Cohorts per analyte", y = "Analytes", fill = NULL) +
  theme_minimal(base_size = 8) + theme(legend.position = "bottom")

cc <- class_cont[Continent != "Unassigned"]
cc[, Continent := factor(Continent, rev(cont_levels))]
pC <- ggplot(cc, aes(participants, reorder(Category, participants, sum), fill = Continent)) +
  geom_col() +
  scale_fill_manual(values = cont_pal, breaks = cont_levels) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.05)), labels = scales::comma) +
  labs(title = "Participants by exposure class and continent", x = "Participants", y = NULL, fill = NULL) +
  theme_minimal(base_size = 8) + theme(legend.position = "bottom")

fig <- pA | (pB / pC)

# Panel titles/subtitles and in-plot labels removed 2026-11-03 (author: easier to read);
# panels carry A/B/C tags only.
fig <- fig + plot_layout(widths = c(1.3, 1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT, "SuppFig_Analyte_Cohort_Overlap.pdf"), fig, width = 13, height = 16)
ggsave(file.path(OUT, "SuppFig_Analyte_Cohort_Overlap.png"), fig, width = 13, height = 16, dpi = 300)

# ---- 5. console summary ------------------------------------------------------
cat("\nExposure analytes:", length(exposure_codes), " classes:", uniqueN(H$Category),
    " cohorts:", uniqueN(H$Std), "\n")
cat("Participants with no TEDDY/ZIP site assignment:", n_unplaced, "\n\n")
print(thresholds)
cat("\nPer class:\n"); print(class_table[, .(Category, analytes, analytes_ge2_cohorts,
                                              analytes_ge2_continents, analytes_US_and_nonUS,
                                              cohorts_total, participants_total, pct_participants_nonUS)])
cat("\nNon-US cohorts:\n")
print(cohort_breadth[!grepl("^North America \\(US\\)$", continents), .(Std, n_analytes, n_classes, classes, continents)])
