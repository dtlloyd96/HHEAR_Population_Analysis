# R4_5_fingerprints.R
# Reviewer 4 -- biomarker source fingerprints: within-person ratios of related
# chemicals that indicate the source, computed in every cohort that measured
# both members of a pair (no questionnaire needed). Ratios within one person,
# one sample basis and one lab cancel urine dilution and much of the lab /
# protocol difference, so they compare more fairly across locations than levels.
#
# Fingerprints (participant-level, from participant median concentrations):
#   PAH  log10(1-naphthol / 2-naphthol)          higher with tobacco smoke (and naphthalene)
#        log10(1-OH-pyrene / (1- + 2-naphthol))   pyrene share: diet / other combustion vs tobacco
#   PHT  DEHTP share = (MEHHTP + MECPTP) / (MEHHTP + MECPTP + MEHHP + MEOHP + MECPP + MEHP)
#        replacement of DEHP by DEHTP (market / regulatory shift)
#   BIS  BPS share = BPS / (BPA + BPS)            BPA phase-out and replacement
#   PAR  log10(methylparaben / propylparaben)    food / pharma (methyl) vs cosmetics (propyl)
#   OP   diethyl share = (DEP + DETP + DEDP) / (all six DAPs)   chlorpyrifos-type vs malathion-type OPs
# Values: to_ng_per_mL; < LOD -> LOD/sqrt(2). A fingerprint needs >= 50%
# detection of each member in that cohort-location, and >= 20 participants.
# Validation: PAH fingerprints vs reported smoking within DCHS, SICAS, ELEMENT,
# PARENTs (same smoking definitions as R4_5_source_meta_analysis.R).
# Shares are not molar-corrected (metabolite molecular weights are similar
# within each family); they are relative indices for comparison, not exact fractions.
# Writes Output/Reviewer_Round2/R4_5_determinants/fingerprints/. Run from project root.

suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(patchwork)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
source("Coding/Paper1_Submission2/harmonization_helpers.R")
out <- "Output/Reviewer_Round2/R4_5_determinants/fingerprints"; dir.create(out, recursive = TRUE, showWarnings = FALSE)
t <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = t); TD <- t$Targeted_Data_All
s <- new.env(); load("Input/Harmonized_Datasets/Harmonized_SDOH.RData", envir = s); SD <- s$SDOH_Data_All
loc_map <- SD %>% distinct(Study, Child_PID, Mat_PID, Location, Country) %>%
  mutate(Study = sub("_SDOH$", "_Targeted", Study),
         loc_id = case_when(Study == "TEDDY_Targeted" ~ paste0("TEDDY_", Location),
                            Study == "ZIP_Targeted" & !is.na(Country) ~ paste0("ZIP_", Country), TRUE ~ Location))
us_locs <- c("Atlanta","Sacramento","Denver","Baltimore","SanFrancisco","Boston","Los_Angeles","Syracuse","New_York","Michigan","California",
             "Washington_State","New_Hampshire","Massachusetts","NorthEast_US","Midwestern_US","United_States","TEDDY_United_States")

NEED <- c("NAP1", "NAP2", "PYR1", "MEHHTP", "MECPTP", "MEHHP", "MEOHP", "MECPP", "MEHP", "BPA", "BPS", "MEPB", "PRPB",
          "DMP", "DMTP", "DMDP", "DEP", "DETP", "DEDP")
P <- TD %>% mutate(Analyte_Code = ifelse(Study == "VDAART_Targeted" & Analyte_Code == "PYR", "PYR1", Analyte_Code)) %>%
  filter(Analyte_Code %in% NEED, Matrix_Class == "Urine") %>%
  mutate(id = coalesce(Mat_PID, Child_PID), panel = ifelse(!is.na(Mat_PID), "Mat", "Child"),
         conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
         det = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
         val = ifelse(det, conc, lod / sqrt(2)), yr = suppressWarnings(as.numeric(Sample_Collection_Year))) %>% filter(!is.na(val), val > 0)
P <- bind_rows(
  P %>% filter(!Study %in% c("TEDDY_Targeted", "ZIP_Targeted")) %>% left_join(loc_map %>% distinct(Study, loc_id), by = "Study"),
  P %>% filter(Study == "TEDDY_Targeted") %>% left_join(loc_map %>% filter(Study == "TEDDY_Targeted") %>% distinct(Child_PID, loc_id), by = "Child_PID"),
  P %>% filter(Study == "ZIP_Targeted") %>% left_join(loc_map %>% filter(Study == "ZIP_Targeted", !is.na(Mat_PID)) %>% distinct(Mat_PID, loc_id), by = "Mat_PID")
) %>% filter(!is.na(loc_id), loc_id != "Central/SouthAmerica")   # ZIP women with no country recorded
W <- P %>% group_by(Study, loc_id, panel, id, Analyte_Code) %>% summarise(v = median(val), det = mean(det), yr = median(yr, na.rm = TRUE), .groups = "drop")
DF <- W %>% group_by(Study, loc_id, Analyte_Code) %>% summarise(df = mean(det), .groups = "drop")
Wv <- W %>% select(Study, loc_id, panel, id, Analyte_Code, v) %>% pivot_wider(names_from = Analyte_Code, values_from = v)
yrs <- W %>% group_by(Study, panel, id) %>% summarise(yr = suppressWarnings(median(yr, na.rm = TRUE)), .groups = "drop")
Wv <- Wv %>% left_join(yrs, by = c("Study", "panel", "id"))
ok <- function(st, lo, members) all(sapply(members, function(m) { z <- DF$df[DF$Study == st & DF$loc_id == lo & DF$Analyte_Code == m]; length(z) && z >= 0.5 }))

sumc <- function(df, cols) { cols <- intersect(cols, names(df)); if (!length(cols)) return(rep(NA_real_, nrow(df))); x <- as.matrix(df[, cols]); ifelse(rowSums(is.na(x)) == length(cols), NA, rowSums(x, na.rm = TRUE)) }
FP <- list(
  pah_nap_ratio = list(label = "1-naphthol / 2-naphthol (log10)", members = c("NAP1", "NAP2"), f = function(d) log10(d$NAP1 / d$NAP2)),
  pah_pyrene_share = list(label = "1-OH-pyrene / naphthols (log10)", members = c("NAP1", "NAP2", "PYR1"), f = function(d) log10(d$PYR1 / (d$NAP1 + d$NAP2))),
  dehtp_share = list(label = "DEHTP share of DEHP + DEHTP", members = c("MEHHTP", "MECPTP", "MEHHP", "MECPP"),
                     f = function(d) { a <- sumc(d, c("MEHHTP", "MECPTP")); b <- sumc(d, c("MEHHP", "MEOHP", "MECPP", "MEHP")); a / (a + b) }),
  bps_share = list(label = "BPS share of BPA + BPS", members = c("BPA", "BPS"), f = function(d) d$BPS / (d$BPA + d$BPS)),
  paraben_ratio = list(label = "Methyl- / propylparaben (log10)", members = c("MEPB", "PRPB"), f = function(d) log10(d$MEPB / d$PRPB)),
  diethyl_share = list(label = "Diethyl share of OP metabolites", members = c("DEP", "DMP"),
                       f = function(d) { a <- sumc(d, c("DEP", "DETP", "DEDP")); b <- sumc(d, c("DMP", "DMTP", "DMDP")); a / (a + b) }))
IND <- bind_rows(lapply(names(FP), function(k) {
  fp <- FP[[k]]
  bind_rows(lapply(split(Wv, list(Wv$Study, Wv$loc_id), drop = TRUE), function(d) {
    if (!all(fp$members %in% names(d)) || !ok(d$Study[1], d$loc_id[1], fp$members)) return(NULL)
    v <- fp$f(d); keep <- is.finite(v); if (sum(keep) < 20) return(NULL)
    tibble(fingerprint = k, label = fp$label, Study = d$Study[keep], loc_id = d$loc_id[keep], panel = d$panel[keep], id = d$id[keep], yr = d$yr[keep], value = v[keep])
  }))
}))
write.csv(IND, file.path(out, "fingerprints_participant.csv"), row.names = FALSE)
LOC <- IND %>% group_by(fingerprint, label, Study, loc_id) %>%
  summarise(n = n(), median = median(value), q25 = quantile(value, .25), q75 = quantile(value, .75), year = median(yr, na.rm = TRUE), .groups = "drop") %>%
  mutate(region = ifelse(loc_id %in% us_locs, "US", ifelse(grepl("^TEDDY_(Finland|Sweden|Germany)", loc_id), "Europe", "Latin America & Africa")),
         cohort = sub("_Targeted", "", Study))
write.csv(LOC, file.path(out, "fingerprints_by_location.csv"), row.names = FALSE)

# ---- validation: PAH fingerprints vs smoking ------------------------------------
smk_sdoh <- SD %>% mutate(sm = suppressWarnings(as.numeric(coalesce(as.character(Mat_Smoking), as.character(Mat_Smoker))))) %>% filter(!is.na(sm)) %>%
  mutate(x = as.integer(sm >= 1), Study = sub("_SDOH$", "_Targeted", Study))
smk <- bind_rows(smk_sdoh %>% filter(!is.na(Mat_PID)) %>% transmute(Study, panel = "Mat", id = as.character(Mat_PID), x),
                 smk_sdoh %>% filter(!is.na(Child_PID)) %>% transmute(Study, panel = "Child", id = as.character(Child_PID), x))
e <- new.env(); load("Input/DCHS_00496_Data_All.RData", envir = e); C <- e[[ls(e)[1]]]$DCHS_EPI
smk <- bind_rows(smk %>% filter(Study != "DCHS_Targeted"),
                 C %>% transmute(Study = "DCHS_Targeted", panel = "Mat", id = as.character(mother_pid), x = smoke_during_pregnancy),
                 C %>% transmute(Study = "DCHS_Targeted", panel = "Child", id = as.character(child_pid), x = smoke_during_pregnancy)) %>% distinct(Study, panel, id, .keep_all = TRUE)
VAL <- IND %>% filter(fingerprint %in% c("pah_nap_ratio", "pah_pyrene_share")) %>% inner_join(smk, by = c("Study", "panel", "id")) %>%
  group_by(fingerprint, Study, panel) %>% filter(n_distinct(x) == 2, min(table(x)) >= 5) %>%
  summarise(n = n(), n_smoke = sum(x), diff = coef(lm(value ~ x))[2], lo = confint(lm(value ~ x))[2, 1], hi = confint(lm(value ~ x))[2, 2],
            p = summary(lm(value ~ x))$coefficients[2, 4], .groups = "drop") %>%
  group_by(fingerprint, Study) %>% arrange(desc(panel == "Mat"), desc(n)) %>% slice(1) %>% ungroup() %>% mutate(cohort = sub("_Targeted", "", Study))  # mother panel preferred (smoking is maternal)
write.csv(VAL, file.path(out, "pah_fingerprint_smoking_validation.csv"), row.names = FALSE)

# ---- time trends of replacement shares ------------------------------------------
TR <- IND %>% filter(fingerprint %in% c("dehtp_share", "bps_share"), !is.na(yr)) %>% group_by(fingerprint) %>%
  summarise(n = n(), cohorts = n_distinct(Study), slope_per_year = coef(lm(value ~ yr + Study))["yr"],
            p = summary(lm(value ~ yr + Study))$coefficients["yr", 4], .groups = "drop")
write.csv(TR, file.path(out, "replacement_share_time_trends.csv"), row.names = FALSE)
print(LOC %>% select(fingerprint, cohort, loc_id, n, median, year) %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% arrange(fingerprint, median) %>% as.data.frame())
print(VAL %>% mutate(across(where(is.numeric), ~ signif(.x, 3))) %>% as.data.frame()); print(TR)

# ---- figure ------------------------------------------------------------------------
cols <- c(US = "#3b6fb6", Europe = "#2e8b57", `Latin America & Africa` = "#c0392b")
pl <- function(k, title) {
  z <- LOC %>% filter(fingerprint == k) %>% mutate(lab = paste0(gsub("_", " ", sub("^(TEDDY|ZIP)_", "", loc_id)), " (", cohort, ")"))
  ggplot(z, aes(median, reorder(lab, median), colour = region)) + geom_linerange(aes(xmin = q25, xmax = q75)) + geom_point(size = 2) +
    scale_colour_manual(values = cols, name = NULL, drop = FALSE) + labs(title = title, x = z$label[1], y = NULL) + theme_bw(base_size = 8) + theme(legend.position = "none")
}
trend <- IND %>% filter(fingerprint == "dehtp_share", !is.na(yr)) %>% group_by(Study, loc_id, yr) %>% summarise(m = median(value), n = n(), .groups = "drop") %>% filter(n >= 10)
pT <- ggplot(trend, aes(yr, m, colour = sub("_Targeted", "", Study))) + geom_line() + geom_point(aes(size = n)) + scale_size(range = c(0.8, 3), guide = "none") +
  labs(title = "F. DEHTP share by sampling year", x = "Sampling year", y = "DEHTP share (median)", colour = NULL) + theme_bw(base_size = 8) + theme(legend.position = "right")
g <- (pl("pah_nap_ratio", "A. PAH: 1-naphthol / 2-naphthol") | pl("pah_pyrene_share", "B. PAH: pyrene vs naphthols")) /
     (pl("dehtp_share", "C. Phthalates: DEHTP share") | pl("bps_share", "D. Bisphenols: BPS share")) /
     (pl("diethyl_share", "E. OP pesticides: diethyl share") | pT) +
  plot_annotation(caption = "Points: location median; bars: IQR. Colour: US (blue), Europe (green), Latin America & Africa (red).")
ggsave(file.path(out, "Fig_Fingerprints.png"), g, width = 11, height = 11, dpi = 200)
ggsave(file.path(out, "Fig_Fingerprints.pdf"), g, width = 11, height = 11, device = cairo_pdf)
