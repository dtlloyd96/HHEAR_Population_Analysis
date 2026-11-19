# R4_5_geospatial_scatter.R
# Reviewer 4 -- pre-specified location-level consistency checks between
# biomarker location medians and the geospatial layers built by
# R4_5_geospatial_extract.R. Descriptive/ecological: one point per location
# (cohort), Spearman rho with an exact permutation p-value. Run from project root.
#
# Pre-specified pairs (PLAN_R4_Geographic_Determinants.md 5.3; author approved
# 2026-11-04 that these are the reported result, heatmap descriptive):
#   urinary PAH (NAP1, NAP2, PYR1)   ~ PM2.5, NO2 (traffic; replaces PKU-PAH
#                                      emissions, not publicly downloadable), population density
#   urinary Pb, Cd                   ~ PM2.5
#   urinary As                       ~ P(groundwater As > 10 ug/L)
# Pesticide pairs wait for PEST-CHEMGRIDS (needs Earthdata login).
#
# Location median: per location, within one matrix, one value per participant
# (participant median over repeat samples), <LOD set to LOD/sqrt(2), units via
# to_ng_per_mL(). Rows whose units to_ng_per_mL() passes through unconverted
# (pmol/mL, nM, kg/m3, ng/g) are DROPPED here because absolute levels are
# compared across cohorts. Locations with detection frequency < 50% are kept
# but flagged (hollow points). Urine dilution basis (ug/g creatinine, SG 1.020,
# none) is shown by point shape; a sensitivity drops the unadjusted ("none") rows.

suppressMessages({library(dplyr); library(ggplot2); library(tidyr)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
source("Coding/Paper1_Submission2/harmonization_helpers.R")
out <- "Output/Reviewer_Round2/R4_5_determinants/geospatial"

L  <- read.csv(file.path(out, "location_layers_time_matched.csv"))
yw <- read.csv(file.path(out, "location_study_year_weights.csv"))
e <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = e)
s <- new.env(); load("Input/Harmonized_Datasets/Harmonized_SDOH.RData",   envir = s)
TD <- e$Targeted_Data_All; SD <- s$SDOH_Data_All

loc_map <- SD %>% distinct(Study, Child_PID, Mat_PID, Location, Country) %>%
  mutate(Study = sub("_SDOH$", "_Targeted", Study),
         loc_id = case_when(Study == "TEDDY_Targeted" ~ paste0("TEDDY_", Location),
                            Study == "ZIP_Targeted" & !is.na(Country) ~ paste0("ZIP_", Country),
                            TRUE ~ Location))

pairs <- tribble(
  ~analyte, ~matrix, ~layer,
  "NAP1", "Urine", "PM25_ugm3",
  "NAP2", "Urine", "PM25_ugm3",
  "PYR1", "Urine", "PM25_ugm3",
  "NAP1", "Urine", "NO2_ppb",
  "NAP2", "Urine", "NO2_ppb",
  "PYR1", "Urine", "NO2_ppb",
  "NAP1", "Urine", "pop_density_km2",
  "NAP2", "Urine", "pop_density_km2",
  "PYR1", "Urine", "pop_density_km2",
  "Pb",   "Urine", "PM25_ugm3",
  "Cd",   "Urine", "PM25_ugm3",
  "As",   "Urine", "GW_As_prob_gt10"
)
pairs <- pairs %>% filter(layer %in% names(L))

ok_units <- c("ng/mL", "ng/ml", "ng/L", "pg/mL", "pg/ml", "ug/dL", "mg/L")
# VDAART reports 1-hydroxypyrene as "PYR" (same metabolite as PYR1 elsewhere);
# remapped here for this analysis only (author approved 2026-11-04). The
# canonical Harmonized_Targeted data and B4 outputs are NOT changed. VDAART's
# "NAP" is combined 1+2-hydroxynaphthalene and is NOT mapped to NAP1/NAP2.
TD <- TD %>% mutate(Analyte_Code = ifelse(Study == "VDAART_Targeted" & Analyte_Code == "PYR", "PYR1", Analyte_Code))
d <- TD %>% filter(Analyte_Code %in% unique(pairs$analyte), Matrix_Class %in% unique(pairs$matrix)) %>%
  mutate(units_ok = is.na(Units) | Units %in% ok_units)
dropped_units <- d %>% filter(!units_ok) %>% count(Study, Analyte_Code, Units)
d <- d %>% filter(units_ok, !is.na(Concentration) | !is.na(LOD)) %>%
  mutate(conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
         det  = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
         val  = ifelse(det, conc, lod / sqrt(2)),
         pid  = coalesce(Mat_PID, Child_PID)) %>% filter(!is.na(val), val > 0)

# attach location (TEDDY by child, ZIP by mother, others by study)
d <- bind_rows(
  d %>% filter(!Study %in% c("TEDDY_Targeted", "ZIP_Targeted")) %>% left_join(loc_map %>% distinct(Study, loc_id), by = "Study"),
  d %>% filter(Study == "TEDDY_Targeted") %>% left_join(loc_map %>% filter(Study == "TEDDY_Targeted") %>% distinct(Child_PID, loc_id), by = "Child_PID"),
  d %>% filter(Study == "ZIP_Targeted") %>% left_join(loc_map %>% filter(Study == "ZIP_Targeted", !is.na(Mat_PID)) %>% distinct(Mat_PID, loc_id), by = "Mat_PID")
) %>% filter(!is.na(loc_id))

loc_med <- function(dd) dd %>%
  group_by(Analyte_Code, Matrix_Class, loc_id, Study, pid) %>%
  summarise(v = median(val), det = mean(det), basis = first(Dilution_Adjustment), .groups = "drop") %>%
  group_by(Analyte_Code, Matrix_Class, loc_id) %>%
  summarise(median = median(v), df = mean(det), n = n(), studies = paste(sort(unique(std_study(Study))), collapse = "+"),
            basis = paste(sort(unique(basis)), collapse = "/"), .groups = "drop")
M_all  <- loc_med(d) %>% mutate(set = "all")
M_adj  <- loc_med(d %>% filter(Matrix_Class != "Urine" | Dilution_Adjustment != "none")) %>% mutate(set = "dilution_adjusted_only")
# Location-analytes with < 20% detection carry no information (e.g. PARENTs NAP1:
# every concentration missing, LOD 1000 ng/L, so the "median" is LOD/sqrt(2)) -- dropped.
M <- bind_rows(M_all, M_adj) %>% filter(n >= 20, df >= 0.2)
write.csv(M, file.path(out, "biomarker_location_medians.csv"), row.names = FALSE)

perm_spearman <- function(x, y, B = 20000) {
  ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]; n <- length(x)
  if (n < 4) return(c(rho = NA, p = NA, n = n))
  r0 <- suppressWarnings(cor(x, y, method = "spearman"))
  set.seed(1); rp <- replicate(B, suppressWarnings(cor(x, sample(y), method = "spearman")))
  c(rho = r0, p = (1 + sum(abs(rp) >= abs(r0) - 1e-12)) / (B + 1), n = n)
}

res <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  p <- pairs[i, ]
  bind_rows(lapply(c("all", "dilution_adjusted_only"), function(st) {
    z <- M %>% filter(set == st, Analyte_Code == p$analyte, Matrix_Class == p$matrix) %>% inner_join(L, by = "loc_id")
    r <- perm_spearman(log10(z[[p$layer]] + ifelse(p$layer == "GW_As_prob_gt10", 0.001, 0)), log10(z$median))
    tibble(analyte = p$analyte, matrix = p$matrix, layer = p$layer, set = st,
           n_locations = r[["n"]], rho = r[["rho"]], perm_p = r[["p"]])
  }))
}))
write.csv(res, file.path(out, "prespecified_pairs_spearman.csv"), row.names = FALSE)
write.csv(dropped_units, file.path(out, "rows_dropped_unconvertible_units.csv"), row.names = FALSE)
print(as.data.frame(res), digits = 3)

# ---- Figure: one panel per pre-specified pair (main set = all) --------------
lab_layer <- c(PM25_ugm3 = "PM2.5, pop-weighted (ug/m3)", NO2_ppb = "NO2, pop-weighted (ppb)", pop_density_km2 = "Population density (per km2)",
               GW_As_prob_gt10 = "P(groundwater As > 10 ug/L), pop-weighted")
fig_d <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) {
  p <- pairs[i, ]
  M %>% filter(set == "all", Analyte_Code == p$analyte, Matrix_Class == p$matrix) %>% inner_join(L, by = "loc_id") %>%
    transmute(panel = paste0(p$analyte, " (", tolower(p$matrix), ") vs ", lab_layer[[p$layer]]),
              analyte = p$analyte, layer = p$layer, x = .data[[p$layer]], y = median, loc_id, studies,
              basis = sub("/.*", " (mixed)", basis), low_df = df < 0.5, assumed = years_assumed)
})) %>% left_join(res %>% filter(set == "all") %>% transmute(analyte, layer,
              lab = sprintf("rho = %.2f, perm p = %.2g, n = %d", rho, perm_p, n_locations)), by = c("analyte", "layer"))
fig_d$panel <- factor(fig_d$panel, levels = unique(fig_d$panel))

g <- ggplot(fig_d, aes(x, y)) +
  geom_point(aes(shape = basis, fill = ifelse(low_df, NA, "detected")), size = 2.4, colour = "grey20") +
  ggrepel::geom_text_repel(aes(label = gsub("_", " ", loc_id)), size = 2.3, colour = "grey30", max.overlaps = 30, seed = 1) +
  geom_text(data = distinct(fig_d, panel, lab), aes(label = lab), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 2.8, inherit.aes = FALSE) +
  scale_x_log10() + scale_y_log10() +
  scale_shape_manual(values = c(21, 22, 24, 23, 25), name = "Urine dilution basis") +
  scale_fill_manual(values = c(detected = "grey55"), na.value = "white", guide = "none") +
  facet_wrap(~panel, scales = "free", ncol = 3) +
  labs(x = "Location-level layer value (matched to sampling years)", y = "Location median biomarker (ng/mL or ug/g creatinine)",
       caption = "One point per location. Hollow = detection frequency < 50%. Ecological, n = locations; consistency check, not effect estimate.") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom", strip.text = element_text(size = 7.5))
ggsave(file.path(out, "SuppFig_Geospatial_Prespecified_Pairs.pdf"), g, width = 11, height = 3.4 * ceiling(nlevels(fig_d$panel) / 3) + 0.8)
ggsave(file.path(out, "SuppFig_Geospatial_Prespecified_Pairs.png"), g, width = 11, height = 3.4 * ceiling(nlevels(fig_d$panel) / 3) + 0.8, dpi = 200)

# ---- Descriptive heatmap: locations x layers (z-scored) ---------------------
hm <- L %>% select(loc_id, type, any_of(c("PM25_ugm3", "NO2_ppb", "pop_density_km2", "polluting_fuel_pct", "GW_As_prob_gt10"))) %>%
  mutate(pop_density_km2 = log10(pop_density_km2)) %>%
  pivot_longer(-c(loc_id, type), names_to = "layer") %>% group_by(layer) %>% mutate(z = as.numeric(scale(value))) %>% ungroup()
ord <- L %>% arrange(PM25_ugm3) %>% pull(loc_id)
gh <- ggplot(hm, aes(layer, factor(loc_id, levels = ord), fill = z)) + geom_tile(colour = "white") +
  geom_text(aes(label = signif(ifelse(layer == "pop_density_km2", 10^value, value), 2)), size = 2.3) +
  scale_fill_gradient2(low = "#3b6fb6", mid = "white", high = "#c0392b", name = "z") +
  scale_x_discrete(labels = c(PM25_ugm3 = "PM2.5", NO2_ppb = "NO2", pop_density_km2 = "Pop. density", polluting_fuel_pct = "% polluting fuel", GW_As_prob_gt10 = "P(GW As>10)")) +
  labs(x = NULL, y = NULL, caption = "Population-weighted, matched to sampling years where time-varying. Descriptive.") +
  theme_minimal(base_size = 8)
ggsave(file.path(out, "SuppFig_Geospatial_Layer_Heatmap.pdf"), gh, width = 5.2, height = 7, device = cairo_pdf)
ggsave(file.path(out, "SuppFig_Geospatial_Layer_Heatmap.png"), gh, width = 5.2, height = 7, dpi = 200)

# ---- Check: does time-matching matter? (justifies single-year NO2, 2015) -----
# PM2.5 is available for every year, so compare the time-matched PM2.5 with a
# fixed-2015 PM2.5: rank agreement across locations, and the pre-specified PAH
# correlations under each version.
pm <- read.csv(file.path(out, "layer_PM25_location_year.csv"))
L15 <- L %>% select(-PM25_ugm3) %>% left_join(pm %>% filter(year == 2015) %>% select(loc_id, PM25_ugm3 = value), by = "loc_id")
tm_check <- tibble(check = "PM2.5 time-matched vs fixed-2015, Spearman across locations",
                   value = cor(L$PM25_ugm3, L15$PM25_ugm3[match(L$loc_id, L15$loc_id)], method = "spearman", use = "complete.obs"))
pah_cmp <- bind_rows(lapply(c("NAP1", "NAP2", "PYR1"), function(a) {
  z <- M %>% filter(set == "all", Analyte_Code == a, Matrix_Class == "Urine")
  m1 <- z %>% inner_join(L,   by = "loc_id"); m2 <- z %>% inner_join(L15, by = "loc_id")
  tibble(analyte = a,
         rho_time_matched = perm_spearman(log10(m1$PM25_ugm3), log10(m1$median))[["rho"]],
         rho_fixed_2015   = perm_spearman(log10(m2$PM25_ugm3), log10(m2$median))[["rho"]], n = nrow(m1))
}))
write.csv(bind_rows(tm_check %>% mutate(analyte = NA), pah_cmp %>% mutate(check = "PAH ~ PM2.5 rho, time-matched vs fixed 2015")),
          file.path(out, "check_time_matching_PM25.csv"), row.names = FALSE)
print(tm_check); print(pah_cmp)
