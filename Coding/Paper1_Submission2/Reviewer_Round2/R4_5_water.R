# R4_5_water.R
# Reviewer 4, determinants of geography: drinking water.
#   Step 1  Arsenic
#     1a MEBCO (Bangladesh): individual urinary As vs measured household
#        drinking-water As (`wateras_c1`, units assumed ug/L; the harmonized
#        column name says "Microliter", read as a labelling slip).
#     1b Location level: P(groundwater As > 10 ug/L) x reliance on unmanaged
#        groundwater (Podgorski & Berg 2020, groundwater_utilization_by_country.csv,
#        rural/urban, from WHO/UNICEF JMP / FAO AQUASTAT). Each ~1 km cell is
#        urban if GHS-POP density > 300/km2, else rural (following the paper's
#        affected-population method); exposure = sum(pop * P * reliance) / sum(pop).
#        The US is not in that table: ASSUMED 0.14 (share of US population on
#        private domestic wells, USGS / Johnson et al. 2019), for all US units,
#        urban and rural alike -- verify, and refine by state later. Tobago,
#        Suriname, Thailand, Uganda, Puerto Rico: no value (NA).
#   Step 2  PFAS in US public drinking water: EPA UCMR3 (2013-2015), PFOS/PFOA/
#        PFHxS/PFNA/PFHpA/PFBS at MRLs 0.01-0.09 ug/L. State-level share of
#        public water systems with >=1 detection of PFOS or PFOA (and of any of
#        the six). UCMR3 has no population-served field, so it is unweighted;
#        large systems (Size L, a census) shown separately. City/region units
#        take their state; multi-state NorthEast_US / Midwestern_US average their
#        states; "United_States" takes the national value. Blood PFOS/PFOA location
#        medians, plus an individual mixed model adjusted for sampling year
#        (PFAS fell steeply over 2000-2020, and locations differ in era).
#   Step 3  ZIP (Brazil, Guatemala, Nicaragua, Puerto Rico; one lab): household
#        primary water source (ZIP_EPI_1$tri_prim_water_src; raw codes mapped to
#        the ZIP codebook via the harmonized column: 1 = water utility,
#        2 = private well, 7 = bottled/delivered, 3-6 = other) vs urinary
#        pesticide metabolites. Does water source explain the between-country
#        differences? Compare country effects with and without water source.
# Raw downloads (cache, outside Drive): ~/HHEAR_geo_cache/as_gw/data/,
# ~/HHEAR_geo_cache/ucmr3/. Writes Output/Reviewer_Round2/R4_5_determinants/water/.
# Run from the project root after R4_5_geospatial_extract.R.

suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(sf); library(terra); library(lme4); library(patchwork)})
sf_use_s2(FALSE)
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
source("Coding/Paper1_Submission2/harmonization_helpers.R")
cache <- Sys.getenv("HHEAR_GEO_CACHE", path.expand("~/HHEAR_geo_cache"))
geo <- "Output/Reviewer_Round2/R4_5_determinants/geospatial"
out <- "Output/Reviewer_Round2/R4_5_determinants/water"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
U  <- st_read(file.path(geo, "location_units.gpkg"), quiet = TRUE)
L  <- read.csv(file.path(geo, "location_layers_time_matched.csv"))
e <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = e); TD <- e$Targeted_Data_All
s <- new.env(); load("Input/Harmonized_Datasets/Harmonized_SDOH.RData", envir = s)
loc_map <- s$SDOH_Data_All %>% distinct(Study, Child_PID, Mat_PID, Location, Country) %>%
  mutate(Study = sub("_SDOH$", "_Targeted", Study),
         loc_id = case_when(Study == "TEDDY_Targeted" ~ paste0("TEDDY_", Location),
                            Study == "ZIP_Targeted" & !is.na(Country) ~ paste0("ZIP_", Country), TRUE ~ Location))
us_units <- c("Atlanta","Sacramento","Denver","Baltimore","SanFrancisco","Boston","Los_Angeles","Syracuse","New_York",
              "Michigan","California","Washington_State","New_Hampshire","Massachusetts","NorthEast_US","Midwestern_US",
              "United_States","TEDDY_United_States")
perm_rho <- function(x, y, B = 20000) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 4) return(c(rho = NA, p = NA, n = length(x)))
  r0 <- cor(x, y, method = "spearman"); set.seed(1); rp <- replicate(B, cor(x, sample(y), method = "spearman"))
  c(rho = r0, p = (1 + sum(abs(rp) >= abs(r0) - 1e-12)) / (B + 1), n = length(x)) }
fmt_r <- function(r) sprintf("%.2f (p = %.2g, n = %d)", r[["rho"]], r[["p"]], as.integer(r[["n"]]))
ok_units <- c("ng/mL", "ng/ml", "ng/L", "pg/mL", "pg/ml", "ug/dL", "mg/L")
prep <- function(dd) dd %>% filter(is.na(Units) | Units %in% ok_units) %>%
  mutate(conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
         det = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
         val = ifelse(det, conc, lod / sqrt(2)), pid = coalesce(Mat_PID, Child_PID)) %>% filter(!is.na(val), val > 0)
attach_loc <- function(dd) bind_rows(
  dd %>% filter(!Study %in% c("TEDDY_Targeted", "ZIP_Targeted")) %>% left_join(loc_map %>% distinct(Study, loc_id), by = "Study"),
  dd %>% filter(Study == "TEDDY_Targeted") %>% left_join(loc_map %>% filter(Study == "TEDDY_Targeted") %>% distinct(Child_PID, loc_id), by = "Child_PID"),
  dd %>% filter(Study == "ZIP_Targeted") %>% left_join(loc_map %>% filter(Study == "ZIP_Targeted", !is.na(Mat_PID)) %>% distinct(Mat_PID, loc_id), by = "Mat_PID")
) %>% filter(!is.na(loc_id))
loc_medians <- function(dd) dd %>% group_by(Analyte_Code, loc_id, Study, pid) %>%
  summarise(v = median(val), det = mean(det), basis = first(Dilution_Adjustment), .groups = "drop") %>%
  group_by(Analyte_Code, loc_id) %>%
  summarise(median = median(v), df = mean(det), n = n(), studies = paste(sort(unique(std_study(Study))), collapse = "+"),
            basis = paste(sort(unique(basis)), collapse = "/"), .groups = "drop") %>% filter(n >= 20)

# ============================ STEP 1a: MEBCO ==================================
# Urinary As (unadjusted, ug/L) and creatinine (mg/dL) from the raw
# MEBCO_Targeted table (all 193 participants); water As from MEBCO_Epi.
# Urinary As three ways: unadjusted, creatinine-divided (= canonical ug/g,
# As / CRE * 100; verified identical ranks), and
#   unadjusted with log creatinine as a covariate (Barr et al. 2005) --
#   creatinine division is known to distort As in undernourished Bangladeshi
#   women because creatinine itself varies with nutrition / As methylation.
m <- new.env(); load("Input/MEBCO_34_Data_All.RData", envir = m); ML <- m[[ls(m)[1]]]
mb <- ML$MEBCO_Targeted %>% transmute(pid = as.character(PID), uAs_raw = As, CRE) %>%
  inner_join(ML$MEBCO_Epi %>% transmute(pid = as.character(PID), wAs = wateras_c1) %>% distinct(pid, .keep_all = TRUE), by = "pid") %>%
  filter(!is.na(wAs), !is.na(uAs_raw), uAs_raw > 0) %>% distinct(pid, .keep_all = TRUE) %>%
  mutate(CRE = ifelse(!is.na(CRE) & CRE > 0, CRE, NA), uAs = uAs_raw / CRE * 100,
         w_cat = cut(wAs, c(-Inf, 10, 50, Inf), labels = c("<10 µg/L (WHO)", "10–50 µg/L", ">50 µg/L (BD standard)")))
one <- function(fm, label, dd = mb) { f <- lm(fm, dd); ci <- confint(f)["log10(wAs)", ]
  tibble(urine_measure = label, n = nobs(f), slope_log10 = coef(f)["log10(wAs)"], lo = ci[1], hi = ci[2],
         fold_per_10x_water = 10^coef(f)["log10(wAs)"], fold_lo = 10^ci[1], fold_hi = 10^ci[2], r2 = summary(f)$r.squared) }
mbc <- mb %>% filter(!is.na(CRE))
s1a <- bind_rows(one(log10(uAs_raw) ~ log10(wAs), "Unadjusted (µg/L), all with water + urine As"),
                 one(log10(uAs_raw) ~ log10(wAs), "Unadjusted (µg/L), creatinine subset", mbc),
                 one(log10(uAs) ~ log10(wAs), "Creatinine-divided (µg/g, canonical), creatinine subset", mbc),
                 one(log10(uAs_raw) ~ log10(wAs) + log10(CRE), "Unadjusted + log creatinine covariate, creatinine subset", mbc)) %>%
  mutate(spearman = c(cor(mb$uAs_raw, mb$wAs, method = "spearman"), cor(mbc$uAs_raw, mbc$wAs, method = "spearman"),
                      cor(mbc$uAs, mbc$wAs, method = "spearman"), NA),
         pct_over_WHO10 = mean(mb$wAs > 10), pct_over_BD50 = mean(mb$wAs > 50))
if (nrow(mbc) == nrow(mb)) s1a <- s1a %>% filter(urine_measure != "Unadjusted (µg/L), creatinine subset") %>%
  mutate(urine_measure = sub(", (all with water \\+ urine As|creatinine subset)$", "", urine_measure))
c1a <- mb %>% group_by(w_cat) %>% summarise(n = n(), median_water = median(wAs), median_urine_raw = median(uAs_raw),
                                            median_urine_ugg = median(uAs, na.rm = TRUE), median_creatinine = median(CRE, na.rm = TRUE), .groups = "drop")
write.csv(s1a, file.path(out, "1a_MEBCO_water_vs_urine_As.csv"), row.names = FALSE)
write.csv(c1a, file.path(out, "1a_MEBCO_by_water_category.csv"), row.names = FALSE)

# ============================ STEP 1b: groundwater reliance ===================
gw <- read.csv(file.path(cache, "as_gw/data/data/groundwater_utilization_by_country.csv"), check.names = FALSE)
names(gw) <- c("country", "rural", "urban")
cname <- c(Accra = "Ghana", Mexico_City = "Mexico", South_Africa = "South Africa", Ecuador = "Ecuador", Chile = "Chile",
           Bangladesh = "Bangladesh", ZIP_Brazil = "Brazil", ZIP_Guatemala = "Guatemala", ZIP_Nicaragua = "Nicaragua",
           TEDDY_Finland = "Finland", TEDDY_Germany = "Germany", TEDDY_Sweden = "Sweden",
           Suriname = "Suriname", Thailand = "Thailand", Uganda = "Uganda", Tobago = "Trinidad and Tobago", `ZIP_Puerto Rico` = "Puerto Rico")
rel <- tibble(loc_id = U$loc_id) %>% mutate(country = ifelse(loc_id %in% us_units, "United States", cname[loc_id])) %>%
  left_join(gw, by = "country") %>%
  mutate(rel_source = case_when(country == "United States" ~ "ASSUMED: USGS private domestic wells (~14%)",
                                !is.na(rural) ~ "Podgorski & Berg 2020 (JMP/AQUASTAT)", TRUE ~ "not available"),
         rural = ifelse(country == "United States", 0.14, rural), urban = ifelse(country == "United States", 0.14, urban))
pop <- rast(list.files(file.path(cache, "ghsl"), pattern = "E2010.*\\.tif$", full.names = TRUE)[1])
asr <- rast(list.files(file.path(cache, "as_gw"), pattern = "probability_arsenic_gt_10ppb\\.tif$", recursive = TRUE, full.names = TRUE)[1])
gw_exp <- bind_rows(lapply(seq_len(nrow(U)), function(i) {
  g <- vect(U[i, ]); r <- rel[i, ]
  pc <- mask(crop(pop, g, snap = "out"), g)
  ac <- resample(crop(asr, ext(pc) + 0.02, snap = "out"), pc, method = "bilinear")
  dens <- pc / cellSize(pc, unit = "km")
  pv <- values(pc, mat = FALSE); av <- values(ac, mat = FALSE); dv <- values(dens, mat = FALSE)
  ok <- is.finite(pv) & is.finite(av) & pv > 0
  urb <- dv[ok] > 300
  relv <- ifelse(urb, r$urban, r$rural)
  tibble(loc_id = r$loc_id, P_As_popw = sum(pv[ok] * av[ok]) / sum(pv[ok]),
         urban_pop_share = sum(pv[ok][urb]) / sum(pv[ok]),
         gw_reliance_popw = if (is.na(r$rural)) NA_real_ else sum(pv[ok] * relv) / sum(pv[ok]),
         P_drink_As_gt10 = if (is.na(r$rural)) NA_real_ else sum(pv[ok] * av[ok] * relv) / sum(pv[ok]),
         rel_source = r$rel_source)
}))
write.csv(gw_exp, file.path(out, "1b_groundwater_As_exposure_by_location.csv"), row.names = FALSE)

as_med <- TD %>% filter(Analyte_Code == "As", Matrix_Class == "Urine") %>% prep() %>% attach_loc() %>% loc_medians()
as_med_adj <- TD %>% filter(Analyte_Code == "As", Matrix_Class == "Urine", Dilution_Adjustment != "none") %>% prep() %>% attach_loc() %>% loc_medians()
z1 <- as_med %>% inner_join(gw_exp, by = "loc_id"); z1a <- as_med_adj %>% inner_join(gw_exp, by = "loc_id")
r1b <- tribble(~layer, ~set, ~res,
  "P(GW As > 10) only", "all", fmt_r(perm_rho(z1$P_As_popw, z1$median)),
  "P(drinking GW As > 10) = P x reliance", "all", fmt_r(perm_rho(z1$P_drink_As_gt10, z1$median)),
  "P(GW As > 10) only", "dilution-adjusted", fmt_r(perm_rho(z1a$P_As_popw, z1a$median)),
  "P(drinking GW As > 10) = P x reliance", "dilution-adjusted", fmt_r(perm_rho(z1a$P_drink_As_gt10, z1a$median)),
  "P x reliance, US units excluded (assumed reliance)", "all", fmt_r(perm_rho(z1$P_drink_As_gt10[!z1$loc_id %in% us_units], z1$median[!z1$loc_id %in% us_units])),
  "P x reliance, Bangladesh excluded", "all", fmt_r(perm_rho(z1$P_drink_As_gt10[z1$loc_id != "Bangladesh"], z1$median[z1$loc_id != "Bangladesh"])))
write.csv(r1b, file.path(out, "1b_As_location_correlations.csv"), row.names = FALSE)
write.csv(z1, file.path(out, "1b_As_location_medians_with_exposure.csv"), row.names = FALSE)

# ============================ STEP 2: UCMR3 PFAS =============================
uc <- data.table::fread(file.path(cache, "ucmr3/UCMR3_All.txt"), sep = "\t", quote = "", header = TRUE, colClasses = "character", encoding = "Latin-1", data.table = FALSE) %>%  # fread: read.delim + latin1 read only the header in a non-UTF-8 locale
  filter(Contaminant %in% c("PFOS", "PFOA", "PFHxS", "PFNA", "PFHpA", "PFBS"))
pws <- uc %>% group_by(State, PWSID, Size) %>%
  summarise(det_PFOS_PFOA = any(AnalyticalResultsSign == "=" & Contaminant %in% c("PFOS", "PFOA")),
            det_PFOS = any(AnalyticalResultsSign == "=" & Contaminant == "PFOS"),
            det_PFOA = any(AnalyticalResultsSign == "=" & Contaminant == "PFOA"),
            det_any6 = any(AnalyticalResultsSign == "="), .groups = "drop")
st <- bind_rows(pws %>% mutate(grp = State), pws %>% mutate(grp = "US")) %>% group_by(grp) %>%
  summarise(n_pws = n(), n_large = sum(Size == "L"), pct_PFOS_PFOA = 100 * mean(det_PFOS_PFOA), pct_PFOS = 100 * mean(det_PFOS),
            pct_PFOA = 100 * mean(det_PFOA), pct_any6 = 100 * mean(det_any6),
            pct_PFOS_PFOA_large = 100 * mean(det_PFOS_PFOA[Size == "L"]), .groups = "drop")
write.csv(st, file.path(out, "2_UCMR3_state_detection.csv"), row.names = FALSE)
st_of <- list(Atlanta = "GA", Sacramento = "CA", Denver = "CO", Baltimore = "MD", SanFrancisco = "CA", Boston = "MA",
              Los_Angeles = "CA", Syracuse = "NY", New_York = "NY", Michigan = "MI", California = "CA", Washington_State = "WA",
              New_Hampshire = "NH", Massachusetts = "MA", NorthEast_US = c("CT", "MA", "RI", "NH", "VT", "ME"),
              Midwestern_US = c("IL", "IN", "IA", "KS", "MI", "MN", "MO", "NE", "ND", "OH", "SD", "WI"),
              United_States = "US", TEDDY_United_States = c("CO", "GA", "FL", "WA"))
uc_loc <- bind_rows(lapply(names(st_of), function(l) st %>% filter(grp %in% st_of[[l]]) %>%
  summarise(across(starts_with("pct_"), mean)) %>% mutate(loc_id = l, states = paste(st_of[[l]], collapse = ","))))
write.csv(uc_loc, file.path(out, "2_UCMR3_by_location.csv"), row.names = FALSE)

pf <- TD %>% filter(Analyte_Code %in% c("PFOS", "PFOA"), Matrix_Class == "Blood") %>% prep() %>% attach_loc() %>%
  filter(loc_id %in% us_units)
pf_med <- loc_medians(pf) %>% left_join(pf %>% group_by(Analyte_Code, loc_id) %>%
  summarise(year_mid = median(Sample_Collection_Year, na.rm = TRUE), year_min = min(Sample_Collection_Year, na.rm = TRUE),
            year_max = max(Sample_Collection_Year, na.rm = TRUE), .groups = "drop"), by = c("Analyte_Code", "loc_id")) %>%
  left_join(uc_loc, by = "loc_id")
write.csv(pf_med, file.path(out, "2_PFAS_location_medians_UCMR3.csv"), row.names = FALSE)
# Exclude NYUWHS (1985-90): sampled ~25 years before UCMR3 and before the PFOS phase-out.
pf_med2 <- pf_med %>% filter(is.finite(year_max), year_max >= 2000)  # also drops ATLAA (no sampling years)
r2 <- bind_rows(lapply(c("PFOS", "PFOA"), function(a) { z <- pf_med2 %>% filter(Analyte_Code == a)
  tibble(analyte = a, metric = c("% PWS with PFOS/PFOA detected", "% PWS with own-analyte detected", "% large PWS with PFOS/PFOA detected"),
         res = c(fmt_r(perm_rho(z$pct_PFOS_PFOA, z$median)), fmt_r(perm_rho(z[[paste0("pct_", a)]], z$median)),
                 fmt_r(perm_rho(z$pct_PFOS_PFOA_large, z$median)))) }))
# individual level, adjusted for sampling year (PFAS secular decline) and cohort-era
pfi <- pf %>% filter(!is.na(Sample_Collection_Year), Sample_Collection_Year >= 2000) %>% left_join(uc_loc %>% select(loc_id, pct_PFOS_PFOA), by = "loc_id") %>%
  mutate(lval = log10(val), yr = Sample_Collection_Year - 2010, uc10 = pct_PFOS_PFOA / 10)
m2 <- bind_rows(lapply(c("PFOS", "PFOA"), function(a) {
  dd <- pfi %>% filter(Analyte_Code == a)
  fm <- lmer(lval ~ uc10 + yr + (1 | loc_id) + (1 | pid), dd, control = lmerControl(calc.derivs = FALSE))
  cf <- summary(fm)$coefficients; ci <- suppressMessages(confint(fm, method = "Wald"))
  tibble(analyte = a, term = c("UCMR3 % PWS detected (per +10 pct-pts)", "Sampling year (per year)"),
         fold = 10^cf[c("uc10", "yr"), 1], lo = 10^ci[c("uc10", "yr"), 1], hi = 10^ci[c("uc10", "yr"), 2],
         n_participants = n_distinct(dd$pid), n_locations = n_distinct(dd$loc_id))
}))
write.csv(r2, file.path(out, "2_PFAS_location_correlations.csv"), row.names = FALSE)
write.csv(m2, file.path(out, "2_PFAS_individual_model.csv"), row.names = FALSE)

# ============================ STEP 3: ZIP water source ========================
z <- new.env(); load("Input/ZIP_1839_Data_All.RData", envir = z); ZE <- z[[ls(z)[1]]]$ZIP_EPI_1
ws <- ZE %>% transmute(Mat_PID = as.character(MOTHER_PID), country = tri_country,
  water = case_when(tri_prim_water_src == 1 ~ "Water utility", tri_prim_water_src == 2 ~ "Private well",
                    tri_prim_water_src == 7 ~ "Bottled/delivered", tri_prim_water_src %in% 3:6 ~ "Other (rain/surface/tank/community well)")) %>%
  filter(!is.na(Mat_PID)) %>% distinct(Mat_PID, .keep_all = TRUE)
write.csv(ws %>% count(country, water) %>% pivot_wider(names_from = water, values_from = n, values_fill = 0),
          file.path(out, "3_ZIP_water_source_by_country.csv"), row.names = FALSE)
zp <- TD %>% filter(Study == "ZIP_Targeted") %>% prep() %>% inner_join(ws, by = "Mat_PID") %>%
  filter(!is.na(country), !is.na(water), water != "Other (rain/surface/tank/community well)") %>%
  mutate(lval = log10(val), water = relevel(factor(water), "Water utility"), country = relevel(factor(country), "Puerto Rico"))
zp_an <- zp %>% group_by(Analyte_Code) %>% summarise(df = mean(det), .groups = "drop") %>% filter(df >= 0.5) %>% pull(Analyte_Code)
m3 <- bind_rows(lapply(zp_an, function(a) {
  dd <- zp %>% filter(Analyte_Code == a)
  f0 <- lmer(lval ~ country + (1 | Mat_PID), dd, REML = FALSE, control = lmerControl(calc.derivs = FALSE))
  f1 <- lmer(lval ~ country + water + (1 | Mat_PID), dd, REML = FALSE, control = lmerControl(calc.derivs = FALSE))
  c0 <- fixef(f0); c1 <- fixef(f1); ci1 <- suppressMessages(confint(f1, method = "Wald"))
  ct <- grep("^country", names(c0), value = TRUE); wt <- grep("^water", names(c1), value = TRUE)
  bind_rows(
    tibble(analyte = a, term = sub("water", "", wt), type = "water source vs utility (within country)",
           fold = 10^c1[wt], lo = 10^ci1[wt, 1], hi = 10^ci1[wt, 2]),
    tibble(analyte = a, term = sub("country", "", ct), type = "country vs Puerto Rico",
           fold = 10^c0[ct], fold_adj_water = 10^c1[ct], pct_attenuation = 100 * (1 - abs(c1[ct]) / abs(c0[ct]))),
    tibble(analyte = a, term = "LRT: water source adds to country", type = "test", p = anova(f0, f1)$`Pr(>Chisq)`[2])) %>%
    mutate(n_participants = n_distinct(dd$Mat_PID), detection = mean(dd$det))
}))
write.csv(m3, file.path(out, "3_ZIP_water_source_models.csv"), row.names = FALSE)

# ============================ print + figure =================================
cat("\n== 1a MEBCO ==\n"); print(s1a, width = Inf); print(c1a)
cat("\n== 1b groundwater As exposure ==\n"); print(z1 %>% select(loc_id, median, P_As_popw, gw_reliance_popw, P_drink_As_gt10, rel_source) %>% arrange(desc(median)) %>% as.data.frame(), digits = 3); print(r1b)
cat("\n== 2 UCMR3 ==\n"); print(uc_loc %>% select(loc_id, states, pct_PFOS_PFOA, pct_PFOS_PFOA_large) %>% as.data.frame(), digits = 3)
print(pf_med %>% select(Analyte_Code, loc_id, studies, median, year_min, year_max, pct_PFOS_PFOA) %>% as.data.frame(), digits = 3); print(r2); print(m2)
cat("\n== 3 ZIP ==\n"); print(read.csv(file.path(out, "3_ZIP_water_source_by_country.csv"))); print(m3 %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())

grp <- function(l) case_when(l %in% us_units ~ "US", grepl("^TEDDY_(Finland|Sweden|Germany)", l) ~ "Europe", TRUE ~ "Latin America, Africa & Asia")
cols <- c(US = "#3b6fb6", Europe = "#2e8b57", `Latin America, Africa & Asia` = "#c0392b")
pA <- ggplot(mb, aes(wAs, uAs_raw)) + geom_point(alpha = 0.5, size = 1.4, colour = "#c0392b") +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey20") + geom_vline(xintercept = c(10, 50), linetype = 2, colour = "grey50") +
  scale_x_log10() + scale_y_log10() +
  annotate("text", x = 0.3, y = max(mb$uAs_raw), hjust = 0, vjust = 1, size = 2.6,
           label = sprintf("Unadjusted: ρ = %.2f, ×%.2f per 10× water (n = %d)\nCreatinine-divided (canonical µg/g, n = %d): ρ = %.2f\nCreatinine as covariate: ×%.2f (95%% CI %.2f–%.2f)",
                           s1a$spearman[1], s1a$fold_per_10x_water[1], s1a$n[1], s1a$n[nrow(s1a) - 1], s1a$spearman[nrow(s1a) - 1], s1a$fold_per_10x_water[nrow(s1a)], s1a$fold_lo[nrow(s1a)], s1a$fold_hi[nrow(s1a)])) +
  labs(title = "A. MEBCO (Bangladesh): individual drinking-water vs urinary As", x = "Household drinking-water As (µg/L, assumed)",
       y = "Urinary As, unadjusted (µg/L)", caption = "Dashed: WHO 10 µg/L and Bangladesh 50 µg/L standards") + theme_bw(base_size = 9)
zb <- z1 %>% select(loc_id, median, df, `P(groundwater As > 10)` = P_As_popw, `P(drinking groundwater As > 10)` = P_drink_As_gt10) %>%
  pivot_longer(starts_with("P("), names_to = "layer", values_to = "x") %>% filter(is.finite(x)) %>% mutate(group = grp(loc_id))
labB <- tibble(layer = c("P(groundwater As > 10)", "P(drinking groundwater As > 10)"), lab = paste("ρ =", sub(" \\(.*", "", r1b$res[1:2])))
pB <- ggplot(zb, aes(x + 1e-3, median)) + geom_point(aes(colour = group), size = 2.3) +
  ggrepel::geom_text_repel(aes(label = gsub("_", " ", loc_id)), size = 2.1, colour = "grey30", seed = 1) +
  geom_text(data = labB, aes(label = lab), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.7, inherit.aes = FALSE) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values = cols, name = NULL) + facet_wrap(~layer, scales = "free_x") +
  labs(title = "B. Location medians: groundwater As vs groundwater As actually drunk", x = "Population-weighted probability (log)",
       y = "Urinary As location median") + theme_bw(base_size = 9) + theme(legend.position = "bottom")
pC <- ggplot(pf_med2 %>% mutate(lab = paste0(gsub("_", " ", loc_id), ifelse(is.finite(year_min), paste0(" (", year_min, "–", year_max, ")"), ""))), aes(pct_PFOS_PFOA, median)) +
  geom_point(colour = "#3b6fb6", size = 2.3) + ggrepel::geom_text_repel(aes(label = lab), size = 2.1, colour = "grey30", seed = 1) +
  facet_wrap(~Analyte_Code, scales = "free_y") + scale_y_log10() +
  geom_text(data = r2 %>% filter(grepl("^% PWS with PFOS/PFOA", metric)) %>% transmute(Analyte_Code = analyte, lab = paste("ρ =", sub(" \\(.*", "", res))),
            aes(label = lab), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.4, size = 2.7, inherit.aes = FALSE) +
  labs(title = "C. US blood PFAS vs UCMR3 drinking-water detections (2013–15)", x = "% public water systems with PFOS or PFOA detected (state)",
       y = "Blood location median (ng/mL)") + theme_bw(base_size = 9)
pD <- ggplot(m3 %>% filter(type == "water source vs utility (within country)"), aes(fold, analyte, colour = term)) +
  geom_vline(xintercept = 1, linetype = 2) + geom_pointrange(aes(xmin = lo, xmax = hi), position = position_dodge(width = 0.5), size = 0.25) +
  scale_x_log10() + scale_colour_manual(values = c("#e08e0b", "#7d3c98"), name = "vs water utility") +
  labs(title = "D. ZIP: urinary pesticide metabolites by household water source (within country)", x = "Fold difference (log)", y = NULL) +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
g <- (pA | pB) / (pC | pD) + plot_layout(widths = c(1, 1.4))
ggsave(file.path(out, "SuppFig_Water_Determinants.png"), g, width = 13, height = 9.5, dpi = 200)
ggsave(file.path(out, "SuppFig_Water_Determinants.pdf"), g, width = 13, height = 9.5, device = cairo_pdf)

# ---- single panels, for reports that split cohort-data and linkage results ----
ggsave(file.path(out, "Panel_MEBCO_water_As.png"), pA, width = 6.5, height = 4.6, dpi = 200)
ggsave(file.path(out, "Panel_groundwater_As.png"), pB, width = 7.5, height = 4.6, dpi = 200)
ggsave(file.path(out, "Panel_UCMR3_PFAS.png"), pC, width = 7, height = 4.4, dpi = 200)
ggsave(file.path(out, "Panel_ZIP_water_source.png"), pD, width = 6.5, height = 4.4, dpi = 200)
