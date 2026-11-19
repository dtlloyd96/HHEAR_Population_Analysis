# R4_5_diet_pesticide_lead.R
# Reviewer 4, determinants of geography -- country-level, time-matched layers for
# three of R4's named determinants:
#   A. Agricultural pesticide use  -> urinary pesticide metabolites
#      FAOSTAT Inputs_Pesticides_Use: "Insecticides" agricultural use (tonnes a.i.)
#      / FAOSTAT Inputs_LandUse "Cropland" (1000 ha) = kg a.i. per ha cropland.
#      The organophosphate / pyrethroid breakdown is missing for the US, Mexico,
#      South Africa, Chile, Uganda and pre-2021 Brazil, so total insecticides are used.
#   B. Diet -> urinary arsenic (and urinary mercury, descriptive)
#      FAOSTAT Food Balance Sheets, "Food supply quantity (kg/capita/yr)":
#      Fish, Seafood (item 2960) and Rice and products (2805). Historic FBS
#      (1961-2013) and new FBS (2010-) are combined, new values preferred.
#      Also a two-determinant model with drinking-water As (R4_5_water.R).
#   C. Regulation -> urinary lead
#      Leaded-petrol ban year in road vehicles (Our World in Data, "Lead petrol
#      ban", compiled from UNEP and other sources, 2021); exposure = years since
#      ban at sampling (sample-weighted). Suriname has no date.
#   (US housing age was planned but the Census API now requires a key.)
# Assignment: US sub-national units -> USA; TEDDY sites and ZIP countries -> own
# country; Tobago -> Trinidad and Tobago. Time matching = sample-weighted mean
# over sampling years (location_study_year_weights.csv), nearest available year.
# Location medians: same rules as R4_5_geospatial_scatter.R (to_ng_per_mL, <LOD
# -> LOD/sqrt(2), participant median, >= 20 participants, >= 20% detected),
# all urine and dilution-adjusted-only. Spearman with permutation p (20,000).
# Writes Output/Reviewer_Round2/R4_5_determinants/diet_pesticide_lead/. Run from project root.

suppressMessages({library(dplyr); library(tidyr); library(readr); library(ggplot2); library(lme4); library(patchwork)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
source("Coding/Paper1_Submission2/harmonization_helpers.R")
cache <- Sys.getenv("HHEAR_GEO_CACHE", path.expand("~/HHEAR_geo_cache"))
geo <- "Output/Reviewer_Round2/R4_5_determinants/geospatial"
out <- "Output/Reviewer_Round2/R4_5_determinants/diet_pesticide_lead"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
yw <- read.csv(file.path(geo, "location_study_year_weights.csv"))
L  <- read.csv(file.path(geo, "location_layers_time_matched.csv"))

us_units <- c("Atlanta","Sacramento","Denver","Baltimore","SanFrancisco","Boston","Los_Angeles","Syracuse","New_York",
              "Michigan","California","Washington_State","New_Hampshire","Massachusetts","NorthEast_US","Midwestern_US",
              "United_States","TEDDY_United_States")
fao_name <- c(Accra = "Ghana", Mexico_City = "Mexico", Tobago = "Trinidad and Tobago", Suriname = "Suriname",
              South_Africa = "South Africa", Ecuador = "Ecuador", Chile = "Chile", Uganda = "Uganda", Bangladesh = "Bangladesh",
              Thailand = "Thailand", ZIP_Brazil = "Brazil", ZIP_Guatemala = "Guatemala", ZIP_Nicaragua = "Nicaragua",
              `ZIP_Puerto Rico` = "Puerto Rico", TEDDY_Finland = "Finland", TEDDY_Germany = "Germany", TEDDY_Sweden = "Sweden")
iso3 <- c(Ghana = "GHA", Mexico = "MEX", `Trinidad and Tobago` = "TTO", Suriname = "SUR", `South Africa` = "ZAF", Ecuador = "ECU",
          Chile = "CHL", Uganda = "UGA", Bangladesh = "BGD", Thailand = "THA", Brazil = "BRA", Guatemala = "GTM", Nicaragua = "NIC",
          `Puerto Rico` = "PRI", Finland = "FIN", Germany = "DEU", Sweden = "SWE", `United States of America` = "USA")
loc_country <- tibble(loc_id = L$loc_id, country = ifelse(L$loc_id %in% us_units, "United States of America", fao_name[L$loc_id]))

wide_to_long <- function(df) df %>% pivot_longer(matches("^Y[0-9]{4}$"), names_to = "year", values_to = "value") %>%
  mutate(year = as.integer(sub("Y", "", year))) %>% filter(!is.na(value))
time_match <- function(series, name) {   # series: country, year, value
  yrs <- yw %>% group_by(loc_id, year) %>% summarise(n = sum(n), .groups = "drop") %>% left_join(loc_country, by = "loc_id")
  bind_rows(lapply(split(yrs, yrs$loc_id), function(z) {
    s <- series %>% filter(country == z$country[1])
    if (!nrow(s)) return(tibble(loc_id = z$loc_id[1], value = NA_real_, mean_gap = NA_real_))
    near <- sapply(z$year, function(y) s$year[which.min(abs(s$year - y))])
    v <- s$value[match(near, s$year)]
    tibble(loc_id = z$loc_id[1], value = weighted.mean(v, z$n), mean_gap = weighted.mean(abs(near - z$year), z$n))
  })) %>% rename(!!name := value, !!paste0(name, "_year_gap") := mean_gap)
}

# ---- layers -----------------------------------------------------------------
# Unzipped FAOSTAT folders may be deleted to save space; re-extract from the zips if needed.
for (z in c("Inputs_Pesticides_Use_E_All_Data", "Inputs_LandUse_E_All_Data", "FoodBalanceSheets_E_All_Data", "FoodBalanceSheetsHistoric_E_All_Data")) {
  if (!dir.exists(file.path(cache, "faostat", z))) unzip(file.path(cache, "faostat", paste0(z, ".zip")), exdir = file.path(cache, "faostat", z))
}
P <- read_csv(file.path(cache, "faostat/Inputs_Pesticides_Use_E_All_Data/Inputs_Pesticides_Use_E_All_Data_NOFLAG.csv"), show_col_types = FALSE)
LU <- read_csv(file.path(cache, "faostat/Inputs_LandUse_E_All_Data/Inputs_LandUse_E_All_Data_NOFLAG.csv"), show_col_types = FALSE)
ins <- P %>% filter(`Item Code` == 1309, Element == "Agricultural Use") %>% wide_to_long() %>% transmute(country = Area, year, t = value)
crop <- LU %>% filter(Item == "Cropland", Element == "Area") %>% wide_to_long() %>% transmute(country = Area, year, kha = value, unit = Unit)
stopifnot(all(crop$unit == "1000 ha"))
insha <- ins %>% inner_join(crop, by = c("country", "year")) %>% transmute(country, year, value = t / kha)  # t per 1000 ha = kg/ha
tot <- P %>% filter(`Item Code` == 1357, Element == "Use per area of cropland") %>% wide_to_long() %>% transmute(country = Area, year, value)

FB_new <- read_csv(file.path(cache, "faostat/FoodBalanceSheets_E_All_Data/FoodBalanceSheets_E_All_Data_NOFLAG.csv"), show_col_types = FALSE)
FB_old <- read_csv(file.path(cache, "faostat/FoodBalanceSheetsHistoric_E_All_Data/FoodBalanceSheetsHistoric_E_All_Data_NOFLAG.csv"), show_col_types = FALSE)
fbs <- function(code) {
  nw <- FB_new %>% filter(`Item Code` == code, `Element Code` == 645) %>% wide_to_long() %>% transmute(country = Area, year, value, src = "new")
  od <- FB_old %>% filter(`Item Code` == code, `Element Code` == 645) %>% wide_to_long() %>% transmute(country = Area, year, value, src = "historic")
  bind_rows(nw, od %>% anti_join(nw, by = c("country", "year"))) %>% select(country, year, value)
}
fish <- fbs(2960); rice <- fbs(2805)

ow <- read_csv(file.path(cache, "lead/lead-petrol-ban.csv"), show_col_types = FALSE)
names(ow)[4] <- "banned"
ban <- ow %>% filter(banned == "Yes") %>% group_by(Code) %>% summarise(ban_year = min(Year), .groups = "drop")
yrs_since <- yw %>% group_by(loc_id, year) %>% summarise(n = sum(n), .groups = "drop") %>% left_join(loc_country, by = "loc_id") %>%
  mutate(Code = iso3[country]) %>% left_join(ban, by = "Code") %>% group_by(loc_id) %>%
  summarise(lead_ban_year = first(ban_year), years_since_lead_ban = if (all(is.na(ban_year))) NA_real_ else weighted.mean(year - ban_year, n), .groups = "drop")

LY <- loc_country %>% left_join(time_match(insha, "insecticide_kg_ha"), by = "loc_id") %>%
  left_join(time_match(tot, "pesticide_total_kg_ha"), by = "loc_id") %>%
  left_join(time_match(fish, "fish_kg_cap"), by = "loc_id") %>% left_join(time_match(rice, "rice_kg_cap"), by = "loc_id") %>%
  left_join(yrs_since, by = "loc_id")
write.csv(LY, file.path(out, "layers_by_location.csv"), row.names = FALSE)
# also append to the shared location layer table
Lall <- L %>% select(-any_of(setdiff(names(LY), "loc_id"))) %>% left_join(LY, by = "loc_id")
write.csv(Lall, file.path(geo, "location_layers_time_matched.csv"), row.names = FALSE)

# ---- biomarker location medians ---------------------------------------------
e <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = e); TD <- e$Targeted_Data_All
s <- new.env(); load("Input/Harmonized_Datasets/Harmonized_SDOH.RData", envir = s)
loc_map <- s$SDOH_Data_All %>% distinct(Study, Child_PID, Mat_PID, Location, Country) %>%
  mutate(Study = sub("_SDOH$", "_Targeted", Study),
         loc_id = case_when(Study == "TEDDY_Targeted" ~ paste0("TEDDY_", Location),
                            Study == "ZIP_Targeted" & !is.na(Country) ~ paste0("ZIP_", Country), TRUE ~ Location))
ok_units <- c("ng/mL", "ng/ml", "ng/L", "pg/mL", "pg/ml", "ug/dL", "mg/L")
PEST <- c("TCP", "PBA", "FPBA", "PNP", "DMP", "DEP", "DMTP", "DETP")
ANALYTES <- c(PEST, "As", "Hg", "Pb")
d <- TD %>% filter(Analyte_Code %in% ANALYTES, Matrix_Class == "Urine", is.na(Units) | Units %in% ok_units) %>%
  mutate(conc = to_ng_per_mL(Concentration, Units), lod = to_ng_per_mL(LOD, Units),
         det = ifelse(is.na(Detected), !is.na(conc) & (is.na(lod) | conc >= lod), Detected),
         val = ifelse(det, conc, lod / sqrt(2)), pid = paste(Study, coalesce(Mat_PID, Child_PID))) %>% filter(!is.na(val), val > 0)
d <- bind_rows(
  d %>% filter(!Study %in% c("TEDDY_Targeted", "ZIP_Targeted")) %>% left_join(loc_map %>% distinct(Study, loc_id), by = "Study"),
  d %>% filter(Study == "TEDDY_Targeted") %>% left_join(loc_map %>% filter(Study == "TEDDY_Targeted") %>% distinct(Child_PID, loc_id), by = "Child_PID"),
  d %>% filter(Study == "ZIP_Targeted") %>% left_join(loc_map %>% filter(Study == "ZIP_Targeted", !is.na(Mat_PID)) %>% distinct(Mat_PID, loc_id), by = "Mat_PID")
) %>% filter(!is.na(loc_id)) %>% mutate(adjusted = Dilution_Adjustment != "none")
d <- d %>% group_by(loc_id, Analyte_Code) %>% filter(mean(det) >= 0.2) %>% ungroup()
med <- function(dd) dd %>% group_by(Analyte_Code, loc_id, pid) %>% summarise(v = median(val), det = mean(det), .groups = "drop") %>%
  group_by(Analyte_Code, loc_id) %>% summarise(median = median(v), df = mean(det), n = n(), .groups = "drop") %>% filter(n >= 20)
M <- bind_rows(med(d) %>% mutate(set = "all"), med(d %>% filter(adjusted)) %>% mutate(set = "dilution_adjusted_only")) %>%
  left_join(d %>% group_by(loc_id) %>% summarise(studies = paste(sort(unique(std_study(Study))), collapse = "+")), by = "loc_id") %>%
  left_join(Lall %>% select(loc_id, insecticide_kg_ha, pesticide_total_kg_ha, fish_kg_cap, rice_kg_cap, years_since_lead_ban, P_As = GW_As_prob_gt10), by = "loc_id")
gwx <- read.csv("Output/Reviewer_Round2/R4_5_determinants/water/1b_groundwater_As_exposure_by_location.csv") %>% select(loc_id, P_drink_As_gt10)
M <- M %>% left_join(gwx, by = "loc_id")
write.csv(M, file.path(out, "location_medians_with_layers.csv"), row.names = FALSE)

perm_rho <- function(x, y, B = 20000) { ok <- is.finite(x) & is.finite(y); x <- x[ok]; y <- y[ok]
  if (length(x) < 4) return(c(rho = NA, p = NA, n = length(x)))
  r0 <- cor(x, y, method = "spearman"); set.seed(1); rp <- replicate(B, cor(x, sample(y), method = "spearman"))
  c(rho = r0, p = (1 + sum(abs(rp) >= abs(r0) - 1e-12)) / (B + 1), n = length(x)) }
pairs <- bind_rows(
  expand_grid(analyte = PEST, layer = c("insecticide_kg_ha", "pesticide_total_kg_ha"), block = "A. Pesticide use"),
  expand_grid(analyte = "As", layer = c("fish_kg_cap", "rice_kg_cap"), block = "B. Diet"),
  expand_grid(analyte = "Hg", layer = c("fish_kg_cap"), block = "B. Diet"),
  expand_grid(analyte = "Pb", layer = c("years_since_lead_ban"), block = "C. Lead regulation"))
R <- bind_rows(lapply(seq_len(nrow(pairs)), function(i) { p <- pairs[i, ]
  bind_rows(lapply(c("all", "dilution_adjusted_only"), function(st) {
    z <- M %>% filter(Analyte_Code == p$analyte, set == st)
    r <- perm_rho(z[[p$layer]], z$median)
    tibble(block = p$block, analyte = p$analyte, layer = p$layer, set = st, n_locations = r[["n"]], rho = r[["rho"]], perm_p = r[["p"]],
           locations = paste(z$loc_id[is.finite(z[[p$layer]])], collapse = ";"))
  })) }))
# arsenic beyond Bangladesh (water-dominated) and with both determinants
zAs <- M %>% filter(Analyte_Code == "As", set == "all")
R <- bind_rows(R,
  tibble(block = "B. Diet", analyte = "As", layer = "fish_kg_cap (Bangladesh excluded)", set = "all",
         n_locations = sum(zAs$loc_id != "Bangladesh" & is.finite(zAs$fish_kg_cap)),
         rho = perm_rho(zAs$fish_kg_cap[zAs$loc_id != "Bangladesh"], zAs$median[zAs$loc_id != "Bangladesh"])[["rho"]],
         perm_p = perm_rho(zAs$fish_kg_cap[zAs$loc_id != "Bangladesh"], zAs$median[zAs$loc_id != "Bangladesh"])[["p"]]))
# robustness: drop the most influential anchors
rob <- function(a, lay, drop, label) { z <- M %>% filter(Analyte_Code == a, set == "all", !loc_id %in% drop); r <- perm_rho(z[[lay]], z$median)
  tibble(block = "robustness", analyte = a, layer = paste0(lay, " (", label, ")"), set = "all", n_locations = r[["n"]], rho = r[["rho"]], perm_p = r[["p"]]) }
R <- bind_rows(R,
  rob("PBA", "insecticide_kg_ha", c("TEDDY_Finland", "TEDDY_Sweden"), "Finland & Sweden excluded"),
  rob("PBA", "insecticide_kg_ha", c("TEDDY_Finland", "TEDDY_Sweden", "TEDDY_Germany", "TEDDY_United_States"), "TEDDY excluded"),
  rob("Pb", "years_since_lead_ban", "Bangladesh", "Bangladesh excluded"),
  rob("TCP", "insecticide_kg_ha", "Baltimore", "Baltimore excluded"))
write.csv(R, file.path(out, "location_correlations.csv"), row.names = FALSE)

# arsenic: two determinants (location level, log-linear), descriptive with n ~ 10
fitAs <- lm(log10(median) ~ log10(fish_kg_cap) + log10(P_drink_As_gt10), data = zAs %>% filter(is.finite(P_drink_As_gt10), is.finite(fish_kg_cap)))
asm <- broom_like <- data.frame(term = names(coef(fitAs)), estimate = coef(fitAs), confint(fitAs), n = nobs(fitAs), r2 = summary(fitAs)$r.squared)
write.csv(asm, file.path(out, "As_two_determinant_model.csv"), row.names = FALSE)

# individual-level between-location models for the pesticide pairs with >= 5 locations
ind <- bind_rows(lapply(c(PEST, "Pb"), function(a) {
  lay <- if (a == "Pb") "years_since_lead_ban" else "insecticide_kg_ha"
  dd <- d %>% filter(Analyte_Code == a) %>% left_join(Lall %>% select(loc_id, x = all_of(lay)), by = "loc_id") %>% filter(is.finite(x)) %>%
    mutate(lval = log10(val), xs = if (a == "Pb") x / 10 else log2(x))
  if (n_distinct(dd$loc_id) < 5) return(NULL)
  if (a == "Pb") dd_noB <- dd %>% filter(loc_id != "Bangladesh")
  fml <- if (n_distinct(dd$Dilution_Adjustment) > 1) lval ~ xs + Dilution_Adjustment + (1 | loc_id) + (1 | pid) else lval ~ xs + (1 | loc_id) + (1 | pid)
  m <- lmer(fml, dd, control = lmerControl(calc.derivs = FALSE)); ci <- suppressMessages(confint(m, method = "Wald"))["xs", ]
  res <- tibble(analyte = a, layer = lay, unit = ifelse(a == "Pb", "per 10 years since ban", "per doubling of insecticide kg/ha"),
         fold = 10^fixef(m)["xs"], lo = 10^ci[1], hi = 10^ci[2], n_locations = n_distinct(dd$loc_id), n_participants = n_distinct(dd$pid), subset = "all")
  if (a == "Pb") { m2 <- lmer(fml, dd_noB, control = lmerControl(calc.derivs = FALSE)); c2 <- suppressMessages(confint(m2, method = "Wald"))["xs", ]
    res <- bind_rows(res, tibble(analyte = a, layer = lay, unit = "per 10 years since ban", fold = 10^fixef(m2)["xs"], lo = 10^c2[1], hi = 10^c2[2],
                                 n_locations = n_distinct(dd_noB$loc_id), n_participants = n_distinct(dd_noB$pid), subset = "Bangladesh excluded")) }
  if (a == "PBA") { dn <- dd %>% filter(!loc_id %in% c("TEDDY_Finland", "TEDDY_Sweden")); m3 <- lmer(fml, dn, control = lmerControl(calc.derivs = FALSE))
    c3 <- suppressMessages(confint(m3, method = "Wald"))["xs", ]
    res <- bind_rows(res, tibble(analyte = a, layer = lay, unit = "per doubling of insecticide kg/ha", fold = 10^fixef(m3)["xs"], lo = 10^c3[1], hi = 10^c3[2],
                                 n_locations = n_distinct(dn$loc_id), n_participants = n_distinct(dn$pid), subset = "Finland & Sweden excluded")) }
  res
}))
write.csv(ind, file.path(out, "individual_between_location_models.csv"), row.names = FALSE)

print(LY %>% filter(!is.na(country)) %>% mutate(across(where(is.numeric), ~ round(.x, 2))) %>% as.data.frame())
print(R %>% select(-locations) %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())
print(asm); print(ind %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())

# ---- figure ------------------------------------------------------------------
grp <- function(l) case_when(l %in% us_units ~ "US", grepl("^TEDDY_(Finland|Sweden|Germany)", l) ~ "Europe", TRUE ~ "Latin America, Africa & Asia")
cols <- c(US = "#3b6fb6", Europe = "#2e8b57", `Latin America, Africa & Asia` = "#c0392b")
panel <- function(a, lay, xlab, title, logx = TRUE, set = "all") {
  z <- M %>% filter(Analyte_Code == a, set == !!set, is.finite(.data[[lay]])) %>%
    mutate(group = factor(grp(loc_id), levels = names(cols)), lowdf = factor(df < 0.5, levels = c(FALSE, TRUE)))
  r <- R %>% filter(analyte == a, layer == lay, set == !!set)
  g <- ggplot(z, aes(.data[[lay]], median)) + geom_point(aes(colour = group, shape = lowdf), size = 2.3) +
    ggrepel::geom_text_repel(aes(label = gsub("_", " ", loc_id)), size = 2, colour = "grey30", seed = 1, max.overlaps = 30) +
    scale_y_log10() + scale_colour_manual(values = cols, name = NULL, drop = FALSE) +
    scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("≥50% detected", "<50% detected"), name = NULL, drop = FALSE) +
    labs(title = title, subtitle = sprintf("ρ = %.2f (perm p = %.2g, n = %d locations)", r$rho, r$perm_p, r$n_locations), x = xlab, y = paste(a, "location median")) + theme_bw(base_size = 8.5) + theme(legend.position = "bottom")
  if (logx) g + scale_x_log10() else g
}
g <- (panel("TCP", "insecticide_kg_ha", "Insecticide use (kg a.i./ha cropland)", "A. TCPy (chlorpyrifos) vs national insecticide use") |
      panel("PBA", "insecticide_kg_ha", "Insecticide use (kg a.i./ha cropland)", "B. 3-PBA (pyrethroids) vs insecticide use") |
      panel("DMP", "insecticide_kg_ha", "Insecticide use (kg a.i./ha cropland)", "C. DMP (OP) vs insecticide use")) /
     (panel("As", "fish_kg_cap", "Fish & seafood supply (kg/person/yr)", "D. Urinary As vs fish & seafood supply") |
      panel("As", "rice_kg_cap", "Rice supply (kg/person/yr)", "E. Urinary As vs rice supply") |
      panel("Pb", "years_since_lead_ban", "Years since leaded-petrol ban, at sampling", "F. Urinary Pb vs time since leaded-petrol ban", logx = FALSE)) +
  plot_layout(guides = "collect") & theme(legend.position = "bottom")
ggsave(file.path(out, "SuppFig_Diet_Pesticide_Lead.png"), g, width = 13, height = 9, dpi = 200)
ggsave(file.path(out, "SuppFig_Diet_Pesticide_Lead.pdf"), g, width = 13, height = 9, device = cairo_pdf)
