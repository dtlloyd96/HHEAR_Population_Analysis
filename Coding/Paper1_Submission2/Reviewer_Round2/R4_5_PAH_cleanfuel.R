# R4_5_PAH_cleanfuel.R
# Reviewer 4, determinants of geography: household solid-fuel use as the
# mechanism behind the between-location PAH ~ PM2.5 association.
#
# Data: WHO Household Energy Database estimates of "access to clean fuels and
# technologies for cooking (% of population)", SDG 7.1.2, retrieved from the
# World Bank WDI API (series EG.CFT.ACCS.ZS national, .UR.ZS urban, .RU.ZS
# rural; source = WHO) because the WHO GHO API was down (503) on 2026-11-04.
# Raw JSON cached in ~/HHEAR_geo_cache/cleanfuel/. Coverage 2000-2023.
# Exposure used: % relying on POLLUTING fuels = 100 - clean.
#
# Assignment: non-US city units (Accra, Mexico City) take the URBAN series;
# all other units take the NATIONAL series (US sub-national units -> USA;
# Tobago -> Trinidad & Tobago; TEDDY sites -> their country). Puerto Rico has
# no estimate (ZIP has no PAH, so no effect here). Time-matched as a
# sample-weighted mean over sampling years; years < 2000 use 2000 (flagged).
#
# Analyses (urinary NAP1, NAP2, PYR1; same data rules as R4_5_PAH_PM25.R):
#   1. Location level: Spearman of location median with % polluting fuel; and
#      PM2.5 vs % polluting fuel across the PAH locations (collinearity).
#   2. Individual level, between-location terms only (location-level covariates):
#        M1: lval ~ log10 PM2.5_loc            + basis + (1|loc) + (1|pid)
#        M2: lval ~ polluting_fuel_loc (per 10%) + basis + (1|loc) + (1|pid)
#        M3: lval ~ both
#      Attenuation of the PM2.5 term from M1 to M3 = how much of the PM2.5
#      association is carried by household solid fuel. ~10 locations, so M3
#      is a descriptive decomposition, not a causal separation.
# Run after R4_5_PAH_PM25.R, from the project root.

suppressMessages({library(dplyr); library(tidyr); library(ggplot2); library(lme4); library(jsonlite); library(patchwork)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
geo <- "Output/Reviewer_Round2/R4_5_determinants/geospatial"; out <- file.path(geo, "pah_pm25")
cache <- file.path(Sys.getenv("HHEAR_GEO_CACHE", path.expand("~/HHEAR_geo_cache")), "cleanfuel")

rd_wb <- function(code) {
  j <- fromJSON(file.path(cache, paste0(code, ".json")))[[2]]
  tibble(iso3 = j$countryiso3code, year = as.integer(j$date), clean = j$value) %>% filter(!is.na(clean))
}
CF <- bind_rows(rd_wb("EG.CFT.ACCS.ZS") %>% mutate(series = "national"),
                rd_wb("EG.CFT.ACCS.UR.ZS") %>% mutate(series = "urban"),
                rd_wb("EG.CFT.ACCS.RU.ZS") %>% mutate(series = "rural"))
write.csv(CF, file.path(geo, "layer_cleanfuel_country_year.csv"), row.names = FALSE)

us_units <- c("Atlanta","Sacramento","Denver","Baltimore","SanFrancisco","Boston","Los_Angeles","Syracuse","New_York",
              "Michigan","California","Washington_State","New_Hampshire","Massachusetts","NorthEast_US","Midwestern_US",
              "United_States","TEDDY_United_States")
iso <- c(Accra = "GHA", Mexico_City = "MEX", Tobago = "TTO", Suriname = "SUR", South_Africa = "ZAF", Ecuador = "ECU",
         Chile = "CHL", Uganda = "UGA", Bangladesh = "BGD", Thailand = "THA", ZIP_Brazil = "BRA", ZIP_Guatemala = "GTM",
         ZIP_Nicaragua = "NIC", `ZIP_Puerto Rico` = "PRI", TEDDY_Finland = "FIN", TEDDY_Germany = "DEU", TEDDY_Sweden = "SWE")
L  <- read.csv(file.path(geo, "location_layers_time_matched.csv"))
yw <- read.csv(file.path(geo, "location_study_year_weights.csv"))
asg <- tibble(loc_id = L$loc_id) %>%
  mutate(iso3 = ifelse(loc_id %in% us_units, "USA", iso[loc_id]),
         series = ifelse(loc_id %in% c("Accra", "Mexico_City"), "urban", "national"))
cf_loc <- yw %>% group_by(loc_id, year) %>% summarise(n = sum(n), .groups = "drop") %>%
  left_join(asg, by = "loc_id") %>% mutate(cf_year = pmin(pmax(year, 2000L), 2023L)) %>%
  left_join(CF, by = c("iso3", "series", "cf_year" = "year")) %>%
  group_by(loc_id, iso3, series) %>%
  summarise(clean_fuel_pct = if (all(is.na(clean))) NA_real_ else weighted.mean(clean, n, na.rm = TRUE),
            frac_years_before_2000 = sum(n[year < 2000]) / sum(n), .groups = "drop") %>%
  mutate(polluting_fuel_pct = 100 - clean_fuel_pct)
L <- L %>% select(-any_of(c("iso3", "cleanfuel_series", "clean_fuel_pct", "polluting_fuel_pct", "cleanfuel_frac_years_before_2000"))) %>%
  left_join(cf_loc %>% rename(cleanfuel_series = series, cleanfuel_frac_years_before_2000 = frac_years_before_2000), by = "loc_id")
write.csv(L, file.path(geo, "location_layers_time_matched.csv"), row.names = FALSE)

# ---- 1. location level ------------------------------------------------------
M <- read.csv(file.path(geo, "biomarker_location_medians.csv")) %>% filter(Analyte_Code %in% c("NAP1", "NAP2", "PYR1"))
Z <- M %>% inner_join(L %>% select(loc_id, PM25_ugm3, NO2_ppb, polluting_fuel_pct), by = "loc_id") %>% filter(!is.na(polluting_fuel_pct))
perm_rho <- function(x, y, B = 20000) { r0 <- cor(x, y, method = "spearman"); set.seed(1)
  rp <- replicate(B, cor(x, sample(y), method = "spearman")); c(rho = r0, p = (1 + sum(abs(rp) >= abs(r0) - 1e-12)) / (B + 1)) }
loc_tab <- Z %>% group_by(analyte = Analyte_Code, set) %>%
  summarise(n = n(), rho_fuel = perm_rho(polluting_fuel_pct, median)[["rho"]], p_fuel = perm_rho(polluting_fuel_pct, median)[["p"]],
            rho_pm = cor(PM25_ugm3, median, method = "spearman"),
            rho_pm_vs_fuel = cor(PM25_ugm3, polluting_fuel_pct, method = "spearman"), .groups = "drop")
write.csv(loc_tab, file.path(out, "D_cleanfuel_location_level.csv"), row.names = FALSE)

# Robustness: polluting-fuel share is ~0-1% at every high-income location, so
# Spearman ranks among those values are noise. Tie everything < 2% at 0, and
# drop Accra (the most influential location); compare with PM2.5 without Accra.
pr <- function(x, y) { r <- perm_rho(x, y); sprintf("%.2f (p = %.2g, n = %d)", r[["rho"]], r[["p"]], length(x)) }
rob <- Z %>% mutate(fuel_tied = ifelse(polluting_fuel_pct < 2, 0, polluting_fuel_pct), noA = loc_id != "Accra") %>%
  group_by(analyte = Analyte_Code, set) %>%
  summarise(fuel = pr(polluting_fuel_pct, median), fuel_tied_lt2pct = pr(fuel_tied, median),
            fuel_drop_Accra = pr(polluting_fuel_pct[noA], median[noA]), fuel_tied_drop_Accra = pr(fuel_tied[noA], median[noA]),
            PM25 = pr(PM25_ugm3, median), PM25_drop_Accra = pr(PM25_ugm3[noA], median[noA]), .groups = "drop")
write.csv(rob, file.path(out, "D_cleanfuel_robustness.csv"), row.names = FALSE)

# ---- 2. individual level: PM2.5 vs solid fuel --------------------------------
d <- readRDS(file.path(out, "pah_sample_level.rds")) %>%
  left_join(L %>% select(loc_id, PM25_ugm3, polluting_fuel_pct), by = "loc_id") %>%
  filter(!is.na(polluting_fuel_pct)) %>% mutate(lpm_loc = log10(PM25_ugm3), fuel10 = polluting_fuel_pct / 10)
fit3 <- function(dd, label) {
  base <- if (length(unique(dd$Dilution_Adjustment)) > 1) " + Dilution_Adjustment + (1|loc_id) + (1|pid)" else " + (1|loc_id) + (1|pid)"
  ms <- list(M1_PM25 = paste("lval ~ lpm_loc", base), M2_fuel = paste("lval ~ fuel10", base), M3_both = paste("lval ~ lpm_loc + fuel10", base))
  bind_rows(lapply(names(ms), function(k) {
    m <- lmer(as.formula(ms[[k]]), data = dd, REML = FALSE, control = lmerControl(calc.derivs = FALSE))
    cf <- summary(m)$coefficients; ci <- suppressMessages(confint(m, method = "Wald"))
    tr <- intersect(c("lpm_loc", "fuel10"), rownames(cf))
    tibble(analyte = unique(dd$Analyte_Code), subset = label, model = k, term = tr, beta = cf[tr, 1], lo = ci[tr, 1], hi = ci[tr, 2],
           AIC = AIC(m), sd_location = as.data.frame(VarCorr(m)) %>% filter(grp == "loc_id") %>% pull(sdcor),
           n_locations = n_distinct(dd$loc_id), n_participants = n_distinct(dd$pid))
  }))
}
D <- bind_rows(lapply(c("NAP1", "NAP2", "PYR1"), function(a) bind_rows(
  fit3(d %>% filter(Analyte_Code == a), "All locations"),
  fit3(d %>% filter(Analyte_Code == a, adjusted), "Dilution-adjusted urine only")))) %>%
  mutate(effect = ifelse(term == "lpm_loc", 10^(beta * log10(2)), 10^beta),
         effect_lo = ifelse(term == "lpm_loc", 10^(lo * log10(2)), 10^lo),
         effect_hi = ifelse(term == "lpm_loc", 10^(hi * log10(2)), 10^hi),
         effect_unit = ifelse(term == "lpm_loc", "fold per PM2.5 doubling", "fold per +10 pct-points polluting fuel"))
att <- D %>% filter(term == "lpm_loc") %>% select(analyte, subset, model, beta) %>% pivot_wider(names_from = model, values_from = beta) %>%
  mutate(pct_PM25_attenuation = 100 * (1 - M3_both / M1_PM25))
write.csv(D, file.path(out, "D_cleanfuel_individual_models.csv"), row.names = FALSE)
write.csv(att, file.path(out, "D_PM25_attenuation_by_fuel.csv"), row.names = FALSE)

print(cf_loc %>% filter(loc_id %in% unique(Z$loc_id)) %>% arrange(desc(polluting_fuel_pct)) %>% as.data.frame(), digits = 3)
print(loc_tab %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())
print(D %>% select(analyte, subset, model, term, effect, effect_lo, effect_hi, AIC, sd_location) %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())
print(att %>% mutate(across(where(is.numeric), ~ round(.x, 3))) %>% as.data.frame())

# ---- Figure -----------------------------------------------------------------
grp <- function(l) case_when(l %in% us_units ~ "US", grepl("^TEDDY_(Finland|Sweden|Germany)", l) ~ "Europe", TRUE ~ "Latin America & Africa")
cols <- c(US = "#3b6fb6", Europe = "#2e8b57", `Latin America & Africa` = "#c0392b")
LL <- L %>% filter(!is.na(polluting_fuel_pct)) %>% mutate(group = grp(loc_id), pah = loc_id %in% unique(Z$loc_id))
pE <- ggplot(LL, aes(polluting_fuel_pct + 0.5, PM25_ugm3, colour = group)) +
  geom_point(aes(shape = pah), size = 2.4) +
  ggrepel::geom_text_repel(aes(label = gsub("_", " ", loc_id)), size = 2.1, colour = "grey30", seed = 2, max.overlaps = 40) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = c(`TRUE` = 16, `FALSE` = 1), labels = c(`TRUE` = "Has urinary PAH", `FALSE` = "No PAH data"), name = NULL) +
  labs(title = "D. Household polluting-fuel use vs ambient PM2.5 (all 35 locations)",
       x = "% population relying on polluting cooking fuels (+0.5, log scale)", y = "PM2.5 (µg/m³)") +
  theme_bw(base_size = 9) + theme(legend.position = "bottom")
labF <- loc_tab %>% filter(set == "all") %>% transmute(Analyte_Code = analyte, lab = sprintf("ρ = %.2f (perm p = %.2g)", rho_fuel, p_fuel))
pF <- ggplot(Z %>% filter(set == "all") %>% mutate(group = grp(loc_id)), aes(polluting_fuel_pct + 0.5, median)) +
  geom_smooth(method = "lm", formula = y ~ x, colour = "grey30", fill = "grey85", linewidth = 0.5) +
  geom_point(aes(colour = group, shape = df < 0.5), size = 2.3) +
  ggrepel::geom_text_repel(aes(label = gsub("_", " ", loc_id)), size = 2.1, colour = "grey30", seed = 1, max.overlaps = 30) +
  geom_text(data = labF, aes(label = lab), x = -Inf, y = Inf, hjust = -0.05, vjust = 1.4, size = 2.5, inherit.aes = FALSE) +
  scale_x_log10() + scale_y_log10() + scale_colour_manual(values = cols, name = NULL) +
  scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 1), labels = c("≥50% detected", "<50% detected"), name = NULL) +
  facet_wrap(~Analyte_Code, scales = "free_y") +
  labs(title = "E. Urinary PAH location medians vs polluting-fuel use", x = "% polluting cooking fuels (+0.5, log scale)",
       y = "Location median (ng/mL or µg/g creatinine)") + theme_bw(base_size = 9) + theme(legend.position = "bottom")
fd <- D %>% mutate(model = factor(model, levels = rev(c("M1_PM25", "M2_fuel", "M3_both")),
                                  labels = rev(c("PM2.5 alone", "Polluting fuel alone", "Both"))),
                   term = ifelse(term == "lpm_loc", "PM2.5 (per doubling)", "Polluting fuel (per +10 pct-pts)"))
pG <- ggplot(fd, aes(effect, model, colour = term)) + geom_vline(xintercept = 1, linetype = 2) +
  geom_pointrange(aes(xmin = effect_lo, xmax = effect_hi), position = position_dodge(width = 0.5), size = 0.25) +
  scale_x_log10() + facet_grid(subset ~ analyte) + scale_colour_manual(values = c("#7d3c98", "#555555"), name = NULL) +
  labs(title = "F. Individual-level models: does solid-fuel use account for the PM2.5 association?",
       x = "Fold change in urinary PAH (log scale)", y = NULL) + theme_bw(base_size = 9) + theme(legend.position = "bottom")
g <- (pE | pF) / pG + plot_layout(widths = c(1, 2), heights = c(1, 1.1))
ggsave(file.path(out, "SuppFig_PAH_CleanFuel.png"), g, width = 13, height = 10, dpi = 200)
ggsave(file.path(out, "SuppFig_PAH_CleanFuel.pdf"), g, width = 13, height = 10, device = cairo_pdf)
