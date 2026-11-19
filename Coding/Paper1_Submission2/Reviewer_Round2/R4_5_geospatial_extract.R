# R4_5_geospatial_extract.R
# Reviewer 4 (determinants of geographic differences) -- geospatial exposure layer.
# Builds one spatial unit per HHEAR location, then extracts population-weighted
# environmental layers matched to each location's sampling years.
#
# Raw rasters live OUTSIDE Google Drive in ~/HHEAR_geo_cache (not synced):
#   pm25/  ACAG V6.GL.02.04 annual PM2.5, 0.1 deg, 1998-2023 (CC BY 4.0)
#   ghsl/  GHS-POP R2023A 30 arcsec, epochs 1990/2000/2010/2020 (weights)
#   as_gw/ Podgorski & Berg 2020 probability groundwater As > 10 ug/L (static)
# Outputs: Output/Reviewer_Round2/R4_5_determinants/geospatial/
#
# Spatial unit rule (PLAN_R4_Geographic_Determinants.md 5.2):
#   city label -> 30 km buffer; state/region -> state polygon(s);
#   country / multi-site "United_States" -> country polygon;
#   TEDDY -> 50 km buffers around published clinical centres, unioned per country.

suppressMessages({library(dplyr); library(sf); library(terra); library(exactextractr); library(rnaturalearth)})
sf_use_s2(FALSE)

root  <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
cache <- Sys.getenv("HHEAR_GEO_CACHE", path.expand("~/HHEAR_geo_cache"))
out   <- file.path(root, "Output/Reviewer_Round2/R4_5_determinants/geospatial")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

# ---- 1. Location -> unit definitions ---------------------------------------
# One row per location label used in the analyses (Harmonized_SDOH `Location`,
# TEDDY split by country, ZIP split by `Country`).
pt <- function(lon, lat) list(c(lon, lat))
units <- tibble::tribble(
  ~loc_id,               ~type,     ~spec,
  "Atlanta",             "city",    pt(-84.388, 33.749),
  "Sacramento",          "city",    pt(-121.494, 38.582),
  "Denver",              "city",    pt(-104.990, 39.739),
  "Baltimore",           "city",    pt(-76.612, 39.290),
  "Mexico_City",         "city",    pt(-99.133, 19.433),
  "SanFrancisco",        "city",    pt(-122.419, 37.775),
  "Accra",               "city",    pt(-0.187, 5.604),
  "Boston",              "city",    pt(-71.059, 42.360),
  "Los_Angeles",         "city",    pt(-118.244, 34.052),
  "Syracuse",            "city",    pt(-76.147, 43.048),
  "New_York",            "city",    pt(-74.006, 40.713),
  "Tobago",              "city",    pt(-60.700, 11.230),   # island; buffer covers it
  "Michigan",            "state",   list("michigan"),
  "California",          "state",   list("california"),
  "Washington_State",    "state",   list("washington"),
  "New_Hampshire",       "state",   list("new hampshire"),
  "Massachusetts",       "state",   list("massachusetts"),
  "NorthEast_US",        "state",   list(c("connecticut","massachusetts","rhode island","new hampshire","vermont","maine")),
  "Midwestern_US",       "state",   list(c("illinois","indiana","iowa","kansas","michigan","minnesota","missouri","nebraska","north dakota","ohio","south dakota","wisconsin")),
  "United_States",       "country", list("United States of America"),
  "Suriname",            "country", list("Suriname"),
  "South_Africa",        "country", list("South Africa"),
  "Ecuador",             "country", list("Ecuador"),
  "Chile",               "country", list("Chile"),
  "Uganda",              "country", list("Uganda"),
  "Bangladesh",          "country", list("Bangladesh"),
  "Thailand",            "country", list("Thailand"),
  "ZIP_Brazil",          "country", list("Brazil"),
  "ZIP_Guatemala",       "country", list("Guatemala"),
  "ZIP_Nicaragua",       "country", list("Nicaragua"),
  "ZIP_Puerto Rico",     "country", list("Puerto Rico"),
  # TEDDY clinical centres (published): Denver CO, Augusta GA, Gainesville FL,
  # Seattle WA; Turku, Oulu, Tampere FI; Munich DE; Malmo/Skane SE.
  "TEDDY_United_States", "teddy",   list(rbind(c(-104.990,39.739), c(-81.975,33.474), c(-82.325,29.652), c(-122.332,47.606))),
  "TEDDY_Finland",       "teddy",   list(rbind(c(22.267,60.452), c(25.473,65.012), c(23.761,61.498))),
  "TEDDY_Germany",       "teddy",   list(rbind(c(11.582,48.135))),
  "TEDDY_Sweden",        "teddy",   list(rbind(c(13.004,55.605)))
)

countries <- ne_countries(scale = 50, returnclass = "sf") %>% select(name = admin)
states    <- st_as_sf(maps::map("state", fill = TRUE, plot = FALSE)) %>% rename(name = ID)
st_crs(states) <- 4326

buf <- function(m, km) st_union(st_buffer(st_transform(st_sfc(lapply(seq_len(nrow(m)), function(i) st_point(m[i, ])), crs = 4326), 3857), km * 1000)) %>% st_transform(4326)
make_geom <- function(type, spec) {
  s <- spec[[1]]
  switch(type,
    city    = buf(matrix(s, nrow = 1), 30),
    teddy   = buf(s, 50),
    state   = st_union(states[states$name %in% s, ]),
    country = { g <- countries[countries$name %in% s, ]; stopifnot(nrow(g) == 1); st_union(g) })
}
geoms <- do.call(c, Map(make_geom, units$type, units$spec))
U <- st_sf(units %>% select(loc_id, type), geometry = st_make_valid(geoms))
stopifnot(all(!st_is_empty(U)))
st_write(U, file.path(out, "location_units.gpkg"), delete_dsn = TRUE, quiet = TRUE)

# ---- 2. Sampling-year weights per location ---------------------------------
e <- new.env(); load(file.path(root, "Input/Harmonized_Datasets/Harmonized_Targeted.RData"), envir = e)
s <- new.env(); load(file.path(root, "Input/Harmonized_Datasets/Harmonized_SDOH.RData"),   envir = s)
TD <- e$Targeted_Data_All; SD <- s$SDOH_Data_All

# Study -> loc_id, using the same labels as the location analyses. MADRES is
# labelled California in SDOH but is a Los Angeles cohort; it stays California
# here so units match the B4 location factor exactly.
loc_map <- SD %>% distinct(Study, Child_PID, Mat_PID, Location, Country) %>%
  mutate(Study = sub("_SDOH$", "_Targeted", Study),
         loc_id = case_when(Study == "TEDDY_Targeted" ~ paste0("TEDDY_", Location),
                            Study == "ZIP_Targeted" & !is.na(Country) ~ paste0("ZIP_", Country),
                            TRUE ~ Location))

# Studies with no Sample_Collection_Year: windows verified against cohort papers
# (2026-11-05; register/verified_dates.md). Each entry is year -> weight.
#   ATLAA   enrollment serum Mar 2014-May 2018 (Chang 2021); HHEAR data are blood/DBS only
#   ELEMENT three periods (HHEAR listing; Perng cohort profile), weighted by the
#           raw file's sample counts: pregnancy 1997-2004 (n=132), early
#           adolescence 2008-2012 (240), late adolescence 2015-2018 (515)
#   MEBCO   enrollment Jan 2008-Jun 2011 (Ahmed 2020)
#   SWSD    SAWASDEE, recruitment Jul 2017-Jun 2019 (Baumert 2022)
#   ZIP     enrollment 2016-2018; fallback only (mothers' visit years used, below)
even <- function(y) setNames(rep(1 / length(y), length(y)), y)
assumed_years <- list(ATLAA_Targeted = even(2014:2018),
                      ELEMENT_Targeted = c(even(1997:2004) * 132, even(2008:2012) * 240, even(2015:2018) * 515) / 887,
                      MEBCO_Targeted = even(2008:2011), SWSD_Targeted = even(2017:2019), ZIP_Targeted = even(2016:2018))

samp <- TD %>% filter(!is.na(Study)) %>% distinct(Study, Child_PID, Mat_PID, SID, Sample_Collection_Year)

# ZIP's sample_collection_year is blank, but its raw EPI_3 table has each
# mother's visit years (visit_date_yr_1..7; 2016-2018). Use the mother's median
# visit year (2026-11-05), so ZIP is no longer assumed.
z <- new.env(); load(file.path(root, "Input/ZIP_1839_Data_All.RData"), envir = z)
zip_yr <- z[[ls(z)[1]]]$ZIP_EPI_3 %>% select(Mat_PID = MOTHER_PID, starts_with("visit_date_yr_")) %>%
  tidyr::pivot_longer(-Mat_PID, values_to = "yr") %>% filter(yr > 2000, yr < 2030) %>%
  group_by(Mat_PID = as.character(Mat_PID)) %>% summarise(zip_year = round(median(yr)), .groups = "drop")
samp <- samp %>% left_join(zip_yr, by = "Mat_PID") %>%
  mutate(Sample_Collection_Year = ifelse(Study == "ZIP_Targeted" & is.na(Sample_Collection_Year), zip_year, Sample_Collection_Year)) %>%
  select(-zip_year)
samp_loc <- bind_rows(
  samp %>% filter(!Study %in% c("TEDDY_Targeted", "ZIP_Targeted")) %>%
    left_join(loc_map %>% distinct(Study, loc_id), by = "Study"),
  samp %>% filter(Study == "TEDDY_Targeted") %>%
    left_join(loc_map %>% filter(Study == "TEDDY_Targeted") %>% distinct(Child_PID, loc_id), by = "Child_PID"),
  samp %>% filter(Study == "ZIP_Targeted") %>%
    left_join(loc_map %>% filter(Study == "ZIP_Targeted", !is.na(Mat_PID)) %>% distinct(Mat_PID, loc_id), by = "Mat_PID")
) %>% filter(!is.na(loc_id), loc_id %in% U$loc_id)

yw_obs <- samp_loc %>% filter(!is.na(Sample_Collection_Year)) %>%
  count(loc_id, Study, year = Sample_Collection_Year, name = "n") %>% mutate(year_source = "observed")
yw_ass <- samp_loc %>% filter(is.na(Sample_Collection_Year), Study %in% names(assumed_years)) %>%
  count(loc_id, Study, name = "n_total") %>%
  rowwise() %>% reframe(loc_id, Study, year = as.numeric(names(assumed_years[[Study]])), n = n_total * unname(assumed_years[[Study]])) %>%
  mutate(year_source = "assumed")
yw <- bind_rows(yw_obs, yw_ass)
write.csv(yw, file.path(out, "location_study_year_weights.csv"), row.names = FALSE)
cat("Locations without any sampling years:", setdiff(U$loc_id, yw$loc_id), "\n")

# ---- 3. Extraction ----------------------------------------------------------
pop_epochs <- c(1990, 2000, 2010, 2020)
pop_file <- function(ep) {
  d <- file.path(cache, "ghsl"); tif <- list.files(d, pattern = paste0("E", ep, ".*\\.tif$"), full.names = TRUE, recursive = TRUE)
  if (!length(tif)) { unzip(file.path(d, sprintf("GHS_POP_E%d_GLOBE_R2023A_4326_30ss_V1_0.zip", ep)), exdir = d)
    tif <- list.files(d, pattern = paste0("E", ep, ".*\\.tif$"), full.names = TRUE, recursive = TRUE) }
  rast(tif[1])
}
pops <- setNames(lapply(pop_epochs, pop_file), pop_epochs)
nearest_epoch <- function(y) as.character(pop_epochs[which.min(abs(pop_epochs - y))])

# GHS-POP's 30" grid is offset half a cell from the layer grids, so per unit the
# population is crop-resampled onto the layer's own grid ("sum" when the layer is
# coarser, e.g. PM2.5 0.1 deg) and used as weights there.
pw_mean <- function(r, ep) sapply(seq_len(nrow(U)), function(i) {
  g  <- U[i, ]; bb <- ext(vect(g)) + 0.2
  rc <- crop(r, bb, snap = "out")
  pc <- crop(pops[[ep]], ext(rc) + 0.02, snap = "out")
  w  <- resample(pc, rc, method = if (res(rc)[1] > 1.5 * res(pc)[1]) "sum" else "bilinear")
  exact_extract(rc, g, "weighted_mean", weights = w, progress = FALSE)
})
# Population count and density per unit (2020 epoch, plus epoch-matched below)
area_km2 <- as.numeric(st_area(st_transform(U, 6933))) / 1e6

pm_years <- 1998:2023
# All PM2.5 years share one 0.1-deg grid, so the per-unit population weights on
# that grid are computed once per GHS-POP epoch and reused for every year.
pm_template <- rast(file.path(cache, "pm25", "V6GL0204_PM25_0p10_2015.nc"))
pm_w <- lapply(setNames(as.character(pop_epochs), pop_epochs), function(ep) lapply(seq_len(nrow(U)), function(i) {
  rc <- crop(pm_template, ext(vect(U[i, ])) + 0.2, snap = "out")
  resample(crop(pops[[ep]], ext(rc) + 0.02, snap = "out"), rc, method = "sum")
}))
pm_long <- bind_rows(lapply(pm_years, function(y) {
  r <- rast(file.path(cache, "pm25", sprintf("V6GL0204_PM25_0p10_%d.nc", y))); ep <- nearest_epoch(y)
  v <- sapply(seq_len(nrow(U)), function(i) {
    w <- pm_w[[ep]][[i]]; rc <- crop(r, w); rc[rc < 0] <- NA
    exact_extract(rc, U[i, ], "weighted_mean", weights = w, progress = FALSE)
  })
  tibble(loc_id = U$loc_id, year = y, layer = "PM25_ugm3", value = v)
}))

pop_long <- bind_rows(lapply(pop_epochs, function(ep) {
  tot <- exact_extract(pops[[as.character(ep)]], U, "sum", progress = FALSE)
  tibble(loc_id = U$loc_id, epoch = ep, pop_total = tot, pop_density_km2 = tot / area_km2)
}))

as_tif <- list.files(file.path(cache, "as_gw"), pattern = "probability_arsenic_gt_10ppb\\.tif$", full.names = TRUE, recursive = TRUE)
as_static <- if (length(as_tif)) {
  r <- rast(as_tif[1])
  tibble(loc_id = U$loc_id, layer = "GW_As_prob_gt10", value = pw_mean(r, "2020"))
} else tibble()

write.csv(pm_long,  file.path(out, "layer_PM25_location_year.csv"), row.names = FALSE)
write.csv(pop_long, file.path(out, "layer_population_location_epoch.csv"), row.names = FALSE)
if (nrow(as_static)) write.csv(as_static, file.path(out, "layer_GW_As_location_static.csv"), row.names = FALSE)

# ---- 4. Time-matched location x layer table ---------------------------------
# PM2.5: weighted mean over the location's sampling years (weights = samples/yr).
# Years outside 1998-2023 take the nearest available year and are flagged.
yl <- yw %>% group_by(loc_id, year) %>% summarise(n = sum(n), assumed = any(year_source == "assumed"), .groups = "drop") %>%
  mutate(pm_year = pmin(pmax(year, min(pm_years)), max(pm_years)), pm_out_of_range = year != pm_year,
         epoch = pop_epochs[sapply(year, function(y) which.min(abs(pop_epochs - y)))])
matched <- yl %>%
  left_join(pm_long %>% select(loc_id, pm_year = year, PM25 = value), by = c("loc_id", "pm_year")) %>%
  left_join(pop_long, by = c("loc_id", "epoch")) %>%
  group_by(loc_id) %>%
  summarise(year_min = min(year), year_max = max(year), n_samples = sum(n),
            PM25_ugm3 = weighted.mean(PM25, n), pop_density_km2 = weighted.mean(pop_density_km2, n),
            frac_years_pm_out_of_range = sum(n[pm_out_of_range]) / sum(n),
            years_assumed = any(assumed), .groups = "drop") %>%
  left_join(U %>% st_drop_geometry(), by = "loc_id")
if (nrow(as_static)) matched <- matched %>% left_join(as_static %>% select(loc_id, GW_As_prob_gt10 = value), by = "loc_id")
write.csv(matched, file.path(out, "location_layers_time_matched.csv"), row.names = FALSE)
print(as.data.frame(matched), digits = 3)
