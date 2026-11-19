# R4_5_geospatial_no2.R
# Adds surface NO2 (Anenberg et al. 2022, Lancet Planet Health; Mohegh & Anenberg
# figshare 12968114; 1 km annual, ppb) to the location layer table built by
# R4_5_geospatial_extract.R. Each raw year is 2.8 GB, so years are processed one
# at a time from ~/HHEAR_geo_cache/no2/ and the per-location values are appended
# to layer_NO2_location_year.csv; the raw .nc can then be deleted.
# Time matching: each location's sampling-year weights are mapped to the NEAREST
# NO2 year available in that CSV (flagged by the year gap).
# Run from project root, after R4_5_geospatial_extract.R.

suppressMessages({library(dplyr); library(sf); library(terra); library(exactextractr)})
sf_use_s2(FALSE)
root  <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
cache <- Sys.getenv("HHEAR_GEO_CACHE", path.expand("~/HHEAR_geo_cache"))
out   <- file.path(root, "Output/Reviewer_Round2/R4_5_determinants/geospatial")
U     <- st_read(file.path(out, "location_units.gpkg"), quiet = TRUE)
csv   <- file.path(out, "layer_NO2_location_year.csv")

pop_epochs <- c(1990, 2000, 2010, 2020)
pops <- setNames(lapply(pop_epochs, function(ep) rast(list.files(file.path(cache, "ghsl"), pattern = paste0("E", ep, ".*\\.tif$"), full.names = TRUE)[1])), pop_epochs)
pw_mean <- function(r, ep) sapply(seq_len(nrow(U)), function(i) {
  g <- U[i, ]; rc <- crop(r, ext(vect(g)) + 0.2, snap = "out")
  w <- resample(crop(pops[[ep]], ext(rc) + 0.02, snap = "out"), rc, method = "bilinear")
  exact_extract(rc, g, "weighted_mean", weights = w, progress = FALSE)
})

prev <- if (file.exists(csv)) read.csv(csv) else data.frame()
for (f in list.files(file.path(cache, "no2"), pattern = "_final_1km\\.nc$", full.names = TRUE)) {
  y <- as.integer(sub(".*NO2_(\\d{4})_.*", "\\1", basename(f)))
  if (nrow(prev) && y %in% prev$year) next
  r <- rast(f); if (nlyr(r) > 1) r <- r[[1]]
  ep <- as.character(pop_epochs[which.min(abs(pop_epochs - y))])
  prev <- bind_rows(prev, tibble(loc_id = U$loc_id, year = y, layer = "NO2_ppb", value = pw_mean(r, ep)))
  cat("NO2", y, "done\n")
}
write.csv(prev, csv, row.names = FALSE)

# Time-match onto the location table
yw <- read.csv(file.path(out, "location_study_year_weights.csv"))
no2_years <- sort(unique(prev$year))
m <- yw %>% group_by(loc_id, year) %>% summarise(n = sum(n), .groups = "drop") %>%
  mutate(no2_year = sapply(year, function(y) no2_years[which.min(abs(no2_years - y))]), gap = abs(year - no2_year)) %>%
  left_join(prev %>% select(loc_id, no2_year = year, NO2 = value), by = c("loc_id", "no2_year")) %>%
  group_by(loc_id) %>% summarise(NO2_ppb = weighted.mean(NO2, n), NO2_mean_year_gap = weighted.mean(gap, n), .groups = "drop")
L <- read.csv(file.path(out, "location_layers_time_matched.csv")) %>% select(-any_of(c("NO2_ppb", "NO2_mean_year_gap"))) %>% left_join(m, by = "loc_id")
write.csv(L, file.path(out, "location_layers_time_matched.csv"), row.names = FALSE)
cat("NO2 years available:", no2_years, "\n"); print(L %>% select(loc_id, PM25_ugm3, NO2_ppb, NO2_mean_year_gap) %>% as.data.frame(), digits = 3)
