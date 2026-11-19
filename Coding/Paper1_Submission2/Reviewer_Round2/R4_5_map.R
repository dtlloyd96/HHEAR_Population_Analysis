# R4_5_map.R
# Reviewer 4 -- publication map: HHEAR locations over the three environmental
# layers with the clearest determinant signals, each location shaded by its
# biomarker median.
#   A  Ambient PM2.5 (ACAG V6.GL.02.04, 2015, aggregated to 0.5 deg) + urinary 1-hydroxypyrene
#   B  % relying on polluting cooking fuels (WHO via World Bank, 2015, national) + urinary 1-hydroxypyrene
#   C  P(groundwater As > 10 ug/L) (Podgorski & Berg 2020, aggregated to 0.5 deg) + urinary As
#   D  National insecticide use (FAOSTAT, kg a.i./ha cropland, 2015) + urinary 3-PBA
# Location points sit at each unit's interior point; locations without that
# biomarker are small grey dots. Layer = fill scale; biomarker = point colour scale. Run after the other R4_5 scripts, from project root.

suppressMessages({library(dplyr); library(ggplot2); library(sf); library(terra); library(rnaturalearth); library(jsonlite); library(readr); library(patchwork); library(tidyr)})
sf_use_s2(FALSE)
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
cache <- Sys.getenv("HHEAR_GEO_CACHE", path.expand("~/HHEAR_geo_cache"))
geo <- "Output/Reviewer_Round2/R4_5_determinants/geospatial"
out <- "Output/Reviewer_Round2/R4_5_determinants/map"; dir.create(out, showWarnings = FALSE)
XL <- c(-128, 100); YL <- c(-38, 68)

world <- ne_countries(scale = 50, returnclass = "sf")
U <- st_read(file.path(geo, "location_units.gpkg"), quiet = TRUE)
pts <- suppressWarnings(st_point_on_surface(U)) %>% mutate(lon = st_coordinates(.)[, 1], lat = st_coordinates(.)[, 2]) %>% st_drop_geometry()
M1 <- read.csv(file.path(geo, "biomarker_location_medians.csv")) %>% filter(set == "all")
M2 <- read.csv("Output/Reviewer_Round2/R4_5_determinants/diet_pesticide_lead/location_medians_with_layers.csv") %>% filter(set == "all")
bio <- function(a, src) src %>% filter(Analyte_Code == a) %>% select(loc_id, median)

rast_df <- function(r, fact) { r <- aggregate(r, fact = fact, fun = "mean", na.rm = TRUE); r <- crop(r, ext(XL[1], XL[2], YL[1], YL[2]))
  df <- as.data.frame(r, xy = TRUE, na.rm = TRUE); names(df)[3] <- "v"; df }
pm <- rast(file.path(cache, "pm25/V6GL0204_PM25_0p10_2015.nc")); pm[pm < 0] <- NA
pm_df <- rast_df(pm, 5)
asr <- rast(list.files(file.path(cache, "as_gw"), pattern = "probability_arsenic_gt_10ppb\\.tif$", recursive = TRUE, full.names = TRUE)[1])
as_df <- rast_df(asr, 60)

cf <- fromJSON(file.path(cache, "cleanfuel/EG.CFT.ACCS.ZS_all_2015.json"))[[2]]
cf <- tibble(iso_a3 = cf$countryiso3code, polluting = 100 - cf$value) %>% filter(!is.na(polluting), nchar(iso_a3) == 3)
P <- read_csv(file.path(cache, "faostat/Inputs_Pesticides_Use_E_All_Data/Inputs_Pesticides_Use_E_All_Data_NOFLAG.csv"), show_col_types = FALSE)
LU <- read_csv(file.path(cache, "faostat/Inputs_LandUse_E_All_Data/Inputs_LandUse_E_All_Data_NOFLAG.csv"), show_col_types = FALSE)
ins <- P %>% filter(`Item Code` == 1309, Element == "Agricultural Use") %>% select(m49 = `Area Code (M49)`, t = Y2015) %>%
  inner_join(LU %>% filter(Item == "Cropland", Element == "Area") %>% select(m49 = `Area Code (M49)`, kha = Y2015), by = "m49") %>%
  transmute(m49 = as.integer(sub("^'", "", m49)), insecticide = t / kha) %>% filter(is.finite(insecticide))
wcf  <- world %>% left_join(cf, by = c("iso_a3_eh" = "iso_a3"))
wins <- world %>% mutate(m49 = suppressWarnings(as.integer(un_a3))) %>% left_join(ins, by = "m49")

base_theme <- theme_void(base_size = 11) + theme(legend.position = "bottom", legend.key.width = unit(1.4, "cm"), legend.key.height = unit(0.3, "cm"),
                                                  plot.title = element_text(face = "bold", size = 12, colour = "#2B2F33", margin = margin(b = 4)),
                                                  legend.title = element_text(size = 9, colour = "#5A6068"), legend.text = element_text(size = 8.5),
                                                  plot.margin = margin(6, 8, 6, 8))
add_points <- function(g, b, lab) {
  p <- pts %>% left_join(b, by = "loc_id")
  g + geom_point(data = p %>% filter(is.na(median)), aes(lon, lat), colour = "grey45", size = 0.7) +
    geom_point(data = p %>% filter(!is.na(median)), aes(lon, lat), shape = 16, size = 4.6, colour = "white") +
    geom_point(data = p %>% filter(!is.na(median)), aes(lon, lat, colour = median), shape = 16, size = 3.6) +
    ggrepel::geom_text_repel(data = p %>% filter(!is.na(median)), aes(lon, lat, label = gsub("_", " ", sub("^ZIP_", "", sub("^TEDDY_United_States$", "TEDDY US sites", sub("^TEDDY_(Finland|Sweden|Germany)$", "\\1", loc_id))))),
                             size = 2.8, colour = "grey10", seed = 3, max.overlaps = 40, min.segment.length = 0.1, box.padding = 0.25,
                             bg.color = "white", bg.r = 0.12) +
    scale_colour_viridis_c(option = "magma", direction = -1, trans = "log10", name = lab,
                           guide = guide_colourbar(title.position = "top", order = 2))
}
mk <- function(layer, fillscale, title) ggplot() + layer + fillscale +
  coord_sf(xlim = XL, ylim = YL, expand = FALSE) + labs(title = title) + base_theme


gA <- add_points(mk(list(geom_raster(data = pm_df, aes(x, y, fill = v)), geom_sf(data = world, fill = NA, colour = "grey40", linewidth = 0.1)),
                    scale_fill_gradientn(colours = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5", "#08306b"), trans = "sqrt", limits = c(0, 80), oob = scales::squish,
                                         name = "PM2.5, 2015 (µg/m³)", guide = guide_colourbar(title.position = "top", order = 1)),
                    "A  Ambient PM2.5 and urinary 1-hydroxypyrene"), bio("PYR1", M1), "1-hydroxypyrene median")
gB <- add_points(mk(geom_sf(data = wcf, aes(fill = polluting), colour = "grey40", linewidth = 0.1),
                    scale_fill_gradientn(colours = c("#fff7ec", "#fdd49e", "#fc8d59", "#d7301f", "#7f0000"), limits = c(0, 100), na.value = "grey90",
                                         name = "% relying on polluting cooking fuels, 2015", guide = guide_colourbar(title.position = "top", order = 1)),
                    "B  Household solid-fuel use and urinary 1-hydroxypyrene"), bio("PYR1", M1), "1-hydroxypyrene median")
gC <- add_points(mk(list(geom_raster(data = as_df, aes(x, y, fill = v)), geom_sf(data = world, fill = NA, colour = "grey40", linewidth = 0.1)),
                    scale_fill_gradientn(colours = c("#f7fcf5", "#c7e9c0", "#74c476", "#238b45", "#00441b"), limits = c(0, 1),
                                         name = "P(groundwater As > 10 µg/L)", guide = guide_colourbar(title.position = "top", order = 1)),
                    "C  Groundwater arsenic and urinary arsenic"), bio("As", M2), "Urinary As median")
gD <- add_points(mk(geom_sf(data = wins, aes(fill = insecticide), colour = "grey40", linewidth = 0.1),
                    scale_fill_gradientn(colours = c("#fcfbfd", "#dadaeb", "#9e9ac8", "#6a51a3", "#3f007d"), trans = "log10", limits = c(0.01, 10), oob = scales::squish,
                                         na.value = "grey90", name = "Insecticide use, 2015 (kg a.i./ha cropland)", guide = guide_colourbar(title.position = "top", order = 1)),
                    "D  Agricultural insecticide use and urinary 3-PBA (pyrethroids)"), bio("PBA", M2), "3-PBA median")
g <- (gA | gB) / (gC | gD) + plot_annotation(caption = "Points: HHEAR locations (unit interior point); filled = location median of the paired biomarker (log scale), grey = biomarker not measured there.\nUS multi-site cohorts are shown at the US interior point. Layers: ACAG V6.GL.02.04; WHO/World Bank SDG 7.1.2; Podgorski & Berg 2020; FAOSTAT.") &
  theme(plot.caption = element_text(size = 9, hjust = 0, colour = "#5A6068"))
ggsave(file.path(out, "Fig_Geospatial_Determinants_Map.png"), g, width = 15, height = 9.5, dpi = 220, bg = "white")
ggsave(file.path(out, "Fig_Geospatial_Determinants_Map.pdf"), g, width = 15, height = 9.5, device = cairo_pdf, bg = "white")
# Panels (without letters) for Figure S2, which stacks the maps above the scatter panels.
strip_letter <- function(x) { x$labels$title <- sub("^[A-D]  ", "", x$labels$title); x }
saveRDS(lapply(list(A = gA, B = gB, C = gC, D = gD), strip_letter), file.path(out, "map_panels.rds"))
cat("map written\n")
