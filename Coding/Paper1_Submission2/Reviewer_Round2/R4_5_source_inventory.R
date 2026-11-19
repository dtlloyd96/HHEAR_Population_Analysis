# R4_5_source_inventory.R
# Reviewer 4 -- systematic inventory of exposure-source questions in every
# HHEAR cohort's raw study file (all non-chemical tables), matched to the
# chemical classes the cohort measured. Identifies which source -> chemical
# pairs can be tested within cohorts, and which recur across cohorts.
#
# A variable is kept if its name (or harmonized ontology name) matches a source
# domain pattern, it has >= 30 non-missing values and >= 2 distinct values, and
# it is not a biomarker / LOD / outcome column. Domains and the chemical classes
# each one plausibly drives:
#   Tobacco           smok, cig, tobacco, nicotine, vape, ecig          -> Smoking, PAH, TraceElements (Cd)
#   Combustion/fuel   stove, fuel, cook, kerosene, wood, gas, heating,  -> PAH
#                     fire, incense, candle
#   Traffic/outdoor   traffic, road, highway, distroad, NO2, PM         -> PAH
#   Diet              fish, seafood, rice, food, canned, packag, meat,  -> TraceElements (As, Hg), Phthalates,
#                     dairy, fast food, milk, eat, diet, plastic, bottle   Phenols, PFAS
#   Water             water, well, tap                                  -> TraceElements, PFAS
#   Pesticides/pests  pest, insect, repellent, flea, tick, mosquito,     -> Pesticides
#                     spray, garden, lawn
#   Personal care     lotion, hair, cosmetic, nail, perfume, fragrance,  -> Phthalates, Parabens, Phenols
#                     soap, shampoo, sunscreen, makeup
#   Housing           paint, renovat, floor, carpet, home type, house,   -> TraceElements (Pb), Phthalates,
#                     dwelling, built, mold, furniture                     FlameRetardants
#   Occupation        occupation, job, work, employ, farm, agricult      -> any
# Writes Output/Reviewer_Round2/R4_5_determinants/source_inventory/. Run from project root.

suppressMessages({library(dplyr); library(tidyr); library(ggplot2)})
root <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
setwd(root)
out <- "Output/Reviewer_Round2/R4_5_determinants/source_inventory"; dir.create(out, recursive = TRUE, showWarnings = FALSE)

dom <- list(
  "Tobacco"          = "smok|cig|tobacco|nicotin|vape|ecig|e_cig",
  "Combustion/fuel"  = "stove|fuel|cook|kerosene|paraffin|woodstove|firewood|wood_?burn|charcoal|heatsource|heat_source|heating|fireplace|incense|candle|fossil|biomass",
  "Traffic/outdoor"  = "traffic|road|highway|freeway|no2|pm10|pm25|pm2_5|benzene|toluene|ambient",
  "Diet"             = "fish|seafood|shellfish|tuna|rice|canned|packag|foodpack|fastfood|fast_food|meat|dairy|(^|_)diet|(^|_)eat|food|beverage|bottle|plastic|microwave|juice|soda",
  "Water"            = "water|(^|_)wells?($|_)|(^|_)tap($|_)",
  "Pesticides/pests" = "pest|insect|repellent|flea|tick|mosq|roach|garden|lawn|herbicid|fumig",
  "Personal care"    = "lotion|(^|_)hair|hairspray|cosmetic|nail|perfume|fragran|soap|shampoo|sunscreen|makeup|deodor|personal.?care",
  "Housing"          = "paint|renovat|floor|carpet|hometype|home_type|house_type|housetype|dwelling|built|mold|mould|furniture|hardwood|residential|vinyl",
  "Occupation"       = "occupat|job|work|employ|farm|agricult")
excl <- "household-income|household_income|income|number_in|number_children|members-in|adults-in|steatosis|tanner|pubic|bacteria\\.|threat|treat|out-of-seat|heatprotector|assistance|allergy|concentration|limit.of|_lod|lod_|weight|wellbeing|well.being|wellness|breath|asthma|wheez|bmi|height|age_|score|depress|stress|_comment|creat|specific.gravity|^sg$|nitric|allergen|endotoxin|apgar|gestation|birth|pid$|_id$|identifier|study-id|date|_yr$|year|visit|time|duration|hours|blood|urine|serum|plasma|milk-concentration|saliva"
chem_driven <- list("Tobacco" = c("Smoking", "PAH", "TraceElements"), "Combustion/fuel" = c("PAH"), "Traffic/outdoor" = c("PAH"),
                    "Diet" = c("TraceElements", "Phthalates", "Phenols", "PFAS"), "Water" = c("TraceElements", "PFAS"),
                    "Pesticides/pests" = c("Pesticides"), "Personal care" = c("Phthalates", "Parabens", "Phenols"),
                    "Housing" = c("TraceElements", "Phthalates", "FlameRetardants"), "Occupation" = c("Pesticides", "TraceElements", "PAH"))

CW <- read.csv("Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv")
t <- new.env(); load("Input/Harmonized_Datasets/Harmonized_Targeted.RData", envir = t)
cls <- t$Targeted_Data_All %>% mutate(cohort = sub("_.*$", "", Study)) %>% distinct(cohort, Analyte_Code) %>%
  left_join(CW %>% distinct(Analyte_Code = Code, Category), by = "Analyte_Code") %>% filter(!is.na(Category)) %>%
  group_by(cohort) %>% summarise(classes = list(sort(unique(Category))), .groups = "drop")

files <- list.files("Input", pattern = "_Data_All\\.RData$", full.names = TRUE)
inv <- bind_rows(lapply(files, function(f) {
  cohort <- sub("_.*$", "", basename(f)); e <- new.env(); load(f, envir = e); L <- e[[ls(e)[1]]]
  cb <- NULL; for (m in names(L)) if (is.data.frame(L[[m]]) && all(c("code", "value") %in% names(L[[m]]))) cb <- L[[m]]
  bind_rows(lapply(names(L), function(m) {
    d <- L[[m]]; if (!is.data.frame(d) || grepl("Targeted|Codebook|CodeBook|Olink|HILIC|META|Metab|Untarg|Map", m)) return(NULL)
    bind_rows(lapply(names(d), function(v) {
      vl <- tolower(v); if (grepl(excl, vl)) return(NULL)
      hit <- names(dom)[sapply(dom, function(p) grepl(p, vl))]; if (!length(hit)) return(NULL)
      x <- d[[v]]; x <- x[!is.na(x) & as.character(x) != "" & !as.character(x) %in% c("-99", "-9", "-999", "99", "999", "NA")]
      nn <- length(x); nl <- length(unique(x)); if (nn < 30 || nl < 2) return(NULL)
      tb <- sort(table(x), decreasing = TRUE)
      labs <- if (!is.null(cb) && nl <= 12) paste(head(paste0(names(tb), ifelse(is.na(match(names(tb), as.character(cb$code))), "", paste0(" [", cb$value[match(names(tb), as.character(cb$code))], "]")), "=", tb), 6), collapse = "; ")
              else if (nl <= 12) paste(head(paste0(names(tb), "=", tb), 6), collapse = "; ") else paste0("continuous / ", nl, " values")
      tibble(cohort = cohort, table = m, variable = v, domain = paste(hit, collapse = "+"), n_nonmissing = nn, n_levels = nl, values = labs)
    }))
  }))
}))
inv <- inv %>% left_join(cls, by = "cohort") %>% rowwise() %>%
  mutate(primary_domain = strsplit(domain, "\\+")[[1]][1],
         matched_classes = paste(intersect(classes, unlist(chem_driven[strsplit(domain, "\\+")[[1]]])), collapse = ", "),
         testable = matched_classes != "") %>% ungroup() %>% select(-classes)
write.csv(inv, file.path(out, "source_question_inventory.csv"), row.names = FALSE)

summ <- inv %>% filter(testable) %>% group_by(cohort, primary_domain) %>% summarise(n_vars = n(), .groups = "drop") %>%
  pivot_wider(names_from = primary_domain, values_from = n_vars, values_fill = 0)
write.csv(summ, file.path(out, "inventory_by_cohort_domain.csv"), row.names = FALSE)
rec <- inv %>% filter(testable) %>% distinct(cohort, primary_domain, matched_classes) %>%
  separate_rows(matched_classes, sep = ", ") %>% count(primary_domain, matched_classes, name = "n_cohorts") %>% arrange(desc(n_cohorts))
write.csv(rec, file.path(out, "recurring_domain_class_pairs.csv"), row.names = FALSE)
cat("cohorts scanned:", length(files), "| variables kept:", nrow(inv), "| testable:", sum(inv$testable), "\n")
print(as.data.frame(summ)); print(as.data.frame(rec))

pd <- inv %>% filter(testable) %>% count(cohort, primary_domain)
g <- ggplot(pd, aes(primary_domain, cohort, fill = n)) + geom_tile(colour = "white") + geom_text(aes(label = n), size = 2.6) +
  scale_fill_gradient(low = "#fde0c5", high = "#c0392b", name = "Questions") +
  labs(title = "Source questions available in each cohort (matched to a chemical class it measured)", x = NULL, y = NULL) +
  theme_minimal(base_size = 9) + theme(axis.text.x = element_text(angle = 30, hjust = 1))
ggsave(file.path(out, "Fig_Source_Inventory.png"), g, width = 8, height = 0.25 * n_distinct(pd$cohort) + 1.6, dpi = 200)
