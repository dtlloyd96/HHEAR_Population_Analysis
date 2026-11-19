# R4_5_supp_tables.py -- supplementary tables for the geographic-determinants
# Methods ("Determinants of geographic differences"), built from the existing
# R4_5 outputs; nothing is re-analysed. Writes to Output/Hero_Numbers/, where the
# Word builds read supplementary tables:
#   SuppTable_SourceInventory.tsv   table S10  questionnaire source items by cohort
#   SuppTable_RegionalLayers.tsv    table S11  regional data layers
#   SuppTable_LocationLinkage.tsv   table S12  locations, areas, sampling years, layer values
#   DataS1_SourceItems.csv          data S1    every retained questionnaire item
# Run from the project root:  python3 Coding/Paper1_Submission2/Reviewer_Round2/R4_5_supp_tables.py
import csv
from collections import Counter, defaultdict

R = "Output/Reviewer_Round2/R4_5_determinants/"
OUT = "Output/Hero_Numbers/"
CLASS = {"TraceElements": "trace elements", "FlameRetardants": "flame retardants", "PAH": "PAHs",
         "PCB": "PCBs", "VOCs": "VOCs", "PFAS": "PFAS", "Phthalates": "phthalates", "Phenols": "phenols",
         "Parabens": "parabens", "Pesticides": "pesticides", "Smoking": "smoking"}
DOMAINS = ["Tobacco", "Combustion/fuel", "Traffic/outdoor", "Diet", "Water", "Pesticides/pests",
           "Personal care", "Housing", "Occupation"]


def write_tsv(path, rows):
    with open(path, "w") as f:
        for r in rows:
            f.write("\t".join(str(c) for c in r) + "\n")


# ---- table S10 and data S1: questionnaire source items ----------------------
# retained items, less the exclusions stated in the Methods: DISCOVER diet and
# product-use items (unlabelled response codes) and MADRES rice intake (237 of 241 yes)
items = [r for r in csv.DictReader(open(R + "source_inventory/source_question_inventory.csv"))
         if r["testable"] == "TRUE"
         and not (r["cohort"] == "DISCOVER" and r["primary_domain"] in ("Diet", "Personal care"))
         and not (r["cohort"] == "MADRES" and "rice" in r["variable"].lower())]
by = defaultdict(list)
for r in items:
    by[r["cohort"]].append(r)
rows = [["Cohort", "Items", "Items by domain", "Chemical classes paired"]]
for coh in sorted(by, key=str.lower):
    dom = Counter(r["primary_domain"] for r in by[coh])
    cls = sorted({CLASS.get(c.strip(), c.strip()) for r in by[coh] for c in r["matched_classes"].split(";")
                  for c in c.split(",") if c.strip()})
    rows.append([coh, len(by[coh]), "; ".join(f"{d} {dom[d]}" for d in DOMAINS if dom[d]), ", ".join(cls)])
rows.append(["Total", len(items), "; ".join(f"{d} {sum(1 for r in items if r['primary_domain'] == d)}"
                                              for d in DOMAINS if any(r["primary_domain"] == d for r in items)), "–"])
write_tsv(OUT + "SuppTable_SourceInventory.tsv", rows)
with open(OUT + "DataS1_SourceItems.csv", "w", newline="") as f:
    w = csv.writer(f)
    w.writerow(["cohort", "source_table", "variable", "domain", "n_responses", "n_levels", "response_values", "chemical_classes"])
    for r in sorted(items, key=lambda r: (r["cohort"].lower(), r["primary_domain"], r["variable"])):
        w.writerow([r["cohort"], r["table"], r["variable"], r["primary_domain"], r["n_nonmissing"], r["n_levels"],
                    r["values"], r["matched_classes"]])

# ---- table S11: regional data layers (from the extraction scripts) ----------
write_tsv(OUT + "SuppTable_RegionalLayers.tsv", [
    ["Layer", "Source", "Resolution", "Years", "Linkage to locations"],
    ["Ambient PM2.5 (µg/m³)", "ACAG V6.GL.02.04 (58)", "0.1°", "1998–2023, annual",
     "Population-weighted area mean; matched to sampling years"],
    ["Surface NO2 (ppb)", "Anenberg et al. (59)", "1 km", "2015",
     "Population-weighted area mean; single year for all locations"],
    ["Groundwater arsenic", "Podgorski and Berg (60); probability that groundwater As exceeds 10 µg/L", "1 km", "Static",
     "Population-weighted area mean; for drinking-water exposure, weighted by national groundwater reliance"],
    ["Population density (persons/km²) and weights", "GHS-POP R2023A (65)", "30 arc-seconds (~1 km)",
     "1990, 2000, 2010, 2020", "Nearest epoch to each sampling year"],
    ["Polluting cooking fuels (% of population)", "WHO Household Energy Database, SDG 7.1.2 (61)", "National (urban or total)",
     "Annual", "Urban series for city areas, total otherwise; matched to sampling years"],
    ["PFAS in public drinking water (% of systems with a detection)", "US EPA UCMR 3 (62); PFOS, PFOA and four other PFAS",
     "US state", "2013–2015", "State value for US locations; multi-state regions averaged; US-wide value for national cohorts"],
    ["Agricultural insecticide use (kg/ha cropland)", "FAOSTAT Pesticides Use and Land Use (63)", "National", "Annual",
     "Matched to sampling years (nearest available year)"],
    ["Fish and seafood, and rice supply (kg/capita/year)", "FAOSTAT Food Balances, current and historic (63)",
     "National", "Annual", "Matched to sampling years (nearest available year)"],
    ["Year leaded petrol was banned", "Our World in Data, from UN Environment Programme (64)", "National", "Static",
     "Years since ban at the median sampling year"],
])

# ---- table S12: locations, areas, sampling years and layer values -----------
AREA = {"city": "City (30 km circle)", "state": "State or region", "country": "Country",
        "teddy": "Clinical centres (50 km circles)"}
NAME = {"Mexico_City": "Mexico City", "SanFrancisco": "San Francisco", "Los_Angeles": "Los Angeles",
        "New_York": "New York", "New_Hampshire": "New Hampshire", "NorthEast_US": "Northeastern US",
        "Midwestern_US": "Midwestern US", "United_States": "United States (national)",
        "Washington_State": "Washington State", "South_Africa": "South Africa",
        "TEDDY_United_States": "TEDDY United States", "TEDDY_Finland": "TEDDY Finland",
        "TEDDY_Germany": "TEDDY Germany", "TEDDY_Sweden": "TEDDY Sweden", "ZIP_Brazil": "Brazil (ZIP)",
        "ZIP_Guatemala": "Guatemala (ZIP)", "ZIP_Nicaragua": "Nicaragua (ZIP)", "ZIP_Puerto Rico": "Puerto Rico (ZIP)"}
types = {}
try:  # area type per location, from the unit definitions written by R4_5_geospatial_extract.R
    import subprocess
    out = subprocess.run(["Rscript", "-e", "suppressMessages(library(sf));g<-st_read('" + R +
                          "geospatial/location_units.gpkg',quiet=TRUE);write.csv(sf::st_drop_geometry(g),stdout(),row.names=FALSE)"],
                         capture_output=True, text=True, check=True).stdout
    types = {r["loc_id"]: r["type"] for r in csv.DictReader(out.splitlines())}
except Exception as e:
    raise SystemExit(f"could not read location_units.gpkg: {e}")
yw = list(csv.DictReader(open(R + "geospatial/location_study_year_weights.csv")))
coh, pub = defaultdict(set), defaultdict(set)
for r in yw:
    c = r["Study"].replace("_Targeted", "").replace("NHBCS_Cot", "NHBCS")
    coh[r["loc_id"]].add(c)
    if c in ("ATLAA", "ELEMENT", "MEBCO", "SWSD"):   # windows from published cohort descriptions
        pub[r["loc_id"]].add(c)
L = {r["loc_id"]: r for r in csv.DictReader(open(R + "geospatial/location_layers_time_matched.csv"))}


def f(x, d=1):
    try:
        return f"{float(x):,.{d}f}"
    except (TypeError, ValueError):
        return "–"


rows = [["Location", "Area", "Cohorts", "Sampling years", "Samples", "PM2.5 (µg/m³)", "NO2 (ppb)",
         "P(groundwater As >10 µg/L)", "Polluting cooking fuel (%)", "Insecticide (kg/ha)",
         "Fish (kg/capita)", "Rice (kg/capita)", "Leaded petrol ban"]]
for loc in sorted(L, key=lambda k: (["city", "state", "country", "teddy"].index(types.get(k, "country")), NAME.get(k, k))):
    r = L[loc]
    cs = ", ".join(c + ("*" if c in pub[loc] else "") for c in sorted(coh[loc], key=str.lower))
    yrs = r["year_min"] if r["year_min"] == r["year_max"] else f"{r['year_min']}–{r['year_max']}"
    rows.append([NAME.get(loc, loc), AREA[types[loc]], cs, yrs, f(r["n_samples"], 0), f(r["PM25_ugm3"]),
                 f(r["NO2_ppb"]), f(r["GW_As_prob_gt10"], 2), f(r["polluting_fuel_pct"], 0),
                 f(r["insecticide_kg_ha"], 2), f(r["fish_kg_cap"]), f(r["rice_kg_cap"]),
                 r["lead_ban_year"] if r["lead_ban_year"] not in ("", "NA") else "–"])
write_tsv(OUT + "SuppTable_LocationLinkage.tsv", rows)
print("S10:", len(by), "cohorts,", len(items), "items; S12:", len(rows) - 1, "locations")
