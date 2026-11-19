#!/usr/bin/env python3

import csv
import re
from collections import defaultdict
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
HARMONIZED_TARGETED = ROOT / "Input" / "Harmonized_Datasets" / "Harmonized_Targeted.csv"
CURRENT_FULL_TABLE = ROOT / "Output" / "Paper1_Output_Submission2" / "Full_Chem_Table.csv"
CHEM_NAME_MAPPING = ROOT / "HHEAR_Man1_App" / "Chem_Name_Mapping.csv"
OUTPUT_MAPPING = ROOT / "Input" / "Harmonized_Datasets" / "Targeted_Analyte_Category_Mapping.csv"
OUTPUT_REPORT = ROOT / "Output" / "Paper1_Output_Submission2" / "Targeted_Analyte_Currently_Unmapped.csv"


CATEGORY_BY_GROUP = {
    "BDC": "Bisdithiocarbamates",
    "BLIPID": "Lipids",
    "BTE": "TraceElements",
    "INFLAM": "Inflammation",
    "OXID": "OxidativeStress",
    "PBFR": "FlameRetardants",
    "SOC": "Pesticides",
    "SPCB": "PCB",
    "SPFAS": "PFAS",
    "STRESS": "Stress",
    "UAS": "TraceElements",
    "UDILUTE": "UrineDilution",
    "UOPFR": "FlameRetardants",
    "UPAH": "PAH",
    "UPB": "Parabens",
    "UPEST": "Pesticides",
    "UPEST2": "Pesticides",
    "UPEST22": "Pesticides",
    "UPHENOL": "Phenols",
    "UPHTH": "Phthalates",
    "USMOKE": "Smoking",
    "UTE": "TraceElements",
    "VOCs": "VOCs",
}

EXPOSURE_CATEGORIES = {
    "TraceElements",
    "Phenols",
    "Phthalates",
    "PAH",
    "PCB",
    "Parabens",
    "Smoking",
    "Pesticides",
    "VOCs",
    "FlameRetardants",
    "PFAS",
    "Bisdithiocarbamates",
}

CODE_CATEGORY_OVERRIDES = {
    # Ethylenethiourea and propylenethiourea are the urinary metabolites of the
    # ethylene-bis-dithiocarbamate fungicides (maneb, mancozeb, zineb). The
    # Chemical_Group code BDC maps them to "Bisdithiocarbamates", which is not
    # one of the 11 exposure classes the manuscript reports, so the class filter
    # silently dropped them from the exposure count. They are pesticide exposures
    # and belong with the pesticides. (2026-08-10)
    "ETU": "Pesticides",
    "PTU": "Pesticides",
    "BUP": "Parabens",
    "BUPB": "Parabens",
    "BZPB": "Parabens",
    "ETP": "Parabens",
    "ETPB": "Parabens",
    "HEPB": "Parabens",
    "MEPB": "Parabens",
    "PRP": "Parabens",
    "PRPB": "Parabens",
    "MEP": "Phthalates",
    "PHEN4": "PAH",
    "PHEN9": "PAH",
    "PHEN2PHEN3": "PAH",
    # ATLAA dried-blood-spot PFAS (added 2026-08). Serum is the canonical PFAS
    # matrix; DBS is retained in the harmonized data under its own codes but is
    # deliberately given NO category so it is excluded from category/theme-level
    # analyses (targeted_theme_codes() drops rows with an empty Category).
    # Detection in DBS is much lower (PFNA_DBS ~15%, PFOA_DBS ~40%) and it is a
    # different biospecimen, so it should not sit alongside serum in PFAS summaries.
    "PFHXS_DBS": "",
    "PFOS_DBS": "",
    "PFOA_DBS": "",
    "PFNA_DBS": "",
    # NHBCS breast-milk PFAS (added 2026-08). Serum is the canonical PFAS matrix;
    # these are retained under their own _BR codes so they are never pooled with
    # serum, and given NO category so they stay out of category/theme analyses.
    # Milk:serum PFOS ratio is ~1%, so mixing the two would be meaningless.
    "PFBS_BR": "", "PFDA_BR": "", "PFHP_BR": "", "PFHX_BR": "", "PFHp_BR": "",
    "PFNA_BR": "", "PFOA_BR": "", "PFOS_BR": "", "PFUN_BR": "",
}

MANUAL_OVERRIDES = {
    "ACE": {"Full_Name": "Acetochlor mercapturate", "CAS": ""},
    # --- New cohorts added 2026-08 (GRAPHS PAH/VOC, SWSD neonicotinoids) ---
    # Category auto-assigned from Chemical_Group via CATEGORY_BY_GROUP
    # (UPEST2->Pesticides, UPAH->PAH, VOCs->VOCs). NOTE: SWSD acetamiprid uses
    # code ACMP because ACE is already Acetochlor mercapturate (above).
    "ACMP": {"Full_Name": "Acetamiprid", "CAS": "135410-20-7"},
    "DMCLO": {"Full_Name": "N-desmethyl-clothianidin", "CAS": ""},
    "IMI": {"Full_Name": "Imidacloprid", "CAS": "138261-41-3"},
    "NDMA": {"Full_Name": "N-desmethyl-acetamiprid", "CAS": ""},
    "OFIMI": {"Full_Name": "Imidacloprid-olefin", "CAS": ""},
    "FLUO3": {"Full_Name": "3-hydroxyfluorene", "CAS": "6344-67-8"},
    "PHEN3": {"Full_Name": "3-hydroxyphenanthrene", "CAS": "605-87-8"},
    "HEMA2": {"Full_Name": "N-Acetyl-S-(2-hydroxyethyl)-L-cysteine", "CAS": "15060-26-1"},
    "AsB": {"Full_Name": "Arsenobetaine", "CAS": ""},
    "AsC": {"Full_Name": "Arsenocholine", "CAS": ""},
    "BCETP": {"Full_Name": "BCETP", "CAS": "3040-56-0"},
    "BCPP": {"Full_Name": "BCPP", "CAS": "789440-10-4"},
    "BDE1": {"Full_Name": "BDE-1", "CAS": "7025-06-1"},
    "BDE2": {"Full_Name": "BDE-2", "CAS": "6876-00-2"},
    "BDE4": {"Full_Name": "BDE-4", "CAS": "51452-87-0"},
    "BDE6": {"Full_Name": "BDE-6", "CAS": "147217-72-9"},
    "BDE8": {"Full_Name": "BDE-8", "CAS": "147217-71-8"},
    "BDE9": {"Full_Name": "BDE-9", "CAS": "337513-66-3"},
    "Be": {"Full_Name": "Beryllium", "CAS": "7440-41-7"},
    "BUP": {"Full_Name": "Butyl paraben", "CAS": "94-26-8"},
    "CHOL": {"Full_Name": "Cholesterol", "CAS": "57-88-5"},
    "CLO": {"Full_Name": "Clothianidin", "CAS": "210880-92-5"},
    "CORTE": {"Full_Name": "Cortisone", "CAS": "53-06-5"},
    "CORTL": {"Full_Name": "Cortisol", "CAS": "50-23-7"},
    "COT": {"Full_Name": "Cotinine", "CAS": "486-56-6"},
    "CRP": {"Full_Name": "C-reactive protein", "CAS": "9007-41-4"},
    "CYMA": {"Full_Name": "N-Acetyl-S-(2-cyanoethyl)-L-cysteine", "CAS": "74514-75-3"},
    "DBUP": {"Full_Name": "DBUP", "CAS": "107-66-4"},
    "DBZP": {"Full_Name": "DBZP", "CAS": "1623-08-1"},
    "DINF": {"Full_Name": "Dinotefuran", "CAS": "165252-70-0"},
    "DMA": {"Full_Name": "Dimethylarsinic acid", "CAS": "75-60-5"},
    "DOCP": {"Full_Name": "DOCP", "CAS": "27378-80-9"},
    "DPCP": {"Full_Name": "DPCP", "CAS": "843-24-3"},
    "ETP": {"Full_Name": "Ethyl paraben", "CAS": "120-47-8"},
    "ETU": {"Full_Name": "Ethylenethiourea", "CAS": "96-45-7"},
    "F2A12I": {"Full_Name": "F2-isoprostane (12-isoprostane isomer)", "CAS": ""},
    "F2A23D": {"Full_Name": "F2-isoprostane (2,3-dinor metabolite)", "CAS": ""},
    "F2A8IP": {"Full_Name": "8-isoprostane", "CAS": ""},
    "F2AVI": {"Full_Name": "F2-isoprostane VI", "CAS": ""},
    "FLUO": {"Full_Name": "Hydroxyfluorene", "CAS": "2443-58-5_6344-67-8"},
    "GMCSF": {"Full_Name": "Granulocyte macrophage Colony stimulating factor", "CAS": "83869-56-1"},
    "HAPTOGLOBIN": {"Full_Name": "Haptoglobin", "CAS": ""},
    "HCOTT": {"Full_Name": "3'-hydroxycotinine", "CAS": "86-56-6"},
    "HMPMA": {
        "Full_Name": "3-Hydroxypropyl methyl mercapturic acids",
        "CAS": "",
    },
    "HPMA2": {"Full_Name": "N-Acetyl-S-(2-hydroxypropyl)-L-cysteine", "CAS": "923-43-3"},
    "iAs": {"Full_Name": "Inorganic arsenic", "CAS": ""},
    "IMI5OH": {"Full_Name": "5-Hydroxy imidacloprid", "CAS": ""},
    "MMA": {"Full_Name": "Monomethylarsonic acid", "CAS": ""},
    "MBP": {"Full_Name": "mono-n-butylphthalate (BzBP, DnBP)", "CAS": "34-74-2"},
    "MBZ": {"Full_Name": "mono-benzyl phthalate", "CAS": "2528-16-7"},
    "MCI": {"Full_Name": "mono-carboxy isononyl phthalate", "CAS": ""},
    "MCMH": {"Full_Name": "mono-2-(carboxymethyl) hexyl phthalate", "CAS": "82975-93-7"},
    "MCO": {"Full_Name": "mono-carboxy isooctyl phthalate", "CAS": "898544-09-7"},
    "MCOMOP": {"Full_Name": "mono (7-COOH-2-methyloctyl) phthalate", "CAS": ""},
    "MCOP": {"Full_Name": "mono-carboxy octyl phthalate", "CAS": ""},
    "MCP": {"Full_Name": "mono (3-carboxypropyl) phthalate (multiple)", "CAS": "66851-46-5"},
    "MEC": {"Full_Name": "mono-(2-ethyl-5-carboxypentyl) phthalate", "CAS": "40809-41-4"},
    "MEH": {"Full_Name": "mono 2-ethyl hexyl phthalate", "CAS": "4376-20-9"},
    "MEHH": {"Full_Name": "mono (2-ethyl-5-hydroxyhexyl) phthalate", "CAS": "40321-99-1"},
    "MEO": {"Full_Name": "mono (2-ethyl-5-oxohexyl) phthalate", "CAS": "40321-98-0"},
    "MEOH": {"Full_Name": "mono (2-ethyl-5-oxohexyl) phthalate", "CAS": "40321-98-0"},
    "MHN": {
        "Full_Name": "Cyclohexane-1,2-dicarboxylic acid mono hydroxyisononyl ester (DINCH)",
        "CAS": "1637562-52-7",
    },
    "MIB": {"Full_Name": "mono-isobutyl phthalate (DiBP)", "CAS": "30833-53-5"},
    "MIN": {"Full_Name": "mono-isononyl phthalate", "CAS": "106610-61-1"},
    "MNB": {"Full_Name": "mono-n-butylphthalate (BzBP, DnBP)", "CAS": "34-74-2"},
    "NAP": {"Full_Name": "Hydroxynaphthalene", "CAS": "90-15-3_135-19-3"},
    "OS8O": {"Full_Name": "8-OH-deoxyguanosine", "CAS": "88847-89-6"},
    "OSMD": {"Full_Name": "Malonodialdehyde", "CAS": "542-78-9"},
    "PBB1": {"Full_Name": "PBB-1", "CAS": "2052-07-5"},
    "PHE": {"Full_Name": "Hydroxyphenanthrene", "CAS": ""},
    "PFDS": {"Full_Name": "Perfluorodecanesulfonate", "CAS": "126105-34-8"},
    "PFHP": {"Full_Name": "Perfluoroheptanoic acid", "CAS": "375-85-9"},
    "PFHX": {"Full_Name": "Perfluorohexanoic acid", "CAS": "307-24-4"},
    "PFHp": {"Full_Name": "Perfluoroheptanoic acid", "CAS": "375-85-9"},
    "PFUN": {"Full_Name": "Perfluoroundecanoic acid", "CAS": "2058-94-8"},
    "PHEN": {"Full_Name": "Hydroxyphenanthrene", "CAS": ""},
    "PHET": {"Full_Name": "Phenanthrene tetrol", "CAS": ""},
    "PHEN2PHEN3": {"Full_Name": "2-hydroxyphenanthrene + 3-hydroxyphenanthrene", "CAS": "605-55-0_605-87-8"},
    "PRP": {"Full_Name": "Propyl paraben", "CAS": "94-13-3"},
    "PTU": {"Full_Name": "Propylenethiourea", "CAS": "96-45-7"},
    "PYR": {"Full_Name": "Hydroxypyrene", "CAS": "5315-79-7"},
    "Pt": {"Full_Name": "Platinum", "CAS": "7440-06-4"},
    "SG": {"Full_Name": "Specific gravity", "CAS": ""},
    "SUF": {"Full_Name": "Sulfoxaflor", "CAS": "946578-00-3"},
    "TBBA": {"Full_Name": "TBBA", "CAS": "27581-13-1"},
    "TG": {"Full_Name": "Triglycerides", "CAS": ""},
    "TIMP2": {"Full_Name": "tissue inhibitor of metalloproteinase 2", "CAS": "127497-59-0"},
    "THI": {"Full_Name": "Thiamethoxam", "CAS": "153719-23-4"},
    "U": {"Full_Name": "Uranium", "CAS": "7440-61-1"},
    "W": {"Full_Name": "Tungsten", "CAS": "7440-33-7"},

    # --- CAS zero-padding repaired 2026-08-11 ---------------------------------
    # These arrived zero-padded on both the first and last group ("0598-02-07"),
    # the signature of a spreadsheet coercing a CAS number to a date. Stripping
    # the padding makes each one satisfy the CAS check digit, which is what
    # confirms the recovered value rather than a guess. Verified with
    # scratchpad/cas_check.py.
    "DEDP": {"Full_Name": "Diethyl dithiophosphate", "CAS": "298-06-6"},
    "DEP":  {"Full_Name": "Diethyl phosphate", "CAS": "598-02-7"},
    "DMTP": {"Full_Name": "Dimethyl thiophosphate", "CAS": "59401-04-6"},
    "PNP":  {"Full_Name": "4-Nitrophenol", "CAS": "100-02-7"},
}

# --- THE PCB CAS COLUMN WAS SHIFTED BY ONE ROW, repaired 2026-08-11 ----------
# 31 of the 43 PCB rows carried the CAS belonging to the NEXT congener row:
# PCB167 held PCB170's number, PCB170 held PCB172_PCB192's, and so on down the
# block -- the signature of the column being pasted one row out of alignment.
# Every value is now rebuilt from the congener number itself, because PubChem
# resolves "PCB 99" directly. Do NOT repair these by hand from IUPAC
# substitution patterns: an earlier attempt to do that put the wrong number on
# PCB89, and the congener-number lookup caught it.
# Re-run scratchpad/pcb_fix.py after any crosswalk rebuild.

# --- CAS repairs confirmed against PubChem, 2026-08-11 -----------------------
# Every component below satisfies the CAS check digit and was resolved by name
# through the PubChem PUG REST API. PCB99 was already correct: the collisions
# were the OTHER entries carrying its number.
CAS_CONFIRMED = {
    "MHXP":         "24539-57-9",                       # was 245539-57-9, duplicated digit
    "PCB85_PCB120": "65510-45-4_68194-12-7",            # was PCB99's number
    "PCB74_PCB61":  "32690-93-0_33284-53-6",            # was PCB90/101/89's string
    "PCB90_101_89": "68194-07-0_37680-73-2_55215-17-3",  # third component was wrong
    "DCP":          "120-83-2_583-78-8",                # separator normalised to '_'
    "MBP":          "131-70-4",                          # 34-74-2 is a valid alias; this is primary
    "MNB":          "131-70-4",
    "MNBP":         "131-70-4",
}

# --- CAS left deliberately BLANK ---------------------------------------------
# Policy: an unconfirmable CAS is left empty rather than carried as a number
# that cannot be traced to a registry entry. Re-populate only from a lookup.
CAS_BLANK_UNCONFIRMED = {
    "TDCCA": "Stored 55401-03-06. The correct value is in the 55701-xx-x family, but the "
             "cis/trans assignment could not be settled: querying the trans compound returns "
             "55701-03-6 among its synonyms while reverse-querying that CAS returns the cis "
             "structure. TDCCA is specifically the trans isomer.",
    "MCIOP": "Stored 895844-09-7, which fails the check digit. No PubChem match under any "
             "mono-carboxy-isooctyl synonym. NB 'monocarboxyoctyl phthalate' resolves to "
             "1923895-92-4, and this table also carries a separate MCOP code with no CAS -- "
             "the two may be the same DINP metabolite under different vendor names, which "
             "would mean one analyte is split across two codes.",
    "PTU":   "Stored 96-45-7, which PubChem confirms is ethylenethiourea (ETU). The correct "
             "value depends on the isomer: 'propylenethiourea' returns 2055-46-1, "
             "'4-methyl-2-imidazolidinethione' returns 2122-19-2.",
}

SUSPICIOUS_CODE_PATTERN = re.compile(r"^[A-Za-z0-9_]+$")


def read_csv_rows(path, encoding="utf-8"):
    with path.open("r", encoding=encoding, newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv_rows(path, rows, fieldnames):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def looks_like_code(value):
    return bool(SUSPICIOUS_CODE_PATTERN.fullmatch(value)) and len(value) <= 20


def normalize_chem_name_mapping():
    rows = read_csv_rows(CHEM_NAME_MAPPING)
    normalized = []
    for row in rows:
        code = row["Code"].strip()
        full_name = row["Full_Name"].strip()
        cas = row["CAS"].strip()
        if ((" " in code) or ("(" in code) or ("-" in code) or len(code) > 20) and looks_like_code(full_name):
            code, full_name = full_name, code
        normalized.append({"Code": code, "Full_Name": full_name, "CAS": cas})
    return normalized


def collect_target_codes():
    codes = {}
    with HARMONIZED_TARGETED.open("r", encoding="latin1", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            code = row["Analyte_Code"].strip()
            group = row["Chemical_Group"].strip()
            if not code or code == "NA":
                continue
            codes.setdefault(code, set()).add(group)
    return codes


def choose_category(code, groups):
    if code in CODE_CATEGORY_OVERRIDES:
        return CODE_CATEGORY_OVERRIDES[code]

    categories = []
    for group in sorted(groups):
        if group in CATEGORY_BY_GROUP:
            categories.append(CATEGORY_BY_GROUP[group])

    unique_categories = []
    for category in categories:
        if category not in unique_categories:
            unique_categories.append(category)

    if not unique_categories:
        return ""
    if len(unique_categories) == 1:
        return unique_categories[0]

    priority = [
        "Parabens",
        "Phthalates",
        "Phenols",
        "FlameRetardants",
        "Pesticides",
        "TraceElements",
        "PFAS",
        "PAH",
        "Smoking",
        "VOCs",
    ]
    for category in priority:
        if category in unique_categories:
            return category
    return unique_categories[0]


def first_non_empty(*values):
    for value in values:
        if value is not None and str(value).strip():
            return str(value).strip()
    return ""


def build_mapping():
    target_codes = collect_target_codes()
    current_full_table = {
        row["Code"].strip(): row for row in read_csv_rows(CURRENT_FULL_TABLE)
    }
    chem_name_rows = normalize_chem_name_mapping()

    chem_name_lookup = defaultdict(list)
    for row in chem_name_rows:
        chem_name_lookup[row["Code"]].append(row)

    rows = []
    report_rows = []

    for code in sorted(target_codes):
        groups = sorted(group for group in target_codes[code] if group and group != "NA")
        current_row = current_full_table.get(code, {})
        chem_rows = chem_name_lookup.get(code, [])
        chem_row = chem_rows[0] if chem_rows else {}
        manual_row = MANUAL_OVERRIDES.get(code, {})

        full_name = first_non_empty(
            manual_row.get("Full_Name"),
            current_row.get("Full_Name"),
            chem_row.get("Full_Name"),
        )
        cas = first_non_empty(
            manual_row.get("CAS"),
            current_row.get("CAS"),
            chem_row.get("CAS"),
        )
        category = choose_category(code, target_codes[code])

        needs_review = "FALSE"
        if (not current_row and code not in MANUAL_OVERRIDES and not chem_rows) or full_name == code or not cas:
            needs_review = "TRUE"

        name_source = "manual" if manual_row.get("Full_Name") else "current_full_table" if current_row.get("Full_Name") else "chem_name_mapping" if chem_row.get("Full_Name") else "fallback_code"
        cas_source = "manual" if manual_row.get("CAS") else "current_full_table" if current_row.get("CAS") else "chem_name_mapping" if chem_row.get("CAS") else "missing"

        if not full_name:
            full_name = code
            name_source = "fallback_code"
            needs_review = "TRUE"

        row = {
            "Code": code,
            "Full_Name": full_name,
            "CAS": cas,
            "Chemical_Group": "|".join(groups),
            "Category": category,
            "Name_Source": name_source,
            "CAS_Source": cas_source,
            "Needs_Review": needs_review,
        }
        rows.append(row)

        old_category = current_row.get("Category", "").strip()
        if code not in current_full_table or old_category != category:
            if category in EXPOSURE_CATEGORIES:
                report_rows.append(
                    {
                        "Code": code,
                        "Full_Name": full_name,
                        "Chemical_Group": "|".join(groups),
                        "Old_Category": old_category,
                        "New_Category": category,
                        "Reason": "missing_from_current_table" if code not in current_full_table else "category_changed",
                    }
                )

    return rows, report_rows


def main():
    rows, report_rows = build_mapping()
    write_csv_rows(
        OUTPUT_MAPPING,
        rows,
        [
            "Code",
            "Full_Name",
            "CAS",
            "Chemical_Group",
            "Category",
            "Name_Source",
            "CAS_Source",
            "Needs_Review",
        ],
    )
    write_csv_rows(
        OUTPUT_REPORT,
        report_rows,
        ["Code", "Full_Name", "Chemical_Group", "Old_Category", "New_Category", "Reason"],
    )
    print(f"Wrote {len(rows)} rows to {OUTPUT_MAPPING}")
    print(f"Wrote {len(report_rows)} rows to {OUTPUT_REPORT}")


if __name__ == "__main__":
    main()
