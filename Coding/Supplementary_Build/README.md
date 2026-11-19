# Supplementary materials build (step 9)

Builds the supplementary materials Word file in the journal's format from
generated inputs, so the document cannot disagree with the analysis. Run it
through `bash Coding/rebuild_supplement.sh`, which first rasterizes the
figure PDFs (`export_supplementary_figs.sh`), or directly:

```bash
cd Coding/Supplementary_Build
npm install
HHEAR_SCIENCE_TEMPLATE=/path/to/Science_Supplementary_Materials_Word_template.docx node build_science.js
```

| File | Role |
|---|---|
| `build_science.js` | assembles the document: Materials and Methods, Supplementary Text, figs. S1–S22 with captions, tables S1–S12 |
| `science_style.js` | paragraph and run styles taken from the journal template |
| `round2_methods.js`, `round2_supp_text.js` | read the Methods and Supplementary Text (Markdown) |
| `captions.json` | figure captions S1–S22 |
| `csvparse.js` | RFC 4180 CSV parser for the tables (a regex split once dropped columns from table S2) |
| `validate_cas.py`, `cas_lookup.py`, `pcb_fix.py` | CAS number curation for table S2: check digits, PubChem lookup, PCB congener → CAS |
| `package.json`, `package-lock.json` | pin the `docx` package |

## Inputs (not in the repository)

- `Output/Supplementary/figs/S1.png` … `S22.png`, from `export_supplementary_figs.sh`
- `Output/Hero_Numbers/*.csv`, `*.tsv`: the tables (`Hero_Numbers_Supp.Rmd`, `R4_5_supp_tables.py`, curated tables)
- `Output/Paper1_Submission3/Text/`: Methods and Supplementary Text
- the journal's Word template (`HHEAR_SCIENCE_TEMPLATE`)

`HHEAR_ROOT`, `HHEAR_SUPP_FIGS` and `HHEAR_SUPP_OUT` override the root, figure
and output folders. The build stops if an input folder is missing rather than
writing a document with gaps.
