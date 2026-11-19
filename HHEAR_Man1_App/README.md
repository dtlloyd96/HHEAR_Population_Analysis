# HHEAR exposome explorer (Shiny)

Live at <https://chiragjp.shinyapps.io/HHEAR_Exposome_Explore/>.

The app reads a data bundle (`.RData`/`.csv` files in this folder) that is
not in the repository. Build it from the analysis outputs with

```bash
bash Coding/build_hhear_man1_app_bundle.sh
```

as the last step of the pipeline (several figure scripts also write into this
folder, and the builder overwrites them), then run
`shiny::runApp("HHEAR_Man1_App")`.
