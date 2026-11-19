# install_packages.R -- install the R packages used anywhere in this repository.
#   Rscript Coding/install_packages.R
cran <- c(
  "Amelia", "FactoMineR", "MatchIt", "RColorBrewer", "arrow", "broom", "bslib",
  "caret", "data.table", "e1071", "exactextractr", "factoextra", "flextable",
  "geosphere", "getopt", "ggalluvial", "ggcorrplot", "ggdist", "gghalves",
  "ggpubr", "ggrepel", "ggridges", "ggsci", "ggspatial", "ggstance", "ggtext",
  "ggthemes", "glmnet", "glue", "gridExtra", "here", "janitor", "jsonlite",
  "knitr", "lme4", "maps", "maptiles", "metafor", "missMDA", "nhanesA",
  "officer", "paletteer", "patchwork", "plotly", "png", "progress", "quantreg",
  "randomForest", "readxl", "rmarkdown", "rnaturalearth", "rnaturalearthdata",
  "scatterpie", "scico", "sf", "shiny", "shinyWidgets", "stars", "survey",
  "survival", "svglite", "terra", "tidyverse", "xgboost")
missing <- setdiff(cran, rownames(installed.packages()))
if (length(missing)) install.packages(missing)
if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install("ComplexHeatmap")
}
