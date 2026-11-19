read_targeted_analyte_mapping <- function(path = "Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv") {
  readr::read_csv(path, show_col_types = FALSE)
}

targeted_theme_codes <- function(mapping = NULL,
                                 path = "Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv") {
  if (is.null(mapping)) {
    mapping <- read_targeted_analyte_mapping(path)
  }

  filtered_mapping <- dplyr::distinct(
    dplyr::filter(mapping, !is.na(Category), Category != ""),
    Code,
    Category
  )

  lapply(split(filtered_mapping, filtered_mapping$Category), function(df) unique(df$Code))
}

add_targeted_category <- function(data,
                                  code_col = "Outcome",
                                  mapping = NULL,
                                  path = "Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv") {
  if (is.null(mapping)) {
    mapping <- read_targeted_analyte_mapping(path)
  }

  lookup <- dplyr::distinct(dplyr::select(mapping, Code, Category))

  dplyr::left_join(data, lookup, by = stats::setNames("Code", code_col))
}
