# Chirag J Patel
# 12/22/2022
# download nhanes tables using nhanesA

library(nhanesA)
library(tidyverse)
library(survey)
library(progress)
library(getopt)

# spec <- matrix(c(
#   'component', 'c', 1, "character",
#   'year', 'y', 1, "character"
# ), byrow=TRUE, ncol=4);
# opt <- getopt(spec)
opt <- list(component = "LAB", year = NULL)

years <- c(1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017)
year_str <- 'all'
component <- opt$component
if (!is.null(opt$year)) {
  years <- opt$year
  year_str <- as.character(years)
}

download_component_table <- function(componentname, year) {
  cat(sprintf("Processing year: %s\n", year))
  
  table_names <- nhanesTables(componentname, year, namesonly = TRUE)
  
  # Optional: filter table names
  table_names <- table_names[!grepl("PAX|SPXRAW", table_names)]
  
  tabs <- vector("list", length(table_names))
  
  for (i in seq_along(table_names)) {
    table_name <- table_names[i]
    cat(sprintf("  Downloading: %s\n", table_name))
    
    result <- tryCatch(
      nhanes(table_name),
      error = function(e) {
        cat(sprintf("  ⚠️ Failed to download: %s\n", table_name))
        tibble()  # preserve length with an empty tibble
      }
    )
    
    tabs[[i]] <- result
  }
  
  # Set names only after tabs are completely filled
  names(tabs) <- table_names
  return(tabs)
}


download_variable_names <- function(componentname, year) {
  table_names <- nhanesTables(componentname, year) 
  names(table_names) <- c("table_name", "table_description")
  map(table_names$table_name,~nhanesTableVars(componentname, .x, details=TRUE, nchar=50))
}

progress_bar_nh <- function(bar_length) {
  progress_bar$new(clear=FALSE, width=100, total = bar_length, format = "  downloading [:bar] :percent eta: :eta elapsed: :elapsed")
}

## years: 1999-2018
years <- c(1999, 2001, 2003, 2005, 2007, 2009, 2011, 2013, 2015, 2017)
component <- "LAB"

tables <- map(years, ~download_component_table(component, .x))
names(tables) <- years
#dictionary <- map(years, ~download_variable_names(component, .x)) %>% bind_rows() # index maps to each year  
# data dictionary
# years

# dat <- list(
#   tables = tables,
#   #dictionary = dictionary,
#   years = years
# )

names(tables) <- years

NHANES_Exposure_Data <- tables

#save(NHANES_Demo_Data, file = "../Input/NHANES_Demo_Data.RData")




