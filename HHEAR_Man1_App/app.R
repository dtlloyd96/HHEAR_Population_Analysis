# app.R
library(tidyverse)
library(shiny)
library(ggplot2)
library(dplyr)
library(plotly)
library(paletteer)
library(ggstance)
library(rlang)
library(shinyWidgets)
library(rnaturalearth)
library(rnaturalearthdata)
library(sf)
library(bslib)
library(ggridges)
load("Child_Plot_Data_Shiny.RData")
load("Child_Interval_Data_Shiny.RData")
load("MatPlot_Data_Shiny.RData")
load("Mat_Interval_Data_Shiny.RData")
load("Harmonized_SDOH.RData")
load("child_list_Targeted.RData")
load("Child_Reference_Intervals.RData")
load("Child_SDOH_Plot_Dfs.RData")
load("Mat_SDOH_Plot_Dfs.RData")
load("Ridge_Plot_Data_Shiny.RData")
load("Mat_Child_Data.RData")
load("Matrix_Lookup.RData")
Locations <- read_csv("Study_Lat_Longs.csv")
Chem_Name_Mapping <- read_csv("Chem_Name_Mapping.csv")

# US site labels, as they appear in Harmonized_SDOH$Location. Anything not listed
# here is treated as non-US, so a new US cohort must be added or its participants
# are silently pooled into "Worldwide". Missing entries here were why Los Angeles
# (MADRES, PARENTs) and Atlanta (ATLAA) were counted as non-US.
US_LOCATIONS <- c(
  "United_States", "NorthEast_US", "Midwestern_US", "New_Hampshire",
  "Washington_State", "Denver", "Baltimore", "Sacramento", "Southern_California",
  "Syracuse", "Michigan", "California", "SanFrancisco", "Massachusetts",
  "Atlanta", "Boston", "Los_Angeles", "New_York"
)

sdoh_locations <- unique(as.character(SDOH_Data_All$Location))
unclassified_us_candidates <- setdiff(US_LOCATIONS, sdoh_locations)
if (length(unclassified_us_candidates)) {
  message("US_LOCATIONS entries absent from the data: ",
          paste(unclassified_us_candidates, collapse = ", "))
}

Dyad_df <- plot_data_vi2 %>%
  mutate(Analyte_Code_chr = as.character(Analyte_Code)) %>%
  left_join(Matrix_Lookup %>% select(Study, Analyte_Code, Unit),
            by = c("Study", "Analyte_Code_chr" = "Analyte_Code")) %>%
  select(-Analyte_Code_chr) %>%
  left_join(.,Chem_Name_Mapping, by = c("Analyte_Code" = "Code")) %>%
  mutate(Full_Name= ifelse(Analyte_Code == "COTT","Cotinine",Full_Name)) %>%
  mutate(Full_Name= ifelse(Analyte_Code == "HCOTT","H-Cotinine",Full_Name)) %>%
  mutate(Full_Name = ifelse(Analyte_Code == "HPMA","N-Acetyl-S-(3-hydroxypropyl)-L-cysteine",Full_Name))

# Dyad_df %>% 
#   dplyr::filter(Analyte_Code == "PHEN9") %>% 
#   View
#length(unique(Chem_Name_Mapping$Code))

DF_plots <- c(Plot_Dfs_Child,Plot_Dfs_Mat)

# Tissue used to be inferred from the location string ("Kampala|Pabna" -> Blood,
# everything else -> Urine). That was wrong for MEBCO, which is Pabna but urine,
# and for every blood cohort not named in the pattern -- ATLAA, GOALS and NYUWHS
# among them. Matrix_Class from the harmonized data is authoritative and is
# unique within study x analyte, so join it instead. NHANES rows carry no HHEAR
# study, so they fall back to the commonest matrix for that analyte.
analyte_matrix_fallback <- Matrix_Lookup %>%
  filter(!is.na(Matrix_Class)) %>%  # unit-only rows (no matrix recorded) must not win the vote
  count(Analyte_Code, Matrix_Class) %>%
  group_by(Analyte_Code) %>%
  slice_max(n, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(Analyte_Code, Fallback_Matrix = Matrix_Class)

add_tissue <- function(df, analyte_col = "curr_analyte") {
  out <- df %>%
    left_join(Matrix_Lookup, by = c("Study" = "Study", stats::setNames("Analyte_Code", analyte_col))) %>%
    left_join(analyte_matrix_fallback, by = stats::setNames("Analyte_Code", analyte_col)) %>%
    mutate(Tissue = dplyr::coalesce(Matrix_Class, Fallback_Matrix, "Unknown")) %>%
    select(-Fallback_Matrix) %>%
    # Units (2026-11-19): Matrix_Lookup carries the display unit per study x
    # analyte. NHANES has no HHEAR study row; its urine is not creatinine-adjusted.
    mutate(
      Unit = dplyr::case_when(
        !is.na(Unit) ~ Unit,
        Study == "NHANES" & Tissue == "Urine" ~ "ng/mL, unadjusted",
        Study == "NHANES" ~ "ng/mL",
        TRUE ~ "unit not recorded"
      ),
      Measure = paste0(Tissue, ": ", Unit)
    )
  n_unknown <- sum(out$Tissue == "Unknown")
  if (n_unknown > 0) message("Rows with no specimen matrix: ", n_unknown, " of ", nrow(out))
  out
}

plot_data <- add_tissue(Ridge_Plot_Data)

# Axis-title unit for whatever is on screen: one unit -> show it; several (urine
# cohorts differ in dilution handling) -> point to the legend, which shows
# Measure = matrix + unit per point.
unit_suffix <- function(units) {
  u <- sort(unique(stats::na.omit(as.character(units))))
  if (length(u) == 1) paste0(" (", u, ")") else ""
}
# ggplot's default shape scale stops at 6 values and silently drops the rest; a
# urinary analyte can span more matrix/unit combinations than that.
measure_shapes <- scale_shape_manual(values = c(16, 17, 15, 3, 4, 8, 1, 2, 0, 5, 6, 18, 7, 9, 10, 12, 13, 14, 11, 19, 20, 21, 22, 23, 24, 25))

# "Units in this view" table shown under each plot (2026-11-19 polish pass).
# One row per study (and location where there is one); NHANES last.
units_table <- function(df) {
  cols <- intersect(c("Study", "Location", "Tissue", "Unit"), names(df))
  if (!nrow(df) || !all(c("Study", "Unit") %in% cols)) return(NULL)
  u <- df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(cols), as.character)) %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(cols))) %>%
    dplyr::arrange(Study == "NHANES", Study)
  header <- c(Study = "Study", Location = "Location", Tissue = "Matrix", Unit = "Unit")[cols]
  mixed <- dplyr::n_distinct(u$Unit) > 1
  tags$details(
    class = "mt-2", open = if (mixed) NA else NULL,
    tags$summary(class = "text-muted small",
                 if (mixed) "Units differ between cohorts in this view; compare urinary values with care."
                 else paste("Units in this view:", u$Unit[1])),
    tags$table(
      class = "table table-sm small mt-1", style = "max-width: 720px;",
      tags$thead(tags$tr(lapply(header, tags$th))),
      tags$tbody(lapply(seq_len(nrow(u)), function(i) tags$tr(lapply(cols, function(cn) tags$td(u[[cn]][i])))))
    )
  )
}

conc_axis_title <- function(units, log1p = FALSE) {
  paste0(if (isTRUE(log1p)) "log1p(Concentration)" else "Concentration", unit_suffix(units))
}

mypal <- paletteer_d("ggsci::springfield_simpsons")

# sort() added 2026-08-12: this drives the study picker, and it took CSV row
# order, so the seven cohorts appended to Study_Lat_Longs.csv (the three new
# ones plus GOALS/NYUWHS/PROGRESS/THS) landed after Z instead of in sequence.
mystudies <- sort(unique(Locations$Study))

ExWAS_SDOH_Exposures <- SDOH_Data_All %>% 
  mutate(Edu_Any = case_when(!is.na(Max_Edu) ~ Max_Edu,
                             is.na(Max_Edu) ~ Mat_Edu,
                             !is.na(Mat_Edu) ~ Mat_Edu)) %>% 
  mutate(Mat_Child_Comb_Race = case_when(!is.na(Child_Race) ~ Child_Race,
                                         is.na(Child_Race) ~ Mat_Race,
                                         !is.na(Mat_Race) ~ Mat_Race)) %>% 
  mutate(Mat_Age = as.numeric(Mat_Age)) %>% 
  dplyr::mutate(Location = relevel(as.factor(Location),ref = "United_States")) %>% 
  mutate(Location_Country = ifelse(Location %in% US_LOCATIONS,
                                   "United_States", "Worldwide")) %>%
  mutate(Mat_Edu_Binary = case_when(
    Mat_Edu %in% c("No_Education",
                   "Primary_Education",
                   "Primary_School",
                   "High_School") ~ "No_College",
    Mat_Edu %in% c("Some_College",
                   "Four_Year_College_Degree",
                   "Some_Grad_School") ~ "Some_College",
    Mat_Edu == "Not_Reported" ~ NA
  )) %>% 
  mutate(Mat_Income_Binary = case_when(
    Mat_Income %in% c("0_25k", "25_49k") ~ "Under_49k",
    Mat_Income %in% c("50_74k", "75_99k", "100_124k", 
                      "Over_100k", "Over_125k", "Over_49k") ~ "Over_49k",
    TRUE ~ NA
  )) %>% 
  mutate(Study =  str_remove(Study, '[_]\\w+|:')) %>% 
  mutate(Mat_Age = as.numeric(Mat_Age))

# ExWAS_SDOH_Exposures %>% 
#   dplyr::filter(Study == "CTS") %>% 
#   View
ExWAS_SDOH_Exposures$Mat_Edu <- ExWAS_SDOH_Exposures$Mat_Edu %>%
  fct_recode(
    "Primary"       = "Primary_School",
    "Primary"       = "Primary_Education",
    "High School"   = "High_School",
    "Some College"  = "Some_College",
    "College"       = "Four_Year_College_Degree",
    "Grad School"   = "Some_Grad_School",
    "Other"         = "Other",
    "None"          = "No_Education"
  ) %>%
  fct_relevel("None", "Primary", "High School", "Some College",
              "College", "Grad School", "Other")



Ref_Vars <- c("Mat_Age","Mat_Edu","Mat_Income", "Max_Edu","Child_Age","Mat_Race","Child_Race","Child_Sex")



Child_SDOH <- ExWAS_SDOH_Exposures %>% 
  dplyr::filter(!is.na(Child_PID) ) %>% 
  dplyr::select(Study,Child_PID,Mat_PID,all_of(Ref_Vars)) 


Mat_SDOH <- ExWAS_SDOH_Exposures %>% 
  dplyr::filter(!is.na(Mat_PID) ) %>% 
  dplyr::select(Study,Child_PID,Mat_PID,all_of(Ref_Vars)) 

ExWAS_SDOH_Exposures <- ExWAS_SDOH_Exposures %>% 
  dplyr::select(Study,Child_PID,Mat_PID,all_of(Ref_Vars)) 

Location_Plotdata <- bind_rows(Location_Plotdata_Child %>% 
                                 mutate(Population = "Child"),
                               Location_Plotdata_Mat%>% 
                                 mutate(Population = "Adult")) %>% 
  #left_join(.,Chem_Name_Mapping, by = c("curr_analyte" = "Code")) %>%
  add_tissue()


interval_units <- Location_Plotdata %>%
  dplyr::filter(Study != "NHANES") %>%
  dplyr::distinct(Population, curr_analyte, Unit) %>%
  dplyr::group_by(Population, curr_analyte) %>%
  dplyr::summarise(Units_Contributing = paste(sort(unique(Unit)), collapse = "; "), .groups = "drop")

Location_interval_data <- bind_rows(Location_interval_data_Child %>% 
                                 mutate(Population = "Child"),
                                 Location_interval_data_Mat%>% 
                                   mutate(Population = "Adult"))

Location_LOD_Data <- bind_rows(
  get0("Location_LOD_Data_Child", ifnotfound = tibble()) %>% mutate(Population = "Child"),
  get0("Location_LOD_Data_Mat", ifnotfound = tibble()) %>% mutate(Population = "Adult")
)

mypal <- c("#FED439FF" ,'#709AE1FF', '#8A9197FF' ,'#D2AF81FF','#FD7446FF', '#D5E4A2FF', '#197EC0FF')
location_colors <- c(
  # --- National anchors ---
  "NHANES-USA"      = "#C51B7D",
  "United_States"   = "#3C3C3C",
  
  # --- U.S. East Coast / Northeast (cool deep blues) ---
  "Boston"          = "#08306B",
  "Northeastern US" = "#08458A",
  "Massachusetts"   = "#08519C",
  "New_Hampshire"   = "#2171B5",
  "New_York"        = "#2171C1",
  "Baltimore"       = "#6BAED6",
  "Michigan"        = "#C6DBEF",  
  "Syracuse"        = "#9ECAE1",
  
  # --- U.S. West Coast (distinct teals/greens) ---
  "California"       = "#00441B",
  "Los_Angeles"      = "#1B7837",
  "Sacramento"       = "#5AAE61",
  "SanFrancisco" =     "#5AAE79",
  "Yakima Valley,Washington" = "#A6DBA0",
  
  # --- U.S. Mountain/Plains/South (olive) ---
  "Denver"           = "#B8DE29",
  "Midwestern_US"    = "#8FB31F",
  "Atlanta"          = "#6B8E23",
  
  # --- Caribbean (aqua) ---
  "Northern Karst,Puerto Rico"      = "#00A7C2",
  
  # --- Latin America (warm reds/oranges) ---
  "Mexico_City"      = "#D7301F",
  "Suchitepéquez,Guatemala"        = "#FC8D59",
  "Managua,Nicaragua"        = "#FDBB84",
  "Suriname"         = "#FDD49E",
  "Bahia, Brazil"           = "#BD0026",
  "Santiago,Chile"            = "#F46D43",
  "Pedro Moncayo, Ecuador"          = "#FDAE61",
  "Tobago,Trinidad and Tobago"           = "#FDAA68",
  
  # --- Europe (greens) ---
  "Turku,Finland" = "#B197FC",  # lavender violet
  "Munich,Germany" = "#C4A7E7",  # muted lilac
  "Malmö,Sweden"  = "#D7BCE8" ,  # pale orchid
  
  
  # --- Africa (earth tones) ---
  "Cape Town,South Africa"     = "#8C510A",
  "Kampala,Uganda"           = "#BF812D",
  "Accra,Ghana"              = "#DFC27D",

  # --- Southeast Asia ---
  "Chiang Mai,Thailand"      = "#E7298A",
  
  # --- South Asia ---
  "Pabna,Bangladesh"       = "#A50F15",
  
  # --- Missing ---
  "NA"               = "#9E9E9E"
)

preferred_analyte_code <- function(codes, preferred = "As") {
  if (preferred %in% codes) preferred else if (length(codes)) codes[1] else NULL
}

preferred_analyte_name <- function(names_vec, preferred = "Arsenic") {
  if (preferred %in% names_vec) preferred else if (length(names_vec)) names_vec[1] else NULL
}



# ui <- fluidPage(
#   titlePanel("Exposome Explorer"),
#   tabsetPanel(
#     id = "tabs",
#     tabPanel(
#       title = "Map • Study Locations",
#       sidebarLayout(
#         sidebarPanel(
#           helpText("Map of public HHEAR study locations."),
#           checkboxInput("show_world", "Show world outline", TRUE),
#           checkboxInput("show_points", "Show study sites", TRUE)
#         ),
#         mainPanel(
#           plotlyOutput("map_plot", height = "600px", width = "100%")
#         )
#       )
#     ),
#     tabPanel(
#       title = "Profiles • By Study",
#       sidebarLayout(
#         sidebarPanel(
#           selectInput("study_pick", "Choose study:", choices = mystudies, selected = mystudies[1]),
#           hr(),
#           checkboxInput("drop_all_na", "Drop columns that are all NA", TRUE),
#           checkboxInput("drop_constant", "Drop columns with a single unique value", TRUE),
#           sliderInput("min_non_na", "Min. non-NA fraction to keep a column",
#                       min = 0, max = 1, value = 0, step = 0.05),
#           sliderInput("max_levels", "Max categories to plot (categorical)",
#                       min = 2, max = 50, value = 20, step = 1),
#           checkboxInput("use_plotly_profiles", "Interactive (Plotly)", FALSE)
#         ),
#         mainPanel(
#           h3("Categorical columns"),
#           uiOutput("cat_plots_ui"),
#           hr(),
#           h3("Numeric columns"),
#           uiOutput("num_plots_ui")
#         )
#       )
#     ),
#   tabPanel(
#       title = "Location Intervals",
#       sidebarLayout(
#         sidebarPanel(
#           selectInput("population", "Population:", choices = c("Child", "Adult"), selected = "Child"),
#           uiOutput("analyte_picker"),
#           hr(),
#           checkboxInput("free_y", "Free Y scale (by analyte)", FALSE),
#           checkboxInput("use_log1p", "Use log1p scale (x-axis)", TRUE),
#           sliderInput("xmax", "X-axis max", min = 2, max = 12, value = 7, step = 0.5),
#           helpText("Tip: zoom & pan with your mouse; double-click to reset.")
#         ),
#         mainPanel(
#           plotlyOutput("child_plot", height = "400px")
#         )
#       )
#     ),
#     # ---- New Quantile Explorer Tab ----
#     tabPanel(
#       title = "Demographic Intervals",
#       sidebarLayout(
#         sidebarPanel(
#           selectInput("df_choice", "Choose dataset:", choices = names(DF_plots)),
#           uiOutput("group_picker"),
#           uiOutput("outcome_picker"),
#           helpText("Select dataset and grouping variable to explore quantiles.")
#         ),
#         mainPanel(
#           plotlyOutput("quantile_plot", height = "400px")
#         )
#       )
#     ),
#     
#     tabPanel("Summary", tags$p("Add tables, KPIs, or overview plots here.")),
#     tabPanel("About",   tags$p("Brief description of the app, data sources, and methods."))
#   )
# )


ui <- fluidPage(
  theme = bs_theme(
    version   = 5,
    bootswatch = "flatly",          # try "cosmo", "minty", "yeti" if you prefer
    base_font = font_google("Inter"),
    heading_font = font_google("Inter"),
    code_font = font_google("JetBrains Mono")
  ),
  
  tags$head(
    # Subtle global styling / spacing
    tags$style(HTML("
      .container-fluid { max-width: 1400px; }
      .sidebar-panel { background: #ffffff; border-radius: 14px; padding: 16px 16px 6px; box-shadow: 0 6px 18px rgba(0,0,0,0.06); }
      .main-card { background: #ffffff; border-radius: 14px; padding: 18px; box-shadow: 0 6px 18px rgba(0,0,0,0.06); }
      .tab-content { padding-top: 12px; }
      .helptext, .form-text, .text-muted { color: #6c757d !important; }
      h2, h3 { letter-spacing: 0.2px; }
      .split-title { display:flex; align-items:center; gap:12px; }
      .split-title .subtitle { font-size: 0.95rem; color:#6c757d; margin-top:2px; }
      .hr-soft { border: 0; height: 1px; background: linear-gradient(90deg, rgba(0,0,0,0), rgba(0,0,0,.08), rgba(0,0,0,0)); margin: 14px 0; }
    "))
  ),
  
  # -------- App Header --------
  div(class = "my-3",
      div(class = "split-title",
          tags$h2(class = "mb-0", "HHEAR Explorer"),
          tags$span(class = "subtitle", "Interactive views of HHEAR datasets")
      ),
      tags$hr(class = "hr-soft")
  ),
  
  tabsetPanel(
    id = "tabs",
    
    tabPanel(
      title = "Landing Page",
      fluidRow(
        column(
          width = 12,
          tags$div(
            style = "max-width: 820px; margin: 24px auto;",
            
            tags$h2(style = "font-weight: 700;", "HHEAR Chemical Exposome Explorer"),
            
            tags$p(
              class = "lead",
              "This application allows interactive exploration of chemical exposure 
           distributions and demographic-specific exposure intervals across HHEAR studies."
            ),
            
            tags$hr(),
            
            tags$h4("What this app does"),
            tags$p(
              "Use this platform to visualize exposure variability, compare exposure 
           intervals across demographic groups, and explore how exposures differ 
           by location or population."
            ),
            
            tags$h4("How to use it"),
            tags$ul(
              tags$li("Select a population (Child/Adult)"),
              tags$li("Choose a chemical feature of the exposome to explore."),
              tags$li("Navigate between tabs to view intervals, distributions, and summaries.")
            ),
            
            tags$h4("Interpreting the plots"),
            tags$p(
              # Units paragraph added 2026-11-19 (DRAFT -- author to edit).
              "Exposure Interval plots display the 25–75% predicted exposure range and 5–95% 
           predicted range, with the predicted median values highlighted. Each point represents an individual HHEAR participant at 
           that location. For HHEAR studies, all individual values are displayed; for NHANES, a random subset of 500 individuals is shown.
           Distribution plots 
           show the full spread of exposure measurements within each group for individuals matched on age and sex."
            ),
            tags$h4("Units"),
            tags$p(
              "Blood, serum and plasma concentrations are shown in ng/mL. Urinary concentrations follow each
           cohort's dilution handling: creatinine-corrected cohorts are in µg/g creatinine, others are
           specific-gravity-adjusted or unadjusted ng/mL, and NHANES urine is unadjusted. Point shapes mark the specimen matrix;
           hover text and the table under each plot give every cohort's unit. Compare urinary values from
           cohorts on different scales with care."
            ),
            tags$h4("Data Access"),
            tags$p(
              "Underlying data can be downloaded from",
              tags$a("Human Health Exposure Analysis Resource (HHEAR) Website.",
                     href = "https://hheardatacenter.mssm.edu/",
                     target = "_blank"),
              "While exposure interval data can be downloaded from the data download tab on this application.",
              "Code can be downloaded from",
              tags$a("github.",
                     href = "https://github.com/dtlloyd96/HHEAR_Population_Analysis",
                     target = "_blank"),
              # DRAFT 2026-11-19 -- author to edit; replace with the citation once accepted.
              tags$br(),
              "Manuscript: Population structure of the exposome shapes epidemiological inference",
              "(Science, under revision; citation to be added on acceptance)."
              
            ),
            tags$hr(),
            
            tags$p(
              class = "text-muted",
              tags$small(
                "Note: Some analytes may not be available for all populations. 
             Sparse strata may be hidden depending on missingness thresholds."
              )
            )
          )
        )
      )
    ),
    
    # ---- Map • Study Locations ----
    tabPanel(
      title = "HHEAR Study Locations",
      sidebarLayout(
        sidebarPanel(class = "sidebar-panel",
                     helpText("Map of public HHEAR study locations."),
                     checkboxInput("show_world", "Show world outline", TRUE),
                     checkboxInput("show_points", "Show study sites", TRUE),
                     tags$small(class = "text-muted",
                                "Tip: pan/zoom in the plot; reset with double-click.")
        ),
        mainPanel(
          div(class = "main-card",
              tags$h3(class = "mb-3", "Study Locations"),
              plotlyOutput("map_plot", height = "600px", width = "100%")
          )
        )
      )
    ),
    
    # ---- Profiles • By Study ----
    tabPanel(
      title = "Demographic Profiles By Study",
      sidebarLayout(
        sidebarPanel(class = "sidebar-panel",
                     selectInput("study_pick", "Choose study:",
                                 choices = mystudies, selected = mystudies[1]),
                     tags$hr(class = "hr-soft"),
                     checkboxInput("drop_all_na", "Drop columns that are all NA", TRUE),
                     checkboxInput("drop_constant", "Drop columns with a single unique value", TRUE),
                     sliderInput("min_non_na", "Min. non-NA fraction to keep a column",
                                 min = 0, max = 1, value = 0, step = 0.05),
                     sliderInput("max_levels", "Max categories to plot (categorical)",
                                 min = 2, max = 50, value = 20, step = 1),
                     checkboxInput("use_plotly_profiles", "Interactive (Plotly)", FALSE),
                     tags$small(class = "text-muted",
                                "Use filters to control which variables are summarized.")
        ),
        mainPanel(
          div(class = "main-card",
              tags$h3("Categorical columns"),
              uiOutput("cat_plots_ui"),
              tags$hr(class = "hr-soft"),
              tags$h3("Numeric columns"),
              uiOutput("num_plots_ui"),
              tags$p(
                class = "mt-3 text-muted",
                "These plots are bar plots and box plots of the categorical and numeric demographic features by study "
              )
          )
        )
      )
    ),
    
    # ---- Location Intervals ----
    tabPanel(
      title = "By Location Exposure Intervals",
      sidebarLayout(
        sidebarPanel(class = "sidebar-panel",
                     selectInput("population", "Population:",
                                 choices = c("Child", "Adult"), selected = "Child"),
                     uiOutput("analyte_picker"),
                     tags$hr(class = "hr-soft"),
                     #checkboxInput("free_y", "Free Y scale (by analyte)", FALSE),
                     checkboxInput("use_log1p", "Use log1p scale (x-axis)", TRUE),
                     #sliderInput("xmax", "X-axis max", min = 2, max = 12, value = 7, step = 0.5),
                     helpText("Zoom & pan with your mouse; double-click to reset.")
        ),
        mainPanel(
          div(class = "main-card",
              tags$h3(class = "mb-3", "Interval Distributions by Location"),
              plotlyOutput("child_plot", height = "560px"),
              tags$p(
                class = "mt-3 text-muted",
                "Points are individual exposure levels, coloured by location (nearest city). Dark bars show the
                pooled predicted 25–75% interval, lighter bars the predicted 5–95% interval, and the large dot the
                predicted median. Open circles mark each cohort's limit of detection."
              ),
              uiOutput("loc_units")
          )
        )
      )
    ),
    
    tabPanel(
      title = "QC By Study Intervals",
      sidebarLayout(
        sidebarPanel(class = "sidebar-panel",
                     selectInput("qc_population", "Population:",
                                 choices = c("Child", "Adult"), selected = "Child"),
                     uiOutput("qc_analyte_picker"),
                     tags$hr(class = "hr-soft"),
                     checkboxInput("qc_use_log1p", "Use log1p scale (x-axis)", TRUE),
                     helpText("Quick QC view using observed within-study quantiles and study-median LOD.")
        ),
        mainPanel(
          div(class = "main-card",
              tags$h3(class = "mb-3", "Observed Intervals by Study"),
              plotlyOutput("study_qc_plot", height = "500px"),
              tags$p(
                class = "mt-3 text-muted",
                "Points are individual observations by study. Dark bars show observed 25–75% ranges, lighter bars show observed 5–95% ranges, the filled dot is the observed median, and the open circle is the study median LOD when available."
              ),
              uiOutput("qc_units")
          )
        )
      )
    ),
    
    # ---- Demographic Intervals ----
    tabPanel(
      title = "By Demographic Exposure Intervals",
      sidebarLayout(
        sidebarPanel(class = "sidebar-panel",
                     uiOutput("study_picker"),
                     uiOutput("group_picker"),
                     uiOutput("outcome_picker"),
                     helpText("Select study and demographic grouping to explore quantiles.")
        ),
        mainPanel(
          div(class = "main-card",
              tags$h3(class = "mb-3", "Exposure Intervals"),
              plotlyOutput("quantile_plot", height = "400px"),
              tags$p(
                class = "mt-3 text-muted",
                "Points are the log concentration of individual exposure levels, colored by demographic feature.
                Dark bars show the predicted 25–75% interval, lighter bars show the predicted 5–95% interval, 
                and large dot represent the predicted median exposure.
                Grouping column is the demographic feature the  exposure concentrations are visualized by."
              )
          )
        )
      )
    ),
    
    tabPanel(
      title = "Matched Strata Concentrations",
      sidebarLayout(
        sidebarPanel(class = "sidebar-panel",
                     uiOutput("ridge_analyte_picker"),
                     helpText("Compares locations within shared age/sex strata, so differences are not driven by who was sampled.")
        ),
        mainPanel(
          div(class = "main-card",
              tags$h3(class = "mb-3", "Exposure Concentrations Across Matched Strata"),
              plotOutput("conc_ridge_plot", height = "640px"),
              tags$p(
                class = "mt-3 text-muted",
                "Each ridge is the concentration distribution in one age (and, for children, sex) stratum, coloured by
                location (nearest city). Strata are on the y axis; the x axis is log1p concentration. Line type marks the
                specimen matrix."
              ),
              uiOutput("ridge_units")
              #plotlyOutput("conc_ridge_plot", height = "500px")
          )
        )
      )
    ),
    
   
    # Mat Child Dyads 
    tabPanel(
      title = "Maternal Child Dyads",
      sidebarLayout(
        sidebarPanel(
          class = "sidebar-panel",
          
          selectizeInput(
            "dyad_analyte",
            "Choose analyte:",
            choices = NULL,
            multiple = FALSE,
            options = list(placeholder = "Start typing an analyte…")
          ),
          
          checkboxInput("dyad_log10", "Log10 scale", TRUE),
          checkboxInput("dyad_show_1to1", "Show 1:1 line", TRUE),
          
          tags$hr(),
          tags$small(
            class = "text-muted",
            "Each point represents a mother–child dyad for the selected analyte."
          )
        ),
        
        # --- UI: swap your mainPanel card to include TWO plots (scatter + violin) ---
        mainPanel(
          div(
            class = "main-card",
            tags$h3(class = "mb-3", "Maternal vs Child Concentrations"),
            plotlyOutput("dyad_scatter", height = "520px"),
            tags$div(style = "margin-top: 6px;", textOutput("dyad_note")),
            uiOutput("dyad_units")
            
            # tags$hr(style = "margin: 12px 0;"),
            # 
            # tags$h4(class = "mb-2", "Child–Mother Log Ratios Across Exposures"),
            # plotlyOutput("dyad_violin", height = "520px")
          )
        )
      )
    ),
    tabPanel(
      title = "Downloads",
      value = "downloads",
      
      sidebarLayout(
        sidebarPanel(
          tags$h4("Download data"),
          tags$p(class = "text-muted",
                 "Choose a dataset and download it as a CSV file."),
          
          selectInput(
            "download_dataset",
            "Choose dataset:",
            choices = c(
              #"Exposure plot data"     = "Location_Plotdata",
              "Exposure interval data" = "Location_interval_data"
            )
          ),
          
          #checkboxInput("download_filtered", "Download filtered (current selections)", TRUE),
          
          downloadButton("download_csv", "Download CSV"),
          
          tags$hr(),
          tags$small(class = "text-muted",
                     "Download HHEAR interval data")
        ),
        
        mainPanel(
          tags$h3("Preview"),
          tableOutput("download_preview")
        )
      )
    ),

    # ---- Summary ----
    # Summary tab hidden 2026-11-19 (polish pass): it only held placeholders
    # ("Placeholder", "Paper Information"). Restore by uncommenting once it has content.
    # tabPanel(
    #   "Summary",
    #   div(class = "main-card",
    #       tags$h3("Overview"),
    #       tags$p(""),
    #       tags$hr(class = "hr-soft"),
    #       fluidRow(
    #         column(
    #           width = 6,
    #           div(class = "main-card",
    #               style = "box-shadow:none; background:#f8f9fa;",
    #               tags$h4("Key Metrics"),
    #               tags$p(class = "text-muted", "Placeholder ")
    #           )
    #         ),
    #         column(
    #           width = 6,
    #           div(class = "main-card",
    #               style = "box-shadow:none; background:#f8f9fa;",
    #               tags$h4("Notes"),
    #               tags$p(class = "text-muted", "Paper Information")
    #           )
    #         )
    #       )
    #   )
    # ),

    # ---- About ----
    tabPanel(
      "About",
      div(class = "main-card",
          # DRAFT 2026-11-19 -- author to edit.
          tags$h3("About this app"),
          tags$p("The HHEAR Chemical Exposome Explorer accompanies the manuscript
                  \"Population structure of the exposome shapes epidemiological inference\"
                  (Science, under revision). It shows targeted chemical exposure measurements
                  from 38 cohorts of the Human Health Exposure Analysis Resource (HHEAR),
                  analysed separately for children and for adults (mothers in birth cohorts)."),
          tags$h4("Data"),
          tags$ul(
            tags$li("Targeted analyte concentrations harmonized across HHEAR cohorts and converted to ng/mL wherever the reported unit allows; mass-per-mass (ng/g) and molar (pmol/mL) values keep their own units."),
            tags$li("Blood, serum and plasma are not dilution-adjusted. Urine is creatinine-corrected (µg/g creatinine), specific-gravity-adjusted or unadjusted, depending on the cohort; the unit of each value is shown in the plots and downloads."),
            tags$li("Each cohort is placed at its sampling location; TEDDY and ZIP participants are assigned to their country sites."),
            tags$li("NHANES (all available survey cycles) is shown for comparison; a random subset of 500 participants is plotted per chemical. NHANES urine is not creatinine-adjusted.")
          ),
          tags$h4("Methods"),
          tags$ul(
            tags$li("Exposure intervals are the 5th–95th and 25th–75th percentiles, with the median, estimated by censored quantile regression adjusted for age and sex; values below the limit of detection are treated as censored, not imputed."),
            tags$li("Intervals are estimated for the pooled HHEAR population and by location, study and demographic group, separately for each specimen matrix."),
            tags$li("Matched-strata plots compare locations within shared age (and, for children, sex) strata.")
          ),
          tags$h4("Data access"),
          tags$p(HTML(paste0(
            'Participant-level data are available from the ',
            '<a href="https://hheardatacenter.mssm.edu/" target="_blank">HHEAR Data Center</a>; ',
            'analysis code is on ',
            '<a href="https://github.com/dtlloyd96/HHEAR_Population_Analysis" target="_blank">GitHub</a>.')))
      )
    )
  )
)
server <- function(input, output, session) {
  
  ##########DEMOS#################
  
  safe_renderPlot <- function(expr, msg = "Plot not available for this variable.") {
    renderPlot({
      out <- tryCatch(
        expr,
        error = function(e) e
      )
      
      if (inherits(out, "error")) {
        validate(need(FALSE, msg))
        return(invisible(NULL))
      }
      
      out
    })
  }
  
  safe_renderPlotly <- function(expr, msg = "Plot not available for this variable.") {
    plotly::renderPlotly({
      out <- tryCatch(
        expr,
        error = function(e) e
      )
      
      if (inherits(out, "error")) {
        validate(need(FALSE, msg))
        return(NULL)
      }
      
      out
    })
  }
  
  sanitize_id <- function(x) gsub("[^A-Za-z0-9_]", "_", x)
  
  study_df_re <- reactive({
    ExWAS_SDOH_Exposures %>% 
    dplyr::filter(Study == input$study_pick) %>% 
      dplyr::select(-Mat_PID,-Child_PID)
  })
  
  filtered_df <- reactive({
    df <- study_df_re()
    req(nrow(df) > 0)
    
    non_na_frac <- sapply(df, function(x) mean(!is.na(x)))
    keep <- rep(TRUE, ncol(df))
    if (isTRUE(input$drop_all_na)) keep <- keep & (non_na_frac > 0)
    if (isTRUE(input$drop_constant)) {
      const_col <- sapply(df, function(x) dplyr::n_distinct(x, na.rm = TRUE) <= 1)
      keep <- keep & !const_col
    }
    keep <- keep & (non_na_frac >= input$min_non_na)
    
    cols_keep <- names(df)[keep]
    df_kept <- dplyr::select(df, tidyselect::all_of(cols_keep))
    
    if (ncol(df_kept) > 0) {
      df_kept <- df_kept[rowSums(is.na(df_kept)) < ncol(df_kept), , drop = FALSE]
    }
    
    df_kept
  })
  
  
  categorical_cols <- reactive({
    df <- filtered_df()
    names(df)[vapply(df, function(x) is.factor(x) || is.character(x), logical(1))]
  })
  
  numeric_cols <- reactive({
    df <- filtered_df()
    names(df)[vapply(df, is.numeric, logical(1))]
  })
  
  # 4) Plot builders
  make_barplot <- function(df, col, max_levels = 20) {
    # collapse rare levels if there are too many
    vec <- df[[col]]
    vec <- if (is.factor(vec)) vec else factor(vec)
    # Drop NA for plotting
    tab <- as.data.frame(table(vec, useNA = "no"), stringsAsFactors = FALSE)
    names(tab) <- c("level", "n")
    tab <- tab[order(tab$n, decreasing = TRUE), ]
    
    if (nrow(tab) > max_levels) {
      top <- tab[1:(max_levels - 1), ]
      other_n <- sum(tab$n[(max_levels):nrow(tab)])
      tab <- rbind(top, data.frame(level = "Other", n = other_n))
    }
    mypal <- paletteer_d("ggsci::springfield_simpsons")
    
   # n_colors <- nrow(tab)
    #pal <- mypal[n_colors]
    
   #ggplot(tab, aes(x = reorder(level, n), y = n, fill = level)) +
    ggplot(tab, aes(x = reorder(level, n), y = n)) +
      geom_col() +
      #scale_fill_manual(values = mypal) +
      coord_flip() +
      labs(x = col, y = "Count") +
      theme_classic(base_size = 12) +
      theme(axis.title.y = element_text(face = "bold"))
  }
  
  make_boxplot <- function(df, col) {
    ggplot(df, aes(x = "", y = .data[[col]])) +
      geom_boxplot(outlier.alpha = 0.5) +
      labs(x = col, y = NULL) +
      theme_classic(base_size = 12) +
      theme(axis.title.x = element_text(face = "bold"),
            axis.text.x  = element_blank(),
            axis.ticks.x = element_blank())
  }
  
  # 5) Dynamic UI containers (flex layout grid)
  output$cat_plots_ui <- renderUI({
    cols <- categorical_cols()
    if (length(cols) == 0) return(tags$p("No categorical columns to display."))
    
    div(
      style = "display:flex; flex-wrap:wrap; gap:16px;",
      lapply(cols, function(col) {
        div(
          style = "flex: 1 1 45%; min-width: 380px;",
          if (isTRUE(input$use_plotly_profiles)) {
            plotlyOutput(paste0("cat_", sanitize_id(col)), height = "350px")
          } else {
            plotOutput(paste0("cat_", sanitize_id(col)), height = "350px")
          }
        )
      })
    )
  })
  
  output$num_plots_ui <- renderUI({
    cols <- numeric_cols()
    if (length(cols) == 0) return(tags$p("No numeric columns to display."))
    
    div(
      style = "display:flex; flex-wrap:wrap; gap:16px;",
      lapply(cols, function(col) {
        div(
          style = "flex: 1 1 30%; min-width: 320px;",
          if (isTRUE(input$use_plotly_profiles)) {
            plotlyOutput(paste0("num_", sanitize_id(col)), height = "300px")
          } else {
            plotOutput(paste0("num_", sanitize_id(col)), height = "300px")
          }
        )
      })
    )
  })
  
  # 6) Render all plots reactively
  observe({
    df <- filtered_df()
    req(nrow(df) > 0)
    
    # Categorical
    for (col in categorical_cols()) {
      local({
        ccol <- col
        pid  <- paste0("cat_", sanitize_id(ccol))
        if (isTRUE(input$use_plotly_profiles)) {
          output[[pid]] <- renderPlotly({
            p <- make_barplot(df, ccol, max_levels = input$max_levels)
            ggplotly(p, tooltip = c("x", "y")) %>%
              layout(margin = list(l = 60, r = 20, t = 30, b = 40))
          })
        } else {
          output[[pid]] <- renderPlot({
            make_barplot(df, ccol, max_levels = input$max_levels)
          })
        }
      })
    }
    
    # Numeric
    for (col in numeric_cols()) {
      local({
        ncol_ <- col
        pid   <- paste0("num_", sanitize_id(ncol_))
        if (isTRUE(input$use_plotly_profiles)) {
          output[[pid]] <-  safe_renderPlotly({
            p <- make_boxplot(df, ncol_)
            plotly::ggplotly(p, tooltip = c("y")) %>%
              plotly::layout(margin = list(l = 60, r = 20, t = 30, b = 40))
          })
        } else {
          output[[pid]] <- safe_renderPlot({
            make_boxplot(df, ncol_)
          })
        }
      })
    }
  })
  

  #########MAP#######
  
  output$map_plot <- renderPlotly({

    lat_long_df <- Locations
    
    # convert to sf object
    my_sf <- Locations %>%
      sf::st_as_sf(coords = c("Longitude", "Latitude")) %>%
      sf::st_set_crs(4326)
    
    # World boundaries
    world <- rnaturalearth::ne_countries(scale = "medium", returnclass = "sf")
    
    # ---- Create ggplot ----
    p <- ggplot() +
      {
        if (input$show_world)
          geom_sf(data = world, fill = "lightyellow", color = "gray80")
      } +
      geom_sf(
        data = my_sf,
        aes(text = paste0(
          "<b>Study:</b> ", Study, "<br>",
          "<b>Location:</b> ", Location, "<br>",
          "<b>Latitude:</b> ", round(st_coordinates(geometry)[,2], 3), "<br>",
          "<b>Longitude:</b> ", round(st_coordinates(geometry)[,1], 3)
        )),
        size = input$point_size,
        color = "dodgerblue"
      ) +
      ggtitle("Map of Public HHEAR Studies") +
      theme_classic() +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
        axis.text = element_blank(),
        axis.ticks = element_blank()
      )
    
    # ---- Convert to plotly ----
    plotly::ggplotly(p, tooltip = "text") %>%
      layout(
        dragmode = "pan",
        hoverlabel = list(bgcolor = "white", font = list(size = 12)),
        margin = list(l = 0, r = 0, t = 60, b = 0)
      )    
   
  })
  
  
  
  
  output$debug_vals <- renderPrint({
    list(
      populations_in_pts = table(Location_Plotdata$Population, useNA = "ifany"),
      populations_in_int = table(Location_interval_data$Population, useNA = "ifany"),
      n_analytes_child = length(unique(Location_Plotdata$curr_analyte[Location_Plotdata$Population=="Child"])),
      n_analytes_adult = length(unique(Location_Plotdata$curr_analyte[Location_Plotdata$Population=="Adult"]))
    )
  })
  
  output$analyte_picker <- renderUI({
    validate(
      need(exists("Location_Plotdata"), "Missing object: Location_Plotdata"),
      need(exists("Location_interval_data"), "Missing object: Location_interval_data"),
      need("Population"   %in% names(Location_Plotdata), "Location_Plotdata missing Population"),
      need("curr_analyte" %in% names(Location_Plotdata), "Location_Plotdata missing curr_analyte"),
      need("Population"   %in% names(Location_interval_data), "Location_interval_data missing Population"),
      need("curr_analyte" %in% names(Location_interval_data), "Location_interval_data missing curr_analyte")
    )
    
  
    
    pop <- input$population %||% "Child"
    
    analytes <- sort(unique(c(
      Location_Plotdata        |> filter(Population == pop, !is.na(curr_analyte))        |> pull(curr_analyte),
      Location_interval_data   |> filter(Population == pop, !is.na(curr_analyte))        |> pull(curr_analyte)
    )))
    
    # fallback: if none found for chosen population, offer all analytes (prevents empty UI)
    if (!length(analytes)) {
      analytes <- sort(unique(c(
        Location_Plotdata        |> filter(!is.na(curr_analyte))      |> pull(curr_analyte),
        Location_interval_data   |> filter(!is.na(curr_analyte))      |> pull(curr_analyte)
      )))
    }
    
    display_names <- Chem_Name_Mapping |>
      filter(
        Code %in% analytes,
        !is.na(Full_Name),
        Full_Name != ""
      ) |>
      distinct(Code, Full_Name)
    
    choices_named <- setNames(
      display_names$Code,
      display_names$Full_Name
    )
    
    # missing_codes <- setdiff(analytes, display_names$Code)
    # if (length(missing_codes)) {
    #   choices_named <- c(
    #     choices_named,
    #     setNames(missing_codes, missing_codes)
    #   )
    # }
    # 
    selectizeInput(
      "analyte", "Choose Chemical Exposure:",
      choices  = choices_named,
      selected = preferred_analyte_code(analytes),
      options  = list(
        placeholder = "Start typing an analyte…"
      )
    )
    # selectizeInput(
    #   "analyte", "Choose Chemical Exposure:",
    #   choices  = analytes,
    #   selected = if (length(analytes)) analytes[1] else NULL,
    #   options  = list(placeholder = if (length(analytes)) "Start typing an analyte…" else "No analytes found for this population")
    # )
  })
  
  
  # Helper for `%||%`
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  .sdoh_color_scale <- reactive({
    if (exists("mypal") && is.atomic(mypal)) {
      ggplot2::scale_color_manual(values = mypal)
    } else {
      ggplot2::scale_color_discrete()
    }
  })
  # ---- Filtered data (one analyte at a time) ----
  pts_data <- reactive({
    req(input$population, input$analyte)
    Location_Plotdata %>%
      dplyr::filter(Population == input$population) %>% 
      dplyr::filter(curr_analyte == input$analyte)
  })
  
  
  int_data <- reactive({
    Location_interval_data %>%
      dplyr::filter(Population == input$population) %>% 
      dplyr::filter(curr_analyte == input$analyte)
    
  })

  lod_data <- reactive({
    if (!exists("Location_LOD_Data") || !nrow(Location_LOD_Data)) {
      return(tibble())
    }
    Location_LOD_Data %>%
      dplyr::filter(Population == input$population) %>%
      dplyr::filter(curr_analyte == input$analyte)
  })
  
  name_data <- reactive({
    chem_namedf <- Chem_Name_Mapping %>%
      #dplyr::filter(Population == input$population) %>% 
      dplyr::filter(Code == input$analyte)
    
    chem_name <- unique(chem_namedf$Full_Name)
    chem_name <- chem_name[1]
  })
  
  # ---- Build ggplot ----
  plot_obj <- reactive({
    req(nrow(pts_data()) > 0, nrow(int_data()) > 0)
    
    # X transform function
    xf <- if (isTRUE(input$use_log1p)) function(x) log1p(x) else identity
    
    # Safely supply a color scale (manual if location_colors exists, else discrete)
    color_layer <-
      if (exists("location_colors") && is.vector(location_colors)) {
        scale_color_manual(values = location_colors)
      } else {
        scale_color_discrete()
      }
    
    p <- ggplot() +
      geom_jitter(
        data = pts_data(),
        aes(x = xf(.data$plot_value), y = .data$Location, color = .data$Location, shape = .data$Tissue,
            text = paste0(.data$Study, " | ", .data$Location, "<br>", signif(.data$plot_value, 3), " ", .data$Unit,
                          " (", .data$Tissue, ")")),
        size = 0.6, alpha = 0.9, width = 0, height = 0.12
      ) +
      color_layer + measure_shapes +
      # 25–75% interval
      geom_errorbarh(
        data = int_data(),
        aes(y = .data$n_unique_locations,
            xmin = xf(.data$Tau25),
            xmax = xf(.data$Tau75)),
        height = int_data()$level_height, color = "#1A5354", linewidth = 1
      ) +
      # 5–95% interval
      geom_errorbarh(
        data = int_data(),
        aes(y = .data$n_unique_locations,
            xmin = xf(.data$Tau05),
            xmax = xf(.data$Tau95)),
        height = int_data()$level_height / 2, color = "#1A5354", linewidth = 1, alpha = 0.6
      ) +
      # Median point
      geom_point(
        data = int_data(),
        aes(y = .data$n_unique_locations, x = xf(.data$Tau50)),
        color = "#1A5354", size = 3
      ) +
      geom_point(
        data = lod_data(),
        aes(x = xf(.data$LOD_value), y = .data$Location),
        inherit.aes = FALSE,
        shape = 21,
        color = "#111111",
        fill = "white",
        stroke = 0.8,
        size = 2.2,
        position = position_jitter(width = 0, height = 0.12, seed = 123)
      ) +
      #coord_cartesian(xlim = c(0, input$xmax), clip = "on") +
      labs(
        x = conc_axis_title(pts_data()$Unit, input$use_log1p),
        y = "Location", shape = "Matrix",
        #title = input$analyte
        title = name_data()
        
      ) +
      theme_classic(base_size = 12) +
      theme(
        # Title block
        plot.title = element_text(size = 20, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 12),
        
        # Facet strips
        strip.background = element_rect(fill = "#F2F2F2", color = NA),
        strip.text = element_text(size = 14, face = "bold"),
        
        # Gridlines to aid reading across
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        
        # Axes
        axis.title.x = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        
        # Legend
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        #legend.text = element_markdown(size = 10),
        legend.box = "vertical",
        legend.key.height = unit(10, "mm"),
        legend.key.width  = unit(8, "mm"),
        
        # Panel spacing
        panel.spacing = unit(10, "pt")
      ) +
      guides(
        color = guide_legend(
          nrow = 4,
          byrow = TRUE,
          label.hjust = 0,
          #label.theme = ggtext::element_markdown(size = 14, lineheight = 1.15),
          override.aes = list(size = 3.2, alpha = 1)
        ),
        size = "none"
      ) 
      
      #coord_cartesian(clip = "off")
    
    # Note: we're showing a single analyte (no facet) for clarity.
    # If you ever want to keep a facet (nrow = 1), use:
    # p + facet_wrap(~ curr_analyte, scales = if (input$free_y) "free_y" else "fixed", nrow = 1)
    
    p
  })
  
  # ---- Render interactive plot ----
  output$loc_units <- renderUI(units_table(pts_data()))

  output$child_plot <- renderPlotly({
    ggplotly(plot_obj(), tooltip = "text") %>%
      layout(
        legend = list(orientation = "v", x = 1.02, y = 1, xanchor = "left"),
        margin = list(l = 60, r = 180, t = 40, b = 60)
      )
  })  
  
  output$qc_analyte_picker <- renderUI({
    pop <- input$qc_population %||% "Child"
    analytes <- sort(unique(c(
      Location_Plotdata |> filter(Population == pop, !is.na(curr_analyte)) |> pull(curr_analyte),
      Location_interval_data |> filter(Population == pop, !is.na(curr_analyte)) |> pull(curr_analyte)
    )))
    
    display_names <- Chem_Name_Mapping |>
      filter(Code %in% analytes, !is.na(Full_Name), Full_Name != "") |>
      distinct(Code, Full_Name)
    
    choices_named <- setNames(display_names$Code, display_names$Full_Name)
    
    selectizeInput(
      "qc_analyte", "Choose Chemical Exposure:",
      choices = choices_named,
      selected = preferred_analyte_code(analytes),
      options = list(placeholder = "Start typing an analyte…")
    )
  })
  
  qc_name_data <- reactive({
    chem_namedf <- Chem_Name_Mapping %>%
      dplyr::filter(Code == input$qc_analyte)
    
    chem_name <- unique(chem_namedf$Full_Name)
    chem_name[1]
  })
  
  study_qc_points <- reactive({
    req(input$qc_population, input$qc_analyte)
    Location_Plotdata %>%
      filter(Population == input$qc_population, curr_analyte == input$qc_analyte) %>%
      filter(!is.na(plot_value), !is.na(Study))
  })
  
  study_qc_intervals <- reactive({
    pts <- study_qc_points()
    req(nrow(pts) > 0)
    
    pts %>%
      group_by(Study) %>%
      summarise(
        Tau05 = quantile(plot_value, 0.05, na.rm = TRUE, names = FALSE),
        Tau25 = quantile(plot_value, 0.25, na.rm = TRUE, names = FALSE),
        Tau50 = quantile(plot_value, 0.50, na.rm = TRUE, names = FALSE),
        Tau75 = quantile(plot_value, 0.75, na.rm = TRUE, names = FALSE),
        Tau95 = quantile(plot_value, 0.95, na.rm = TRUE, names = FALSE),
        n_obs = n(),
        .groups = "drop"
      ) %>%
      arrange(Tau50) %>%
      mutate(Study = factor(Study, levels = Study))
  })
  
  study_qc_lod <- reactive({
    req(input$qc_population, input$qc_analyte)
    lod_df <- get0("Location_LOD_Data", ifnotfound = tibble())
    if (!nrow(lod_df)) {
      return(tibble())
    }
    int_df <- study_qc_intervals()
    lod_df %>%
      filter(Population == input$qc_population, curr_analyte == input$qc_analyte) %>%
      group_by(Study) %>%
      summarise(LOD_value = median(LOD_value, na.rm = TRUE), .groups = "drop") %>%
      mutate(Study = factor(Study, levels = levels(int_df$Study)))
  })
  
  output$qc_units <- renderUI(units_table(dplyr::select(study_qc_points(), Study, Tissue, Unit)))

  output$study_qc_plot <- renderPlotly({
    pts <- study_qc_points()
    ints <- study_qc_intervals()
    lods <- study_qc_lod()
    req(nrow(pts) > 0, nrow(ints) > 0)
    
    xf <- if (isTRUE(input$qc_use_log1p)) function(x) log1p(x) else identity
    pts <- pts %>% mutate(Study = factor(Study, levels = levels(ints$Study)))
    
    p <- ggplot() +
      geom_jitter(
        data = pts,
        aes(x = xf(plot_value), y = Study, color = Study, shape = Tissue,
            text = paste0(Study, "<br>", signif(plot_value, 3), " ", Unit, " (", Tissue, ")")),
        size = 0.6, alpha = 0.55, width = 0, height = 0.12
      ) + measure_shapes +
      geom_errorbarh(
        data = ints,
        aes(y = Study, xmin = xf(Tau25), xmax = xf(Tau75)),
        height = 0.28, color = "#1A5354", linewidth = 1
      ) +
      geom_errorbarh(
        data = ints,
        aes(y = Study, xmin = xf(Tau05), xmax = xf(Tau95)),
        height = 0.14, color = "#1A5354", linewidth = 1, alpha = 0.6
      ) +
      geom_point(
        data = ints,
        aes(y = Study, x = xf(Tau50)),
        color = "#1A5354", size = 3
      ) +
      geom_point(
        data = lods,
        aes(x = xf(LOD_value), y = Study),
        inherit.aes = FALSE,
        shape = 21,
        color = "#111111",
        fill = "white",
        stroke = 0.8,
        size = 2.2
      ) +
      labs(
        x = conc_axis_title(pts$Unit, input$qc_use_log1p),
        y = "Study", shape = "Matrix",
        title = paste0(qc_name_data(), " by Study")
      ) +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 18, face = "bold", hjust = 0),
        axis.title.x = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10, face = "bold"),
        legend.position = "right",
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_blank()
      ) +
      guides(color = guide_legend(ncol = 1), size = "none")
    
    ggplotly(p, tooltip = "text") %>%
      layout(
        legend = list(orientation = "v", x = 1.02, y = 1, xanchor = "left"),
        margin = list(l = 80, r = 180, t = 40, b = 60)
      )
  })
  
#############SDOH Plots############
  dataset_index <- reactive({
    df_names <- names(DF_plots)
    tibble(df_name = df_names) %>%
      mutate(
        name_parts = strsplit(df_name, "__", fixed = TRUE),
        study = vapply(name_parts, function(x) x[1], character(1)),
        demographic = vapply(name_parts, function(x) if (length(x) >= 2) x[2] else "Unspecified", character(1))
      ) %>%
      select(-name_parts)
  })
  
  output$study_picker <- renderUI({
    idx <- dataset_index()
    studies <- sort(unique(idx$study))
    req(length(studies) > 0)
    selectInput("study_choice", "Choose study:", choices = studies, selected = studies[1])
  })
  
  DF_plot_re <- reactive({
    req(input$study_choice, input$group_col)
    idx <- dataset_index() %>%
      filter(study == input$study_choice, demographic == input$group_col)
    
    req(nrow(idx) > 0)
    dfs <- lapply(idx$df_name, function(nm) DF_plots[[nm]])
    
    if (length(dfs) == 1) {
      return(dfs[[1]])
    }
    
    dplyr::bind_rows(dfs)
  })
  
  # --- Dynamically generate grouping picker ---
  output$group_picker <- renderUI({
    req(input$study_choice)
    idx <- dataset_index() %>%
      filter(study == input$study_choice)
    
    demographics <- sort(unique(idx$demographic))
    labels <- gsub("_", " ", demographics)
    choices_named <- stats::setNames(demographics, labels)
    
    selectInput("group_col", "Choose demographic:", choices = choices_named, selected = demographics[1])
  })
  
  # --- Palette logic ---
  palette_for_groups <- reactive({
    df <- DF_plot_re()
    req(df, input$group_col)
    groups <- unique(df[[input$group_col]])
    n_groups <- length(groups)
    
    if (exists("mypal", inherits = TRUE) && length(mypal) > 0) {
      base_pal <- unname(mypal)
      pal <- if (length(base_pal) >= n_groups) {
        base_pal[seq_len(n_groups)]
      } else {
        grDevices::colorRampPalette(base_pal)(n_groups)
      }
      names(pal) <- groups
      return(pal)
    }
    
    # backup palette only if mypal is unavailable
    pal <- grDevices::hcl.colors(n_groups, palette = "Set 3")
    names(pal) <- groups
    pal
  })
  
  output$outcome_picker <- renderUI({
    df <- DF_plot_re()
    req(df, "Outcome" %in% names(df))
    choices <- sort(unique(df$Outcome))
    
    display_names <- Chem_Name_Mapping |>
      filter(
        Code %in% choices,
        !is.na(Full_Name),
        Full_Name != ""
      ) |>
      distinct(Code, Full_Name)
    
    choices_named <- setNames(
      display_names$Code,
      display_names$Full_Name
    )
    
    # missing_codes <- setdiff(analytes, display_names$Code)
    # if (length(missing_codes)) {
    #   choices_named <- c(
    #     choices_named,
    #     setNames(missing_codes, missing_codes)
    #   )
    # }
    # 
    selectizeInput(
      "outcome", "Choose Chemical Exposure:",
      choices  = choices_named,
      selected = preferred_analyte_code(choices),
      options  = list(
        placeholder = "Start typing an analyte…"
      )
    )
    
    #selectInput("outcome", "Exposure:", choices = choices, selected = choices[1])
  })
  
  
  # --- Plot logic ---
  output$quantile_plot <- renderPlotly({
    DF_plot <- DF_plot_re()
   # print(DF_plot)
    req(DF_plot, input$group_col, input$outcome, "Outcome" %in% names(DF_plot), "Value" %in% names(DF_plot))
    group_col <- input$group_col
    group_palette <- palette_for_groups()
    
    n_groups   <- length(unique(DF_plot[[group_col]]))
    print(n_groups)
    error_y    <- (n_groups + 1) / 2
    bar_height <- n_groups + 1
    
    DF_one <- dplyr::filter(DF_plot, Outcome == input$outcome)
    

      chem_namedf <- Chem_Name_Mapping %>%
        #dplyr::filter(Population == input$population) %>% 
        dplyr::filter(Code == input$outcome)
      
      chem_name <- unique(chem_namedf$Full_Name)
      chem_name <- chem_name[1]
  
    
    Intervals <- DF_one %>% 
      dplyr::select(Outcome,Tau05,Tau25,Tau50,Tau75,Tau95)
    
    Intervals_y <- Intervals %>%
      mutate(error_y = error_y, bar_height = bar_height) %>%
      filter(Outcome %in% DF_plot$Outcome)
    
    validate(
      need(nrow(DF_plot) > 0, "No rows in the selected dataset."),
      need(nrow(Intervals_y) > 0, "No matching outcomes found in Intervals.")
    )
    #print(mypal)

    p <- ggplot() +
      geom_jitter(
        data = DF_one,
        aes(x = log1p(Value), y = !!sym(group_col), color = !!sym(group_col)), size = 1.5, alpha = 0.5
      ) +
      scale_color_manual(values = group_palette) +
      # 25–75
      geom_errorbarh(
        data = Intervals_y,
        aes(y = error_y, xmin = log1p(Tau25), xmax = log1p(Tau75)),
        inherit.aes = FALSE, height = bar_height, color = "#1A5354", linewidth = 1
      ) +
      # 5–95
      geom_errorbarh(
        data = Intervals_y,
        aes(y = error_y, xmin = log1p(Tau05), xmax = log1p(Tau95)),
        inherit.aes = FALSE, height = bar_height / 2, color = "#1A5354", linewidth = 1, alpha = 0.6
      ) +
      # median
      geom_point(
        data = Intervals_y,
        aes(y = error_y, x = log1p(Tau50)),
        inherit.aes = FALSE, color = "#1A5354", size = 6
      ) +
      facet_wrap(~ Outcome, scales = "free_x",nrow = 1) +
      theme_classic() +
      theme_classic(base_size = 12) +
      theme(
        plot.title = element_text(size = 20, face = "bold", hjust = 0),
        plot.subtitle = element_text(size = 12),
        strip.background = element_rect(fill = "#F2F2F2", color = NA),
        strip.text = element_text(size = 14, face = "bold"),
        panel.grid.major.x = element_line(color = "grey90", linewidth = 0.4),
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size = 13, face = "bold"),
        axis.text.x = element_text(size = 8, face = "bold"),
        axis.text.y = element_blank(),
        axis.ticks = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        legend.box = "vertical",
        legend.key.height = unit(5, "mm"),
        legend.key.width  = unit(8, "mm"),
        panel.spacing = unit(10, "pt")
      ) +
      guides(
        color = guide_legend(
          nrow = 1,
          byrow = TRUE,
          override.aes = list(size = 3.2, alpha = 1)
        )
      ) +
      labs(title = chem_name, y = gsub("_", " ", group_col), color = gsub("_", " ", group_col)) +
      coord_cartesian(clip = "off")  +
      labs(x = conc_axis_title(
        Matrix_Lookup$Unit[Matrix_Lookup$Study == input$study_choice & Matrix_Lookup$Analyte_Code == input$outcome],
        log1p = TRUE))
    
    ggplotly(p) %>% layout(autosize = "FALSE") 
  })
  
  ######match plots#####
  # ---- Strata ridgeline plot (filtered by curr_analyte) ----
  output$ridge_analyte_picker <- renderUI({
    validate(
      need(exists("plot_data"), "Missing object: plot_data"),
      need("curr_analyte" %in% names(plot_data), "plot_data missing column: curr_analyte")
    )
    
    analytes <- sort(unique(plot_data$curr_analyte[!is.na(plot_data$curr_analyte)]))
    
    display_names <- Chem_Name_Mapping |>
      filter(
        Code %in% analytes,
        !is.na(Full_Name),
        Full_Name != ""
      ) |>
      distinct(Code, Full_Name)
    
    choices_named <- setNames(
      display_names$Code,
      display_names$Full_Name
    )
    
    # missing_codes <- setdiff(analytes, display_names$Code)
    # if (length(missing_codes)) {
    #   choices_named <- c(
    #     choices_named,
    #     setNames(missing_codes, missing_codes)
    #   )
    # }
    # 
    selectizeInput(
      "ridge_analyte", "Choose Chemical Exposure:",
      choices  = choices_named,
      selected = preferred_analyte_code(analytes),
      options  = list(
        placeholder = "Start typing an analyte…"
      )
    )
    
    
    # selectizeInput(
    #   "ridge_analyte", "Choose Chemical Exposure:",
    #   choices = analytes,
    #   selected = if (length(analytes)) analytes[1] else NULL,
    #   options = list(placeholder = "Start typing an analyte…")
    # )
  })
  
  ridge_df <- reactive({
    validate(
      need(exists("plot_data"), "Missing object: plot_data"),
      need(all(c("plot_value","Strata","Location","facet_label","curr_analyte") %in% names(plot_data)),
           "plot_data must contain: plot_value, Strata, Location, facet_label, curr_analyte")
    )
    req(input$ridge_analyte)
    
    plot_data %>%
      dplyr::filter(curr_analyte == input$ridge_analyte) %>%   # <-- the only requested change
      dplyr::filter(!is.na(plot_value), !is.na(Strata), !is.na(Location))
  })
  
  strata_chem_name <- reactive({
    
    chem_namedf <- Chem_Name_Mapping %>%
      #dplyr::filter(Population == input$population) %>% 
      dplyr::filter(Code == input$ridge_analyte)
    
    chem_name <- unique(chem_namedf$Full_Name)
    chem_name <- chem_name[1]
    
  })
  

  # output$conc_ridge_plot <- renderPlotly({
  #   df <- ridge_df()
  #   req(nrow(df) > 0)
  #   
  #   # Use your existing location_colors if available; otherwise fall back
  #   fill_scale <-
  #     if (exists("location_colors") && is.vector(location_colors)) {
  #       scale_fill_manual(values = location_colors,
  #                         guide = guide_legend(title = NULL, byrow = TRUE))
  #     } else {
  #       scale_fill_discrete(guide = guide_legend(title = NULL, byrow = TRUE))
  #     }
  #   
  #   p <- ggplot(
  #     df,
  #     aes(
  #       x = log1p(plot_value),
  #       y = forcats::fct_rev(Strata),
  #       fill = factor(Location)
  #     )
  #   ) +
  #     ggridges::geom_density_ridges(
  #       scale = 1.6,
  #       rel_min_height = 0.01,
  #       size = 0.2,
  #       color = "white",
  #       alpha = 0.9
  #     ) +
  #     fill_scale +
  #     labs(
  #       x = "Concentration",
  #       y = "Age–Sex Stratum",
  #       title = "Exposure Concentrations Across Matched Strata"
  #     ) +
  #     facet_wrap(~ facet_label, scales = "free", ncol = 2, strip.position = "top") +
  #     theme_minimal(base_size = 12) +
  #     theme(
  #       panel.grid = element_blank(),
  #       strip.background = element_rect(fill = "grey95", color = NA),
  #       strip.text = element_text(size = 11, face = "bold"),
  #       axis.text.y = element_text(size = 10, color = "grey20", face = "bold"),
  #       axis.text.x = element_text(size = 9),
  #       plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
  #       legend.position = "bottom",
  #       panel.spacing = unit(1, "lines")
  #     )
  #   
  #   ggplotly(p) %>% layout(legend = list(orientation = "h"))
  # })
  # 
  
  output$ridge_units <- renderUI(units_table(ridge_df()))

  output$conc_ridge_plot <- renderPlot({
    df <- ridge_df()
    req(nrow(df) > 0)
    
    p <- ggplot(df, aes(x = log1p(plot_value),
                        y = forcats::fct_rev(Strata),
                        fill = factor(Location),
                        linetype = Tissue)) +
      ggridges::geom_density_ridges(
        scale = 1.6,
        rel_min_height = 0.01,
        size = 0.2,
        color = "white",
        alpha = 0.9
      ) +
      scale_fill_manual(values = location_colors, guide = guide_legend(title = NULL, byrow = TRUE, ncol = 4)) +
      scale_linetype_manual(name = "Matrix", values = c("dashed", "solid")) +
      #facet_wrap(~ facet_label, scales = "free", ncol = 2, strip.position = "top") +
      facet_wrap(~ facet_label,scales = "free_y", ncol = 2, strip.position = "top") +
      theme_minimal(base_size = 12) +
      theme(
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "grey95", color = NA),
        strip.text = element_text(size = 11, face = "bold"),
        axis.text.y = element_text(size = 10, color = "grey20", face = "bold"),
        axis.text.x = element_text(size = 9),
        plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
        legend.position = "bottom",
        legend.box = "vertical",
        legend.text = element_text(size = 9),
        panel.spacing = unit(1, "lines")
      ) +
      labs(x = conc_axis_title(df$Unit, log1p = TRUE), y = "Age–Sex Stratum",
           title = paste0(strata_chem_name(),""),
)
    
    p
  }, res = 110)
  
  
  #Data Download
  
  # 1) Put your downloadable objects in one place
  download_registry <- reactive({
    list(
      #Location_Plotdata        = Location_Plotdata,
      Location_interval_data   = Location_interval_data %>% 
        dplyr::select(-bar_height,-Flag,-error_y,-n_unique_locations,-level_height) %>%
        # Pooled intervals span cohorts; list every unit that fed them (2026-11-19).
        dplyr::left_join(interval_units, by = c("Population", "curr_analyte")),
      Matched_Data     = plot_data
    )
  })
  
  # 2) Optional filtering logic based on current inputs
  download_data <- reactive({
    reg <- download_registry()
    key <- req(input$download_dataset)
    
    df <- reg[[key]]
    validate(need(!is.null(df), "Selected dataset not available."))
    
    if (!isTRUE(input$download_filtered)) return(df)

    
    df
  })
  
  # 3) Preview table
  output$download_preview <- renderTable({
    df <- download_data()
    utils::head(df, 20)
  })
  
  # 4) CSV download handler
  output$download_csv <- downloadHandler(
    filename = function() {
      paste0(
        input$download_dataset,
        if (isTRUE(input$download_filtered)) "_filtered" else "_full",
        "_",
        format(Sys.Date(), "%Y-%m-%d"),
        ".csv"
      )
    },
    content = function(file) {
      df <- download_data()
      
      # Prefer readr if available (faster + consistent), else base utils
      if (requireNamespace("readr", quietly = TRUE)) {
        readr::write_csv(df, file)
      } else {
        utils::write.csv(df, file, row.names = FALSE)
      }
    }
  )

  
  
  #Dyad #######
  
  # Populate analyte dropdown
  observe({
    req(exists("Dyad_df"))
    
    choices <- Dyad_df |>
      dplyr::distinct(Full_Name) |>
      dplyr::arrange(Full_Name) |>
      dplyr::pull(Full_Name)


    # choices <- Dyad_df |>
    #   dplyr::distinct(Analyte_Code) |>
    #   dplyr::arrange(Analyte_Code) |>
    #   dplyr::pull(Analyte_Code)
    
  

    updateSelectizeInput(
      session,
      "dyad_analyte",
      choices = choices,
      selected = preferred_analyte_name(choices),
      server = TRUE
    )
  })
  
  
  # Reactive filtered data
  dyad_plot_df <- reactive({
    req(input$dyad_analyte, exists("Dyad_df"))
    
    Dyad_df 
    
    Dyad_df |>
      dplyr::filter(
        Full_Name == input$dyad_analyte,
        !is.na(Concentration_Mat),
        !is.na(Concentration_Child)
      )
  })
  
  mypal <- c("#FED439FF", "#709AE1FF")
  # Scatter plot
  output$dyad_scatter <- plotly::renderPlotly({
    dfp <- dyad_plot_df()
    req(nrow(dfp) > 0)
    
    dfp <- dfp %>% 
      dplyr::filter(Concentration_Mat >0)
    eps <- 1e-8
    
    x_vals <- if (isTRUE(input$dyad_log10)) log10(dfp$Concentration_Mat + eps) else dfp$Concentration_Mat
    y_vals <- if (isTRUE(input$dyad_log10)) log10(dfp$Concentration_Child + eps) else dfp$Concentration_Child
    
    plot_dat <- dfp |>
      dplyr::mutate(x = x_vals, y = y_vals)
    
    p <- ggplot(plot_dat, aes(x = x, y = y,color = Study)) +
      geom_point(alpha = 0.6, size = 2) +
      scale_color_manual(values =mypal) +
      labs(
        x = paste0(if (isTRUE(input$dyad_log10)) "Maternal concentration (log10)" else "Maternal concentration", unit_suffix(plot_dat$Unit)),
        y = paste0(if (isTRUE(input$dyad_log10)) "Child concentration (log10)" else "Child concentration", unit_suffix(plot_dat$Unit)),
        title = paste0("Analyte: ", input$dyad_analyte)
      ) +
      theme_classic(base_size = 13)
    
    if (isTRUE(input$dyad_show_1to1)) {
      rng <- range(c(plot_dat$x, plot_dat$y), finite = TRUE)
      p <- p +
        geom_abline(intercept = 0, slope = 1, linetype = "dashed") +
        coord_cartesian(xlim = rng, ylim = rng)
    }
    
    plotly::ggplotly(p, tooltip = c("x", "y")) %>%
      plotly::layout(margin = list(l = 60, r = 20, t = 50, b = 60))
  })
  
  
  
  output$dyad_units <- renderUI(units_table(dplyr::select(dyad_plot_df(), Study, Unit)))

  # N count text
  output$dyad_note <- renderText({
    dfp <- dyad_plot_df()
    paste0("N dyads: ", format(nrow(dfp), big.mark = ","))
  })
 
  # --- SERVER: build a per-analyte log-ratio violin plot under the scatter ---
  # Assumes Dyad_df has: Analyte_Code, Concentration_Mat, Concentration_Child
  # and you have Chem_Name_Mapping with Code -> Category (or swap Category source)
  
  dyad_ratio_df <- reactive({
    req(exists("Dyad_df"))
    
    Dyad_df 
    # |>
    #   dplyr::filter(
    #     Full_Name == input$dyad_analyte,
    #     !is.na(Concentration_Mat),
    #     !is.na(Concentration_Child)
    #   )
  })
  
  output$dyad_violin <- plotly::renderPlotly({
    dfp <- dyad_ratio_df()
    req(nrow(dfp) > 0)
    
    # palette (keeps your style; uses a discrete palette without hardcoding)
    cats <- sort(unique(dfp$Category))
    mypal <- c("#FED439FF", "#709AE1FF")
    
    
    vi_plot <- ggplot2::ggplot(
      dfp,
      #ggplot2::aes(x = log_ratio, y = Full_Name, color = Study, fill = Study)
      ggplot2::aes(x = log_ratio, y = Analyte_Code)
      
    ) +
      geom_violin(
        linewidth = 0.9,
        alpha = 1
      ) +
      ggplot2::geom_boxplot(
        width = 0.4,
        outlier.size = 0.5,
        outlier.alpha = 0.6,
        fatten = 0.7
      ) +
      ggplot2::geom_vline(
        xintercept = 0,
        linetype = "dashed",
        color = "grey40",
        linewidth = 0.5
      ) +
      ggplot2::labs(
        x = "log(Child / Mother)",
        y = NULL,
        title = "Child–Mother Log Ratios Across Exposures"
      ) +
      ggplot2::scale_fill_manual(values = mypal) +
      ggplot2::scale_color_manual(values = mypal) +
      ggplot2::scale_x_continuous(
        limits = c(-5, 5),
        breaks = c(-5, -3, -1, 0, 1, 3, 5),
        labels = c(
          "148× Mother",
          "20× Mother",
          "2.7× Mother",
          "Equal",
          "2.7× Child",
          "20× Child",
          "148× Child"
        )
      ) +
      ggplot2::theme_classic(base_size = 16) +
      ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", hjust = 0.5, size = 14),
        axis.text.y = ggplot2::element_text(size = 13),
        axis.ticks.y = ggplot2::element_blank(),
        legend.position = "bottom",
        axis.text.x = ggplot2::element_text(angle = -45)
      ) +
      ggplot2::coord_cartesian(clip = "off") +
      ggplot2::coord_flip()
    
    plotly::ggplotly(vi_plot, tooltip = c("x", "y", "colour")) %>%
      plotly::layout(margin = list(l = 80, r = 20, t = 60, b = 80))
  })
  
  
}

shinyApp(ui, server)
