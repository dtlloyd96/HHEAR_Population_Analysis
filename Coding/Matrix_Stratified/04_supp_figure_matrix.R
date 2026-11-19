## ============================================================================
## 04_supp_figure_matrix.R  --  supplementary figure for the matrix confound
##
## Two panels, in the order the argument runs:
##   A. the MECHANISM  - the same analyte differs by orders of magnitude between
##      blood and urine, so medium is a huge latent effect.
##   B. the CONSEQUENCE - location contribution pooled vs within matrix. Every
##      single-matrix class lands exactly on its pooled value, which is the
##      built-in control: stratification moves only what it should.
##
## Form: paired dot plot (dumbbell) both panels -- the data's job is comparing
## two values within a group, which is what a dumbbell is for.
## Colour: two categorical hues only (Okabe-Ito blue/orange, CVD-validated,
## worst adjacent dE 29.2 protan). "Published pooled" is a BASELINE, not a third
## category, so it is a grey reference mark rather than a hue.
## 2026-11-04: titles, subtitles and number labels removed for readability
## (author request); pooled shown as a legend-keyed grey point.
## Writes only to Output/Matrix_Stratified/.
## ============================================================================
suppressPackageStartupMessages({
  library(tidyverse); library(data.table); library(patchwork)
})

ROOT <- Sys.getenv("HHEAR_ROOT", getwd())   # run from the repository root
OUT  <- file.path(ROOT, "Output", "Matrix_Stratified")

URINE <- "#0072B2"   # blue
BLOOD <- "#E69F00"   # orange
REF   <- "grey45"    # pooled baseline: reference, not a category

base_theme <- theme_classic(base_size = 10) +
  theme(axis.text        = element_text(size = 8.5, colour = "grey15"),
        axis.title       = element_text(size = 9.5),
        plot.title       = element_text(size = 11, face = "bold"),
        plot.subtitle    = element_text(size = 8.5, colour = "grey35"),
        panel.grid.major.y = element_line(colour = "grey93", linewidth = 0.25),
        legend.position  = "bottom",
        legend.title     = element_blank(),
        legend.text      = element_text(size = 8.5))

## ---------------------------------------------------------------- Panel A ---
load(file.path(ROOT, "Input/Harmonized_Datasets/Harmonized_Targeted.RData"))
CW <- fread(file.path(ROOT, "Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv"))
H <- as.data.table(Targeted_Data_All)
H <- H[!is.na(Matrix_Class) & !is.na(Concentration) & Concentration > 0 &
         Matrix_Class %in% c("Blood", "Urine")]
## Restrict to the class the argument is about. Without this the panel also
## picks up inflammation markers (IL6, IL1B, TNFA) and Mg, which are measured in
## both media but are not the trace elements the subtitle claims.
H <- H[Analyte_Code %in% CW[Category == "TraceElements", Code]]

med <- H[, .(med = median(Concentration)), by = .(Analyte_Code, Matrix_Class)]
both <- dcast(med, Analyte_Code ~ Matrix_Class, value.var = "med")[
  !is.na(Blood) & !is.na(Urine)]
both[, fold := Blood / Urine]
setorder(both, -fold)
both[, Analyte_Code := factor(Analyte_Code, levels = rev(Analyte_Code))]

pa_long <- melt(both[, .(Analyte_Code, Blood, Urine)], id.vars = "Analyte_Code",
                variable.name = "Matrix", value.name = "med")

pA <- ggplot(both, aes(y = Analyte_Code)) +
  geom_segment(aes(x = pmin(Blood, Urine), xend = pmax(Blood, Urine),
                   yend = Analyte_Code), colour = "grey75", linewidth = 0.5) +
  # legend comes from panel B, which carries all three keys
  geom_point(data = pa_long, aes(x = med, colour = Matrix), size = 2.4, show.legend = FALSE) +
  scale_x_log10(labels = function(x) trimws(formatC(x, format = "fg", digits = 2)),
                expand = expansion(mult = c(0.06, 0.06))) +
  scale_colour_manual(values = c(Urine = URINE, Blood = BLOOD)) +
  labs(title = "Trace-element concentration by sample matrix",
       x = "Median concentration (log scale)", y = NULL) +
  base_theme

## ---------------------------------------------------------------- Panel B ---
w <- fread(file.path(OUT, "location_contribution_by_matrix.csv"))
w <- w[!is.na(published)]
## exposure classes only (Inflammation is an assay panel, not one of the 11)
w <- w[Category != "Inflammation"]
setorder(w, published)
w[, Category := factor(Category, levels = Category)]
w[, n_matrix := (!is.na(Urine)) + (!is.na(Blood))]

pb_long <- melt(w[, .(Category, Pooled = published, Urine, Blood)], id.vars = "Category",
                variable.name = "Matrix", value.name = "R2")[!is.na(R2)]
pb_long[, Matrix := factor(Matrix, c("Pooled", "Urine", "Blood"))]

pB <- ggplot(w, aes(y = Category)) +
  geom_segment(aes(x = pmin(published, Urine, Blood, na.rm = TRUE),
                   xend = pmax(published, Urine, Blood, na.rm = TRUE),
                   yend = Category), colour = "grey80", linewidth = 0.5) +
  geom_point(data = pb_long, aes(x = R2, colour = Matrix), size = 2.6) +
  scale_x_continuous(limits = c(0, NA), expand = expansion(mult = c(0.02, 0.05))) +
  labs(title = "Location contribution, pooled vs within matrix",
       x = "Location contribution (%)", y = NULL) +
  base_theme

## one shared legend, fixed order, so identity is never colour-alone and the two
## panels cannot disagree about which hue is which
fig <- pA + pB +
  plot_layout(widths = c(1, 1.15), guides = "collect") +
  plot_annotation(tag_levels = "A") &
  scale_colour_manual(values = c(Pooled = REF, Urine = URINE, Blood = BLOOD),
                      breaks = c("Pooled", "Urine", "Blood"),
                      labels = c(Pooled = "Pooled (all matrices)", Urine = "Urine", Blood = "Blood")) &
  theme(legend.position = "bottom")

ggsave(file.path(OUT, "SuppFig_Matrix_Confound.pdf"), fig, width = 11, height = 5.2)
ggsave(file.path(OUT, "SuppFig_Matrix_Confound.png"), fig, width = 11, height = 5.2, dpi = 200)
cat("written: Output/Matrix_Stratified/SuppFig_Matrix_Confound.{pdf,png}\n")
print(both[, .(Analyte_Code, Blood = signif(Blood,3), Urine = signif(Urine,3), fold = round(fold,1))])
