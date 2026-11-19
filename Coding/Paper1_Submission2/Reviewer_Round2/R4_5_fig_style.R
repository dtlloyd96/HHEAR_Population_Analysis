# R4_5_fig_style.R
# Shared look for every Reviewer-4 (R4_5) figure: one theme, one palette, one
# save function, and small helpers for location labels and correlation labels.
# Sourced by R4_5_report_figures.R, R4_5_supp_figures.R and R4_5_map.R.

suppressMessages({library(ggplot2); library(ggrepel); library(patchwork)})

pal <- list(red = "#B23A3A", blue = "#2F5D8C", orange = "#D98A2B", green = "#3C8D5A", purple = "#6A4C93",
            grey = "#8C9196", dark = "#2B2F33", mid = "#5A6068", light = "#E3E5E8")
pal_region  <- c("US" = pal$blue, "Non-US" = pal$red)
pal_verdict <- c("FDR < 0.05" = pal$red, "p < 0.05" = pal$orange, "No association" = pal$grey, "Opposite direction" = pal$blue)
pal_urine   <- c("Dilution-adjusted" = 16, "Unadjusted" = 21)

theme_r4 <- function(base = 10) {
  theme_minimal(base_size = base) +
    theme(panel.grid.minor = element_blank(),
          panel.grid.major = element_line(colour = "#ECEDEF", linewidth = 0.35),
          panel.border = element_rect(fill = NA, colour = "#C9CDD2", linewidth = 0.4),
          plot.title = element_text(face = "bold", size = base + 1, colour = pal$dark, margin = margin(b = 2)),
          plot.subtitle = element_text(size = base - 1, colour = pal$mid, margin = margin(b = 6)),
          plot.title.position = "plot", plot.caption.position = "plot",
          plot.caption = element_text(size = base - 2, colour = pal$mid, hjust = 0),
          axis.title = element_text(size = base - 0.5, colour = "#3A3F45"),
          axis.text = element_text(size = base - 1.5, colour = "#4A4F55"),
          strip.text = element_text(face = "bold", size = base - 1, colour = pal$dark, hjust = 0, margin = margin(b = 3, t = 3)),
          strip.background = element_blank(),
          legend.position = "bottom", legend.title = element_text(size = base - 1, colour = pal$mid),
          legend.text = element_text(size = base - 1.5), legend.key.size = unit(0.4, "cm"),
          legend.margin = margin(t = -2), plot.tag = element_text(face = "bold", size = base + 3),
          plot.margin = margin(6, 10, 6, 6))
}
theme_set(theme_r4())
update_geom_defaults("text_repel", list(size = 2.7, colour = "#2B2F33"))

save_fig <- function(g, dir, name, w, h) {
  dir.create(dir, showWarnings = FALSE, recursive = TRUE)
  ggsave(file.path(dir, paste0(name, ".png")), g, width = w, height = h, dpi = 300, bg = "white")
  ggsave(file.path(dir, paste0(name, ".pdf")), g, width = w, height = h, device = cairo_pdf, bg = "white")
}

lab_loc <- function(x) gsub("_", " ", x)
an_lab <- c(NAP1 = "1-naphthol", NAP2 = "2-naphthol", PYR1 = "1-hydroxypyrene", PHEN1 = "1-OH-phenanthrene", NNAL = "NNAL",
            COTT = "cotinine", As = "arsenic", Hg = "mercury", Pb = "lead", Cd = "cadmium", MBZP = "MBzP", MEHHP = "MEHHP", MEP = "MEP",
            PRPB = "propylparaben", MEPB = "methylparaben", PBA = "3-PBA", TCP = "TCPy", DEP = "DEP", DETP = "DETP", DMP = "DMP",
            DMTP = "DMTP", PNP = "4-nitrophenol", PFOS = "PFOS", PFOA = "PFOA")
alab <- function(x) ifelse(is.na(an_lab[x]), x, an_lab[x])

# Forest-plot layers: reference line at 1, interval, point (log x).
forest_layers <- function(...) list(geom_vline(xintercept = 1, linetype = 2, colour = "#9AA0A6", linewidth = 0.4),
                                    geom_linerange(aes(xmin = lo, xmax = hi), linewidth = 0.6, ...),
                                    geom_point(size = 2, ...))
