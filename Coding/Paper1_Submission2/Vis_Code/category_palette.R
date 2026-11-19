## ============================================================================
## category_palette.R  -- single source of truth for exposure-category colors
##
## Values are the existing IGV-style palette already used by vi_plot; nothing is
## recolored, only made consistent. Keys include BOTH spelling variants
## (e.g. PCB/PCBs, Pesticide/Pesticides, VOC/VOCs, TraceElements/Trace_Elements,
## FlameRetardants/Flame_Retardant) so a category renders identically across
## figures regardless of how its name happens to be spelled in that dataset.
##
## Usage in a figure:
##   source("category_palette.R")
##   p + scale_fill_manual(values = CATEGORY_PALETTE, na.value = "grey80") +
##       scale_colour_manual(values = CATEGORY_PALETTE, na.value = "grey80")
## ============================================================================

CATEGORY_PALETTE <- c(
  "FlameRetardants" = "#FED439FF",  # yellow
  "Flame_Retardant" = "#FED439FF",
  "PAH"             = "#709AE1FF",  # blue
  "Parabens"        = "#8A9197FF",  # gray
  "PCB"             = "#D2AF81FF",  # tan
  "PCBs"            = "#D2AF81FF",
  "Pesticide"       = "#FD7446FF",  # orange
  "Pesticides"      = "#FD7446FF",
  "PFAS"            = "#D5E4A2FF",  # light green
  "Phenols"         = "#197EC0FF",  # deep blue
  "Phthalates"      = "#F05C3BFF",  # red-orange
  "Smoking"         = "#46732EFF",  # dark green
  "TraceElements"   = "#71D0F5FF",  # light blue
  "Trace_Elements"  = "#71D0F5FF",
  "VOC"             = "#370335FF",  # dark purple
  "VOCs"            = "#370335FF"
)
