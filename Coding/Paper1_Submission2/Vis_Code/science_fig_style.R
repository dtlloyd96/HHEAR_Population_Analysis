# science_fig_style.R -- Science figure specifications, applied to existing ggplot /
# patchwork objects (aeh0807 revision checklist, "Guide to Preparing Figures"):
#   * one file per figure, vector: PDF (cairo, fonts embedded as live text) + SVG
#   * width 1 column = 3.6 in (9 cm) or 2 columns = 7.25 in (18.3 cm); figure plus
#     legend must fit the page box, so main-figure height is kept <= ~8 in
#   * sans serif (Helvetica), at most three type sizes: 6 pt (ticks, legend text,
#     in-plot labels), 7 pt (axis titles, strip and panel titles), 10 pt bold panel
#     letters (A, B, C)
#   * minimum line weight 0.28 pt at final size (ggplot linewidth >= 0.1)
#   * raster elements >= 300 ppi at final size; RGB
# Sourced by Science_Figures.R.
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })

SCI_W1 <- 3.6; SCI_W2 <- 7.25
SCI_FONT <- "Helvetica"
SCI_SMALL <- 6; SCI_MID <- 7; SCI_TAG <- 10
SCI_MIN_LW <- 0.1                         # mm; 0.1 * .pt = 0.28 pt

theme_science <- function() theme(
  text            = element_text(family = SCI_FONT, size = SCI_MID, colour = "black"),
  axis.text       = element_text(size = SCI_SMALL, colour = "black"),
  axis.title      = element_text(size = SCI_MID),
  strip.text      = element_text(size = SCI_MID, margin = margin(1.5, 1.5, 1.5, 1.5)),
  plot.title      = element_text(size = SCI_MID, face = "bold", margin = margin(0, 0, 2, 0)),
  plot.subtitle   = element_text(size = SCI_SMALL, margin = margin(0, 0, 2, 0)),
  plot.caption    = element_text(size = SCI_SMALL),
  legend.text     = element_text(size = SCI_SMALL),
  legend.title    = element_text(size = SCI_MID),
  legend.key.size = unit(8, "pt"),
  legend.key.height = unit(7, "pt"), legend.key.width = unit(11, "pt"),   # colour bars are 5 x width
  legend.key.spacing.x = unit(2, "pt"), legend.key.spacing.y = unit(1, "pt"),
  legend.spacing  = unit(2, "pt"),
  legend.margin   = margin(1, 1, 1, 1),
  legend.box.spacing = unit(3, "pt"),
  plot.margin     = margin(3, 3, 3, 3),
  plot.tag        = element_text(family = SCI_FONT, size = SCI_TAG, face = "bold"),
  plot.tag.position = c(0, 1)
)

# Scale the geoms that were sized for a 15-20 in canvas: text to 6 pt, points and
# lines by `pt_scale`, never below the minimum line weight.
sci_fix_layers <- function(p, pt_scale = 0.5, text_pt = SCI_SMALL) {
  if (!inherits(p, "ggplot")) return(p)
  # a saved plot can carry two scales for one aesthetic (a palette added on top of
  # the original); keep the last, which is the one ggplot would use
  sc <- p$scales$scales
  if (length(sc) > 1) {
    keep <- rep(TRUE, length(sc))
    for (i in seq_along(sc)) for (j in seq_along(sc)) if (j > i && length(intersect(sc[[i]]$aesthetics, sc[[j]]$aesthetics))) keep[i] <- FALSE
    p$scales <- p$scales$clone(); p$scales$scales <- sc[keep]
  }
  for (i in seq_along(p$layers)) {
    l <- p$layers[[i]]
    g <- class(l$geom)[1]
    ap <- l$aes_params
    if (grepl("Text|Label|RichText", g)) {
      ap$size <- text_pt / .pt
    } else {
      if (!is.null(ap$size) && is.numeric(ap$size)) ap$size <- max(ap$size * pt_scale, 0.2)
      if (!is.null(ap$linewidth) && is.numeric(ap$linewidth)) ap$linewidth <- max(ap$linewidth * pt_scale, SCI_MIN_LW)
      if (!is.null(ap$stroke) && is.numeric(ap$stroke)) ap$stroke <- max(ap$stroke * pt_scale, 0.2)
    }
    l$aes_params <- ap
    p$layers[[i]] <- l
  }
  # mapped sizes (e.g. point area by number of cohorts) shrink by the same factor
  for (i in seq_along(p$scales$scales)) {
    sc <- p$scales$scales[[i]]
    if (any(c("size", "linewidth") %in% sc$aesthetics) && is.function(sc$palette)) {
      # ggproto resolves its parent lazily by name, so bind it in its own frame
      p$scales$scales[[i]] <- local({
        parent <- sc; pal <- sc$palette; k <- pt_scale
        ggproto(NULL, parent, palette = function(x) pal(x) * k)
      })
    }
  }
  p
}

# Apply to every plot inside a patchwork (recursively).
sci_walk <- function(p, f) {
  if (inherits(p, "patchwork")) {
    p$patches$plots <- lapply(p$patches$plots, sci_walk, f = f)
  }
  f(p)
}

# Theme for one plot. The panels set sizes on leaf elements (axis.text.y = 9,
# strip.text = 14, ...), which beat anything set on a parent, so every text
# element the plot defines is re-sized: titles, axis titles, strips and legend
# titles 7 pt; everything else 6 pt; tags 10 pt. Elements drawn with
# ggtext::element_markdown keep that class.
sci_size_for <- function(nm) {
  if (grepl("^plot\\.tag", nm)) return(SCI_TAG)
  if (grepl("^(plot\\.title|axis\\.title|strip\\.text|legend\\.title)", nm)) return(SCI_MID)
  if (nm == "text" || nm == "title") return(SCI_MID)
  SCI_SMALL
}
sci_theme_for <- function(q) {
  th <- theme_science()
  # anything the plot switched off (maps: axes, titles) stays off
  for (nm in names(q$theme)) if (inherits(q$theme[[nm]], "element_blank") ||
                                 inherits(q$theme[[nm]], "ggplot2::element_blank")) th[[nm]] <- element_blank()
  for (nm in union(names(q$theme), names(th))) {
    old <- q$theme[[nm]]; cur <- th[[nm]]
    if (inherits(old, "element_blank") || inherits(old, "ggplot2::element_blank")) next
    is_txt <- function(x) inherits(x, "element_text") || inherits(x, "ggplot2::element_text")
    if (!(is_txt(old) || is_txt(cur))) next
    face <- if (!is.null(cur) && is_txt(cur)) cur$face else NULL
    if (is.null(face) && is_txt(old)) face <- old$face
    if (inherits(old, "element_markdown")) {
      th[[nm]] <- ggtext::element_markdown(family = SCI_FONT, size = sci_size_for(nm), face = face)
    } else {
      th[[nm]] <- element_text(family = SCI_FONT, size = sci_size_for(nm), face = face,
                               margin = if (!is.null(cur) && is_txt(cur)) cur$margin else NULL)
    }
  }
  th
}

# Layers, mapped sizes and theme for every plot in a ggplot/patchwork. After
# assembling, add sci_tags() so the panel letters are 10 pt bold.
sci_prepare <- function(p, pt_scale = 0.5) sci_walk(p, function(q) {
  if (!inherits(q, "ggplot")) return(q)
  sci_fix_layers(q, pt_scale) + sci_theme_for(q) + labs(caption = NULL)
})
# tags sit in their own margin so they never overprint a panel title
sci_tags <- function() theme(plot.tag = element_text(family = SCI_FONT, size = SCI_TAG, face = "bold"),
                             plot.tag.location = "margin", plot.tag.position = "topleft")
# dense multi-panel figures: titles and axis titles at the small size too (still
# only two text sizes plus the 10-pt tags)
sci_dense <- function() theme(plot.title = element_text(size = SCI_SMALL, face = "bold"),
                              axis.title = element_text(size = SCI_SMALL),
                              strip.text = element_text(size = SCI_SMALL),
                              legend.title = element_text(size = SCI_SMALL))

# Save one figure: <name>.pdf (cairo, embedded fonts), optionally <name>.svg, plus
# a 300-ppi PNG preview. Returns the paths invisibly.
sci_save <- function(p, name, width = SCI_W2, height, dir, svg = TRUE) {
  stopifnot(width <= SCI_W2 + 1e-6, height <= 9.7)
  dir.create(dir, recursive = TRUE, showWarnings = FALSE)
  pdf_f <- file.path(dir, paste0(name, ".pdf")); svg_f <- file.path(dir, paste0(name, ".svg"))
  png_f <- file.path(dir, paste0(name, "_preview.png"))
  # cairo_pdf straight over an existing file on Google Drive can leave it 0 bytes;
  # write locally, then copy over
  tmp <- tempfile(fileext = ".pdf")
  ggsave(tmp, p, width = width, height = height, units = "in", device = cairo_pdf)
  if (file.exists(pdf_f)) file.remove(pdf_f)
  stopifnot(file.copy(tmp, pdf_f), file.size(pdf_f) > 0)
  # same local-then-copy route for the SVG and preview (a 0-byte Figure2.svg on 2026-11-18)
  put <- function(tmp, dest) { if (file.exists(dest)) file.remove(dest); stopifnot(file.copy(tmp, dest), file.size(dest) > 0) }
  if (svg) { t2 <- tempfile(fileext = ".svg"); ggsave(t2, p, width = width, height = height, units = "in", device = svglite::svglite); put(t2, svg_f) }
  t3 <- tempfile(fileext = ".png"); ggsave(t3, p, width = width, height = height, units = "in", dpi = 300, bg = "white"); put(t3, png_f)
  cat(sprintf("%-40s %.2f x %.2f in\n", name, width, height))
  invisible(c(pdf_f, svg_f, png_f))
}
