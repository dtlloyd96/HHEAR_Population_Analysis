# capture_plots.R -- run a figure script unchanged and keep the plot objects it saves,
# so Science_Supp_Figures.R can re-export them to Science specifications without
# editing each script. Masks, in the global environment:
#   ggsave()            saves `plot` as <basename>.rds before writing the file
#   pdf(), cairo_pdf()  remember the file being opened
#   print() of a ggplot / patchwork while such a device is open saves the object
# Objects go to Output/Paper1_Submission3/Science_Figures/objects/<basename>.rds
# (a second print to the same device is kept as <basename>__2.rds).
#
#   Rscript Coding/Paper1_Submission2/Vis_Code/capture_plots.R <script.R | report.Rmd>
.CAP_DIR <- normalizePath(file.path(getwd(), "Output/Paper1_Submission3/Science_Figures/objects"), mustWork = FALSE)
dir.create(.CAP_DIR, recursive = TRUE, showWarnings = FALSE)
.cap_state <- new.env()
.cap_save <- function(obj, file) {
  if (is.null(file) || !(inherits(obj, "ggplot") || inherits(obj, "patchwork"))) return(invisible())
  key <- tools::file_path_sans_ext(basename(file))
  n <- (if (is.null(.cap_state[[key]])) 0 else .cap_state[[key]]) + 1
  .cap_state[[key]] <- n
  out <- file.path(.CAP_DIR, paste0(key, if (n > 1) paste0("__", n) else "", ".rds"))
  saveRDS(obj, out); message("[captured] ", basename(out))
  # Plots built by a sourced script refer to helpers and data in the global
  # environment, which saveRDS does not carry; keep the small ones next to the plot.
  g <- mget(setdiff(ls(globalenv(), all.names = FALSE), c("ggsave", "pdf", "cairo_pdf", "dev.off",
            "print.patchwork", "print.ggplot2::ggplot")), envir = globalenv())
  g <- g[vapply(g, function(x) as.numeric(object.size(x)) < 50e6, logical(1))]
  g$.cap_theme <- ggplot2::theme_get()        # scripts that rely on theme_set()
  saveRDS(g, sub("\\.rds$", "__env.rds", out))
}
suppressPackageStartupMessages({ library(ggplot2); library(patchwork) })
# Capture only: figure files are not re-drawn (the originals are already current,
# and drawing the 20-38 in originals is what makes a run slow). Set
# CAPTURE_DRAW=1 to draw as well.
.cap_draw <- nzchar(Sys.getenv("CAPTURE_DRAW"))
ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
  .cap_save(plot, filename); if (.cap_draw) ggplot2::ggsave(filename, plot, ...)
  invisible(filename)
}
pdf <- function(file = NULL, ...) { .cap_state$.file <- file; if (.cap_draw) grDevices::pdf(file, ...) else grDevices::pdf(NULL) }
cairo_pdf <- function(filename = NULL, ...) { .cap_state$.file <- filename; if (.cap_draw) grDevices::cairo_pdf(filename, ...) else grDevices::pdf(NULL) }
dev.off <- function(...) { .cap_state$.file <- NULL; grDevices::dev.off(...) }
# Printing goes through the S3 methods registered by ggplot2 and patchwork, which a
# global-environment definition does not reliably mask, so replace the registered
# methods themselves (the originals are kept and called).
local({
  for (cls in c("patchwork", "ggplot2::ggplot")) {
    ns  <- if (cls == "patchwork") "patchwork" else "ggplot2"
    orig <- getS3method("print", cls, envir = emptyenv())
    wrap <- local({ o <- orig; function(x, ...) {
      .cap_save(x, .cap_state$.file)
      if (.cap_draw || is.null(.cap_state$.file)) o(x, ...) else invisible(x)
    } })
    registerS3method("print", cls, wrap, envir = asNamespace(ns))
  }
})

if (sys.nframe() == 0L) {
  f <- commandArgs(TRUE)[1]
  if (grepl("\\.Rmd$", f)) {
    # knitr's autoprint goes around the print masks, so run the report's code at top
    # level instead, printing visible values as the chunks would
    rf <- file.path(tempdir(), sub("\\.Rmd$", ".R", basename(f)))
    knitr::purl(f, output = rf, quiet = TRUE, documentation = 0)
    owd <- setwd(dirname(f))
    source(rf, print.eval = TRUE, echo = FALSE, max.deparse.length = 0)
    setwd(owd)
  } else source(f, chdir = FALSE)
}
