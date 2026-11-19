## Per-analyte spread of the location contribution within each class, current
## pipeline vs the ug/g (creatinine x100) build. Reads the per-analyte table
## written by compare_location_unadjusted.R. Each analyte is averaged over the
## child and adult panels first so it counts once.
## Writes Output/Dilution_Unadjusted/location_range_by_class.csv and
## location_by_analyte_wide.csv. Run from the project root.
suppressPackageStartupMessages(library(data.table))
OUT <- "Output/Dilution_Unadjusted"
pa  <- fread(file.path(OUT, "location_contribution_by_analyte.csv"))
pa  <- pa[source %in% c("published", "ugg") & Category %in%
          c("FlameRetardants","PAH","Parabens","PCB","Pesticides","PFAS",
            "Phenols","Phthalates","Smoking","TraceElements","VOCs")]
an  <- pa[, .(R2 = 100 * mean(R2, na.rm = TRUE), panels = paste(sort(unique(panel)), collapse = "+")),
          by = .(Category, Analyte_Code, source)]
wide <- dcast(an, Category + Analyte_Code + panels ~ source, value.var = "R2")
setnames(wide, c("published", "ugg"), c("current_pct", "ugg_pct"))
wide[, change_pts := ugg_pct - current_pct]
fwrite(wide[order(Category, -ugg_pct)], file.path(OUT, "location_by_analyte_wide.csv"))

fmt <- function(x) sprintf("%.1f", x)
rng <- wide[, .(
  analytes      = .N,
  current_mean  = fmt(mean(current_pct)),
  current_median_range = sprintf("%s (%s-%s)", fmt(median(current_pct)), fmt(min(current_pct)), fmt(max(current_pct))),
  ugg_mean      = fmt(mean(ugg_pct)),
  ugg_median_IQR = sprintf("%s (%s-%s)", fmt(median(ugg_pct)), fmt(quantile(ugg_pct, .25)), fmt(quantile(ugg_pct, .75))),
  ugg_range     = sprintf("%s-%s", fmt(min(ugg_pct)), fmt(max(ugg_pct))),
  ugg_top       = paste(head(Analyte_Code[order(-ugg_pct)], 3), collapse = ", "),
  n_ge20_ugg    = sum(ugg_pct >= 20)), by = Category][order(-as.numeric(ugg_mean))]
fwrite(rng, file.path(OUT, "location_range_by_class.csv"))
print(rng)
cat("\nTop 15 analytes, ug/g:\n"); print(wide[order(-ugg_pct)][1:15, .(Category, Analyte_Code, panels, current = round(current_pct,1), ugg = round(ugg_pct,1))])
