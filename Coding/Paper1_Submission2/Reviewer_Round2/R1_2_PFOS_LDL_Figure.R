# R1 comment 2, option (2): an exposure-outcome relationship that appears only
# after accounting for population structure -- PFOS and LDL cholesterol in
# children.
#
# Found in Cohort_Regime_Diagnostic (cohort_regime_pairs.csv): cohort-adjusted
# (within-cohort, study fixed effects) slope +4.55, t = 3.5, all cohorts
# positive; naive pooled slope +0.63, t = 0.7. Rebuilt here from the same
# inputs and the same prep as cohort_regime_diagnostic.R (exposure ng/mL ->
# log1p, mean per study x child; outcome from Harmonized_Outcomes.csv;
# unadjusted), and asserted to reproduce those numbers.
#
# A  the data: one colour per cohort, each cohort's own fitted line, cohort
#    means (large points), and the naive pooled line (black, dashed)
# B  the estimates: each cohort's slope, the cohort-adjusted estimate and the
#    naive pooled estimate, with 95% CIs
# Writes Output/Reviewer_Round2/R1_2_PFOS_LDL/. Run from the project root.
suppressPackageStartupMessages({ library(data.table); library(ggplot2); library(patchwork) })
source("Coding/Paper1_Submission2/harmonization_helpers.R")
OUT <- "Output/Reviewer_Round2/R1_2_PFOS_LDL"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
IN  <- "Health_Outcome_Analysis/Input"

tg <- fread(file.path(IN, "Harmonized_Targeted.csv"), encoding = "Latin-1",
            colClasses = list(character = c("Child_PID", "Mat_PID")))
tg <- tg[Analyte_Code == "PFOS" & !is.na(Child_PID)]
tg[, `:=`(Study = std_study(Study), conc = to_ng_per_mL(Concentration, Units))]
ex <- tg[!is.na(conc) & conc >= 0, .(x = mean(log1p(conc))), by = .(Study, PID = Child_PID)]

oc <- fread(file.path(IN, "Harmonized_Outcomes.csv"), colClasses = list(character = c("Child_PID", "Mat_PID")))
oc <- oc[!is.na(Child_PID), .(Study = std_study(Study), PID = Child_PID,
                              y = suppressWarnings(as.numeric(Cholesterol_LDL)))][!is.na(y)]
d <- merge(ex, oc, by = c("Study", "PID"))
d <- d[, if (.N >= 10 && uniqueN(x) > 1) .SD, by = Study]

# estimates
per <- d[, { f <- lm(y ~ x); s <- summary(f)$coefficients["x", ]
             .(n = .N, est = s[1], se = s[2]) }, by = Study]
fw <- lm(y ~ x + Study, data = d); sw <- summary(fw)$coefficients["x", ]
fn <- lm(y ~ x, data = d);         sn <- summary(fn)$coefficients["x", ]
# within-between decomposition: y ~ Bw (x - cohort mean) + Bb cohort mean. Without
# covariates Bw equals the cohort-adjusted slope, so panel C labels that row "= Bw"
# rather than drawing it twice; asserted here.
d[, xbar := mean(x), by = Study]
fwb <- lm(y ~ I(x - xbar) + xbar, data = d); swb <- summary(fwb)$coefficients["I(x - xbar)", ]
stopifnot(abs(swb[1] - sw[1]) < 1e-8)
# cohort means with 95% CIs for PFOS (x) and LDL (y)
means <- d[, .(x = mean(x), y = mean(y), n = .N, se_x = sd(x) / sqrt(.N), se_y = sd(y) / sqrt(.N)), by = Study]
means[, `:=`(x_lo = x - 1.96 * se_x, x_hi = x + 1.96 * se_x, y_lo = y - 1.96 * se_y, y_hi = y + 1.96 * se_y)]
cat(sprintf("cohorts %d, N %d | within %.2f (t %.2f) | naive %.2f (t %.2f)\n",
            nrow(per), nrow(d), sw[1], sw[3], sn[1], sn[3]))
stopifnot(abs(sw[1] - 4.55) < 0.05, abs(sn[1] - 0.63) < 0.05)   # = cohort_regime_pairs.csv

est <- rbind(
  per[, .(label = sprintf("%s (n = %d)", Study, n), est, lo = est - 1.96 * se, hi = est + 1.96 * se, type = "Single cohort")],
  data.table(label = "Cohort-adjusted (= Bw, all cohorts)", est = sw[1], lo = sw[1] - 1.96 * sw[2], hi = sw[1] + 1.96 * sw[2], type = "Cohort-adjusted"),
  data.table(label = "Naive pooled (all cohorts)",    est = sn[1], lo = sn[1] - 1.96 * sn[2], hi = sn[1] + 1.96 * sn[2], type = "Naive pooled"))
est[, label := factor(label, rev(label))]
fwrite(est, file.path(OUT, "pfos_ldl_estimates.csv"))
fwrite(means, file.path(OUT, "pfos_ldl_cohort_means.csv"))

pal <- setNames(c("#4C72B0", "#DD8452", "#55A868", "#8172B3")[seq_len(nrow(per))], per$Study)
xr <- range(d$x)
# naive pooled fit, and the between-cohort trend through the cohort means
# (n-weighted) -- the two reference lines that explain the cancellation
fb <- lm(y ~ x, data = means, weights = n)
ref_lines <- rbind(
  data.table(x = xr, y = coef(fn)[1] + coef(fn)[2] * xr, line = "Naive pooled fit"),
  data.table(x = range(means$x), y = coef(fb)[1] + coef(fb)[2] * range(means$x), line = "Between-cohort trend (cohort means)"))
cat(sprintf("between-cohort slope (n-weighted, cohort means) %.2f\n", coef(fb)[2]))

pA <- ggplot(d, aes(x, y, colour = Study)) +
  geom_point(alpha = 0.45, size = 1.5, stroke = 0) +
  geom_smooth(method = "lm", formula = y ~ x, se = FALSE, linewidth = 1.1) +
  geom_line(data = ref_lines, aes(x, y, linetype = line, colour = line), inherit.aes = FALSE,
            linewidth = 1.3) +
  geom_errorbar(data = means, aes(x = x, ymin = y_lo, ymax = y_hi), inherit.aes = FALSE,
                colour = "black", width = 0.12, linewidth = 0.8) +
  geom_errorbar(data = means, aes(y = y, xmin = x_lo, xmax = x_hi), inherit.aes = FALSE,
                colour = "black", width = 2.5, linewidth = 0.8, orientation = "y") +
  geom_point(data = means, aes(x, y, fill = Study, shape = "Cohort mean"), inherit.aes = FALSE,
             size = 4.5, colour = "black", stroke = 0.8) +
  scale_colour_manual(values = c(pal, "Naive pooled fit" = "black",
                                 "Between-cohort trend (cohort means)" = "grey55"),
                      breaks = names(pal), name = "Cohort (own fit)") +
  scale_fill_manual(values = pal, guide = "none") +
  scale_linetype_manual(values = c("Naive pooled fit" = "22", "Between-cohort trend (cohort means)" = "solid"), name = NULL) +
  scale_shape_manual(values = c("Cohort mean" = 21), name = NULL) +
  guides(colour = guide_legend(order = 1, override.aes = list(alpha = 1, size = 2)),
         linetype = guide_legend(order = 2, override.aes = list(linewidth = 1.1,
                                 colour = c("grey55", "black"))),
         shape = guide_legend(order = 3, override.aes = list(fill = "grey70"))) +
  coord_cartesian(ylim = quantile(d$y, c(0.05, 0.95))) +
  labs(title = "Child LDL cholesterol vs PFOS, by cohort",
       x = "PFOS, log(1 + ng/mL)", y = "LDL cholesterol (mg/dL)") +
  theme_classic(base_size = 11) +
  theme(legend.position = "bottom", legend.box = "vertical", legend.margin = margin(0, 0, 0, 0))

est[, col := fifelse(type == "Single cohort", pal[sub(" .*$", "", as.character(label))],
              fifelse(type == "Cohort-adjusted", "#C44E52", "black"))]
pB <- ggplot(est, aes(est, label)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey60") +
  geom_errorbar(aes(xmin = lo, xmax = hi, colour = col), width = 0.15, linewidth = 0.8, orientation = "y") +
  geom_point(aes(colour = col), size = 3) +
  scale_colour_identity() +
  labs(title = "Slope estimates", x = "Change in LDL (mg/dL) per unit log(1 + PFOS)", y = NULL) +
  theme_classic(base_size = 11)

# panels reused in the combined Figure S17 (R1_2_S17_assemble.R)
saveRDS(list(pA = pA, pB = pB), file.path(OUT, "S17_panels_pfos_ldl.rds"))
fig <- pA + pB + plot_layout(widths = c(1.3, 1)) + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT, "Fig_PFOS_LDL_population_structure.pdf"), fig, width = 11, height = 5.6)
ggsave(file.path(OUT, "Fig_PFOS_LDL_population_structure.png"), fig, width = 11, height = 5.6, dpi = 200)
print(est); print(means)
