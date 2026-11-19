# C2b -- calendar-year trends compared ACROSS COUNTRIES
# (Round 2: Reviewer 4, comment 4 -- "comparing them across countries with
#  different regulatory histories")
#
# Part 1 -- TEDDY by country. TEDDY sampled Finland, Sweden and the US (and 9
#   German children, too few) under one protocol, one lab and the same 2012-2017
#   window, so differences in trend between countries are not lab or era
#   artifacts. Censored median regression of log concentration on year, within
#   each country: 100*(exp(b)-1) is % change per year. Heterogeneity: the
#   year x country interaction in a pooled TEDDY model (US reference), tested
#   with a bootstrap Wald test.
# Part 2 -- the same chemical across cohorts in different countries, from the
#   C2 within-cohort slopes (different labs and windows; descriptive only).
#
# Log scale makes the slope invariant to the ug/g creatinine x100.
# Writes Output/Reviewer_Round2/C2b_cross_country/. Run from the project root.
suppressPackageStartupMessages({
  library(data.table); library(dplyr); library(quantreg); library(survival); library(ggplot2); library(patchwork)
})
source("Coding/Paper1_Submission2/harmonization_helpers.R")
source("Coding/Paper1_Submission2/censored_quantile_helpers.R")
set.seed(20261104)
OUT <- "Output/Reviewer_Round2/C2b_cross_country"; dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
CW  <- fread("Input/Harmonized_Datasets/Targeted_Analyte_Category_Mapping.csv")
CLASSES_11 <- c("FlameRetardants","PAH","Parabens","PCB","Pesticides","PFAS",
                "Phenols","Phthalates","Smoking","TraceElements","VOCs")
MIN_N <- 30; R_BOOT <- 200

# ---- Part 1: TEDDY by country --------------------------------------------------
load("Input/Harmonized_Datasets/Harmonized_Targeted.RData")
T <- as.data.table(Targeted_Data_All)[Study == "TEDDY_Targeted" & Analyte_Code %in% CW[Category %in% CLASSES_11, Code]]
yr <- unique(T[!is.na(Sample_Collection_Year), .(Child_PID, Sample_Collection_Year)])[
  , .(Year = mean(Sample_Collection_Year)), by = Child_PID]
P <- as.data.table(build_censored_targeted_panel(as.data.frame(T), "Child_PID", to_ng_per_mL))
e <- new.env(); load("Input/Harmonized_Datasets/Harmonized_SDOH.RData", envir = e)
S <- as.data.table(e$SDOH_Data_All)[Study == "TEDDY_SDOH", .(Child_PID, Country = as.character(Location))]
P <- merge(merge(P, yr, by = "Child_PID"), S, by = "Child_PID")
P <- P[!is.na(Exposure_Obs) & Exposure_Obs > 0 & !is.na(Year) & Country != "Germany"]
P[, `:=`(y = log(Exposure_Obs), d = as.numeric(Detected %in% TRUE))]
P[, Country := factor(Country, c("United_States", "Sweden", "Finland"))]

fit_slope <- function(d) {
  if (nrow(d) < MIN_N || uniqueN(d$Year) < 3) return(NULL)
  f <- tryCatch(if (all(d$d == 1)) rq(y ~ Year, tau = 0.5, data = d)
                else crq(Surv(y, d, type = "left") ~ Year, data = d, method = "Portnoy"),
                error = function(e) NULL)
  if (is.null(f)) return(NULL)
  sm <- tryCatch(if (inherits(f, "rq")) summary(f, se = "boot", R = R_BOOT)$coefficients["Year", 1:2]
                 else { s <- summary(f, taus = 0.5, R = R_BOOT)[[1]]$coefficients; s["Year", c(1, 4)] },
                 error = function(e) NULL)
  if (is.null(sm)) return(NULL)
  b <- unname(sm[1]); se <- unname(sm[2])
  data.table(n = nrow(d), pct_cens = round(100 * mean(d$d == 0), 1), b = b, se = se,
             pct_per_year = 100 * (exp(b) - 1), lo = 100 * (exp(b - 1.96 * se) - 1),
             hi = 100 * (exp(b + 1.96 * se) - 1), p = 2 * pnorm(-abs(b / se)))
}
by_country <- P[, fit_slope(.SD), by = .(Analyte_Code, Country)]

# heterogeneity: year x country interaction, bootstrap Wald (uncensored rq on
# Exposure_Obs when censoring is light; analytes with >20% censored are skipped)
het <- P[, {
  if (mean(d == 0) > 0.2 || uniqueN(Country) < 2 || .N < 3 * MIN_N) NULL else {
    f <- tryCatch(rq(y ~ Year * Country, tau = 0.5, data = .SD), error = function(e) NULL)
    if (is.null(f)) NULL else {
      B <- tryCatch(boot.rq(model.matrix(f), .SD$y, tau = 0.5, R = R_BOOT)$B, error = function(e) NULL)
      ix <- grep("^Year:", names(coef(f)))
      if (is.null(B) || !length(ix)) NULL else {
        b <- coef(f)[ix]; V <- cov(B[, ix, drop = FALSE])
        W <- tryCatch(drop(t(b) %*% solve(V) %*% b), error = function(e) NA_real_)
        list(n = .N, wald = W, df = length(ix), p_heterogeneity = pchisq(W, length(ix), lower.tail = FALSE))
      }
    }
  }
}, by = Analyte_Code]
by_country <- merge(by_country, het[, .(Analyte_Code, p_heterogeneity)], by = "Analyte_Code", all.x = TRUE)
by_country <- merge(by_country, unique(CW[, .(Analyte_Code = Code, Category)]), by = "Analyte_Code")
fwrite(by_country[order(Category, Analyte_Code, Country)], file.path(OUT, "teddy_year_trend_by_country.csv"))

# ---- Part 2: same chemical, cohorts in different countries ----------------------
m <- fread("Output/Reviewer_Round2/C2/year_slopes_median.csv")
country <- c(DCHS = "South Africa", GOC = "Chile", TEDDY = "TEDDY (FI/SE/US)", CCREOH = "Suriname")
m[, Country := fifelse(Cohort %in% names(country), country[Cohort], "US")]
multi <- m[, .(k = uniqueN(Country)), by = Analyte_Code][k >= 2, Analyte_Code]
cross <- m[Analyte_Code %in% multi & is.finite(pct_per_year),
           .(Analyte_Code, Category, Cohort, Panel, Country, years = paste0(first_year, "-", last_year),
             pct_per_year = round(pct_per_year, 1), pct_lo = round(pct_lo, 1), pct_hi = round(pct_hi, 1),
             p = signif(P_Value, 2), direction)]
fwrite(cross[order(Category, Analyte_Code, Country)], file.path(OUT, "cross_cohort_country_slopes.csv"))
agree <- cross[direction != "No significant trend", .(sig_countries = uniqueN(Country),
                 dirs = paste(sort(unique(direction)), collapse = " / ")), by = .(Category, Analyte_Code)]
fwrite(agree[order(Category, Analyte_Code)], file.path(OUT, "cross_cohort_direction_agreement.csv"))

# ---- figure ---------------------------------------------------------------------
SHOW <- c("BPA", "BPS", "TCS", "BP3", "MEPB", "PRPB", "MEHP", "MEHHP", "MECPP", "MBZP", "MNBP", "MIBP",
          "MEP", "PYR1", "NAP1", "NAP2", "Pb", "Cd", "Mn")
cap <- function(x, lim = 40) pmax(pmin(x, lim), -lim)
pal <- c(United_States = "#4C72B0", Sweden = "#DD8452", Finland = "#55A868")
a <- by_country[Analyte_Code %in% SHOW]
a[, Analyte_Code := factor(Analyte_Code, rev(SHOW))]
a[, `:=`(x_c = cap(pct_per_year), lo_c = cap(lo), hi_c = cap(hi))]
pA <- ggplot(a, aes(x_c, Analyte_Code, colour = Country)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo_c, xmax = hi_c), height = 0, alpha = 0.6,
                 position = position_dodge(width = 0.6)) +
  geom_point(position = position_dodge(width = 0.6), size = 1.8) +
  scale_colour_manual(values = pal, labels = c(United_States = "United States")) +
  labs(title = "TEDDY, 2012-2017, by country", x = "% change per year", y = NULL, colour = NULL) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")

cpal <- c("US" = "#4C72B0", "TEDDY (FI/SE/US)" = "#8172B3", "Chile" = "#DD8452",
          "South Africa" = "#C44E52", "Suriname" = "#937860")
b <- cross[Analyte_Code %in% SHOW]
b[, Analyte_Code := factor(Analyte_Code, rev(SHOW))]
b[, `:=`(x_c = cap(pct_per_year), lo_c = cap(pct_lo), hi_c = cap(pct_hi))]
pB <- ggplot(b, aes(x_c, Analyte_Code, colour = Country)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbarh(aes(xmin = lo_c, xmax = hi_c), height = 0, alpha = 0.5,
                 position = position_dodge(width = 0.7)) +
  geom_point(position = position_dodge(width = 0.7), size = 1.6) +
  scale_colour_manual(values = cpal) +
  labs(title = "Same chemical, cohorts in different countries", x = "% change per year", y = NULL, colour = NULL) +
  theme_minimal(base_size = 9) + theme(legend.position = "bottom")

# panel C of Figure S8 (assembled by C2_S8_assemble.R): eight key chemicals.
# The cross-cohort panel (pB here) is no longer in S8 (2026-11-04); its result
# is reported in the text and in cross_cohort_country_slopes.csv.
S8C <- c("BPA", "BPS", "MEHHP", "MECPP", "MBZP", "BP3", "PYR1", "Pb")
ac <- a[Analyte_Code %in% S8C]
ac[, Analyte_Code := factor(as.character(Analyte_Code), rev(S8C))]
pC_s8 <- ggplot(ac, aes(x_c, Analyte_Code, colour = Country)) +
  geom_vline(xintercept = 0, linetype = 2, colour = "grey50") +
  geom_errorbar(aes(xmin = lo_c, xmax = hi_c), width = 0, alpha = 0.7,
                position = position_dodge(width = 0.6), orientation = "y") +
  geom_point(position = position_dodge(width = 0.6), size = 2.4) +
  scale_colour_manual(values = pal, labels = c(United_States = "United States")) +
  labs(title = "TEDDY 2012-2017, by country", x = "% change per year", y = NULL, colour = NULL) +
  theme_minimal(base_size = 10) + theme(legend.position = "bottom")
saveRDS(list(pC = pC_s8), file.path(OUT, "S8_panels_CD.rds"))
fig <- pA + pB + plot_annotation(tag_levels = "A")
ggsave(file.path(OUT, "SuppFig_Cross_Country_Year_Trends.pdf"), fig, width = 11, height = 7)
ggsave(file.path(OUT, "SuppFig_Cross_Country_Year_Trends.png"), fig, width = 11, height = 7, dpi = 200)

# ---- console ----------------------------------------------------------------------
options(width = 200)
w <- dcast(by_country, Category + Analyte_Code + p_heterogeneity ~ Country,
           value.var = "pct_per_year", fun.aggregate = function(x) round(x[1], 1))
cat("TEDDY % change per year by country (median), p for year x country:\n")
print(w[order(p_heterogeneity)][, p_heterogeneity := signif(p_heterogeneity, 2)][])
cat("\nTEDDY: analytes with significant (p<0.05) trends, by country and direction:\n")
print(by_country[p < 0.05, .(n = .N, analytes = paste(sort(Analyte_Code), collapse = ",")),
                 by = .(Country, dir = fifelse(b < 0, "down", "up"))])
cat("\nCross-cohort: chemicals with significant trends in >=2 countries:\n")
print(agree[sig_countries >= 2])
