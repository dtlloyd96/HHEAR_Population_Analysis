## ============================================================================
## harmonization_helpers.R
##
## Single source of truth for the three harmonization steps that were previously
## reimplemented in every script, with divergent behaviour:
##
##   1. to_ng_per_mL()             - unit conversion    (was 6 copies)
##   2. std_study()                - cohort key         (was 3 definitions)
##   3. standardize_binary_outcome() - binary coding    (was an alphabetical guess)
##
## Added 2026-08-11 in response to the code audit (CODE_AUDIT_2026-08-11.md
## sections 1.1, 1.2, 1.4, 3.5). Source this file instead of redefining these
## functions:
##
##   source("harmonization_helpers.R")            # from Coding/Paper1_Submission2/
##   source("../harmonization_helpers.R")         # from Vis_Code/
##
## ============================================================================

## ---------------------------------------------------------------------------
## 1. Unit conversion
## ---------------------------------------------------------------------------
##
## REVERTED 2026-08-11 at the author's request. This is the ORIGINAL conversion,
## byte-for-byte in behaviour with the per-file copies it replaced: unrecognised
## units fall through the final `TRUE ~ value` branch and are treated as ng/mL.
##
## The audit finding is NOT fixed here and still stands (see
## CODE_AUDIT_2026-08-11.md section 1.1): mg/dL, ug/mL, ng/g, pmol/mL, nM and
## kg/m3 -- 39,475 rows, 5.9% -- pass through unconverted. Left in place
## deliberately so the units question can be settled separately.

to_ng_per_mL <- function(value, units) {
  dplyr::case_when(
    is.na(units)                  ~ value,          # assume already ng/mL if missing
    units == "ng/mL"              ~ value,
    units == "ug/dL"              ~ value * 10,     # 1 ug/dL = 10 ng/mL
    units %in% c("ng/L")          ~ value / 1000,   # 1 ng/L = 0.001 ng/mL
    units %in% c("pg/mL", "pg/ml")~ value / 1000,   # 1000 pg/mL = 1 ng/mL
    units == "mg/L"               ~ value * 1000,   # 1 mg/L = 1000 ng/mL
    # ng/g RESTORED to passthrough 2026-08-11 (author's decision: keep VIVA in).
    # It had briefly been mapped to NA, which removed VIVA entirely -- VIVA reports
    # nothing but ng/g, all 28,460 rows across 20 trace elements and 1,423
    # participants. That drop cost only 0.9 R2 points on the trace-element
    # category, so VIVA is retained.
    #
    # The caveat stands and is NOT fixed: ng/g is mass-per-mass, a different
    # sample matrix from everyone else's fluid concentrations, and passing it
    # through means VIVA's values are treated as if they were ng/mL. Harmless for
    # within-analyte, within-cohort work such as the Location quantile models;
    # NOT valid anywhere absolute concentrations are compared across cohorts
    # (Figure 3 Panel C's fold-change, the exposure-burden percentiles).
    # Same caveat applies to pmol/mL, pmol/ml and nM, also left as passthrough.
    # See CODE_AUDIT_2026-08-11.md section 1.1.
    TRUE                          ~ value           # unknown unit: pass through
  )
}

## ---------------------------------------------------------------------------
## 2. Cohort key
## ---------------------------------------------------------------------------
##
## Three different definitions were in use:
##   sub("_.*$", "", x)            -> "NHBCS"      (ExWAS, regime diagnostic)
##   str_remove(x, "[_]\\w+|:")    -> "NHBCS"      (power arrays, SDOH intervals)
##   str_remove(x, "_Targeted$")   -> "NHBCS_Cot"  (Figure 3, combined figure)
##
## The third caused NHBCS_Cot's concentrations to fail to join against ExWAS
## results keyed on NHBCS, so that cohort was silently dropped from Figure 3's
## Panel C fold-change.
##
## Canonical choice: strip the FIRST underscore-delimited suffix onwards, i.e.
## the majority behaviour. NHBCS_Cot (the cotinine sub-assay) collapses into
## NHBCS, which is correct for cohort identity -- it is the same cohort, a
## different assay batch.

#' Canonical cohort key. "GOALS_Targeted" -> "GOALS", "NHBCS_Cot_Targeted" -> "NHBCS".
std_study <- function(x) {
  sub("_.*$", "", as.character(x))
}

## ---------------------------------------------------------------------------
## 3. Binary outcome coding
## ---------------------------------------------------------------------------
##
## The previous implementation decided which level was the "case" by sorting the
## levels as CHARACTER and taking the second one. That is an arbitrary rule that
## silently inverted outcomes:
##
##   c(0, 1)          -> 1 is the case   (correct)
##   c("No","Yes")    -> "Yes" is case   (correct by luck)
##   c("Case","Control") -> "Control" is case   (WRONG)
##   c(1, 2)          -> 2 is the case   (WRONG where 1=Yes, 2=No)
##
## Verified instance: Preterm_Birth is coded 0/1 in CEHC and PARENTs but 1/2 in
## CCREOH, so CCREOH's models were estimating the effect on TERM birth while the
## other cohorts estimated PRETERM. Those betas then entered Figure 3, whose
## headline result is cross-cohort sign discordance.
##
## Policy now: recognise only unambiguous codings, plus an explicit registry for
## known conventions. Anything else is a hard error rather than a guess.

# Explicit level -> case mappings for codings that are NOT self-evident.
# Format: list(levels_sorted_as_character = case_level_as_character)
#
# 1/2 is the standard "1 = Yes, 2 = No" questionnaire convention. Applied here
# for CCREOH's Preterm_Birth, where it is consistent with prevalence (81/649 =
# 12.5% coded 1, a plausible preterm rate; 87.5% coded 2).
#
# !! CONFIRM AGAINST THE CCREOH CODEBOOK BEFORE PUBLICATION. If CCREOH in fact
# !! codes 1 = No, invert this entry. The audit could establish that the two
# !! codings are inconsistent but not, from data alone, which direction is right.
BINARY_CODING_REGISTRY <- list(
  "1|2" = "1"   # 1 = Yes/affected, 2 = No/unaffected
)

#' Map a two-level outcome to 0/1 with an explicit, auditable rule.
#'
#' @param y     the outcome vector
#' @param name  outcome name, used in error messages
#' @return integer 0/1 of the same length as the non-missing input
standardize_binary_outcome <- function(y, name = "") {
  y <- y[!is.na(y)]
  uniq <- unique(y)
  if (length(uniq) != 2) {
    stop(sprintf("%s: expected exactly 2 levels for a binary outcome, got %d.",
                 name, length(uniq)))
  }

  # Unambiguous: already 0/1
  if (setequal(as.character(uniq), c("0", "1"))) {
    return(as.integer(as.numeric(as.character(y)) == 1))
  }
  # Unambiguous: logical
  if (is.logical(y)) {
    return(as.integer(y))
  }
  # Unambiguous: conventional yes/no style words
  lower <- tolower(trimws(as.character(uniq)))
  YES <- c("yes", "y", "true", "t", "positive", "pos", "case", "affected", "present")
  NO  <- c("no", "n", "false", "f", "negative", "neg", "control", "unaffected", "absent")
  if (sum(lower %in% YES) == 1 && sum(lower %in% NO) == 1) {
    case_level <- as.character(uniq)[lower %in% YES]
    return(as.integer(as.character(y) == case_level))
  }

  # Registered convention
  key <- paste(sort(as.character(uniq)), collapse = "|")
  if (!is.null(BINARY_CODING_REGISTRY[[key]])) {
    case_level <- BINARY_CODING_REGISTRY[[key]]
    return(as.integer(as.character(y) == case_level))
  }

  stop(sprintf(
    paste0("%s: ambiguous binary coding with levels {%s}. Refusing to guess which ",
           "level is the case -- the previous alphabetical rule silently inverted ",
           "outcomes. Add an entry to BINARY_CODING_REGISTRY in ",
           "harmonization_helpers.R after checking the cohort codebook."),
    name, paste(sort(as.character(uniq)), collapse = ", ")))
}

## ---------------------------------------------------------------------------
## 4. Shared assertions
## ---------------------------------------------------------------------------

#' Collapse a covariate table to one row per participant.
#'
#' Harmonized_SDOH.csv is NOT one row per participant. Some cohorts contribute
#' repeated measures: GOALS is longitudinal and stores one row per visit, so a
#' single Child_PID appears with Child_Age 11.2, 12.1, 13.3, 14.3 and a maternal
#' age that increments alongside it. Measured on the child panel:
#'
#'   11,259 rows / 9,663 unique Child_PID
#'      548 exact duplicate rows
#'    1,320 rows across 272 PIDs genuinely CONFLICT (different covariate values)
#'          -- 1,312 of them in GOALS, the rest in AIRWEIGHS, CHARGE, IMBEDD
#'
#' A cross-sectional ExWAS that leaves these in enters those participants into
#' each model several times, which is pseudo-replication regardless of the
#' sample-year join fixed separately (audit section 1.3).
#'
#' !! POLICY DECISION -- NEEDS AUTHOR SIGN-OFF.
#' !! The default here ("earliest") takes the row with the smallest age, which is
#' !! deterministic and documented but is NOT the same as matching each
#' !! participant's covariates to the visit their biospecimen came from. If the
#' !! exposure's Sample_Collection_Year can be mapped to a visit, that pairing is
#' !! the scientifically correct choice and should replace this.
#' !! The previous behaviour was to keep every row (silent pseudo-replication);
#' !! elsewhere in the codebase the convention is distinct(.keep_all = TRUE),
#' !! which takes whichever row happens to sort first.
#'
#' @param df      covariate table
#' @param id_col  participant id column
#' @param age_col column used to order visits, if present
#' @param policy  "earliest" (default), "latest", or "first_row" (legacy)
dedupe_participant_rows <- function(df, id_col, age_col = NULL,
                                    policy = c("earliest", "latest", "first_row"),
                                    label = id_col) {
  policy <- match.arg(policy)
  n0 <- nrow(df)

  df <- dplyr::distinct(df)                       # drop exact repeats first
  n1 <- nrow(df)

  dup_ids <- unique(df[[id_col]][duplicated(df[[id_col]])])
  if (length(dup_ids) == 0) {
    if (n0 != n1) {
      message(sprintf("[dedupe: %s] removed %d exact duplicate row(s); one row per participant.",
                      label, n0 - n1))
    }
    return(df)
  }

  if (!is.null(age_col) && age_col %in% names(df) && policy != "first_row") {
    ord <- suppressWarnings(as.numeric(df[[age_col]]))
    # participants missing the ordering variable sort last, so a row WITH an age
    # is preferred over one without
    ord[is.na(ord)] <- if (policy == "earliest") Inf else -Inf
    df <- df[order(df[[id_col]], if (policy == "earliest") ord else -ord), , drop = FALSE]
  } else {
    df <- df[order(df[[id_col]]), , drop = FALSE]
  }

  out <- df[!duplicated(df[[id_col]]), , drop = FALSE]

  affected <- if ("Study" %in% names(df)) {
    tab <- table(df$Study[df[[id_col]] %in% dup_ids])
    paste(sprintf("%s(%d)", names(tab), as.integer(tab)), collapse = ", ")
  } else "unknown"

  warning(sprintf(
    paste0("[dedupe: %s] %d participant(s) had multiple covariate rows (%d rows -> %d). ",
           "Kept the '%s' row per participant. Cohorts affected: %s. ",
           "This is a POLICY choice -- see dedupe_participant_rows() in harmonization_helpers.R."),
    label, length(dup_ids), n1, nrow(out), policy, affected), call. = FALSE)
  out
}

#' Fail loudly if a join multiplied rows.
#'
#' The ExWAS covariate join silently inflated the child SDOH table by 58.8%
#' (GOALS +305%) because sample-collection-year was joined on PID alone while
#' 2,355 PIDs carry more than one year.
assert_no_row_multiplication <- function(before, after, what = "join") {
  if (nrow(after) > nrow(before)) {
    stop(sprintf(
      "%s multiplied rows: %d -> %d (+%.1f%%). Aggregate the right-hand table to one row per key first.",
      what, nrow(before), nrow(after), 100 * (nrow(after) / nrow(before) - 1)))
  }
  invisible(after)
}

#' Apply FDR across a whole scan, not within one model's coefficient table.
#'
#' The figure scripts called p.adjust() INSIDE their per-study / per-analyte
#' loops, so the correction ran over the handful of coefficients in a single
#' model rather than across the thousands of analyte x cohort x tau tests in the
#' scan -- under-correcting by orders of magnitude
#' (CODE_AUDIT_2026-08-11.md section 2.8).
#'
#' @param df    bound results for the WHOLE scan (all studies, all analytes)
#' @param p_col p-value column
#' @param by    optional character vector of columns defining separate test
#'   families. `by = "Tau"` treats each quantile as its own family, which is
#'   defensible since each tau answers a different question; the default (NULL)
#'   corrects across everything, which is the most conservative reading.
#' @param alpha threshold used for the Sig_10 indicator
apply_scan_fdr <- function(df, p_col = "P_Value", by = NULL, alpha = 0.10,
                           what = "scan") {
  # WARN, don't stop (2026-08-11). This used to call assert_has_inference(),
  # which errors when every p-value is NA. That is too aggressive here: a figure
  # script such as Location_Compare.Rmd loads ALL cached result sets and runs
  # them through this function, even though only one covariate set reaches the
  # plotted panel. A single stale set with no inference then halted the whole
  # figure.
  #
  # A loud warning still surfaces the condition -- which is the point, since the
  # original failure mode was figures silently emptying -- while letting the sets
  # that do have inference through. Adj_P and Sig_10 come back NA for the ones
  # that do not, which is the honest answer. Call assert_has_inference() directly
  # where a hard stop is wanted.
  # An empty or column-less result set is a legitimate input here: some cached
  # covariate sets produced no estimable models at all (e.g.
  # Mat_Location_Study_Quantile_Results). Hand it back with the columns the
  # callers expect rather than erroring, so one dead set cannot halt a figure.
  if (nrow(df) == 0 || !p_col %in% names(df)) {
    warning(sprintf(
      "[fdr: %s] %s -- returning Adj_P/Sig_10 as NA. Nothing from this set is significance-filtered.",
      what, if (nrow(df) == 0) "result set is EMPTY" else sprintf("column '%s' is missing", p_col)),
      call. = FALSE)
    df[["Adj_P"]] <- NA_real_
    df[["Sig_10"]] <- NA_integer_
    return(df)
  }
  if (all(is.na(df[[p_col]]))) {
    warning(sprintf(
      paste0("[fdr: %s] every value of '%s' is NA -- Adj_P and Sig_10 will be NA. ",
             "This result set was produced without a usable bootstrap (summary_R < 25), ",
             "or its fits did not converge. Anything plotted from it is not significance-filtered."),
      what, p_col), call. = FALSE)
  }

  adjust_block <- function(d) {
    d[["Adj_P"]] <- stats::p.adjust(d[[p_col]], method = "fdr")
    d[["Sig_10"]] <- ifelse(is.na(d[["Adj_P"]]), NA_integer_,
                            as.integer(d[["Adj_P"]] <= alpha))
    d
  }

  if (is.null(by)) {
    out <- adjust_block(df)
    n_fam <- 1L
  } else {
    by <- intersect(by, names(df))
    if (length(by) == 0) {
      out <- adjust_block(df); n_fam <- 1L
    } else {
      out <- df %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(by))) %>%
        dplyr::group_modify(~ adjust_block(.x)) %>%
        dplyr::ungroup()
      n_fam <- dplyr::n_distinct(df[by])
    }
  }

  message(sprintf(
    "[fdr: %s] %d tests in %d famil%s; %d significant at Adj_P <= %.2f (raw p <= %.2f: %d)",
    what, nrow(out), n_fam, ifelse(n_fam == 1, "y", "ies"),
    sum(out$Sig_10 %in% 1L), alpha, alpha, sum(out[[p_col]] <= alpha, na.rm = TRUE)))
  out
}

#' Fail loudly if a results table carries no usable inference.
#'
#' run_censored_quantile_reg() returns NA standard errors and p-values whenever
#' summary_R < 25 (it falls back to fallback_crq_summary). Downstream,
#' p.adjust(NA) -> NA and ifelse(NA <= 0.10, 1, 0) -> NA, so a figure silently
#' empties instead of erroring. Call this after every regression sweep.
assert_has_inference <- function(df, col = "P_Value", what = "results") {
  if (!col %in% names(df)) {
    stop(sprintf("%s: column '%s' is missing entirely.", what, col))
  }
  if (all(is.na(df[[col]]))) {
    stop(sprintf(
      paste0("%s: every value of '%s' is NA. This is what happens when ",
             "run_censored_quantile_reg() is called without summary_R (the default ",
             "of 2 is below min_boot_R = 25, so inference falls back to NA). Pass ",
             "summary_R explicitly."),
      what, col))
  }
  invisible(df)
}
