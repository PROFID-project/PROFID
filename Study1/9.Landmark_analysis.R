###############################################################################
# Landmark Approach Sensitivity Analysis
#
# Purpose:
#   Sensitivity analysis to the primary time-dependent Cox model, addressing
#   immortal time bias via landmark restriction rather than counting process
#   formatting. Evaluates whether inappropriate ICD shock (FIS) occurring
#   prior to each landmark is associated with all-cause mortality.
#
# Design:
#   - Full pooled cohort (HELIOS, EU-CERT-ICD, ISRAEL-ICD, PROSE-ICD)
#   - Landmark time points: 6, 12, and 24 months post-ICD implantation
#   - At each landmark: cohort restricted to patients alive beyond landmark;
#     exposure (FIS) fixed as any inappropriate shock prior to or at landmark;
#     follow-up measured from landmark time
#   - Multiple imputation: 20 imputations (MICE), pooled via Rubin's rules
#
# Model:
#   Surv(t_landmark, event) ~ FIS_L + Age_10 + LVEF_5 + eGFR_10 +
#     QRS_log1p + bin_beta_blockers + bin_diabetes + bin_stroke_tia +
#     bin_af_atrial_flutter + bin_sex_male + NYHA_grp + strata(DB)
#
#   - NYHA dichotomised (I-II vs III-IV) as covariate for model stability
#     within reduced landmark risk sets
#   - strata(DB) allows cohort-specific baseline hazards
#   - QRS log1p-transformed to address right skew
#
# Outputs:
#   - Forest plots (PNG + PDF) per landmark
#   - HR tables (PDF + HTML) per landmark
#   - Landmark risk table: N at risk, exposed/unexposed, deaths by group
#   - Consistency comparison table: primary Cox vs all landmark estimates
#   - Summary table across all landmarks
###############################################################################


# Required packages ---------------------------------------------------------
pkgs <- c(
  "tidyverse","dplyr" ,"survival", "data.table", "mice", "naniar", "grid","gridExtra",
  "openxlsx", "readxl", "gt", "ggplot2"
)

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))



# =============================================================================

# 0) INPUTS & SETTINGS

# =============================================================================

INFILE  <- "T:/Study_1/Imputed_data/mice_full_object1.rds"

OUTDIR  <- "T:/Study_1/Supplementary_data/Sensitivity"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
msg <- function(...) cat(sprintf(...), "\n")
# =============================================================================

# 1) LOAD IMPUTED OBJECT + VERIFY KEY VARIABLES

# =============================================================================

tt_full <- readRDS(INFILE)

m_full  <- tt_full$m

msg("Full cohort MICE object loaded: m = %d imputations", m_full)

# Quick verification on imputation 1
dt_ref <- as.data.table(complete(tt_full, 1))

required_vars <- c(
  "Time_death_days", "Status_death",
  "Status_FIS",      "Time_FIS_days",
  "NYHA",            "DB",
  "Age",             "LVEF",          "eGFR",
  "QRS",                        # raw QRS — log1p derived inside build_landmark_dt()
  "bin_beta_blockers",           "bin_diabetes",
  "bin_stroke_tia",              "bin_af_atrial_flutter",
  "bin_sex_male"
)
missing_vars <- setdiff(required_vars, names(dt_ref))
if (length(missing_vars) > 0) {
  stop(sprintf("Required variables missing: %s", paste(missing_vars, collapse=", ")))
} else {
  cat("All required variables present.\n")
}
cat("QRS range:", range(as.numeric(dt_ref$QRS), na.rm=TRUE), "\n")
cat("QRS NAs:  ", sum(is.na(dt_ref$QRS)), "\n\n")
rm(dt_ref)

# =============================================================================
# 1) CONSTANTS
# =============================================================================
VAR_TDEATH <- "Time_death_days"
VAR_DEATH  <- "Status_death"
VAR_FIS    <- "Status_FIS"
VAR_TFIS   <- "Time_FIS_days"
VAR_NYHA   <- "NYHA"

LANDMARKS_DAYS <- setNames(c(183, 365, 730), c("6m","12m","24m"))
MIN_EVENTS <- 3

COVARS <- c(
  "Age_10",
  "LVEF_5",
  "eGFR_10",
  "QRS_log1p",
  "bin_beta_blockers",
  "bin_diabetes",
  "bin_stroke_tia",
  "bin_af_atrial_flutter",
  "bin_sex_male"
)

PRETTY_LABELS <- c(
  "FIS_L"                 = "Inappropriate ICD shock by landmark (Yes vs No)",
  "Age_10"                = "Age (per 10 years)",
  "LVEF_5"                = "LVEF (per 5%)",
  "eGFR_10"               = "eGFR (per 10 units)",
  "QRS_log1p"             = "QRS duration (log-transformed(per unit))",
  "bin_beta_blockers"     = "Beta blockers (Yes vs No)",
  "bin_diabetes"          = "Diabetes (Yes vs No)",
  "bin_stroke_tia"        = "Stroke/TIA (Yes vs No)",
  "bin_af_atrial_flutter" = "AF/Atrial flutter (Yes vs No)",
  "bin_sex_male"          = "Sex (Male vs Female)",
  "NYHA_grpIII_IV"        = "NYHA class III-IV (vs I-II)"
)

# NYHA binary covariate; strata(DB) for cohort heterogeneity
RHS <- paste(c("FIS_L", COVARS, "NYHA_grp", "strata(DB)"), collapse = " + ")
FML <- as.formula(paste("Surv(t_landmark, event) ~", RHS))
msg("Model formula: %s\n", deparse(FML))

# =============================================================================
# 2) HELPERS
# =============================================================================
html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&",  "&amp;",  x, fixed = TRUE)
  x <- gsub("<",  "&lt;",   x, fixed = TRUE)
  x <- gsub(">",  "&gt;",   x, fixed = TRUE)
  x <- gsub("\"", "&quot;", x, fixed = TRUE)
  x
}

save_html_table <- function(dt, file, title = NULL, digits = 4) {
  dt <- as.data.table(dt)
  if (nrow(dt) == 0) dt <- data.table(Note="(Empty table)")
  for (j in names(dt)) if (is.numeric(dt[[j]])) dt[[j]] <- round(dt[[j]], digits)
  th   <- paste0("<tr>", paste(sprintf("<th>%s</th>", html_escape(names(dt))), collapse=""), "</tr>")
  rows <- apply(dt, 1, function(r)
    paste0("<tr>", paste(sprintf("<td>%s</td>", html_escape(r)), collapse=""), "</tr>"))
  rows <- paste(rows, collapse="\n")
  ttl  <- if (!is.null(title))
    sprintf("<h2 style='font-family:Arial;'>%s</h2>", html_escape(title)) else ""
  html <- sprintf(
    "<html><head><meta charset='utf-8'></head><body>%s<table border='1' cellpadding='6' cellspacing='0' style='border-collapse:collapse;font-family:Arial;font-size:12px;'>%s\n%s</table></body></html>",
    ttl, th, rows)
  writeLines(html, con=file)
  invisible(file)
}

# Binary NYHA — NA rows NOT dropped (fixes empty dataset bug from original)
make_nyha_grp <- function(x) {
  z <- toupper(trimws(as.character(x)))
  grp <- rep(NA_character_, length(z))
  grp[z %in% c("I","II","1","2")]   <- "I_II"
  grp[z %in% c("III","IV","3","4")] <- "III_IV"
  factor(grp, levels=c("I_II","III_IV"))
}

# DB factor — actual codes in data: CERT, HELS, ISRL, PRSE
make_db_grp <- function(x) {
  z <- toupper(trimws(as.character(x)))
  factor(z, levels=c("CERT","HELS","ISRL","PRSE"))
}

# Unchanged from original
derive_FIS_L <- function(status_fis, time_fis, L_days) {
  s <- suppressWarnings(as.integer(status_fis))
  t <- suppressWarnings(as.numeric(time_fis))
  s[is.na(s) & !is.na(t)] <- 1L
  s[is.na(s) & is.na(t)]  <- 0L
  ambig <- (s == 1L & is.na(t))
  out   <- as.integer(s == 1L & !is.na(t) & t <= L_days)
  out[ambig] <- NA_integer_
  out
}

# build_landmark_dt — structure unchanged from original
# Key changes: NYHA NA rows not dropped; DB added; QRS_log1p derived from QRS
build_landmark_dt <- function(dt, L_days) {
  dt <- as.data.table(copy(dt))
  
  dt <- dt[is.finite(get(VAR_TDEATH))]
  if (nrow(dt) == 0) return(data.table())
  
  dt <- dt[get(VAR_TDEATH) > L_days]
  if (nrow(dt) == 0) return(data.table())
  
  dt[, FIS_L := derive_FIS_L(get(VAR_FIS), get(VAR_TFIS), L_days)]
  dt <- dt[!is.na(FIS_L)]
  if (nrow(dt) == 0) return(data.table())
  
  dt[, t_landmark := as.numeric(get(VAR_TDEATH)) - L_days]
  dt[, event      := as.integer(get(VAR_DEATH) == 1L)]
  dt <- dt[is.finite(t_landmark) & t_landmark > 0]
  if (nrow(dt) == 0) return(data.table())
  
  # NYHA binary covariate — NA rows NOT dropped
  dt[, NYHA_grp := make_nyha_grp(get(VAR_NYHA))]
  
  # DB — only drop rows where DB itself is missing
  dt[, DB := make_db_grp(DB)]
  dt <- dt[!is.na(DB)]
  if (nrow(dt) == 0) return(data.table())
  
  # Rescale continuous covariates
  dt[, Age_10    := as.numeric(Age)  / 10]
  dt[, LVEF_5    := as.numeric(LVEF) / 5]
  dt[, eGFR_10   := as.numeric(eGFR) / 10]
  # QRS_log1p derived here from raw QRS — not in mids object
  # QRS not imputed so identical across all 20 completions
  dt[, QRS_log1p := log1p(as.numeric(QRS))]
  
  keep <- c("t_landmark","event","FIS_L", COVARS, "NYHA_grp","DB")
  dt[, ..keep]
}

# Unchanged from original
tidy_train_pooled <- function(pooled_summary_dt) {
  dt <- as.data.table(copy(pooled_summary_dt))
  dt[, term  := gsub("^(.*?)(?:1)$", "\\1", as.character(term))]
  dt[, label := as.character(term)]
  dt[term %in% names(PRETTY_LABELS), label := PRETTY_LABELS[term]]
  out <- dt[, .(
    term, label,
    HR      = estimate,
    CI_low  = `2.5 %`,
    CI_high = `97.5 %`,
    p_value = p.value
  )]
  out[, HR_CI := sprintf("%.2f (%.2f\u2013%.2f)", HR, CI_low, CI_high)]
  out[, p_fmt := ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))]
  ord_terms <- c("FIS_L", COVARS, "NYHA_grpIII_IV")
  out[, ord  := match(term, ord_terms)]
  setorder(out, ord)
  out[, ord := NULL]
  out[]
}

# =============================================================================
# 3) OUTPUT — unchanged from original
# =============================================================================
save_hr_table_pdf <- function(tidy_dt, file, title) {
  tab <- tidy_dt[, .(Variable=label, `Hazard ratio (95% CI)`=HR_CI, `p-value`=p_fmt)]
  tg  <- gridExtra::tableGrob(
    tab, rows=NULL,
    theme=gridExtra::ttheme_minimal(
      base_size=10,
      core    = list(fg_params=list(hjust=0, x=0.02)),
      colhead = list(fg_params=list(fontface="bold"))
    )
  )
  pdf(file, width=8, height=max(3, 0.33*nrow(tab)+1))
  gridExtra::grid.arrange(tg, top=grid::textGrob(title, gp=grid::gpar(fontface="bold", fontsize=12)))
  dev.off()
}

save_forest_plot_pub <- function(tidy_dt, lm_label) {
  dt <- copy(tidy_dt)
  CLIP_MAX <- 10
  dt[, CI_low_plot  := pmax(CI_low,  1/CLIP_MAX)]
  dt[, CI_high_plot := pmin(CI_high, CLIP_MAX)]
  dt[, label_wrapped := sapply(label, function(s) paste(strwrap(s, width=40), collapse="\n"))]
  dt[, label_wrapped := factor(label_wrapped, levels=rev(unique(label_wrapped)))]
  
  p <- ggplot(dt, aes(x=HR, y=label_wrapped)) +
    geom_vline(xintercept=1, linetype=2, linewidth=0.6) +
    geom_errorbarh(aes(xmin=CI_low_plot, xmax=CI_high_plot), height=0.18, linewidth=0.7) +
    geom_point(size=2.2) +
    scale_x_log10(limits=c(1/CLIP_MAX, CLIP_MAX)) +
    labs(title=sprintf("Landmark %s (full cohort)", lm_label),
         x="Hazard ratio (log scale)", y=NULL) +
    theme_minimal(base_size=12) +
    theme(plot.title=element_text(face="bold"), panel.grid.minor=element_blank())
  
  png_file <- file.path(OUTDIR, sprintf("forest_%s.png", lm_label))
  pdf_file <- file.path(OUTDIR, sprintf("forest_%s.pdf", lm_label))
  ggsave(png_file, p, width=8.2, height=5.3, dpi=300)
  ggsave(pdf_file, p, width=8.2, height=5.3)
  list(png=png_file, pdf=pdf_file)
}

# =============================================================================
# 4) FIT — unchanged from original fit_train(), uses tt_full
# =============================================================================
fit_landmark <- function(L_days, lm_label) {
  fits  <- vector("list", m_full)
  n_fit <- 0
  
  for (i in 1:m_full) {
    d <- build_landmark_dt(mice::complete(tt_full, i), L_days)
    if (nrow(d) == 0) next
    if (sum(d$event, na.rm=TRUE) < MIN_EVENTS) next
    
    f <- try(coxph(FML, data=d, ties="efron"), silent=TRUE)
    if (!inherits(f, "try-error")) {
      fits[[i]] <- f
      n_fit <- n_fit + 1
    }
  }
  
  fits_ok <- fits[!sapply(fits, is.null)]
  if (length(fits_ok) < 2) return(list(ok=FALSE, n_fit=n_fit))
  
  mira_obj        <- list(call=match.call(), analyses=fits_ok)
  class(mira_obj) <- "mira"
  pooled          <- mice::pool(mira_obj)
  raw_dt          <- as.data.table(summary(pooled, conf.int=TRUE, exponentiate=TRUE))
  tidy_dt         <- tidy_train_pooled(raw_dt)
  
  # Return results — saving done outside so save errors never kill model results
  list(ok=TRUE, n_fit=n_fit, tidy=tidy_dt)
}

# =============================================================================
# 5) MAIN LOOP — unchanged from original
# =============================================================================
msg("\nLandmark analysis — full cohort (m = %d)", m_full)
msg("Model: %s", deparse(FML))
msg("Output folder: %s\n", normalizePath(OUTDIR, winslash="/", mustWork=FALSE))

all_res <- list()

for (lm_label in names(LANDMARKS_DAYS)) {
  L <- unname(LANDMARKS_DAYS[[lm_label]])
  msg("=== Landmark %s (L=%d days) ===", lm_label, L)
  
  res <- tryCatch(
    fit_landmark(L, lm_label),
    error=function(e) { msg("ERROR: %s", e$message); list(ok=FALSE, n_fit=0) }
  )
  
  if (isTRUE(res$ok)) {
    fis <- res$tidy[term=="FIS_L"]
    msg("Models pooled : %d / %d",       res$n_fit, m_full)
    msg("FIS_L result  : HR %s  p = %s", fis$HR_CI, fis$p_fmt)
    
    # Save outputs separately — errors here never kill model results
    tryCatch({
      pdf_tab <- file.path(OUTDIR, sprintf("HR_table_%s.pdf", lm_label))
      save_hr_table_pdf(res$tidy, pdf_tab,
                        sprintf("Pooled Cox — full cohort, Landmark %s", lm_label))
      msg("HR table PDF  : %s", pdf_tab)
    }, error=function(e) msg("WARNING: HR table save failed — %s", e$message))
    
    tryCatch({
      fp <- save_forest_plot_pub(res$tidy, lm_label)
      msg("Forest plot   : %s", fp$pdf)
    }, error=function(e) msg("WARNING: Forest plot save failed — %s", e$message))
    
  } else {
    msg("Pooling not possible — %d models fitted", res$n_fit)
  }
  
  d_counts  <- build_landmark_dt(mice::complete(tt_full, 1), L)
  n_risk    <- if (nrow(d_counts) > 0) nrow(d_counts)                         else NA_integer_
  n_exposed <- if (nrow(d_counts) > 0) sum(d_counts$FIS_L == 1L, na.rm=TRUE) else NA_integer_
  n_unexposed    <- if (nrow(d_counts) > 0) sum(d_counts$FIS_L == 0L, na.rm=TRUE)                      else NA_integer_
  n_deaths  <- if (nrow(d_counts) > 0) sum(d_counts$event == 1L, na.rm=TRUE) else NA_integer_
  n_deaths_exp   <- if (nrow(d_counts) > 0) sum(d_counts$event == 1L & d_counts$FIS_L == 1L, na.rm=TRUE) else NA_integer_
  
  n_deaths_unexp <- if (nrow(d_counts) > 0) sum(d_counts$event == 1L & d_counts$FIS_L == 0L, na.rm=TRUE) else NA_integer_
  
  
  all_res[[lm_label]] <- list(
    landmark    = lm_label,
    L_days      = L,
    ok          = isTRUE(res$ok),
    n_fit       = res$n_fit,
    n_risk      = n_risk,
    n_exposed   = n_exposed,
    n_unexposed    = n_unexposed,
    n_deaths    = n_deaths,
    n_deaths_exp   = n_deaths_exp,
    n_deaths_unexp = n_deaths_unexp,
    FIS_HR_CI   = if (isTRUE(res$ok)) res$tidy[term=="FIS_L", HR_CI]   else NA_character_,
    FIS_p       = if (isTRUE(res$ok)) res$tidy[term=="FIS_L", p_fmt]   else NA_character_,
    FIS_HR      = if (isTRUE(res$ok)) res$tidy[term=="FIS_L", HR]      else NA_real_,
    FIS_CI_low  = if (isTRUE(res$ok)) res$tidy[term=="FIS_L", CI_low]  else NA_real_,
    FIS_CI_high = if (isTRUE(res$ok)) res$tidy[term=="FIS_L", CI_high] else NA_real_,
    FIS_p_raw   = if (isTRUE(res$ok)) res$tidy[term=="FIS_L", p_value] else NA_real_
  )
  msg("")
}

# =============================================================================
# 6) SUMMARY TABLE — unchanged from original
# =============================================================================
summary_dt <- rbindlist(lapply(names(all_res), function(lbl) {
  x <- all_res[[lbl]]
  data.table(
    Landmark       = x$landmark,
    L_days         = x$L_days,
    Pooling_ok     = x$ok,
    Models_fitted  = x$n_fit,
    N_at_risk      = x$n_risk,
    N_exposed_FIS  = x$n_exposed,
    N_deaths       = x$n_deaths,
    FIS_HR_95CI    = x$FIS_HR_CI,
    FIS_p_value    = x$FIS_p,
    Stratification = "NYHA_grp (I-II vs III-IV) covariate + strata(DB)"
  )
}), fill=TRUE)

save_html_table(summary_dt,
                file.path(OUTDIR, "SUMMARY_Landmark_fullcohort.html"),
                title="Landmark analysis summary — full cohort")
fwrite(summary_dt, file.path(OUTDIR, "SUMMARY_Landmark_fullcohort.csv"))

msg("Summary saved: %s", file.path(OUTDIR, "SUMMARY_Landmark_fullcohort.html"))
msg("DONE.")


# =============================================================================
# 7) LANDMARK RISK TABLE
#    Counts at each landmark: N at risk, exposed vs unexposed, deaths in each group
#    Addresses reviewer request for deaths among exposed to contextualise HR
# =============================================================================
msg("\n=== Building landmark risk table ===")

risk_table <- rbindlist(lapply(names(all_res), function(lbl) {
  x <- all_res[[lbl]]
  data.table(
    Landmark              = sprintf("%s (%d days)", lbl, x$L_days),
    N_at_risk             = x$n_risk,
    N_exposed_FIS         = x$n_exposed,
    N_unexposed           = x$n_unexposed,
    Deaths_exposed        = x$n_deaths_exp,
    Deaths_unexposed      = x$n_deaths_unexp,
    Death_pct_exposed     = ifelse(!is.na(x$n_deaths_exp) & !is.na(x$n_exposed) & x$n_exposed > 0,
                                   sprintf("%.1f%%", x$n_deaths_exp / x$n_exposed * 100),
                                   NA_character_),
    Death_pct_unexposed   = ifelse(!is.na(x$n_deaths_unexp) & !is.na(x$n_unexposed) & x$n_unexposed > 0,
                                   sprintf("%.1f%%", x$n_deaths_unexp / x$n_unexposed * 100),
                                   NA_character_)
  )
}), fill = TRUE)

setnames(risk_table,
         old = c("Landmark", "N_at_risk", "N_exposed_FIS", "N_unexposed",
                 "Deaths_exposed", "Deaths_unexposed",
                 "Death_pct_exposed", "Death_pct_unexposed"),
         new = c("Landmark", "Patients at risk", "Inappropriate shock, n", "No inappropriate shock, n",
                 "Deaths (shock), n", "Deaths (no shock), n",
                 "Mortality, shock (%)", "Mortality, no shock (%)")
)

save_html_table(risk_table,
                file.path(OUTDIR, "LANDMARK_risk_table.html"),
                title = "Landmark risk table: exposed vs unexposed patients at each landmark time point")
fwrite(risk_table, file.path(OUTDIR, "LANDMARK_risk_table.csv"))
msg("Landmark risk table saved: %s", file.path(OUTDIR, "LANDMARK_risk_table.html"))

# =============================================================================
# 8) CONSISTENCY COMPARISON TABLE: Primary Cox vs Landmark analyses
#    Reads primary model FIS result from CSV
#    Produces side-by-side comparison of FIS HR across all analytical approaches
# =============================================================================
msg("\n=== Building consistency comparison table ===")

INFILE_PRIMARY <- "T:/Study_1/Supplementary_data/Primary/07_primary_model_strata_DB.csv"

primary_row <- tryCatch({
  pm <- fread(INFILE_PRIMARY)
  
  # Identify FIS term — handle common column name variants
  term_col <- intersect(c("term", "Term", "variable", "Variable"), names(pm))[1]
  hr_col   <- intersect(c("HR", "hr", "estimate", "Estimate"),     names(pm))[1]
  lo_col   <- intersect(c("lower", "CI_low", "2.5 %", "ci_low"),   names(pm))[1]
  hi_col   <- intersect(c("upper", "CI_high", "97.5 %", "ci_high"),names(pm))[1]
  p_col    <- intersect(c("p.value", "p_value", "pvalue", "P"),    names(pm))[1]
  
  cat("Primary model columns detected:
")
  cat(sprintf("  term=%s  HR=%s  lower=%s  upper=%s  p=%s\n",
              term_col, hr_col, lo_col, hi_col, p_col))
  print(pm)
  
  # Find FIS row
  fis <- pm[get(term_col) %in% c("FIS_td", "FIS_L", "Status_FIS", "inappropriate_shock")]
  if (nrow(fis) == 0) {
    stop(sprintf("FIS term not found in primary model CSV. Terms present: %s",
                 paste(pm[[term_col]], collapse=", ")))
  }
  
  hr  <- as.numeric(fis[[hr_col]][1])
  lo  <- as.numeric(fis[[lo_col]][1])
  hi  <- as.numeric(fis[[hi_col]][1])
  pv  <- as.numeric(fis[[p_col]][1])
  
  data.table(
    Analysis      = "Primary: time-dependent Cox (full follow-up)",
    N_at_risk     = NA_integer_,
    N_exposed_FIS = NA_integer_,
    N_deaths      = NA_integer_,
    HR            = round(hr, 2),
    CI_low        = round(lo, 2),
    CI_high       = round(hi, 2),
    HR_95CI       = sprintf("%.2f (%.2f–%.2f)", hr, lo, hi),
    p_value       = ifelse(pv < 0.001, "<0.001", sprintf("%.3f", pv)),
    Immortal_time = "Addressed via start-stop counting process",
    NYHA_handling = "4-level factor (I, II, III, IV vs I)",
    Cohort        = "strata(DB)"
  )
}, error = function(e) {
  msg("WARNING: Could not read primary model CSV — %s", e$message)
  msg("         Check path: %s", INFILE_PRIMARY)
  msg("         Primary row will appear as NA in comparison table.")
  data.table(
    Analysis      = "Primary: time-dependent Cox (full follow-up)",
    N_at_risk     = NA_integer_,
    N_exposed_FIS = NA_integer_,
    N_deaths      = NA_integer_,
    HR            = NA_real_,
    CI_low        = NA_real_,
    CI_high       = NA_real_,
    HR_95CI       = NA_character_,
    p_value       = NA_character_,
    Immortal_time = "Addressed via start-stop counting process",
    NYHA_handling = "4-level factor (I, II, III, IV vs I)",
    Cohort        = "strata(DB)"
  )
})

landmark_rows <- rbindlist(lapply(names(all_res), function(lbl) {
  x     <- all_res[[lbl]]
  L_str <- sprintf("%s (%d days)", lbl, x$L_days)
  data.table(
    Analysis      = sprintf("Sensitivity: landmark Cox (%s)", L_str),
    N_at_risk     = x$n_risk,
    N_exposed_FIS = x$n_exposed,
    N_deaths      = x$n_deaths,
    HR            = if (!is.na(x$FIS_HR))      round(x$FIS_HR,      2) else NA_real_,
    CI_low        = if (!is.na(x$FIS_CI_low))  round(x$FIS_CI_low,  2) else NA_real_,
    CI_high       = if (!is.na(x$FIS_CI_high)) round(x$FIS_CI_high, 2) else NA_real_,
    HR_95CI       = x$FIS_HR_CI,
    p_value       = x$FIS_p,
    Immortal_time = sprintf("Addressed via landmark restriction at %s", L_str),
    NYHA_handling = "Binary covariate (I-II vs III-IV)",
    Cohort        = "strata(DB)"
  )
}), fill = TRUE)

comparison_dt <- rbind(primary_row, landmark_rows, fill = TRUE)

save_html_table(
  comparison_dt,
  file.path(OUTDIR, "COMPARISON_primary_vs_landmark.html"),
  title = paste0(
    "Consistency of inappropriate ICD shock exposure effect across analytical approaches — ",
    "Primary time-dependent Cox vs landmark sensitivity analyses"
  )
)
fwrite(comparison_dt, file.path(OUTDIR, "COMPARISON_primary_vs_landmark.csv"))
msg("Comparison table saved: %s", file.path(OUTDIR, "COMPARISON_primary_vs_landmark.html"))

# Console print
msg("\n=== CONSISTENCY CHECK: FIS HR across all analyses ===")
msg("%-52s  %-22s  %s", "Analysis", "HR (95% CI)", "p-value")
msg("%s", strrep("-", 85))
for (i in seq_len(nrow(comparison_dt))) {
  msg("%-52s  %-22s  %s",
      comparison_dt$Analysis[i],
      comparison_dt$HR_95CI[i],
      comparison_dt$p_value[i])
}

msg("\nDONE.")