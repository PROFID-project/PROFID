###############################################################################
# LANDMARK COMPETING-RISK SENSITIVITY (6/12/24m) WITH MICE — Fine & Gray (sHR)
# NYHA EXCLUDED (to avoid sparse strata + ensure pooling works)
#
# TRAIN: pooled Fine–Gray sHRs across imputations (Rubin rules on log(sHR))
# CIF curves (TRAIN + TEST) by exposure for:
#   - Status=1 (Appropriate shock)
#   - Status=2 (Death)
#
# Notes:
# - Survival_time is in MONTHS -> converted to DAYS (30.4375 days/month)
# - Time_FIS_days may be in months or days; heuristic conversion included
# - No CSV saved. Only PDF tables + PNG/PDF plots.
# - Continuous covariates rescaled for interpretability:
#     Age_10  = Age / 10   -> sHR per 10 years
#     LVEF_5  = LVEF / 5   -> sHR per 5% LVEF
#     eGFR_10 = eGFR / 10  -> sHR per 10 eGFR units
###############################################################################

suppressPackageStartupMessages({
  library(mice)
  library(data.table)
  library(ggplot2)
  library(grid)
  library(gtable)
  library(gridExtra)
  library(cmprsk)  # crr, cuminc
})

# -------------------------- USER SETTINGS ------------------------------------

PATH_TRAIN <- "T:/Imputed_data/mice_train_object3.rds"
PATH_TEST  <- "T:/Imputed_data/mice_test_object3.rds"

OUTDIR <- "T:/Study_1/Sensitivity_analysis"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# Which imputations to use for CIF figures (standard: 1)
CIF_IMP_TRAIN <- 1
CIF_IMP_TEST  <- 1

DAYS_PER_MONTH <- 30.4375

# Landmarks in months -> converted to days
LANDMARKS_MONTHS <- setNames(c(6, 12, 24), c("6m","12m","24m"))
LANDMARKS_DAYS   <- round(LANDMARKS_MONTHS * DAYS_PER_MONTH)

MIN_MODELS_FOR_POOL <- 2
MIN_EVENTS_CAUSE1   <- 2

# Variable names
VAR_TIME   <- "Survival_time"   # months (confirmed)
VAR_STATUS <- "Status"          # 0=alive/censored, 1=appropriate shock, 2=death
VAR_FIS    <- "Status_FIS"
VAR_TFIS   <- "Time_FIS_days"

# Covariates (raw names in your data)
COVARS_RAW <- c(
  "Age","LVEF","eGFR",
  "bin_beta_blockers","bin_diabetes","bin_stroke_tia",
  "bin_af_atrial_flutter","bin_sex_male"
)

# Pretty labels (MODEL terms, not raw variable names)
PRETTY_LABELS <- c(
  "FIS_L"                 = "Inappropriate ICD therapy by landmark (Yes vs No)",
  "Age_10"                = "Age (per 10 years)",
  "LVEF_5"                = "LVEF (per 5%)",
  "eGFR_10"               = "eGFR (per 10 units)",
  "bin_beta_blockers1"     = "Beta blockers (Yes vs No)",
  "bin_diabetes1"          = "Diabetes (Yes vs No)",
  "bin_stroke_tia1"        = "Stroke/TIA (Yes vs No)",
  "bin_af_atrial_flutter1" = "AF/Atrial flutter (Yes vs No)",
  "bin_sex_male"          = "Sex (Male vs Female)"
)

# -------------------------- LOAD DATA ----------------------------------------

tt_train <- readRDS(PATH_TRAIN)
tt_test  <- readRDS(PATH_TEST)
stopifnot(inherits(tt_train, "mids"), inherits(tt_test, "mids"))

m_train <- tt_train$m
m_test  <- tt_test$m

cat(sprintf("\nLoaded: tt_train (m=%d), tt_test (m=%d)\n", m_train, m_test))
cat(sprintf("Output folder: %s\n\n", normalizePath(OUTDIR, winslash="/")))

# -------------------------- HELPERS ------------------------------------------

# Heuristic: if Time_FIS_days looks like months (<=120 at 95th percentile), convert to days.
tfis_to_days <- function(tfis_raw) {
  q95 <- suppressWarnings(as.numeric(stats::quantile(tfis_raw, 0.95, na.rm = TRUE)))
  if (is.finite(q95) && q95 <= 120) return(as.numeric(tfis_raw) * DAYS_PER_MONTH)
  as.numeric(tfis_raw)
}

derive_FIS_L <- function(status_fis, time_fis_days, L_days) {
  s <- suppressWarnings(as.integer(status_fis))
  t <- suppressWarnings(as.numeric(time_fis_days))
  
  # both missing => assume no inappropriate therapy recorded
  s[is.na(s) & is.na(t)] <- 0L
  # status missing but time exists => assume event recorded
  s[is.na(s) & !is.na(t)] <- 1L
  
  # ambiguous: status says event but time missing => drop
  ambig <- (s == 1L & is.na(t))
  
  out <- as.integer(s == 1L & !is.na(t) & t <= L_days)
  out[ambig] <- NA_integer_
  out
}

# Build landmark dataset in (days since landmark), with rescaled covariates
build_landmark_dt <- function(df, L_days) {
  dt <- as.data.table(copy(df))
  
  need <- unique(c(VAR_TIME, VAR_STATUS, VAR_FIS, VAR_TFIS, COVARS_RAW))
  miss <- setdiff(need, names(dt))
  if (length(miss)) stop("Missing variables: ", paste(miss, collapse=", "))
  
  # Survival_time is months -> days
  dt[, time_days := as.numeric(get(VAR_TIME)) * DAYS_PER_MONTH]
  dt[, tfis_days := tfis_to_days(get(VAR_TFIS))]
  
  # survive beyond landmark
  dt <- dt[time_days > L_days]
  if (nrow(dt) == 0) return(data.table())
  
  # exposure at landmark
  dt[, FIS_L := derive_FIS_L(get(VAR_FIS), tfis_days, L_days)]
  dt <- dt[!is.na(FIS_L)]
  
  # restart clock at landmark
  dt[, t_landmark := time_days - L_days]
  dt[, Status_fg  := as.integer(get(VAR_STATUS))]  # 0/1/2
  
  # Rescale continuous predictors for interpretable sHRs
  dt[, Age_10  := as.numeric(Age)  / 10]
  dt[, LVEF_5  := as.numeric(LVEF) / 5]
  dt[, eGFR_10 := as.numeric(eGFR) / 10]
  
  keep <- c(
    "t_landmark","Status_fg","FIS_L",
    "Age_10","LVEF_5","eGFR_10",
    "bin_beta_blockers","bin_diabetes","bin_stroke_tia",
    "bin_af_atrial_flutter","bin_sex_male"
  )
  dt[, ..keep]
}

# Fine–Gray via cmprsk::crr; return coef + var diagonal (fit has $var, not vcov())
fit_crr_return <- function(dt) {
  if (sum(dt$Status_fg == 1L, na.rm = TRUE) < MIN_EVENTS_CAUSE1) return(NULL)
  
  X <- model.matrix(
    ~ FIS_L + Age_10 + LVEF_5 + eGFR_10 +
      bin_beta_blockers + bin_diabetes + bin_stroke_tia +
      bin_af_atrial_flutter + bin_sex_male,
    data = dt
  )
  if (ncol(X) >= 1 && colnames(X)[1] == "(Intercept)") X <- X[, -1, drop=FALSE]
  
  fit <- try(cmprsk::crr(
    ftime    = dt$t_landmark,
    fstatus  = dt$Status_fg,
    cov1     = X,
    failcode = 1,
    cencode  = 0
  ), silent = TRUE)
  if (inherits(fit, "try-error")) return(NULL)
  
  b <- fit$coef
  vdiag <- diag(fit$var)
  names(vdiag) <- names(b)
  list(coef = b, vdiag = vdiag)
}

# Rubin pooling for a single coefficient:
# q_i = log(sHR), u_i = var(log(sHR))
pool_rubin <- function(q, u) {
  m <- length(q)
  qbar <- mean(q)
  ubar <- mean(u)
  b <- stats::var(q)
  Tvar <- ubar + (1 + 1/m) * b
  se <- sqrt(Tvar)
  list(est = qbar, se = se)
}

make_pooled_table <- function(coef_list, vdiag_list) {
  terms <- names(coef_list[[1]])
  
  out <- rbindlist(lapply(terms, function(term) {
    q <- vapply(coef_list, function(b) unname(b[term]), numeric(1))
    u <- vapply(vdiag_list, function(v) unname(v[term]), numeric(1))
    if (any(!is.finite(q)) || any(!is.finite(u))) return(NULL)
    
    pu <- pool_rubin(q, u)
    logHR <- pu$est; se <- pu$se
    
    HR <- exp(logHR)
    lo <- exp(logHR - 1.96*se)
    hi <- exp(logHR + 1.96*se)
    
    z <- logHR / se
    p <- 2 * stats::pnorm(abs(z), lower.tail = FALSE)
    
    data.table(term = term, HR = HR, CI_low = lo, CI_high = hi, p_value = p)
  }), fill = TRUE)
  
  if (nrow(out) == 0) return(data.table())
  
  out[, label := term]
  out[term %in% names(PRETTY_LABELS), label := PRETTY_LABELS[term]]
  out[, HR_CI := sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)]
  out[, p_fmt := ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))]
  out[]
}

save_table_pdf <- function(tidy_dt, file, title) {
  tab <- tidy_dt[, .(Variable = label, `sHR (95% CI)` = HR_CI, `p-value` = p_fmt)]
  tg <- tableGrob(
    tab, rows = NULL,
    theme = ttheme_minimal(
      base_size = 11,
      core = list(fg_params = list(hjust = 0, x = 0.02)),
      colhead = list(fg_params = list(fontface = "bold"))
    )
  )
  tg <- gtable_add_rows(tg, heights = grobHeight(textGrob(title)) + unit(3, "mm"), pos = 0)
  tg <- gtable_add_grob(tg, textGrob(title, gp = gpar(fontface="bold", fontsize=13)),
                        t = 1, l = 1, r = ncol(tg))
  
  grDevices::pdf(file, width = 9, height = max(3.2, 0.35*nrow(tab) + 1.2), useDingbats = FALSE)
  grid::grid.draw(tg)
  grDevices::dev.off()
}

save_forest_plot <- function(tidy_dt, title, out_prefix) {
  dt <- copy(tidy_dt)
  
  CLIP_MAX <- 10
  dt[, CI_low_plot  := pmax(CI_low,  1/CLIP_MAX)]
  dt[, CI_high_plot := pmin(CI_high, CLIP_MAX)]
  dt[, label_wrapped := sapply(label, function(s) paste(strwrap(s, width = 46), collapse = "\n"))]
  dt[, label_wrapped := factor(label_wrapped, levels = rev(unique(label_wrapped)))]
  
  p <- ggplot(dt, aes(x = HR, y = label_wrapped)) +
    geom_vline(xintercept = 1, linetype = 2, linewidth = 0.8, color = "grey40") +
    geom_errorbarh(aes(xmin = CI_low_plot, xmax = CI_high_plot),
                   height = 0.18, linewidth = 0.9, color = "grey20") +
    geom_point(size = 2.8, color = "black") +
    scale_x_log10(limits = c(1/CLIP_MAX, CLIP_MAX)) +
    labs(title = title, x = "Subdistribution hazard ratio (log scale)", y = NULL) +
    theme_minimal(base_size = 13) +
    theme(
      plot.title = element_text(face="bold"),
      axis.text.y = element_text(size = 11),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(OUTDIR, paste0(out_prefix, ".png")), p, width = 9.2, height = 5.8, dpi = 350)
  ggsave(file.path(OUTDIR, paste0(out_prefix, ".pdf")), p, width = 9.2, height = 5.8)
}

# CIF plot (colored + renamed legend)
make_cif_plot <- function(dt, title) {
  dt <- as.data.table(copy(dt))
  dt <- dt[is.finite(t_landmark) & t_landmark >= 0 & !is.na(Status_fg) & !is.na(FIS_L)]
  if (nrow(dt) == 0) stop("CIF: empty dataset.")
  
  ci <- cmprsk::cuminc(ftime = dt$t_landmark, fstatus = dt$Status_fg, group = dt$FIS_L)
  nm <- names(ci)
  keep <- nm[grepl("^[01] [12]$", nm)]  # "<group> <cause>"
  
  long <- rbindlist(lapply(keep, function(k) {
    parts <- strsplit(k, " ")[[1]]
    grp <- as.integer(parts[1])
    cause <- as.integer(parts[2])
    obj <- ci[[k]]
    data.table(time = obj$time, cif = obj$est, FIS_L = grp, cause = cause)
  }), fill = TRUE)
  
  long[, Group := ifelse(FIS_L == 1,
                         "Inappropriate ICD therapy: Yes",
                         "Inappropriate ICD therapy: No")]
  long[, Event := ifelse(cause == 1, "Appropriate shock (Status=1)", "Death (Status=2)")]
  
  cols <- c("Inappropriate ICD therapy: No"  = "#1f77b4",
            "Inappropriate ICD therapy: Yes" = "#d62728")
  
  ggplot(long, aes(x = time, y = cif, color = Group)) +
    geom_step(linewidth = 1.05) +
    facet_wrap(~ Event, ncol = 1, scales = "free_y") +
    scale_color_manual(values = cols) +
    labs(
      title = title,
      x = "Time since landmark (days)",
      y = "Cumulative incidence (CIF)",
      color = NULL
    ) +
    theme_minimal(base_size = 14) +
    theme(
      plot.title = element_text(face="bold"),
      legend.position = "right",
      panel.grid.minor = element_blank()
    )
}

save_cif_figures <- function(L_days, lm_label) {
  d_tr <- build_landmark_dt(complete(tt_train, CIF_IMP_TRAIN), L_days)
  if (nrow(d_tr) > 0) {
    p <- make_cif_plot(d_tr, sprintf("CIF (TRAIN) — Landmark %s", lm_label))
    ggsave(file.path(OUTDIR, sprintf("TRAIN_CIF_%s.png", lm_label)), p, width = 7.8, height = 7.8, dpi = 350)
    ggsave(file.path(OUTDIR, sprintf("TRAIN_CIF_%s.pdf", lm_label)), p, width = 7.8, height = 7.8)
  }
  
  d_te <- build_landmark_dt(complete(tt_test, CIF_IMP_TEST), L_days)
  if (nrow(d_te) > 0) {
    p <- make_cif_plot(d_te, sprintf("CIF (TEST) — Landmark %s", lm_label))
    ggsave(file.path(OUTDIR, sprintf("TEST_CIF_%s.png", lm_label)), p, width = 7.8, height = 7.8, dpi = 350)
    ggsave(file.path(OUTDIR, sprintf("TEST_CIF_%s.pdf", lm_label)), p, width = 7.8, height = 7.8)
  }
}

# -------------------------- TRAIN POOLING (NO NYHA) --------------------------

pool_fg_train <- function(L_days, lm_label) {
  coef_list <- list()
  vdiag_list <- list()
  
  for (i in 1:m_train) {
    d <- build_landmark_dt(complete(tt_train, i), L_days)
    if (nrow(d) == 0) next
    
    fit_obj <- fit_crr_return(d)
    if (is.null(fit_obj)) next
    
    coef_list[[length(coef_list)+1]] <- fit_obj$coef
    vdiag_list[[length(vdiag_list)+1]] <- fit_obj$vdiag
  }
  
  if (length(coef_list) < MIN_MODELS_FOR_POOL) {
    cat("  TRAIN: pooling not feasible (too few successful models). Fitting imp=1 only.\n")
    
    d1 <- build_landmark_dt(complete(tt_train, 1), L_days)
    fit1 <- fit_crr_return(d1)
    
    if (is.null(fit1)) {
      cat("  TRAIN imp=1: Fine–Gray not feasible (too few cause-1 events).\n")
      return(invisible(FALSE))
    }
    
    b <- fit1$coef
    se <- sqrt(fit1$vdiag)
    sh <- exp(b)
    lo <- exp(b - 1.96*se)
    hi <- exp(b + 1.96*se)
    p  <- 2*pnorm(abs(b/se), lower.tail = FALSE)
    
    res <- data.table(term = names(b), HR = sh, CI_low = lo, CI_high = hi, p_value = p)
    res[, label := term]
    res[term %in% names(PRETTY_LABELS), label := PRETTY_LABELS[term]]
    res[, HR_CI := sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)]
    res[, p_fmt := ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))]
    
    cat("  TRAIN imp=1 results:\n")
    print(res[, .(label, HR_CI, p_fmt)], row.names = FALSE)
    return(invisible(TRUE))
  }
  
  pooled_dt <- make_pooled_table(coef_list, vdiag_list)
  
  title_tab <- sprintf("Pooled Fine–Gray (TRAIN) — Landmark %s (NYHA excluded)", lm_label)
  pdf_file <- file.path(OUTDIR, sprintf("TRAIN_sHR_table_%s.pdf", lm_label))
  save_table_pdf(pooled_dt, pdf_file, title_tab)
  
  save_forest_plot(
    pooled_dt,
    title = sprintf("Pooled Fine–Gray sHRs (TRAIN) — %s (NYHA excluded)", lm_label),
    out_prefix = sprintf("TRAIN_forest_%s", lm_label)
  )
  
  fis <- pooled_dt[term == "FIS_L"]
  if (nrow(fis) == 1) cat(sprintf("  TRAIN FIS_L sHR: %s, p=%s\n", fis$HR_CI, fis$p_fmt))
  cat(sprintf("  Saved TRAIN table: %s\n", pdf_file))
  invisible(TRUE)
}

# -------------------------- RUN ----------------------------------------------

for (lm_label in names(LANDMARKS_DAYS)) {
  L_days <- unname(LANDMARKS_DAYS[[lm_label]])
  cat(sprintf("=== Landmark %s (L=%d days; %.0f months) ===\n",
              lm_label, L_days, LANDMARKS_MONTHS[[lm_label]]))
  
  pool_fg_train(L_days, lm_label)
  save_cif_figures(L_days, lm_label)
  cat("  Saved CIF figures.\n\n")
}

cat("Done.\n")






