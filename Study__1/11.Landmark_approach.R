###############################################################################

# LANDMARK ANALYSIS WITH MICE (TRAIN + TEST)

#

# Outputs (per landmark):

#  - TRAIN: pooled HR table (CSV) + forest plot (PNG + PDF)

#  - TEST : C-index, calibration slope, time-dependent AUC (CSV)

#  - FINAL: summary table across landmarks (CSV)

###############################################################################



suppressPackageStartupMessages({
  
  library(mice)
  library(data.table)
  library(survival)
  library(timeROC)
  library(ggplot2)
  library(gridExtra)
  library(grid)
  library(gtable)
})
  




# =============================================================================

# 0) INPUTS & SETTINGS

# =============================================================================

OUTDIR <- "T:/Study_1/Landmark_approach"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

tt_train <- readRDS("T:/Imputed_data/mice_train_object3.rds")

tt_test <- readRDS("T:/Imputed_data/mice_test_object3.rds")

msg <- function(...) cat(sprintf(...), "\n")
m_train <- tt_train$m
m_test  <- tt_test$m

VAR_TDEATH <- "Time_death_days"
VAR_DEATH  <- "Status_death"
VAR_FIS    <- "Status_FIS"
VAR_TFIS   <- "Time_FIS_days"
VAR_NYHA   <- "NYHA"

LANDMARKS_DAYS <- setNames(c(183, 365, 730), c("6m","12m","24m"))

MIN_EVENTS_TRAIN <- 3
MIN_EVENTS_TEST  <- 2

# tdAUC horizons AFTER landmark (days): 3y + 5y
AUC_TIMES  <- c(1095, 1825)
AUC_LABELS <- c("AUC_3y","AUC_5y")

# Calibration settings
CAL_TAU_TARGET <- 1095
CAL_GROUPS     <- 10
CAL_IMP_TRAIN  <- 1
CAL_IMP_TEST   <- 1

# Rescaled covariates for interpretability
COVARS <- c(
  "Age_10",
  "LVEF_5",
  "eGFR_10",
  "bin_beta_blockers",
  "bin_diabetes",
  "bin_stroke_tia",
  "bin_af_atrial_flutter",
  "bin_sex_male"
)

PRETTY_LABELS <- c(
  "FIS_L"                 = "Inappropriate ICD therapy by landmark (Yes vs No)",
  "Age_10"                = "Age (per 10 years)",
  "LVEF_5"                = "LVEF (per 5%)",
  "eGFR_10"               = "eGFR (per 10 units)",
  "bin_beta_blockers"     = "Beta blockers (Yes vs No)",
  "bin_diabetes"          = "Diabetes (Yes vs No)",
  "bin_stroke_tia"        = "Stroke/TIA (Yes vs No)",
  "bin_af_atrial_flutter" = "AF/Atrial flutter (Yes vs No)",
  "bin_sex_male"          = "Sex (Male vs Female)"
)

RHS <- paste(c("FIS_L", COVARS, "strata(NYHA_grp)"), collapse = " + ")
FML <- as.formula(paste("Surv(t_landmark, event) ~", RHS))

# ============================== LOG ==========================================
msg <- function(...) cat(sprintf(...), "\n")

# ============================ HTML HELPERS ===================================

html_escape <- function(x) {
  x <- as.character(x)
  x <- gsub("&", "&amp;", x, fixed = TRUE)
  x <- gsub("<", "&lt;",  x, fixed = TRUE)
  x <- gsub(">", "&gt;",  x, fixed = TRUE)
  x <- gsub("\"","&quot;", x, fixed = TRUE)
  x
}

save_html_table <- function(dt, file, title = NULL, digits = 4) {
  dt <- as.data.table(dt)
  if (nrow(dt) == 0) dt <- data.table(Note="(Empty table)")
  
  for (j in names(dt)) if (is.numeric(dt[[j]])) dt[[j]] <- round(dt[[j]], digits)
  
  th <- paste0("<tr>", paste(sprintf("<th>%s</th>", html_escape(names(dt))), collapse=""), "</tr>")
  rows <- apply(dt, 1, function(r) paste0("<tr>", paste(sprintf("<td>%s</td>", html_escape(r)), collapse=""), "</tr>"))
  rows <- paste(rows, collapse="\n")
  
  ttl <- if (!is.null(title)) sprintf("<h2 style='font-family:Arial;'>%s</h2>", html_escape(title)) else ""
  html <- sprintf(
    "<html><head><meta charset='utf-8'></head><body>%s<table border='1' cellpadding='6' cellspacing='0' style='border-collapse:collapse;font-family:Arial;font-size:12px;'>%s\n%s</table></body></html>",
    ttl, th, rows
  )
  writeLines(html, con=file)
  invisible(file)
}

# ============================ CORE HELPERS ===================================

summarize_metric <- function(x) {
  n <- sum(!is.na(x))
  if (n < 2) return(list(mean=NA_real_, ci=c(NA_real_, NA_real_), n_valid=n))
  mu <- mean(x, na.rm=TRUE); sdv <- sd(x, na.rm=TRUE); se <- sdv/sqrt(n)
  list(mean=mu, ci=c(mu-1.96*se, mu+1.96*se), n_valid=n)
}

make_nyha_grp <- function(x) {
  z <- toupper(trimws(as.character(x)))
  grp <- rep(NA_character_, length(z))
  grp[z %in% c("I","II","1","2")]   <- "I_II"
  grp[z %in% c("III","IV","3","4")] <- "III_IV"
  factor(grp, levels=c("I_II","III_IV"))
}

derive_FIS_L <- function(status_fis, time_fis, L_days) {
  s <- suppressWarnings(as.integer(status_fis))
  t <- suppressWarnings(as.numeric(time_fis))
  # if status missing but time present -> assume status=1 (matches your earlier logic)
  s[is.na(s) & !is.na(t)] <- 1L
  s[is.na(s) & is.na(t)]  <- 0L
  
  ambig <- (s == 1L & is.na(t))
  out <- as.integer(s == 1L & !is.na(t) & t <= L_days)
  out[ambig] <- NA_integer_
  out
}

build_landmark_dt <- function(dt, L_days) {
  dt <- as.data.table(copy(dt))
  
  # must have time-to-death
  dt <- dt[is.finite(get(VAR_TDEATH))]
  if (nrow(dt) == 0) return(data.table())
  
  # alive beyond landmark
  dt <- dt[get(VAR_TDEATH) > L_days]
  if (nrow(dt) == 0) return(data.table())
  
  # exposure at landmark
  dt[, FIS_L := derive_FIS_L(get(VAR_FIS), get(VAR_TFIS), L_days)]
  dt <- dt[!is.na(FIS_L)]
  if (nrow(dt) == 0) return(data.table())
  
  # follow-up and event
  dt[, t_landmark := as.numeric(get(VAR_TDEATH)) - L_days]
  dt[, event := as.integer(get(VAR_DEATH) == 1L)]
  dt <- dt[is.finite(t_landmark) & t_landmark > 0]
  if (nrow(dt) == 0) return(data.table())
  
  # NYHA strata
  dt[, NYHA_grp := make_nyha_grp(get(VAR_NYHA))]
  dt <- dt[!is.na(NYHA_grp)]
  if (nrow(dt) == 0) return(data.table())
  
  # rescale
  dt[, Age_10  := as.numeric(Age)  / 10]
  dt[, LVEF_5  := as.numeric(LVEF) / 5]
  dt[, eGFR_10 := as.numeric(eGFR) / 10]
  
  keep <- c("t_landmark","event","FIS_L", COVARS, "NYHA_grp")
  dt[, ..keep]
}

td_auc_timeROC <- function(time, status, marker, times_days) {
  ok <- is.finite(time) & time > 0 & !is.na(status) & is.finite(marker)
  time <- time[ok]; status <- status[ok]; marker <- marker[ok]
  if (length(times_days) < 1) return(numeric(0))
  if (sum(status == 1L, na.rm=TRUE) < 2) return(rep(NA_real_, length(times_days)))
  
  roc <- try(timeROC::timeROC(
    T=time, delta=status, marker=marker,
    cause=1, weighting="marginal",
    times=times_days, iid=FALSE
  ), silent=TRUE)
  
  if (inherits(roc, "try-error")) return(rep(NA_real_, length(times_days)))
  as.numeric(roc$AUC)
}

tidy_train_pooled <- function(pooled_summary_dt) {
  dt <- as.data.table(copy(pooled_summary_dt))
  dt[, term := gsub("^(.*?)(?:1)$", "\\1", as.character(term))]
  
  dt[, label := as.character(term)]
  dt[term %in% names(PRETTY_LABELS), label := PRETTY_LABELS[term]]
  
  out <- dt[, .(
    term, label,
    HR = estimate,
    CI_low = `2.5 %`,
    CI_high = `97.5 %`,
    p_value = p.value
  )]
  out[, HR_CI := sprintf("%.2f (%.2f–%.2f)", HR, CI_low, CI_high)]
  out[, p_fmt := ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value))]
  
  ord_terms <- c("FIS_L", COVARS)
  out[, ord := match(term, ord_terms)]
  setorder(out, ord)
  out[, ord := NULL]
  out[]
}

# ============================ OUTPUT (TRAIN) =================================

save_hr_table_pdf <- function(tidy_dt, file, title) {
  tab <- tidy_dt[, .(Variable = label, `Hazard ratio (95% CI)` = HR_CI, `p-value` = p_fmt)]
  tg <- tableGrob(
    tab, rows = NULL,
    theme = ttheme_minimal(
      base_size = 10,
      core = list(fg_params = list(hjust = 0, x = 0.02)),
      colhead = list(fg_params = list(fontface = "bold"))
    )
  )
  tg <- gtable_add_rows(tg, heights = grobHeight(textGrob(title)) + unit(3, "mm"), pos = 0)
  tg <- gtable_add_grob(tg, textGrob(title, gp = gpar(fontface="bold", fontsize=12)),
                        t = 1, l = 1, r = ncol(tg))
  
  pdf(file, width = 8, height = max(3, 0.33*nrow(tab) + 1))
  grid::grid.draw(tg)
  dev.off()
}

save_forest_plot_pub <- function(tidy_dt, lm_label) {
  dt <- copy(tidy_dt)
  CLIP_MAX <- 10
  dt[, CI_low_plot  := pmax(CI_low,  1/CLIP_MAX)]
  dt[, CI_high_plot := pmin(CI_high, CLIP_MAX)]
  dt[, label_wrapped := sapply(label, function(s) paste(strwrap(s, width=40), collapse="\n"))]
  dt[, label_wrapped := factor(label_wrapped, levels = rev(unique(label_wrapped)))]
  
  p <- ggplot(dt, aes(x=HR, y=label_wrapped)) +
    geom_vline(xintercept=1, linetype=2, linewidth=0.6) +
    geom_errorbarh(aes(xmin=CI_low_plot, xmax=CI_high_plot), height=0.18, linewidth=0.7) +
    geom_point(size=2.2) +
    scale_x_log10(limits=c(1/CLIP_MAX, CLIP_MAX)) +
    labs(title=sprintf("Pooled hazard ratios (TRAIN) — Landmark %s", lm_label),
         x="Hazard ratio (log scale)", y=NULL) +
    theme_minimal(base_size=12) +
    theme(plot.title = element_text(face="bold"),
          panel.grid.minor = element_blank())
  
  png_file <- file.path(OUTDIR, sprintf("TRAIN_forest_%s.png", lm_label))
  pdf_file <- file.path(OUTDIR, sprintf("TRAIN_forest_%s.pdf", lm_label))
  ggsave(png_file, p, width=8.2, height=5.3, dpi=300)
  ggsave(pdf_file, p, width=8.2, height=5.3)
  list(png=png_file, pdf=pdf_file)
}

fit_train <- function(L_days, lm_label) {
  fits <- vector("list", m_train)
  n_fit <- 0
  
  for (i in 1:m_train) {
    d <- build_landmark_dt(mice::complete(tt_train, i), L_days)
    if (nrow(d) == 0) next
    if (sum(d$event, na.rm=TRUE) < MIN_EVENTS_TRAIN) next
    
    f <- try(coxph(FML, data=d, ties="efron"), silent=TRUE)
    if (!inherits(f, "try-error")) {
      fits[[i]] <- f
      n_fit <- n_fit + 1
    }
  }
  
  fits_ok <- fits[!sapply(fits, is.null)]
  if (length(fits_ok) < 2) {
    return(list(ok=FALSE, n_fit=n_fit))
  }
  
  mira_obj <- list(call=match.call(), analyses=fits_ok)
  class(mira_obj) <- "mira"
  
  pooled <- mice::pool(mira_obj)
  raw_dt <- as.data.table(summary(pooled, conf.int=TRUE, exponentiate=TRUE))
  tidy_dt <- tidy_train_pooled(raw_dt)
  
  pdf_tab <- file.path(OUTDIR, sprintf("TRAIN_HR_table_%s.pdf", lm_label))
  save_hr_table_pdf(tidy_dt, pdf_tab, sprintf("Pooled Cox model (TRAIN) — Landmark %s", lm_label))
  fp <- save_forest_plot_pub(tidy_dt, lm_label)
  
  list(ok=TRUE, n_fit=n_fit, tidy=tidy_dt,
       table_pdf=pdf_tab, forest_png=fp$png, forest_pdf=fp$pdf)
}

# ============================ VALIDATION (TEST) ==============================

validate_test <- function(L_days, lm_label) {
  c_idx <- rep(NA_real_, m_test)
  cal_s <- rep(NA_real_, m_test)
  
  auc_mat <- matrix(NA_real_, nrow=m_test, ncol=length(AUC_TIMES))
  colnames(auc_mat) <- AUC_LABELS
  
  # diagnostics (to explain NA)
  diag <- data.table(imputation=1:m_test,
                     n_test=0L, events_test=0L, max_follow=NA_real_,
                     times_ok="", skipped_reason="")
  
  for (i in 1:m_test) {
    k <- ((i - 1) %% m_train) + 1
    
    d_tr <- build_landmark_dt(mice::complete(tt_train, k), L_days)
    if (nrow(d_tr) == 0 || sum(d_tr$event, na.rm=TRUE) < MIN_EVENTS_TRAIN) {
      diag[i, skipped_reason := "train_empty_or_low_events"]
      next
    }
    
    fit_tr <- try(coxph(FML, data=d_tr, ties="efron", x=TRUE), silent=TRUE)
    if (inherits(fit_tr, "try-error")) {
      diag[i, skipped_reason := "train_fit_error"]
      next
    }
    
    d_te <- build_landmark_dt(mice::complete(tt_test, i), L_days)
    diag[i, `:=`(n_test=nrow(d_te), events_test=sum(d_te$event, na.rm=TRUE))]
    if (nrow(d_te) == 0 || sum(d_te$event, na.rm=TRUE) < MIN_EVENTS_TEST) {
      diag[i, skipped_reason := "test_empty_or_low_events"]
      next
    }
    
    lp <- try(predict(fit_tr, newdata=d_te, type="lp"), silent=TRUE)
    if (inherits(lp, "try-error")) {
      diag[i, skipped_reason := "predict_error"]
      next
    }
    lp <- as.numeric(lp)
    
    cc <- try(survival::concordance(Surv(t_landmark, event) ~ lp, data=d_te), silent=TRUE)
    if (!inherits(cc, "try-error")) c_idx[i] <- as.numeric(cc$concordance)
    
    cal_fit <- try(coxph(Surv(t_landmark, event) ~ lp, data=d_te, ties="efron"), silent=TRUE)
    if (!inherits(cal_fit, "try-error")) cal_s[i] <- as.numeric(coef(cal_fit)[1])
    
    # --- AUC guard (critical) ---
    maxT <- max(d_te$t_landmark, na.rm=TRUE)
    diag[i, max_follow := maxT]
    if (!is.finite(maxT) || maxT <= 0) {
      diag[i, skipped_reason := "bad_followup"]
      next
    }
    
    times_ok <- AUC_TIMES[AUC_TIMES <= maxT]
    diag[i, times_ok := paste(times_ok, collapse=",")]
    
    if (length(times_ok) >= 1) {
      aucs_ok <- td_auc_timeROC(d_te$t_landmark, d_te$event, lp, times_ok)
      auc_mat[i, match(times_ok, AUC_TIMES)] <- aucs_ok
    }
  }
  
  c_sum <- summarize_metric(c_idx)
  s_sum <- summarize_metric(cal_s)
  
  auc_sum <- rbindlist(lapply(seq_along(AUC_LABELS), function(j) {
    s <- summarize_metric(auc_mat[, j])
    data.table(metric=AUC_LABELS[j], mean=s$mean, ci_low=s$ci[1], ci_high=s$ci[2], n_valid=s$n_valid)
  }))
  
  per_imp <- data.table(imputation=1:m_test, c_index=c_idx, cal_slope=cal_s)
  per_imp <- cbind(per_imp, as.data.table(auc_mat))
  
  # HTML outputs
  save_html_table(per_imp,
                  file.path(OUTDIR, sprintf("TEST_metrics_%s.html", lm_label)),
                  title=sprintf("TEST metrics per imputation — Landmark %s", lm_label))
  
  save_html_table(auc_sum,
                  file.path(OUTDIR, sprintf("TEST_tdAUC_summary_%s.html", lm_label)),
                  title=sprintf("TEST time-dependent AUC summary — Landmark %s", lm_label))
  
  save_html_table(diag,
                  file.path(OUTDIR, sprintf("TEST_diagnostics_%s.html", lm_label)),
                  title=sprintf("TEST diagnostics (why NA/skip) — Landmark %s", lm_label),
                  digits=2)
  
  list(c_index=c_sum, cal_slope=s_sum, auc_summary=auc_sum, diagnostics=diag)
}

# ============================ CALIBRATION PLOT ===============================

get_s0_tau <- function(fit, tau) {
  base_sf <- survfit(fit)
  s <- try(summary(base_sf, times=tau)$surv, silent=TRUE)
  if (!inherits(s, "try-error") && length(s) >= 1) {
    s0 <- suppressWarnings(as.numeric(s[length(s)]))
  } else {
    s0 <- NA_real_
  }
  if (!is.finite(s0) || is.na(s0)) s0 <- suppressWarnings(as.numeric(tail(base_sf$surv, 1)))
  if (!is.finite(s0) || is.na(s0)) stop("Calibration: could not obtain baseline survival.")
  s0
}

make_calibration_plot_guaranteed <- function(L_days, lm_label) {
  d_train <- build_landmark_dt(mice::complete(tt_train, CAL_IMP_TRAIN), L_days)
  d_test  <- build_landmark_dt(mice::complete(tt_test,  CAL_IMP_TEST),  L_days)
  
  if (nrow(d_train) == 0 || sum(d_train$event) < 2) stop("Calibration: not enough TRAIN events.")
  if (nrow(d_test)  == 0) stop("Calibration: empty TEST set at this landmark.")
  
  fit <- coxph(FML, data=d_train, ties="efron", x=TRUE)
  lp  <- as.numeric(predict(fit, newdata=d_test, type="lp"))
  
  max_t <- max(d_test$t_landmark, na.rm=TRUE)
  if (!is.finite(max_t) || max_t < 1) stop("Calibration: invalid follow-up on TEST.")
  
  tau <- min(CAL_TAU_TARGET, floor(stats::quantile(d_test$t_landmark, 0.80, na.rm=TRUE)))
  tau <- max(30, min(tau, floor(max_t)))
  
  s0_tau <- get_s0_tau(fit, tau)
  pred_surv <- s0_tau ^ exp(lp)
  pred_risk <- 1 - pred_surv
  
  dt <- as.data.table(d_test)
  dt[, pred_risk := pred_risk]
  
  # grouped KM if feasible, else fallback
  if (sum(dt$event, na.rm=TRUE) >= 2) {
    dt[, r := frank(pred_risk, ties.method="average")]
    dt[, grp := cut(r,
                    breaks=quantile(r, probs=seq(0,1,length.out=CAL_GROUPS+1), na.rm=TRUE),
                    include.lowest=TRUE, labels=FALSE)]
    
    cal_dt <- rbindlist(lapply(sort(unique(dt$grp)), function(g) {
      dd <- dt[grp == g]
      km <- survfit(Surv(t_landmark, event) ~ 1, data=dd)
      s_obs <- summary(km, times=tau)$surv
      if (length(s_obs) == 0) s_obs <- NA_real_
      data.table(group=g, n=nrow(dd), events=sum(dd$event, na.rm=TRUE),
                 pred=mean(dd$pred_risk, na.rm=TRUE),
                 obs=1 - as.numeric(s_obs),
                 method="KM_groups")
    }))
    
    if (nrow(cal_dt[!is.na(obs)]) >= 2) {
      p <- ggplot(cal_dt[!is.na(obs)], aes(x=pred, y=obs)) +
        geom_abline(intercept=0, slope=1, linetype=3, linewidth=0.8) +
        geom_line(linewidth=1) + geom_point(size=2.6) +
        coord_equal(xlim=c(0,1), ylim=c(0,1)) +
        labs(title=sprintf("Calibration plot (TEST) — Landmark %s", lm_label),
             subtitle=sprintf("Observed vs predicted risk by %d days (KM, %d groups)", tau, CAL_GROUPS),
             x=sprintf("Predicted risk by %d days", tau),
             y=sprintf("Observed risk by %d days (KM)", tau)) +
        theme_minimal(base_size=13) +
        theme(plot.title = element_text(face="bold"))
      return(list(plot=p, data=cal_dt, tau_used=tau, method="KM_groups"))
    }
  }
  
  dt[, y_tau := as.integer(event == 1L & t_landmark <= tau)]
  dt2 <- dt[(t_landmark >= tau) | (y_tau == 1L)]
  if (nrow(dt2) < 20) dt2 <- dt
  
  dt2[, r2 := frank(pred_risk, ties.method="average")]
  dt2[, bin := cut(r2,
                   breaks=quantile(r2, probs=seq(0,1,length.out=11), na.rm=TRUE),
                   include.lowest=TRUE, labels=FALSE)]
  
  cal2 <- dt2[, .(n=.N, pred=mean(pred_risk, na.rm=TRUE), obs=mean(y_tau, na.rm=TRUE)),
              by=bin][order(bin)]
  cal2[, method := "Binary_tau"]
  
  p2 <- ggplot(cal2, aes(x=pred, y=obs)) +
    geom_abline(intercept=0, slope=1, linetype=3, linewidth=0.8) +
    geom_line(linewidth=1) + geom_point(size=2.6) +
    coord_equal(xlim=c(0,1), ylim=c(0,1)) +
    labs(title=sprintf("Calibration plot (TEST) — Landmark %s", lm_label),
         subtitle=sprintf("Observed vs predicted risk by %d days (fallback)", tau),
         x=sprintf("Predicted risk by %d days", tau),
         y=sprintf("Observed risk by %d days", tau)) +
    theme_minimal(base_size=13) +
    theme(plot.title = element_text(face="bold"))
  
  list(plot=p2, data=cal2, tau_used=tau, method="Binary_tau")
}

# ============================ RUN ============================================

msg("\nLandmark analysis (Train m=%d, Test m=%d)", m_train, m_test)
msg("Model: %s", deparse(FML))
msg("Output folder: %s\n", normalizePath(OUTDIR, winslash="/", mustWork=FALSE))

all_res <- list()

for (lm_label in names(LANDMARKS_DAYS)) {
  L <- unname(LANDMARKS_DAYS[[lm_label]])
  msg("=== Landmark %s (L=%d days) ===", lm_label, L)
  
  tr <- tryCatch(fit_train(L, lm_label), error=function(e) { msg("TRAIN error: %s", e$message); list(ok=FALSE, n_fit=0) })
  te <- tryCatch(validate_test(L, lm_label), error=function(e) { msg("TEST error: %s", e$message); NULL })
  
  # Calibration plot (don’t kill run)
  cal_obj <- tryCatch(make_calibration_plot_guaranteed(L, lm_label),
                      error=function(e){ msg("Calibration error: %s", e$message); NULL })
  
  if (!is.null(cal_obj)) {
    cal_png  <- file.path(OUTDIR, sprintf("TEST_calibration_%s.png", lm_label))
    cal_pdf  <- file.path(OUTDIR, sprintf("TEST_calibration_%s.pdf", lm_label))
    cal_html <- file.path(OUTDIR, sprintf("TEST_calibration_points_%s.html", lm_label))
    ggsave(cal_png, cal_obj$plot, width=6.8, height=5.2, dpi=300)
    ggsave(cal_pdf, cal_obj$plot, width=6.8, height=5.2)
    save_html_table(as.data.table(cal_obj$data), cal_html,
                    title=sprintf("Calibration points — %s (method=%s, tau=%d)", lm_label, cal_obj$method, cal_obj$tau_used))
  }
  
  if (isTRUE(tr$ok)) {
    fis <- tr$tidy[term=="FIS_L"]
    msg("TRAIN fitted models: %d", tr$n_fit)
    msg("TRAIN FIS_L: %s, p=%s", fis$HR_CI, fis$p_fmt)
    msg("Saved TRAIN: %s | %s", tr$table_pdf, tr$forest_pdf)
  } else {
    msg("TRAIN: pooling not possible (fitted models=%d).", tr$n_fit)
  }
  
  if (!is.null(te)) {
    if (!is.na(te$c_index$mean)) msg("TEST C-index: %.3f (%.3f–%.3f)", te$c_index$mean, te$c_index$ci[1], te$c_index$ci[2]) else msg("TEST C-index: NA")
    if (!is.na(te$cal_slope$mean)) msg("TEST calib slope: %.3f (%.3f–%.3f)", te$cal_slope$mean, te$cal_slope$ci[1], te$cal_slope$ci[2]) else msg("TEST calib slope: NA")
    msg("TEST metrics HTML: %s", file.path(OUTDIR, sprintf("TEST_metrics_%s.html", lm_label)))
    msg("TEST diagnostics HTML: %s", file.path(OUTDIR, sprintf("TEST_diagnostics_%s.html", lm_label)))
  } else {
    msg("TEST: validation failed completely (see console).")
  }
  
  all_res[[lm_label]] <- list(
    L=L, train_ok=isTRUE(tr$ok),
    train_fitted=tr$n_fit,
    cal_tau=ifelse(is.null(cal_obj), NA_integer_, cal_obj$tau_used),
    cal_method=ifelse(is.null(cal_obj), NA_character_, cal_obj$method),
    test_cindex=ifelse(is.null(te), NA_real_, te$c_index$mean),
    test_calslope=ifelse(is.null(te), NA_real_, te$cal_slope$mean)
  )
  
  msg("")
}

summary_dt <- rbindlist(lapply(names(all_res), function(lbl) {
  x <- all_res[[lbl]]
  data.table(
    Landmark=lbl, L_days=x$L,
    Train_ok=x$train_ok,
    Train_fitted_models=x$train_fitted,
    Test_Cindex=ifelse(is.na(x$test_cindex), NA_character_, sprintf("%.3f", x$test_cindex)),
    Test_CalSlope=ifelse(is.na(x$test_calslope), NA_character_, sprintf("%.3f", x$test_calslope)),
    Cal_tau_days=x$cal_tau,
    Cal_method=x$cal_method
  )
}), fill=TRUE)

save_html_table(summary_dt, file.path(OUTDIR, "SUMMARY_Landmark_TrainTest.html"),
                title="SUMMARY — Landmark Train/Test Results")

msg("Saved summary HTML: %s", file.path(OUTDIR, "SUMMARY_Landmark_TrainTest.html"))
msg("DONE.")

