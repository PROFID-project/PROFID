###############################################################################
# SENSITIVITY ANALYSIS 3 (SAP)
# TIME-DEPENDENT PROPENSITY SCORE (MSM-IPTW) FOR FIS_td
# WITH PS MODEL CALIBRATION ON TEST SET
###############################################################################

suppressPackageStartupMessages({
  library(mice)
  library(data.table)
  library(survival)
  library(ggplot2)
  install.packages("cobalt")
  library(cobalt)
  install.packages("webshot2")
  library(webshot2)
  library(grid)
  
})

# -------------------------------------------------------------------------
# 0) Paths / Inputs
# -------------------------------------------------------------------------
OUTDIR <- "T:/FINAL ICD COHORT/Sensitivity analysis"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)



# Provide your mids objects here (recommended) ------------------------------
TRAIN_RDS <- "T:/Imputed_data/mice_train_object3.rds"  # <-- CHANGE
TEST_RDS  <- "T:/Imputed_data/mice_test_object3.rds"  # <-- CHANGE

t1_train <- readRDS(TRAIN_RDS)
t1_test  <- readRDS(TEST_RDS)


m_train <- t1_train$m
m_test  <- t1_test$m

EPS <- 1e-7
msg <- function(...) cat(sprintf(...), "\n")



# -------------------------------------------------------------------------

# 1) Primary model covariates (your final TD-Cox variables)

# -------------------------------------------------------------------------

primary_covars <- c(
  
  "Age",
  
  "LVEF",
  
  "eGFR",
  
  "bin_beta_blockers",
  
  "bin_diabetes",
  
  "bin_stroke_tia",
  
  "bin_af_atrial_flutter",
  
  "bin_sex_male",
  
  "NYHA"  # included in PS; stratified in outcome Cox
  
)



# -------------------------------------------------------------------------

# 2) Sanity checks: ensure needed columns exist

# -------------------------------------------------------------------------

dt_ref <- as.data.table(complete(t1_train, 1))



need_core <- c("ID","Time_death_days","Status_death","Status_FIS","Time_FIS_days")

miss_core <- setdiff(need_core, names(dt_ref))

if (length(miss_core)) stop("Missing CORE columns in train imp1: ", paste(miss_core, collapse=", "))



miss_cov <- setdiff(primary_covars, names(dt_ref))

if (length(miss_cov)) {
  
  msg("Dropping primary covariates not present: %s", paste(miss_cov, collapse=", "))
  
  primary_covars <- setdiff(primary_covars, miss_cov)
  
}



stopifnot(length(primary_covars) >= 3)



# -------------------------------------------------------------------------

# 3) Helpers: standardise types + fix exposure missingness

# -------------------------------------------------------------------------

standardise_types <- function(dt){
  
  dt <- as.data.table(dt)
  
  
  
  # outcome
  
  dt[, Time_death_days := as.numeric(Time_death_days)]
  
  dt[, Status_death    := as.integer(Status_death)]
  
  
  
  # exposure
  
  dt[, Status_FIS    := as.integer(Status_FIS)]
  
  dt[, Time_FIS_days := as.numeric(Time_FIS_days)]
  
  
  
  # bin_* as integer
  
  for (v in grep("^bin_", names(dt), value = TRUE)) dt[, (v) := as.integer(get(v))]
  
  
  
  # numeric covariates (as applicable)
  
  for (v in intersect(names(dt), c("Age","LVEF","eGFR"))) dt[, (v) := as.numeric(get(v))]
  
  
  
  # NYHA: factor for stratification
  
  if ("NYHA" %in% names(dt)) dt[, NYHA := factor(NYHA)]
  
  
  
  dt
  
}



fix_exposure_na <- function(dt){
  
  dt <- as.data.table(dt)
  
  
  
  # A: both missing => assume not exposed
  
  dt[is.na(Status_FIS) & is.na(Time_FIS_days), Status_FIS := 0L]
  
  
  
  # B: time present but status missing => exposed
  
  dt[is.na(Status_FIS) & !is.na(Time_FIS_days), Status_FIS := 1L]
  
  
  
  # C: exposed but time missing => ambiguous => exclude from TD analysis
  
  dt[, excl_fis_ambig := (Status_FIS == 1L & is.na(Time_FIS_days))]
  
  
  
  # If not exposed, time should be NA
  
  dt[Status_FIS == 0L, Time_FIS_days := NA_real_]
  
  
  
  # If exposure time >= death time -> treat as unexposed (no exposed person-time)
  
  dt[!is.na(Time_FIS_days) & !is.na(Time_death_days) & Time_FIS_days >= Time_death_days,
     
     `:=`(Status_FIS = 0L, Time_FIS_days = NA_real_)]
  
  
  
  dt
  
}



# -------------------------------------------------------------------------

# 4) TD dataset via tmerge

# -------------------------------------------------------------------------

make_td_tmerge <- function(dt, vars_base){
  
  dt <- as.data.table(dt)
  
  
  
  # require outcome
  
  dt <- dt[!is.na(Time_death_days) & !is.na(Status_death)]
  
  
  
  # exclude ambiguous exposure timing
  
  dt <- dt[excl_fis_ambig != TRUE]
  
  
  
  # enforce 1 row per ID for baseline
  
  dt <- dt[order(ID)][, .SD[1], by = ID]
  
  
  
  base_cols <- c("ID","Time_death_days","Status_death", vars_base)
  
  data1 <- dt[, ..base_cols]
  
  
  
  d2_cols <- c("ID","Time_death_days","Status_death", vars_base, "Status_FIS","Time_FIS_days")
  
  data2 <- dt[, ..d2_cols]
  
  
  
  td <- tmerge(
    
    data1 = data1,
    
    data2 = data2,
    
    id    = ID,
    
    death = event(Time_death_days, Status_death == 1),
    
    FIS_td = tdc(Time_FIS_days)
    
  )
  
  
  
  td <- as.data.table(td)
  
  td[, tstop := pmax(tstop, EPS)]
  
  td <- td[tstop > tstart]
  
  setorder(td, ID, tstart, tstop)
  
  td
  
}



# -------------------------------------------------------------------------

# 5) Rubin pooling (single coefficient)

# -------------------------------------------------------------------------

pool_1coef <- function(q, u){
  
  m <- length(q)
  
  qbar <- mean(q)
  
  ubar <- mean(u)
  
  b <- stats::var(q)
  
  tvar <- ubar + (1 + 1/m) * b
  
  se <- sqrt(tvar)
  
  z <- qbar / se
  
  p <- 2 * pnorm(abs(z), lower.tail = FALSE)
  
  list(est=qbar, se=se, p=p)
  
}



# -------------------------------------------------------------------------

# 6) Time-dependent PS (initiation) + stabilized IPTW

#    - We model initiation only (switch 0->1) among those not yet exposed

#    - Denominator: primary covars + time

#    - Numerator  : time only (stabilization)

# -------------------------------------------------------------------------

make_initiation_long <- function(td){
  
  td <- as.data.table(td)
  
  setorder(td, ID, tstart, tstop)
  
  
  
  td[, A_start := shift(FIS_td, n=1L, type="lag", fill=0L), by=ID]
  
  td[, A_stop  := FIS_td]
  
  td[, switch  := as.integer(A_start == 0L & A_stop == 1L)]
  
  
  
  ps_dt <- td[A_start == 0L]           # at-risk for initiation
  
  ps_dt[, t_mid := (tstart + tstop)/2] # time term
  
  ps_dt
  
}



compute_sw_one_imp <- function(td, denom_vars, numer_time_only = TRUE, trunc_q=c(0.01,0.99)){
  
  td <- as.data.table(td)
  
  ps_dt <- make_initiation_long(td)
  
  
  
  if (nrow(ps_dt) == 0L) {
    
    td[, sw := 1]
    
    td[, sw_trunc := 1]
    
    return(list(td=td, ps_model_d=NULL, ps_model_n=NULL,
                
                w_info=data.table(sw_mean=1, sw_sd=0, sw_min=1, sw_max=1,
                                  
                                  swt_mean=1, swt_sd=0, swt_min=1, swt_max=1)))
    
  }
  
  
  
  # Denominator PS model: switch ~ X + time
  
  f_d <- as.formula(paste("switch ~", paste(c(denom_vars, "t_mid"), collapse=" + ")))
  
  
  
  # Numerator PS model: time-only stabilization (recommended simple)
  
  f_n <- if (numer_time_only) {
    
    as.formula("switch ~ t_mid")
    
  } else {
    
    # if you want same as denom (not recommended), set numer_time_only=FALSE
    
    as.formula(paste("switch ~", paste(c(denom_vars, "t_mid"), collapse=" + ")))
    
  }
  
  
  
  ps_model_d <- try(glm(f_d, data=ps_dt, family=binomial(),
                        
                        control=glm.control(maxit=50)), silent=TRUE)
  
  ps_model_n <- try(glm(f_n, data=ps_dt, family=binomial(),
                        
                        control=glm.control(maxit=50)), silent=TRUE)
  
  
  
  if (inherits(ps_model_d, "try-error") || inherits(ps_model_n, "try-error")) {
    
    return(NULL)
    
  }
  
  
  
  ps_dt[, p_d := pmin(pmax(predict(ps_model_d, type="response"), EPS), 1 - EPS)]
  
  ps_dt[, p_n := pmin(pmax(predict(ps_model_n, type="response"), EPS), 1 - EPS)]
  
  
  
  # Contribution per interval for initiation hazard
  
  ps_dt[, contrib_d := ifelse(switch==1L, p_d, 1 - p_d)]
  
  ps_dt[, contrib_n := ifelse(switch==1L, p_n, 1 - p_n)]
  
  
  
  # Person-level stabilized weight = product(numer)/product(denom)
  
  sw_id <- ps_dt[, .(
    
    sw = prod(contrib_n) / prod(contrib_d),
    
    n_at_risk_intervals = .N,
    
    n_switch = sum(switch)
    
  ), by=ID]
  
  
  
  td <- merge(td, sw_id[, .(ID, sw)], by="ID", all.x=TRUE)
  
  td[is.na(sw), sw := 1]
  
  
  
  lo <- as.numeric(quantile(td$sw, probs=trunc_q[1], na.rm=TRUE))
  
  hi <- as.numeric(quantile(td$sw, probs=trunc_q[2], na.rm=TRUE))
  
  td[, sw_trunc := pmin(pmax(sw, lo), hi)]
  
  
  
  w_info <- td[, .(
    
    sw_mean = mean(sw, na.rm=TRUE),
    
    sw_sd   = sd(sw, na.rm=TRUE),
    
    sw_min  = min(sw, na.rm=TRUE),
    
    sw_max  = max(sw, na.rm=TRUE),
    
    swt_mean = mean(sw_trunc, na.rm=TRUE),
    
    swt_sd   = sd(sw_trunc, na.rm=TRUE),
    
    swt_min  = min(sw_trunc, na.rm=TRUE),
    
    swt_max  = max(sw_trunc, na.rm=TRUE)
    
  )]
  
  
  
  list(td=td, ps_model_d=ps_model_d, ps_model_n=ps_model_n, w_info=w_info)
  
}



fit_msm_one_imp <- function(mids_obj, i, denom_vars, msm_rhs, trunc_q=c(0.01,0.99)){
  
  dt <- as.data.table(complete(mids_obj, i))
  
  dt <- fix_exposure_na(standardise_types(dt))
  
  td <- make_td_tmerge(dt, vars_base=denom_vars)
  
  
  
  sw_out <- compute_sw_one_imp(td, denom_vars=denom_vars, numer_time_only=TRUE, trunc_q=trunc_q)
  
  if (is.null(sw_out)) return(NULL)
  
  
  
  td_w <- sw_out$td
  
  
  
  # Outcome model: weighted TD-Cox, stratify by NYHA if present
  
  if ("NYHA" %in% names(td_w)) {
    
    f <- as.formula(paste("Surv(tstart, tstop, death) ~", msm_rhs, "+ strata(NYHA)"))
    
  } else {
    
    f <- as.formula(paste("Surv(tstart, tstop, death) ~", msm_rhs))
    
  }
  
  
  
  fit <- try(
    
    coxph(
      
      f,
      
      data    = td_w,
      
      weights = sw_trunc,
      
      robust  = TRUE,
      
      cluster = ID,
      
      ties    = "efron"
      
    ),
    
    silent=TRUE
    
  )
  
  
  
  if (inherits(fit, "try-error")) return(NULL)
  
  
  
  list(fit=fit, sw_out=sw_out)
  
}



pool_fis_coef <- function(fit_objs){
  
  q <- c(); u <- c()
  
  for (x in fit_objs){
    
    fit <- x$fit
    
    if (!("FIS_td" %in% names(coef(fit)))) next
    
    q <- c(q, unname(coef(fit)["FIS_td"]))
    
    u <- c(u, unname(vcov(fit)["FIS_td","FIS_td"]))
    
  }
  
  if (length(q) < 2) return(list(ok=FALSE, est=NA, se=NA, p=NA, n_ok=length(q)))
  
  pooled <- pool_1coef(q, u)
  
  list(ok=TRUE, est=pooled$est, se=pooled$se, p=pooled$p, n_ok=length(q))
  
}



# -------------------------------------------------------------------------

# 7) Model definitions

# -------------------------------------------------------------------------

ps_denom_vars <- primary_covars



# PURE MSM (causal): Surv ~ FIS_td (weighted) + strata(NYHA)

msm_rhs_pure <- "FIS_td"



# Doubly robust: Surv ~ FIS_td + baseline covars (exclude NYHA because stratified)

adj_covs <- setdiff(ps_denom_vars, "NYHA")

msm_rhs_adj <- paste(c("FIS_td", adj_covs), collapse=" + ")



TRUNC_Q <- c(0.01, 0.99)



msg("Denominator covars: %s", paste(ps_denom_vars, collapse=", "))

msg("Truncation: %.1f%% / %.1f%%\n", TRUNC_Q[1]*100, TRUNC_Q[2]*100)



# -------------------------------------------------------------------------

# 8) TRAIN: fit across imputations

# -------------------------------------------------------------------------

msg("=== TRAIN MSM-IPTW: PURE ===")

fit_pure <- list()

wdiag_list <- list()

fail_pure <- integer(0)



for (i in 1:m_train){
  
  if (i %% 5 == 0 || i == 1) msg("  Train imputation %d/%d", i, m_train)
  
  
  
  out <- fit_msm_one_imp(
    
    mids_obj = t1_train,
    
    i        = i,
    
    denom_vars = ps_denom_vars,
    
    msm_rhs  = msm_rhs_pure,
    
    trunc_q  = TRUNC_Q
    
  )
  
  if (is.null(out)) {
    
    fail_pure <- c(fail_pure, i)
    
    next
    
  }
  
  fit_pure[[length(fit_pure)+1]] <- out
  
  wdiag_list[[length(wdiag_list)+1]] <- cbind(imputation=i, out$sw_out$w_info)
  
}



pooled_pure <- pool_fis_coef(fit_pure)

if (!pooled_pure$ok) {
  
  msg("\nERROR: Only %d successful PURE imputations. Failed imputations: %s",
      
      pooled_pure$n_ok, paste(fail_pure, collapse=", "))
  
  stop("Fewer than 2 successful imputations; cannot pool. (Try reducing covariates or check data issues.)")
  
}



msg("\nTRAIN PURE pooled over %d successful imputations:", pooled_pure$n_ok)

msg("  HR = %.3f", exp(pooled_pure$est))

msg("  95%% CI = (%.3f, %.3f)", exp(pooled_pure$est - 1.96*pooled_pure$se), exp(pooled_pure$est + 1.96*pooled_pure$se))

msg("  p = %.4f\n", pooled_pure$p)



msg("=== TRAIN MSM-IPTW: ADJUSTED (doubly robust) ===")

fit_adj <- list()

fail_adj <- integer(0)



for (i in 1:m_train){
  
  if (i %% 5 == 0 || i == 1) msg("  Train imputation %d/%d", i, m_train)
  
  
  
  out <- fit_msm_one_imp(
    
    mids_obj = t1_train,
    
    i        = i,
    
    denom_vars = ps_denom_vars,
    
    msm_rhs  = msm_rhs_adj,
    
    trunc_q  = TRUNC_Q
    
  )
  
  if (is.null(out)) {
    
    fail_adj <- c(fail_adj, i)
    
    next
    
  }
  
  fit_adj[[length(fit_adj)+1]] <- out
  
}



pooled_adj <- pool_fis_coef(fit_adj)

if (!pooled_adj$ok) {
  
  msg("\nERROR: Only %d successful ADJ imputations. Failed imputations: %s",
      
      pooled_adj$n_ok, paste(fail_adj, collapse=", "))
  
  stop("Fewer than 2 successful imputations; cannot pool adjusted model.")
  
}



msg("\nTRAIN ADJ pooled over %d successful imputations:", pooled_adj$n_ok)

msg("  HR = %.3f", exp(pooled_adj$est))

msg("  95%% CI = (%.3f, %.3f)", exp(pooled_adj$est - 1.96*pooled_adj$se), exp(pooled_adj$est + 1.96*pooled_adj$se))

msg("  p = %.4f\n", pooled_adj$p)



# Comparison table (in-memory)

comp <- data.table(
  
  model = c("SENS3 MSM-IPTW pure (weighted + strata(NYHA))",
            
            "SENS3 MSM-IPTW adjusted (weighted + covars + strata(NYHA))"),
  
  HR_FIS = c(exp(pooled_pure$est), exp(pooled_adj$est)),
  
  CI_lower = c(exp(pooled_pure$est - 1.96*pooled_pure$se),
               
               exp(pooled_adj$est - 1.96*pooled_adj$se)),
  
  CI_upper = c(exp(pooled_pure$est + 1.96*pooled_pure$se),
               
               exp(pooled_adj$est + 1.96*pooled_adj$se)),
  
  p_FIS  = c(pooled_pure$p, pooled_adj$p)
  
)



# Weight diagnostics (console)

wdiag_dt <- rbindlist(wdiag_list, fill=TRUE)

wdiag_summary <- wdiag_dt[, .(
  
  Mean_SW  = mean(sw_mean, na.rm=TRUE),
  
  SD_SW    = mean(sw_sd, na.rm=TRUE),
  
  Min_SW   = min(sw_min, na.rm=TRUE),
  
  Max_SW   = max(sw_max, na.rm=TRUE),
  
  Mean_SWT = mean(swt_mean, na.rm=TRUE),
  
  SD_SWT   = mean(swt_sd, na.rm=TRUE),
  
  Min_SWT  = min(swt_min, na.rm=TRUE),
  
  Max_SWT  = max(swt_max, na.rm=TRUE)
  
)]

msg("Weight diagnostics (train, across successful imputations):")

print(wdiag_summary)



# -------------------------------------------------------------------------

# 9) Balance diagnostics (TRAIN imp1) — FIXED + shows unweighted & weighted

# -------------------------------------------------------------------------

msg("\n=== Balance diagnostics (TRAIN imp1) ===")



dt1 <- fix_exposure_na(standardise_types(as.data.table(complete(t1_train, 1))))

td1 <- make_td_tmerge(dt1, ps_denom_vars)

sw1 <- compute_sw_one_imp(td1, denom_vars=ps_denom_vars, numer_time_only=TRUE, trunc_q=TRUNC_Q)



if (is.null(sw1)) stop("Failed to compute weights on TRAIN imp1 (PS model did not fit).")



# baseline one row per ID

base_cols <- c("ID", ps_denom_vars, "Status_FIS")

base_dt <- dt1[, ..base_cols]

base_dt <- base_dt[order(ID)][, .SD[1], by = ID]

base_dt[, ever_FIS := as.integer(Status_FIS == 1L)]



# attach weights at ID-level (weights are in td1, but sw1$td has sw_trunc per row; take 1 per ID)

w_id <- unique(sw1$td[, .(ID, sw_trunc)])

base_dt <- merge(base_dt, w_id, by="ID", all.x=TRUE)

base_dt[is.na(sw_trunc), sw_trunc := 1]



# prepare balance vars:

# - cobalt handles factors fine; keep NYHA factor

balance_vars <- ps_denom_vars



X_unw <- as.data.frame(base_dt[, ..balance_vars])

treat <- base_dt$ever_FIS

wts   <- base_dt$sw_trunc



# Unweighted balance

bal_unw <- cobalt::bal.tab(
  
  x = X_unw,
  
  treat = treat,
  
  s.d.denom = "pooled",
  
  abs = TRUE
  
)



# Weighted balance

bal_w <- cobalt::bal.tab(
  
  x = X_unw,
  
  treat = treat,
  
  weights = wts,
  
  method = "weighting",
  
  s.d.denom = "pooled",
  
  abs = TRUE
  
)



# Build a clean summary table

bal_sum <- data.table(
  
  covariate = rownames(bal_w$Balance),
  
  SMD_unweighted = bal_unw$Balance[rownames(bal_w$Balance), "Diff.Un"],
  
  SMD_weighted   = bal_w$Balance[, "Diff.Adj"]
  
)



msg("Balance summary (TRAIN imp1):")

msg("  Unweighted: %d/%d covariates with |SMD| > 0.1",
    
    sum(abs(bal_sum$SMD_unweighted) > 0.1, na.rm=TRUE), nrow(bal_sum))

msg("  Weighted:   %d/%d covariates with |SMD| > 0.1",
    
    sum(abs(bal_sum$SMD_weighted) > 0.1, na.rm=TRUE), nrow(bal_sum))



msg("\nTop covariates by |SMD| AFTER weighting:")

print(bal_sum[order(-abs(SMD_weighted))][1:min(15, .N)], row.names = FALSE)

# -------------------------------------------------------------------------

# 9b) Pretty balance table (gt) — keep ALL variables, but label NYHA separately

# -------------------------------------------------------------------------


bal_tbl <- copy(bal_sum)

bal_tbl[, abs_smd := abs(SMD_weighted)]

bal_tbl[, Balance := ifelse(abs_smd < 0.10, "Yes", "No")]



# Label NYHA rows (shown for transparency, but controlled via strata in outcome)

bal_tbl[, Group := ifelse(grepl("^NYHA", covariate),
                          
                          "Covariate controlled via stratification",
                          
                          "Covariates included in PS model")]



# Optional: nicer variable labels (edit as needed)

label_map <- c(
  
  Age = "Age (years)",
  
  LVEF = "LVEF (%)",
  
  eGFR = "eGFR (mL/min/1.73m²)",
  
  bin_diabetes_2 = "Diabetes",
  
  bin_af_atrial_flutter_2 = "Atrial fibrillation/flutter",
  
  bin_beta_blockers_2 = "Beta-blocker use",
  
  bin_sex_male = "Sex",
  
  bin_stroke_tia_2 = "Prior stroke/TIA",
  
  NYHA_I = "NYHA class I",
  
  NYHA_II = "NYHA class II",
  
  NYHA_III = "NYHA class III",
  
  NYHA_IV = "NYHA class IV"
  
)



bal_tbl[, Covariate := ifelse(covariate %in% names(label_map),
                              
                              unname(label_map[covariate]),
                              
                              covariate)]



bal_tbl <- bal_tbl[, .(
  
  Group, Covariate, SMD_weighted, Balance
  
)][order(Group, -abs(SMD_weighted))]



gt_balance <- bal_tbl |>
  
  gt(groupname_col = "Group") |>
  
  fmt_number(columns = SMD_weighted, decimals = 2) |>
  
  cols_label(
    
    Covariate = "Covariate",
    
    SMD_weighted = "SMD After Weighting",
    
    Balance = "Balance Achieved"
    
  ) |>
  
  tab_header(
    
    title = "Covariate Balance After IPTW (Training Dataset)"
    
  ) |>
  
  tab_source_note(
    
    source_note = "Balance defined as |SMD| < 0.10. NYHA is controlled via stratification in the outcome model and shown here for transparency."
    
  )



# Save as PDF (simple, robust)

gtsave(gt_balance, filename = file.path(OUTDIR, "SENS3_balance_table.pdf"))

# Save balance table as HTML (robust, no browser dependency)

gtsave(
  
  gt_balance,
  
  filename = file.path(OUTDIR, "SENS3_balance_table.html")
  
)

# ESS (approx)

ess <- sum(wts)^2 / sum(wts^2)

msg("\nEffective sample size (TRAIN imp1): %.1f (%.1f%% of %d)\n",
    
    ess, 100*ess/nrow(base_dt), nrow(base_dt))



# Love plot object (ggplot)

p_love <- love.plot(
  
  bal_w,
  
  threshold = 0.10,
  
  abs = TRUE,
  
  var.order = "unadjusted",
  
  title = "SENS3 MSM-IPTW: Covariate balance (TRAIN imp1)"
  
)



ggsave(file.path(OUTDIR, "SENS3_loveplot_imp1.pdf"), p_love, width=10, height=7)



# Weight histogram (ggplot)

p_w <- ggplot(base_dt, aes(x = sw_trunc)) +
  
  geom_histogram(bins = 50) +
  
  labs(title="SENS3: Truncated stabilized weights (TRAIN imp1)",
       
       x="sw_trunc", y="Count") +
  
  theme_minimal()



ggsave(file.path(OUTDIR, "SENS3_weights_hist_imp1.png"), p_w, width=8, height=5)



# -------------------------------------------------------------------------

# 10) PS MODEL CALIBRATION (TEST) — using training denominator PS model

# -------------------------------------------------------------------------

msg("=== PS Model Calibration (TEST) ===")



# use the PS denominator model from TRAIN imp1 (most reproducible)

ps_model_d_train1 <- sw1$ps_model_d

if (is.null(ps_model_d_train1)) {
  
  msg("WARNING: No train1 PS model found; skipping PS calibration.")
  
  do_calib <- FALSE
  
} else {
  
  do_calib <- TRUE
  
}



ps_calib_pooled <- NULL

calib_slope <- NA_real_

calib_intercept <- NA_real_

p_calib <- NULL



if (do_calib) {
  
  ps_calib_list <- list()
  
  
  
  for (i in 1:m_test) {
    
    dt_te <- fix_exposure_na(standardise_types(as.data.table(complete(t1_test, i))))
    
    td_te <- make_td_tmerge(dt_te, ps_denom_vars)
    
    ps_te <- make_initiation_long(td_te)
    
    if (nrow(ps_te) == 0) next
    
    
    
    # predict initiation hazard under denom PS model
    
    ps_te[, ps_pred := tryCatch(
      
      predict(ps_model_d_train1, newdata = ps_te, type = "response"),
      
      error = function(e) rep(NA_real_, .N)
      
    )]
    
    ps_te <- ps_te[!is.na(ps_pred)]
    
    if (nrow(ps_te) < 50) next
    
    
    
    breaks <- quantile(ps_te$ps_pred, probs = seq(0, 1, 0.1), na.rm = TRUE)
    
    if (length(unique(breaks)) < 11) next
    
    
    
    ps_te[, ps_decile := cut(ps_pred, breaks = breaks, include.lowest = TRUE, labels = FALSE)]
    
    
    
    calib <- ps_te[!is.na(ps_decile), .(
      
      predicted_ps  = mean(ps_pred, na.rm = TRUE),
      
      observed_rate = mean(switch,  na.rm = TRUE),
      
      n = .N
      
    ), by = ps_decile]
    
    
    
    calib[, imputation := i]
    
    ps_calib_list[[length(ps_calib_list) + 1]] <- calib
    
  }
  
  
  
  if (length(ps_calib_list) > 0) {
    
    ps_calib_all <- rbindlist(ps_calib_list, fill = TRUE)
    
    
    
    ps_calib_pooled <- ps_calib_all[, .(
      
      predicted_ps  = mean(predicted_ps,  na.rm = TRUE),
      
      observed_rate = mean(observed_rate, na.rm = TRUE),
      
      n = mean(n, na.rm = TRUE)
      
    ), by = ps_decile][order(ps_decile)]
    
    
    
    calib_lm <- lm(observed_rate ~ predicted_ps, data = ps_calib_pooled)
    
    calib_slope <- unname(coef(calib_lm)[2])
    
    calib_intercept <- unname(coef(calib_lm)[1])
    
    
    
    msg("PS Calibration (TEST):")
    
    msg("  Slope: %.3f (ideal=1.0)", calib_slope)
    
    msg("  Intercept: %.3f (ideal=0.0)\n", calib_intercept)
    
    
    
    p_calib <- ggplot(ps_calib_pooled, aes(x = predicted_ps, y = observed_rate)) +
      
      geom_point(aes(size = n), alpha = 0.7) +
      
      geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      
      geom_smooth(method = "lm", se = FALSE, linetype = "dotted") +
      
      labs(
        
        title = "Propensity Score Calibration (TEST)",
        
        subtitle = sprintf("Slope = %.3f, Intercept = %.3f", calib_slope, calib_intercept),
        
        x = "Mean predicted initiation probability (deciles)",
        
        y = "Observed initiation rate (deciles)",
        
        size = "N (person-intervals)"
        
      ) +
      
      theme_minimal()
    
    
    
    ggsave(file.path(OUTDIR, "SENS3_ps_calibration_test.png"), p_calib, width = 8, height = 6)
    
  } else {
    
    msg("WARNING: No valid calibration data generated for TEST.\n")
    
  }
  
}



# -------------------------------------------------------------------------

# 11) TEST discrimination (C-index) using TRAIN PURE fits (cycle)

# -------------------------------------------------------------------------

msg("=== TEST C-index for Outcome Model (discrimination) ===")



cidx <- rep(NA_real_, m_test)

n_fit_pure <- length(fit_pure)



for (i in 1:m_test){
  
  k <- ((i - 1) %% n_fit_pure) + 1
  
  fit_tr <- fit_pure[[k]]$fit
  
  
  
  dt_te <- fix_exposure_na(standardise_types(as.data.table(complete(t1_test, i))))
  
  td_te <- make_td_tmerge(dt_te, ps_denom_vars)
  
  
  
  lp <- tryCatch(predict(fit_tr, newdata=td_te, type="lp"), error=function(e) rep(NA_real_, nrow(td_te)))
  
  ok <- !is.na(lp)
  
  if (sum(ok) < 10) next
  
  
  
  cidx[i] <- survConcordance(Surv(td_te$tstart[ok], td_te$tstop[ok], td_te$death[ok]) ~ lp[ok])$concordance
  
}



msg("TEST Outcome Model C-index: mean=%.4f sd=%.4f min=%.4f max=%.4f\n",
    
    mean(cidx, na.rm=TRUE), sd(cidx, na.rm=TRUE),
    
    min(cidx, na.rm=TRUE), max(cidx, na.rm=TRUE))



p_cidx <- ggplot(data.table(cindex=cidx), aes(x=cindex)) +
  
  geom_histogram(bins=20) +
  
  labs(title="Outcome model C-index across test imputations", x="C-index", y="Count") +
  
  theme_minimal()



# -------------------------------------------------------------------------

# 12) PDF REPORT (all outputs printed)

# -------------------------------------------------------------------------

report_path <- file.path(OUTDIR, "SENS3_COMPLETE_REPORT_PRIMARYCOVARS.pdf")

msg("=== Creating PDF Report: %s ===", report_path)



pdf(report_path, width = 11, height = 8.5)



# Page 1: Title + key results

par(mar=c(0,0,0,0))

plot.new()

text(0.5, 0.92, "SENSITIVITY ANALYSIS 3: MSM-IPTW RESULTS", cex=2, font=2)

text(0.5, 0.84, "Time-Dependent Propensity Score (Initiation) + Stabilized IPTW", cex=1.3)

text(0.5, 0.78, paste("Generated:", Sys.time()), cex=1)

text(0.5, 0.72, sprintf("Train imputations: %d | Test imputations: %d", m_train, m_test), cex=1)



text(0.5, 0.60, "KEY RESULTS", cex=1.5, font=2)

text(0.08, 0.52, sprintf("Pure MSM (weighted + strata(NYHA)): HR = %.3f, p = %.4f",
                         
                         exp(pooled_pure$est), pooled_pure$p), adj=0, cex=1.2)

text(0.08, 0.46, sprintf("95%% CI: (%.3f – %.3f)",
                         
                         exp(pooled_pure$est - 1.96*pooled_pure$se),
                         
                         exp(pooled_pure$est + 1.96*pooled_pure$se)), adj=0, cex=1.2)



text(0.08, 0.36, sprintf("Adjusted MSM (weighted + covars + strata(NYHA)): HR = %.3f, p = %.4f",
                         
                         exp(pooled_adj$est), pooled_adj$p), adj=0, cex=1.2)

text(0.08, 0.30, sprintf("95%% CI: (%.3f – %.3f)",
                         
                         exp(pooled_adj$est - 1.96*pooled_adj$se),
                         
                         exp(pooled_adj$est + 1.96*pooled_adj$se)), adj=0, cex=1.2)



# Page 2: Model comparison table

plot.new()

text(0.5, 0.95, "MODEL COMPARISON", cex=1.8, font=2)



comp_display <- comp[, .(
  
  Model = model,
  
  HR    = sprintf("%.3f", HR_FIS),
  
  `95% CI` = sprintf("(%.3f – %.3f)", CI_lower, CI_upper),
  
  `P-value`= sprintf("%.4f", p_FIS)
  
)]



y_start <- 0.82

y_step  <- 0.12



text(0.05, y_start, "Model", font=2, adj=0, cex=1.05)

text(0.58, y_start, "HR", font=2, adj=0, cex=1.05)

text(0.70, y_start, "95% CI", font=2, adj=0, cex=1.05)

text(0.88, y_start, "P-value", font=2, adj=0, cex=1.05)

segments(0.05, y_start-0.03, 0.95, y_start-0.03, lwd=2)



for (r in 1:nrow(comp_display)){
  
  y <- y_start - r*y_step
  
  text(0.05, y, comp_display$Model[r], adj=0, cex=0.95)
  
  text(0.58, y, comp_display$HR[r], adj=0, cex=0.95)
  
  text(0.70, y, comp_display$`95% CI`[r], adj=0, cex=0.95)
  
  text(0.88, y, comp_display$`P-value`[r], adj=0, cex=0.95)
  
}



# Page 3: Weight diagnostics summary

plot.new()

text(0.5, 0.95, "WEIGHT DIAGNOSTICS (TRAIN)", cex=1.8, font=2)



text(0.10, 0.82, "Stabilized Weights (before truncation)", font=2, adj=0, cex=1.2)

text(0.12, 0.74, sprintf("Mean: %.3f", wdiag_summary$Mean_SW), adj=0, cex=1.1)

text(0.55, 0.74, sprintf("SD: %.3f", wdiag_summary$SD_SW), adj=0, cex=1.1)

text(0.12, 0.67, sprintf("Min: %.3f", wdiag_summary$Min_SW), adj=0, cex=1.1)

text(0.55, 0.67, sprintf("Max: %.3f", wdiag_summary$Max_SW), adj=0, cex=1.1)



text(0.10, 0.52, sprintf("Truncated Weights (%.0f–%.0f percentile)", TRUNC_Q[1]*100, TRUNC_Q[2]*100),
     
     font=2, adj=0, cex=1.2)

text(0.12, 0.44, sprintf("Mean: %.3f", wdiag_summary$Mean_SWT), adj=0, cex=1.1)

text(0.55, 0.44, sprintf("SD: %.3f", wdiag_summary$SD_SWT), adj=0, cex=1.1)

text(0.12, 0.37, sprintf("Min: %.3f", wdiag_summary$Min_SWT), adj=0, cex=1.1)

text(0.55, 0.37, sprintf("Max: %.3f", wdiag_summary$Max_SWT), adj=0, cex=1.1)



text(0.10, 0.22, sprintf("Effective Sample Size (TRAIN imp1): %.1f (%.1f%%)",
                         
                         ess, 100*ess/nrow(base_dt)),
     
     adj=0, cex=1.2, font=2)



# Page 4: Love plot

print(p_love)



# Page 5: Weight histogram

print(p_w)



# Page 6: PS Calibration plot or note

if (!is.null(p_calib)) {
  
  print(p_calib)
  
} else {
  
  plot.new()
  
  text(0.5, 0.6, "PS Calibration plot not available (insufficient decile data).", cex=1.3)
  
}



# Page 7: C-index distribution plot

print(p_cidx)



# Page 8: Methods summary

plot.new()

text(0.5, 0.95, "METHODS SUMMARY", cex=1.8, font=2)



y <- 0.86

text(0.08, y, "Time-dependent propensity score:", font=2, adj=0, cex=1.1); y <- y-0.06

text(0.10, y, "• Pooled logistic regression for initiation (switch 0→1) among not-yet-exposed", adj=0, cex=1.0); y <- y-0.05

text(0.10, y, paste0("• Denominator covariates: ", paste(ps_denom_vars, collapse=", "), " + time"), adj=0, cex=1.0); y <- y-0.05

text(0.10, y, "• Numerator: time-only stabilization", adj=0, cex=1.0); y <- y-0.08



text(0.08, y, "Weights:", font=2, adj=0, cex=1.1); y <- y-0.06

text(0.10, y, sprintf("• Stabilized IPTW, truncated at %.0fth and %.0fth percentiles",
                      
                      TRUNC_Q[1]*100, TRUNC_Q[2]*100), adj=0, cex=1.0); y <- y-0.08



text(0.08, y, "Outcome model:", font=2, adj=0, cex=1.1); y <- y-0.06

text(0.10, y, "• Time-dependent Cox model weighted by IPTW", adj=0, cex=1.0); y <- y-0.05

text(0.10, y, "• Stratified by NYHA", adj=0, cex=1.0); y <- y-0.05

text(0.10, y, "• Robust variance + clustering by patient ID", adj=0, cex=1.0); y <- y-0.08



text(0.08, y, "Diagnostics:", font=2, adj=0, cex=1.1); y <- y-0.06

text(0.10, y, "• Covariate balance via standardized mean differences (unweighted vs weighted)", adj=0, cex=1.0); y <- y-0.05

text(0.10, y, "• Outcome discrimination via C-index on test imputations", adj=0, cex=1.0)



dev.off()



msg("PDF report saved: %s", report_path)

msg("Also saved: loveplot PDF + weights hist PNG + PS calibration PNG (if available).")

msg("====================================")

msg("SENS3 COMPLETE. Outputs in: %s", OUTDIR)

msg("====================================")

