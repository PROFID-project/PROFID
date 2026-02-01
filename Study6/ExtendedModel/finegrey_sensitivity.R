## finegrey_sensitivity_90mo.R
## Fine–Gray subdistribution hazards for SCD vs non-cardiac death
## Extended model sensitivity analysis, restricted to 90-month follow-up

library(survival)
library(mice)
library(dplyr)
library(splines)
library(cmprsk)

setwd("T:/Dokumente/PROFID/Study6")

## 1) Load extended mids object
imp_ext <- readRDS("mice_imputed_data_extended.RDS")
m <- imp_ext$m

## 2) Follow-up horizon (months)
horizon <- 90

## 3) Covariate set (match UPDATED extended Cox model, excluding BMI)
## NOTE: COPD_cat and Cancer removed to stay consistent with updated model.
vars_extended <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute", "Cholesterol", "HDL", "LDL", "Triglycerides",
  "Stroke_TIA", "ICD_status"
)

## 4) Define BMI spline knots (same scheme as main analysis)
q <- quantile(
  imp_ext$data$BMI,
  probs = c(.05,.10,.35,.65,.90,.95),
  na.rm = TRUE
)
BND <- as.numeric(q[c(1,6)])          # 5th & 95th percentiles
K4  <- as.numeric(q[c(2,3,4,5)])      # 10th, 35th, 65th, 90th percentiles

## 5) Apply 90-month horizon for competing risks
apply_horizon_fg <- function(dat) {
  # Assumed coding for competing risks:
  #   0 = censored
  #   1 = SCD (event of interest)
  #   2 = non-cardiac death (competing event)
  #
  # If your dataset uses different coding, recode Status_fg here.
  dat$Status_fg <- dat$Status
  
  dat$ftime_h   <- pmin(dat$Survival_time, horizon)
  dat$fstatus_h <- ifelse(dat$Survival_time > horizon, 0L, dat$Status_fg)
  dat$fstatus_h <- as.integer(dat$fstatus_h)
  
  dat
}

## 6) Helper: completed data for imputation i + horizon vars
make_dat_i <- function(i) {
  dat_i <- complete(imp_ext, i)
  dat_i <- apply_horizon_fg(dat_i)
  dat_i
}

## 7) Fit Fine–Gray model in selected imputations
max_imps <- 3
imps_to_fit <- 1:min(max_imps, m)

coef_list <- vector("list", length(imps_to_fit))
vcov_list <- vector("list", length(imps_to_fit))

for (k in seq_along(imps_to_fit)) {
  i <- imps_to_fit[k]
  cat("Fitting Fine–Gray model, imputation", i, "of", length(imps_to_fit), "\n")
  flush.console()
  
  dat_i <- make_dat_i(i)
  
  # Build design matrix: BMI as RCS + all extended covariates
  mm_i <- model.matrix(
    ~ ns(BMI, knots = K4, Boundary.knots = BND) +
      Age + Sex +
      Diabetes + Hypertension + Smoking + MI_history +
      LVEF +
      eGFR + Haemoglobin +
      ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
      Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
      Stroke_TIA + ICD_status,
    data = dat_i
  )
  
  # Drop intercept
  mm_i <- mm_i[, -1, drop = FALSE]
  
  # Fine–Gray model (SCD vs non-cardiac death), restricted to horizon
  fg_i <- crr(
    ftime   = dat_i$ftime_h,
    fstatus = dat_i$fstatus_h,
    cov1    = mm_i
  )
  
  coef_list[[k]] <- fg_i$coef
  vcov_list[[k]] <- fg_i$var
}

## 8) Pool Fine–Gray coefficients across fitted imputations (Rubin-style)
pool_crr <- function(coef_list, vcov_list) {
  # Drop any NULLs defensively
  keep <- which(!vapply(coef_list, is.null, logical(1)) &
                  !vapply(vcov_list, is.null, logical(1)))
  coef_list <- coef_list[keep]
  vcov_list <- vcov_list[keep]
  
  m <- length(coef_list)
  if (m == 0) stop("No fitted imputations found in coef_list/vcov_list.")
  
  # Stack coefficients: p x m
  Q <- do.call(cbind, coef_list)
  rownames(Q) <- names(coef_list[[1]])
  
  # Average coefficients
  qbar <- rowMeans(Q)
  
  # Within-imputation variance (average vcov)
  ubar <- Reduce("+", vcov_list) / m
  
  # Between-imputation variance
  if (m > 1) {
    Q_centered <- sweep(Q, 1, qbar, "-")
    bmat <- (Q_centered %*% t(Q_centered)) / (m - 1)
  } else {
    bmat <- matrix(0, nrow = length(qbar), ncol = length(qbar),
                   dimnames = list(names(qbar), names(qbar)))
  }
  
  # Total variance
  Tmat <- ubar + (1 + 1/m) * bmat
  
  se    <- sqrt(diag(Tmat))
  zval  <- qbar / se
  pval  <- 2 * pnorm(-abs(zval))
  lower <- qbar - 1.96 * se
  upper <- qbar + 1.96 * se
  
  res <- data.frame(
    term      = names(qbar),
    estimate  = qbar,
    se        = se,
    z         = zval,
    p.value   = pval,
    conf.low  = lower,
    conf.high = upper,
    row.names = NULL
  )
  
  list(summary = res, coef = qbar, vcov = Tmat)
}

fg_pooled <- pool_crr(coef_list, vcov_list)
fg_summary <- fg_pooled$summary

## 9) Save pooled results
write.csv(
  fg_summary,
  "FineGray_SCD_vs_noncardiac_extended_pooled_results_90mo.csv",
  row.names = FALSE
)

saveRDS(
  list(
    coef_list = coef_list,
    vcov_list = vcov_list,
    pooled    = fg_pooled,
    horizon   = horizon,
    knots     = list(K4 = K4, BND = BND),
    imputations_fitted = imps_to_fit
  ),
  "FineGray_SCD_vs_noncardiac_extended_allfits_90mo.RDS"
)

cat("Fine–Gray competing risks analysis (90-month horizon) completed and saved.\n")
