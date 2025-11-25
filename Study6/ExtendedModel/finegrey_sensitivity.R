## sensitivity_FG_BMI_extended.R
## Fine–Gray subdistribution hazards for SCD vs non-cardiac death
## using extended imputed dataset (mice_imputed_data_extended.RDS)

library(survival)
library(mice)
library(dplyr)
library(splines)
library(cmprsk)

setwd("T:/Dokumente/PROFID/Study6")

## 1. Load extended mids object (already contains ICD_status, Cancer, COPD_cat, etc.)

imp_ext <- readRDS("mice_imputed_data_extended.RDS")
m       <- imp_ext$m

# Quick sanity check
imp_ext
names(imp_ext$data)[1:30]

## 2. Define covariate set (same as extended Cox model, excluding BMI)

vars_extended <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute", "Cholesterol", "HDL", "LDL", "Triglycerides",
  "COPD_cat", "Cancer", "Stroke_TIA", "ICD_status"
)

## 3. Define BMI spline knots (use same scheme as main analysis)

q   <- quantile(imp_ext$data$BMI,
                probs = c(.05,.10,.35,.65,.90,.95),
                na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])          # 5th & 95th
K4  <- as.numeric(q[c(2,3,4,5)])      # 10, 35, 65, 90

## 4. Helper to get completed data for imputation i

make_dat_i <- function(i) {
  dat_i <- complete(imp_ext, i)
  
  # Status coding assumption:
  #   0 = censored
  #   1 = SCD (event of interest)
  #   2 = non-cardiac death (competing event)
  # If that’s already true, we can use Status directly as fstatus.
  # If not, recode here.
  dat_i$Status_fg <- dat_i$Status
  
  dat_i
}

##  Fit Fine–Gray model in each imputation

coef_list <- vector("list", m)
vcov_list <- vector("list", m)

max_imps <- 3

for (i in 1:max_imps) {
  cat("Fitting Fine–Gray model, imputation", i, "of", m, "\n")
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
      COPD_cat + Cancer + Stroke_TIA + ICD_status,
    data = dat_i
  )
  
  # Drop intercept
  mm_i <- mm_i[, -1, drop = FALSE]
  
  # Fine–Gray subdistribution hazards model:
  # fstatus: 0=censor, 1= SCD (of interest), 2 = non-cardiac death (competing)
  fg_i <- crr(
    ftime   = dat_i$Survival_time,
    fstatus = dat_i$Status_fg,
    cov1    = mm_i
  )
  
  coef_list[[i]] <- fg_i$coef
  vcov_list[[i]] <- fg_i$var
}

##  Pool Fine–Gray coefficients across imputations (Rubin-style)

pool_crr <- function(coef_list, vcov_list) {
  m <- length(coef_list)
  
  # stack coefs: p x m
  Q <- do.call(cbind, coef_list)
  rownames(Q) <- names(coef_list[[1]])
  
  # average coefficients
  qbar <- rowMeans(Q)
  
  # within-imputation variance (average of vcov)
  ubar <- Reduce("+", vcov_list) / m
  
  # between-imputation variance
  Q_centered <- sweep(Q, 1, qbar, "-")
  bmat <- (Q_centered %*% t(Q_centered)) / (m - 1)
  
  # total variance
  Tmat <- ubar + (1 + 1/m) * bmat
  
  se   <- sqrt(diag(Tmat))
  zval <- qbar / se
  pval <- 2 * pnorm(-abs(zval))
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

## 7. Save pooled results

fg_summary <- fg_pooled$summary

write.csv(
  fg_summary,
  "FineGray_SCD_vs_noncardiac_extended_pooled_results.csv",
  row.names = FALSE
)

saveRDS(
  list(
    coef_list = coef_list,
    vcov_list = vcov_list,
    pooled    = fg_pooled
  ),
  "FineGray_SCD_vs_noncardiac_extended_allfits.RDS"
)

cat("Fine–Gray competing risks analysis completed and saved.\n")
