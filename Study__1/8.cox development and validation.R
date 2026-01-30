################################################################################
# COX PROPORTIONAL HAZARDS MODEL — DEVELOPMENT AND INTERNAL VALIDATION
# Time-dependent exposure analysis with inappropriate ICD therapy
# UPDATED: Proper tmerge implementation, correlation pruning, and validation
#
# Analysis overview:
#   - Outcome: All-cause mortality
#   - Exposure: Time-dependent inappropriate ICD therapy (FIS)
#   - Time scale: Days from study entry to death or censoring
#
# Model development (TRAIN set):
#   - Candidate covariates defined a priori based on clinical relevance
#   - Correlation-based pruning to reduce multicollinearity (|r| > 0.70)
#   - Full multivariable Cox model fitted with time-dependent exposure
#   - Backward selection applied to non-forced covariates
#   - Proportional hazards assumption assessed; NYHA stratified if violated
#   - interaction testing
#   - Prespecified subgroup analyses  

# Internal validation (TEST set):
#   - Patient-level train/test split to avoid information leakage
#   - Discrimination assessed using Harrell’s C-index and time-dependent AUC
#   - Calibration assessed using calibration slope at prespecified horizons
#
################################################################################


# Required packages ---------------------------------------------------------
pkgs <- c(
  "tidyverse","dplyr" ,"survival", "data.table", "mice", "naniar", "grid","gridExtra",
  "openxlsx", "readxl", "gt", "ggplot2", "timeROC","prodlim","riskRegression"
)

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))


OUTDIR <- "T:/Study_1/Primary analysis"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

t_train <- readRDS("T:/Imputed_data/mice_train_object3.rds")
t_test <- readRDS("T:/Imputed_data/mice_test_object3.rds") 
m_train <- t_train$m
m_test  <- t_test$m
names(complete(t_test, 1))
d1 = complete(t_train,1)
summary(d1$Time_death_days)

EPS <- 1e-7
msg <- function(...) cat(sprintf(...), "\n")


# -------------------------------------------------------------------------
# 1) Variable lists
# -------------------------------------------------------------------------
forced_vars <- c(
  "Age",
  "bin_sex_male",
  "LVEF",
  "NYHA",
  "bin_diabetes",
  "eGFR",
  "bin_beta_blockers",
  "bin_af_atrial_flutter"
)

candidate_vars <- c(
  "bin_hypertension","bin_smoking","Haemoglobin","SBP","BMI","HR",
   "QRS_log1p",
  "bin_hf","bin_stroke_tia","bin_av_block","bin_lbbb",
  "bin_ace_inhibitor","bin_ace_inhibitor_arb","bin_diuretics",
  "bin_anti_arrhythmic_iii","bin_anti_coagulant","bin_anti_platelet",
  "bin_lipid_lowering","bin_anti_diabetic_oral","bin_digitalis_glycosides",
  "bin_pci","HasMRI"
)

all_covars <- unique(c(forced_vars, candidate_vars))

# -------------------------------------------------------------------------
# 2) Sanity checks 
# -------------------------------------------------------------------------
dt_ref <- as.data.table(complete(t_train, 1))

is.ordered(dt_ref$NYHA)
getOption("contrasts")

# learn NYHA levels from train imp1 to keep stable everywhere
nyha_levels <- if ("NYHA" %in% names(dt_ref)) levels(factor(dt_ref$NYHA)) else NULL
#cvd_levels  <- c("1","2","3","4")

# -------------------------------------------------------------------------

# 3) Helpers (minimal, robust)

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
  
  if ("HasMRI" %in% names(dt)) dt[, HasMRI := as.integer(HasMRI)]
  
  
  
  # numeric covariates if present
  
  num_candidates <- c("Age","LVEF","eGFR","Haemoglobin","SBP","BMI","HR","CRP_log1p","QRS_log1p")
  
  for (v in intersect(names(dt), num_candidates)) dt[, (v) := as.numeric(get(v))]
  
  
  
  # factors
  
  #if (!is.null(nyha_levels) && "NYHA" %in% names(dt)) {
    
   # dt[, NYHA := factor(NYHA, levels = nyha_levels)]
    
  #}
  
  if (!is.null(nyha_levels) && "NYHA" %in% names(dt)) {
    
    dt[, NYHA := factor(as.character(NYHA), levels = nyha_levels, ordered = FALSE)]
    
    contrasts(dt$NYHA) <- contr.treatment(nlevels(dt$NYHA))
    
  }
  #if ("CVD_risk_region" %in% names(dt)) {
    
    #dt[, CVD_risk_region := factor(as.character(CVD_risk_region), levels = cvd_levels)]
       
  #}
  
  
  
  dt
  
}



# -------------------------------------------------------------------------
# 3a) Missing Data Handling for Time-Dependent Exposure Variables
# -------------------------------------------------------------------------
fix_exposure_na <- function(dt){
  dt <- as.data.table(dt)
  
  # A: both missing => assume not exposed
  dt[is.na(Status_FIS) & is.na(Time_FIS_days), Status_FIS := 0L]
  
  # B: time present but status missing => exposed
  dt[is.na(Status_FIS) & !is.na(Time_FIS_days), Status_FIS := 1L]
  
  # C: exposed but time missing => ambiguous => exclude from TD
  dt[, excl_fis_ambig := (Status_FIS == 1L & is.na(Time_FIS_days))]
  
  # if not exposed, time should be NA
  dt[Status_FIS == 0L, Time_FIS_days := NA_real_]
  
  # TD rule: exposure at/after end of follow-up => no exposed person-time => treat unexposed
  dt[!is.na(Time_FIS_days) & !is.na(Time_death_days) & Time_FIS_days >= Time_death_days,
     `:=`(Status_FIS = 0L, Time_FIS_days = NA_real_)]
  
  dt
}

# -------------------------------------------------------------------------
# 3b) Restructuring Data to Long (Start–Stop) Format Using tmerge
# -------------------------------------------------------------------------
make_td_tmerge <- function(dt, vars_base){
  dt <- as.data.table(dt)
  
  # require outcome
  dt <- dt[!is.na(Time_death_days) & !is.na(Status_death)]
  
  # exclude ambiguous exposure timing
  dt <- dt[excl_fis_ambig != TRUE]
  
  # enforce 1 row per ID (baseline row)
  dt <- dt[order(ID)][, .SD[1], by = ID]
  
  # ========== KEY FIX: PROPER data1 and data2 structure ==========
  
  # data1: ID + outcome time/status + baseline covariates
  base_cols <- c("ID", "Time_death_days", "Status_death", vars_base)
  data1 <- dt[, ..base_cols]
  
  # data2: SAME as data1 + BOTH Status_FIS and Time_FIS_days
  d2_cols <- c("ID", "Time_death_days", "Status_death", vars_base, 
               "Status_FIS", "Time_FIS_days")
  data2 <- dt[, ..d2_cols]
  
  # Create time-dependent dataset using tmerge
  td <- tmerge(
    data1 = data1,                                    # Baseline data
    data2 = data2,                                    # Baseline + exposure vars
    id    = ID,                                       # ID variable
    death = event(Time_death_days, Status_death == 1), # Event indicator
    FIS_td = tdc(Time_FIS_days)                       # Time-dependent covariate
  )
  
  # Clean up and format
  td <- as.data.table(td)
  td[, tstop := pmax(tstop, EPS)]        # Ensure positive tstop
  td <- td[tstop > tstart]                # Remove invalid intervals
  setorder(td, ID, tstart, tstop)         # Sort by ID and time
  
  return(td)
}

# -------------------------------------------------------------------------
# 3c) Enhanced diagnostic (train imp1) - show wide vs. long conversion
# -------------------------------------------------------------------------
msg("\n=== DIAGNOSTIC: Wide to Long Conversion (train imp1) ===")

dt0 <- fix_exposure_na(standardise_types(copy(dt_ref)))
td0 <- make_td_tmerge(dt0, all_covars)

msg("Wide format: %d patients", uniqueN(dt0$ID))
msg("Long format: %d rows, %d unique patients, %d deaths",
    nrow(td0), uniqueN(td0$ID), sum(td0$death))

# Example 1: Exposed patient (FIS occurred)
ex_exposed <- dt0[Status_FIS == 1L & !is.na(Time_FIS_days), ID][1]
if (!is.na(ex_exposed)) {
  msg("\n=== EXAMPLE 1: Exposed Patient (ID = %s) ===", ex_exposed)
  msg("WIDE FORMAT:")
  print(dt0[ID == ex_exposed, .(ID, Status_FIS, Time_FIS_days, Status_death, Time_death_days)])
  msg("\nLONG FORMAT (time-split intervals):")
  print(td0[ID == ex_exposed, .(ID, tstart, tstop, death, FIS_td)])
  msg("Interpretation: Patient split into 2 intervals - unexposed (FIS_td=0) before Time_FIS_days, exposed (FIS_td=1) after")
}

# Example 2: Unexposed patient (no FIS)
ex_unexposed <- dt0[Status_FIS == 0L, ID][1]
if (!is.na(ex_unexposed)) {
  msg("\n=== EXAMPLE 2: Unexposed Patient (ID = %s) ===", ex_unexposed)
  msg("WIDE FORMAT:")
  print(dt0[ID == ex_unexposed, .(ID, Status_FIS, Time_FIS_days, Status_death, Time_death_days)])
  msg("\nLONG FORMAT (single interval):")
  print(td0[ID == ex_unexposed, .(ID, tstart, tstop, death, FIS_td)])
  msg("Interpretation: Patient has 1 interval - never exposed (FIS_td=0) throughout follow-up")
}

# Summary statistics
fis_summary <- dt0[, .(
  n_total = .N,
  n_exposed = sum(Status_FIS == 1L, na.rm=TRUE),
  n_unexposed = sum(Status_FIS == 0L, na.rm=TRUE),
  n_ambiguous = sum(excl_fis_ambig, na.rm=TRUE)
)]
msg("\nExposure Summary:")
print(fis_summary)


# -------------------------------------------------------------------------
# 4) POOLING FUNCTIONS
# -------------------------------------------------------------------------

# Pool single coefficient (binary/continuous variables)
pool_single_coef <- function(coef_vector, variance_vector) {
  m <- length(coef_vector)
  coef_pooled <- mean(coef_vector)
  var_within <- mean(variance_vector)
  var_between <- var(coef_vector)
  var_total <- var_within + (1 + 1/m) * var_between
  se_pooled <- sqrt(var_total)
  z_stat <- coef_pooled / se_pooled
  p_value <- 2 * pnorm(abs(z_stat), lower.tail = FALSE)
  
  list(estimate = coef_pooled, se = se_pooled, p = p_value)
}

# Pool multiple coefficients (factor variables) using Wald test
pool_multiple_coefs <- function(coef_matrix, vcov_list) {
  m <- nrow(coef_matrix)
  k <- ncol(coef_matrix)
  coef_pooled <- colMeans(coef_matrix)
  vcov_within <- Reduce("+", vcov_list) / m
  vcov_between <- cov(coef_matrix)
  vcov_total <- vcov_within + (1 + 1/m) * vcov_between
  wald_stat <- as.numeric(t(coef_pooled) %*% solve(vcov_total) %*% coef_pooled)
  p_value <- pchisq(wald_stat, df = k, lower.tail = FALSE)
  
  list(wald = wald_stat, df = k, p = p_value)
}

# Smart pooling - automatically chooses method
pool_variable <- function(mids_obj, vars_base, var_name) {
  formula_rhs <- paste("FIS_td +", var_name)
  other_vars <- setdiff(vars_base, var_name)
  if (length(other_vars) > 0) {
    formula_rhs <- paste(formula_rhs, "+", paste(other_vars, collapse=" + "))
  }
  
  coef_list <- list()
  vcov_list <- list()
  
  for (i in 1:mids_obj$m) {
    dt <- as.data.table(complete(mids_obj, i))
    dt <- fix_exposure_na(standardise_types(dt))
    all_vars_needed <- unique(c(vars_base, var_name))
    td <- make_td_tmerge(dt, all_vars_needed)
    
    if (nrow(td) == 0) next
    
    f <- as.formula(paste("Surv(tstart, tstop, death) ~", formula_rhs))
    fit <- try(coxph(f, data=td, ties="efron"), silent=TRUE)
    
    if (inherits(fit, "try-error")) next
    
    all_terms <- names(coef(fit))
    var_terms <- if (var_name %in% all_terms) {
      var_name
    } else {
      grep(paste0("^", var_name), all_terms, value=TRUE)
    }
    
    if (length(var_terms) == 0) next
    
    coef_list[[length(coef_list) + 1]] <- coef(fit)[var_terms]
    vcov_list[[length(vcov_list) + 1]] <- vcov(fit)[var_terms, var_terms, drop=FALSE]
  }
  
  if (length(coef_list) < 2) {
    return(list(p = NA_real_, type = "failed"))
  }
  
  n_coefs <- length(coef_list[[1]])
  
  if (n_coefs == 1) {
    coef_vec <- sapply(coef_list, function(x) as.numeric(x))
    var_vec <- sapply(vcov_list, function(x) as.numeric(x))
    result <- pool_single_coef(coef_vec, var_vec)
    return(list(p = result$p, type = "single"))
  } else {
    coef_matrix <- do.call(rbind, lapply(coef_list, as.numeric))
    result <- pool_multiple_coefs(coef_matrix, vcov_list)
    return(list(p = result$p, df = result$df, type = "factor"))
  }
}

# Fit model and pool all coefficients
fit_and_pool_all <- function(mids_obj, vars_base, formula_rhs) {
  fit_list <- list()
  
  for (i in 1:mids_obj$m) {
    dt <- as.data.table(complete(mids_obj, i))
    dt <- fix_exposure_na(standardise_types(dt))
    td <- make_td_tmerge(dt, vars_base)
    
    if (nrow(td) == 0) next
    
    f <- as.formula(paste("Surv(tstart, tstop, death) ~", formula_rhs))
    fit <- try(coxph(f, data=td, ties="efron"), silent=TRUE)
    
    if (!inherits(fit, "try-error")) {
      fit_list[[length(fit_list) + 1]] <- fit
    }
  }
  
  if (length(fit_list) < 2) return(NULL)
  
  all_terms <- names(coef(fit_list[[1]]))
  
  results <- lapply(all_terms, function(term) {
    coef_vec <- sapply(fit_list, function(f) coef(f)[term])
    var_vec <- sapply(fit_list, function(f) vcov(f)[term, term])
    pooled <- pool_single_coef(coef_vec, var_vec)
    
    data.table(
      term = term,
      estimate = pooled$estimate,
      se = pooled$se,
      HR = exp(pooled$estimate),
      lower = exp(pooled$estimate - 1.96 * pooled$se),
      upper = exp(pooled$estimate + 1.96 * pooled$se),
      p.value = pooled$p
    )
  })
  
  rbindlist(results)
}


# 5) STEP 1: Screening 
# -------------------------------------------------------------------------

P_THRESH <- 0.05

# Crude FIS model
msg("Fitting crude FIS_td model...")
fis_crude <- fit_and_pool_all(t_train, character(0), "FIS_td")

if (is.null(fis_crude) || !("FIS_td" %in% fis_crude$term)) {
  stop("Crude FIS_td model failed!")
}

msg("\n*** CRUDE FIS_td RESULT ***")
msg("HR: %.3f (95%% CI: %.3f-%.3f)", 
    fis_crude[term=="FIS_td", HR],
    fis_crude[term=="FIS_td", lower],
    fis_crude[term=="FIS_td", upper])
msg("P-value: %.4f\n", fis_crude[term=="FIS_td", p.value])

# Save crude result
gtsave(fis_crude, file.path(OUTDIR, "00_crude_FIS_model.html"))
gtsave(
  gt::gt(fis_crude),
  filename = file.path(OUTDIR, "00_crude_FIS_model.html"),
  inline_css = TRUE
)

# Screen each covariate
msg("Screening %d covariates...", length(all_covars))

screen_results <- rbindlist(lapply(seq_along(all_covars), function(idx) {
  x <- all_covars[idx]
  
  if (idx %% 5 == 0) msg("  Progress: %d/%d", idx, length(all_covars))
  
  # Fit model: Surv ~ FIS_td + X
  pooled_model <- fit_and_pool_all(t_train, x, paste("FIS_td +", x))
  
  if (is.null(pooled_model)) {
    return(data.table(
      covariate = x,
      HR = NA_real_,
      CI_lower = NA_real_,
      CI_upper = NA_real_,
      p_value = NA_real_,
      HR_FIS_crude = fis_crude[term=="FIS_td", HR],
      HR_FIS_adj = NA_real_,
      confounding_pct = NA_real_,
      test_type = "failed",
      keep = FALSE
    ))
  }
  
  # Get FIS_td from adjusted model
  fis_row <- pooled_model[term == "FIS_td", ]
  
  # Get covariate results
  x_rows <- pooled_model[term != "FIS_td", ]
  
  # Calculate confounding
  hr_crude_log <- log(fis_crude[term=="FIS_td", HR])
  hr_adj_log <- log(fis_row$HR)
  confound_pct <- if (is.finite(hr_crude_log) && hr_crude_log != 0) {
    100 * (hr_adj_log - hr_crude_log) / abs(hr_crude_log)
  } else NA_real_
  
  if (nrow(x_rows) == 0) {
    x_hr <- NA_real_
    x_ci_low <- NA_real_
    x_ci_high <- NA_real_
    x_p <- NA_real_
    test_type <- "error"
  } else if (nrow(x_rows) == 1) {
    x_hr <- x_rows$HR
    x_ci_low <- x_rows$lower
    x_ci_high <- x_rows$upper
    x_p <- x_rows$p.value
    test_type <- "single"
  } else {
    result_wald <- pool_variable(t_train, x, x)
    x_hr <- NA_real_
    x_ci_low <- NA_real_
    x_ci_high <- NA_real_
    x_p <- result_wald$p
    test_type <- "factor"
  }
  
  data.table(
    covariate = x,
    HR = x_hr,
    CI_lower = x_ci_low,
    CI_upper = x_ci_high,
    p_value = x_p,
    HR_FIS_crude = fis_crude[term=="FIS_td", HR],
    HR_FIS_adj = fis_row$HR,
    confounding_pct = confound_pct,
    test_type = test_type,
    keep = !is.na(x_p) & x_p < P_THRESH
  )
}))

setorder(screen_results, p_value)
#fwrite(screen_results, file.path(OUTDIR, "01_screening_summary3.csv"))
gtsave(
  gt::gt(screen_results),
  filename = file.path(OUTDIR, "01_screening_summary.html"),
  inline_css = TRUE
)

# Save factor details
factor_vars <- screen_results[test_type == "factor", covariate]
if (length(factor_vars) > 0) {
  factor_details <- rbindlist(lapply(factor_vars, function(x) {
    pooled <- fit_and_pool_all(t_train, x, paste("FIS_td +", x))
    if (is.null(pooled)) return(NULL)
    x_rows <- pooled[term != "FIS_td", ]
    x_rows[, covariate := x]
    x_rows
  }))
  fwrite(factor_details, file.path(OUTDIR, "01_screening_factor_details1.csv"))
}

# Print summary
msg("\n=== Screening Summary (Threshold: p < %.2f) ===", P_THRESH)
summary_print <- screen_results[, .(
  Covariate = covariate,
  HR = ifelse(is.na(HR), "see_factor", sprintf("%.3f", HR)),
  CI_95 = ifelse(is.na(CI_lower), "-", sprintf("%.3f-%.3f", CI_lower, CI_upper)),
  P_value = sprintf("%.2e", p_value),
  Type = test_type,
  Keep = keep
)]
print(summary_print[1:min(20, nrow(summary_print)), ])
gtsave(
  gt::gt(summary_print),
  filename = file.path(OUTDIR, "01_screening_summary_formatted.html"),
  inline_css = TRUE
)
screened_vars <- screen_results[keep==TRUE, covariate]
msg("\nScreened-in (%d variables): %s", 
    length(screened_vars),
    paste(screened_vars, collapse=", "))

full_vars <- unique(c(forced_vars, screened_vars))
msg("Total for next step (forced + screened): %d variables", length(full_vars))


# -------------------------------------------------------------------------

# 5) STEP 2: Correlation pruning (|r| > 0.70) ===")

# -------------------------------------------------------------------------



dt_num <- standardise_types(as.data.table(complete(t_train, 1)))

num_vars <- intersect(full_vars, names(dt_num)[sapply(dt_num, is.numeric)])



remove_cor <- character()

if (length(num_vars) >= 2) {
  
  cor_mat <- cor(dt_num[, ..num_vars], use="pairwise.complete.obs")
  
  cor_mat[lower.tri(cor_mat, diag=TRUE)] <- NA
  
  
  
  high_cor <- which(abs(cor_mat) > 0.70, arr.ind=TRUE)
  
  
  
  if (nrow(high_cor) > 0) {
    
    pairs_df <- data.table(
      
      var1 = rownames(cor_mat)[high_cor[,1]],
      
      var2 = colnames(cor_mat)[high_cor[,2]],
      
      correlation = cor_mat[high_cor]
      
    )
    
    
    
    # Add absolute correlation column for sorting
    
    pairs_df[, abs_correlation := abs(correlation)]
    
    setorder(pairs_df, -abs_correlation)  # Sort by absolute correlation (descending)
    
    gtsave(
      gt::gt(pairs_df),
      filename = file.path(OUTDIR, "02_high_correlations.html"),
      inline_css = TRUE
    )
    
    fwrite(pairs_df, file.path(OUTDIR, "02_high_correlations.csv"))
    
    msg("Found %d high correlation pairs", nrow(pairs_df))
    
    
    
    for (i in 1:nrow(pairs_df)) {
      
      v1 <- pairs_df$var1[i]
      
      v2 <- pairs_df$var2[i]
      
      
      
      if (v1 %in% remove_cor || v2 %in% remove_cor) next
      
      
      
      if (v1 %in% forced_vars) {
        
        remove_cor <- c(remove_cor, v2)
        
        msg("  Keeping %s (forced), removing %s (r=%.3f)", v1, v2, pairs_df$correlation[i])
        
      } else if (v2 %in% forced_vars) {
        
        remove_cor <- c(remove_cor, v1)
        
        msg("  Keeping %s (forced), removing %s (r=%.3f)", v2, v1, pairs_df$correlation[i])
        
      } else {
        
        p1 <- screen_results[covariate==v1, p_value]
        
        p2 <- screen_results[covariate==v2, p_value]
        
        if (p1 > p2) {
          
          remove_cor <- c(remove_cor, v1)
          
          msg("  Removing %s (p=%.3f), keeping %s (p=%.3f)", v1, p1, v2, p2)
          
        } else {
          
          remove_cor <- c(remove_cor, v2)
          
          msg("  Removing %s (p=%.3f), keeping %s (p=%.3f)", v2, p2, v1, p1)
          
        }
        
      }
      
    }
    
  } else {
    
    msg("No high correlations found")
    
  }
  
}



model_vars <- setdiff(full_vars, unique(remove_cor))

msg("\nAfter correlation pruning: %d variables", length(model_vars))

msg("Variables: %s", paste(model_vars, collapse=", "))



# -------------------------------------------------------------------------
# 5) STEP 3: Full model
# -------------------------------------------------------------------------

rhs_full <- paste("FIS_td +", paste(model_vars, collapse=" + "))
full_pooled <- fit_and_pool_all(t_train, model_vars, rhs_full)

if (is.null(full_pooled)) stop("Full model failed!")

setorder(full_pooled, p.value)
fwrite(full_pooled, file.path(OUTDIR, "03_full_model.csv"))

msg("\n*** FIS_td in FULL MODEL ***")
msg("HR: %.3f (95%% CI: %.3f-%.3f)", 
    full_pooled[term=="FIS_td", HR],
    full_pooled[term=="FIS_td", lower],
    full_pooled[term=="FIS_td", upper])
msg("P-value: %.4f", full_pooled[term=="FIS_td", p.value])

msg("\nTop 5 predictors (by p-value):")
print(full_pooled[1:min(5, nrow(full_pooled)), .(term, HR, lower, upper, p.value)])


#-------------------------------------------------------------------------
# 5) STEP 4: Backward selection (p < 0.05)
# -------------------------------------------------------------------------

vars_current <- model_vars

repeat {
  removable <- setdiff(vars_current, forced_vars)
  if (length(removable) == 0) {
    msg("All remaining variables are forced - stopping")
    break
  }
  
  var_pvals <- sapply(removable, function(v) {
    result <- pool_variable(t_train, vars_current, v)
    result$p
  })
  
  max_p <- max(var_pvals, na.rm=TRUE)
  if (!is.finite(max_p) || max_p < 0.05) {
    msg("All p-values < 0.05 - stopping")
    break
  }
  
  drop_var <- names(which.max(var_pvals))
  msg("  Dropping %s (p=%.4f)", drop_var, max_p)
  vars_current <- setdiff(vars_current, drop_var)
}

final_vars <- vars_current
msg("\nFinal model: %d variables", length(final_vars))
msg("Variables: %s", paste(final_vars, collapse=", "))



# -------------------------------------------------------------------------
# 5) STEP 5: Final model results
# -------------------------------------------------------------------------

rhs_final <- paste("FIS_td +", paste(final_vars, collapse=" + "))
final_pooled <- fit_and_pool_all(t_train, final_vars, rhs_final)

setorder(final_pooled, p.value)
gtsave(
  gt::gt(final_pooled),
  filename = file.path(OUTDIR, "04_final_model1.html"),
  inline_css = TRUE
)
fwrite(final_pooled, file.path(OUTDIR, "04_final_model1.csv"))

msg("\n*** FINAL FIS_td RESULT ***")
msg("HR: %.3f (95%% CI: %.3f-%.3f)", 
    final_pooled[term=="FIS_td", HR],
    final_pooled[term=="FIS_td", lower],
    final_pooled[term=="FIS_td", upper])
msg("P-value: %.4f", final_pooled[term=="FIS_td", p.value])

msg("\nAll final model coefficients:")
print(final_pooled[, .(term, HR, lower, upper, p.value)])

# TEMPORARY row renaming (presentation-only) + save as HTML
# (does NOT modify final_pooled)

gtsave(
  gt::gt(final_pooled) |>
    gt::text_transform(
      locations = gt::cells_body(columns = term),
      fn = function(x) {
        dplyr::recode(
          x,
          "LVEF"                  = "LVEF",
          "Age"                   = "Age",
          "eGFR"                  = "eGFR",
          "bin_beta_blockers"     = "Beta-blockers",
          "bin_diabetes"          = "Diabetes",
          "bin_stroke_tia"        = "Stroke/TIA",
          "bin_af_atrial_flutter" = "AF/Atrial flutter",
          "bin_sex_male"          = "Sex",
          "FIS_td"                = "Inappropriate ICD therapy (time-dependent)",
          "NYHA2"                 = "NYHA II",
          "NYHA3"                 = "NYHA III",
          "NYHA4"                 = "NYHA IV",
          .default = x
        )
      }
    ) |>
    gt::cols_label(
      term     = "Covariate",
      estimate = "β (log HR)",
      se       = "SE",
      HR       = "HR",
      lower    = "95% CI (Lower)",
      upper    = "95% CI (Upper)",
      p.value  = "P-value"
    ) |>
    gt::fmt_number(columns = c(estimate, se, HR, lower, upper), decimals = 2) |>
    gt::fmt_scientific(columns = p.value, decimals = 2),
  filename = file.path(OUTDIR, "04_final_model1_pretty.html"),
  inline_css = TRUE
)


# -------------------------------------------------------------------------
# 5) STEP 6: Assumptions (train imp1)
# -------------------------------------------------------------------------

dt1 <- fix_exposure_na(standardise_types(as.data.table(complete(t_train, 1))))
td1 <- make_td_tmerge(dt1, final_vars)
fit1 <- coxph(as.formula(paste("Surv(tstart,tstop,death) ~", rhs_final)), 
              data=td1, ties="efron")

# PH assumption
ph <- cox.zph(fit1)
sink(file.path(OUTDIR, "05_ph_test.txt"))
print(ph)
sink()

pdf(file.path(OUTDIR, "05_ph_plots.pdf"), width=12, height=10)
plot(ph)
dev.off()

ph_violations <- which(ph$table[,"p"] < 0.05)
if (length(ph_violations) > 0) {
  msg("WARNING: PH assumption violated for: %s", 
      paste(rownames(ph$table)[ph_violations], collapse=", "))
} else {
  msg("PH assumption satisfied for all variables")
}

# Linearity
cont_vars <- c("Age","LVEF","eGFR","Haemoglobin","SBP","BMI","HR","CRP_log1p","QRS_log1p")
cont_in_model <- intersect(final_vars, cont_vars)

if (length(cont_in_model) > 0) {
  fit_null <- coxph(Surv(tstart,tstop,death) ~ 1, data=td1)
  mart <- residuals(fit_null, type="martingale")
  
  pdf(file.path(OUTDIR, "05_linearity_plots.pdf"), width=12, height=10)
  par(mfrow=c(ceiling(length(cont_in_model)/2), 2))
  for (v in cont_in_model) {
    x <- td1[[v]]
    ok <- is.finite(x) & is.finite(mart)
    plot(x[ok], mart[ok], pch=16, cex=0.5, xlab=v, ylab="Martingale residual",
         main=paste("Linearity check:", v))
    lines(lowess(x[ok], mart[ok]), lwd=2, col="red")
    abline(h=0, lty=2)
  }
  dev.off()
  msg("Linearity plots saved for %d continuous variables", length(cont_in_model))
}


# -------------------------------------------------------------------------

# 6) Cox Model with NYHA Stratification (PH Assumption)

# -------------------------------------------------------------------------



# Keep ORIGINAL final model objects (do not modify)

final_vars_orig   <- final_vars

rhs_final_orig    <- rhs_final

final_pooled_orig <- final_pooled



# Create STRATIFIED copy (PH-fix)

final_vars_strat <- final_vars_orig



# Remove NYHA from covariate list because it will be handled via strata(NYHA)

if ("NYHA" %in% final_vars_strat) {
  
  final_vars_strat <- setdiff(final_vars_strat, "NYHA")
  
} else {
  
  msg("NOTE: NYHA not in final_vars; stratification will still be added via strata(NYHA).")
  
}



# RHS formula for stratified model

rhs_final_strat <- paste(
  
  "FIS_td +",
  
  paste(final_vars_strat, collapse = " + "),
  
  "+ strata(NYHA)"
  
)



# IMPORTANT FIX:

# vars_base is used to BUILD the tmerge dataset; it MUST include NYHA because RHS uses strata(NYHA)

vars_for_td_strat <- unique(c(final_vars_strat, "NYHA"))



final_pooled_strat <- fit_and_pool_all(
  
  mids_obj    = t_train,
  
  vars_base   = vars_for_td_strat,
  
  formula_rhs = rhs_final_strat
  
)



if (is.null(final_pooled_strat)) stop("NYHA-stratified copy model failed!")



setorder(final_pooled_strat, p.value)

fwrite(final_pooled_strat, file.path(OUTDIR, "04_final_model_NYHA_stratified_COPY.csv"))

msg("\n=== NYHA-STRATIFIED COPY: ALL MODEL COEFFICIENTS ===")

print(final_pooled_strat[, .(term, HR, lower, upper, p.value)])

msg("\n*** NYHA-STRATIFIED COPY: FIS_td RESULT ***")

msg("HR: %.3f (95%% CI: %.3f-%.3f)",
    
    final_pooled_strat[term=="FIS_td", HR],
    
    final_pooled_strat[term=="FIS_td", lower],
    
    final_pooled_strat[term=="FIS_td", upper])

msg("P-value: %.4f", final_pooled_strat[term=="FIS_td", p.value])

gtsave(
  gt::gt(final_pooled_strat) |>
    gt::text_transform(
      locations = gt::cells_body(columns = term),
      fn = function(x) {
        dplyr::recode(
          x,
          "LVEF"                  = "LVEF",
          "Age"                   = "Age",
          "eGFR"                  = "eGFR",
          "bin_beta_blockers"     = "Beta-blockers",
          "bin_diabetes"          = "Diabetes",
          "bin_stroke_tia"        = "Stroke/TIA",
          "bin_af_atrial_flutter" = "AF/Atrial flutter",
          "bin_sex_male"          = "Sex",
          "FIS_td"                = "Inappropriate ICD therapy (time-dependent)",
          .default = x
        )
      }
    ) |>
    gt::cols_label(
      term     = "Covariate",
      estimate = "β (log HR)",
      se       = "SE",
      HR       = "HR",
      lower    = "95% CI (Lower)",
      upper    = "95% CI (Upper)",
      p.value  = "P-value"
    ) |>
    gt::fmt_number(
      columns = c(estimate, se, HR, lower, upper),
      decimals = 2
    ) |>
    gt::fmt_scientific(
      columns = p.value,
      decimals = 2
    ) |>
    gt::tab_source_note(
      source_note = "NYHA functional class included as a stratification factor in the Cox model."
    ),
  filename = file.path(OUTDIR, "04_final_model_NYHA_stratified.html"),
  inline_css = TRUE
)


# -------------------------------------------------------------------------
# 7)Power and Feasibility Checks for Time-Dependent FIS Exposure
# -------------------------------------------------------------------------

dt_check <- fix_exposure_na(standardise_types(as.data.table(complete(t_train, 1))))
td_check <- make_td_tmerge(dt_check, final_vars)

msg("\nExposure Distribution:")
msg("  Total person-intervals: %d", nrow(td_check))
msg("  Intervals with FIS_td=0: %d (%.1f%%)", 
    sum(td_check$FIS_td == 0), 100*mean(td_check$FIS_td == 0))
msg("  Intervals with FIS_td=1: %d (%.1f%%)", 
    sum(td_check$FIS_td == 1), 100*mean(td_check$FIS_td == 1))

msg("\nEvent Distribution:")
msg("  Deaths during FIS_td=0: %d", sum(td_check$death==TRUE & td_check$FIS_td==0))
msg("  Deaths during FIS_td=1: %d", sum(td_check$death==TRUE & td_check$FIS_td==1))

# Calculate crude event rate by exposure
rate_unexposed <- sum(td_check$death==TRUE & td_check$FIS_td==0) / 
  sum(td_check[FIS_td==0, tstop - tstart])
rate_exposed <- sum(td_check$death==TRUE & td_check$FIS_td==1) / 
  sum(td_check[FIS_td==1, tstop - tstart])

msg("\nCrude Event Rates:")
msg("  Unexposed: %.4f per person-day", rate_unexposed)
msg("  Exposed: %.4f per person-day", rate_exposed)
msg("  Rate ratio: %.2f", rate_exposed / rate_unexposed)

# Rule of thumb: need ~10 events per covariate
n_covariates <- length(final_vars)
events_total <- sum(td_check$death == TRUE)
events_per_covariate <- events_total / n_covariates

msg("\nSample Size Assessment:")
msg("  Total deaths: %d", events_total)
msg("  Number of covariates: %d", n_covariates)
msg("  Events per covariate: %.1f", events_per_covariate)

if (events_per_covariate < 10) {
  msg("  WARNING: <10 events per covariate (rule of thumb)")
  msg("  → Model may be overfitted")
  msg("  → Wide confidence intervals expected")
} else {
  msg("  ✓ Adequate events per covariate (≥10)")
}

# Check FIS-specific power
deaths_during_fis <- sum(td_check$death==TRUE & td_check$FIS_td==1)
if (deaths_during_fis < 10) {
  msg("\nCRITICAL: Only %d deaths during FIS exposure", deaths_during_fis)
  msg("  → Insufficient power to detect FIS effect")
  msg("  → Wide CI and non-significance expected")
  msg("  → Consider: longer follow-up, larger sample, or pooled studies")
}

# --- TEST SET ---
msg("\n2. TEST SET:")
dt_test <- fix_exposure_na(standardise_types(as.data.table(complete(t_test, 1))))
td_test <- make_td_tmerge(dt_test, character(0))

n_test <- uniqueN(dt_test$ID)
fis_test <- sum(dt_test$Status_FIS==1, na.rm=TRUE)
deaths_test <- sum(dt_test$Status_death==1, na.rm=TRUE)
deaths_fis_test <- sum(td_test$death==TRUE & td_test$FIS_td==1)

msg("   Total patients: %d", n_test)
msg("   Patients with FIS: %d (%.1f%%)", fis_test, 100*fis_test/n_test)
msg("   Total deaths: %d (%.1f%%)", deaths_test, 100*deaths_test/n_test)
msg("   Deaths during FIS exposure: %d", deaths_fis_test)


# -------------------------------------------------------------------------

# 8) INTERACTION TESTING IN FINAL MODEL (NYHA STRATIFIED) 

# -------------------------------------------------------------------------

# Interaction candidates = final baseline covariates only (exclude NYHA itself)

# NOTE: final_vars_strat is the final covariate list with NYHA removed (because it is stratified)

interaction_vars <- final_vars_strat

interaction_vars <- setdiff(interaction_vars, "NYHA")

interaction_vars <- interaction_vars[interaction_vars %in% names(dt_ref)]  # safety



msg("Interaction candidates (final-only): %d", length(interaction_vars))

msg("Variables: %s", paste(interaction_vars, collapse=", "))



# Helper: add main effect of xvar to RHS if it isn't already present

add_main_effect_if_missing <- function(rhs, xvar){
  
  present <- grepl(paste0("(^|\\s|\\+)\\s*", xvar, "\\s*($|\\s|\\+)"), rhs)
  
  if (present) return(rhs)
  
  paste(rhs, "+", xvar)
  
}



# Pooled overall p-value for interaction FIS_td:xvar (Wald pooled; factors handled as multi-df)

pool_interaction_p <- function(mids_obj, vars_for_td_backbone, rhs_backbone, xvar){
  
  
  
  coef_list <- list()
  
  vcov_list <- list()
  
  
  
  for (i in 1:mids_obj$m) {
    
    dt <- as.data.table(complete(mids_obj, i))
    
    dt <- fix_exposure_na(standardise_types(dt))
    
    
    
    if (!(xvar %in% names(dt))) next
    
    
    
    # Build TD dataset including backbone vars + xvar
    
    td_vars <- unique(c(vars_for_td_backbone, xvar))
    
    td <- make_td_tmerge(dt, td_vars)
    
    if (nrow(td) == 0) next
    
    
    
    # Ensure main effect of xvar exists (proper interaction test)
    
    rhs_use <- add_main_effect_if_missing(rhs_backbone, xvar)
    
    
    
    # Fit interaction model
    
    f_int <- as.formula(paste0(
      
      "Surv(tstart,tstop,death) ~ ",
      
      rhs_use,
      
      " + FIS_td:", xvar
      
    ))
    
    
    
    fit <- try(coxph(f_int, data=td, ties="efron"), silent=TRUE)
    
    if (inherits(fit, "try-error")) next
    
    
    
    terms <- names(coef(fit))
    
    
    
    # Find interaction coefficients (handles factors: FIS_td:xvarLEVEL)
    
    int_terms <- grep(paste0("^FIS_td:", xvar), terms, value=TRUE)
    
    if (length(int_terms) == 0) int_terms <- grep(paste0(xvar, ":FIS_td"), terms, value=TRUE)
    
    if (length(int_terms) == 0) next
    
    
    
    coef_list[[length(coef_list)+1]] <- coef(fit)[int_terms]
    
    vcov_list[[length(vcov_list)+1]] <- vcov(fit)[int_terms, int_terms, drop=FALSE]
    
  }
  
  
  
  if (length(coef_list) < 2) {
    
    return(list(p=NA_real_, df=NA_integer_, type="failed", n_imp=length(coef_list)))
    
  }
  
  
  
  k <- length(coef_list[[1]])
  
  
  
  # If factor levels drop in some imputations, term counts can differ
  
  if (any(vapply(coef_list, length, integer(1)) != k)) {
    
    return(list(p=NA_real_, df=NA_integer_, type="inconsistent_terms", n_imp=length(coef_list)))
    
  }
  
  
  
  if (k == 1) {
    
    coef_vec <- vapply(coef_list, function(z) as.numeric(z), numeric(1))
    
    var_vec  <- vapply(vcov_list, function(v) as.numeric(v), numeric(1))
    
    pooled <- pool_single_coef(coef_vec, var_vec)
    
    return(list(p=pooled$p, df=1L, type="single", n_imp=length(coef_list)))
    
  } else {
    
    coef_mat <- do.call(rbind, lapply(coef_list, as.numeric))
    
    pooled <- pool_multiple_coefs(coef_mat, vcov_list)
    
    return(list(p=pooled$p, df=pooled$df, type="factor", n_imp=length(coef_list)))
    
  }
  
}



# Backbone from STEP 6 (must exist)

vars_for_td_backbone <- vars_for_td_strat   # includes NYHA so strata(NYHA) works

rhs_backbone <- rhs_final_strat



# Run interaction tests (final-only)

interaction_tab_finalonly <- rbindlist(lapply(seq_along(interaction_vars), function(j){
  
  x <- interaction_vars[j]
  
  if (j %% 5 == 0) msg("  Progress: %d/%d", j, length(interaction_vars))
  
  
  
  res <- pool_interaction_p(
    
    mids_obj = t_train,
    
    vars_for_td_backbone = vars_for_td_backbone,
    
    rhs_backbone = rhs_backbone,
    
    xvar = x
    
  )
  
  
  
  data.table(
    
    covariate = x,
    
    p_interaction = res$p,
    
    df = res$df,
    
    test_type = res$type,
    
    n_imputations_used = res$n_imp
    
  )
  
}))



setorder(interaction_tab_finalonly, p_interaction)

fwrite(interaction_tab_finalonly, file.path(OUTDIR, "07_interaction_tests_FIS_td_FINALONLY_STRATNYHA.csv"))



msg("\nSaved: 07_interaction_tests_FIS_td_FINALONLY_STRATNYHA.csv")

msg("Interaction results:")

print(interaction_tab_finalonly)
# Save presentation-ready interaction tests (HTML) — temporary conversion + labels

gtsave(
  
  gt::gt(interaction_tab_finalonly) |>
    
    gt::text_transform(
      
      locations = gt::cells_body(columns = covariate),
      
      fn = function(x) {
        
        dplyr::recode(
          
          x,
          
          "bin_stroke_tia"        = "Stroke / TIA",
          
          "bin_sex_male"          = "Sex",
          
          "eGFR"                  = "eGFR",
          
          "Age"                   = "Age",
          
          "bin_af_atrial_flutter" = "AF / Atrial flutter",
          
          "LVEF"                  = "LVEF",
          
          "bin_beta_blockers"     = "Beta-blockers",
          
          "bin_diabetes"          = "Diabetes",
          
          .default = x
          
        )
        
      }
      
    ),
  
  filename = file.path(OUTDIR, "07_interaction_tests_FIS_td_FINALONLY_STRATNYHA.html"),
  
  inline_css = TRUE
  
)




# -------------------------------------------------------------------------

#("STEP 9: SUBGROUP ANALYSIS- SENSITIVITY TEST)

# -------------------------------------------------------------------------



# -------------------------------------------------------------------------

# 9A) Ensure subgroup variables exist consistently

# -------------------------------------------------------------------------

add_subgroup_vars <- function(dt){
  
  dt <- as.data.table(dt)
  
  
  
  # Age_group (use existing if present; else create)
  
  if (!("Age_group" %in% names(dt)) && "Age" %in% names(dt)) {
    
    dt[, Age_group := fifelse(Age < 65, "<65",
                              
                              fifelse(Age <= 75, "65-75", ">75"))]
    
  }
  
  if ("Age_group" %in% names(dt)) {
    
    dt[, Age_group := factor(Age_group, levels=c("<65","65-75",">75"))]
    
  }
  
  
  
  # LVEF categories
  
  if (!("LVEF_cat" %in% names(dt)) && "LVEF" %in% names(dt)) {
    
    dt[, LVEF_cat := fifelse(LVEF < 30, "<30",
                             
                             fifelse(LVEF <= 35, "30-35", ">35"))]
    
  }
  
  if ("LVEF_cat" %in% names(dt)) {
    
    dt[, LVEF_cat := factor(LVEF_cat, levels=c("<30","30-35",">35"))]
    
  }
  
  
  
  # NYHA binary (I-II vs III-IV)
  
  if (!("NYHA_bin" %in% names(dt)) && "NYHA" %in% names(dt)) {
    
    nyha_num <- suppressWarnings(as.integer(as.character(dt$NYHA)))
    
    if (all(is.na(nyha_num))) {
      
      nyha_chr <- toupper(as.character(dt$NYHA))
      
      nyha_num <- fifelse(nyha_chr %in% c("I","1"), 1L,
                          
                          fifelse(nyha_chr %in% c("II","2"), 2L,
                                  
                                  fifelse(nyha_chr %in% c("III","3"), 3L,
                                          
                                          fifelse(nyha_chr %in% c("IV","4"), 4L, NA_integer_))))
      
    }
    
    dt[, NYHA_bin := fifelse(nyha_num %in% c(1,2), "I-II",
                             
                             fifelse(nyha_num %in% c(3,4), "III-IV", NA_character_))]
    
    dt[, NYHA_bin := factor(NYHA_bin, levels=c("I-II","III-IV"))]
    
  }
  
  
  
  dt
  
}



# -------------------------------------------------------------------------

# 9B) Subgroup pooled HR(FIS_td) — suppress strata(NYHA) inside NYHA subgroups

# -------------------------------------------------------------------------

pool_fis_in_subset <- function(mids_obj, subset_var, subset_level,
                               
                               vars_for_td_backbone, rhs_backbone) {
  
  
  
  drop_strata_nyha <- function(rhs) {
    
    rhs <- gsub("\\+\\s*strata\\(NYHA\\)", "", rhs)
    
    rhs <- gsub("strata\\(NYHA\\)\\s*\\+\\s*", "", rhs)
    
    rhs <- gsub("strata\\(NYHA\\)", "", rhs)
    
    rhs <- gsub("\\s+", " ", rhs)
    
    trimws(rhs)
    
  }
  
  
  
  rhs_use <- rhs_backbone
  
  if (subset_var %in% c("NYHA_bin", "NYHA")) rhs_use <- drop_strata_nyha(rhs_backbone)
  
  
  
  coefs <- c(); vars <- c()
  
  
  
  for (i in 1:mids_obj$m) {
    
    dt <- as.data.table(complete(mids_obj, i))
    
    dt <- fix_exposure_na(standardise_types(dt))
    
    dt <- add_subgroup_vars(dt)
    
    
    
    if (!(subset_var %in% names(dt))) next
    
    
    
    dsub <- dt[get(subset_var) == subset_level]
    
    if (nrow(dsub) == 0) next
    
    
    
    td <- make_td_tmerge(dsub, unique(c(vars_for_td_backbone, subset_var)))
    
    if (nrow(td) == 0) next
    
    
    
    fit <- try(
      
      coxph(as.formula(paste("Surv(tstart,tstop,death) ~", rhs_use)),
            
            data = td, ties = "efron"),
      
      silent = TRUE
      
    )
    
    if (inherits(fit, "try-error")) next
    
    if (!("FIS_td" %in% names(coef(fit)))) next
    
    
    
    coefs <- c(coefs, as.numeric(coef(fit)["FIS_td"]))
    
    vars  <- c(vars,  as.numeric(vcov(fit)["FIS_td","FIS_td"]))
    
  }
  
  
  
  if (length(coefs) < 2) {
    
    return(data.table(
      
      subgroup_var = subset_var,
      
      level        = as.character(subset_level),
      
      n_imp_used   = length(coefs),
      
      HR           = NA_real_, lower = NA_real_, upper = NA_real_, p.value = NA_real_
      
    ))
    
  }
  
  
  
  pooled <- pool_single_coef(coefs, vars)
  
  
  
  data.table(
    
    subgroup_var = subset_var,
    
    level        = as.character(subset_level),
    
    n_imp_used   = length(coefs),
    
    HR           = exp(pooled$estimate),
    
    lower        = exp(pooled$estimate - 1.96 * pooled$se),
    
    upper        = exp(pooled$estimate + 1.96 * pooled$se),
    
    p.value      = pooled$p
    
  )
  
}



# -------------------------------------------------------------------------

# 9C) Run SAP subgroup analyses (train, pooled across imputations)

# -------------------------------------------------------------------------

sap_subgroups <- list(
  
  Age_group = c("<65","65-75",">75"),
  
  LVEF_cat  = c("<30","30-35",">35"),
  
  bin_sex_male = c(0,1),
  
  bin_af_atrial_flutter = c(0,1),
  
  bin_diabetes = c(0,1),
  
  NYHA_bin = c("I-II","III-IV")
  
)





# - vars_for_td_strat includes NYHA so strata(NYHA) works

# - rhs_final_strat is "FIS_td + <final covars> + strata(NYHA)"

stopifnot(exists("vars_for_td_strat"), exists("rhs_final_strat"))



vars_for_td_backbone <- vars_for_td_strat

rhs_backbone <- rhs_final_strat



msg("Running subgroup models (pooled across imputations)...")



subgroup_results <- rbindlist(lapply(names(sap_subgroups), function(v){
  
  levs <- sap_subgroups[[v]]
  
  rbindlist(lapply(levs, function(L){
    
    pool_fis_in_subset(
      
      mids_obj = t_train,
      
      subset_var = v,
      
      subset_level = L,
      
      vars_for_td_backbone = vars_for_td_backbone,
      
      rhs_backbone = rhs_backbone
      
    )
    
  }))
  
}), fill = TRUE)



# Save raw subgroup results (as before)

fwrite(subgroup_results, file.path(OUTDIR, "08_subgroup_HR_FIS_td.csv"))

msg("Saved: 08_subgroup_HR_FIS_td.csv")



# -------------------------------------------------------------------------

# 9D) Publication-style table image: Characteristic | HR | 95% CI | p-value

# -------------------------------------------------------------------------


tab <- copy(subgroup_results)



# ---- Labels ----

var_labels <- c(
  
  "Age_group" = "Age group (years)",
  
  "LVEF_cat"  = "LVEF category (%)",
  
  "bin_sex_male" = "Sex",
  
  "bin_af_atrial_flutter" = "Atrial fibrillation/flutter",
  
  "bin_diabetes" = "Diabetes",
  
  "NYHA_bin"  = "NYHA class"
  
)



format_level <- function(v, lvl) {
  
  lvl <- as.character(lvl)
  
  if (v == "bin_sex_male") return(ifelse(lvl == "1", "Male", "Female"))
  
  if (v == "bin_af_atrial_flutter") return(ifelse(lvl == "1", "Yes", "No"))
  
  if (v == "bin_diabetes") return(ifelse(lvl == "1", "Yes", "No"))
  
  if (v == "NYHA_bin") return(ifelse(lvl == "I-II", "I–II", "III–IV"))
  
  lvl
  
}



tab[, Subgroup := var_labels[subgroup_var]]

tab[, Level := mapply(format_level, subgroup_var, level)]



# ---- Keep only estimable rows ----

tab <- tab[!is.na(HR) & is.finite(HR) & !is.na(lower) & !is.na(upper) & !is.na(p.value)]



# ---- Ordering (SAP order) ----

sap_order <- c("Age_group","LVEF_cat","bin_sex_male",
               
               "bin_af_atrial_flutter","bin_diabetes","NYHA_bin")

tab[, subgroup_var := factor(subgroup_var, levels = sap_order)]

setorder(tab, subgroup_var, Level)



# ---- Build Characteristic column in the style you want ----

# e.g. "Sex: Male", "Age group (years): <65"

tab[, Characteristic := paste0(Subgroup, ": ", Level)]



# ---- Format columns ----

fmt_num <- function(x, d=2) formatC(x, digits=d, format="f")

fmt_p <- function(p) ifelse(p < 0.001, "<0.001", formatC(p, digits=3, format="f"))



out_df <- data.table(
  
  Characteristic = tab$Characteristic,
  
  HR            = fmt_num(tab$HR, 2),
  
  `95% CI`      = paste0(fmt_num(tab$lower,2), ", ", fmt_num(tab$upper,2)),
  
  `p-value`     = fmt_p(tab$p.value)
  
)



# Optional: wrap long characteristic text

wrap_text <- function(x, width=38) vapply(x, function(s) paste(strwrap(s, width=width), collapse="\n"), "")

out_df[, Characteristic := wrap_text(Characteristic, 42)]



df <- as.data.frame(out_df)



# ---- Main table grob (header styling) ----

table_theme <- ttheme_minimal(
  
  base_size = 11,
  
  core = list(fg_params = list(hjust = 0, x = 0.02)),
  
  colhead = list(
    
    fg_params = list(fontface = "bold", hjust = 0, x = 0.02),
    
    bg_params = list(fill = "#F2F2F2", col = NA)
    
  ),
  
  padding = unit(c(4, 6), "mm")
  
)



tg <- tableGrob(df, rows = NULL, theme = table_theme)



# ---- Footnote ----

foot <- textGrob(
  
  "Abbreviations: CI = Confidence Interval, HR = Hazard Ratio",
  
  x = 0, hjust = 0,
  
  gp = gpar(fontsize = 10)
  
)



# Combine table + footnote (stacked)

g <- arrangeGrob(tg, foot, ncol = 1, heights = unit.c(unit(1, "npc") - unit(1.2, "lines"), unit(1.2, "lines")))



# ---- Save PNG ----

png(file.path(OUTDIR, "08_subgroup_table_FIS_td_pubstyle.png"),
    
    width = 2400, height = 1500, res = 300)

grid.newpage(); grid.draw(g)

dev.off()



# ---- Save PDF ----

pdf(file.path(OUTDIR, "08_subgroup_table_FIS_td_pubstyle.pdf"),
    
    width = 10.5, height = 6.8)

grid.newpage(); grid.draw(g)

dev.off()



msg("Saved: 08_subgroup_table_FIS_td_pubstyle1.png")

msg("Saved: 08_subgroup_table_FIS_td_pubstyle1.pdf")


# -------------------------------------------------------------------------

# 9E) Publication-ready FOREST PLOT (drop NA rows)

# -------------------------------------------------------------------------

fp <- copy(tab)

fp <- fp[!is.na(HR) & is.finite(HR) & !is.na(lower) & !is.na(upper)]

if (nrow(fp) == 0) stop("No estimable subgroup results to plot.")



fp[, plot_label := paste0(Subgroup, " — ", Level)]

fp[, plot_label := factor(plot_label, levels = rev(unique(plot_label)))]



p <- ggplot(fp, aes(x = HR, y = plot_label)) +
  
  geom_point(size = 2) +
  
  geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2) +
  
  geom_vline(xintercept = 1, linetype = 2) +
  
  scale_x_log10() +
  
  labs(
    
    x = "Hazard ratio for all-cause mortality (Inappropriate shock vs no Inappropriate shock)",
    
    y = NULL,
    
    title = "Subgroup analyses"
    
  ) +
  
  theme_bw(base_size = 11) +
  
  theme(
    
    plot.title = element_text(face = "bold"),
    
    axis.text.y = element_text(size = 10)
    
  )



ggsave(file.path(OUTDIR, "08_forest_plot_subgroups_FIS_td_publication1.pdf"),
       
       p, width = 9.5, height = 6.5)

msg("Saved: 08_forest_plot_subgroups_FIS_td_publication1.pdf")





# -------------------------------------------------------------------------

# 10) INTERNAL VALIDATION

# -------------------------------------------------------------------------


TIMES       <- c(365, 1095, 1825)   # 3y, 5y, 10y
TIME_LABELS <- c("1 years", "3 years", "5 years")

MIN_EVENTS_TEST_PLOT <- 5   # minimum events by tau to draw KM calibration plot
CAL_GROUPS <- 10            # use 5 if unstable
CAL_TAU    <- 1095          # target horizon for calibration plot (days)

# --------------------------- HELPERS -----------------------------------------

# Summary helper for vectors
summ_vec <- function(x) {
  n <- sum(!is.na(x))
  if (n < 2) return(list(mean = NA_real_, sd = NA_real_, n = n))
  list(mean = mean(x, na.rm = TRUE), sd = sd(x, na.rm = TRUE), n = n)
}

# TD AUC using timeROC on collapsed (1 row per ID)
td_auc_timeROC <- function(T, delta, marker, times_days) {
  ok <- is.finite(T) & T > 0 & !is.na(delta) & is.finite(marker)
  T <- T[ok]; delta <- delta[ok]; marker <- marker[ok]
  if (sum(delta == 1L, na.rm = TRUE) < 2) return(rep(NA_real_, length(times_days)))
  roc <- try(timeROC(T = T, delta = delta, marker = marker,
                     cause = 1, weighting = "marginal",
                     times = times_days, iid = FALSE), silent = TRUE)
  if (inherits(roc, "try-error")) return(rep(NA_real_, length(times_days)))
  as.numeric(roc$AUC)
}

# Robust KM calibration plot for TD-Cox (start/stop)
tdcox_calibration_km <- function(fit_tr, td_te,
                                 tau,
                                 id_var = "ID",
                                 tstart_var = "tstart",
                                 tstop_var  = "tstop",
                                 death_var  = "death",
                                 groups = 10) {
  
  td_te <- as.data.table(td_te)
  
  # Collapse to 1 row per ID (use last interval row)
  setorderv(td_te, c(id_var, tstop_var))
  td1 <- td_te[, .SD[.N], by = id_var]
  
  # Subject-level end time/event
  end_dt <- td_te[, .(
    fu_time  = max(get(tstop_var), na.rm = TRUE),
    fu_event = as.integer(max(get(death_var), na.rm = TRUE) == 1L)
  ), by = id_var]
  
  td1 <- merge(td1, end_dt, by = id_var, all.x = TRUE)
  td1 <- td1[is.finite(fu_time) & fu_time > 0]
  if (nrow(td1) < 30) return(list(ok = FALSE, reason = "Too few TEST subjects after collapse."))
  
  # Choose tau not beyond follow-up
  tau_use <- min(tau, floor(stats::quantile(td1$fu_time, 0.80, na.rm = TRUE)))
  tau_use <- max(30, tau_use)
  
  # LP on collapsed TEST
  lp <- try(predict(fit_tr, newdata = td1, type = "lp"), silent = TRUE)
  if (inherits(lp, "try-error")) return(list(ok = FALSE, reason = "predict(lp) failed on collapsed TEST."))
  lp <- as.numeric(lp)
  
  # Baseline cumulative hazard at tau_use
  bh <- basehaz(fit_tr, centered = FALSE)
  bh <- bh[order(bh$time), ]
  if (nrow(bh) == 0) return(list(ok = FALSE, reason = "basehaz() empty."))
  idx <- max(which(bh$time <= tau_use))
  H0_tau <- if (!is.finite(idx)) 0 else as.numeric(bh$hazard[idx])
  
  # Predicted risk at tau_use
  pred_surv <- exp(-H0_tau * exp(lp))
  pred_risk <- 1 - pred_surv
  td1[, pred_risk := pred_risk]
  
  # Observed by tau_use: censor at tau_use
  td1[, time_tau  := pmin(fu_time, tau_use)]
  td1[, event_tau := as.integer(fu_event == 1L & fu_time <= tau_use)]
  
  if (sum(td1$event_tau, na.rm = TRUE) < MIN_EVENTS_TEST_PLOT) {
    return(list(ok = FALSE, reason = sprintf("Too few events by tau=%d for KM groups.", tau_use)))
  }
  
  # Quantile groups
  td1[, r := frank(pred_risk, ties.method = "average")]
  br <- unique(quantile(td1$r, probs = seq(0, 1, length.out = groups + 1), na.rm = TRUE))
  if (length(br) < 3) return(list(ok = FALSE, reason = "Not enough unique predicted values."))
  
  td1[, grp := cut(r, breaks = br, include.lowest = TRUE, labels = FALSE)]
  
  cal_dt <- rbindlist(lapply(sort(unique(td1$grp)), function(g) {
    dd <- td1[grp == g]
    if (nrow(dd) < 5) return(NULL)
    
    km <- survfit(Surv(time_tau, event_tau) ~ 1, data = dd)
    s_obs <- summary(km, times = tau_use)$surv
    if (length(s_obs) == 0) s_obs <- NA_real_
    
    data.table(
      group = g,
      n = nrow(dd),
      events = sum(dd$event_tau, na.rm = TRUE),
      pred = mean(dd$pred_risk, na.rm = TRUE),
      obs  = 1 - as.numeric(s_obs),
      tau_used = tau_use
    )
  }), fill = TRUE)
  
  cal_dt <- cal_dt[!is.na(obs)]
  if (nrow(cal_dt) < 2) return(list(ok = FALSE, reason = "Not enough KM calibration points."))
  
  p <- ggplot(cal_dt, aes(x = pred, y = obs)) +
    geom_abline(intercept = 0, slope = 1, linetype = 3, linewidth = 0.8) +
    geom_line(linewidth = 1) +
    geom_point(size = 2.6) +
    coord_equal(xlim = c(0, 1), ylim = c(0, 1)) +
    labs(
      title = sprintf("Calibration (TEST) — TD-Cox — tau = %d days", cal_dt$tau_used[1]),
      subtitle = sprintf("Observed vs predicted risk by %d days (KM, %d groups)",
                         cal_dt$tau_used[1], length(unique(cal_dt$group))),
      x = sprintf("Predicted risk by %d days", cal_dt$tau_used[1]),
      y = sprintf("Observed risk by %d days (KM)", cal_dt$tau_used[1])
    ) +
    theme_minimal(base_size = 13) +
    theme(plot.title = element_text(face = "bold"))
  
  list(ok = TRUE, plot = p, data = cal_dt, tau_used = cal_dt$tau_used[1])
}

# --------------------------- RUN VALIDATION ----------------------------------

cat("\n=== TD-COX Internal validation on TEST set ===\n")
cat(sprintf("Train imputations: %d | Test imputations: %d\n", m_train, m_test))
cat(sprintf("Output folder: %s\n\n", normalizePath(OUTDIR, winslash="/")))

cindex_vec <- rep(NA_real_, m_test)
calib_slope_vec <- rep(NA_real_, m_test)
calib_intercept_vec <- rep(NA_real_, m_test)

first_fit   <- NULL
first_td_te <- NULL

for (i in 1:m_test) {
  k <- ((i - 1) %% m_train) + 1
  
  # ---------------- Train (imputation k) ----------------
  dt_tr <- fix_exposure_na(standardise_types(as.data.table(complete(t_train, k))))
  td_tr <- make_td_tmerge(dt_tr, vars_for_td_strat)
  
  fit_tr <- try(
    coxph(
      as.formula(paste0("Surv(tstart, tstop, death) ~ ", rhs_final_strat)),
      data = td_tr, ties = "efron"
    ),
    silent = TRUE
  )
  if (inherits(fit_tr, "try-error")) next
  
  # ---------------- Test (imputation i) ----------------
  dt_te <- fix_exposure_na(standardise_types(as.data.table(complete(t_test, i))))
  td_te <- make_td_tmerge(dt_te, vars_for_td_strat)
  
  if (is.null(first_fit)) {
    first_fit   <- fit_tr
    first_td_te <- td_te
  }
  
  # ---------------- C-index ----------------
  lp <- try(predict(fit_tr, newdata = td_te, type = "lp"), silent = TRUE)
  if (inherits(lp, "try-error")) next
  
  cindex_vec[i] <- survConcordance(Surv(td_te$tstart, td_te$tstop, td_te$death) ~ lp)$concordance
  
  # ---------------- Calibration slope + intercept (simple recalibration) ----
  fit_calib <- try(coxph(Surv(tstart, tstop, death) ~ lp, data = td_te), silent = TRUE)
  if (!inherits(fit_calib, "try-error")) {
    calib_slope_vec[i] <- as.numeric(coef(fit_calib)[1])
    
    obs_rate  <- mean(td_te$death, na.rm = TRUE)
    pred_rate <- mean(exp(lp), na.rm = TRUE)
    if (is.finite(obs_rate) && is.finite(pred_rate) && obs_rate > 0 && pred_rate > 0) {
      calib_intercept_vec[i] <- log(obs_rate) - log(pred_rate)
    }
  }
  
  if (i %% 5 == 0 || i == 1) cat(sprintf("  processed test imputation %d/%d\n", i, m_test))
}

# Save metrics per imputation
metrics_dt <- data.table(
  imputation = 1:m_test,
  c_index = cindex_vec,
  calib_slope = calib_slope_vec,
  calib_intercept = calib_intercept_vec
)
fwrite(metrics_dt, file.path(OUTDIR, "TEST_metrics_tdcox1.csv"))

# Print summary
c_sum <- summ_vec(cindex_vec)
s_sum <- summ_vec(calib_slope_vec)
i_sum <- summ_vec(calib_intercept_vec)

cat("\n--- SUMMARY (TEST across imputations) ---\n")
cat(sprintf("C-index: mean=%.4f, sd=%.4f (n=%d)\n", c_sum$mean, c_sum$sd, c_sum$n))
cat(sprintf("Calib slope: mean=%.4f, sd=%.4f (n=%d)\n", s_sum$mean, s_sum$sd, s_sum$n))
cat(sprintf("Calib intercept: mean=%.4f, sd=%.4f (n=%d)\n\n", i_sum$mean, i_sum$sd, i_sum$n))

# ---------------- Time-dependent AUC (first imputation only) ------------------
auc_results <- NULL
if (!is.null(first_fit) && !is.null(first_td_te)) {
  
  td_te1 <- as.data.table(first_td_te)
  setorder(td_te1, ID, tstop)
  
  # collapse to 1 row per ID at last tstop for AUC calculation
  td_col <- td_te1[, .SD[.N], by = ID]
  td_col$marker <- as.numeric(predict(first_fit, newdata = td_col, type = "lp"))
  
  max_time <- max(td_col$tstop, na.rm = TRUE)
  valid_times  <- TIMES[TIMES <= max_time]
  valid_labels <- TIME_LABELS[TIMES <= max_time]
  
  if (length(valid_times) > 0) {
    auc <- try(timeROC(
      T = td_col$tstop, delta = td_col$death,
      marker = td_col$marker, cause = 1,
      times = valid_times, iid = FALSE
    ), silent = TRUE)
    
    if (!inherits(auc, "try-error")) {
      auc_results <- data.table(
        Time_Label = valid_labels,
        Time_Days  = valid_times,
        AUC        = as.numeric(auc$AUC)
      )
      fwrite(auc_results, file.path(OUTDIR, "TEST_tdAUC_tdcox.csv"))
      cat("--- Time-dependent AUC (first imputation) ---\n")
      print(auc_results)
      cat("\n")
    }
  }
}

# ---------------- KM CALIBRATION PLOT (first imputation) ---------------------
if (!is.null(first_fit) && !is.null(first_td_te)) {
  
  cal_obj <- tdcox_calibration_km(
    fit_tr = first_fit,
    td_te  = first_td_te,
    tau    = CAL_TAU,
    id_var = "ID",
    tstart_var = "tstart",
    tstop_var  = "tstop",
    death_var  = "death",
    groups = CAL_GROUPS
  )
  
  if (isTRUE(cal_obj$ok)) {
    cal_png <- file.path(OUTDIR, sprintf("TEST_calibration_tdcox_tau%d1.png", cal_obj$tau_used))
    cal_pdf <- file.path(OUTDIR, sprintf("TEST_calibration_tdcox_tau%d1.pdf", cal_obj$tau_used))
    cal_csv <- file.path(OUTDIR, sprintf("TEST_calibration_points_tdcox_tau%d1.csv", cal_obj$tau_used))
    
    ggsave(cal_png, cal_obj$plot, width = 6.8, height = 5.2, dpi = 300)
    ggsave(cal_pdf, cal_obj$plot, width = 6.8, height = 5.2)
    fwrite(cal_obj$data, cal_csv)
    
    cat(sprintf("Calibration plot saved: %s | %s\n", cal_png, cal_pdf))
  } else {
    cat("Calibration plot skipped: ", cal_obj$reason, "\n")
    cat("Tip: try groups=5 or tau=365/730 if events are low.\n")
  }
}

cat("\nDone. Outputs saved in: ", normalizePath(OUTDIR, winslash="/"), "\n")

gtsave(
  gt::gt(metrics_dt),
  filename = file.path(OUTDIR, "TEST_metrics_tdcox1.html"),
  inline_css = TRUE
)

gtsave(
  gt::gt(auc_results),
  filename = file.path(OUTDIR, "TEST_tdAUC_tdcox.html"),
  inline_css = TRUE
)

