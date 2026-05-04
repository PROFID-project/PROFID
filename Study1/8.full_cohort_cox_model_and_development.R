###############################################################################
# COX PROPORTIONAL HAZARDS MODEL — PRIMARY ANALYSIS
# Time-dependent exposure analysis with inappropriate ICD shock
#
# Analysis overview:
#   - Outcome:    All-cause mortality
#   - Exposure:   Time-dependent inappropriate ICD shock (FIS)
#   - Time scale: Days from ICD implantation to death or censoring
#
# Model development (full cohort):
#   - Candidate covariates defined a priori based on clinical relevance
#   - Correlation-based pruning to reduce multicollinearity (|r| > 0.70)
#   - Full multivariable Cox model fitted with time-dependent exposure
#   - Backward selection applied to non-forced covariates
#   - Proportional hazards assumption assessed via Schoenfeld residuals
#   - Linearity assessed via Martingale residuals for continuous covariates
#   - Cohort (DB) incorporated as stratification variable (strata(DB))
#   - Formal interaction testing for all model covariates
#   - Prespecified subgroup analyses with forest plot
#
# Pooling: Rubin's rules across 20 MICE imputations
###############################################################################





# Required packages ---------------------------------------------------------
pkgs <- c(
  "tidyverse","dplyr" ,"survival", "data.table", "mice", "naniar", "grid","gridExtra",
  "openxlsx", "readxl", "gt", "ggplot2", "timeROC","prodlim","riskRegression"
)

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))


OUTDIR <- "T:/Study_1/Supplementary_data/Primary"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

t_full <- readRDS("T:/Study_1/Imputed_data/mice_full_object1.rds")
m_full  <- t_full$m

# Derive QRS_log1p across all imputations
for (i in 1:m_full) {
  comp_i <- complete(t_full, i)
  missing_rows <- which(is.na(t_full$data$QRS))
  t_full$imp[["QRS_log1p"]] <- if (is.null(t_full$imp[["QRS_log1p"]])) {
    as.data.frame(matrix(NA, nrow = nrow(t_full$data), ncol = m_full))
  } else t_full$imp[["QRS_log1p"]]
  t_full$imp[["QRS_log1p"]][missing_rows, i] <- log1p(comp_i$QRS[missing_rows])
}

dt_ref <- as.data.table(complete(t_full, 1))

# Verify
cat("QRS_log1p exists:", "QRS_log1p" %in% names(dt_ref), "\n")
cat("QRS_log1p NAs:", sum(is.na(dt_ref$QRS_log1p)), "\n")
summary
head(dt_ref)
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
  
  "bin_hypertension", "bin_smoking", "Haemoglobin", "SBP", "BMI", "HR",
  
  "QRS_log1p",
  
  "bin_hf", "bin_stroke_tia", "bin_av_block", "bin_lbbb",
  
  "bin_ace_inhibitor", "bin_ace_inhibitor_arb", "bin_diuretics",
  
  "bin_anti_arrhythmic_iii", "bin_anti_coagulant", "bin_anti_platelet",
  
  "bin_lipid_lowering", "bin_anti_diabetic_oral", "bin_digitalis_glycosides",
  
  "bin_pci", "HasMRI"
  
)



all_covars <- unique(c(forced_vars, candidate_vars))



# Learn NYHA levels from full cohort imp1

nyha_levels <- if ("NYHA" %in% names(dt_ref)) levels(factor(dt_ref$NYHA)) else NULL

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
  
  num_candidates <- c("Age","LVEF","eGFR","Haemoglobin","SBP","BMI","HR","CRP_log1p")
  
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
  # NYHA_bin — binary collapse — created HERE before tmerge
  # Used ONLY in interaction testing and subgroup analyses
  # Primary model uses 4-level NYHA — NOT affected
  # FIXED: use string labels not integer conversion
  if ("NYHA" %in% names(dt) && !("NYHA_bin" %in% names(dt))) {
    nyha_chr <- as.character(dt$NYHA)
    dt[, NYHA_bin := fifelse(nyha_chr %in% c("I","II"), "I-II",
                             fifelse(nyha_chr %in% c("III","IV"), "III-IV", NA_character_))]
    dt[, NYHA_bin := factor(NYHA_bin, levels=c("I-II","III-IV"))]
  }
  # DB as factor (cohort indicator for stratification)
  if ("DB" %in% names(dt)) {
    dt[, DB := factor(as.character(DB))]
  }
  
  # Derive QRS_log1p if QRS exists but QRS_log1p not present
  if ("QRS" %in% names(dt) && !"QRS_log1p" %in% names(dt)) {
    dt[, QRS_log1p := log1p(as.numeric(QRS))]
  }
  
  dt
  
}


# -------------------------------------------------------------------------
# 3b) Restructuring Data to Long (Start–Stop) Format Using tmerge
# -------------------------------------------------------------------------
make_td_tmerge <- function(dt, vars_base){
  dt <- as.data.table(dt)
  
  # require outcome
  dt <- dt[!is.na(Time_death_days) & !is.na(Status_death)]
  
  
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

dt0 <- standardise_types(copy(dt_ref))
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
  n_unexposed = sum(Status_FIS == 0L, na.rm=TRUE)
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
    dt <- standardise_types(dt)
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
    dt <- standardise_types(dt)
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
fis_crude <- fit_and_pool_all(t_full, character(0), "FIS_td")

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
  pooled_model <- fit_and_pool_all(t_full, x, paste("FIS_td +", x))
  
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
    result_wald <- pool_variable(t_full, x, x)
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
    pooled <- fit_and_pool_all(t_full, x, paste("FIS_td +", x))
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



dt_num <- standardise_types(as.data.table(complete(t_full, 1)))

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
full_pooled <- fit_and_pool_all(t_full, model_vars, rhs_full)

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
    result <- pool_variable(t_full, vars_current, v)
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
final_pooled <- fit_and_pool_all(t_full, final_vars, rhs_final)

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
          "QRS_loglp"             = "QRS duration (ms)" ,
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
# 5) STEP 6: Assumptions 
# -------------------------------------------------------------------------

dt1 <- standardise_types(as.data.table(complete(t_full, 1)))
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


################################################################################

# COX MODEL — STEP 7 ONWARDS

# Continuation after PH assumption + linearity checks (Step 6 / STEP 6 in

# your existing script). Paste this immediately after your linearity block.

#

# Assumes the following objects exist from earlier steps:

#   t_full          — mids object (full cohort, 20 imputations)

#   final_vars      — character vector of selected covariates after backward selection

#   OUTDIR          — output directory path

#   msg()           — logging helper function

#   standardise_types(), make_td_tmerge()  — helper functions

#   fit_and_pool_all(), pool_single_coef(),

#   pool_multiple_coefs(), pool_variable()  — pooling functions

################################################################################





# -------------------------------------------------------------------------

# 7) PRIMARY MODEL: NYHA as covariate + strata(DB)

#

#    Rationale:

#    - PH assumption checked in Step 6 — not violated for any covariate

#    - NYHA not violating PH (Schoenfeld p=0.12, GLOBAL p=0.38)

#      => retained as regular covariate, HR reported in Table 2

#    - DB (cohort) stratified to account for between-cohort heterogeneity

#      in baseline hazard (enrollment era, device programming, follow-up

#      intensity) — no HR estimated for DB in primary model

# -------------------------------------------------------------------------



msg("\n=== STEP 7: PRIMARY MODEL: NYHA covariate + strata(DB) ===")



# DB must be in vars_base so tmerge can pass it to strata(DB)

vars_for_td_primary <- unique(c(final_vars, "DB"))



rhs_primary <- paste(
  
  "FIS_td +",
  
  paste(final_vars, collapse = " + "),
  
  "+ strata(DB)"
  
)



msg("Primary model formula: Surv(tstart,tstop,death) ~ %s", rhs_primary)



final_pooled_primary <- fit_and_pool_all(
  
  mids_obj    = t_full,
  
  vars_base   = vars_for_td_primary,
  
  formula_rhs = rhs_primary
  
)



if (is.null(final_pooled_primary)) stop("Primary model (strata DB) failed!")



setorder(final_pooled_primary, p.value)

fwrite(final_pooled_primary, file.path(OUTDIR, "07_primary_model_strata_DB.csv"))



msg("\n*** PRIMARY MODEL: FIS_td RESULT ***")

msg("HR: %.3f (95%% CI: %.3f-%.3f)",
    
    final_pooled_primary[term=="FIS_td", HR],
    
    final_pooled_primary[term=="FIS_td", lower],
    
    final_pooled_primary[term=="FIS_td", upper])

msg("P-value: %.4f", final_pooled_primary[term=="FIS_td", p.value])



msg("\nAll coefficients:")

print(final_pooled_primary[, .(term, HR, lower, upper, p.value)])



# Save formatted HTML table (Table 2)

gtsave(
  
  gt::gt(final_pooled_primary) |>
    
    gt::text_transform(
      
      locations = gt::cells_body(columns = term),
      
      fn = function(x) {
        
        dplyr::recode(x,
                      
                      "FIS_td"                = "Inappropriate ICD therapy (time-dependent)",
                      
                      "Age"                   = "Age (years)",
                      
                      "LVEF"                  = "LVEF (%)",
                      
                      "eGFR"                  = "eGFR (mL/min/1.73m²)",
                      "QRS_loglp"             = "QRS duration (ms)" ,
                      "bin_sex_male"          = "Sex (male)",
                      
                      "bin_diabetes"          = "Diabetes mellitus",
                      
                      "bin_beta_blockers"     = "Beta-blocker therapy",
                      
                      "bin_af_atrial_flutter" = "Atrial fibrillation/flutter",
                      
                      "bin_stroke_tia"        = "Prior stroke/TIA",
                      
                      "NYHA2"                 = "NYHA class II",
                      
                      "NYHA3"                 = "NYHA class III",
                      
                      "NYHA4"                 = "NYHA class IV",
                      
                      .default = x
                      
        )
        
      }
      
    ) |>
    
    gt::cols_label(
      
      term     = "Covariate",
      
      estimate = "\u03b2 (log HR)",
      
      se       = "SE",
      
      HR       = "HR",
      
      lower    = "95% CI (Lower)",
      
      upper    = "95% CI (Upper)",
      
      p.value  = "P-value"
      
    ) |>
    
    gt::fmt_number(columns = c(estimate, se, HR, lower, upper), decimals = 2) |>
    
    gt::fmt_scientific(columns = p.value, decimals = 2) |>
    
    gt::tab_source_note(
      
      source_note = paste(
        
        "Cohort (DB) included as stratification variable — baseline hazard allowed to vary by cohort.",
        
        "NYHA class I is the reference category.",
        
        "Inappropriate ICD therapy modelled as time-dependent exposure.",
        
        sep = " "
        
      )
      
    ),
  
  filename   = file.path(OUTDIR, "07_primary_model_strata_DB.html"),
  
  inline_css = TRUE
  
)



msg("Saved: 07_primary_model_strata_DB.csv/.html")










# -------------------------------------------------------------------------

# 8) Power and Feasibility Checks

#    Diagnostic sanity check on the long-format data (imp1 only)

#    Documents: exposure distribution, deaths during FIS exposure,

#    events per parameter (EPP) — context for interpreting null result

# -------------------------------------------------------------------------



msg("\n=== STEP 8: Power and Feasibility Checks ===")



dt_check <- standardise_types(as.data.table(complete(t_full, 1)))

td_check <- make_td_tmerge(dt_check, final_vars)



n_intervals_unexposed <- sum(td_check$FIS_td == 0)

n_intervals_exposed   <- sum(td_check$FIS_td == 1)

deaths_unexposed      <- sum(td_check$death == TRUE & td_check$FIS_td == 0)

deaths_exposed        <- sum(td_check$death == TRUE & td_check$FIS_td == 1)

total_deaths          <- sum(td_check$death == TRUE)

n_params              <- length(final_vars) + 1   # +1 for FIS_td



msg("Exposure distribution (long format, imp1):")

msg("  Intervals FIS_td=0 (unexposed): %d (%.1f%%)",
    
    n_intervals_unexposed, 100 * n_intervals_unexposed / nrow(td_check))

msg("  Intervals FIS_td=1 (exposed):   %d (%.1f%%)",
    
    n_intervals_exposed, 100 * n_intervals_exposed / nrow(td_check))

msg("  Deaths during FIS_td=0: %d", deaths_unexposed)

msg("  Deaths during FIS_td=1: %d", deaths_exposed)

msg("  Total deaths: %d", total_deaths)

msg("  Parameters in model (incl. FIS_td): %d", n_params)

msg("  Events per parameter (EPP): %.1f", total_deaths / n_params)



if (deaths_exposed < 10) {
  
  msg("\nCRITICAL: Only %d deaths during FIS exposure.", deaths_exposed)
  
  msg("         Wide confidence intervals expected — interpret HR with caution.")
  
}

if (total_deaths / n_params < 10) {
  
  msg("WARNING: EPP = %.1f (below recommended threshold of 10).", total_deaths / n_params)
  
}


# =========================================================================

# STEP 9: INTERACTION TESTING

# =========================================================================



msg("\n=== STEP 9: INTERACTION TESTING ===")



# Helper: add main effect to rhs if not already present

add_main_effect_if_missing <- function(rhs, xvar) {
  
  if (grepl(paste0("\\b", xvar, "\\b"), rhs)) return(rhs)
  
  paste(rhs, "+", xvar)
  
}



# Pool interaction p-value across imputations (Rubin's rules)

pool_interaction_p <- function(mids_obj, vars_for_td_backbone, rhs_backbone, xvar) {
  
  coef_list <- list()
  
  vcov_list <- list()
  
  
  
  for (i in 1:mids_obj$m) {
    
    dt <- as.data.table(complete(mids_obj, i))
    
    dt <- standardise_types(dt)
    
    
    
    # NYHA_bin: replace 4-level NYHA with binary in backbone
    
    # to avoid collinearity — gsub places it as main effect already
    
    if (xvar == "NYHA_bin") {
      
      nyha_chr <- toupper(as.character(dt$NYHA))
      
      dt[, NYHA_bin := fifelse(nyha_chr %in% c("I","1"),   "I-II",
                               
                               fifelse(nyha_chr %in% c("II","2"),  "I-II",
                                       
                                       fifelse(nyha_chr %in% c("III","3"), "III-IV",
                                               
                                               fifelse(nyha_chr %in% c("IV","4"),  "III-IV",
                                                       
                                                       NA_character_))))]
      
      dt[, NYHA_bin := factor(NYHA_bin, levels=c("I-II","III-IV"))]
      
      td_vars_use <- unique(c(setdiff(vars_for_td_backbone, "NYHA"), "NYHA_bin"))
      
      rhs_use     <- gsub("\\bNYHA\\b", "NYHA_bin", rhs_backbone)
      
      f_int       <- as.formula(paste0("Surv(tstart,tstop,death) ~ ",
                                       
                                       rhs_use, " + FIS_td:NYHA_bin"))
      
    } else {
      
      td_vars_use <- unique(c(vars_for_td_backbone, xvar))
      
      rhs_use     <- add_main_effect_if_missing(rhs_backbone, xvar)
      
      f_int       <- as.formula(paste0("Surv(tstart,tstop,death) ~ ",
                                       
                                       rhs_use, " + FIS_td:", xvar))
      
    }
    
    
    
    if (!(xvar %in% names(dt))) next
    
    td <- make_td_tmerge(dt, td_vars_use)
    
    if (nrow(td) == 0) next
    
    
    
    fit <- try(coxph(f_int, data=td, ties="efron"), silent=TRUE)
    
    if (inherits(fit, "try-error")) next
    
    
    
    terms     <- names(coef(fit))
    
    int_terms <- grep(paste0("^FIS_td:", xvar), terms, value=TRUE)
    
    if (length(int_terms) == 0)
      
      int_terms <- grep(paste0(xvar, ":FIS_td"), terms, value=TRUE)
    
    if (length(int_terms) == 0) next
    
    
    
    coef_list[[length(coef_list)+1]] <- coef(fit)[int_terms]
    
    vcov_list[[length(vcov_list)+1]] <- vcov(fit)[int_terms, int_terms, drop=FALSE]
    
  }
  
  
  
  if (length(coef_list) < 2)
    
    return(list(p=NA_real_, df=NA_integer_, type="failed",
                
                n_imp=length(coef_list)))
  
  
  
  k <- length(coef_list[[1]])
  
  if (any(vapply(coef_list, length, integer(1)) != k))
    
    return(list(p=NA_real_, df=NA_integer_, type="inconsistent_terms",
                
                n_imp=length(coef_list)))
  
  
  
  if (k == 1) {
    
    pooled <- pool_single_coef(
      
      vapply(coef_list, as.numeric, numeric(1)),
      
      vapply(vcov_list, as.numeric, numeric(1))
      
    )
    
    return(list(p=pooled$p, df=1L, type="single", n_imp=length(coef_list)))
    
  } else {
    
    pooled <- pool_multiple_coefs(
      
      do.call(rbind, lapply(coef_list, as.numeric)), vcov_list
      
    )
    
    return(list(p=pooled$p, df=pooled$df, type="factor",
                
                n_imp=length(coef_list)))
    
  }
  
}



# Interaction candidates: all final vars except 4-level NYHA, add NYHA_bin

interaction_vars <- setdiff(final_vars, "NYHA")

interaction_vars <- interaction_vars[
  
  interaction_vars %in% names(as.data.table(complete(t_full, 1)))
  
]

interaction_vars <- c(interaction_vars, "NYHA_bin")



msg("Interaction candidates (%d): %s",
    
    length(interaction_vars), paste(interaction_vars, collapse=", "))



interaction_tab <- rbindlist(lapply(seq_along(interaction_vars), function(j) {
  
  x <- interaction_vars[j]
  
  if (j %% 3 == 0) msg("  Progress: %d/%d", j, length(interaction_vars))
  
  res <- pool_interaction_p(
    
    mids_obj             = t_full,
    
    vars_for_td_backbone = vars_for_td_primary,
    
    rhs_backbone         = rhs_primary,
    
    xvar                 = x
    
  )
  
  data.table(
    
    covariate          = x,
    
    p_interaction      = res$p,
    
    df                 = res$df,
    
    test_type          = res$type,
    
    n_imputations_used = res$n_imp
    
  )
  
}))



setorder(interaction_tab, p_interaction)

fwrite(interaction_tab, file.path(OUTDIR, "08_interaction_tests.csv"))

gtsave(gt::gt(interaction_tab),
       
       filename   = file.path(OUTDIR, "08_interaction_tests.html"),
       
       inline_css = TRUE)



msg("\nInteraction results:")

print(interaction_tab)

msg("Saved: 08_interaction_tests.csv/.html")





# =========================================================================

# STEP 10: SUBGROUP ANALYSES

# =========================================================================



msg("\n=== STEP 10: SUBGROUP ANALYSES ===")



# -------------------------------------------------------------------------

# Helper: create subgroup variables in wide data before tmerge

# -------------------------------------------------------------------------

make_subgroup_vars <- function(dt) {
  
  dt <- as.data.table(dt)
  
  
  
  if (!("Age_group" %in% names(dt)) && "Age" %in% names(dt))
    
    dt[, Age_group := fifelse(Age < 65, "<65",
                              
                              fifelse(Age <= 75, "65-75", ">75"))]
  
  if ("Age_group" %in% names(dt))
    
    dt[, Age_group := factor(Age_group, levels=c("<65","65-75",">75"))]
  
  
  
  if (!("LVEF_cat" %in% names(dt)) && "LVEF" %in% names(dt))
    
    dt[, LVEF_cat := fifelse(LVEF < 30, "<30",
                             
                             fifelse(LVEF <= 35, "30-35", ">35"))]
  
  if ("LVEF_cat" %in% names(dt))
    
    dt[, LVEF_cat := factor(LVEF_cat, levels=c("<30","30-35",">35"))]
  
  
  
  # NYHA binary — string matching (NYHA stored as "I","II","III","IV")
  
  if (!("NYHA_bin" %in% names(dt)) && "NYHA" %in% names(dt)) {
    
    nyha_chr <- toupper(as.character(dt$NYHA))
    
    dt[, NYHA_bin := fifelse(nyha_chr %in% c("I","1"),   "I-II",
                             
                             fifelse(nyha_chr %in% c("II","2"),  "I-II",
                                     
                                     fifelse(nyha_chr %in% c("III","3"), "III-IV",
                                             
                                             fifelse(nyha_chr %in% c("IV","4"),  "III-IV",
                                                     
                                                     NA_character_))))]
    
    dt[, NYHA_bin := factor(NYHA_bin, levels=c("I-II","III-IV"))]
    
  }
  
  
  
  dt
  
}



# -------------------------------------------------------------------------

# Helper: fit Cox model within one subgroup level, pool across imputations

# -------------------------------------------------------------------------

pool_fis_in_subset <- function(mids_obj, subset_var, subset_level,
                               
                               vars_for_td_backbone, rhs_backbone) {
  
  
  
  rhs_use <- rhs_backbone
  
  if (subset_var == "DB") {
    
    rhs_use <- gsub("\\+\\s*strata\\(DB\\)", "", rhs_use)
    
    rhs_use <- gsub("strata\\(DB\\)\\s*\\+\\s*", "", rhs_use)
    
    rhs_use <- gsub("strata\\(DB\\)", "", rhs_use)
    
    rhs_use <- trimws(gsub("\\s+", " ", rhs_use))
    
  }
  
  
  
  coefs     <- c()
  
  vars      <- c()
  
  n_total   <- NA_integer_
  
  n_exposed <- NA_integer_
  
  n_events  <- NA_integer_
  
  
  
  for (i in 1:mids_obj$m) {
    
    dt <- as.data.table(complete(mids_obj, i))
    
    dt <- standardise_types(dt)
    
    dt <- make_subgroup_vars(dt)
    
    
    
    if (!(subset_var %in% names(dt))) next
    
    dsub <- dt[get(subset_var) == subset_level]
    
    if (nrow(dsub) == 0) next
    
    
    
    if (i == 1) {
      
      n_total   <- nrow(dsub)
      
      n_exposed <- sum(dsub$Status_FIS   == 1L, na.rm=TRUE)
      
      n_events  <- sum(dsub$Status_death == 1L, na.rm=TRUE)
      
    }
    
    
    
    td <- make_td_tmerge(dsub, unique(c(vars_for_td_backbone, subset_var)))
    
    if (nrow(td) == 0) next
    
    
    
    fit <- try(
      
      coxph(as.formula(paste("Surv(tstart,tstop,death) ~", rhs_use)),
            
            data=td, ties="efron"),
      
      silent=TRUE
      
    )
    
    if (inherits(fit, "try-error")) next
    
    if (!("FIS_td" %in% names(coef(fit)))) next
    
    
    
    coefs <- c(coefs, as.numeric(coef(fit)["FIS_td"]))
    
    vars  <- c(vars,  as.numeric(vcov(fit)["FIS_td","FIS_td"]))
    
  }
  
  
  
  if (length(coefs) < 2) {
    
    return(data.table(
      
      subgroup_var=subset_var, level=as.character(subset_level),
      
      n_total=n_total, n_exposed=n_exposed, n_events=n_events,
      
      n_imp_used=length(coefs),
      
      HR=NA_real_, lower=NA_real_, upper=NA_real_, p.value=NA_real_
      
    ))
    
  }
  
  
  
  pooled <- pool_single_coef(coefs, vars)
  
  data.table(
    
    subgroup_var = subset_var,
    
    level        = as.character(subset_level),
    
    n_total      = n_total,
    
    n_exposed    = n_exposed,
    
    n_events     = n_events,
    
    n_imp_used   = length(coefs),
    
    HR           = exp(pooled$estimate),
    
    lower        = exp(pooled$estimate - 1.96 * pooled$se),
    
    upper        = exp(pooled$estimate + 1.96 * pooled$se),
    
    p.value      = pooled$p
    
  )
  
}



# -------------------------------------------------------------------------

# SAP subgroups — run ALL levels (both Male and Female etc.)

# Filtering for plot done separately below

# -------------------------------------------------------------------------

sap_subgroups <- list(
  
  Age_group             = c("<65","65-75",">75"),
  
  LVEF_cat              = c("<30","30-35",">35"),
  
  bin_sex_male          = c(0, 1),
  
  bin_af_atrial_flutter = c(0, 1),
  
  bin_diabetes          = c(0, 1),
  
  NYHA_bin              = c("I-II","III-IV")
  
)



msg("Running subgroup analyses...")



subgroup_results <- rbindlist(lapply(names(sap_subgroups), function(v) {
  
  rbindlist(lapply(sap_subgroups[[v]], function(L) {
    
    pool_fis_in_subset(
      
      mids_obj             = t_full,
      
      subset_var           = v,
      
      subset_level         = L,
      
      vars_for_td_backbone = vars_for_td_primary,
      
      rhs_backbone         = rhs_primary
      
    )
    
  }))
  
}), fill=TRUE)



fwrite(subgroup_results, file.path(OUTDIR, "09_subgroup_HR_FIS_td.csv"))

msg("Saved: 09_subgroup_HR_FIS_td.csv")



# -------------------------------------------------------------------------

# Labels

# -------------------------------------------------------------------------

var_labels <- c(
  
  "Age_group"             = "Age group (years)",
  
  "LVEF_cat"              = "LVEF category (%)",
  
  "bin_sex_male"          = "Sex",
  
  "bin_af_atrial_flutter" = "Atrial fibrillation/flutter",
  
  "bin_diabetes"          = "Diabetes mellitus",
  
  "NYHA_bin"              = "NYHA class"
  
)



format_level <- function(v, lvl) {
  
  lvl <- as.character(lvl)
  
  if (v == "bin_sex_male")          return(ifelse(lvl=="1", "Male", "Female"))
  
  if (v == "bin_af_atrial_flutter") return(ifelse(lvl=="1", "Yes", "No"))
  
  if (v == "bin_diabetes")          return(ifelse(lvl=="1", "Yes", "No"))
  
  if (v == "NYHA_bin")              return(ifelse(lvl=="I-II", "I\u2013II", "III\u2013IV"))
  
  lvl
  
}



tab <- copy(subgroup_results)

tab[, Subgroup := var_labels[subgroup_var]]

tab[, Level    := mapply(format_level, subgroup_var, level)]

tab <- tab[!is.na(HR) & is.finite(HR) & !is.na(lower) & !is.na(upper) & !is.na(p.value)]



sap_order <- c("Age_group","LVEF_cat","bin_sex_male",
               
               "bin_af_atrial_flutter","bin_diabetes","NYHA_bin")

tab[, subgroup_var := factor(subgroup_var, levels=sap_order)]

setorder(tab, subgroup_var, Level)



# Attach p-interaction from Step 9

pint_map <- data.table(
  
  subgroup_var = names(sap_subgroups),
  
  covariate    = c("Age_group","LVEF_cat","bin_sex_male",
                   
                   "bin_af_atrial_flutter","bin_diabetes","NYHA_bin")
  
)

pint_vals <- merge(pint_map,
                   
                   interaction_tab[, .(covariate, p_interaction)],
                   
                   by="covariate", all.x=TRUE)

tab <- merge(tab, pint_vals[, .(subgroup_var, p_interaction)],
             
             by="subgroup_var", all.x=TRUE)



# =========================================================================

# ACADEMIC TABLE — all subgroup levels including Female

# p-interaction shown once per subgroup variable

# =========================================================================



fmt2 <- function(x) formatC(x, digits=2, format="f")

fmt3 <- function(x) ifelse(x < 0.001, "<0.001", formatC(x, digits=3, format="f"))



tbl_full <- rbindlist(lapply(unique(tab$subgroup_var), function(sv) {
  
  
  
  pv     <- tab$p_interaction[tab$subgroup_var == sv][1]
  
  pv_fmt <- ifelse(is.na(pv), "", fmt3(pv))
  
  
  
  # Bold header row with p-interaction
  
  hdr <- data.frame(
    
    Subgroup        = var_labels[as.character(sv)],
    
    N               = "",
    
    `Exposed (FIS)` = "",
    
    Deaths          = "",
    
    HR              = "",
    
    `95% CI`        = "",
    
    `p-value`       = "",
    
    `p-interaction` = pv_fmt,
    
    is_hdr          = TRUE,
    
    check.names=FALSE, stringsAsFactors=FALSE
    
  )
  
  
  
  # Data rows — indented
  
  grp <- tab[subgroup_var == sv]
  
  dat <- data.frame(
    
    Subgroup        = paste0("  ", grp$Level),
    
    N               = ifelse(is.na(grp$n_total),   "", as.character(grp$n_total)),
    
    `Exposed (FIS)` = ifelse(is.na(grp$n_exposed), "", as.character(grp$n_exposed)),
    
    Deaths          = ifelse(is.na(grp$n_events),  "", as.character(grp$n_events)),
    
    HR              = fmt2(grp$HR),
    
    `95% CI`        = paste0(fmt2(grp$lower), "\u2013", fmt2(grp$upper)),
    
    `p-value`       = fmt3(grp$p.value),
    
    `p-interaction` = "",
    
    is_hdr          = FALSE,
    
    check.names=FALSE, stringsAsFactors=FALSE
    
  )
  
  
  
  rbind(hdr, dat)
  
}), fill=TRUE)



is_hdr_tbl  <- tbl_full$is_hdr

tbl_display <- as.data.frame(tbl_full[, !("is_hdr"), with=FALSE])



tbl_theme <- ttheme_minimal(
  
  base_size = 11,
  
  core = list(
    
    fg_params = list(
      
      hjust    = 0, x = 0.02,
      
      fontface = ifelse(is_hdr_tbl, "bold", "plain")
      
    ),
    
    bg_params = list(
      
      fill = ifelse(is_hdr_tbl, "grey90",
                    
                    rep(c("white","grey97"), length.out=nrow(tbl_full))),
      
      col  = NA
      
    )
    
  ),
  
  colhead = list(
    
    fg_params = list(fontface="bold", hjust=0, x=0.02),
    
    bg_params = list(fill="#DDEEFF", col=NA)
    
  ),
  
  padding = unit(c(4,6), "mm")
  
)



tg_tbl <- tableGrob(tbl_display, rows=NULL, theme=tbl_theme)



foot_tbl <- textGrob(
  
  paste(
    
    "HR = Hazard Ratio; CI = Confidence Interval; FIS = patients with inappropriate ICD shock; Deaths = all-cause deaths in subgroup.",
    
    "Each HR reflects the within-subgroup association between inappropriate ICD shock and all-cause mortality.",
    
    "p-interaction tests whether the effect of shock differs across subgroup levels (Step 9); reported once per subgroup variable.",
    
    sep="\n"
    
  ),
  
  x=0, hjust=0, gp=gpar(fontsize=9, col="grey30")
  
)



g_tbl <- arrangeGrob(tg_tbl, foot_tbl, ncol=1,
                     
                     heights=unit.c(unit(1,"npc")-unit(2.5,"lines"), unit(2.5,"lines")))



png(file.path(OUTDIR, "09_subgroup_table.png"),
    
    width=3600, height=nrow(tbl_full)*220+500, res=300)

grid.newpage(); grid.draw(g_tbl); dev.off()



pdf(file.path(OUTDIR, "09_subgroup_table.pdf"),
    
    width=14, height=max(7, nrow(tbl_full)*0.42+1.5))

grid.newpage(); grid.draw(g_tbl); dev.off()



msg("Saved: 09_subgroup_table.png/.pdf")





# =========================================================================

# STEP 10: ALL OUTPUTS

# 1. Academic table  — all subgroup levels, grey headers, supplementary

# 2. Forest plot     — clean, Male only for sex, Figure C in manuscript

# 3. HR side table   — matches forest plot rows exactly, sits beside it

# =========================================================================



# -------------------------------------------------------------------------

# Labels and formatting helpers

# -------------------------------------------------------------------------

var_labels <- c(
  
  "Age_group"             = "Age group (years)",
  
  "LVEF_cat"              = "LVEF category (%)",
  
  "bin_sex_male"          = "Sex",
  
  "bin_af_atrial_flutter" = "Atrial fibrillation/flutter",
  
  "bin_diabetes"          = "Diabetes mellitus",
  
  "NYHA_bin"              = "NYHA class"
  
)



format_level <- function(v, lvl) {
  
  lvl <- as.character(lvl)
  
  if (v == "bin_sex_male")          return(ifelse(lvl=="1", "Male", "Female"))
  
  if (v == "bin_af_atrial_flutter") return(ifelse(lvl=="1", "Yes", "No"))
  
  if (v == "bin_diabetes")          return(ifelse(lvl=="1", "Yes", "No"))
  
  if (v == "NYHA_bin")              return(ifelse(lvl=="I-II", "I\u2013II", "III\u2013IV"))
  
  lvl
  
}



fmt2 <- function(x) formatC(x, digits=2, format="f")

fmt3 <- function(x) ifelse(x < 0.001, "<0.001", formatC(x, digits=3, format="f"))



# -------------------------------------------------------------------------

# Build tab — labelled, ordered, p-interaction attached

# -------------------------------------------------------------------------

tab <- copy(subgroup_results)

tab[, Subgroup := var_labels[subgroup_var]]

tab[, Level    := mapply(format_level, subgroup_var, level)]

tab <- tab[!is.na(HR) & is.finite(HR) & !is.na(lower) & !is.na(upper) & !is.na(p.value)]



sap_order <- c("Age_group","LVEF_cat","NYHA_bin",
               
               "bin_af_atrial_flutter","bin_diabetes","bin_sex_male")

tab[, subgroup_var := factor(subgroup_var, levels=sap_order)]

setorder(tab, subgroup_var, Level)



# Attach p-interaction from Step 9

pint_map <- data.table(
  
  subgroup_var = c("Age_group","LVEF_cat","NYHA_bin",
                   
                   "bin_af_atrial_flutter","bin_diabetes","bin_sex_male"),
  
  covariate    = c("Age_group","LVEF_cat","NYHA_bin",
                   
                   "bin_af_atrial_flutter","bin_diabetes","bin_sex_male")
  
)

pint_vals <- merge(pint_map,
                   
                   interaction_tab[, .(covariate, p_interaction)],
                   
                   by="covariate", all.x=TRUE)

tab <- merge(tab, pint_vals[, .(subgroup_var, p_interaction)],
             
             by="subgroup_var", all.x=TRUE)

tab[, subgroup_var := factor(subgroup_var, levels=sap_order)]

setorder(tab, subgroup_var, Level)



# =========================================================================

# OUTPUT 1: ACADEMIC TABLE

# All subgroup levels (including Female) — for supplementary

# Columns: Subgroup / N / Exposed (FIS) / Deaths / HR / 95% CI / p-value

# Grey column header, grey bold subgroup header rows, alternating data rows

# =========================================================================



tbl_full <- rbindlist(lapply(unique(tab$subgroup_var), function(sv) {
  
  
  
  # Bold header row
  
  hdr <- data.frame(
    
    Subgroup        = var_labels[as.character(sv)],
    
    N               = "",
    
    `Exposed (FIS)` = "",
    
    Deaths          = "",
    
    HR              = "",
    
    `95% CI`        = "",
    
    `p-value`       = "",
    
    is_hdr          = TRUE,
    
    check.names=FALSE, stringsAsFactors=FALSE
    
  )
  
  
  
  # Data rows — indented
  
  grp <- tab[subgroup_var == sv]
  
  dat <- data.frame(
    
    Subgroup        = paste0("  ", grp$Level),
    
    N               = ifelse(is.na(grp$n_total),   "", as.character(grp$n_total)),
    
    `Exposed (FIS)` = ifelse(is.na(grp$n_exposed), "", as.character(grp$n_exposed)),
    
    Deaths          = ifelse(is.na(grp$n_events),  "", as.character(grp$n_events)),
    
    HR              = fmt2(grp$HR),
    
    `95% CI`        = paste0(fmt2(grp$lower), "\u2013", fmt2(grp$upper)),
    
    `p-value`       = fmt3(grp$p.value),
    
    is_hdr          = FALSE,
    
    check.names=FALSE, stringsAsFactors=FALSE
    
  )
  
  
  
  rbind(hdr, dat)
  
}), fill=TRUE)



is_hdr_tbl  <- tbl_full$is_hdr

tbl_display <- as.data.frame(tbl_full[, !("is_hdr"), with=FALSE])



tbl_theme <- ttheme_minimal(
  
  base_size = 11,
  
  core = list(
    
    fg_params = list(
      
      hjust    = 0, x = 0.02,
      
      fontface = ifelse(is_hdr_tbl, "bold", "plain")
      
    ),
    
    bg_params = list(
      
      fill = ifelse(is_hdr_tbl, "grey88",
                    
                    rep(c("white","grey97"), length.out=nrow(tbl_full))),
      
      col  = NA
      
    )
    
  ),
  
  colhead = list(
    
    fg_params = list(fontface="bold", hjust=0, x=0.02),
    
    bg_params = list(fill="grey78", col=NA)
    
  ),
  
  padding = unit(c(4, 6), "mm")
  
)



tg_tbl <- tableGrob(tbl_display, rows=NULL, theme=tbl_theme)



foot_tbl <- textGrob(
  
  paste(
    
    "HR = Hazard Ratio; CI = Confidence Interval; FIS = patients with inappropriate ICD shock; Deaths = all-cause deaths in subgroup.",
    
    "Each HR reflects the within-subgroup association between inappropriate ICD shock and all-cause mortality.",
    
    sep="\n"
    
  ),
  
  x=0, hjust=0, gp=gpar(fontsize=9, col="grey30")
  
)



g_tbl <- arrangeGrob(tg_tbl, foot_tbl, ncol=1,
                     
                     heights=unit.c(unit(1,"npc")-unit(2,"lines"), unit(2,"lines")))



png(file.path(OUTDIR, "09_subgroup_table.png"),
    
    width=3200, height=nrow(tbl_full)*220+400, res=300)

grid.newpage(); grid.draw(g_tbl); dev.off()



pdf(file.path(OUTDIR, "09_subgroup_table.pdf"),
    
    width=12, height=max(6, nrow(tbl_full)*0.42+1.5))

grid.newpage(); grid.draw(g_tbl); dev.off()



msg("Saved: 09_subgroup_table.png/.pdf")



# =========================================================================

# Prepare fp — filtered dataset for forest plot + HR side table

# Male only for sex, same row order for both outputs

# =========================================================================



fp <- copy(tab)

fp <- fp[!is.na(HR) & is.finite(HR) & !is.na(lower) & !is.na(upper)]

fp <- fp[!(subgroup_var == "bin_sex_male" & Level == "Female")]



fp[, plot_label := paste0(Subgroup, " \u2014 ", Level)]

fp[, plot_label := factor(plot_label, levels=rev(unique(plot_label)))]



plot_h <- max(6, nrow(fp) * 0.55 + 2)



# =========================================================================

# OUTPUT 2: FOREST PLOT ONLY

# Clean, no HR table beside it — Figure C in manuscript

# =========================================================================



fp_plot <- ggplot(fp, aes(y=plot_label)) +
  
  geom_vline(xintercept=1, linetype=2, colour="grey40", linewidth=0.5) +
  
  geom_errorbarh(aes(xmin=lower, xmax=upper), height=0.25, linewidth=0.6) +
  
  geom_point(aes(x=HR), size=2.5, shape=16) +
  
  scale_x_log10(
    
    limits = c(0.1, 10),
    
    breaks = c(0.3, 1, 3),
    
    labels = c("0.3", "1", "3")
    
  ) +
  
  labs(
    
    x       = "Hazard ratio for all-cause mortality\n(Inappropriate shock vs no Inappropriate shock)",
    
    y       = NULL,
    
    title   = "Subgroup analyses",
    
    caption = paste(
      
      "Each HR reflects the within-subgroup association between inappropriate ICD shock and all-cause mortality.",
      
      "Sex: HR shown for male patients only (largest subgroup). Female subgroup result available in supplementary table.",
      
      "These are not comparisons between subgroup levels; no reference category applies.",
      
      sep="\n"
      
    )
    
  ) +
  
  theme_bw(base_size=11) +
  
  theme(
    
    plot.title       = element_text(face="bold", size=12, hjust=0.5),
    
    plot.caption     = element_text(size=8, colour="grey40", hjust=0),
    
    axis.text.y      = element_text(size=10),
    
    panel.grid.minor = element_blank(),
    
    plot.margin      = margin(5, 10, 10, 5, "mm")
    
  )



ggsave(file.path(OUTDIR, "14_forest_plot_subgroups.png"),
       
       fp_plot, width=10, height=plot_h, dpi=300)

ggsave(file.path(OUTDIR, "14_forest_plot_subgroups.pdf"),
       
       fp_plot, width=10, height=plot_h)

msg("Saved: 14_forest_plot_subgroups.png/.pdf")



# =========================================================================

# OUTPUT 3: HR SIDE TABLE

# Matches forest plot rows exactly — same order, same labels, Male only

# Columns: Subgroup / HR / 95% CI

# Sits beside forest plot in manuscript

# =========================================================================



fp_ordered <- fp[order(-as.integer(plot_label))]



tbl_hr <- data.frame(
  
  Subgroup = as.character(fp_ordered$plot_label),
  
  HR       = sprintf("%.2f", fp_ordered$HR),
  
  `95% CI` = sprintf("%.2f\u2013%.2f", fp_ordered$lower, fp_ordered$upper),
  
  check.names=FALSE, stringsAsFactors=FALSE
  
)



tbl_hr_theme <- ttheme_minimal(
  
  base_size = 11,
  
  core = list(
    
    fg_params = list(hjust=0, x=0.05),
    
    bg_params = list(
      
      fill = rep(c("white","grey97"), length.out=nrow(tbl_hr)),
      
      col  = NA
      
    )
    
  ),
  
  colhead = list(
    
    fg_params = list(fontface="bold", hjust=0, x=0.05),
    
    bg_params = list(fill="grey78", col=NA)
    
  ),
  
  padding = unit(c(4, 6), "mm")
  
)



tg_hr <- tableGrob(tbl_hr, rows=NULL, theme=tbl_hr_theme)



png(file.path(OUTDIR, "14_forest_HR_sidetable.png"),
    
    width=1800, height=nrow(tbl_hr)*220+200, res=300)

grid.newpage(); grid.draw(tg_hr); dev.off()



pdf(file.path(OUTDIR, "14_forest_HR_sidetable.pdf"),
    
    width=6, height=max(4, nrow(tbl_hr)*0.42+0.5))

grid.newpage(); grid.draw(tg_hr); dev.off()



msg("Saved: 14_forest_HR_sidetable.png/.pdf")

