## sensitivity_FP_BMI_extended.R

library(survival)
library(dplyr)
library(mice)

setwd("T:/Dokumente/PROFID/Study6")

# 1. Load original mids object (UNMODIFIED)
imp0   <- readRDS("mice_imputed_data.RDS")  # fresh copy
preimp <- read.csv("combined_BMI_outcomefiltered.csv")

# Column names in preimp for extra vars
icd_col    <- "ICD_status"
cancer_col <- "Cancer"
COPD_col   <- "COPD"

# Covariates in extended model (excluding BMI, which we model as FP)
vars_extended <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute", "Cholesterol", "HDL", "LDL", "Triglycerides",
  "COPD_cat", "Cancer", "Stroke_TIA", "ICD_status"
)

## ---- Fractional polynomial models: loop over imputations safely ----

fp_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

make_single_fp_term <- function(power, var = "BMI") {
  if (power == 0) {
    sprintf("log(%s)", var)
  } else {
    sprintf("I(%s^%s)", var, power)
  }
}

fp_bmi_terms <- function(p, q = NULL, var = "BMI") {
  if (is.null(q)) {
    # FP1
    make_single_fp_term(p, var)
  } else if (p == q) {
    # FP2 with repeated power
    base <- make_single_fp_term(p, var)
    if (p == 0) {
      extra <- sprintf("I(log(%s)^2)", var)
    } else {
      extra <- sprintf("I(%s^%s * log(%s))", var, p, var)
    }
    paste(base, extra, sep = " + ")
  } else {
    # FP2 with different powers
    term1 <- make_single_fp_term(p, var)
    term2 <- make_single_fp_term(q, var)
    paste(term1, term2, sep = " + ")
  }
}

# Build FP spec list (FP1 + FP2)
fp_specs <- list()
for (p in fp_powers) {
  nm <- paste0("FP1_", p)
  fp_specs[[nm]] <- list(type = "FP1", p = p, q = NULL)
}
for (i in seq_along(fp_powers)) {
  for (j in i:length(fp_powers)) {
    p <- fp_powers[i]
    q <- fp_powers[j]
    nm <- paste0("FP2_", p, "_", q)
    fp_specs[[nm]] <- list(type = "FP2", p = p, q = q)
  }
}

spec_names <- names(fp_specs)

# Containers for fits and AICs, named by spec
fp_fits     <- vector("list", length(spec_names))
fp_aic_list <- vector("list", length(spec_names))
names(fp_fits)     <- spec_names
names(fp_aic_list) <- spec_names

# Main loop: over FP specs, then imputations
for (s in seq_along(spec_names)) {
  
  nm   <- spec_names[s]
  spec <- fp_specs[[nm]]
  
  bmi_part <- fp_bmi_terms(spec$p, spec$q, var = "BMI")
  rhs  <- paste(bmi_part, paste(vars_extended, collapse = " + "), sep = " + ")
  fmla <- as.formula(paste("Surv(Survival_time, Status_cs1) ~", rhs))
  
  cat("Fitting", nm, "...\n")
  
  fits_nm <- vector("list", imp0$m)
  aic_nm  <- numeric(imp0$m)
  
  for (i in 1:imp0$m) {
    cat("   Imputation", i, "of", imp0$m, "\n")
    flush.console()
    
    #  Completed data from untouched mids
    dat_i <- complete(imp0, i)
    
    #  Sanity check: same number of rows as preimp
    stopifnot(nrow(dat_i) == nrow(preimp))
    
    #  Add ICD_status, Cancer, COPD_cat from preimp *by row index*
    dat_i$ICD_status <- preimp[[icd_col]]
    
    dat_i$Cancer <- ifelse(
      is.na(preimp[[cancer_col]]),
      0L,
      as.integer(preimp[[cancer_col]])
    )
    
    dat_i$COPD_cat <- factor(
      dplyr::case_when(
        is.na(preimp[[COPD_col]])      ~ "missing",
        preimp[[COPD_col]] == 1        ~ "Yes",
        TRUE                           ~ "No"
      ),
      levels = c("No", "Yes", "missing")
    )
    
    # Cause-specific status
    dat_i$Status_cs1 <- ifelse(dat_i$Status == 1, 1L, 0L)
    
    #  Fit model in this imputation
    fits_nm[[i]] <- coxph(fmla, data = dat_i, x = TRUE, y = TRUE)
    aic_nm[i]    <- AIC(fits_nm[[i]])
  }
  
  # store results!!!
  fp_fits[[nm]]     <- fits_nm
  fp_aic_list[[nm]] <- aic_nm
}

## ---- Summarise AIC ----

# which specs actually have valid AIC vectors?
valid_specs <- spec_names[
  sapply(fp_aic_list, function(x) is.numeric(x) && length(x) > 0)
]

if (length(valid_specs) == 0) {
  stop("No valid AIC vectors in fp_aic_list – check that the FP loop ran successfully.")
}

fp_aic_summary <- do.call(
  rbind,
  lapply(valid_specs, function(nm) {
    vec <- fp_aic_list[[nm]]
    data.frame(
      spec      = nm,
      type      = ifelse(substr(nm, 1, 3) == "FP1", "FP1", "FP2"),
      mean_AIC  = mean(vec, na.rm = TRUE),
      sd_AIC    = sd(vec,   na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)

fp_aic_summary <- fp_aic_summary[order(fp_aic_summary$mean_AIC), ]

write.csv(fp_aic_summary,
          "FP_AIC_summary_extended.csv",
          row.names = FALSE)

# Per-imputation AIC wide table (for valid specs only)
aic_mat <- sapply(valid_specs, function(nm) fp_aic_list[[nm]])
aic_df  <- data.frame(
  imputation = seq_len(nrow(aic_mat)),
  aic_mat
)
write.csv(aic_df,
          "FP_AIC_per_imputation_extended.csv",
          row.names = FALSE)

print(head(fp_aic_summary, 10))

# Best specs
best_overall <- fp_aic_summary$spec[which.min(fp_aic_summary$mean_AIC)]
best_fp1     <- fp_aic_summary$spec[fp_aic_summary$type == "FP1"][
  which.min(fp_aic_summary$mean_AIC[fp_aic_summary$type == "FP1"])
]
best_fp2     <- fp_aic_summary$spec[fp_aic_summary$type == "FP2"][
  which.min(fp_aic_summary$mean_AIC[fp_aic_summary$type == "FP2"])
]

cat("Best overall FP spec:", best_overall, "\n")
cat("Best FP1 spec:",       best_fp1,     "\n")
cat("Best FP2 spec:",       best_fp2,     "\n")

# Save fits for best overall model
saveRDS(fp_fits[[best_overall]],
        paste0("fit_list_cs1_extended_", best_overall, ".RDS"))
