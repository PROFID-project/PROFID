install.packages("dplyr")
install.packages("mice")

library(survival)
library(dplyr)
library(splines)
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

## ---- Fractional polynomial models: FP1 + diagonal FP2 only ----

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
    # FP2 with repeated power (diagonal)
    base <- make_single_fp_term(p, var)
    if (p == 0) {
      extra <- sprintf("I(log(%s)^2)", var)
    } else {
      extra <- sprintf("I(%s^%s * log(%s))", var, p, var)
    }
    paste(base, extra, sep = " + ")
  } else {
    stop("This script is set up for diagonal FP2 only (p = q).")
  }
}

# Build FP spec list
fp_specs <- list()

# FP1: all powers
for (p in fp_powers) {
  nm <- paste0("FP1_", p)
  fp_specs[[nm]] <- list(type = "FP1", p = p, q = NULL)
}

# FP2: only repeated powers (p, p)  -> "diagonal" FP2s
for (p in fp_powers) {
  nm <- paste0("FP2_", p, "_", p)
  fp_specs[[nm]] <- list(type = "FP2", p = p, q = p)
}

spec_names <- names(fp_specs)

# (keep the rest of the script – the big loop, AIC summary, etc. – exactly as in the fixed version I sent)


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
s
if (length(valid_specs) == 0) {
  stop("No valid AIC vectors in fp_aic_list – check that the FP loop ran successfully.")
}

fp_aic_summary <- do.call(
  rbind,
  lapply(valid_specs, function(nm) {
    vec <- fp_aic_list[[nm]]
    data.frame(
      spec      = nm,
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


#  Load pre-imputation data (for BMI distribution)
preimp <- read.csv("combined_BMI_outcomefiltered.csv")

#  Load FP AIC summary and find best spec
fp_aic_summary <- read.csv("FP_AIC_summary_extended.csv")

best_overall <- fp_aic_summary$spec[which.min(fp_aic_summary$mean_AIC)]
cat("Best overall FP spec (from CSV):", best_overall, "\n")

# Load the saved fits for that spec
fit_file  <- paste0("fit_list_cs1_extended_", best_overall, ".RDS")
mods_best <- readRDS(fit_file)   # list of coxph objects, one per imputation
m_best    <- length(mods_best)

#  Parse the FP spec name, e.g. "FP1_0" or "FP2_-1_-1"
make_fp_numeric_terms <- function(bmi, p, q = NULL) {
  if (is.null(q)) {
    # FP1
    if (p == 0) {
      return(log(bmi))
    } else {
      return(bmi^p)
    }
  } else {
    # FP2 (diagonal p = q)
    if (p == 0) {
      return(c(log(bmi), log(bmi)^2))
    } else {
      return(c(bmi^p, bmi^p * log(bmi)))
    }
  }
}

parts <- strsplit(best_overall, "_")[[1]]

if (parts[1] == "FP1") {
  p_best <- as.numeric(parts[2])
  q_best <- NULL
} else if (parts[1] == "FP2") {
  p_best <- as.numeric(parts[2])
  q_best <- as.numeric(parts[3])
} else {
  stop("Unrecognised best_overall spec format: ", best_overall)
}

cat("Using powers p =", p_best,
    ifelse(is.null(q_best), "", paste("and q =", q_best)),
    "\n")

#  BMI grid (e.g. 5th–95th percentile range)
q_bmi    <- quantile(preimp$BMI, probs = c(0.05, 0.95), na.rm = TRUE)
bmi_grid <- seq(q_bmi[1], q_bmi[2], length.out = 300)

# Build FP term matrix for each BMI
n_terms <- if (is.null(q_best)) 1 else 2

Lmat <- t(vapply(
  bmi_grid,
  function(b) as.numeric(make_fp_numeric_terms(b, p_best, q_best)),
  numeric(n_terms)
))

# Reference at BMI = 25
ref_terms <- as.numeric(make_fp_numeric_terms(25, p_best, q_best))

# Contrast matrix vs ref BMI = 25
Lcontr <- sweep(Lmat, 2, ref_terms, FUN = "-")

#  Identify BMI-related coefficient positions in the model
coef_names <- names(coef(mods_best[[1]]))
sidx <- grep("BMI", coef_names)

if (length(sidx) != ncol(Lcontr)) {
  stop("Number of BMI FP coefficients (", length(sidx),
       ") does not match contrast matrix columns (", ncol(Lcontr), ").")
}

#  Rubin-style pooling of log-HR for each BMI value
pool_contrast <- function(Lrow) {
  qi <- vapply(mods_best, function(fm) {
    as.numeric(Lrow %*% coef(fm)[sidx])
  }, numeric(1))
  
  ui <- vapply(mods_best, function(fm) {
    V <- vcov(fm)[sidx, sidx, drop = FALSE]
    as.numeric(Lrow %*% V %*% Lrow)
  }, numeric(1))
  
  qbar <- mean(qi)
  ubar <- mean(ui)
  bvar <- if (length(qi) > 1) var(qi) else 0
  Tvar <- ubar + (1 + 1/m_best) * bvar
  
  c(logHR = qbar, se = sqrt(Tvar))
}

res   <- t(apply(Lcontr, 1, pool_contrast))
HR    <- exp(res[, "logHR"])
LCL   <- exp(res[, "logHR"] - 1.96 * res[, "se"])
UCL   <- exp(res[, "logHR"] + 1.96 * res[, "se"])

fp_curve <- data.frame(
  BMI = bmi_grid,
  HR  = HR,
  LCL = LCL,
  UCL = UCL
)

# 8. Save curve and plot
curve_csv_name <- paste0("FP_BMI_curve_", best_overall, "_extended.csv")
curve_pdf_name <- paste0("FP_BMI_curve_", best_overall, "_extended.pdf")

write.csv(fp_curve, curve_csv_name, row.names = FALSE)

pdf(curve_pdf_name, width = 7, height = 5)
plot(fp_curve$BMI, fp_curve$HR, type = "l",
     xlab = "BMI (kg/m²)",
     ylab = "Hazard ratio (ref = BMI 25)",
     ylim = range(c(fp_curve$LCL, fp_curve$UCL)),
     lwd = 2,
     main = paste("Fractional polynomial model:", best_overall))
lines(fp_curve$BMI, fp_curve$LCL, lty = 2)
lines(fp_curve$BMI, fp_curve$UCL, lty = 2)
abline(h = 1,  lty = 3, col = "gray60")
abline(v = 25, lty = 3, col = "gray60")
dev.off()

cat("FP BMI curve saved to:\n  ", curve_csv_name, "\n  ", curve_pdf_name, "\n")

