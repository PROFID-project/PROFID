library(survival)
library(dplyr)
library(mice)

setwd("T:/Dokumente/PROFID/Study6")

# ---- Load imputed object (use the same one you used for the extended model) ----
# If your extended model uses imp2, use imp2 here.
imp2 <- readRDS("mice_imputed_data_extended.RDS")

# Pre-imputation/raw data (only needed if you have variables not in the mids object)
preimp <- read.csv("combined_BMI_outcomefiltered.csv")

# If ICD_status isn't in the imputed datasets and you need to attach it:
icd_col <- "ICD_status"

# ---- Updated extended covariates (EXCLUDING BMI, which is FP) ----
# (Removed COPD_cat + Cancer to match your updated extended model)
vars_extended <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute", "Cholesterol", "HDL", "LDL", "Triglycerides",
  "Stroke_TIA", "ICD_status"
)

# ---- 90-month horizon variables ----
horizon <- 90

make_90mo <- function(d) {
  tt <- d$Survival_time
  ss <- d$Status
  ss <- suppressWarnings(as.integer(as.character(ss)))
  
  d$Survival_time_h <- pmin(tt, horizon)
  d$Status_cs1_h    <- ifelse(!is.na(ss) & ss == 1L & tt <= horizon, 1L, 0L)
  d
}

# ---- Fractional polynomial models: FP1 + diagonal FP2 only ----
fp_powers <- c(-2, -1, -0.5, 0, 0.5, 1, 2, 3)

make_single_fp_term <- function(power, var = "BMI") {
  if (power == 0) sprintf("log(%s)", var) else sprintf("I(%s^%s)", var, power)
}

fp_bmi_terms <- function(p, q = NULL, var = "BMI") {
  if (is.null(q)) {
    make_single_fp_term(p, var)                     # FP1
  } else if (p == q) {
    base <- make_single_fp_term(p, var)             # FP2 diagonal
    extra <- if (p == 0) sprintf("I(log(%s)^2)", var) else sprintf("I(%s^%s * log(%s))", var, p, var)
    paste(base, extra, sep = " + ")
  } else {
    stop("This script is set up for diagonal FP2 only (p = q).")
  }
}

# Build FP spec list
fp_specs <- list()
for (p in fp_powers) fp_specs[[paste0("FP1_", p)]] <- list(type = "FP1", p = p, q = NULL)
for (p in fp_powers) fp_specs[[paste0("FP2_", p, "_", p)]] <- list(type = "FP2", p = p, q = p)

spec_names <- names(fp_specs)

# Containers
fp_fits     <- vector("list", length(spec_names)); names(fp_fits) <- spec_names
fp_aic_list <- vector("list", length(spec_names)); names(fp_aic_list) <- spec_names

# ---- Main loop over FP specs and imputations ----
for (nm in spec_names) {
  
  spec <- fp_specs[[nm]]
  
  bmi_part <- fp_bmi_terms(spec$p, spec$q, var = "BMI")
  rhs  <- paste(bmi_part, paste(vars_extended, collapse = " + "), sep = " + ")
  fmla <- as.formula(paste("Surv(Survival_time_h, Status_cs1_h) ~", rhs))
  
  cat("Fitting", nm, "...\n")
  
  fits_nm <- vector("list", imp2$m)
  aic_nm  <- numeric(imp2$m)
  
  for (i in 1:imp2$m) {
    cat("   Imputation", i, "of", imp2$m, "\n")
    flush.console()
    
    dat_i <- complete(imp2, i)
    
    # If ICD_status needs to be brought in from preimp:
    if (!("ICD_status" %in% names(dat_i)) && icd_col %in% names(preimp)) {
      stopifnot(nrow(dat_i) == nrow(preimp))
      dat_i$ICD_status <- preimp[[icd_col]]
    }
    
    # Create 90-month horizon outcome vars
    dat_i <- make_90mo(dat_i)
    
    # Fit and store AIC
    fits_nm[[i]] <- coxph(fmla, data = dat_i, x = TRUE, y = TRUE)
    aic_nm[i]    <- AIC(fits_nm[[i]])
  }
  
  fp_fits[[nm]]     <- fits_nm
  fp_aic_list[[nm]] <- aic_nm
}

# ---- Summarise AIC ----
valid_specs <- spec_names[sapply(fp_aic_list, function(x) is.numeric(x) && length(x) > 0)]
if (length(valid_specs) == 0) stop("No valid AIC vectors in fp_aic_list.")

fp_aic_summary <- do.call(
  rbind,
  lapply(valid_specs, function(nm) {
    vec <- fp_aic_list[[nm]]
    data.frame(
      spec     = nm,
      type     = ifelse(grepl("^FP1_", nm), "FP1", "FP2"),
      mean_AIC = mean(vec, na.rm = TRUE),
      sd_AIC   = sd(vec,   na.rm = TRUE),
      stringsAsFactors = FALSE
    )
  })
)
fp_aic_summary <- fp_aic_summary[order(fp_aic_summary$mean_AIC), ]

write.csv(fp_aic_summary, "FP_AIC_summary_extended_90mo.csv", row.names = FALSE)

# Best specs
best_overall <- fp_aic_summary$spec[which.min(fp_aic_summary$mean_AIC)]
best_fp1 <- fp_aic_summary$spec[fp_aic_summary$type == "FP1"][which.min(fp_aic_summary$mean_AIC[fp_aic_summary$type == "FP1"])]
best_fp2 <- fp_aic_summary$spec[fp_aic_summary$type == "FP2"][which.min(fp_aic_summary$mean_AIC[fp_aic_summary$type == "FP2"])]

cat("Best overall FP spec:", best_overall, "\n")
cat("Best FP1 spec:", best_fp1, "\n")
cat("Best FP2 spec:", best_fp2, "\n")

# Save fits for best overall model
saveRDS(fp_fits[[best_overall]], paste0("fit_list_cs1_extended_90mo_", best_overall, ".RDS"))

# ---- Build pooled FP BMI curve for best model (ref BMI=25) ----
mods_best <- fp_fits[[best_overall]]
m_best <- length(mods_best)

make_fp_numeric_terms <- function(bmi, p, q = NULL) {
  if (is.null(q)) {
    if (p == 0) log(bmi) else bmi^p
  } else {
    if (p == 0) c(log(bmi), log(bmi)^2) else c(bmi^p, bmi^p * log(bmi))
  }
}

parts <- strsplit(best_overall, "_")[[1]]
if (parts[1] == "FP1") {
  p_best <- as.numeric(parts[2]); q_best <- NULL
} else {
  p_best <- as.numeric(parts[2]); q_best <- as.numeric(parts[3])
}

q_bmi <- quantile(preimp$BMI, probs = c(0.05, 0.95), na.rm = TRUE)
bmi_grid <- seq(q_bmi[1], q_bmi[2], length.out = 300)

n_terms <- if (is.null(q_best)) 1 else 2
Lmat <- t(vapply(bmi_grid, function(b) as.numeric(make_fp_numeric_terms(b, p_best, q_best)), numeric(n_terms)))

ref_terms <- as.numeric(make_fp_numeric_terms(25, p_best, q_best))
Lcontr <- sweep(Lmat, 2, ref_terms, FUN = "-")

coef_names <- names(coef(mods_best[[1]]))
sidx <- grep("BMI", coef_names)

if (length(sidx) != ncol(Lcontr)) {
  stop("Number of BMI FP coefficients (", length(sidx),
       ") does not match contrast columns (", ncol(Lcontr), ").")
}

pool_contrast <- function(Lrow) {
  qi <- vapply(mods_best, function(fm) as.numeric(Lrow %*% coef(fm)[sidx]), numeric(1))
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

res <- t(apply(Lcontr, 1, pool_contrast))
HR  <- exp(res[, "logHR"])
LCL <- exp(res[, "logHR"] - 1.96 * res[, "se"])
UCL <- exp(res[, "logHR"] + 1.96 * res[, "se"])

fp_curve <- data.frame(BMI = bmi_grid, HR = HR, LCL = LCL, UCL = UCL)

curve_csv <- paste0("FP_BMI_curve_", best_overall, "_extended_90mo.csv")
curve_pdf <- paste0("FP_BMI_curve_", best_overall, "_extended_90mo.pdf")

write.csv(fp_curve, curve_csv, row.names = FALSE)

pdf(curve_pdf, width = 7, height = 5)
plot(fp_curve$BMI, fp_curve$HR, type = "l", lwd = 2,
     xlab = "BMI (kg/m²)",
     ylab = "Hazard ratio (ref = BMI 25)",
     ylim = range(c(fp_curve$LCL, fp_curve$UCL)),
     main = paste("Fractional polynomial model:", best_overall))
lines(fp_curve$BMI, fp_curve$LCL, lty = 2)
lines(fp_curve$BMI, fp_curve$UCL, lty = 2)
abline(h = 1,  lty = 3, col = "gray60")
abline(v = 25, lty = 3, col = "gray60")
dev.off()

cat("Saved:\n", curve_csv, "\n", curve_pdf, "\n")
