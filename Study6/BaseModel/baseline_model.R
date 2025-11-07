install.packages("survival")
install.packages("survminer")
install.packages("dplyr")
install.packages("svglite")
install.packages("rms")
install.packages("mice")

library(survival)
library(survminer)
library(dplyr)
library(rms)
library(mice)
library(splines)
library(dplyr)
library(broom)

setwd("T:/Dokumente/PROFID/Study6")

# Load the imputed object
imp <- readRDS("mice_imputed_data.RDS")

# Quick checks
imp
names(imp$data)
imp$method[imp$method != ""]   # which vars were imputed

# Variables in the base model
vars_base <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute"
)

# Fixed knots from observed BMI in the original data inside imp$data 
q   <- quantile(imp$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])              # boundary knots 5th & 95th
K4  <- as.numeric(q[c(2,3,4,5)])           # inner knots 10,35,65,90

# Sensitivity options (3rd and 5th knots)
K3 <- as.numeric(quantile(imp$data$BMI, probs = c(.10,.50,.90), na.rm = TRUE))
K5 <- as.numeric(quantile(imp$data$BMI, probs = c(.05,.275,.50,.725,.95), na.rm = TRUE))[2:4]

# Recode Status: event of interest = 1, everything else (0 or 2) = 0
imp$data$Status_cs1 <- ifelse(imp$data$Status == 1, 1L, 0L)

# Quick sanity check
table(imp$data$Status, imp$data$Status_cs1, useNA = "ifany")

# Refit cause-specific Cox (y=TRUE so we can compute C-index reliably)
fit_list_cs1 <- with(
  imp,
  coxph(
    Surv(Survival_time, Status_cs1) ~
      ns(BMI, knots = K4, Boundary.knots = BND) +
      Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
      LVEF  + eGFR + Haemoglobin +
      ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
      Revascularisation_acute,
    x = TRUE, y = TRUE
  )
)

pool_fit_cs1 <- pool(fit_list_cs1)
summary(pool_fit_cs1)


# Save results 
saveRDS(pool_fit_cs1, "cox_RCS_pooled_fit.RDS")


# Output the summary above
pooled_tidy <- broom::tidy(pool_fit_cs1, conf.int = TRUE)
write.csv(pooled_tidy,
          "T:/Dokumente/PROFID/Study6/cox_RCS_pooled_results.csv",
          row.names = FALSE)


#  model comparison/AIC across imputations
mean(sapply(fit_list_cs1$analyses, AIC))
# Pooled AIC: average AIC across imputations
aic_vec <- sapply(fit_list_cs1$analyses, AIC)
mean(aic_vec); sd(aic_vec)


# Pooled BMI RCS curve for cause-specific model (Status_cs1)
mods <- fit_list_cs1$analyses
m    <- length(mods)

# Use the same knots/boundaries defined earlier
bmi_min <- max(min(imp$data$BMI, na.rm = TRUE), BND[1])
bmi_max <- min(max(imp$data$BMI, na.rm = TRUE), BND[2])
bmi_grid <- seq(bmi_min, bmi_max, length.out = 400)

Xbmi <- ns(bmi_grid, knots = K4, Boundary.knots = BND)
Xref <- ns(25,        knots = K4, Boundary.knots = BND)
Lmat <- Xbmi - matrix(rep(Xref, each = nrow(Xbmi)), ncol = ncol(Xbmi))  # contrasts vs BMI=25

# indices of BMI spline terms (same across imputations)
coef_names <- names(coef(mods[[1]]))
sidx <- grep("^ns\\(BMI", coef_names)

# pool a single contrast row using Rubin's rules
pool_contrast <- function(Lrow){
  qi <- vapply(mods, function(fm) as.numeric(Lrow %*% coef(fm)[sidx]), numeric(1))
  ui <- vapply(mods, function(fm) {
    as.numeric(Lrow %*% vcov(fm)[sidx, sidx, drop=FALSE] %*% Lrow)
  }, numeric(1))
  qbar <- mean(qi); ubar <- mean(ui); bvar <- var(qi)
  Tvar <- ubar + (1 + 1/m) * bvar
  c(logHR = qbar, se = sqrt(Tvar))
}

res   <- t(apply(Lmat, 1, pool_contrast))
HR    <- exp(res[, "logHR"])
LCL   <- exp(res[, "logHR"] - 1.96*res[, "se"])
UCL   <- exp(res[, "logHR"] + 1.96*res[, "se"])

# export CSV
curve_cs1 <- data.frame(BMI = bmi_grid, HR = HR, LCL = LCL, UCL = UCL)
write.csv(curve_cs1, "cox_RCS_BMI_curve_cs1.csv", row.names = FALSE)

# export plot
pdf("cox_RCS_BMI_curve_cs1.pdf", width = 7, height = 5)
plot(bmi_grid, HR, type = "l", xlab = "BMI (kg/m²)", ylab = "Hazard ratio (ref = 25)",
     ylim = range(c(LCL, UCL)), lwd = 2)
lines(bmi_grid, LCL, lty = 2); lines(bmi_grid, UCL, lty = 2)
abline(h = 1, lty = 3); abline(v = 25, lty = 3)
dev.off()

# Harrell's C per imputation
cindex_per_imp <- sapply(fit_list_cs1$analyses, function(fm){
  y  <- fm$y
  lp <- predict(fm, type="lp")
  survival::survConcordance(Surv(y[,1], y[,2]) ~ lp)$concordance
})

cindex_summary <- data.frame(
  Mean_Cindex   = mean(cindex_per_imp),
  SD_Cindex     = sd(cindex_per_imp),
  N_Imputations = length(cindex_per_imp)
)

write.csv(data.frame(Imputation=seq_along(cindex_per_imp), C_index=cindex_per_imp),
          "C_index_per_imputation_cs1.csv", row.names=FALSE)
write.csv(cindex_summary, "C_index_summary_cs1.csv", row.names=FALSE)

print(cindex_summary)

# Save pooled summaries (coef table) and fits
saveRDS(fit_list_cs1, "fit_list_cs1.RDS")             # all per-imputation fits
saveRDS(pool_fit_cs1, "pool_fit_cs1.RDS")             # pooled mipo object
write.csv(broom::tidy(pool_fit_cs1, conf.int = TRUE),
          "pooled_results_cs1.csv", row.names = FALSE)



# Proportional Hazards tests per imputation
ph_tbl <- lapply(fit_list_cs1$analyses, function(fm) {
  as.data.frame(cox.zph(fm)$table)   # rows = covariates + "GLOBAL"
})

# Collect p-values into a single wide table
ph_pvals <- Reduce(function(a, b) cbind(a["p"], b["p"]),
                   ph_tbl)
colnames(ph_pvals) <- paste0("imp", seq_len(ncol(ph_pvals)))
ph_pvals$term <- rownames(ph_tbl[[1]])
ph_pvals <- ph_pvals[ph_pvals$term != "GLOBAL", ]     # drop global row for per-covariate summary
row.names(ph_pvals) <- NULL

# Proportion of imputations with p < 0.05 (violation rate) for each covariate
viol_rate <- data.frame(
  term = ph_pvals$term,
  prop_p_lt_0.05 = apply(ph_pvals[ , startsWith(names(ph_pvals), "imp")], 1,
                         function(x) mean(x < 0.05, na.rm = TRUE))
)
viol_rate <- viol_rate[order(-viol_rate$prop_p_lt_0.05), ]

# Export p-values and summary
write.csv(ph_pvals,  "PH_pvalues_per_imputation.csv", row.names = FALSE)
write.csv(viol_rate, "PH_violation_rate.csv",         row.names = FALSE)

# Quick on-screen view (largest violation rates first)
viol_rate



# assumes: vars_base, BND, K3, K5, fit_list_cs1 already defined
# and imp$data$Status_cs1 exists (1 = event of interest, 0 = otherwise)


# Fit 3-knot and 5-knot models (cause-specific)
fit_k3 <- with(imp, {
  coxph(Surv(Survival_time, Status_cs1) ~ 
          ns(BMI, knots = K3, Boundary.knots = BND) + 
          Age + Sex + Diabetes + Hypertension + Smoking + MI_history + 
          LVEF + eGFR + Haemoglobin + ACE_inhibitor_ARB + 
          Beta_blockers + Lipid_lowering + Revascularisation_acute,
        x = TRUE, y = TRUE)
})

fit_k5 <- with(imp, {
  coxph(Surv(Survival_time, Status_cs1) ~ 
          ns(BMI, knots = K5, Boundary.knots = BND) + 
          Age + Sex + Diabetes + Hypertension + Smoking + MI_history + 
          LVEF + eGFR + Haemoglobin + ACE_inhibitor_ARB + 
          Beta_blockers + Lipid_lowering + Revascularisation_acute,
        x = TRUE, y = TRUE)
})


# AIC per imputation
aic_k3 <- sapply(fit_k3$analyses, AIC)
aic_k4 <- sapply(fit_list_cs1$analyses, AIC)  # your main 4-knot model
aic_k5 <- sapply(fit_k5$analyses, AIC)

# Summary + export
aic_summary <- data.frame(
  model    = c("K3","K4","K5"),
  mean_AIC = c(mean(aic_k3), mean(aic_k4), mean(aic_k5)),
  sd_AIC   = c(sd(aic_k3),   sd(aic_k4),   sd(aic_k5))
)

write.csv(
  data.frame(imp = seq_along(aic_k3), AIC_K3 = aic_k3, AIC_K4 = aic_k4, AIC_K5 = aic_k5),
  "AIC_per_imputation_K3_K4_K5.csv",
  row.names = FALSE
)
write.csv(aic_summary, "AIC_summary_K3_K4_K5.csv", row.names = FALSE)

aic_summary

# backward selection with p values on 1 imputed dataset 
# Build helpers
build_formula <- function(vars){
  as.formula(paste0(
    "Surv(Survival_time, Status_cs1) ~ ",
    "ns(BMI, knots = K4, Boundary.knots = BND) + ",
    paste(vars, collapse = " + ")
  ))
}

# Work on the first completed dataset
d1 <- mice::complete(imp, 1)

# Full model (base set)
f_full <- build_formula(vars_base)
fit1_full <- coxph(f_full, data = d1, x = TRUE, y = TRUE)

# Report variables with p<0.10 in the full model (per term, LR test)
full_drop <- drop1(fit1_full, test = "Chisq")
full_keep_010 <- rownames(full_drop)[which(full_drop$"Pr(>Chi)" < 0.10)]
full_keep_010 <- setdiff(full_keep_010, rownames(full_drop)[grepl("^ns\\(BMI", rownames(full_drop))])
full_keep_010
# (this is for reporting, start backward from the full model)

# Backward elimination to p<0.05 retention
keep <- vars_base
repeat {
  f <- build_formula(keep)
  fit <- coxph(f, data = d1, x = TRUE, y = TRUE)
  dr <- drop1(fit, test = "Chisq")
  # remove BMI spline row from consideration
  dr <- dr[!grepl("^ns\\(BMI", rownames(dr)), , drop = FALSE]
  # if nothing to drop, or all p<0.05, stop
  if (nrow(dr) == 0 || all(is.na(dr$"Pr(>Chi)") | dr$"Pr(>Chi)" < 0.05)) break
  # drop the worst (largest p)
  worst <- rownames(dr)[which.max(dr$"Pr(>Chi)")]
  if (max(dr$"Pr(>Chi)", na.rm = TRUE) < 0.05) break
  keep <- setdiff(keep, worst)
  if (length(keep) == 0) break
}

keep

# make sure Status_cs1 is inside imp$data
if (!"Status_cs1" %in% names(imp$data)) {
  imp$data$Status_cs1 <- ifelse(imp$data$Status == 1, 1L, 0L)
}

# safety: Survival_time must be present too
stopifnot("Survival_time" %in% names(imp$data))

# refit pooled simplified model with the 'keep' variables
fit_list_cs1_simpl <- with(imp, {
  rhs  <- paste0("ns(BMI, knots = K4, Boundary.knots = BND) + ",
                 paste(keep, collapse = " + "))
  form <- as.formula(paste("Surv(Survival_time, Status_cs1) ~", rhs))
  coxph(form, x = TRUE, y = TRUE)
})

pool_fit_cs1_simpl <- pool(fit_list_cs1_simpl)
summary(pool_fit_cs1_simpl)


# Export pooled coefficients for simplified model
write.csv(broom::tidy(pool_fit_cs1_simpl, conf.int = TRUE),
          "pooled_results_base_simplified.csv", row.names = FALSE)

# Compare AIC (mean across imputations)
aic_full  <- sapply(fit_list_cs1$analyses,       AIC)
aic_simpl <- sapply(fit_list_cs1_simpl$analyses, AIC)
data.frame(model = c("Full base", "Simplified base"),
           mean_AIC = c(mean(aic_full), mean(aic_simpl)),
           sd_AIC   = c(sd(aic_full),   sd(aic_simpl)))
# lower mean_AIC is better

# Compare discrimination
cidx_full <- sapply(fit_list_cs1$analyses, function(fm){
  y <- fm$y; lp <- predict(fm, type="lp")
  survival::survConcordance(Surv(y[,1], y[,2]) ~ lp)$concordance
})
cidx_simpl <- sapply(fit_list_cs1_simpl$analyses, function(fm){
  y <- fm$y; lp <- predict(fm, type="lp")
  survival::survConcordance(Surv(y[,1], y[,2]) ~ lp)$concordance
})
data.frame(model=c("Full base","Simplified base"),
           mean_C = c(mean(cidx_full), mean(cidx_simpl)),
           sd_C   = c(sd(cidx_full),   sd(cidx_simpl)))


# Extract and average AIC per model across imputations
aic_full <- sapply(fit_list_cs1$analyses, AIC)
aic_simpl <- sapply(fit_list_cs1_simpl$analyses, AIC)

# Summarise results
aic_summary <- data.frame(
  model = c("Full base", "Simplified base"),
  mean_AIC = c(mean(aic_full), mean(aic_simpl)),
  sd_AIC   = c(sd(aic_full), sd(aic_simpl))
)

print(aic_summary)
