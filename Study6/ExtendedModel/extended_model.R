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

library(dplyr)

#  Read pre-imputation data; keep only id + ICD_status
preimp <- read.csv("combined_BMI_outcomefiltered.csv")
 
# Load the imputed object
imp <- readRDS("mice_imputed_data.RDS")

# Make sure the columns are named consistently
# (adjust these two names if your file uses different ones)
id_col <- "ID"
icd_col <- "ICD_status"
cancer_col <- "Cancer"
stroke_col <- "Stroke_TIA"

stopifnot(id_col %in% names(preimp), icd_col %in% names(preimp))

cs_map <- preimp %>%
  transmute(
    ID = .data[[id_col]],
    
    # ICD status as before
    ICD_status = .data[[icd_col]],
    
    # Cancer: treat missing as "no" (0)
    Cancer = if_else(is.na(.data[[cancer_col]]),
                     0L,
                     as.integer(.data[[cancer_col]])),
    
    # Stroke: keep original 0/1, but also build a 3-level factor
    Stroke_TIA_raw = .data[[stroke_col]],
    
    Stroke_cat = case_when(
      is.na(.data[[stroke_col]]) ~ "missing",
      .data[[stroke_col]] == 1   ~ "yes",
      TRUE                       ~ "no"
    )
  ) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(
    Stroke_cat = factor(Stroke_cat,
                        levels = c("no", "yes", "missing"))
  )


#  Ensure imp$data has ID; if not, copy from the same preimp (row-aligned)
if (!"ID" %in% names(imp$data)) {
  stopifnot(nrow(preimp) == nrow(imp$data))
  imp$data$ID <- preimp[[id_col]]
}

# Merge into imp$data (no imputation needed)
imp$data <- imp$data %>%
  left_join(cs_map %>% select(ID, ICD_status, Cancer, Stroke_cat),
            by = "ID")


# Sanity checks
cat("N missing ICD_status after join:", sum(is.na(imp$data$ICD_status)), "\n")
cat("Distinct ICD_status values:", paste(sort(unique(imp$data$ICD_status)), collapse=", "), "\n")
cat("Cancer table:\n"); print(table(imp$data$Cancer, useNA = "ifany"))
cat("Stroke_cat table:\n"); print(table(imp$data$Stroke_cat, useNA = "ifany"))


# Quick checks
imp
names(imp$data)
imp$method[imp$method != ""]   # which vars were imputed

# Variables in the extended model
vars_extended <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute", "Cholesterol", "HDL", "LDL", "Triglycerides", "Stroke_TIA", "ICD_status"
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
      Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
      Cancer + Stroke_cat + ICD_status,
    x = TRUE, y = TRUE
  )
)

pool_fit_cs1 <- pool(fit_list_cs1)
summary(pool_fit_cs1)


# Save results 
saveRDS(pool_fit_cs1, "cox_RCS_pooled_fit_extended.RDS")


# Output the summary above
pooled_tidy <- broom::tidy(pool_fit_cs1, conf.int = TRUE)
write.csv(pooled_tidy,
          "T:/Dokumente/PROFID/Study6/cox_RCS_pooled_results_extended.csv",
          row.names = FALSE)


#  model comparison/AIC across imputations
mean(sapply(fit_list_cs1$analyses, AIC))


# Pooled AIC: average AIC across imputations
aic_vec <- sapply(fit_list_cs1$analyses, AIC)
mean(aic_vec); sd(aic_vec)


# Pooled BMI RCS curve for cause-specific model (Status_cs1)
mods <- fit_list_cs1$analyses
m    <- length(mods)

# Use the same knots/boundaries you defined earlier
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
write.csv(curve_cs1, "cox_RCS_BMI_curve_cs1_extended.csv", row.names = FALSE)

# export plot
pdf("cox_RCS_BMI_curve_cs1_extended.pdf", width = 7, height = 5)
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
          "C_index_per_imputation_cs1_extended.csv", row.names=FALSE)
write.csv(cindex_summary, "C_index_summary_cs1_extended.csv", row.names=FALSE)

print(cindex_summary)

# Save pooled summaries (coef table) and fits
saveRDS(fit_list_cs1, "fit_list_cs1_extended.RDS")             # all per-imputation fits
saveRDS(pool_fit_cs1, "pool_fit_cs1_extended.RDS")             # pooled mipo object
write.csv(broom::tidy(pool_fit_cs1, conf.int = TRUE),
          "pooled_results_cs1_extended.csv", row.names = FALSE)



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
write.csv(ph_pvals,  "PH_pvalues_per_imputation_extended.csv", row.names = FALSE)
write.csv(viol_rate, "PH_violation_rate_extended.csv",         row.names = FALSE)

# Quick on-screen view (largest violation rates first)
viol_rate

# assumes: vars_extended, BND, K3, K5, fit_list_cs1 already defined
# and imp$data$Status_cs1 exists (1 = event of interest, 0 = otherwise)

# Fit 3-knot and 5-knot models (cause-specific)
fit_k3 <- with(imp, {
  coxph(
    Surv(Survival_time, Status_cs1) ~
      ns(BMI, knots = K4, Boundary.knots = BND) +
      Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
      LVEF  + eGFR + Haemoglobin +
      ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
      Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
      Cancer + Stroke_cat + ICD_status,
    x = TRUE, y = TRUE
  )
})

fit_k5 <- with(imp, {
  coxph(
    Surv(Survival_time, Status_cs1) ~
      ns(BMI, knots = K4, Boundary.knots = BND) +
      Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
      LVEF  + eGFR + Haemoglobin +
      ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
      Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
      Cancer + Stroke_cat + ICD_status,
    x = TRUE, y = TRUE
  )
})


# AIC per imputation
aic_k3 <- sapply(fit_k3$analyses, AIC)
aic_k4 <- sapply(fit_list_cs1$analyses, AIC)  #  main 4-knot model
aic_k5 <- sapply(fit_k5$analyses, AIC)

# Summary + export
aic_summary <- data.frame(
  model    = c("K3","K4","K5"),
  mean_AIC = c(mean(aic_k3), mean(aic_k4), mean(aic_k5)),
  sd_AIC   = c(sd(aic_k3),   sd(aic_k4),   sd(aic_k5))
)

write.csv(
  data.frame(imp = seq_along(aic_k3), AIC_K3 = aic_k3, AIC_K4 = aic_k4, AIC_K5 = aic_k5),
  "AIC_per_imputation_K3_K4_K5_extended.csv",
  row.names = FALSE
)
write.csv(aic_summary, "AIC_summary_K3_K4_K5_extended.csv", row.names = FALSE)

aic_summary

# compare base and extended model visually 
# Load curves
curve_base <- read.csv("cox_RCS_BMI_curve_cs1.csv")
curve_ext  <- read.csv("cox_RCS_BMI_curve_cs1_extended.csv")

# Export high-quality PDF
pdf("BMI_curve_comparison_base_vs_extended_shaded.pdf", width = 7, height = 5)

# Base plot setup
plot(curve_base$BMI, curve_base$HR, type = "n",
     ylim = range(c(curve_base$LCL, curve_base$UCL,
                    curve_ext$LCL, curve_ext$UCL)),
     xlab = "BMI (kg/m²)",
     ylab = "Hazard ratio (ref = 25)",
     main = "Restricted cubic spline: Base vs Extended models")

# Add shaded 95% CI for base model (blue)
polygon(c(curve_base$BMI, rev(curve_base$BMI)),
        c(curve_base$LCL, rev(curve_base$UCL)),
        col = adjustcolor("darkblue", alpha.f = 0.15), border = NA)

# Add shaded 95% CI for extended model (red)
polygon(c(curve_ext$BMI, rev(curve_ext$BMI)),
        c(curve_ext$LCL, rev(curve_ext$UCL)),
        col = adjustcolor("firebrick", alpha.f = 0.15), border = NA)

# Add the mean HR curves
lines(curve_base$BMI, curve_base$HR, col = "darkblue", lwd = 2)
lines(curve_ext$BMI, curve_ext$HR, col = "firebrick", lwd = 2)

# Reference lines
abline(h = 1, lty = 3, col = "gray50")
abline(v = 25, lty = 3, col = "gray50")

# Legend
legend("topleft", legend = c("Base model", "Extended model"),
       col = c("darkblue", "firebrick"), lwd = 2, bty = "n")

dev.off()
