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

#  Read pre-imputation data; keep only id + ICD_status
preimp <- read.csv("combined_BMI_outcomefiltered.csv")
 
# Load the imputed object
imp <- readRDS("mice_imputed_data.RDS")

# Make sure the columns are named consistently
id_col <- "ID"
icd_col <- "ICD_status"
cancer_col <- "Cancer"
COPD_col <- "COPD"
stopifnot(
  id_col %in% names(preimp),
  icd_col %in% names(preimp),
  cancer_col %in% names(preimp),
  COPD_col %in% names(preimp)
)

cs_map <- preimp %>%
  transmute(
    ID = .data[[id_col]],
    
    ICD_status = .data[[icd_col]],
    
    # Cancer: treat missing as "no" (0)
    Cancer = case_when(
      .data[[cancer_col]] %in% c("Yes", 1L) ~ 1L,
      is.na(.data[[cancer_col]])            ~ 0L,
      TRUE                                  ~ 0L
    ),
    
    # COPD: original is No/Yes (or 0/1); make 3-level factor
    COPD_raw = .data[[COPD_col]],
    COPD_cat = case_when(
      is.na(.data[[COPD_col]])            ~ "Missing",
      .data[[COPD_col]] %in% c("Yes", 1L) ~ "Yes",
      TRUE                                ~ "No"
    )
  ) %>%
  distinct(ID, .keep_all = TRUE) %>%
  mutate(
    COPD_cat = factor(COPD_cat,
                      levels = c("No", "Yes", "Missing"))
  )

imp_long <- complete(imp, "long", include = TRUE)  # has .imp and .id

## Make sure ID is present in imp_long; if not, map by row order
if (!"ID" %in% names(imp_long)) {
  stopifnot(nrow(preimp) == length(unique(imp_long$.id)))
  imp_long <- imp_long %>%
    group_by(.imp) %>%
    arrange(.id) %>%
    mutate(ID = preimp[[id_col]]) %>%
    ungroup()
}

## Merge in ICD_status, Cancer, COPD_cat
imp_long <- imp_long %>%
  left_join(cs_map %>% select(ID, ICD_status, Cancer, COPD_cat),
            by = "ID") %>%
  mutate(
    # Recode Status_cs1 in the long data
    Status_cs1 = ifelse(Status == 1, 1L, 0L)
  )

## Back to a well-formed mids object
imp2 <- as.mids(imp_long)



# Sanity checks
cat("N missing ICD_status after join:", sum(is.na(imp2$data$ICD_status)), "\n")
cat("Distinct ICD_status values:", paste(sort(unique(imp2$data$ICD_status)), collapse=", "), "\n")
cat("Cancer table:\n"); print(table(imp2$data$Cancer, useNA = "ifany"))
cat("COPD_cat table:\n"); print(table(imp2$data$COPD_cat, useNA = "ifany"))

saveRDS(imp2, "mice_imputed_data_extended.RDS")


# Quick checks
imp2
names(imp2$data)
imp2$method[imp2$method != ""]   # which vars were imputed

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
q   <- quantile(imp2$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])              # boundary knots 5th & 95th
K4  <- as.numeric(q[c(2,3,4,5)])           # inner knots 10,35,65,90

# Sensitivity options (3rd and 5th knots)
K3 <- as.numeric(quantile(imp2$data$BMI, probs = c(.10,.50,.90), na.rm = TRUE))
K5 <- as.numeric(quantile(imp2$data$BMI, probs = c(.05,.275,.50,.725,.95), na.rm = TRUE))[2:4]

# Recode Status: event of interest = 1, everything else (0 or 2) = 0
imp2$data$Status_cs1 <- ifelse(imp2$data$Status == 1, 1L, 0L)

# Quick sanity check
table(imp2$data$Status, imp2$data$Status_cs1, useNA = "ifany")

horizon <- 90

imp2$data <- imp2$data %>%
  mutate(
    Survival_time_h = pmin(Survival_time, horizon),
    Status_cs1_h    = ifelse(Survival_time <= horizon & Status == 1, 1L, 0L)
  )


# Refit cause-specific Cox (y=TRUE so we can compute C-index reliably)
fit_list_cs1 <- with(
  imp2,
  coxph(
    Surv(Survival_time_h, Status_cs1_h) ~
      ns(BMI, knots = K4, Boundary.knots = BND) +
      Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
      LVEF  + eGFR + Haemoglobin +
      ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
      Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
       + Stroke_TIA + ICD_status,
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

# Use the same knots/boundaries defined earlier
bmi_min <- max(min(imp2$data$BMI, na.rm = TRUE), BND[1])
bmi_max <- min(max(imp2$data$BMI, na.rm = TRUE), BND[2])
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
plot(bmi_grid, HR, type = "l", xlab = "BMI (kg/m²)", ylab = "Hazard ratio (reference BMI 25 kg/m²)",
     ylim = range(c(LCL, UCL)), lwd = 2, col = "black")
     
     lines(bmi_grid, LCL,
           lty = 2,
           lwd = 1.5,
           col = "grey60")
     
     lines(bmi_grid, UCL,
           lty = 2,
           lwd = 1.5,
           col = "grey60")

#lines(bmi_grid, LCL, lty = 2); lines(bmi_grid, UCL, lty = 2)
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
colnames(ph_pvals) <- paste0("imp2", seq_len(ncol(ph_pvals)))
ph_pvals$term <- rownames(ph_tbl[[1]])
ph_pvals <- ph_pvals[ph_pvals$term != "GLOBAL", ]     # drop global row for per-covariate summary
row.names(ph_pvals) <- NULL

# Proportion of imputations with p < 0.05 (violation rate) for each covariate
viol_rate <- data.frame(
  term = ph_pvals$term,
  prop_p_lt_0.05 = apply(ph_pvals[ , startsWith(names(ph_pvals), "imp2")], 1,
                         function(x) mean(x < 0.05, na.rm = TRUE))
)
viol_rate <- viol_rate[order(-viol_rate$prop_p_lt_0.05), ]

# Export p-values and summary
write.csv(ph_pvals,  "PH_pvalues_per_imputation_extended.csv", row.names = FALSE)
write.csv(viol_rate, "PH_violation_rate_extended.csv",         row.names = FALSE)

# Quick on-screen view (largest violation rates first)
viol_rate


# Fit 3-knot and 5-knot models (cause-specific)
fit_k3 <- with(imp2, {
  coxph(
    Surv(Survival_time_h, Status_cs1_h) ~
      ns(BMI, knots = K3, Boundary.knots = BND) +
      Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
      LVEF  + eGFR + Haemoglobin +
      ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
      Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
      + Stroke_TIA + ICD_status,
    x = TRUE, y = TRUE
  )
})

fit_k5 <- with(imp2, {
  coxph(
    Surv(Survival_time_h, Status_cs1_h) ~
      ns(BMI, knots = K5, Boundary.knots = BND) +
      Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
      LVEF  + eGFR + Haemoglobin +
      ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
      Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
      + Stroke_TIA + ICD_status,
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
  data.frame(imp2 = seq_along(aic_k3), AIC_K3 = aic_k3, AIC_K4 = aic_k4, AIC_K5 = aic_k5),
  "AIC_per_imputation_K3_K4_K5_extended.csv",
  row.names = FALSE
)
write.csv(aic_summary, "AIC_summary_K3_K4_K5_extended.csv", row.names = FALSE)

aic_summary

# compare base and extended model visually 
# Load curves
curve_base <- read.csv("cox_RCS_BMI_curve_cs1.csv")
curve_ext  <- read.csv("cox_RCS_BMI_curve_cs1_extended.csv")

##  Find minimum HR for each model 
i_base <- which.min(curve_base$HR)
opt_base <- with(curve_base, data.frame(
  BMI = BMI[i_base],
  HR  = HR[i_base],
  LCL = LCL[i_base],
  UCL = UCL[i_base]
))

i_ext <- which.min(curve_ext$HR)
opt_ext <- with(curve_ext, data.frame(
  BMI = BMI[i_ext],
  HR  = HR[i_ext],
  LCL = LCL[i_ext],
  UCL = UCL[i_ext]
))

opt_base
opt_ext   # print these to see the “optimal” BMI and HR + CI

# Find the minimum hazard ratio
HRmin_ext <- min(curve_ext$HR, na.rm = TRUE)

## Plot 
pdf("BMI_curve_comparison_base_vs_extended_shaded.pdf", width = 7, height = 5)

# Base plot setup (suppress default x-axis so we can make it granular)
plot(curve_base$BMI, curve_base$HR, type = "n",
     ylim = range(c(curve_base$LCL, curve_base$UCL,
                    curve_ext$LCL,  curve_ext$UCL)),
     xlab = "BMI (kg/m²)",
     ylab = "Hazard ratio (ref = 25)",
     main = "Restricted cubic spline: Base vs Extended models",
     xaxt = "n")

# Custom, more granular x-axis
axis(1, at = seq(18, 36, by = 1))           # ticks every 1 kg/m²
# Optional minor ticks:
# axis(1, at = seq(18, 36, by = 0.5), tcl = -0.2, labels = FALSE)

# Shaded 95% CI for base model (blue)
polygon(c(curve_base$BMI, rev(curve_base$BMI)),
        c(curve_base$LCL, rev(curve_base$UCL)),
        col = adjustcolor("darkblue", alpha.f = 0.15), border = NA)

# Shaded 95% CI for extended model (red)
polygon(c(curve_ext$BMI, rev(curve_ext$BMI)),
        c(curve_ext$LCL, rev(curve_ext$UCL)),
        col = adjustcolor("firebrick", alpha.f = 0.15), border = NA)

# Mean HR curves
lines(curve_base$BMI, curve_base$HR, col = "darkblue",  lwd = 2)
lines(curve_ext$BMI,  curve_ext$HR,  col = "firebrick", lwd = 2)

# Reference lines
abline(h = 1,  lty = 3, col = "gray50")
abline(v = 25, lty = 3, col = "gray80")

## Mark the minimum-risk BMI for each model

# 1) Points at the minimum
points(opt_base$BMI, opt_base$HR, pch = 19, col = "darkblue")
points(opt_ext$BMI,  opt_ext$HR,  pch = 19, col = "firebrick")

# 2) Vertical dashed lines at each minimum
abline(v = opt_base$BMI, lty = 3, col = adjustcolor("darkblue", alpha.f = 0.6))
abline(v = opt_ext$BMI,  lty = 3, col = adjustcolor("firebrick", alpha.f = 0.6))

# 3) Error bars = CI for HR at that minimum BMI
arrows(opt_base$BMI, opt_base$LCL,
       opt_base$BMI, opt_base$UCL,
       angle = 90, code = 3, length = 0.05, col = "darkblue")

arrows(opt_ext$BMI, opt_ext$LCL,
       opt_ext$BMI, opt_ext$UCL,
       angle = 90, code = 3, length = 0.05, col = "firebrick")

# Legend
legend("topleft",
       legend = c("Base model", "Extended model"),
       col    = c("darkblue", "firebrick"),
       lwd    = 2,
       pch    = 19,
       bty    = "n")

dev.off()


# comparing knot versions 
get_pooled_spline_curve <- function(fit_list, knots, BND,
                                    bmi_ref = 25,
                                    ngrid = 400,
                                    label) {
  
  mods <- fit_list$analyses
  m    <- length(mods)
  
  bmi_min  <- max(min(imp2$data$BMI, na.rm = TRUE), BND[1])
  bmi_max  <- min(max(imp2$data$BMI, na.rm = TRUE), BND[2])
  bmi_grid <- seq(bmi_min, bmi_max, length.out = ngrid)
  
  Xbmi <- ns(bmi_grid, knots = knots, Boundary.knots = BND)
  Xref <- ns(bmi_ref,  knots = knots, Boundary.knots = BND)
  
  Lmat <- Xbmi - matrix(rep(Xref, each = nrow(Xbmi)),
                        ncol = ncol(Xbmi))
  
  coef_names <- names(coef(mods[[1]]))
  sidx <- grep("^ns\\(BMI", coef_names)
  
  pool_contrast <- function(Lrow){
    qi <- vapply(mods, function(fm)
      as.numeric(Lrow %*% coef(fm)[sidx]), numeric(1))
    
    ui <- vapply(mods, function(fm)
      as.numeric(Lrow %*% vcov(fm)[sidx, sidx, drop=FALSE] %*% Lrow),
      numeric(1))
    
    qbar <- mean(qi)
    ubar <- mean(ui)
    bvar <- var(qi)
    Tvar <- ubar + (1 + 1/m) * bvar
    
    c(logHR = qbar, se = sqrt(Tvar))
  }
  
  res <- t(apply(Lmat, 1, pool_contrast))
  
  data.frame(
    BMI   = bmi_grid,
    HR    = exp(res[, "logHR"]),
    LCL   = exp(res[, "logHR"] - 1.96 * res[, "se"]),
    UCL   = exp(res[, "logHR"] + 1.96 * res[, "se"]),
    model = label
  )
}

curve_k3 <- get_pooled_spline_curve(
  fit_list = fit_k3,
  knots    = K3,
  BND      = BND,
  label    = "3 internal knots"
)

curve_k4 <- get_pooled_spline_curve(
  fit_list = fit_list_cs1,   # your main model
  knots    = K4,
  BND      = BND,
  label    = "4 internal knots"
)

curve_k5 <- get_pooled_spline_curve(
  fit_list = fit_k5,
  knots    = K5,
  BND      = BND,
  label    = "5 internal knots"
)

curves_all <- rbind(curve_k3, curve_k4, curve_k5)

pdf("BMI_RCS_knot_sensitivity_cs1.pdf", width = 7, height = 5)

ylim_all <- range(curves_all$LCL, curves_all$UCL)

plot(curve_k4$BMI, curve_k4$HR, type = "n",
     ylim = ylim_all,
     xlab = "BMI (kg/m²)",
     ylab = "Hazard ratio (ref = 25)",
     main = "")

# HR curves
# Okabe–Ito colour-blind safe palette
col_k3 <- "#0072B2"  # blue
col_k4 <- "black"
col_k5 <- "#D55E00"  # orange

lines(curve_k3$BMI, curve_k3$HR, col = col_k3, lwd = 2, lty = 2)
lines(curve_k4$BMI, curve_k4$HR, col = col_k4, lwd = 2, lty = 1)
lines(curve_k5$BMI, curve_k5$HR, col = col_k5, lwd = 2, lty = 2)


# Reference lines
abline(h = 1,  lty = 3, col = "grey60")
abline(v = 25, lty = 3, col = "grey80")

legend("topleft",
       legend = c("3 internal knots", "4 internal knots", "5 internal knots"),
       col    = c(col_k3, "black", col_k5),
       lty    = c(2, 1, 2),
       lwd    = 2,
       bty    = "n")

dev.off()
