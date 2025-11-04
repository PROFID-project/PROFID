install.packages("survival")
install.packages("survminer")
install.packages("dplyr")
install.packages("svglite")
install.packages("rms")

library(survival)
library(survminer)
library(dplyr)
library(rms)
library(mice)
library(splines)
library(dplyr)
library(broom)

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
mean(sapply(fit_list$analyses, AIC))


# Pooled AIC: average AIC across imputations
aic_vec <- sapply(fit_list$analyses, AIC)
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



# Proportional Hazards checks 
ph_tbl <- map(fit_list$analyses, ~cox.zph(.x)$table[, "p"])
# quick summary: proportion of imputations with p<0.05 per covariate
ph_flag <- sapply(names(ph_tbl[[1]]), function(v){
  mean(sapply(ph_tbl, function(pv) ifelse(is.na(pv[v]), NA, pv[v] < 0.05)), na.rm=TRUE)
})
sort(ph_flag, decreasing=TRUE)

# Sensitivity: 3 and 5 knots 
fit_k3 <- with(imp, coxph(
  Surv(Survival_time, Status) ~ ns(BMI, knots=K3, Boundary.knots=BND) + ... , x=TRUE))
fit_k5 <- with(imp, coxph(
  Surv(Survival_time, Status) ~ ns(BMI, knots=K5, Boundary.knots=BND) + ... , x=TRUE))

mean(sapply(fit_k3$analyses, AIC))
mean(sapply(fit_k5$analyses, AIC))



