
library(survival)
library(dplyr)
library(ggplot2)
library(survminer)
library(patchwork)

# Set working directory
setwd("T:/Dokumente/PROFID/Study6")

# Load objects ----
# Imputed mids object with extended variables
imp2         <- readRDS("mice_imputed_data_extended.RDS")

# List of Cox models across imputations (extended BMI spline model)
fit_list_cs1 <- readRDS("fit_list_cs1_extended.RDS")


# Recompute knots
q   <- quantile(imp2$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])              # boundary knots 5th & 95th
K4  <- as.numeric(q[c(2,3,4,5)])          # inner knots 10,35,65,90


# Quick sanity check
length(fit_list_cs1$analyses)
fit_list_cs1$analyses[[1]]
names(imp2$data)[1:20]


# Schoenfeld residual tests per imputation ----

# Run cox.zph on each imputed model
ph_tbl_list <- lapply(fit_list_cs1$analyses, function(fm) {
  as.data.frame(cox.zph(fm)$table)
})

# Combine p-values across imputations
# Each table has rows = terms (including "GLOBAL"), columns = rho, chisq, p
ph_pvals <- Reduce(
  function(a, b) cbind(a["p"], b["p"]),
  ph_tbl_list
)

colnames(ph_pvals) <- paste0("imp", seq_len(ncol(ph_pvals)))
ph_pvals$term <- rownames(ph_tbl_list[[1]])
ph_pvals <- ph_pvals[, c("term", grep("^imp", names(ph_pvals), value = TRUE))]

# Proportion of imputations with p < 0.05 (evidence of PH violation)
viol_rate <- ph_pvals %>%
  filter(term != "GLOBAL") %>%  # drop global if you want per-covariate view
  rowwise() %>%
  mutate(
    prop_p_lt_0.05 = mean(c_across(starts_with("imp")) < 0.05, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(desc(prop_p_lt_0.05))

viol_rate

# Save PH test outputs
write.csv(ph_pvals,
          "PH_pvalues_per_imputation_extended_new.csv",
          row.names = FALSE)

write.csv(viol_rate,
          "PH_violation_rate_extended_new.csv",
          row.names = FALSE)


fm1  <- fit_list_cs1$analyses[[1]]
czph <- cox.zph(fm1)

# Helper to plot all rows that match a variable name
plot_terms <- function(zph_obj, pattern, title = NULL) {
  idx <- grep(pattern, rownames(zph_obj$table))
  if (length(idx) == 0) {
    message("No coefficient rows match: ", pattern)
  } else {
    for (i in idx) {
      plot(zph_obj[i], main = paste0(pattern, " — ", rownames(zph_obj$table)[i]))
    }
  }
}

pdf("PH_selected_covariates_imp1_extended.pdf", width = 8, height = 10)

# BMI (all spline basis terms)
plot_terms(czph, "^ns\\(BMI")

# Continuous variables
plot_terms(czph, "^Age")
plot_terms(czph, "^LVEF")
plot_terms(czph, "^eGFR")

# Binary/factors
plot_terms(czph, "^Sex")
plot_terms(czph, "^ICD_status")

dev.off()


# Log(-log) survival plots for key groups 


dat_all <- imp2$data

# Sex
sf_sex <- survfit(Surv(Survival_time, Status_cs1) ~ Sex, data = dat_all)

p_sex <- ggsurvplot(
  sf_sex,
  fun = "cloglog",
  xlab = "Time",
  ylab = "log(-log(Survival))",
  legend.title = "Sex",
  legend.labs = levels(dat_all$Sex),
  palette = c("#1F77B4", "#FF7F0E")
)$plot +
  ggtitle("log(-log(Survival)) by Sex")

# ICD status
sf_icd <- survfit(Surv(Survival_time, Status_cs1) ~ ICD_status, data = dat_all)

p_icd <- ggsurvplot(
  sf_icd,
  fun = "cloglog",
  xlab = "Time",
  ylab = "log(-log(Survival))",
  legend.title = "ICD status",
  legend.labs = levels(dat_all$ICD_status),
  palette = c("#2CA02C", "#D62728")
)$plot +
  ggtitle("log(-log(Survival)) by ICD status")

#  LVEF grouped (e.g. <35, 35–50, >50)
dat_all <- dat_all %>%
  mutate(
    LVEF_group = cut(
      LVEF,
      breaks = c(-Inf, 35, 50, Inf),
      labels = c("<35%", "35–50%", ">50%")
    )
  )

sf_lvef <- survfit(Surv(Survival_time, Status_cs1) ~ LVEF_group, data = dat_all)

p_lvef <- ggsurvplot(
  sf_lvef,
  fun = "cloglog",
  xlab = "Time",
  ylab = "log(-log(Survival))",
  legend.title = "LVEF group",
  legend.labs = levels(dat_all$LVEF_group),
  palette = c("#9467BD", "#8C564B", "#17BECF")
)$plot +
  ggtitle("log(-log(Survival)) by LVEF group")

# Combine plots into one panel and save
panel_loglog <- (p_sex | p_icd) / p_lvef

ggsave(
  "loglog_PH_checks_extended_imputed.pdf",
  panel_loglog,
  width = 10, height = 8
)


czph <- cox.zph(fm1)

# plot all rows that match a variable name
plot_terms <- function(zph_obj, pattern, title = NULL) {
  idx <- grep(pattern, rownames(zph_obj$table))
  if (length(idx) == 0) {
    message("No coefficient rows match: ", pattern)
  } else {
    for (i in idx) {
      plot(zph_obj[i], main = paste0(pattern, " — ", rownames(zph_obj$table)[i]))
    }
  }
}


terms_all <- rownames(czph$table)
terms_all <- terms_all[terms_all != "GLOBAL"]  # drop global row

fm1  <- fit_list_cs1$analyses[[1]]
czph <- cox.zph(fm1)

# all terms except GLOBAL
idx <- which(rownames(czph$table) != "GLOBAL")

pdf("PH_Schoenfeld_all_covariates_panels_imp1_extended.pdf",
    width = 8, height = 10)

n_rows <- 3
n_cols <- 2
par(mfrow = c(n_rows, n_cols))

for (i in idx) {
  term_name <- rownames(czph$table)[i]
  plot(czph[i], main = term_name)
 
}

dev.off()
par(mfrow = c(1, 1))

# BMI (all spline basis terms)
plot_terms(czph, "^ns\\(BMI")

# Continuous variables
plot_terms(czph, "^Age")
plot_terms(czph, "^LVEF")
plot_terms(czph, "^eGFR")

# Binary/factors
plot_terms(czph, "^Sex")
plot_terms(czph, "^ICD_status")

dev.off()

# Variables with high violation rate
head(viol_rate, 10)


#below is code for the 90 day cutoff 

library(mice)
library(survival)
library(splines)

# Assumes these exist in environment from main extended modelling code:
# K4 (inner knots), BND (boundary knots)


form_ext_90 <- as.formula(
  paste0(
    "Surv(time90, status90) ~ ",
    "ns(BMI, knots = K4, Boundary.knots = BND) + ",
    "Age + Sex + Diabetes + Hypertension + Smoking + MI_history + ",
    "LVEF + eGFR + Haemoglobin + ",
    "ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering + ",
    "Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides + ",
    "COPD_cat + Cancer + Stroke_TIA + ICD_status"
  )
)

fit_ext_90 <- with(imp2, {
  coxph(
    Surv(pmin(Survival_time, 90),
         as.integer(Status_cs1 == 1 & Survival_time <= 90)) ~
      ns(BMI, knots = K4, Boundary.knots = BND) +
      Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
      LVEF + eGFR + Haemoglobin +
      ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
      Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
      COPD_cat + Cancer + Stroke_TIA + ICD_status,
    ties = "efron",
    x = TRUE
  )
})

pool_fit_ext_90 <- pool(fit_ext_90)
summary(pool_fit_ext_90)

# pick one fitted model (e.g., imputation 1)
fit1 <- fit_ext_90$analyses[[1]]

cz <- cox.zph(fit1, transform = "km")
print(cz)


pdf("schoenfeld_panels_BMI_Age_LVEF_eGFR.pdf", width = 8.5, height = 9)

op <- par(mfrow = c(2,2), mar = c(4,4,3,1), oma = c(0,0,3,0))



plot(cz, var = "ns(BMI, knots = K4, Boundary.knots = BND)",
     main = "a) BMI spline", xlab = "Time", ylab = "Scaled Schoenfeld residuals")
plot(cz, var = "Age",
     main = "b) Age", xlab = "Time", ylab = "Scaled Schoenfeld residuals")
plot(cz, var = "LVEF",
     main = "c) LVEF", xlab = "Time", ylab = "Scaled Schoenfeld residuals")
plot(cz, var = "eGFR",
     main = "d) eGFR", xlab = "Time", ylab = "Scaled Schoenfeld residuals")

mtext(sprintf("Global proportional hazards test: p = %.3f",
              cz$table["GLOBAL","p"]),
      outer = TRUE, cex = 1.0)


par(op)
dev.off()

