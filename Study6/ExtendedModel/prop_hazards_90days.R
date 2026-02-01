
library(mice)
library(survival)
library(splines)
library(dplyr)
library(ggplot2)
library(survminer)
library(patchwork)

t0 <- 90

make_90day <- function(d,
                       time_var   = "Survival_time",
                       status_var = "Status_cs1") {
  
  tt <- d[[time_var]]
  ss <- d[[status_var]]s
  s
  # Ensure ss is 0/1 integer (Status_cs1 should already be, but be safe)
  ss <- suppressWarnings(as.integer(as.character(ss)))
  
  d$Survival_time_90 <- pmin(tt, t0)
  d$Status_cs1_90    <- ifelse(!is.na(ss) & ss == 1L & tt <= t0, 1L, 0L)
  
  d
}


# imp2 already loaded in your script
# imp2 <- readRDS("mice_imputed_data_extended.RDS")

form_ext_90 <- Surv(Survival_time_90, Status_cs1_90) ~
  ns(BMI, knots = K4, Boundary.knots = BND) +
  Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
  LVEF + eGFR + Haemoglobin +
  ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
  Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
  COPD_cat + Cancer + Stroke_TIA + ICD_status

fit_ext_90 <- with(imp2, {
  d <- data
  d <- make_90day(data)
  coxph(form_ext_90, data = d, ties = "efron", x = TRUE)
})

# Use the first imputation for diagnostic plots (standard approach)
fm1_90 <- fit_ext_90$analyses[[1]]


czph_90 <- cox.zph(fm1_90, transform = "km")

pdf("PH_Schoenfeld_panels_extended_90day.pdf", width = 8, height = 10)

par(mfrow = c(2,2),
    mar = c(4.5, 4.8, 3, 1),
    oma = c(0, 0, 2, 0))

yl <- "Scaled Schoenfeld residuals"

# BMI spline term name in cox.zph will match how it appears in the model
# You can confirm with: rownames(czph_90$table)
plot(czph_90, var = "ns(BMI, knots = K4, Boundary.knots = BND)",
     resid = TRUE, se = TRUE,
     main = "a) BMI spline", ylab = yl, xlim = c(0, 90))

plot(czph_90, var = "Age",
     resid = TRUE, se = TRUE,
     main = "b) Age", ylab = yl, xlim = c(0, 90))

plot(czph_90, var = "LVEF",
     resid = TRUE, se = TRUE,
     main = "c) LVEF", ylab = yl, xlim = c(0, 90))

plot(czph_90, var = "eGFR",
     resid = TRUE, se = TRUE,
     main = "d) eGFR", ylab = yl, xlim = c(0, 90))

mtext(
  paste0("Global proportional hazards test (90-day censored): p = ",
         format.pval(czph_90$table["GLOBAL","p"], digits = 3, eps = 1e-3)),
  outer = TRUE, cex = 0.95
)

dev.off()
