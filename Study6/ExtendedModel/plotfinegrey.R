## plot_FineGray_BMI_SCD.R
## Uses the already-fitted Fine–Gray model to plot BMI spline curve
## WITHOUT refitting the expensive model.

library(mice)
library(splines)

setwd("T:/Dokumente/PROFID/Study6")

# ---- Load data + model ------------------------------------------------------

# Extended imputed object (used only to get BMI distribution + Sex levels)
imp_ext  <- readRDS("mice_imputed_data_extended.RDS")

# Fine–Gray model fitted previously: BMI spline + Age + Sex
fg_small <- readRDS("FineGray_fg_small_BMI_AgeSex.RDS")

# Take first completed dataset to get Age/Sex structure
dat1 <- complete(imp_ext, 1)

# ---- Recreate knots exactly as in original analysis -------------------------

q <- quantile(
  dat1$BMI,
  probs = c(.05, .10, .35, .65, .90, .95),
  na.rm = TRUE
)

BND <- as.numeric(q[c(1, 6)])      # 5th & 95th
K4  <- as.numeric(q[c(2, 3, 4, 5)])# 10th, 35th, 65th, 90th

# ---- Build prediction grid --------------------------------------------------

bmi_grid <- seq(18, 40, length.out = 400)

pred_df <- data.frame(
  BMI = bmi_grid,
  Age = median(dat1$Age, na.rm = TRUE),
  Sex = levels(dat1$Sex)[1]   # reference sex (e.g. "Female")
)

# Make sure Sex has same levels as in the model
if (is.factor(dat1$Sex)) {
  pred_df$Sex <- factor(pred_df$Sex, levels = levels(dat1$Sex))
}

# Design matrix matching the Fine–Gray model (no intercept)
X <- model.matrix(
  ~ ns(BMI, knots = K4, Boundary.knots = BND) + Age + Sex,
  data = pred_df
)
X <- X[, -1, drop = FALSE]

# ---- Predict subdistribution hazard ratios vs BMI ---------------------------

lp <- as.vector(X %*% fg_small$coef)

# Use BMI 25 as reference
ref_idx <- which.min(abs(bmi_grid - 25))
HR <- exp(lp - lp[ref_idx])

# Optional: save curve data for reproducibility / plotting elsewhere
curve_fg <- data.frame(BMI = bmi_grid, HR = HR)
write.csv(curve_fg, "FineGray_BMI_SCD_curve_data.csv", row.names = FALSE)

# ---- Make and save the plot -------------------------------------------------

pdf("FineGray_BMI_SCD_curve.pdf", width = 7, height = 5)

plot(
  bmi_grid, HR,
  type = "l", lwd = 3,
  xlab = "BMI (kg/m²)",
  ylab = "Subdistribution hazard ratio (ref = BMI 25)",
  main = "Fine–Gray: BMI and SCD risk",
  ylim = range(HR) * c(0.95, 1.05)
)
abline(h = 1,  lty = 3, col = "gray40")
abline(v = 25, lty = 3, col = "gray40")

dev.off()

## compare_Cox_vs_FG_BMI.R
## Overlay cause-specific Cox RCS curve and Fine–Gray subdistribution curve

setwd("T:/Dokumente/PROFID/Study6")

# ---- Load previously saved curves ------------------------------------------

# Cause-specific Cox RCS curve (extended model)
cox_curve <- read.csv("cox_RCS_BMI_curve_cs1_extended.csv")

# Fine–Gray curve (BMI spline + Age + Sex)
fg_curve  <- read.csv("FineGray_BMI_SCD_curve_data.csv")

# Optional: restrict to common BMI range
bmi_min <- max(min(cox_curve$BMI), min(fg_curve$BMI))
bmi_max <- min(max(cox_curve$BMI), max(fg_curve$BMI))

cox_sub <- subset(cox_curve, BMI >= bmi_min & BMI <= bmi_max)
fg_sub  <- subset(fg_curve,  BMI >= bmi_min & BMI <= bmi_max)

# ---- Overlay plot ----------------------------------------------------------

pdf("Cox_vs_FineGray_BMI_overlay.pdf", width = 7, height = 5)

plot(
  cox_sub$BMI, cox_sub$HR,
  type = "l", lwd = 3,
  xlab = "BMI (kg/m²)",
  ylab = "Relative hazard (ref = BMI 25)",
  main = "Cause-specific Cox vs Fine–Gray: BMI and SCD risk",
  ylim = range(c(cox_sub$HR, fg_sub$HR), na.rm = TRUE) * c(0.95, 1.05)
)

# Add Cox 95% CI band (optional, shaded)
polygon(
  c(cox_sub$BMI, rev(cox_sub$BMI)),
  c(cox_sub$LCL, rev(cox_sub$UCL)),
  border = NA,
  col = adjustcolor("gray80", alpha.f = 0.6)
)

# Redraw Cox central curve on top of band
lines(cox_sub$BMI, cox_sub$HR, lwd = 3)

# Add Fine–Gray curve (dashed)
lines(fg_sub$BMI, fg_sub$HR, lwd = 2, lty = 2)

# Reference lines
abline(h = 1,  lty = 3, col = "gray40")
abline(v = 25, lty = 3, col = "gray40")

legend(
  "topleft",
  legend = c("Cause-specific Cox (RCS)", "Fine–Gray (BMI spline, age+sex)"),
  lwd    = c(3, 2),
  lty    = c(1, 2),
  bty    = "n"
)

dev.off()

