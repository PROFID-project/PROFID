## plotfinegrey_90mo.R
## Plot Fine–Gray BMI spline curve and optionally overlay with cause-specific Cox curve (90-month horizon)

library(mice)
library(splines)
library(cmprsk)   # for crr() if refitting small model
library(grDevices)

setwd("T:/Dokumente/PROFID/Study6")

# -----------------------
# Settings
# -----------------------
horizon <- 90

# Input/output files
imp_file   <- "mice_imputed_data_extended.RDS"
fg_model_file <- "FineGray_fg_small_BMI_AgeSex_90mo.RDS"

fg_curve_csv <- "FineGray_BMI_SCD_curve_data_90mo.csv"
fg_curve_pdf <- "FineGray_BMI_SCD_curve_90mo.pdf"

cox_curve_csv <- "cox_RCS_BMI_curve_cs1_extended_90mo.csv"
overlay_pdf   <- "Cox_vs_FineGray_BMI_overlay_90mo.pdf"

# If the small Fine–Gray model file doesn't exist, refit it on imputation 1 and save
refit_if_missing <- TRUE

# -----------------------
# Load imputed data
# -----------------------
imp_ext <- readRDS(imp_file)
dat1 <- complete(imp_ext, 1)

# Assumed competing risks coding:
#   0 = censored
#   1 = SCD (event of interest)
#   2 = non-cardiac death (competing event)
#
# Apply 90-month horizon:
dat1$Status_fg <- dat1$Status
dat1$ftime_h   <- pmin(dat1$Survival_time, horizon)
dat1$fstatus_h <- ifelse(dat1$Survival_time > horizon, 0L, dat1$Status_fg)
dat1$fstatus_h <- as.integer(dat1$fstatus_h)

# -----------------------
# Recreate knots as in main analysis
# -----------------------
q <- quantile(
  dat1$BMI,
  probs = c(.05, .10, .35, .65, .90, .95),
  na.rm = TRUE
)
BND <- as.numeric(q[c(1, 6)])       # 5th & 95th
K4  <- as.numeric(q[c(2, 3, 4, 5)]) # 10th, 35th, 65th, 90th

# -----------------------
# Fit or load the Fine–Gray small model (BMI spline + Age + Sex)
# -----------------------
if (!file.exists(fg_model_file)) {
  if (!refit_if_missing) stop("Fine–Gray model file not found: ", fg_model_file)
  
  message("Fine–Gray model not found. Refitting small model and saving to: ", fg_model_file)
  
  mm <- model.matrix(
    ~ ns(BMI, knots = K4, Boundary.knots = BND) + Age + Sex,
    data = dat1
  )
  mm <- mm[, -1, drop = FALSE]
  
  fg_small <- crr(
    ftime   = dat1$ftime_h,
    fstatus = dat1$fstatus_h,
    cov1    = mm
  )
  
  saveRDS(fg_small, fg_model_file)
} else {
  fg_small <- readRDS(fg_model_file)
}

# -----------------------
# Build prediction grid
# -----------------------
bmi_grid <- seq(18, 40, length.out = 400)

pred_df <- data.frame(
  BMI = bmi_grid,
  Age = median(dat1$Age, na.rm = TRUE),
  Sex = if (is.factor(dat1$Sex)) levels(dat1$Sex)[1] else dat1$Sex[which(!is.na(dat1$Sex))[1]]
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

# Predict subdistribution hazard ratios vs BMI (relative to BMI 25)
lp <- as.vector(X %*% fg_small$coef)

ref_idx <- which.min(abs(bmi_grid - 25))
HR <- exp(lp - lp[ref_idx])

curve_fg <- data.frame(BMI = bmi_grid, HR = HR)
write.csv(curve_fg, fg_curve_csv, row.names = FALSE)

# Plot Fine–Gray curve
pdf(fg_curve_pdf, width = 7, height = 5)
plot(
  bmi_grid, HR,
  type = "l", lwd = 3,
  xlab = "BMI (kg/m²)",
  ylab = "Subdistribution hazard ratio (ref = BMI 25)",
  main = paste0("Fine–Gray: BMI and SCD risk (", horizon, "-month horizon)"),
  ylim = range(HR, na.rm = TRUE) * c(0.95, 1.05)
)
abline(h = 1,  lty = 3, col = "gray40")
abline(v = 25, lty = 3, col = "gray40")
dev.off()

# -----------------------
# Overlay cause-specific Cox curve (if available)
# -----------------------
if (file.exists(cox_curve_csv)) {
  
  cox_curve <- read.csv(cox_curve_csv)
  fg_curve  <- read.csv(fg_curve_csv)
  
  # Restrict to common BMI range
  bmi_min <- max(min(cox_curve$BMI), min(fg_curve$BMI))
  bmi_max <- min(max(cox_curve$BMI), max(fg_curve$BMI))
  
  cox_sub <- subset(cox_curve, BMI >= bmi_min & BMI <= bmi_max)
  fg_sub  <- subset(fg_curve,  BMI >= bmi_min & BMI <= bmi_max)
  
  pdf(overlay_pdf, width = 7, height = 5)
  
  plot(
    cox_sub$BMI, cox_sub$HR,
    type = "l", lwd = 3,
    xlab = "BMI (kg/m²)",
    ylab = "Relative hazard (ref = BMI 25)",
    main = paste0("Cause-specific Cox vs Fine–Gray (", horizon, "-month horizon)"),
    ylim = range(c(cox_sub$HR, fg_sub$HR, cox_sub$LCL, cox_sub$UCL), na.rm = TRUE) * c(0.95, 1.05)
  )
  
  # Cox 95% CI band (if present)
  if (all(c("LCL","UCL") %in% names(cox_sub))) {
    polygon(
      c(cox_sub$BMI, rev(cox_sub$BMI)),
      c(cox_sub$LCL, rev(cox_sub$UCL)),
      border = NA,
      col = adjustcolor("gray80", alpha.f = 0.6)
    )
    lines(cox_sub$BMI, cox_sub$HR, lwd = 3)
  }
  
  # Fine–Gray curve (dashed)
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
} else {
  message("Cox curve CSV not found (skipping overlay): ", cox_curve_csv)
}
