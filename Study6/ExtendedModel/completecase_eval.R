library(survival)
library(splines)
library(dplyr)

setwd("T:/Dokumente/PROFID/Study6")

# ---------------------------
# 1) Load MI object (for knots) and raw data (for complete-case model)
# ---------------------------

# Use the same object that the extended MI analysis used for knots
imp2 <- readRDS("mice_imputed_data_extended.RDS")

# Original (non-imputed) data for complete-case analysis
preimp <- read.csv("combined_BMI_outcomefiltered.csv")

# ---------------------------
# 2) Define knots EXACTLY as in extended_model (7).R
# ---------------------------

q   <- quantile(imp2$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])          # boundary knots: 5th & 95th percentiles
K4  <- as.numeric(q[c(2,3,4,5)])      # inner knots: 10,35,65,90

# ---------------------------
# 3) Define horizon + outcome (match extended MI analysis)
# ---------------------------

horizon <- 90

preimp <- preimp %>%
  mutate(
    Survival_time_h = pmin(Survival_time, horizon),
    Status_cs1_h    = ifelse(Survival_time <= horizon & Status == 1, 1L, 0L)
  )

# ---------------------------
# 4) Variables in the EXTENDED model (match vars_extended in MI script)
# ---------------------------

vars_extended <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute", "Cholesterol", "HDL", "LDL", "Triglycerides",
  "Stroke_TIA", "ICD_status"
)

# Variables needed for complete-case filtering and modeling
vars_needed <- c("BMI", "Survival_time_h", "Status_cs1_h", vars_extended)

# Safety check: do all required columns exist?
missing_cols <- setdiff(vars_needed, names(preimp))
if (length(missing_cols) > 0) {
  stop(paste0("Missing required columns in combined_BMI_outcomefiltered.csv: ",
              paste(missing_cols, collapse = ", ")))
}

# ---------------------------
# 5) Complete-case dataset
# ---------------------------

cca_data <- preimp %>%
  filter(complete.cases(across(all_of(vars_needed))))

cat("Rows in raw data:", nrow(preimp), "\n")
cat("Rows in complete-case data:", nrow(cca_data), "\n")

# ---------------------------
# 6) Fit complete-case extended Cox model (cause-specific, 90-month horizon)
# ---------------------------

cox_cca <- coxph(
  Surv(Survival_time_h, Status_cs1_h) ~
    ns(BMI, knots = K4, Boundary.knots = BND) +
    Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
    LVEF  + eGFR + Haemoglobin +
    ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
    Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
    Stroke_TIA + ICD_status,
  data = cca_data,
  x = TRUE, y = TRUE
)

saveRDS(cox_cca, "cox_RCS_cs1_extended_CCA.RDS")
print(summary(cox_cca))

# ---------------------------
# 7) BMI–HR curve vs BMI=25 with correct contrast-based SE/CI
#    (matches the MI curve logic; avoids se.fit error for contrasts)
# ---------------------------

# BMI grid (stays within the same boundary knots used in MI)
bmi_min  <- max(min(cca_data$BMI, na.rm = TRUE), BND[1])
bmi_max  <- min(max(cca_data$BMI, na.rm = TRUE), BND[2])
bmi_grid <- seq(bmi_min, bmi_max, length.out = 400)

# Spline basis for grid and reference BMI=25
Xbmi <- ns(bmi_grid, knots = K4, Boundary.knots = BND)
Xref <- ns(25,       knots = K4, Boundary.knots = BND)

# Contrast matrix: each row = (basis(BMI) - basis(25))
Lmat <- Xbmi - matrix(rep(Xref, each = nrow(Xbmi)), ncol = ncol(Xbmi))

# Extract spline coefficients + covariance
coef_names <- names(coef(cox_cca))
sidx <- grep("^ns\\(BMI", coef_names)

beta <- coef(cox_cca)[sidx]
V    <- vcov(cox_cca)[sidx, sidx, drop = FALSE]

# Log-HR and SE for contrast
logHR <- as.vector(Lmat %*% beta)
se    <- sqrt(diag(Lmat %*% V %*% t(Lmat)))

curve_cca <- data.frame(
  BMI = bmi_grid,
  HR  = exp(logHR),
  LCL = exp(logHR - 1.96 * se),
  UCL = exp(logHR + 1.96 * se)
)

write.csv(curve_cca, "cox_RCS_BMI_curve_cs1_extended_CCA.csv", row.names = FALSE)

# ---------------------------
# 8) Optional: C-index for the CC model (same horizon)
# ---------------------------

lp_cca <- predict(cox_cca, type = "lp")
ccon   <- survConcordance(Surv(Survival_time_h, Status_cs1_h) ~ lp_cca, data = cca_data)

cat("Complete-case C-index:", ccon$concordance, "\n")
saveRDS(ccon, "cindex_cs1_extended_CCA.RDS")

lp <- predict(cox_cca, type = "lp")

ccon <- survival::survConcordance(
  Surv(Survival_time_h, Status_cs1_h) ~ lp,
  data = cca_data
)

C  <- as.numeric(ccon$concordance)
se <- as.numeric(ccon$std.err)

ci_low  <- C - 1.96 * se
ci_high <- C + 1.96 * se

aic_cc <- AIC(cox_cca)

cat(sprintf("Complete-case C-index: %.3f (95%% CI %.3f–%.3f)\n", C, ci_low, ci_high))
cat(sprintf("Complete-case AIC: %.1f\n", aic_cc))

# --- assumes we already have cox_cca fitted on:
# Surv(Survival_time_h, Status_cs1_h) ~ ns(BMI, knots=K4, Boundary.knots=BND) + ...same covariates...

bmi_min  <- max(min(cca_data$BMI, na.rm = TRUE), BND[1])
bmi_max  <- min(max(cca_data$BMI, na.rm = TRUE), BND[2])
bmi_grid <- seq(bmi_min, bmi_max, length.out = 400)

Xbmi <- ns(bmi_grid, knots = K4, Boundary.knots = BND)
Xref <- ns(25,       knots = K4, Boundary.knots = BND)
Lmat <- Xbmi - matrix(rep(Xref, each = nrow(Xbmi)), ncol = ncol(Xbmi))

coef_names <- names(coef(cox_cca))
sidx <- grep("^ns\\(BMI", coef_names)

beta <- coef(cox_cca)[sidx]
V    <- vcov(cox_cca)[sidx, sidx, drop = FALSE]

logHR <- as.vector(Lmat %*% beta)
se    <- sqrt(diag(Lmat %*% V %*% t(Lmat)))

curve_cca <- data.frame(
  BMI = bmi_grid,
  HR  = exp(logHR),
  LCL = exp(logHR - 1.96 * se),
  UCL = exp(logHR + 1.96 * se)
)

write.csv(curve_cca, "cox_RCS_BMI_curve_cs1_extended_CCA.csv", row.names = FALSE)

# ---------------------------
# Plot MICE vs Complete-case curves and export PDF
# ---------------------------

# Read in the two curve files
curve_mice <- read.csv("cox_RCS_BMI_curve_cs1_extended.csv")
curve_cca  <- read.csv("cox_RCS_BMI_curve_cs1_extended_CCA.csv")

# Basic sanity checks
required_cols <- c("BMI", "HR", "LCL", "UCL")
stopifnot(all(required_cols %in% names(curve_mice)))
stopifnot(all(required_cols %in% names(curve_cca)))

# Set common x/y limits so both curves are visible
xlim <- range(c(curve_mice$BMI, curve_cca$BMI), na.rm = TRUE)
ylim <- range(c(curve_mice$LCL, curve_mice$UCL, curve_cca$LCL, curve_cca$UCL), na.rm = TRUE)

# Optional: cap the y-axis if you want a cleaner plot (edit/remove as needed)
# ylim[2] <- min(ylim[2], 2.5)

pdf("BMI_curve_MICE_vs_CCA_extended.pdf", width = 7.2, height = 5)

par(mar = c(4.6, 5.0, 1.4, 1.2))

# Empty plot canvas
plot(NA,
     xlim = xlim, ylim = ylim,
     xlab = expression(BMI~(kg/m^2)),
     ylab = "Hazard ratio (ref = 25)",
     axes = FALSE)

axis(1)
axis(2, las = 1)
box()

# Reference lines: HR=1 and BMI=25
abline(h = 1, lty = 3, col = "grey60")s
abline(v = 25, lty = 3, col = "grey60")

# --- MICE shaded CI and line (blue) ---
polygon(
  x = c(curve_mice$BMI, rev(curve_mice$BMI)),
  y = c(curve_mice$LCL, rev(curve_mice$UCL)),
  border = NA,
  col = adjustcolor("blue", alpha.f = 0.18)
)
lines(curve_mice$BMI, curve_mice$HR, col = "blue", lwd = 2)

# --- Complete-case shaded CI and line (red) ---
polygon(
  x = c(curve_cca$BMI, rev(curve_cca$BMI)),
  y = c(curve_cca$LCL, rev(curve_cca$UCL)),
  border = NA,
  col = adjustcolor("red", alpha.f = 0.12)
)
lines(curve_cca$BMI, curve_cca$HR, col = "red", lwd = 2)

# Legend
legend("topleft",
       legend = c("MICE (main analysis)", "Complete case"),
       col = c("blue", "red"),
       lwd = 2,
       bty = "n")

dev.off()

