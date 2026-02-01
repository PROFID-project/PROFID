library(survival)
library(mice)
library(dplyr)
library(splines)
library(broom)
library(ggplot2)

setwd("T:/Dokumente/PROFID/Study6")

# ----------------------------
# Load imputed extended mids
# ----------------------------
imp2 <- readRDS("mice_imputed_data_extended.RDS")

# Ensure cause-specific event indicator exists
if (!"Status_cs1" %in% names(imp2$data)) {
  imp2$data$Status_cs1 <- ifelse(imp2$data$Status == 1, 1L, 0L)
}

# ----------------------------
# Define horizon variables (match main analysis)
# ----------------------------
horizon <- 90
imp2$data <- imp2$data %>%
  mutate(
    Survival_time_h = pmin(Survival_time, horizon),
    Status_cs1_h    = ifelse(Survival_time <= horizon & Status == 1, 1L, 0L)
  )

# ----------------------------
# Define knots exactly as in extended model code
# ----------------------------
q   <- quantile(imp2$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])
K4  <- as.numeric(q[c(2,3,4,5)])

# ----------------------------
# Take the same 20k subsample (for memory reasons)
# IMPORTANT: keep the same individuals for std and TVC
# ----------------------------
set.seed(123)
dat1 <- complete(imp2, 1)

n_sub <- 20000
sub_ids <- sample(seq_len(nrow(dat1)), n_sub)
dat_sub <- dat1[sub_ids, ]

# ----------------------------
# Standard Cox on the subsample (MATCH extended model covariates)
# ----------------------------
cox_std_sub <- coxph(
  Surv(Survival_time_h, Status_cs1_h) ~
    ns(BMI, knots = K4, Boundary.knots = BND) +
    Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
    LVEF + eGFR + Haemoglobin +
    ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
    Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
    Stroke_TIA + ICD_status,
  data = dat_sub,
  x = TRUE, y = TRUE
)

saveRDS(cox_std_sub, "cox_standard_subsample_extended.RDS")

# ----------------------------
# TVC Cox on the same subsample
# Add log(time) interactions for variables flagged in PH checks
# (your paper text mentions Age, LVEF, eGFR) :contentReference[oaicite:1]{index=1}
# ----------------------------
cox_tvc_sub <- coxph(
  Surv(Survival_time_h, Status_cs1_h) ~
    ns(BMI, knots = K4, Boundary.knots = BND) +
    Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
    LVEF + eGFR + Haemoglobin +
    ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
    Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
    Stroke_TIA + ICD_status +
    tt(Age) + tt(eGFR) + tt(Haemoglobin),
  data = dat_sub,
  tt = function(x, t, ...) x * log(pmax(t, 1e-3)),  # avoid log(0)
  x = TRUE, y = TRUE
)

saveRDS(cox_tvc_sub, "cox_TVC_subsample_extended.RDS")

# ----------------------------
# Table S9-style: compare BMI spline coefficients
# ----------------------------
tidy_std <- tidy(cox_std_sub, conf.int = TRUE)
tidy_tvc <- tidy(cox_tvc_sub, conf.int = TRUE)

bmi_std <- subset(tidy_std, grepl("^ns\\(BMI", term))
bmi_tvc <- subset(tidy_tvc, grepl("^ns\\(BMI", term))

bmi_compare <- merge(
  bmi_std[, c("term","estimate","conf.low","conf.high")],
  bmi_tvc[, c("term","estimate","conf.low","conf.high")],
  by = "term",
  suffixes = c("_standard","_tvc")
) %>%
  mutate(
    diff = abs(estimate_tvc - estimate_standard)
  )

write.csv(
  bmi_compare,
  "BMI_spline_coefficients_standard_vs_TVC_subsample.csv",
  row.names = FALSE
)

# ----------------------------
# Figure S6-style: BMI curve (standard vs TVC) on same 20k subsample
# ----------------------------
bmi_grid <- seq(BND[1], BND[2], length.out = 400)

Xbmi <- ns(bmi_grid, knots = K4, Boundary.knots = BND)
Xref <- ns(25,       knots = K4, Boundary.knots = BND)
Lmat <- Xbmi - matrix(rep(Xref, each = nrow(Xbmi)), ncol = ncol(Xbmi))

get_curve_single <- function(fit, label) {
  cf   <- coef(fit)
  V    <- vcov(fit)
  sidx <- grep("^ns\\(BMI", names(cf))
  beta <- cf[sidx]
  Vb   <- V[sidx, sidx, drop = FALSE]
  
  logHR <- as.numeric(Lmat %*% beta)
  se    <- sqrt(rowSums((Lmat %*% Vb) * Lmat))
  
  data.frame(
    BMI   = bmi_grid,
    HR    = exp(logHR),
    LCL   = exp(logHR - 1.96 * se),
    UCL   = exp(logHR + 1.96 * se),
    Model = label
  )
}

curve_std <- get_curve_single(cox_std_sub, "Standard Cox (20k subsample)")
curve_tvc <- get_curve_single(cox_tvc_sub, "TVC Cox (20k subsample)")

df_all <- bind_rows(curve_std, curve_tvc)

library(ggplot2)

# Choose which model should be dashed (the one drawn "on top")
top_model <- "TVC Cox (20k subsample)"   # change if you want the other dashed

p <- ggplot(df_all, aes(x = BMI, y = HR, colour = Model, fill = Model)) +
  geom_ribbon(aes(ymin = LCL, ymax = UCL), alpha = 0.15, colour = NA) +
  # Draw solid line first (underneath) by subsetting
  geom_line(
    data = subset(df_all, Model != top_model),
    linewidth = 1.1,
    linetype = "solid"
  ) +
  # Draw dashed line second (on top)
  geom_line(
    data = subset(df_all, Model == top_model),
    linewidth = 1.1,
    linetype = "dashed"
  ) +
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_vline(xintercept = 25, linetype = "dotted") +
  labs(
    x = "BMI (kg/m²)",
    y = "Hazard ratio for SCD (ref = 25 kg/m²)"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave("BMI_spline_standard_vs_TVC_subsample_dashed.pdf", p, width = 8, height = 5)
