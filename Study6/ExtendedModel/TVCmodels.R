library(survival)
library(mice)
library(dplyr)
library(splines)
library(broom)

setwd("T:/Dokumente/PROFID/Study6")

imp2 <- readRDS("mice_imputed_data_extended.RDS")

# ensure Status_cs1 exists
if (!"Status_cs1" %in% names(imp2$data)) {
  imp2$data$Status_cs1 <- ifelse(imp2$data$Status == 1, 1L, 0L)
}

dat1 <- complete(imp2, 1)

set.seed(123)
n_sub <- 20000   # if this still fails, reduce to e.g. 10000 or 5000
sub_ids <- sample(seq_len(nrow(dat1)), n_sub)
dat_sub <- dat1[sub_ids, ]

q   <- quantile(imp2$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])
K4  <- as.numeric(q[c(2,3,4,5)])

cox_std_sub <- coxph(
  Surv(Survival_time, Status_cs1) ~
    ns(BMI, knots = K4, Boundary.knots = BND) +
    Age + LVEF + eGFR +
    Sex + ICD_status + Diabetes + Hypertension + Smoking,
  data = dat_sub
)

summary(cox_std_sub)
saveRDS(cox_std_sub, "cox_standard_subsample_extended.RDS")

cox_tvc_sub <- coxph(
  Surv(Survival_time, Status_cs1) ~
    ns(BMI, knots = K4, Boundary.knots = BND) +
    Age + LVEF + eGFR +
    Sex + ICD_status + Diabetes + Hypertension + Smoking +
    tt(Age) + tt(LVEF) + tt(eGFR),
  data = dat_sub,
  tt = function(x, t, ...) x * log(t)
)

summary(cox_tvc_sub)
saveRDS(cox_tvc_sub, "cox_TVC_subsample_extended.RDS")

tidy_std <- tidy(cox_std_sub, conf.int = TRUE)
tidy_tvc <- tidy(cox_tvc_sub, conf.int = TRUE)

bmi_std <- subset(tidy_std, grepl("^ns\\(BMI", term))
bmi_tvc <- subset(tidy_tvc, grepl("^ns\\(BMI", term))

bmi_compare <- merge(
  bmi_std[, c("term","estimate","conf.low","conf.high")],
  bmi_tvc[, c("term","estimate","conf.low","conf.high")],
  by = "term",
  suffixes = c("_standard","_tvc")
)

bmi_compare
write.csv(bmi_compare,
          "BMI_spline_coefficients_standard_vs_TVC_subsample.csv",
          row.names = FALSE)


rm(list = ls())

library(survival)
library(splines)
library(ggplot2)
library(dplyr)

setwd("T:/Dokumente/PROFID/Study6")

# 1. Load subsample models (same 20k sample)
cox_std_sub <- readRDS("cox_standard_subsample_extended.RDS")
cox_tvc_sub <- readRDS("cox_TVC_subsample_extended.RDS")

# 2. Recreate BMI grid and knots (consistent with main analysis)
imp2 <- readRDS("mice_imputed_data_extended.RDS")
q   <- quantile(imp2$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])
K4  <- as.numeric(q[c(2,3,4,5)])

bmi_grid <- seq(BND[1], BND[2], length.out = 400)

Xbmi <- ns(bmi_grid, knots = K4, Boundary.knots = BND)
Xref <- ns(25,       knots = K4, Boundary.knots = BND)
Lmat <- Xbmi - matrix(rep(Xref, each = nrow(Xbmi)), ncol = ncol(Xbmi))

# Helper to get curve from a Cox model
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
    LCL   = exp(logHR - 1.96*se),
    UCL   = exp(logHR + 1.96*se),
    Model = label
  )
}

curve_std <- get_curve_single(cox_std_sub, "Standard Cox (20k subsample)")
curve_tvc <- get_curve_single(cox_tvc_sub, "TVC Cox (20k subsample)")

# For plotting: CI for standard only (keeps it clean)
df_lines <- bind_rows(
  curve_std %>% select(BMI, HR, Model),
  curve_tvc %>% select(BMI, HR, Model)
)

p <- ggplot() +
  geom_ribbon(
    data = curve_std,
    aes(x = BMI, ymin = LCL, ymax = UCL),
    fill = "steelblue", alpha = 0.15
  ) +
  geom_line(
    data = df_lines,
    aes(x = BMI, y = HR, colour = Model, linetype = Model),
    size = 1.1
  ) +
  geom_hline(yintercept = 1, linetype = "dotted", colour = "grey50") +
  geom_vline(xintercept = 25, linetype = "dotted", colour = "grey70") +
  scale_colour_manual(values = c(
    "Standard Cox (20k subsample)" = "steelblue4",
    "TVC Cox (20k subsample)"      = "firebrick3"
  )) +
  scale_linetype_manual(values = c("solid", "dashed")) +
  labs(
    x = "BMI (kg/m²)",
    y = "Hazard ratio for SCD (ref = 25 kg/m²)",
    title = "Restricted cubic spline for BMI\nStandard vs time-varying coefficients (same 20k subsample)"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "bottom",
    legend.title    = element_blank()
  )

p
ggsave("BMI_spline_standard_vs_TVC_subsample.pdf", p, width = 7, height = 5)
