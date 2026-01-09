#############################################################
# Seasonal Cox & Extended Models for SCD
# Author: Alexia Sampri
#
# Data assumptions (df_cleaned.rds):
# - Time_zero_Ym : Date of follow-up start (Date)
# - Survival_time: follow-up duration in MONTHS
# - Status       : 0 = censored, 1 = SCD (event), 2 = non-SCD death (competing)
# - season       : factor with "Winter" as reference
# - LVEF_std     : numeric LVEF (%) for stratified analyses & adjustment
# - Age          : baseline age (years)
# - Sex          : "Male"/"Female"
#############################################################

## =========================================================
## 0. Setup
## =========================================================

library(dplyr)
library(lubridate)
library(survival)
library(ggplot2)
library(purrr)
library(broom)
library(readr)
library(cmprsk)
library(splines)

# Load cleaned cohort

data_dir   <- "T:/PROFID/data/processed"
output_dir <- "T:/PROFID/output/survival_seasonal"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

setwd(data_dir)

# Load cleaned cohort
df <- readRDS("df_cleaned.rds")


## =========================================================
## 1. Core Data Preparation
## =========================================================

df <- df %>%
  mutate(
    # Time scale: MONTHS (no conversion to years)
    fu_mo = Survival_time,
    
    # Cause-specific event indicator for SCD only
    status_scd = ifelse(Status == 1, 1, 0),
    
    # Month of follow-up start (1–12)
    month_start = month(Time_zero_Ym),
    
    # Ensure season factor with Winter as reference
    season = factor(season,
                    levels = c("Winter", "Spring", "Summer", "Autumn")),
    
    # Continuous seasonality based on month (for main sin/cos models)
    sin_time = sin(2 * pi * month_start / 12),
    cos_time = cos(2 * pi * month_start / 12),
    
    # Day of year (for 365-day sinusoid later)
    doy = yday(Time_zero_Ym),
    
    # Sex factor
    Sex = factor(Sex)
  )

# Remove non-positive follow-up if any
df <- df %>% filter(fu_mo > 0)

# Create Age_cat consistently (even if it already exists, we overwrite)
df <- df %>%
  mutate(
    Age_cat = case_when(
      Age < 55 ~ "<55",
      Age >= 55 & Age < 70 ~ "55-69",
      Age >= 70 ~ ">=70"
    ),
    Age_cat = factor(Age_cat, levels = c("<55", "55-69", ">=70"))
  )

# Create LVEF_cat for stratified models (main 3 categories)
df <- df %>%
  mutate(
    LVEF_cat = case_when(
      LVEF_std <= 35 ~ "≤35%",
      LVEF_std > 35 & LVEF_std <= 50 ~ "36–50%",
      LVEF_std > 50 ~ ">50%"
    ),
    LVEF_cat = factor(LVEF_cat,
                      levels = c("≤35%", "36–50%", ">50%"))
  )

# Drop unused levels to avoid set_names / factor level issues
df <- df %>%
  mutate(
    season   = droplevels(season),
    Age_cat  = droplevels(Age_cat),
    LVEF_cat = droplevels(LVEF_cat)
  )

# Quick sanity checks
cat("Status vs status_scd:\n")
print(table(df$Status, df$status_scd))
cat("\nSeason distribution:\n")
print(table(df$season, useNA = "ifany"))
cat("\nAge_cat distribution:\n")
print(table(df$Age_cat, useNA = "ifany"))
cat("\nLVEF_cat distribution:\n")
print(table(df$LVEF_cat, useNA = "ifany"))


#############################################################
## 2.3 Survival Analysis with Seasonal Components
#############################################################


## =========================================================
## 2.3.1 Base Model:
##       Season as Categorical (Winter = Reference)
## =========================================================

cox_base <- coxph(
  Surv(fu_mo, status_scd) ~ season,
  data = df
)
summary(cox_base)

cox_base_hr <- tidy(cox_base, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::filter(grepl("^season", term))

write_csv(cox_base_hr,
          file.path(output_dir, "cox_base_season_hr.csv"))

p_base <- ggplot(cox_base_hr,
                 aes(x = reorder(term, estimate),
                     y = estimate,
                     ymin = conf.low,
                     ymax = conf.high)) +
  geom_pointrange() +
  geom_hline(yintercept = 1, linetype = "dashed") +
  coord_flip() +
  theme_minimal(base_size = 14) +
  labs(
    x = "Season (vs Winter)",
    y = "Hazard ratio (95% CI)",
    title = "Seasonal Hazard of SCD (Base Cox Model)"
  )

ggsave(file.path(output_dir, "cox_base_season_forest.png"),
       plot = p_base, width = 7, height = 5, dpi = 300)


## =========================================================
## 2.3.2 Continuous Time:
##       Sine & Cosine of Calendar Month
## =========================================================

cox_sincos <- coxph(
  Surv(fu_mo, status_scd) ~ sin_time + cos_time,
  data = df
)
summary(cox_sincos)

cox_sincos_hr <- tidy(cox_sincos, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::filter(term %in% c("sin_time", "cos_time"))

write_csv(cox_sincos_hr,
          file.path(output_dir, "cox_sincos_hr.csv"))

pred_months <- data.frame(
  month_start = 1:12
) %>%
  mutate(
    sin_time = sin(2 * pi * month_start / 12),
    cos_time = cos(2 * pi * month_start / 12)
  )

pred_months$lin_pred <- predict(cox_sincos, newdata = pred_months, type = "lp")
pred_months$rel_hazard <- exp(pred_months$lin_pred - mean(pred_months$lin_pred))

write_csv(pred_months,
          file.path(output_dir, "continuous_seasonal_predicted_hazard_by_month.csv"))

p_cont <- ggplot(pred_months,
                 aes(x = month_start, y = rel_hazard)) +
  geom_line(size = 1.2) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    title = "Predicted Seasonal Hazard of SCD (Continuous Month sin/cos)",
    x = "Month",
    y = "Relative Hazard"
  )

ggsave(file.path(output_dir, "seasonal_continuous_month_sincos.png"),
       plot = p_cont, width = 7, height = 5, dpi = 300)


## =========================================================
## 2.3.3 Time-Varying Coefficients:
##       Seasonal Effects Over Follow-up Duration
## =========================================================

# Define follow-up bands in MONTHS: 0–24, 24–60, >60
df <- df %>%
  mutate(
    fu_band = cut(
      fu_mo,
      breaks = c(0, 24, 60, Inf),
      labels = c("0–24 months", "24–60 months", ">60 months"),
      right = FALSE
    ),
    fu_band = droplevels(fu_band)
  )

table(df$fu_band, useNA = "ifany")

cox_tv <- coxph(
  Surv(fu_mo, status_scd) ~ (sin_time + cos_time) * fu_band,
  data = df
)
summary(cox_tv)

anova(cox_sincos, cox_tv, test = "LRT")

cox_tv_hr <- tidy(cox_tv, exponentiate = TRUE, conf.int = TRUE)
write_csv(cox_tv_hr,
          file.path(output_dir, "cox_timevarying_season_sincos_by_band.csv"))

pred_months_tv <- expand.grid(
  month_start = 1:12,
  fu_band = levels(df$fu_band)
) %>%
  mutate(
    sin_time = sin(2 * pi * month_start / 12),
    cos_time = cos(2 * pi * month_start / 12)
  )

pred_months_tv$lin_pred <- predict(cox_tv, newdata = pred_months_tv, type = "lp")

pred_months_tv <- pred_months_tv %>%
  group_by(fu_band) %>%
  mutate(rel_hazard = exp(lin_pred - mean(lin_pred))) %>%
  ungroup()

write_csv(pred_months_tv,
          file.path(output_dir, "timevarying_seasonal_predicted_hazard_by_band.csv"))

p_tv <- ggplot(pred_months_tv,
               aes(x = month_start, y = rel_hazard, colour = fu_band)) +
  geom_line(size = 1.1) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    title = "Time-Varying Seasonal Hazard of SCD by Follow-up Duration",
    x = "Month",
    y = "Relative Hazard",
    colour = "Follow-up Band"
  )

ggsave(file.path(output_dir, "seasonal_timevarying_by_band.png"),
       plot = p_tv, width = 7, height = 5, dpi = 300)


## =========================================================
## 2.3.4 Stratified Analysis:
##       LVEF Categories & Age Groups
## =========================================================

## ---------- 2.3.4a LVEF strata ----------

df_LVEF <- df %>%
  filter(!is.na(LVEF_cat)) %>%
  mutate(LVEF_cat = droplevels(LVEF_cat))

cox_by_LVEF <- df_LVEF %>%
  group_split(LVEF_cat) %>%
  set_names(levels(df_LVEF$LVEF_cat)) %>%
  map(~ coxph(Surv(fu_mo, status_scd) ~ sin_time + cos_time, data = .x))

cox_LVEF_results <- map_dfr(cox_by_LVEF, tidy, .id = "LVEF_cat")
write_csv(cox_LVEF_results,
          file.path(output_dir, "cox_continuous_season_by_LVEF.csv"))

pred_LVEF <- expand.grid(
  month_start = 1:12,
  LVEF_cat = levels(df_LVEF$LVEF_cat)
) %>%
  mutate(
    sin_time = sin(2 * pi * month_start / 12),
    cos_time = cos(2 * pi * month_start / 12)
  )

pred_LVEF <- pred_LVEF %>%
  group_by(LVEF_cat) %>%
  mutate(
    lin_pred = predict(
      cox_by_LVEF[[unique(LVEF_cat)]],
      newdata = cur_data(),
      type = "lp"
    ),
    rel_hazard = exp(lin_pred - mean(lin_pred))
  ) %>%
  ungroup()

write_csv(pred_LVEF,
          file.path(output_dir, "continuous_seasonal_predicted_by_LVEF.csv"))

p_LVEF <- ggplot(pred_LVEF,
                 aes(x = month_start, y = rel_hazard, colour = LVEF_cat)) +
  geom_line(size = 1.2) +
  scale_x_continuous(breaks = 1:12) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Seasonal Hazard of SCD by LVEF Category",
    x = "Month",
    y = "Relative Hazard",
    colour = "LVEF Category"
  )

ggsave(file.path(output_dir, "seasonal_continuous_by_LVEF.png"),
       plot = p_LVEF, width = 7, height = 5, dpi = 300)


## ---------- 2.3.4b Age strata ----------

df_age <- df %>%
  filter(!is.na(Age_cat)) %>%
  mutate(Age_cat = droplevels(Age_cat))

cox_by_age <- df_age %>%
  group_split(Age_cat) %>%
  set_names(levels(df_age$Age_cat)) %>%
  map(~ coxph(Surv(fu_mo, status_scd) ~ sin_time + cos_time, data = .x))

cox_age_results <- map_dfr(cox_by_age, tidy, .id = "Age_cat")
write_csv(cox_age_results,
          file.path(output_dir, "cox_continuous_season_by_Agecat.csv"))

pred_age <- expand.grid(
  month_start = 1:12,
  Age_cat = levels(df_age$Age_cat)
) %>%
  mutate(
    sin_time = sin(2 * pi * month_start / 12),
    cos_time = cos(2 * pi * month_start / 12)
  )

pred_age <- pred_age %>%
  group_by(Age_cat) %>%
  mutate(
    lin_pred = predict(
      cox_by_age[[unique(Age_cat)]],
      newdata = cur_data(),
      type = "lp"
    ),
    rel_hazard = exp(lin_pred - mean(lin_pred))
  ) %>%
  ungroup()

write_csv(pred_age,
          file.path(output_dir, "continuous_seasonal_predicted_by_Agecat.csv"))

p_age <- ggplot(pred_age,
                aes(x = month_start, y = rel_hazard, colour = Age_cat)) +
  geom_line(size = 1.2) +
  scale_x_continuous(breaks = 1:12) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Seasonal Hazard of SCD by Age Group",
    x = "Month",
    y = "Relative Hazard",
    colour = "Age group"
  )

ggsave(file.path(output_dir, "seasonal_continuous_by_Agecat.png"),
       plot = p_age, width = 7, height = 5, dpi = 300)


#############################################################
## 3. EXTRA MODELS / SENSITIVITY ANALYSES
#############################################################


## =========================================================
## 3.1 Fine–Gray Competing Risks (SCD vs non-SCD death)
## =========================================================

# Unadjusted Fine–Gray with season only TAKES TOO LONG TO RUN
X_fg_unadj <- model.matrix(~ season, data = df)[, -1]

fg_unadj <- crr(
  ftime   = df$fu_mo,
  fstatus = df$Status,  # 1 = SCD, 2 = competing, 0 = censored
  cov1    = X_fg_unadj
)
summary(fg_unadj)

fg_adj <- crr(
  ftime   = df$fu_mo,
  fstatus = df$Status,
  cov1    = cbind(
    season_Spring = ifelse(df$season == "Spring", 1, 0),
    season_Summer = ifelse(df$season == "Summer", 1, 0),
    season_Autumn = ifelse(df$season == "Autumn", 1, 0),
    Age      = df$Age,
    LVEF_std = df$LVEF_std,
    Sex      = as.numeric(df$Sex == "Male")
  )
)
summary(fg_adj)


# Adjusted Fine–Gray model: season + Age + Sex + LVEF_std
covariates <- c("season", "Age", "Sex", "LVEF_std")

X_fg_adj <- model.matrix(
  as.formula(paste("~", paste(covariates, collapse = " + "))),
  data = df
)[, -1]

fg_adj <- crr(
  ftime   = df$fu_mo,
  fstatus = df$Status,
  cov1    = X_fg_adj
)
summary(fg_adj)


## =========================================================
## 3.2 Spline Seasonality (Month as smooth term)
## =========================================================

# Unadjusted spline on month
cox_spline_unadj <- coxph(
  Surv(fu_mo, status_scd) ~ ns(month_start, df = 4),
  data = df
)
summary(cox_spline_unadj)

# Adjusted spline model
cox_spline_adj <- coxph(
  Surv(fu_mo, status_scd) ~ ns(month_start, df = 4) +
    Age + LVEF_std + Sex,
  data = df
)
summary(cox_spline_adj)

# Predict smooth adjusted seasonal curve (1–12 months)
pred_spline <- data.frame(
  month_start = seq(1, 12, by = 0.1),
  Age         = mean(df$Age, na.rm = TRUE),
  LVEF_std    = mean(df$LVEF_std, na.rm = TRUE),
  Sex         = levels(df$Sex)[1]
)

pred_spline$lin_pred <- predict(cox_spline_adj, newdata = pred_spline, type = "lp")
pred_spline$rel_hazard <- exp(pred_spline$lin_pred - mean(pred_spline$lin_pred))

p_spline <- ggplot(pred_spline, aes(x = month_start, y = rel_hazard)) +
  geom_line(size = 1.1) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    title = "Spline-based Seasonal Hazard of SCD (Adjusted)",
    x = "Month",
    y = "Relative hazard"
  )

ggsave(file.path(output_dir, "seasonal_spline_adjusted.png"),
       plot = p_spline, width = 7, height = 5, dpi = 300)


## =========================================================
## 3.3 365-day Sinusoid (Day-of-year resolution)
## =========================================================

df <- df %>%
  mutate(
    sin365 = sin(2 * pi * doy / 365.25),
    cos365 = cos(2 * pi * doy / 365.25)
  )

# Unadjusted 365-day sinusoid
cox_365_unadj <- coxph(
  Surv(fu_mo, status_scd) ~ sin365 + cos365,
  data = df
)
summary(cox_365_unadj)

# Adjusted 365-day sinusoid
cox_365_adj <- coxph(
  Surv(fu_mo, status_scd) ~ sin365 + cos365 + Age + LVEF_std + Sex,
  data = df
)
summary(cox_365_adj)

pred_doy <- data.frame(
  doy      = 1:365,
  Age      = mean(df$Age, na.rm = TRUE),
  LVEF_std = mean(df$LVEF_std, na.rm = TRUE),
  Sex      = levels(df$Sex)[1]
) %>%
  mutate(
    sin365 = sin(2 * pi * doy / 365.25),
    cos365 = cos(2 * pi * doy / 365.25)
  )

pred_doy$lin_pred <- predict(cox_365_adj, newdata = pred_doy, type = "lp")
pred_doy$rel_hazard <- exp(pred_doy$lin_pred - mean(pred_doy$lin_pred))

p_365 <- ggplot(pred_doy, aes(x = doy, y = rel_hazard)) +
  geom_line(size = 1.1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "365-day Sinusoidal Seasonal Hazard of SCD (Adjusted)",
    x = "Day of year",
    y = "Relative hazard"
  )

ggsave(file.path(output_dir, "seasonal_365day_sinusoid_adjusted.png"),
       plot = p_365, width = 7, height = 5, dpi = 300)


## =========================================================
## 3.4 Season × Age Interaction
## =========================================================

df <- df %>%
  mutate(
    Age_cat = droplevels(Age_cat),
    Age10   = Age / 10
  )

# Categorical Age_cat interaction
cox_season_agecat_adj <- coxph(
  Surv(fu_mo, status_scd) ~ season * Age_cat + LVEF_std + Sex,
  data = df
)
summary(cox_season_agecat_adj)

cox_season_agecat_no_int <- coxph(
  Surv(fu_mo, status_scd) ~ season + Age_cat + LVEF_std + Sex,
  data = df
)

anova(cox_season_agecat_no_int, cox_season_agecat_adj, test = "LRT")

# Continuous Age interaction (per 10 years)
cox_season_age_cont <- coxph(
  Surv(fu_mo, status_scd) ~ season * Age10 + LVEF_std + Sex,
  data = df
)
summary(cox_season_age_cont)


## =========================================================
## 3.5 Adjusted Versions of Main 2.3 Models
## =========================================================

# Adjusted base model (season categorical)
cox_base_adj <- coxph(
  Surv(fu_mo, status_scd) ~ season + Age + LVEF_std + Sex,
  data = df
)
summary(cox_base_adj)

# Adjusted month-based sin/cos model
cox_sincos_adj <- coxph(
  Surv(fu_mo, status_scd) ~ sin_time + cos_time + Age + LVEF_std + Sex,
  data = df
)
summary(cox_sincos_adj)

# Adjusted time-varying sin/cos × fu_band
cox_tv_adj <- coxph(
  Surv(fu_mo, status_scd) ~ (sin_time + cos_time) * fu_band +
    Age + LVEF_std + Sex,
  data = df
)
summary(cox_tv_adj)

#############################################################
# End of full seasonal + extended modelling script
#############################################################
