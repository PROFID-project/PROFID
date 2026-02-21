#############################################################
# 2.3 Survival Analysis with Seasonal Components – SCD
# Author: Alexia Sampri
#
# Data assumptions (df_cleaned.rds):
# - Time_zero_Ym : Date of follow-up start (Date)
# - Survival_time: follow-up duration in MONTHS
# - Status       : 0 = censored, 1 = SCD (event), 2 = non-SCD death (competing)
# - season       : factor with "Winter" as reference (already OK)
# - LVEF_std     : numeric LVEF (%) for stratified analyses
# - Age_cat      : age group factor (e.g. "<55", "55-69", ">=70")
#############################################################

library(dplyr)
library(lubridate)
library(survival)
library(ggplot2)
library(purrr)
library(broom)
library(readr)

# =========================================================
# 1. Data preparation
# =========================================================

data_dir   <- "T:/PROFID/data/processed"
output_dir <- "T:/PROFID/output/survival_seasonal"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

setwd(data_dir)

# Load cleaned cohort
df <- readRDS("df_cleaned.rds")

df <- df %>%
  mutate(
    # Time scale: MONTHS (no conversion to years)
    fu_mo = Survival_time,
    
    # Cause-specific event indicator for SCD only
    status_scd = ifelse(Status == 1, 1, 0),
    
    # Month of follow-up start (1–12) for continuous seasonal terms
    month_start = month(Time_zero_Ym),
    
    # If season not already a factor with Winter ref, enforce it
    season = factor(season,
                    levels = c("Winter", "Spring", "Summer", "Autumn")),
    
    # Continuous seasonality: sine/cosine of calendar month
    sin_time = sin(2 * pi * month_start / 12),
    cos_time = cos(2 * pi * month_start / 12)
  )

# Quick sanity checks
table(df$Status, df$status_scd)
table(df$season, useNA = "ifany")
summary(df$fu_mo)

# Remove any non-positive follow-up (if present)
df <- df %>% filter(fu_mo > 0)


# =========================================================
# 2.3.1 Base Model:
#      Season as Categorical (Winter = Reference)
# =========================================================

cox_base <- coxph(
  Surv(fu_mo, status_scd) ~ season,
  data = df
)
summary(cox_base)

# Tidy HRs for seasons
cox_base_hr <- tidy(cox_base, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::filter(grepl("^season", term))

write_csv(cox_base_hr,
          file.path(output_dir, "cox_base_season_hr.csv"))

# Forest plot for base model
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


# =========================================================
# 2.3.2 Continuous Time:
#      Sine & Cosine of Calendar Time
# =========================================================

cox_sincos <- coxph(
  Surv(fu_mo, status_scd) ~ sin_time + cos_time,
  data = df
)
summary(cox_sincos)

# HRs for sin/cos terms
cox_sincos_hr <- tidy(cox_sincos, exponentiate = TRUE, conf.int = TRUE) %>%
  dplyr::filter(term %in% c("sin_time", "cos_time"))

write_csv(cox_sincos_hr,
          file.path(output_dir, "cox_sincos_hr.csv"))

# Predicted relative hazard by month (continuous seasonal curve)
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
    title = "Predicted Seasonal Hazard of SCD (Continuous Time Model)",
    x = "Month",
    y = "Relative Hazard"
  )

ggsave(file.path(output_dir, "seasonal_continuous.png"),
       plot = p_cont, width = 7, height = 5, dpi = 300)


# =========================================================
# 2.3.3 Time-Varying Coefficients:
#      Seasonal Effects Vary Over Follow-up Duration
# =========================================================

# Define follow-up duration bands in MONTHS (e.g. 0–24, 24–60, >60)
df <- df %>%
  mutate(
    fu_band = cut(
      fu_mo,
      breaks = c(0, 24, 60, Inf),
      labels = c("0–24 months", "24–60 months", ">60 months"),
      right = FALSE
    )
  )

table(df$fu_band, useNA = "ifany")

# Model: interaction between continuous season and follow-up band
cox_tv <- coxph(
  Surv(fu_mo, status_scd) ~ (sin_time + cos_time) * fu_band,
  data = df
)
summary(cox_tv)

# Compare to simple sin/cos model (LRT)
anova(cox_sincos, cox_tv, test = "LRT")

cox_tv_hr <- tidy(cox_tv, exponentiate = TRUE, conf.int = TRUE)
write_csv(cox_tv_hr,
          file.path(output_dir, "cox_timevarying_season_sincos_by_band.csv"))

# Predicted seasonal hazard by follow-up band
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


# =========================================================
# 2.3.4 Stratified Analysis:
#      LVEF Categories & Age Groups
# =========================================================

# ---------- 2.3.4a LVEF categories ----------
df <- df %>%
  mutate(
    LVEF_cat = case_when(
      LVEF_std <= 35 ~ "≤35%",
      LVEF_std > 35 & LVEF_std <= 50 ~ "36–50%",
      LVEF_std > 50 ~ ">50%"
    ) %>%
      factor(levels = c("≤35%", "36–50%", ">50%"))
  ) %>%
  filter(!is.na(LVEF_cat))

table(df$LVEF_cat, useNA = "ifany")

cox_by_LVEF <- df %>%
  group_split(LVEF_cat) %>%
  set_names(levels(df$LVEF_cat)) %>%
  map(~ coxph(Surv(fu_mo, status_scd) ~ sin_time + cos_time, data = .x))

cox_LVEF_results <- map_dfr(cox_by_LVEF, tidy, .id = "LVEF_cat")
write_csv(cox_LVEF_results,
          file.path(output_dir, "cox_continuous_season_by_LVEF.csv"))

pred_LVEF <- expand.grid(
  month_start = 1:12,
  LVEF_cat = levels(df$LVEF_cat)
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


# ---------- 2.3.4b Age group strata (Age_cat already in data) ----------
df <- df %>%
  mutate(
    Age_cat = case_when(
      Age < 55 ~ "<55",
      Age >= 55 & Age < 70 ~ "55-69",
      Age >= 70 ~ ">=70"
    ),
    Age_cat = factor(Age_cat, levels = c("<55", "55-69", ">=70"))
  )

df_age <- df %>%
  filter(!is.na(Age_cat)) %>%
  mutate(Age_cat = droplevels(Age_cat))


table(df$Age_cat, useNA = "ifany")

cox_by_age <- df_age %>%
  group_split(Age_cat) %>%
  set_names(levels(df_age$Age_cat)) %>%
  map(~ coxph(Surv(fu_mo, status_scd) ~ sin_time + cos_time, data = .x))

cox_age_results <- map_dfr(cox_by_age, tidy, .id = "Age_cat")
write_csv(cox_age_results,
          file.path(output_dir, "cox_continuous_season_by_Agecat.csv"))

pred_age <- expand.grid(
  month_start = 1:12,
  Age_cat = levels(df$Age_cat)
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
# End of Section 2.3 – Seasonal Cox Models for SCD
#############################################################
