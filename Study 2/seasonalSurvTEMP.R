#############################################################
# Seasonal Cox Proportional Hazards Models for SCD
# Author: Alexia Sampri
# Dataset assumptions:
#   - df_cleaned: individual-level dataset
#   - Time_zero_Ym: date of follow-up start (format: Date or Year-Month)
#   - Survival_time: follow-up duration in months
#   - Status: 0 = censored, 1 = SCD, 2 = non-SCD death
#   - LVEF_std: numeric (%) left ventricular ejection fraction
#   - Age: baseline age
#############################################################

library(dplyr)
library(lubridate)
library(survival)
library(ggplot2)
library(cmprsk)

# =========================================================
# 1. Data preparation
# =========================================================
data_dir   <- "T:/PROFID/data/processed"
output_dir <- "T:/PROFID/output/"
dir.create(output_dir, showWarnings = FALSE)

setwd(data_dir)

# Load cleaned cohort
df <- readRDS("df_cleaned.rds")

df <- df %>%
  mutate(
    # Convert follow-up to years
    fu_yrs = Survival_time,
    
    # Define event indicator for cause-specific Cox (SCD only)
    status_scd = ifelse(Status == 1, 1, 0),
    
    # Extract month of follow-up start (1–12)
    month_start = month(Time_zero_Ym),
    
    # Define season categories
    season = case_when(
      month_start %in% c(12, 1, 2)  ~ "Winter",
      month_start %in% c(3, 4, 5)   ~ "Spring",
      month_start %in% c(6, 7, 8)   ~ "Summer",
      month_start %in% c(9, 10, 11) ~ "Autumn"
    ) %>% factor(levels = c("Winter", "Spring", "Summer", "Autumn")),
    
    # Continuous seasonality terms (for sine/cosine models)
    sin_time = sin(2 * pi * month_start / 12),
    cos_time = cos(2 * pi * month_start / 12)
  )

# =========================================================
# 2. Base Model — Season as Categorical (Winter = Reference)
# =========================================================

cox_base <- coxph(
  Surv(fu_yrs, status_scd) ~ season,
  data = df
)
summary(cox_base)

# Interpretation:
# HR > 1 → higher risk vs. winter
# HR < 1 → lower risk vs. winter


# =========================================================
# 3. Continuous Time Model — Sine & Cosine Transformation
# =========================================================

cox_sincos <- coxph(
  Surv(fu_yrs, status_scd) ~ sin_time + cos_time,
  data = df
)
summary(cox_sincos)

# Predict relative hazard by month
pred_months <- data.frame(
  month = 1:12,
  sin_time = sin(2 * pi * (1:12) / 12),
  cos_time = cos(2 * pi * (1:12) / 12)
)

pred_months$lin_pred <- predict(cox_sincos, newdata = pred_months, type = "lp")
pred_months$rel_hazard <- exp(pred_months$lin_pred - mean(pred_months$lin_pred))

ggplot(pred_months, aes(x = month, y = rel_hazard)) +
  geom_line(size = 1.2, colour = "#0072B2") +
  theme_minimal(base_size = 14) +
  labs(
    title = "Predicted Seasonal Hazard of SCD (Continuous Time Model)",
    x = "Month",
    y = "Relative Hazard"
  )

ggsave("T:/PROFID/output/survival_seasonal/seasonal_continuous.png", plot = p)

# =========================================================
# 4. Time-Varying Coefficients — Seasonal Effect Changes Over Follow-up
# =========================================================

# Define follow-up duration bands
df <- df %>%
  mutate(
    fu_band = cut(fu_yrs, breaks = c(0, 2, 5, Inf),
                  labels = c("0–2 years", "2–5 years", ">5 years"),
                  right = FALSE)
  )

cox_tv <- coxph(
  Surv(fu_yrs, status_scd) ~ (sin_time + cos_time) * fu_band,
  data = df
)
summary(cox_tv)

# Test improvement vs. simple sine/cosine model
anova(cox_sincos, cox_tv, test = "LRT")

# Visualize predicted seasonal hazard by follow-up band
pred_months_tv <- expand.grid(
  month_start = 1:12,
  fu_band = levels(df$fu_band)
)


# Add sine and cosine transformations for each month
pred_months_tv <- pred_months_tv %>%
  mutate(
    sin_time = sin(2 * pi * month / 12),
    cos_time = cos(2 * pi * month / 12)
  )

# ---- Predict linear predictor from Cox model ----
# Now sin_time and cos_time exist inside pred_months_tv
pred_months_tv$lin_pred <- predict(cox_tv, newdata = pred_months_tv, type = "lp")

# Compute relative hazard within each follow-up band (centre around mean)
pred_months_tv <- pred_months_tv %>%
  group_by(fu_band) %>%
  mutate(rel_hazard = exp(lin_pred - mean(lin_pred))) %>%
  ungroup()

ggplot(pred_months_tv, aes(x = month_start, y = rel_hazard, colour = fu_band)) +
  geom_line(size = 1.1) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    title = "Time-Varying Seasonal Hazard of SCD by Follow-up Duration",
    x = "Month",
    y = "Relative Hazard",
    colour = "Follow-up Band"
  )


# =========================================================
# 5. Stratified Analysis — Separate Models by Risk Subgroups
# =========================================================

library(purrr)
library(broom)

# Exclude missing LVEF before grouping
df <- df %>%
  mutate(
    LVEF_cat = case_when(
      LVEF_std <= 35 ~ "≤35%",
      LVEF_std > 35 & LVEF_std <= 50 ~ "36–50%",
      LVEF_std > 50 ~ ">50%"
    ) %>% factor(levels = c("≤35%", "36–50%", ">50%"))
  ) %>%
  filter(!is.na(LVEF_cat))  # 🔹 Remove missing categories


# Fit separate Cox models for each LVEF group
cox_by_LVEF <- df %>%
  group_split(LVEF_cat) %>%
  set_names(levels(df$LVEF_cat)) %>%
  map(~ coxph(Surv(fu_yrs, status_scd) ~ sin_time + cos_time, data = .x))

# Summaries (optional)
cox_LVEF_results <- map_dfr(cox_by_LVEF, tidy, .id = "LVEF_group")
print(cox_LVEF_results)

# ---- Build prediction grid ----
pred_df <- expand.grid(
  month = 1:12,
  LVEF_cat = levels(df$LVEF_cat)
) %>%
  mutate(
    sin_time = sin(2 * pi * month / 12),
    cos_time = cos(2 * pi * month / 12)
  )

# ---- Predict linear predictors safely ----
pred_df <- pred_df %>%
  group_by(LVEF_cat) %>%
  mutate(
    lin_pred = predict(
      cox_by_LVEF[[unique(LVEF_cat)]],
      newdata = cur_data(),   # use only that subgroup’s grid
      type = "lp"
    ),
    rel_hazard = exp(lin_pred - mean(lin_pred))
  ) %>%
  ungroup()

# ---- Plot ----
ggplot(pred_df, aes(x = month, y = rel_hazard, colour = LVEF_cat)) +
  geom_line(size = 1.2) +
  scale_x_continuous(breaks = 1:12) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Seasonal Hazard of SCD by LVEF Category",
    x = "Month",
    y = "Relative Hazard",
    colour = "LVEF Category"
  )


df <- df %>%
  mutate(
    Age_cat = case_when(
      Age < 55 ~ "<55",
      Age >= 55 & Age < 70 ~ "55-69",
      Age >= 70 ~ ">=70"
    ),
    Age_cat = factor(Age_cat, levels = c("<55", "55-69", ">=70"))
  )

# Fit separate Cox models for each LVEF group
cox_by_age <- df %>%
  group_split(Age_cat) %>%
  set_names(levels(df$Age_cat)) %>%
  map(~ coxph(Surv(fu_yrs, status_scd) ~ sin_time + cos_time, data = .x))

# Summaries (optional)
cox_age_results <- map_dfr(cox_by_age, tidy, .id = "Age_cat")
print(cox_age_results)

# ---- Build prediction grid ----
pred_df <- expand.grid(
  month = 1:12,
  Age_cat = levels(df$Age_cat)
) %>%
  mutate(
    sin_time = sin(2 * pi * month / 12),
    cos_time = cos(2 * pi * month / 12)
  )

# ---- Predict linear predictors safely ----
pred_df <- pred_df %>%
  group_by(Age_cat) %>%
  mutate(
    lin_pred = predict(
      cox_by_age[[unique(Age_cat)]],
      newdata = cur_data(),   # use only that subgroup’s grid
      type = "lp"
    ),
    rel_hazard = exp(lin_pred - mean(lin_pred))
  ) %>%
  ungroup()

# ---- Plot ----
ggplot(pred_df, aes(x = month, y = rel_hazard, colour = Age_cat)) +
  geom_line(size = 1.2) +
  scale_x_continuous(breaks = 1:12) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Seasonal Hazard of SCD by Age Group",
    x = "Month",
    y = "Relative Hazard",
    colour = "Age group"
  )


#############################################################
# End of survival section: Seasonal Cox models
# #############################################################
# Seasonal variation in SCD risk differed markedly by left ventricular ejection fraction (LVEF). Individuals with severely reduced LVEF (<35%) showed the strongest and clinically meaningful seasonal pattern, with a pronounced summer peak in predicted hazard (~25–30% above average) and a clear winter trough. Patients with moderately reduced LVEF (36–50%) demonstrated a more attenuated seasonal pattern, whereas those with preserved LVEF (>50%) exhibited minimal variation throughout the year. These findings suggest that structural cardiac vulnerability amplifies sensitivity to environmental and physiological seasonal triggers of SCD.
# 
# Strong seasonal effects were observed in younger individuals (<55), with a pronounced summer peak in SCD hazard (~1.7–1.8) and a winter nadir. The 55–69 group showed moderate seasonality, whereas individuals aged ≥70 demonstrated minimal seasonal variation. This suggests that seasonal environmental triggers play a larger role in SCD risk among younger patients, while structural cardiac disease and age-related comorbidities dominate risk in older patients.

