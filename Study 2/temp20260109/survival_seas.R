library(dplyr)
library(lubridate)

df_filter <- df %>%
  filter(Status %in% c(0, 1))


df_filter <- df_filter %>%
  mutate(
    month = month(Time_zero_Ym),  # extract month (1–12)
    season = case_when(
      month %in% c(12, 1, 2)  ~ "Winter",
      month %in% c(3, 4, 5)   ~ "Spring",
      month %in% c(6, 7, 8)   ~ "Summer",
      month %in% c(9, 10, 11) ~ "Autumn"
    ),
    season = factor(season, levels = c("Winter", "Spring", "Summer", "Autumn"))  # Winter = reference
  )


library(survival)

cox_base <- coxph(Surv(Survival_time, Status) ~ season, data = df_filter)
summary(cox_base)

cox.zph(cox_base)

install.packages('survminer')
library(survminer)

fit_season <- survfit(Surv(Survival_time, Status) ~ season, data = df_filter)

ggsurvplot(
  fit_season,
  data = df,
  risk.table = TRUE,
  ggtheme = theme_minimal(),
  palette = "Dark2",
  title = "Survival by Season of Follow-up Start (Winter as Reference)",
  xlab = "Time (months)",
  ylab = "Survival probability"
)


library(dplyr)
library(lubridate)

df_filter <- df_filter%>%
  mutate(
    month = month(Time_zero_Ym),  # extract month number (1–12)
    sin_time = sin(2 * pi * month / 12),
    cos_time = cos(2 * pi * month / 12)
  )

library(survival)

cox_sin_cos <- coxph(Surv(Survival_time, Status) ~ sin_time + cos_time, data = df_filter)
summary(cox_sin_cos)


library(ggplot2)

pred_months <- data.frame(
  month = 1:12,
  sin_time = sin(2 * pi * (1:12) / 12),
  cos_time = cos(2 * pi * (1:12) / 12)
)

pred_months$lin_pred <- predict(cox_sin_cos, newdata = pred_months, type = "lp")

ggplot(pred_months, aes(x = month, y = exp(lin_pred))) +
  geom_line(size = 1.2, color = "#0072B2") +
  theme_minimal() +
  labs(
    title = "Predicted Seasonal Hazard of SCD (Continuous Time Model)",
    x = "Month",
    y = "Relative Hazard"
  )

ggsave(
  filename = "output/Continuous_Seasonal_Hazard_SCD.png",
  plot = last_plot(),
  width = 8, height = 5, dpi = 300
)


## =========================================================
## Time-Varying Coefficients:
## Seasonal effects that change over follow-up duration
## =========================================================

library(dplyr)
library(lubridate)
library(survival)
library(ggplot2)

# 1. Prepare survival time and SCD event indicator -------------------------

df_tv <- df %>%
  # follow-up in YEARS (Survival_time is in months)
  mutate(
    fu_yrs    = Survival_time / 12,
    # SCD event indicator: 1 = SCD, 0 = censored / other death
    status_scd = ifelse(Status == 1, 1, 0)
  )

# 2. Create seasonal (sin/cos) terms based on START month ------------------

df_tv <- df_tv %>%
  mutate(
    month_start = month(Time_zero_Ym),                  # 1–12
    sin_time    = sin(2 * pi * month_start / 12),
    cos_time    = cos(2 * pi * month_start / 12)
  )

# 3. Define follow-up bands (where seasonality might differ) ---------------

df_tv <- df_tv %>%
  mutate(
    fu_band = cut(
      fu_yrs,
      breaks = c(0, 2, 5, Inf),                         # 0-2, 2-5, >5 years
      labels = c("0–2 years", "2–5 years", ">5 years"),
      right  = FALSE
    )
  )

# Optional: check how many patients per band
table(df_tv$fu_band, useNA = "ifany")


# 4. Base continuous-time seasonal model (no time-varying effect) ---------

cox_base <- coxph(
  Surv(fu_yrs, status_scd) ~ sin_time + cos_time,
  data = df_tv
)

summary(cox_base)


# 5. Time-varying seasonal model (interactions with follow-up band) -------

cox_tv <- coxph(
  Surv(fu_yrs, status_scd) ~ 
    sin_time * fu_band +         # sin seasonal effect varies by fu_band
    cos_time * fu_band,          # cos seasonal effect varies by fu_band
  data = df_tv
)

summary(cox_tv)

# 6. Formal test: does seasonality change with follow-up? ------------------

anova(cox_base, cox_tv, test = "LRT")
# If p < 0.05  -> evidence that seasonal effect differs by follow-up band
# If p >= 0.05 -> no strong evidence; constant seasonal effect is adequate


# 7. Visualise predicted seasonal hazard by follow-up band -----------------

# Create a grid of months (1–12) and follow-up bands
pred_months <- expand.grid(
  month_start = 1:12,
  fu_band     = levels(df_tv$fu_band)
)

pred_months <- pred_months %>%
  mutate(
    sin_time = sin(2 * pi * month_start / 12),
    cos_time = cos(2 * pi * month_start / 12)
  )

# Linear predictor from the time-varying model
pred_months$lin_pred <- predict(
  cox_tv,
  newdata = pred_months,
  type    = "lp"
)

# Convert to relative hazard within each follow-up band (centre at band mean)
pred_months <- pred_months %>%
  group_by(fu_band) %>%
  mutate(rel_hazard = exp(lin_pred - mean(lin_pred))) %>%
  ungroup()

# Plot
p_tv <- ggplot(pred_months,
               aes(x = month_start, y = rel_hazard, colour = fu_band)) +
  geom_line(size = 1.1) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    title = "Predicted Seasonal Hazard of SCD\nTime-Varying Seasonal Effects",
    x     = "Month of Year",
    y     = "Relative Hazard",
    colour = "Follow-up band"
  ) +
  theme_minimal(base_size = 14)

print(p_tv)

# Optional: save figure
ggsave("output/Cox_TimeVarying_Seasonal_Hazard.png",
       p_tv, width = 8, height = 5, dpi = 300)

