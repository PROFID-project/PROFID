#############################################################
# Piecewise Exponential Models: Allow baseline hazard to vary by season and calendar period
# #############################################################

df_pe <- df %>%
  mutate(
    # SCD event indicator (1 = SCD, 0 = censored or other death)
    status_scd = ifelse(Status == 1, 1, 0),
    
    # Convert survival time from MONTHS → YEARS
    fu_yrs = Survival_time / 12,
    
    # Ensure follow-up start date is valid
    start_date = as.Date(Time_zero_Ym)
  ) %>%
  filter(!is.na(fu_yrs), fu_yrs > 0, !is.na(start_date))

df_pe <- df_pe %>%
  mutate(
    start_month = month(start_date),
    start_year  = year(start_date),
    
    season = case_when(
      start_month %in% c(12, 1, 2)  ~ "Winter",
      start_month %in% c(3, 4, 5)   ~ "Spring",
      start_month %in% c(6, 7, 8)   ~ "Summer",
      start_month %in% c(9, 10, 11) ~ "Autumn"
    ),
    season = factor(season, levels = c("Winter","Spring","Summer","Autumn")),
    
    calendar_period = cut(
      start_year,
      breaks = c(1995, 2000, 2005, 2010, 2015, 2021),
      labels = c("1995–1999","2000–2004","2005–2009","2010–2014","2015–2020"),
      include.lowest = TRUE
    )
  ) %>%
  filter(!is.na(season), !is.na(calendar_period))

pe_agg <- df_pe %>%
  group_by(season, calendar_period) %>%
  summarise(
    events      = sum(status_scd, na.rm = TRUE),
    person_time = sum(fu_yrs,      na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rate = events / person_time * 1000   # events per 1000 person-years
  )


fit_pe <- glm(
  events ~ season + calendar_period,
  offset = log(person_time),
  family = poisson(link = "log"),
  data   = pe_aggs
)
summary(fit_pe)

library(broom)
pe_results <- tidy(fit_pe, exponentiate = TRUE, conf.int = TRUE)
pe_results

p_pe <- ggplot(pe_agg, aes(x = calendar_period, y = rate, color = season, group = season)) +
  geom_line(size = 1.2) +
  geom_point(size = 2) +
  labs(
    title = "Crude SCD Rates by Season & Calendar Period (Piecewise Exponential Model)",
    x = "Calendar period",
    y = "Rate per 1,000 person-years",
    color = "Season"
  ) +
  theme_minimal(base_size = 14)

# Save to file
ggsave(
  filename = "output/crude_scd_rates_piecewise.png",
  plot = p_pe,
  width = 10,
  height = 6,
  dpi = 300
)



#############################################################
# Spline-Based Approaches: Use penalised splines to model smooth seasonal hazard functions
# #############################################################


library(dplyr)
library(lubridate)

df_spline <- df_cleaned %>%
  mutate(
    Month_spline = month(Time_zero_Ym)   # numeric 1–12
  )

library(survival)

cox_spline <- coxph(
  Surv(Survival_time, Status == 1) ~ pspline(Month_spline),
  data = df_spline
)
summary(cox_spline)

pred_months <- data.frame(
  Month_spline = 1:12
)

pred_months$lin_pred <- predict(
  cox_spline,
  newdata = pred_months,
  type = "lp"
)

# Convert to hazard ratio relative to mean
pred_months$HR <- exp(pred_months$lin_pred - mean(pred_months$lin_pred))

library(ggplot2)

p_spline <- ggplot(pred_months, aes(x = Month_spline, y = HR)) +
  geom_line(size = 1.3, color = "#0072B2") +
  geom_point(size = 2) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    title = "Seasonal Hazard of SCD (Penalised Spline Model)",
    x = "Month of Year",
    y = "Relative Hazard"
  ) +
  theme_minimal(base_size = 15)

print(p_spline)

dir.create("output", showWarnings = FALSE)

ggsave(
  filename = "output/spline_seasonal_hazard.png",
  plot = p_spline,
  width = 10,
  height = 6,
  dpi = 300
)
