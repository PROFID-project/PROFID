#############################################################
# 3. Flexible Survival Models
#  - Piecewise exponential models
#  - Penalised spline Cox models
#  - Competing risks (cause-specific + Fine–Gray)
#############################################################

# Just in case:
library(survival)
library(dplyr)
library(lubridate)
library(readr)
library(cmprsk)
library(splines)

data_dir   <- "T:/PROFID/data/processed"
output_dir <- "T:/PROFID/output"
flex_dir <- file.path(output_dir, "flexible_models")
dir.create(flex_dir, showWarnings = FALSE, recursive = TRUE)

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
# 3.1 Piecewise Exponential Models
#     Baseline hazard varies by follow-up interval,
#     season, and calendar period
#############################################################

# Define follow-up intervals in MONTHS
cut_points <- c(12, 24, 60)  # 0–12, 12–24, 24–60, >60

# Split data into piecewise intervals
df_pexp <- survSplit(
  Surv(fu_mo, status_scd) ~ .,
  data  = df,
  cut   = cut_points,
  start = "tstart",
  end   = "tstop",
  event = "event_scd"
)

df_pexp <- df_pexp %>%
  mutate(
    interval = cut(
      tstop,
      breaks = c(0, cut_points, Inf),
      labels = c("0–12", "12–24", "24–60", ">60"),
      right = TRUE
    ),
    person_time = tstop - tstart,
    year_start  = year(Time_zero_Ym)
  )

# Define broad calendar periods using tertiles of year
year_breaks <- quantile(df_pexp$year_start,
                        probs = c(0, 1/3, 2/3, 1),
                        na.rm = TRUE)
df_pexp <- df_pexp %>%
  mutate(
    cal_period = cut(
      year_start,
      breaks = year_breaks,
      include.lowest = TRUE,
      labels = c("Early", "Middle", "Late")
    ),
    cal_period = droplevels(cal_period)
  )

# Piecewise exponential model via Poisson regression
pexp_model <- glm(
  event_scd ~ interval + season + cal_period,
  offset = log(person_time),
  family = poisson,
  data   = df_pexp
)
summary(pexp_model)

pexp_coef <- broom::tidy(pexp_model, exponentiate = TRUE, conf.int = TRUE)
write_csv(pexp_coef, file.path(flex_dir, "piecewise_exponential_coefficients.csv"))

# Optional: predicted rate by season & interval (keeping calendar period at "Middle")
newdata_pexp <- expand.grid(
  interval   = levels(df_pexp$interval),
  season     = levels(df_pexp$season),
  cal_period = "Middle",
  person_time = 1
)

newdata_pexp$rate <- predict(pexp_model, newdata = newdata_pexp,
                             type = "response")

write_csv(newdata_pexp,
          file.path(flex_dir, "piecewise_rates_by_interval_season.csv"))

p_pexp <- ggplot(newdata_pexp,
                 aes(x = interval, y = rate, group = season,
                     colour = season)) +
  geom_line(aes(group = season)) +
  geom_point() +
  theme_minimal(base_size = 14) +
  labs(
    title = "Piecewise Exponential Rates of SCD",
    x = "Follow-up interval (months)",
    y = "Estimated rate (per month, cal. period = Middle)",
    colour = "Season at baseline"
  )

ggsave(file.path(flex_dir, "piecewise_exponential_rates.png"),
       plot = p_pexp, width = 7, height = 5, dpi = 300)


#############################################################
# 3.2 Spline-Based Cox Models
#     Penalised splines for smooth seasonal hazard
#############################################################

# Penalised spline of month_start in Cox model
cox_pspline <- coxph(
  Surv(fu_mo, status_scd) ~ pspline(month_start, df = 4),
  data = df
)
summary(cox_pspline)

# Prediction grid for smooth seasonal curve
pred_spline <- data.frame(
  month_start = seq(1, 12, by = 0.1)
)

pred_spline$lp <- predict(cox_pspline,
                          newdata = pred_spline,
                          type = "lp")
pred_spline$rel_hazard <- exp(pred_spline$lp - mean(pred_spline$lp))

write_csv(pred_spline,
          file.path(flex_dir, "pspline_seasonal_curve.csv"))

p_pspline <- ggplot(pred_spline,
                    aes(x = month_start, y = rel_hazard)) +
  geom_line(size = 1.1) +
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = 1:12) +
  labs(
    title = "Penalised Spline Seasonal Hazard of SCD",
    x = "Month",
    y = "Relative hazard"
  )

ggsave(file.path(flex_dir, "pspline_seasonal_hazard.png"),
       plot = p_pspline, width = 7, height = 5, dpi = 300)


#############################################################
# 3.3 Competing Risks:
#     Cause-specific hazards & Fine–Gray model
#############################################################

# Cause-specific Cox for SCD (already have as cox_base / cox_sincos),
# but here we also model non-SCD death as its own endpoint.

df <- df %>%
  mutate(
    status_nonscd = ifelse(Status == 2, 1, 0)
  )

# Cause-specific Cox for non-SCD death
cox_cs_nonscd <- coxph(
  Surv(fu_mo, status_nonscd) ~ season,
  data = df
)
summary(cox_cs_nonscd)

cox_cs_nonscd_hr <- broom::tidy(cox_cs_nonscd,
                                exponentiate = TRUE,
                                conf.int = TRUE)
write_csv(cox_cs_nonscd_hr,
          file.path(flex_dir, "cause_specific_nonscd_season_hr.csv"))

# Fine–Gray model for SCD with simple season effect
# (manual dummies to keep it fast and stable)
fg_cov <- cbind(
  season_Spring = ifelse(df$season == "Spring", 1, 0),
  season_Summer = ifelse(df$season == "Summer", 1, 0),
  season_Autumn = ifelse(df$season == "Autumn", 1, 0)
)

fg_scd <- crr(
  ftime   = df$fu_mo,
  fstatus = df$Status,  # 1 = SCD, 2 = non-SCD death, 0 = censored
  cov1    = fg_cov
)
summary(fg_scd)

# Extract Fine–Gray coefficients
fg_coef <- data.frame(
  term = rownames(fg_scd$coef),
  estimate = fg_scd$coef,
  se = sqrt(diag(fg_scd$var))
)
fg_coef <- fg_coef %>%
  mutate(
    HR = exp(estimate),
    lower = exp(estimate - 1.96 * se),
    upper = exp(estimate + 1.96 * se)
  )

write_csv(fg_coef,
          file.path(flex_dir, "finegray_season_hr.csv"))

# Optional: cumulative incidence by season (SCD vs non-SCD)
ci_season <- cuminc(
  ftime   = df$fu_mo,
  fstatus = df$Status,
  group   = df$season
)

# Quick plot of CIF for SCD by season (grouped)
# (Convert to data.frame for ggplot)
ci_list <- lapply(names(ci_season), function(nm) {
  out <- ci_season[[nm]]
  data.frame(
    time = out$time,
    est  = out$est,
    group = nm
  )
})
ci_df <- bind_rows(ci_list) %>%
  filter(grepl("^1 ", group)) %>%       # event type 1 = SCD
  mutate(
    season = sub("^1 ", "", group)
  )

p_cif <- ggplot(ci_df,
                aes(x = time, y = est, colour = season)) +
  geom_line(size = 1.1) +
  theme_minimal(base_size = 14) +
  labs(
    title = "Cumulative Incidence of SCD by Season (Fine–Gray framework)",
    x = "Follow-up time (months)",
    y = "Cumulative incidence",
    colour = "Season at baseline"
  )

ggsave(file.path(flex_dir, "cif_scd_by_season.png"),
       plot = p_cif, width = 7, height = 5, dpi = 300)


#############################################################
# ADDITIONAL: Fast Fine–Gray model
#  - keeps only season dummies
#  - runs in seconds
#############################################################

message("Running fast Fine–Gray seasonal model...")

# Prepare vectors
ftime   <- df$fu_mo
fstatus <- df$Status  # 1=SCD, 2=non-SCD, 0=censored

# Minimal dummy covariate matrix
fg_cov_fast <- cbind(
  Spring = as.integer(df$season == "Spring"),
  Summer = as.integer(df$season == "Summer"),
  Autumn = as.integer(df$season == "Autumn")
)

# Fine–Gray competing risk for SCD
fg_scd_fast <- crr(
  ftime     = ftime,
  fstatus   = fstatus,
  cov1      = fg_cov_fast,
  failcode  = 1,   # SCD
  cencode   = 0
)

summary(fg_scd_fast)

# Extract clean HR table
fg_coef_fast <- data.frame(
  term = rownames(fg_scd_fast$coef),
  estimate = fg_scd_fast$coef,
  se = sqrt(diag(fg_scd_fast$var))
) %>%
  mutate(
    HR = exp(estimate),
    lower = exp(estimate - 1.96 * se),
    upper = exp(estimate + 1.96 * se)
  )

write_csv(
  fg_coef_fast,
  file.path(flex_dir, "fast_finegray_season_hr.csv")
)

print("Fast Fine–Gray model saved at:")
print(file.path(flex_dir, "fast_finegray_season_hr.csv"))


#############################################################
# End of Flexible Survival Models section
#############################################################
