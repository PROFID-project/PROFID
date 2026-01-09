# ============================================================
# Task: Advanced TS — Fit SARIMA(p,d,q)(P,D,Q)[12] to SCD rates
# Data needed: scd_monthly with columns: Time_zero_Ym (Date), rate
#   - rate should be incidence per 1,000 person-years (already done)
# ============================================================
install.packages("forecast")
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(forecast)   # auto.arima, Arima, BoxCox.lambda, checkresiduals
  library(ggplot2)
})

# 0) Ensure monthly continuity & order ---------------------------------
scd_monthly <- scd_monthly %>%
  arrange(month_year)

scd_sarima<- scd_monthly

# Optional: flag months with too little person-time


library(forecast)   # for ts, auto.arima, forecast, autoplot
library(lubridate)  # for year(), month()

# Extract start year & month
start_year  <- year(min(scd_sarima$month_year))
start_month <- month(min(scd_sarima_df$month_year))

# Create monthly time series of counts
y_ts <- ts(
  scd_sarima_df$events,
  start    = c(start_year, start_month),
  frequency = 12               # s = 12 → monthly seasonality
)

plot(y_ts, main = "Monthly SCD Events", ylab = "Events", xlab = "Year")

# Estimate Box-Cox lambda
lambda <- BoxCox.lambda(y_ts)

# Transform series (this is what we actually model)
y_bc <- BoxCox(y_ts, lambda)

set.seed(123)  # for reproducibility of model search

fit_sarima <- auto.arima(
  y_bc,
  seasonal    = TRUE,   # allow seasonal terms (P,D,Q)_12
  D           = 1,      # force one seasonal difference
  max.D       = 1,
  stepwise    = FALSE,  # broader search (slower but better)
  approximation = FALSE # exact likelihood
)

summary(fit_sarima)

checkresiduals(fit_sarima)

fit_sarima_fix <- Arima(
  y_bc,
  order = c(1, 1, 2),     # add AR(1)
  seasonal = c(1, 0, 0),  # keep seasonal part as before
  include.constant = FALSE
)

summary(fit_sarima_fix)
checkresiduals(fit_sarima_fix)


# 1) Build a monthly ts object -----------------------------------------
start_year  <- year(min(scd_monthly$month_year))
start_month <- month(min(scd_monthly$month_year))
y_ts <- ts(scd_monthly$rate, start = c(start_year, start_month), frequency = 12)


# 3) Fit SARIMA via auto.arima (let it search seasonal space) ----------

# y_ts → your monthly time series (frequency = 12).
# seasonal = TRUE → allow search over seasonal ARIMA (P,D,Q)_12.
# stepwise = FALSE → performs an exhaustive search (slower but more accurate).
# approximation = FALSE → uses full maximum likelihood, not a fast approximation.
# lambda = NULL → no Box-Cox transformation; models raw counts directly.
fit <- auto.arima(
  y_ts,
  seasonal = TRUE,
  stepwise = FALSE,    # more thorough search
  approximation = FALSE,
  lambda = NULL        # we already BoxCox-transformed
)

fit

# 4) Diagnostics --------------------------------------------------------
# Residuals should look like white noise; no seasonality left.
checkresiduals(fit)  # ACF, Ljung-Box p-value; p > 0.05 is good

# 5) Forecast next 12 months ----------------------------------------
h <- 12  # forecast 1 year 
fc_bc <- forecast(fit, h = h)

# Back-transform from Box–Cox to original scale
fc <- InvBoxCox(fc_bc$mean, lambda = lambda)
fc_lower <- InvBoxCox(fc_bc$lower[,"95%"], lambda = lambda)
fc_upper <- InvBoxCox(fc_bc$upper[,"95%"], lambda = lambda)

# Build a tidy data frame for plotting/saving
fc_dates <- seq(max(scd_monthly$month_year) %m+% months(1),
                by = "month", length.out = h)
fc_df <- tibble(
  date  = fc_dates,
  mean  = as.numeric(fc),
  lower = as.numeric(fc_lower),
  upper = as.numeric(fc_upper)
)

# ========================================
# Simple SARIMA Forecast Plot
# ========================================

library(forecast)

# Forecast next 12 months
fc_sarima <- forecast(fit_sarima, h = 12)

# Simple built-in plot
p_sarima<- autoplot(fc_sarima) +
  labs(
    title = "Forecast of Monthly SCD Incidence (SARIMA Model)",
    x = "Time (months)",
    y = "Predicted incidence rate per 100,000 person-years"
  ) +
  theme_minimal(base_size = 14)

# 8) Save outputs -------------------------------------------------------
dir.create("./output", showWarnings = FALSE)
saveRDS(fit, "./output/sarima_fit.rds")
write.csv(fc_df, "./output/sarima_forecast_12m.csv", row.names = FALSE)
ggsave("./output/sarima_forecast_plot.png", p_sarima, width = 9, height = 6, dpi = 300)
# ============================================================
