############################################################
# PROFID – Seasonal Patterns in Sudden Cardiac Death (SCD)
# Time-Series Pipeline (Counts of SCD over Calendar Time)
#
# Author: Dr Alexia Sampri
# Objective: To investigate temporal variation and seasonal
#            patterns in SCD incidence post-MI, using:
#            - Classical decomposition
#            - X-13ARIMA-SEATS
#            - Spectral analysis
#            - SARIMA models
#            - State-space models (KFAS)
############################################################
############################################################
# Install ALL necessary packages for the full analysis
############################################################

install.packages(c(
  "dplyr",
  "tidyr",
  "lubridate",
  "ggplot2",
  "forecast",
  "seasonal",
  "KFAS",
  "zoo",
  "xts",
  "tsibble",
  "fable",
  "fabletools",
  "lubridate",
  "tibble"
))

## 0. Load libraries ---------------------------------------

library(dplyr)
library(lubridate)
library(tidyr)
library(ggplot2)
library(forecast)    # for auto.arima, forecast, autoplot.ts
install.packages("seasonal") # if not installed
install.packages("KFAS")     # if not installed
library(seasonal)
library(KFAS)

setwd("T:/PROFID/data/processed")
# Load the dataset back into R
df <- readRDS("df_cleaned.rds")
############################################################
# PROJECT: Exploring Seasonal Patterns in SCD Incidence Post-MI
# AUTHOR: Dr Alexia Sampri
# TASK: Full Time Series Analysis Pipeline (All Steps in One Script)
# DATE: YYYY-MM-DD
#
# INPUTS:
#   - df with variables:
#         Time_zero_Ym (Date)
#         Survival_time (months)
#         Status (0=censored, 1=SCD, 2=other death)
#
# OUTPUTS:
#   - Clean monthly SCD series
#   - Classical decomposition
#   - X-13ARIMA-SEATS seasonal adjustment
#   - Spectral (Fourier) analysis
#   - SARIMA model + forecast
#   - State-space model (KFAS)
############################################################


############# 0. LOAD PACKAGES ###########################################

library(dplyr)
library(tidyr)
library(lubridate)
library(ggplot2)

# Time series packages
library(forecast)
library(seasonal)
library(KFAS)



############# 1. DATA PREPARATION ########################################

# Ensure date format is correct
df <- df %>%
  mutate(
    Time_zero_Ym = as.Date(Time_zero_Ym),
    surv_months = round(Survival_time)      # convert survival time to full months
  )

# Compute event date for SCD cases
df <- df %>%
  mutate(
    event_date = if_else(
      Status == 1,
      Time_zero_Ym %m+% months(surv_months),
      as.Date(NA)
    )
  )

# Extract monthly indicator
df <- df %>%
  mutate(month_year = floor_date(event_date, "month"))


############# 2. BUILD MONTHLY SCD TIME SERIES ############################

# Keep ONLY SCD events (Status==1)
scd_events <- df %>% filter(Status == 1, !is.na(event_date))

# Count events per month
scd_monthly <- scd_events %>%
  count(month_year, name="events")

# Build continuous monthly grid (no gaps)
full_months <- tibble(
  month_year = seq.Date(
    from = min(scd_events$event_date, na.rm=TRUE),
    to   = max(scd_events$event_date, na.rm=TRUE),
    by = "month"
  )
)

# Merge & fill missing months with 0 events
scd_monthly_full <- full_months %>%
  left_join(scd_monthly, by="month_year") %>%
  mutate(events = replace_na(events, 0))


############# 3. CREATE TS OBJECT #########################################

scd_ts <- ts(
  scd_monthly_full$events,
  start = c(year(min(scd_monthly_full$month_year)),
            month(min(scd_monthly_full$month_year))),
  frequency = 12
)



############# 4. CLASSICAL SEASONAL DECOMPOSITION ##########################

scd_decomp <- decompose(scd_ts, type="additive")

# Plot
plot(scd_decomp, main="Classical Decomposition of SCD Monthly Series")

# Save
ggsave("plots/scd_classical_decomposition.png", width=8, height=6, dpi=300)



############# 5. X-13ARIMA-SEATS (Advanced Adjustment) #####################

options(seas.warn=FALSE)

scd_x13 <- seas(scd_ts)

# Plot
plot(scd_x13, main="X-13ARIMA-SEATS Seasonal Adjustment")

# Save
ggsave("~/output/figures/timeseries/scd_x13_adjustment.png", width=8, height=6, dpi=300)



############# 6. SPECTRAL / FOURIER ANALYSIS ###############################

spec <- spectrum(scd_ts, plot=FALSE)

# Identify dominant frequency
peak_index <- which.max(spec$spec)
dom_freq <- spec$freq[peak_index]

cat("Dominant frequency =", dom_freq, "\n")
cat("Approx cycle length =", 1/dom_freq, "months\n")

# Plot periodogram
png("~/output/figures/scd_spectral_analysis.png", width=800, height=600)
plot(spec, main="Spectral Analysis of SCD Monthly Counts")
abline(v=dom_freq, col="red", lwd=2)
dev.off()



############# 7. SARIMA MODEL #############################################

fit_sarima <- auto.arima(
  scd_ts,
  seasonal=TRUE,
  stepwise=FALSE,
  approximation=FALSE
)

summary(fit_sarima)

# Forecast 12 months
fc <- forecast(fit_sarima, h=12)

autoplot(fc) +
  ggtitle("SARIMA Forecast for SCD Events") +
  theme_minimal()

ggsave("~/output/figures/scd_sarima_forecast.png", width=8, height=6, dpi=300)



############# 8. STATE SPACE MODEL (KFAS) #################################

ss_model <- SSModel(scd_ts ~ SSMtrend(1) + SSMseasonal(12), H = NA)

fit_ss <- fitSSM(ss_model, inits = rep(0.1, 2))

kfas_out <- KFS(fit_ss$model)

# Extract smoothed trend
trend <- kfas_out$alphahat[,1]

png("~/output/figures/scd_kfas_state_space.png", width=800, height=600)
plot(trend, type="l", lwd=2, col="blue",
     main="State-Space Trend Component (KFAS)",
     ylab="Smoothed Trend")
dev.off()



#######################################################################
#                     END OF FULL SCRIPT
#######################################################################
