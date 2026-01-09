############################################################
# Project: PROFID – Exploring Seasonal Patterns in SCD 
#          Incidence Post-Myocardial Infarction
# Author:  Dr Alexia Sampri
#
# Task 2.2 – Time Series Analysis
# --------------------------------
# Data assumptions:
#   df_cleaned.rds contains at least:
#     - Time_zero_Ym   : Date (baseline / follow-up start)
#     - Survival_time  : numeric (months from baseline to event/censoring)
#     - Status         : 0 = alive/censored,
#                        1 = SCD event,
#                        2 = non-SCD death
#
#   SCD events are defined by Status == 1.
#   SCD event date is reconstructed as:
#      event_date = Time_zero_Ym + Survival_time (in months)
#
# Outputs:
#   - Monthly SCD counts
#   - Classical decomposition (trend/seasonal/random)
#   - X-13ARIMA-SEATS seasonal adjustment:
#       * full series
#       * post-2000
#       * "stable" 2005–2018
#   - Spectral (Fourier) analysis:
#       * full series
#       * pre-2000 vs post-2000 comparison
############################################################

## =============================
## 0. Load packages & set paths
## =============================

# Install once if needed:
# install.packages(c("dplyr", "tidyr", "tibble", "lubridate",
#                    "ggplot2", "seasonal", "forecast", "zoo"))

library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)
library(ggplot2)
library(seasonal)
library(forecast)
library(zoo)

data_dir   <- "T:/PROFID/data/processed"
output_dir <- "T:/PROFID/output/"
dir.create(output_dir, showWarnings = FALSE)

setwd(data_dir)

# Load cleaned cohort
df <- readRDS("df_cleaned.rds")


## =========================================================
## 1. Construct SCD event dates and aggregate monthly counts
## =========================================================

# 1.1 Keep SCD events with valid baseline date & survival time
df_scd <- df %>%
  filter(
    Status == 1,
    !is.na(Time_zero_Ym),
    !is.na(Survival_time)
  ) %>%
  mutate(Time_zero_Ym = as.Date(Time_zero_Ym))

# 1.2 Reconstruct event date from baseline + survival time (months)
df_scd <- df_scd %>%
  mutate(
    Survival_time_round = round(Survival_time),     # months
    event_date = Time_zero_Ym %m+% months(Survival_time_round),
    month_year = floor_date(event_date, "month")
  )

# 1.3 Aggregate SCD event counts per month and fill missing months with 0
scd_monthly <- df_scd %>%
  group_by(month_year) %>%
  summarise(events = n(), .groups = "drop") %>%
  complete(
    month_year = seq(min(month_year), max(month_year), by = "month"),
    fill = list(events = 0L)
  ) %>%
  arrange(month_year)

# Quick sanity plot of raw monthly SCD counts
p_raw <- ggplot(scd_monthly, aes(x = month_year, y = events)) +
  geom_line() +
  labs(
    title = "Monthly SCD Event Counts",
    x = "Calendar Month",
    y = "Number of SCD Events"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(output_dir, "Monthly_SCD_events_raw.png"),
  plot = p_raw,
  width = 9, height = 6, dpi = 300
)

# 1.4 Convert to time series object (monthly frequency)
first_month <- min(scd_monthly$month_year)
scd_ts <- ts(
  data = scd_monthly$events,
  start = c(year(first_month), month(first_month)),
  frequency = 12
)


## =========================================================
## 2. Classical Decomposition (moving averages)
##    Separate trend, seasonal, irregular components (overall)
## =========================================================

# 2.1 Apply classical additive decomposition
decomp <- decompose(scd_ts, type = "additive")

# 2.2 Put components into a tidy data frame
scd_decomp_df <- tibble(
  month_year = scd_monthly$month_year,
  observed   = as.numeric(decomp$x),
  trend      = as.numeric(decomp$trend),
  seasonal   = as.numeric(decomp$seasonal),
  random     = as.numeric(decomp$random)
)

scd_decomp_long <- scd_decomp_df %>%
  pivot_longer(
    cols = c(observed, trend, seasonal, random),
    names_to = "component",
    values_to = "value"
  )

# 2.3 ggplot version of decomposition (instead of base plot())
p_decomp <- ggplot(scd_decomp_long,
                   aes(x = month_year, y = value)) +
  geom_line(color = "steelblue") +
  facet_wrap(~ component, ncol = 1, scales = "free_y") +
  labs(
    title = "Classical Decomposition of Monthly SCD Events",
    x = "Calendar Month",
    y = "Value"
  ) +
  theme_minimal(base_size = 12)

ggsave(
  filename = file.path(output_dir, "Classical_Decomposition_SCD.png"),
  plot = p_decomp,
  width = 9, height = 8, dpi = 300
)

# 2.4 Quantify seasonal amplitude (peak vs trough months)
seasonal_df <- scd_decomp_df %>%
  mutate(month = month(month_year, label = TRUE, abbr = TRUE)) %>%
  group_by(month) %>%
  summarise(
    mean_seasonal = mean(seasonal, na.rm = TRUE),
    sd_seasonal   = sd(seasonal, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(desc(mean_seasonal))

peak_month   <- seasonal_df$month[which.max(seasonal_df$mean_seasonal)]
trough_month <- seasonal_df$month[which.min(seasonal_df$mean_seasonal)]
peak_val     <- max(seasonal_df$mean_seasonal, na.rm = TRUE)
trough_val   <- min(seasonal_df$mean_seasonal, na.rm = TRUE)
amplitude    <- peak_val - trough_val
mean_rate    <- mean(scd_decomp_df$observed, na.rm = TRUE)
relative_amplitude <- (amplitude / mean_rate) * 100

cat("\n==== Seasonal amplitude summary (classical decomposition) ====\n")
cat("Peak month:   ", as.character(peak_month),
    " (mean seasonal = ", round(peak_val, 2), ")\n", sep = "")
cat("Trough month: ", as.character(trough_month),
    " (mean seasonal = ", round(trough_val, 2), ")\n", sep = "")
cat("Absolute amplitude: ", round(amplitude, 2), "\n", sep = "")
cat("Relative amplitude: ", round(relative_amplitude, 2),
    "% of mean monthly SCD count\n\n", sep = "")

# 2.5 Plot average seasonal profile
p_seasonal <- ggplot(seasonal_df,
                     aes(x = month, y = mean_seasonal, group = 1)) +
  geom_line(color = "steelblue", linewidth = 1) +
  geom_point(color = "darkred", size = 2) +
  labs(
    title = "Average Seasonal Effect in Monthly SCD Events",
    subtitle = "From classical decomposition",
    x = "Month",
    y = "Seasonal component (deviation from trend)"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(output_dir, "Average_Monthly_Seasonal_Component.png"),
  plot = p_seasonal,
  width = 8, height = 5, dpi = 300
)

saveRDS(scd_decomp_df,  file.path(output_dir, "scd_decomp_df.rds"))
saveRDS(scd_decomp_long,file.path(output_dir, "scd_decomp_long.rds"))
saveRDS(seasonal_df,    file.path(output_dir, "seasonal_df.rds"))
saveRDS(
  list(
    peak_month = peak_month,
    trough_month = trough_month,
    amplitude = amplitude,
    relative_amplitude = relative_amplitude
  ),
  file = file.path(output_dir, "seasonal_amplitude_summary.rds")
)


## =========================================================
## 3. X-13ARIMA-SEATS Seasonal Adjustment – multiple windows
##    Full / post-2000 / stable 2005–2018
## =========================================================

# 3.0 Clean TS if needed (interpolate tiny gaps)
scd_ts_clean <- ts(
  na.interp(scd_ts),
  start = start(scd_ts),
  frequency = 12
)

# 3.1 Build windows
scd_full      <- scd_ts_clean
scd_post2000  <- window(scd_ts_clean, start = c(2000, 1))
scd_stable    <- window(scd_ts_clean, start = c(2005, 1), end = c(2018, 12))

# Add small constant to avoid log(0) issues
scd_full     <- scd_full + 0.1
scd_post2000 <- scd_post2000 + 0.1
scd_stable   <- scd_stable + 0.1

# 3.2 Helper function to run X-13 safely
run_x13_safe <- function(series, label) {
  tryCatch({
    seas(
      x = series,
      transform.function = "none",
      regression.variables = "td",
      outlier.types = "ao",
      outlier.critical = 4.0,
      x11 = "",
      regression.aictest = NULL
    )
  },
  error = function(e) {
    message(paste("⚠️ X-13 failed for", label, ":", e$message))
    return(NULL)
  })
}

x13_full     <- run_x13_safe(scd_full, "Full")
x13_post2000 <- run_x13_safe(scd_post2000, "Post-2000")
x13_stable   <- run_x13_safe(scd_stable, "Stable (2005–2018)")

# 3.3 Helper to extract adjusted series with dates
get_adj_df <- function(model, label) {
  if (is.null(model)) return(NULL)
  adj <- final(model)
  s   <- start(adj)
  dates <- seq(
    from = as.Date(sprintf("%d-%02d-01", s[1], s[2])),
    by   = "month",
    length.out = length(adj)
  )
  tibble(
    Date = dates,
    Adjusted = as.numeric(adj),
    Model = label
  )
}

df_full     <- get_adj_df(x13_full, "Full")
df_post2000 <- get_adj_df(x13_post2000, "Post-2000")
df_stable   <- get_adj_df(x13_stable, "Stable (2005–2018)")

df_x13_all <- bind_rows(df_full, df_post2000, df_stable)

# 3.4 Overlay plot of adjusted series
p_x13_multi <- ggplot(df_x13_all, aes(x = Date, y = Adjusted, color = Model)) +
  geom_line(linewidth = 1) +
  scale_color_manual(values = c(
    "Full" = "grey50",
    "Post-2000" = "steelblue",
    "Stable (2005–2018)" = "darkred"
  )) +
  labs(
    title = "Seasonally Adjusted Monthly SCD Events (X-13ARIMA-SEATS)",
    subtitle = "Comparison of full, post-2000, and 2005–2018 stable period",
    x = "Calendar Month",
    y = "Adjusted Events"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

ggsave(
  filename = file.path(output_dir, "X13_Sensitivity_Comparison.png"),
  plot = p_x13_multi,
  width = 9, height = 6, dpi = 300
)

# Save models + combined DF
saveRDS(x13_full,     file.path(output_dir, "x13_full_model.rds"))
saveRDS(x13_post2000, file.path(output_dir, "x13_post2000_model.rds"))
saveRDS(x13_stable,   file.path(output_dir, "x13_stable_model.rds"))
saveRDS(df_x13_all,   file.path(output_dir, "x13_adjusted_all.rds"))


## =========================================================
## 4. Spectral (Fourier) Analysis – full + pre/post-2000
## =========================================================

# 4.1 Full series spectral analysis
spec_full <- spec.pgram(
  scd_ts_clean,
  log = "no",
  taper = 0.1,
  plot = FALSE
)

spec_full_df <- tibble(
  frequency = spec_full$freq * 12,   # cycles per year
  spectral_density = spec_full$spec,
  period = "Full"
)

peak_full_idx  <- which.max(spec_full_df$spectral_density)
peak_full_freq <- spec_full_df$frequency[peak_full_idx]

cat("Dominant frequency (FULL): ",
    round(peak_full_freq, 3), " cycles/year\n", sep = "")

# 4.2 Pre-2000 and Post-2000 windows
scd_pre2000  <- window(scd_ts_clean, end   = c(1999, 12))
scd_post2000 <- window(scd_ts_clean, start = c(2000, 1))

spec_pre  <- spec.pgram(scd_pre2000,  log = "no", taper = 0.1, plot = FALSE)
spec_post <- spec.pgram(scd_post2000, log = "no", taper = 0.1, plot = FALSE)

spec_pre_df <- tibble(
  frequency = spec_pre$freq * 12,
  spectral_density = spec_pre$spec,
  period = "Pre-2000"
)

spec_post_df <- tibble(
  frequency = spec_post$freq * 12,
  spectral_density = spec_post$spec,
  period = "Post-2000"
)

# Identify dominant frequencies
peak_pre_idx   <- which.max(spec_pre_df$spectral_density)
peak_post_idx  <- which.max(spec_post_df$spectral_density)
peak_pre_freq  <- spec_pre_df$frequency[peak_pre_idx]
peak_post_freq <- spec_post_df$frequency[peak_post_idx]

cat("Dominant frequency (Pre-2000):  ",
    round(peak_pre_freq, 3), " cycles/year\n", sep = "")
cat("Dominant frequency (Post-2000): ",
    round(peak_post_freq, 3), " cycles/year\n", sep = "")

# Peak strength ratio (post vs pre)
peak_strength_ratio <- max(spec_post_df$spectral_density) /
  max(spec_pre_df$spectral_density)

cat("Peak strength ratio (Post / Pre): ",
    round(peak_strength_ratio, 3), "\n", sep = "")

# 4.3 Combine and plot spectral curves for comparison
spec_combined <- bind_rows(spec_pre_df, spec_post_df)

p_spec_comp <- ggplot(spec_combined,
                      aes(x = frequency, y = spectral_density, color = period)) +
  geom_line(linewidth = 1) +
  labs(
    title = "Spectral Comparison of Monthly SCD Events",
    subtitle = "Pre-2000 versus Post-2000",
    x = "Frequency (cycles/year)",
    y = "Spectral Density"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "top")

ggsave(
  filename = file.path(output_dir, "Spectral_Comparison_PrePost2000.png"),
  plot = p_spec_comp,
  width = 8, height = 5, dpi = 300
)

# 4.4 Also save full-series spectral plot with annotated peak
p_spec_full <- ggplot(spec_full_df,
                      aes(x = frequency, y = spectral_density)) +
  geom_line(color = "steelblue") +
  geom_vline(xintercept = peak_full_freq,
             linetype = "dashed", color = "red") +
  annotate(
    "text",
    x = peak_full_freq,
    y = max(spec_full_df$spectral_density) * 0.9,
    label = paste0("Peak ≈ ", round(peak_full_freq, 2), " cycles/year"),
    hjust = 0, color = "red", size = 3.5
  ) +
  labs(
    title = "Spectral Density of Monthly SCD Events (Full Series)",
    x = "Frequency (cycles/year)",
    y = "Spectral Density"
  ) +
  theme_minimal(base_size = 13)

ggsave(
  filename = file.path(output_dir, "Spectral_Density_Full_SCD.png"),
  plot = p_spec_full,
  width = 8, height = 5, dpi = 300
)

# Save spectral data and peaks
saveRDS(spec_full_df, file.path(output_dir, "spectral_full_df.rds"))
saveRDS(spec_pre_df,  file.path(output_dir, "spectral_pre2000_df.rds"))
saveRDS(spec_post_df, file.path(output_dir, "spectral_post2000_df.rds"))
saveRDS(
  list(
    peak_full_freq  = peak_full_freq,
    peak_pre_freq   = peak_pre_freq,
    peak_post_freq  = peak_post_freq,
    peak_strength_ratio = peak_strength_ratio
  ),
  file.path(output_dir, "spectral_peak_summary.rds")
)


## =========================================================
## 5. Save monthly series & workspace snapshot
## =========================================================

write.csv(scd_monthly,
          file.path(output_dir, "scd_monthly_counts.csv"),
          row.names = FALSE)
saveRDS(scd_monthly,
        file.path(output_dir, "scd_monthly_counts.rds"))

save.image(file = file.path(output_dir, "Timeseries_Analysis_Workspace.RData"))

cat("\n✅ Time series analysis (Task 2.2) completed and all outputs saved in:\n",
    output_dir, "\n")
