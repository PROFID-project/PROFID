library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)
library(ggplot2)
install.packages('seasonal')
library(seasonal)

setwd("T:/PROFID/data/raw")
# Load the dataset back into R
df <- readRDS("df_cleaned.rds")
# 
# df <- df %>%
#   filter(Status ==1)
# =====================================================
# Task 2.2 – Time Series Analysis: Classical Decomposition
# Objective: Separate trend, seasonal, and irregular components
# =====================================================

# --- Step 1: Aggregate SCD events and person-time by month ---
# We use Time_zero_Ym (follow-up start) to define calendar months

# Step 1. Aggregate monthly SCD events
scd_monthly <- df %>%
  mutate(month_year = floor_date(as.Date(date_SCD), "month")) %>%
  group_by(month_year) %>%
  summarise(events = sum(Status == 1, na.rm = TRUE)) %>%
  arrange(month_year)

scd_monthly <- scd_monthly %>%
  arrange(month_year) %>%
  mutate(
    year = year(month_year),
    month = month(month_year)
  )    

# Step 2. Convert to time series object
scd_ts <- ts(scd_monthly$events,
             start = c(year(min(scd_monthly$month_year)),
                       month(min(scd_monthly$month_year))),
             frequency = 12)  # monthly frequency

# scd_ts <- ts(
#   start = c(min(scd_monthly$year), min(scd_monthly$month)),
#   frequency = 12
# )

# --- Step 3: Apply classical decomposition ---
# Decompose into trend, seasonal, and random components
decomp <- decompose(scd_ts, type = "additive")  # use "multiplicative" if seasonal changes increase with trend



# --- Step 4: Basic decomposition plot ---

# Save as PNG
png("output/Classical_Decomposition_SCD.png",
    width = 900, height = 700, res = 150)

plot(decomp, col = "steelblue")

dev.off()


# Optional: If we want a smoother version
# We can replace decompose() with stl() (more robust and flexible):
scd_stl <- stl(scd_ts, s.window = "periodic")

png("output/Classical_DecompositionPeriodic_SCD.png",
    width = 900, height = 700, res = 150)

plot(scd_stl, main = "STL Decomposition (Locally Weighted Moving Averages)")
dev.off()

# ======================================================
# Quantify Seasonal Amplitude in SCD Incidence
# ======================================================

library(dplyr)
library(lubridate)

scd_decomp_df <- data.frame(
  date = scd_monthly$month_year,
  observed = decomp$x,
  trend = decomp$trend,
  seasonal = decomp$seasonal,
  random = decomp$random
)

# --- Step 1: Extract seasonal component from the decomposition ---
seasonal_df <- scd_decomp_df %>%
  mutate(month = month(date, label = TRUE, abbr = TRUE)) %>%
  group_by(month) %>%
  summarise(
    mean_seasonal = mean(seasonal, na.rm = TRUE),
    sd_seasonal = sd(seasonal, na.rm = TRUE)
  ) %>%
  arrange(desc(mean_seasonal))

# --- Step 2: Identify peak and trough months ---
peak_month <- seasonal_df$month[which.max(seasonal_df$mean_seasonal)]
trough_month <- seasonal_df$month[which.min(seasonal_df$mean_seasonal)]
peak_val <- max(seasonal_df$mean_seasonal)
trough_val <- min(seasonal_df$mean_seasonal)

# --- Step 3: Calculate amplitude and relative amplitude ---
amplitude <- peak_val - trough_val
mean_rate <- mean(scd_decomp_df$observed, na.rm = TRUE)
relative_amplitude <- (amplitude / mean_rate) * 100

# --- Step 4: Print summary ---
cat("📊 Seasonal amplitude summary:\n")
cat("Peak month: ", peak_month, " (mean seasonal = ", round(peak_val, 2), ")\n", sep = "")
cat("Trough month: ", trough_month, " (mean seasonal = ", round(trough_val, 2), ")\n", sep = "")
cat("Absolute amplitude: ", round(amplitude, 2), "\n", sep = "")
cat("Relative amplitude (% of mean SCD rate): ", round(relative_amplitude, 2), "%\n", sep = "")

# --- Step 5: Optional plot: average seasonal pattern by month ---
aver_seas <- ggplot(seasonal_df, aes(x = month, y = mean_seasonal, group = 1)) +
  geom_line(color = "steelblue", size = 1) +
  geom_point(color = "darkred", size = 2) +
  theme_minimal(base_size = 13) +
  labs(
    title = "Average Monthly Seasonal Effect in SCD Incidence",
    subtitle = "Derived from classical decomposition",
    x = "Month",
    y = "Seasonal component (relative deviation)"
  )

# Save as PNG (high quality)
ggsave(
  filename = "output/Average_Monthly_seasonal.png",
  plot = aver_seas,
  width = 9,       # in inches
  height = 6,
  dpi = 300        # 300 dpi = publication quality
)

saveRDS(seasonal_df, "T:/PROFID/output/seasonal_pattern_summary_classicaldecomp.rds")
saveRDS(data.frame(
  peak_month, trough_month, amplitude, relative_amplitude
), "T:/PROFID/output/seasonal_amplitude_summary_classicaldecomp.rds")


# ==============================================================
# Task: X-13ARIMA-SEATS seasonal adjustment for SCD incidence
# ==============================================================

diffs <- diff(scd_monthly$month_year)
table(diffs)

ggplot(scd_monthly, aes(x = month_year, y = events)) +
  geom_line() +
  labs(title = "Monthly SCD Event Counts",
       x = "Month",
       y = "Number of Events") +
  theme_minimal()

# Save as PNG (high quality)
ggsave(
  filename = "output/Monthly_SCDeventcounts.png",
  plot = aver_seas,
  width = 9,       # in inches
  height = 6,
  dpi = 300        # 300 dpi = publication quality
)

# --- Step 1: Prepare the monthly time series of SCD incidence ---
# Ensure your date_sdc is a proper Date and sorted

scd_monthly <- scd_monthly %>%
  mutate(
    month_year = as.Date(month_year),
    events = as.numeric(events)
  ) %>%
  complete(month_year = seq(min(month_year), max(month_year), by = "month"),
           fill = list(events = 0)) %>%
  arrange(month_year)

# Create the time series object
scd_monthly_ts <- ts(
  scd_monthly$events,
  start = c(year(min(scd_monthly$month_year)), month(min(scd_monthly$month_year))),
  frequency = 12
)

plot(scd_monthly_ts, main = "Monthly SCD Events", ylab = "Events", xlab = "Year")


# --- Step 2: Run X-13ARIMA-SEATS adjustment ---
scd_x13 <- seas(
  x = scd_monthly_ts,
  x11 = "",
  regression.aictest = c("td", "easter"),
  transform.function = "auto"
)

# When I run it Too many outliers or sudden jumps → the model tried to detect >80 “additive outliers” (AOs).
# 
# Zeros in your series → transform.function = "auto" or "log" fails, because log(0) is undefined.
# 
# The series is quite noisy with long zero periods (e.g., before 2000 or after 2018), so X-13 struggles to fit.

# Restrict to stable years
scd_stable <- window(scd_monthly_ts, start = c(2005, 1), end = c(2018, 12))

# Add a small constant. This has no meaningful effect on trends but avoids the log(0) issue.
scd_stable_adj <- scd_stable + 0.1

# Simplify the model

scd_x13 <- seas(
  x = scd_stable_adj,
  transform.function = "none",           # skip log transform
  regression.aictest = NULL,             # disable automatic variable selection
  outlier.types = c("ao"),               # only additive outliers
  outlier.critical = 4.0,                # stricter threshold (fewer detected)
  regression.variables = "td",           # trading-day only (skip Easter)
  x11 = ""                               # use X-11 decomposition
)


# --- Step 3: Inspect results ---
summary(scd_x13)
plot(scd_x13)

# Re enable SEATS( model based seasonal decomposition)

scd_x13_seats <- seas(
  x = scd_stable_adj,
  transform.function = "none",
  regression.variables = "td",
  outlier.types = c("ao"),
  outlier.critical = 4.0
)
plot(scd_x13_seats)

# -------------------------------------------------------------
# STEP 1: Load required package
# -------------------------------------------------------------
library(seasonal)
library(ggplot2)

# -------------------------------------------------------------
# STEP 2: Define time windows
# -------------------------------------------------------------
# Full series
scd_full <- scd_monthly_ts

# Post-2000
scd_post2000 <- window(scd_monthly_ts, start = c(2000, 1))

# Stable years only (2005–2018)
scd_stable <- window(scd_monthly_ts, start = c(2005, 1), end = c(2018, 12))

# Add small constant to avoid log(0) issue
scd_full <- scd_full + 0.1
scd_post2000 <- scd_post2000 + 0.1
scd_stable <- scd_stable + 0.1

# -------------------------------------------------------------
# STEP 3: Define a function to run X-13 safely
# -------------------------------------------------------------
run_x13_safe <- function(series, label) {
  tryCatch({
    seas(
      x = series,
      transform.function = "none",         # no log
      regression.variables = "td",         # trading day correction
      outlier.types = c("ao"),             # additive outliers only
      outlier.critical = 4.0,              # reduce number of detected outliers
      x11 = "",                            # use X-11 decomposition
      regression.aictest = NULL
    )
  },
  error = function(e) {
    message(paste("⚠️ X13 failed for", label, ":", e$message))
    return(NULL)
  })
}

# -------------------------------------------------------------
# STEP 4: Run the models
# -------------------------------------------------------------
x13_full <- run_x13_safe(scd_full, "Full")
x13_post2000 <- run_x13_safe(scd_post2000, "Post-2000")
x13_stable <- run_x13_safe(scd_stable, "Stable (2005–2018)")

# -------------------------------------------------------------
# STEP 5: Extract adjusted series (aligned by time)
# -------------------------------------------------------------
get_aligned_series <- function(model, start_year, end_year) {
  if (is.null(model)) return(NULL)
  
  adj <- final(model)
  start_ts <- start(adj)
  end_ts <- end(adj)
  
  # Convert to Date sequence
  dates <- seq(
    from = as.Date(sprintf("%d-%02d-01", start_ts[1], start_ts[2])),
    by = "month",
    length.out = length(adj)
  )
  
  data.frame(Date = dates, Adjusted = as.numeric(adj))
}

# Extract each adjusted series
df_full      <- get_aligned_series(x13_full, 1995, 2020)
df_post2000  <- get_aligned_series(x13_post2000, 2000, 2020)
df_stable    <- get_aligned_series(x13_stable, 2005, 2018)

# Add labels
df_full$Model <- "Full"
df_post2000$Model <- "Post-2000"
df_stable$Model <- "Stable (2005–2018)"

# -------------------------------------------------------------
# STEP 6: Combine and align all by Date
# -------------------------------------------------------------
df_all <- dplyr::bind_rows(df_full, df_post2000, df_stable)

# Ensure full timeline for comparison
date_seq <- seq(min(df_all$Date), max(df_all$Date), by = "month")
df_all <- dplyr::right_join(df_all, data.frame(Date = date_seq), by = "Date")

# -------------------------------------------------------------
# STEP 7: Plot overlapping adjusted series
# -------------------------------------------------------------
library(ggplot2)

ggplot(df_all, aes(x = Date, y = Adjusted, color = Model)) +
  geom_line(size = 1) +
  scale_color_manual(values = c("Full" = "grey50", "Post-2000" = "steelblue", "Stable (2005–2018)" = "darkred")) +
  labs(
    title = "Comparison of Seasonally Adjusted Monthly SCD Events (X-13ARIMA-SEATS)",
    subtitle = "Aligned across full, post-2000, and stable (2005–2018) periods",
    y = "Adjusted Events",
    x = "Year"
  ) +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top")

# -------------------------------------------------------------
# STEP 7: Save results
# -------------------------------------------------------------
dir.create("output", showWarnings = FALSE)
ggsave("output/X13_Sensitivity_Comparison.png", width = 9, height = 6, dpi = 300)

# ==============================================================
# Task: Spectral Analysis - Identify cyclical patterns (Fourier)
# ==============================================================

library(ggplot2)

# --- Step 1: Prepare your cleaned, monthly time series ---
# We'll use the same scd_monthly_ts (from your X13 analysis)
# Ensure it is a 'ts' object with monthly frequency = 12
# If not, quickly recreate:
scd_monthly_ts <- ts(scd_monthly$events,
                     start = c(1995, 1),
                     frequency = 12)

# --- Step 2: Compute the periodogram using spec.pgram() ---
spec <- spec.pgram(
  scd_monthly_ts,
  log = "no",         # no log-transform (raw power)
  taper = 0.1,        # smooth the edges a bit
  main = "Spectral (Fourier) Analysis of Monthly SCD Events",
  xlab = "Cycles per year",
  ylab = "Spectral Density (Power)"
)

# --- Step 3: Prepare the spectral data for plotting ---
# Convert the output to a data frame to detect peaks
spec_df <- data.frame(
  frequency = spec$freq * 12,   # Convert to cycles per year (since frequency = 1/12 per month)
  spectral_density = spec$spec
)

# Find the frequency with maximum spectral power
peak <- spec_df[which.max(spec_df$spectral_density), ]
cat("Dominant frequency:", round(peak$frequency, 2), "cycles per year\n")

# Interpretation hint:
# - 1.0 = annual cycle
# - 2.0 = biannual (every 6 months)
# - <1.0 = multi-year trend

# ===============================================================
# STEP 4. Plot (ggplot version, optional)
# ===============================================================

library(ggplot2)


ggplot(spec_df, aes(x = frequency, y = spectral_density)) +
  geom_line(color = "steelblue") +
  annotate("point", x = peak$frequency, y = max(spec_df$spec), color = "red", size = 3) +
  annotate("text", x = peak$frequency, y = max(spec_df$spec) * 0.9,
           label = paste0("Peak = ", round(peak$frequency, 3), " cycles/year"),
           color = "red", hjust = 0, size = 3.5) +
  labs(
    title = "Spectral Density of Monthly SCD Events",
    x = "Frequency (cycles/year)",
    y = "Spectral Density"
  ) +
  theme_minimal()

spec_combined <- rbind(
  data.frame(freq = spec_pre$freq, spec = spec_pre$spec, period = "Pre-2000"),
  data.frame(freq = spec_post$freq, spec = spec_post$spec, period = "Post-2000")
)

spectral_plot <- ggplot(spec_combined, aes(x = freq, y = spec, color = period)) +
  geom_line(size = 1) +
  geom_vline(xintercept = c(0.208, 0.05), linetype = "dashed", color = c("blue", "red")) +
  annotate("text", x = 0.208, y = max(spec_combined$spec)*0.8, label = "Pre-2000 peak (0.21)", color = "blue") +
  annotate("text", x = 0.05, y = max(spec_combined$spec)*0.6, label = "Post-2000 peak (0.05)", color = "red") +
  labs(title = "Spectral Comparison of SCD Events", x = "Frequency (cycles/year)", y = "Spectral Density") +
  theme_minimal()

# Save as high-resolution PNG
ggsave("output/Spectral_Comparison_SCD.png", width = 8, height = 5, dpi = 300)


# ==============================================================
# Identify and summarise dominant frequencies from spectral analysis
# ==============================================================
library(forecast)
scd_monthly_clean <- ts(
  na.interp(scd_monthly$events),      # fill small gaps
  frequency = 12,
  start = c(1995, 1)
)

# Stabilize variance
scd_log <- log1p(scd_monthly_clean)    # log(1+x) handles zeros
scd_detrend <- diff(scd_log)           # remove long-term trend

scd_pre2000  <- window(scd_detrend, end = c(1999,12))
scd_post2000 <- window(scd_detrend, start = c(2000,1))

par(mfrow=c(1,2))
spec.pgram(scd_pre2000, spans=c(3,3), main="Pre-2000")
spec.pgram(scd_post2000, spans=c(3,3), main="Post-2000")

time_index <- 1:length(scd_monthly_clean)
fit_harmonic <- lm(
  scd_monthly_clean ~ sin(2*pi*time_index/12) + cos(2*pi*time_index/12)
)
summary(fit_harmonic)


library(zoo)
roll_spec <- rollapply(scd_detrend, width = 120, by = 12, FUN = function(x) {
  s <- spec.pgram(x, plot = FALSE)
  s$freq[which.max(s$spec)] * 12
})


plot(roll_spec, type='l', ylab="Dominant Frequency (cycles/year)",
     main="Evolution of Seasonal Frequency (10-year rolling window)")

par(mfrow = c(1, 2))
spec.pgram(scd_pre2000, spans = c(3,3), main = "Pre-2000: dominant freq ≈ 0.21", log = "no")
spec.pgram(scd_post2000, spans = c(3,3), main = "Post-2000: dominant freq ≈ 0.05", log = "no")

spec_pre <- spec.pgram(scd_pre2000, plot = FALSE)
spec_post <- spec.pgram(scd_post2000, plot = FALSE)

peak_strength_ratio <- max(spec_post$spec) / max(spec_pre$spec)
peak_strength_ratio
# 
# 
# To quantify changes in temporal structure, we compared the peak spectral power between periods.
# The ratio of post-2000 to pre-2000 peak spectral strength was 1.77, indicating that although the dominant frequency shifted from approximately 0.21 to 0.05 cycles per year (from ~5-year to ~20-year oscillations), the post-2000 cyclic component explained nearly 80% more variance than the pre-2000 pattern.
# This suggests a transition from irregular, noisy fluctuations in earlier years to a more stable, coherent long-term temporal trend in SCD incidence after 2000.


# ==============================================================
# SAVE ALL TIME SERIES ANALYSIS OUTPUTS
# ==============================================================
# Folder structure: save in your current project (e.g., "./output/")
# or replace with an absolute path if needed.
dir.create("./output", showWarnings = FALSE)

# 1️⃣ Classical Decomposition Results ----------------------------
saveRDS(scd_decomp_df, file = "./output/scd_decomp.rds")
write.csv(scd_decomp_long, "./output/scd_decomp_long.csv", row.names = FALSE)

# 2️⃣ Seasonal Metrics (amplitude, peak, trough, etc.) ------------
saveRDS(seasonal_df, file = "./output/seasonal_metrics.rds")
write.csv(seasonal_df, "./output/seasonal_metrics.csv", row.names = FALSE)

# 3️⃣ X-13ARIMA-SEATS) ------------------------------
saveRDS(scd_x13, file = "./output/scd_x13_object.rds")


# 4️⃣ Spectral (Fourier) Analysis ---------------------------------
saveRDS(spec_df, file = "./output/spec_object.rds")
write.csv(spec_df, "./output/spec_df.csv", row.names = FALSE)
write.csv(top_peaks, "./output/spec_top_peaks.csv", row.names = FALSE)

# 5️⃣ Figures -----------------------------------------------------
# Save ggplot objects if you created them as variables
ggsave("./output/classical_decomposition_plot.png", plot = last_plot(), width = 9, height = 6, dpi = 300)

# 6️⃣ Aggregated Monthly Series -----------------------------------
write.csv(scd_monthly, "./output/scd_monthly.csv", row.names = FALSE)
saveRDS(scd_monthly, "./output/scd_monthly.rds")

# 7️⃣ Workspace snapshot (everything)
save.image(file = "./output/full_workspace_timeseries.RData")

cat("\n✅ All time series analysis outputs saved successfully in './output/' folder.\n")

