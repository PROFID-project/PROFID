# Load required packages
library(dplyr)
library(readr)
library(lubridate)

df_cleaned <- readRDS("T:/PROFID/data/processed/df_cleaned.rds")

# Read the metadata for research centres
research_centres <- read.csv("T:/PROFID/data/raw/research_centres_UPDATED.txt") %>%
  rename(center_id = ID)

# Read climate data
climate_all <- read_rds("T:/PROFID/data/raw/era5_city_raw/climate_all_centres.rds")

# Merge metadata to climate
climate_merged <- climate_all %>%
  left_join(research_centres, by = "center_id")

# Save as RDS if you plan to load it back into R later
saveRDS(climate_merged, file = "climate_merged.rds")

# Check
head(climate_merged)


# # 1. Filter only patients with SCD events
# df_scd <- df_cleaned %>% filter(Status == 1)
# 
# # # 2. Join only the needed patients with the climate data
# # df_scd_env <- df_scd %>%
# #   left_join(climate_merged, by = "DB")


# STEP 1: Filter SCD patients from ASTN
df_astn <- df_cleaned %>%
  filter(Status == 1, DB == "ASTN")

library(purrr)

# For each patient, get climate data from lag 0 to lag 3
df_scd_lagged <- df_astn %>%
  rowwise() %>%
  mutate(
    lag_days = list(seq(event_date - days(3), event_date, by = "day"))
  ) %>%
  unnest(lag_days) %>%
  rename(date = lag_days)

# Merge with climate data for their center
df_scd_climate <- df_scd_lagged %>%
  left_join(climate_astn, by = c("date", "DB"))

df_scd_lag3_avg <- df_scd_climate %>%
  group_by(ID) %>%
  summarise(
    avg_temp_lag0_3 = mean(temp_mean, na.rm = TRUE),
    avg_pressure_lag0_3 = mean(pressure, na.rm = TRUE),
    avg_sunshine_lag0_3 = mean(sunshine, na.rm = TRUE)
  )


df_scd_final <- df_astn %>%
  left_join(df_scd_lag3_avg, by = "ID")

# Check structure
str(df_scd_final)

# Make sure key variables exist
summary(df_scd_final[c("temp_mean", "pressure", "daylight")])


# subset ASTN patients

library(dplyr)
library(lubridate)


astn <- df_cleaned %>% 
  filter(DB == "ASTN") %>%
  mutate(event_date = as.Date(event_date))

# subset climate data for ASTN and calculate scores before joining

climate_astn <- climate_all %>%
  filter(center_id == "1") %>%
  mutate(date = as.Date(date),
         month = month(date)) %>%
  group_by(month) %>%
  mutate(
    temp_zscore = (temp_mean - mean(temp_mean, na.rm = TRUE)) / sd(temp_mean, na.rm = TRUE),
    heat_extreme = temp_mean > quantile(temp_mean, 0.90, na.rm = TRUE),
    cold_extreme = temp_mean < quantile(temp_mean, 0.10, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(date) %>%
  mutate(pressure_change = pressure - lag(pressure))

library(tidyr)
library(purrr)

# Create lagged weather windows (eg -7 to 0 days before event)

astn_lagged <- astn %>%
  rowwise() %>%
  mutate(lag_dates = list(seq(event_date - 7, event_date, by = "day"))) %>%
  unnest(cols = c(lag_dates)) %>%
  left_join(climate_astn, by = c("lag_dates" = "date"))  # joins climate by date

# Aggregate weather exposures for each individual

astn_weather_summary <- astn_lagged %>%
  group_by(ID) %>%
  summarise(
    mean_temp = mean(temp_mean, na.rm = TRUE),
    max_temp = ifelse(all(is.na(temp_max)), NA, max(temp_max, na.rm = TRUE)),
    mean_temp_zscore = mean(temp_zscore, na.rm = TRUE),
    heat_days = sum(heat_extreme, na.rm = TRUE),
    cold_days = sum(cold_extreme, na.rm = TRUE),
    mean_pressure = mean(pressure, na.rm = TRUE),
    mean_pressure_change = mean(abs(pressure_change), na.rm = TRUE),
    mean_daylight = mean(daylight, na.rm = TRUE),
    mean_sunshine = mean(sunshine, na.rm = TRUE)
  )

# merge with ASTN patient dataset

astn_final <- astn %>%
  left_join(astn_weather_summary, by = "ID")

# Check distributions
summary(astn_final[, c("mean_temp", "heat_days", "cold_days", "mean_pressure_change", "mean_daylight")])

library(survival)

astn_final <- astn_final %>%
  mutate(time = as.numeric(event_date - Time_zero_Ym))  # In days


cox_model <- coxph(Surv(time, Status) ~ mean_temp + heat_days + mean_pressure_change + mean_daylight, data = astn_final)
summary(cox_model)

