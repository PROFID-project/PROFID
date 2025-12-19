library(httr)
library(jsonlite)
library(dplyr)

centers <- read.csv("T:/PROFID/data/raw/research_centres_UPDATED.txt")
View(centers)

get_climate <- function(lat, lon,
                        start = "1994-01-01",
                        end   = "2022-12-31") {
  
  url <- paste0(
    "https://archive-api.open-meteo.com/v1/era5?",
    "latitude=", lat,
    "&longitude=", lon,
    "&start_date=", start,
    "&end_date=", end,
    "&daily=temperature_2m_mean,temperature_2m_min,temperature_2m_max,",
    "surface_pressure_mean,sunshine_duration,daylight_duration",
    "&timezone=UTC"
  )
  
  res <- tryCatch({
    fromJSON(content(GET(url), "text", encoding = "UTF-8"))
  }, error = function(e) NULL)
  
  if (is.null(res) || is.null(res$daily) || is.null(res$daily$time)) {
    return(NULL)
  }
  
  daily <- res$daily
  
  tibble(
    date      = as.Date(daily$time),
    temp_mean = daily$temperature_2m_mean,
    temp_min  = daily$temperature_2m_min,
    temp_max  = daily$temperature_2m_max,
    pressure  = daily$surface_pressure_mean,
    sunshine  = daily$sunshine_duration / 3600,
    daylight  = daily$daylight_duration  / 3600
  )
}


get_climate_safe <- function(lat, lon, id, tries = 3) {
  for (k in 1:tries) {
    res <- get_climate(lat, lon)
    if (!is.null(res)) {
      message("  ✓ success for center ", id)
      return(res)
    }
    message("  ✗ failed attempt ", k, " — retrying after 25 sec...")
    Sys.sleep(25)
  }
  message("  ✗ center ", id, " failed after retries.")
  return(NULL)
}


all_centers <- list()

for (i in 14:nrow(centers)) {
  
  lat <- centers$Latitude[i]
  lon <- centers$Longitude[i]
  id  <- centers$ID[i]
  
  message(">>> Fetching center ", id, " (row ", i, ")")
  
  Sys.sleep(20)   # IMPORTANT → prevent throttling
  
  clim <- get_climate_safe(lat, lon, id)
  
  if (!is.null(clim)) {
    clim$center_id <- id
    saveRDS(clim, paste0("climate_center_", id, ".rds"))
    all_centers[[as.character(id)]] <- clim
  }
}


climate_all <- bind_rows(all_centers)
saveRDS(climate_all, "climate_all_32_centers.rds")

climate_all %>% count(center_id)




library(dplyr)
library(lubridate)

climate_summary <- climate_all %>%
  mutate(
    month = month(date),
    season = case_when(
      month %in% c(12,1,2) ~ "winter",
      month %in% c(3,4,5) ~ "spring",
      month %in% c(6,7,8) ~ "summer",
      month %in% c(9,10,11) ~ "autumn"
    )
  ) %>%
  group_by(center_id) %>%
  summarise(
    
    # Temperature
    mean_annual_temp = mean(temp_mean, na.rm = TRUE),
    mean_winter_temp = mean(temp_mean[season == "winter"], na.rm = TRUE),
    mean_summer_temp = mean(temp_mean[season == "summer"], na.rm = TRUE),
    temp_amp = mean_summer_temp - mean_winter_temp,
    temp_p05 = quantile(temp_mean, 0.05, na.rm = TRUE),
    temp_p95 = quantile(temp_mean, 0.95, na.rm = TRUE),
    
    # Pressure
    mean_pressure = mean(pressure, na.rm = TRUE),
    pressure_sd   = sd(pressure, na.rm = TRUE),
    pressure_p05  = quantile(pressure, 0.05, na.rm = TRUE),
    pressure_p95  = quantile(pressure, 0.95, na.rm = TRUE),
    
    # Daylight
    mean_daylight = mean(daylight, na.rm = TRUE),
    winter_daylight = min(daylight, na.rm = TRUE),
    summer_daylight = max(daylight, na.rm = TRUE),
    daylight_amp = summer_daylight - winter_daylight,
    
    # Sunshine
    mean_sunshine = mean(sunshine, na.rm = TRUE),
    sunshine_sd   = sd(sunshine, na.rm = TRUE),
    
    # Placeholder for air quality (if you add later)
    air_quality = NA
  )

climate_summary <- climate_summary %>%
  rename(center_id = ID)

centers <- centers %>%
  rename(center_id = ID)

climate_centers <- centers %>%
  left_join(climate_summary, by = "center_id")


setdiff(centers$center_id, climate_summary$center_id)

missing_ids <- 32

missing_rows <- which(centers$ID %in% missing_ids)
missing_rows

get_climate_safe <- function(lat, lon, id, tries = 4) {
  for (k in 1:tries) {
    Sys.sleep(25)   # IMPORTANT → prevent rate limiting
    res <- get_climate(lat, lon)
    if (!is.null(res)) {
      message("  ✓ successfully retrieved center ", id)
      return(res)
    }
    message("  ✗ attempt ", k, " failed for center ", id)
  }
  message("  ✗ center ", id, " failed after ", tries, " attempts")
  return(NULL)
}


redownload_list <- list()

for (id in missing_ids) {
  
  idx <- which(centers$ID == id)
  lat <- centers$Latitude[idx]
  lon <- centers$Longitude[idx]
  
  message(">>> Redownloading center ", id, " (row ", idx, ")")
  
  clim <- get_climate_safe(lat, lon, id)
  
  if (!is.null(clim)) {
    clim$center_id <- id
    saveRDS(clim, paste0("climate_center_", id, ".rds"))
    redownload_list[[as.character(id)]] <- clim
  } else {
    message("center ", id, " still empty — IP likely throttled. Reset IP & retry.")
  }
}

# Save as RDS if you plan to load it back into R later
saveRDS(climate_all, file = "climate_all.rds")

# Load necessary packages
library(dplyr)
library(readr)   # For reading RDS files
library(purrr)   # For mapping over files

# Step 1: Set your folder path (update accordingly)
folder_path <- "T:/PROFID/data/raw/era5_city_raw_until2022"  
# Step 2: List all climate RDS files in the folder
climate_files <- list.files(path = folder_path, pattern = "^climate_\\d+\\.rds$", full.names = TRUE)
library(dplyr)
library(purrr)

# path to your folder with RDS files
path <- "T:/PROFID/data/raw/era5_city_raw_until2022"

# list all files matching the naming pattern
files <- list.files(path, pattern = "^climate_center_[0-9]+\\.rds$", full.names = TRUE)

# read and bind
climate_all <- files %>%
  map(readRDS) %>%
  bind_rows()

# save combined file
saveRDS(climate_all, file = file.path(path, "climate_all_centers.rds"))

climate_all
