###############################################
# TASK 2.4 – ENVIRONMENTAL LINKAGE + ANALYSIS
# CENTRE-LEVEL DAILY CLIMATE (ROBUST VERSION)
###############################################

library(dplyr)
library(lubridate)
library(stringi)
library(splines)
library(survival)
library(ggplot2)

# ---------------------------
# 0) FILE PATHS
# ---------------------------
path_df   <- "T:/PROFID/data/processed/df_cleaned.rds"
path_clim <- "T:/PROFID/data/raw/era5_city_raw_until2022/climate_all_centers.rds"
path_cent <- "T:/PROFID/data/raw/research_centres_UPDATED.txt"   # has center_id, Center, DB, etc.
path_cert <- "T:/PROFID/data/raw/centreIDSEUCERT.csv"            # has pat_id, ctr_name

# ---------------------------
# 1) LOAD
# ---------------------------
df_cleaned <- readRDS(path_df)

climate_all <- readRDS(path_clim) %>%
  mutate(
    date = as.Date(date),
    center_id = as.integer(center_id)
  )

centres <- read.csv(path_cent, stringsAsFactors = FALSE)

# IMPORTANT: don't rename blindly. Ensure you truly have center_id.
# If your centres file has column "ID" not "center_id", then uncomment:
centres <- centres %>% rename(center_id = ID)

stopifnot("center_id" %in% names(centres))
stopifnot("Center" %in% names(centres))

centres <- centres %>%
  mutate(
    center_id = as.integer(center_id),
    Center_clean = tolower(stri_trans_general(Center, "Latin-ASCII")) |> trimws()
  )

cert_map <- read.csv(path_cert, stringsAsFactors = FALSE) %>%
  mutate(
    pat_id = trimws(as.character(pat_id)),
    # manual fix for the known Munich encoding issue
    ctr_name = case_when(
      ctr_name %in% c("Technische UniversitÃ¤t MÃ¼nchen", "Technische Universität München") ~
        "Technische Universitaet Muenchen",
      TRUE ~ ctr_name
    ),
    ctr_name_clean = tolower(stri_trans_general(ctr_name, "Latin-ASCII")) |> trimws()
  )

# ---------------------------
# 2) EXCLUSIONS (no centre info)
# ---------------------------
exclude_regs <- c("MDRT", "MDII", "SHFT", "DRVT")

df_base <- df_cleaned %>%
  filter(!DB %in% exclude_regs) %>%
  mutate(
    event_date = as.Date(event_date)
  )

# ---------------------------
# 3) ASSIGN centre_id
#    - For non-CERT: use DB->centre only if DB is single-centre in centres table
#    - For CERT: use patient mapping pat_id -> ctr_name -> center_id
# ---------------------------

# 3A) Determine which DBs are single-centre in the centres table
single_centre_dbs <- centres %>%
  count(DB, name = "n_centres") %>%
  filter(n_centres == 1) %>%
  pull(DB)

centres_single <- centres %>%
  filter(DB %in% single_centre_dbs) %>%
  select(DB, center_id) %>%
  distinct()

# 3B) Start with non-CERT patients (safe DB->center_id)
df_noncert <- df_base %>%
  filter(DB != "CERT") %>%
  left_join(centres_single, by = "DB")

# 3C) CERT patients: build CERT_pat_id from ID (CERT_USB00590 -> USB00590)
df_cert <- df_base %>%
  filter(DB == "CERT") %>%
  mutate(
    CERT_pat_id = sub("^CERT_", "", ID),
    CERT_pat_id = trimws(as.character(CERT_pat_id))
  ) %>%
  left_join(cert_map %>% select(pat_id, ctr_name_clean), by = c("CERT_pat_id" = "pat_id")) %>%
  left_join(centres %>% select(center_id, Center_clean), by = c("ctr_name_clean" = "Center_clean"))

# 3D) Combine back
df_env <- bind_rows(df_noncert, df_cert)

# sanity: centre assignment
cat("\nMissing center_id by DB:\n")
print(
  df_env %>%
    group_by(DB) %>%
    summarise(n = n(), missing_center = sum(is.na(center_id)), .groups="drop") %>%
    arrange(desc(missing_center))
)

# If you still have missing center_id for CERT, inspect:
# df_env %>% filter(DB=="CERT", is.na(center_id)) %>% distinct(ctr_name_clean)

# ---------------------------
# 4) JOIN CLIMATE (centre_id + event_date)
# ---------------------------
df_env <- df_env %>%
  left_join(
    climate_all,
    by = c("center_id", "event_date" = "date")
  )

cat("\nMissing climate (temp_mean) by DB:\n")
print(
  df_env %>%
    group_by(DB) %>%
    summarise(n = n(), missing_temp = sum(is.na(temp_mean)), .groups="drop") %>%
    arrange(desc(missing_temp))
)

# Keep an analysis-ready set: define event and restrict to SCD vs alive (0/1)
analysis_df <- df_env %>%
  filter(Status %in% c(0, 1)) %>%
  mutate(event = ifelse(Status == 1, 1, 0)) %>%
  # only keep rows where climate exists for the main analysis
  filter(!is.na(temp_mean), !is.na(pressure), !is.na(daylight), !is.na(sunshine)) %>%
  # optional: remove nonsense survival times
  filter(!is.na(Survival_time) & Survival_time > 0)

# ---------------------------
# 5) MODELS (main + robust SE)
# ---------------------------

# Main: temperature spline, stratified by DB, robust SE clustered by centre
model_temp <- coxph(
  Surv(Survival_time, event) ~ ns(temp_mean, df = 4) + strata(DB) + cluster(center_id),
  data = analysis_df
)

model_pressure <- coxph(
  Surv(Survival_time, event) ~ scale(pressure) + strata(DB) + cluster(center_id),
  data = analysis_df
)

model_daylight <- coxph(
  Surv(Survival_time, event) ~ scale(daylight) + strata(DB) + cluster(center_id),
  data = analysis_df
)

model_sunshine <- coxph(
  Surv(Survival_time, event) ~ scale(sunshine) + strata(DB) + cluster(center_id),
  data = analysis_df
)

print(summary(model_temp))
print(summary(model_pressure))
print(summary(model_daylight))
print(summary(model_sunshine))

# ---------------------------
# 6) PLOTS (publication-style simple curves)
# ---------------------------

# helper for spline curve from model_temp
t1 <- quantile(analysis_df$temp_mean, 0.01, na.rm = TRUE)
t99 <- quantile(analysis_df$temp_mean, 0.99, na.rm = TRUE)
temp_seq <- seq(t1, t99, length.out = 300)

pred_temp <- data.frame(
  temp_mean = temp_seq,
  DB = analysis_df$DB[1],          # any valid stratum label
  center_id = analysis_df$center_id[1]  # any valid centre for predict()
)

pt <- predict(model_temp, newdata = pred_temp, type = "lp", se.fit = TRUE)
pred_temp$HR <- exp(pt$fit)
pred_temp$HR_low <- exp(pt$fit - 1.96 * pt$se.fit)
pred_temp$HR_high <- exp(pt$fit + 1.96 * pt$se.fit)

p_temp <- ggplot(pred_temp, aes(x = temp_mean, y = HR)) +
  geom_ribbon(aes(ymin = HR_low, ymax = HR_high), alpha = 0.25) +
  geom_line(linewidth = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed") +
  labs(
    title = "Temperature–SCD Risk Curve",
    subtitle = "Cox model stratified by registry (DB) with centre-clustered SEs",
    x = "Temperature (°C) on event/censoring date",
    y = "Hazard Ratio"
  ) +
  theme_minimal(base_size = 14)

print(p_temp)

# ---------------------------
# 7) SAVE
# ---------------------------
saveRDS(df_env,        "df_env_linked_dailycentre.rds")
saveRDS(analysis_df,   "analysis_df_env_dailycentre.rds")
saveRDS(model_temp,    "model_temp_dailycentre.rds")
saveRDS(model_pressure,"model_pressure_dailycentre.rds")
saveRDS(model_daylight,"model_daylight_dailycentre.rds")
saveRDS(model_sunshine,"model_sunshine_dailycentre.rds")
