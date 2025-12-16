###################### LCV — Study 3: Cleaning, Events, Exclusions ######################

rm(list = ls())

library(data.table)
library(readxl)
library(stringr)

###################### 1. PATHS ######################

path_lcv      <- "LCV.xlsx"
path_smallmap <- "metadata/02_small_map.xlsx"

###################### 2. READ RAW DATA ######################

dt_lcv <- as.data.table(read_excel(path_lcv, sheet = "Hopkins_LVSCDstudy"))

###################### 3. RENAME USING SMALL MAP ######################

map_small <- as.data.table(read_excel(path_smallmap, sheet = 1))

map_lcv <- map_small[
  !is.na(LCV) & LCV != "" &
    !is.na(harmonised_name) & harmonised_name != ""
]

rename_map <- map_lcv[, .(original = LCV, new = harmonised_name)]

for (i in seq_len(nrow(rename_map))) {
  old <- rename_map$original[i]
  new <- rename_map$new[i]
  if (old %in% names(dt_lcv)) setnames(dt_lcv, old, new)
}

###################### 4. DEATH TYPE & STATUS ######################

dt_lcv[, death_type :=
         fifelse(arrhdeath == 1, "sudden cardiac death",
                 fifelse(hfdeath == 1, "non-sudden cardiac death",
                         "unknown"))]

dt_lcv[, status := fifelse(death_flag == 1, 1, 0)]

###################### 5. COERCE DURATION VARIABLES ######################

duration_vars <- c(
  "days_to_death",
  "time_to_app_therapy",
  "time_to_inap_therapy",
  "days_to_app_shock",
  "days_to_inapp_shock"
)

for (v in duration_vars) {
  if (v %in% names(dt_lcv)) {
    dt_lcv[, (v) := suppressWarnings(as.numeric(get(v)))]
  }
}

###################### 6. CONVERT YEAR-BASED DURATIONS TO DAYS ######################
# (LCV reports some durations in years)

year_based_vars <- c(
  "time_to_app_therapy",
  "time_to_inap_therapy",
  "days_to_death"
)

for (v in year_based_vars) {
  if (v %in% names(dt_lcv)) {
    dt_lcv[, (v) := get(v) * 365.25]
  }
}

###################### 7. FOLLOW-UP DEFINITION ######################
# Time from ICD implant to event or censoring (PROFID rule)

dt_lcv[, t_followup_days := pmax(
  days_to_death,
  time_to_app_therapy,
  time_to_inap_therapy,
  na.rm = TRUE
)]

dt_lcv[is.infinite(t_followup_days), t_followup_days := NA_real_]
dt_lcv[t_followup_days < 0, t_followup_days := NA_real_]



###################### 8. STUDY-3 AGE EXCLUSION ######################

dt <- copy(dt_lcv)

if ("age_icd" %in% names(dt)) {
  dt <- dt[age_icd >= 18 | is.na(age_icd)]
}

###################### 9. SELECT STUDY-3 EVENT VARIABLES ######################

event_vars <- c(
  "patient_id",
  "status",
  "age_icd",
  "icd_type",
  "icd_implant_date",
  "inapp_shock_flag",
  "inapp_therapy",
  "Inapp_ATP",
  "days_to_inapp_shock",
  "time_to_inap_therapy",
  "app_shock_flag",
  "app_therapy_flag",
  "app_atp_flag",
  "days_to_app_shock",
  "time_to_app_therapy",
  "death_date",
  "death_type",
  "death_flag",
  "days_to_death",
  "last_fu_days",
  "t_followup_days"
)

event_vars <- event_vars[event_vars %in% names(dt)]
dt_lcv_final <- dt[, ..event_vars]

###################### 10. DESCRIPTIVE SNAPSHOT ######################

desc_overview <- data.table(
  n_total        = nrow(dt_lcv_final),
  n_deaths       = sum(dt_lcv_final$death_flag == 1, na.rm = TRUE),
  median_fu_days = median(dt_lcv_final$t_followup_days, na.rm = TRUE),
  mean_fu_days   = mean(dt_lcv_final$t_followup_days, na.rm = TRUE),
  missing_fu     = sum(is.na(dt_lcv_final$t_followup_days))
)

desc_overview

###################### 11. EXPORT ######################

saveRDS(dt_lcv_final, "lcv_events_clean.rds")
