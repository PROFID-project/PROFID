###################### PROSE — Study 3: Cleaning, Events, Exclusions ######################

rm(list = ls())

library(data.table)
library(readxl)
library(lubridate)

###################### 1. PATHS ######################

path_prose    <- "prose.xlsx"
path_smallmap <- "metadata/02_small_map.xlsx"

###################### 2. READ PROSE DATA ######################

dt_prose <- as.data.table(read_excel(path_prose, sheet = "Hopkins_PROSEstudy"))

###################### 3. RENAME USING SMALL MAP ######################

map_small <- as.data.table(read_excel(path_smallmap, sheet = 1))

map_prose <- map_small[
  !is.na(PROSE) & PROSE != "" &
    !is.na(harmonised_name) & harmonised_name != ""
]

rename_map <- map_prose[, .(original = PROSE, new = harmonised_name)]

for (i in seq_len(nrow(rename_map))) {
  old <- rename_map$original[i]
  new <- rename_map$new[i]
  if (old %in% names(dt_prose)) setnames(dt_prose, old, new)
}

###################### 4. FOLLOW-UP DEFINITION ######################
# PROSE has no calendar FU date → use max observed event time

dt_prose[, last_fu_days :=
           pmax(
             days_to_death,
             days_to_app_shock,
             days_to_inapp_shock,
             na.rm = TRUE
           )]

dt_prose[, t_followup_days := last_fu_days]

dt_prose[t_followup_days < 0, t_followup_days := NA]



###################### 5. STUDY-3 EXCLUSIONS ######################

keep_exact <- c("SINGLE - Single", "DUAL - Dual")
dt_prose <- dt_prose[icd_type %in% keep_exact]

if ("age_icd" %in% names(dt_prose)) {
  dt_prose <- dt_prose[age_icd >= 18 | is.na(age_icd)]
}

###################### 6. SELECT STUDY-3 EVENT VARIABLES ######################

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
  "inapp_shock_date",
  "time_to_inap_therapy",
  "time_to_inap_atp",
  "app_shock_flag",
  "app_therapy_flag",
  "app_atp_flag",
  "app_shock_date",
  "time_to_ap_atp",
  "days_to_app_shock",
  "time_to_app_therapy",
  "total_ATP",
  "total_app_shock",
  "total_inapp_ATP",
  "total_inapp_shock",
  "death_date",
  "death_type",
  "death_flag",
  "days_to_death",
  "last_fu_date",
  "last_fu_days",
  "t_followup_days",
  "sudden_cardiac_death_flag"
)

event_vars <- event_vars[event_vars %in% names(dt_prose)]
dt_prose_final <- dt_prose[, ..event_vars]

###################### 7. DESCRIPTIVE SNAPSHOT ######################

desc_overview <- data.table(
  n_total        = nrow(dt_prose_final),
  n_deaths       = sum(dt_prose_final$death_flag %chin% c("yes","YES","Yes"), na.rm = TRUE),
  median_fu_days = median(dt_prose_final$t_followup_days, na.rm = TRUE),
  mean_fu_days   = mean(dt_prose_final$t_followup_days, na.rm = TRUE),
  missing_fu     = sum(is.na(dt_prose_final$t_followup_days))
)

desc_overview

###################### 8. EXPORT ######################

saveRDS(dt_prose_final, "prose_events_clean.rds")
