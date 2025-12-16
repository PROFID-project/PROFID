###################### ISRAEL — Study 3: Cleaning, Events, Exclusions & Descriptives ######################

rm(list = ls())

library(data.table)
library(readxl)
library(stringr)
library(lubridate)

###################### 1. PATHS ######################

path_israel   <- "israeli.csv"
path_smallmap <- "metadata/02_small_map.xlsx"

dt_israel <- fread(path_israel)

###################### 2. RENAME USING SMALL MAP ######################

map_small <- as.data.table(read_excel(path_smallmap, sheet = 1))

map_israel <- map_small[
  !is.na(ISRAEL) & ISRAEL != "" &
    !is.na(harmonised_name) & harmonised_name != ""
]

rename_map <- map_israel[, .(original = ISRAEL, new = harmonised_name)]

for (i in seq_len(nrow(rename_map))) {
  old <- rename_map$original[i]
  new <- rename_map$new[i]
  if (old %in% names(dt_israel)) setnames(dt_israel, old, new)
}

###################### 3. ADD DEATH FLAG ######################

dt_israel[, death_flag :=
            ifelse(is.na(death_date) | death_date == "" | str_trim(death_date) == "",
                   "no", "yes")]

###################### 4. DATE PARSING ######################

dt_israel[, icd_implant_date := suppressWarnings(
  as.IDate(icd_implant_date, format = "%b %d, %Y")
)]
dt_israel[is.na(icd_implant_date),
          icd_implant_date := suppressWarnings(
            as.IDate(icd_implant_date, format = "%d-%b-%y")
          )]

dt_israel[, icd_implant_date := as.IDate(
  ifelse(!is.na(year(icd_implant_date)) & year(icd_implant_date) < 2000,
         paste0("20", format(icd_implant_date, "%y-%m-%d")),
         format(icd_implant_date, "%Y-%m-%d"))
)]

dt_israel[, death_date := suppressWarnings(
  as.IDate(death_date, format = "%d-%b-%y")
)]
dt_israel[is.na(death_date),
          death_date := suppressWarnings(
            as.IDate(death_date, format = "%b %d, %Y")
          )]

dt_israel[, death_date := as.IDate(
  ifelse(!is.na(year(death_date)) & year(death_date) < 2000,
         paste0("20", format(death_date, "%y-%m-%d")),
         format(death_date, "%Y-%m-%d"))
)]

dt_israel[, last_fu_date := suppressWarnings(
  as.IDate(last_fu_date, format = "%b %d, %Y")
)]
dt_israel[is.na(last_fu_date),
          last_fu_date := suppressWarnings(
            as.IDate(last_fu_date, format = "%d-%b-%y")
          )]

dt_israel[, last_fu_date := as.IDate(
  ifelse(!is.na(year(last_fu_date)) & year(last_fu_date) < 2000,
         paste0("20", format(last_fu_date, "%y-%m-%d")),
         format(last_fu_date, "%Y-%m-%d"))
)]

###################### 5. TIME VARIABLES ######################

dt_israel[, days_to_death := as.numeric(death_date - icd_implant_date)]

dt_israel[death_flag == "no",
          days_to_death := as.numeric(last_fu_date - icd_implant_date)]

dt_israel[days_to_death < 0, days_to_death := NA]

dt_israel[, last_fu_days := as.numeric(last_fu_date - icd_implant_date)]

###################### 6. STATUS ######################
# 0 = alive/censored
# 1 = cardiac death
# 2 = non-cardiac death

dt_israel[, status := NA_integer_]

dt_israel[death_flag == "yes" &
            sudden_cardiac_death_flag %chin% c("NO","no","No"),
          status := 2]

dt_israel[death_flag == "yes" &
            sudden_cardiac_death_flag %chin% c("YES","yes","Yes"),
          status := 1]

dt_israel[death_flag == "no", status := 0]

###################### 7. FOLLOW-UP ######################

dt_israel[, t_followup_days :=
            fifelse(death_flag %chin% c("yes","YES","Yes"),
                    days_to_death,
                    last_fu_days)
]

dt_israel[t_followup_days < 0, t_followup_days := NA]

###################### 8. DIAGNOSTICS (UNCHANGED) ######################

dt_israel[
  death_flag == "yes",
  summary(t_followup_days == days_to_death)
]

dt_israel[
  death_flag == "no",
  summary(t_followup_days == last_fu_days)
]

###################### 9. EVENT VARIABLE SELECTION ######################

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

event_vars <- event_vars[event_vars %in% names(dt_israel)]
dt <- dt_israel[, ..event_vars]

###################### 10. STUDY-3 EXCLUSIONS ######################

dt <- dt[!grepl(";", icd_type)]

include_patterns <- c("single","dual","dx","icd_1","icd_2")
exclude_patterns <- c("crt","biv","bi-v","triple","plug","3")

dt <- dt[
  grepl(paste(include_patterns, collapse="|"), tolower(icd_type)) &
    !grepl(paste(exclude_patterns, collapse="|"), tolower(icd_type))
]

if ("age_icd" %in% names(dt)) {
  dt <- dt[age_icd >= 18 | is.na(age_icd)]
}
###################### ISRAEL — Study 3: Cleaning, Events, Exclusions & Descriptives ######################

rm(list = ls())

library(data.table)
library(readxl)
library(stringr)
library(lubridate)

###################### 1. PATHS ######################

path_israel   <- "israeli.csv"
path_smallmap <- "metadata/02_small_map.xlsx"

dt_israel <- fread(path_israel)

###################### 2. RENAME USING SMALL MAP ######################

map_small <- as.data.table(read_excel(path_smallmap, sheet = 1))

map_israel <- map_small[
  !is.na(ISRAEL) & ISRAEL != "" &
    !is.na(harmonised_name) & harmonised_name != ""
]

rename_map <- map_israel[, .(original = ISRAEL, new = harmonised_name)]

for (i in seq_len(nrow(rename_map))) {
  old <- rename_map$original[i]
  new <- rename_map$new[i]
  if (old %in% names(dt_israel)) setnames(dt_israel, old, new)
}

###################### 3. ADD DEATH FLAG ######################

dt_israel[, death_flag :=
            ifelse(is.na(death_date) | death_date == "" | str_trim(death_date) == "",
                   "no", "yes")]

###################### 4. DATE PARSING ######################

dt_israel[, icd_implant_date := suppressWarnings(
  as.IDate(icd_implant_date, format = "%b %d, %Y")
)]
dt_israel[is.na(icd_implant_date),
          icd_implant_date := suppressWarnings(
            as.IDate(icd_implant_date, format = "%d-%b-%y")
          )]

dt_israel[, icd_implant_date := as.IDate(
  ifelse(!is.na(year(icd_implant_date)) & year(icd_implant_date) < 2000,
         paste0("20", format(icd_implant_date, "%y-%m-%d")),
         format(icd_implant_date, "%Y-%m-%d"))
)]

dt_israel[, death_date := suppressWarnings(
  as.IDate(death_date, format = "%d-%b-%y")
)]
dt_israel[is.na(death_date),
          death_date := suppressWarnings(
            as.IDate(death_date, format = "%b %d, %Y")
          )]

dt_israel[, death_date := as.IDate(
  ifelse(!is.na(year(death_date)) & year(death_date) < 2000,
         paste0("20", format(death_date, "%y-%m-%d")),
         format(death_date, "%Y-%m-%d"))
)]

dt_israel[, last_fu_date := suppressWarnings(
  as.IDate(last_fu_date, format = "%b %d, %Y")
)]
dt_israel[is.na(last_fu_date),
          last_fu_date := suppressWarnings(
            as.IDate(last_fu_date, format = "%d-%b-%y")
          )]

dt_israel[, last_fu_date := as.IDate(
  ifelse(!is.na(year(last_fu_date)) & year(last_fu_date) < 2000,
         paste0("20", format(last_fu_date, "%y-%m-%d")),
         format(last_fu_date, "%Y-%m-%d"))
)]

###################### 5. TIME VARIABLES ######################

dt_israel[, days_to_death := as.numeric(death_date - icd_implant_date)]

dt_israel[death_flag == "no",
          days_to_death := as.numeric(last_fu_date - icd_implant_date)]

dt_israel[days_to_death < 0, days_to_death := NA]

dt_israel[, last_fu_days := as.numeric(last_fu_date - icd_implant_date)]

###################### 6. STATUS ######################
# 0 = alive/censored
# 1 = cardiac death
# 2 = non-cardiac death

dt_israel[, status := NA_integer_]

dt_israel[death_flag == "yes" &
            sudden_cardiac_death_flag %chin% c("NO","no","No"),
          status := 2]

dt_israel[death_flag == "yes" &
            sudden_cardiac_death_flag %chin% c("YES","yes","Yes"),
          status := 1]

dt_israel[death_flag == "no", status := 0]

###################### 7. FOLLOW-UP ######################

dt_israel[, t_followup_days :=
            fifelse(death_flag %chin% c("yes","YES","Yes"),
                    days_to_death,
                    last_fu_days)
]

dt_israel[t_followup_days < 0, t_followup_days := NA]

###################### 8. DIAGNOSTICS (UNCHANGED) ######################

dt_israel[
  death_flag == "yes",
  summary(t_followup_days == days_to_death)
]

dt_israel[
  death_flag == "no",
  summary(t_followup_days == last_fu_days)
]

###################### 9. EVENT VARIABLE SELECTION ######################

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

event_vars <- event_vars[event_vars %in% names(dt_israel)]
dt <- dt_israel[, ..event_vars]

###################### 10. STUDY-3 EXCLUSIONS ######################

dt <- dt[!grepl(";", icd_type)]

include_patterns <- c("single","dual","dx","icd_1","icd_2")
exclude_patterns <- c("crt","biv","bi-v","triple","plug","3")

dt <- dt[
  grepl(paste(include_patterns, collapse="|"), tolower(icd_type)) &
    !grepl(paste(exclude_patterns, collapse="|"), tolower(icd_type))
]

if ("age_icd" %in% names(dt)) {
  dt <- dt[age_icd >= 18 | is.na(age_icd)]
}

###################### 11. EXPORT ######################

dt_israel_final <- copy(dt)
saveRDS(dt_israel_final, "israel_events_clean.rds")

###################### DESCRIPTIVE SNAPSHOT ######################

desc_overview <- data.table(
  n_total        = nrow(dt_israel_final),
  n_deaths       = sum(dt_israel_final$death_flag %chin% c("yes","YES","Yes"), na.rm = TRUE),
  n_cardiac      = sum(dt_israel_final$status == 1, na.rm = TRUE),
  n_noncardiac   = sum(dt_israel_final$status == 2, na.rm = TRUE),
  median_fu_days = median(dt_israel_final$t_followup_days, na.rm = TRUE),
  mean_fu_days   = mean(dt_israel_final$t_followup_days, na.rm = TRUE),
  missing_fu     = sum(is.na(dt_israel_final$t_followup_days))
)

desc_overview

########### diagnostic 

## Inappropriate shock vs follow-up
dt_israel_final[
  !is.na(days_to_inapp_shock) & !is.na(last_fu_days),
  summary(days_to_inapp_shock <= last_fu_days)
]

## Inappropriate therapy vs follow-up
dt_israel_final[
  !is.na(time_to_inap_therapy) & !is.na(last_fu_days),
  summary(time_to_inap_therapy <= last_fu_days)
]

## Death vs follow-up
dt_israel_final[
  !is.na(days_to_death) & !is.na(last_fu_days),
  summary(days_to_death <= last_fu_days)
]


# Inappropriate shock violations
dt_israel_final[
  !is.na(days_to_inapp_shock) &
  !is.na(last_fu_days) &
  days_to_inapp_shock > last_fu_days,
  .(patient_id, days_to_inapp_shock, last_fu_days)
]

# Inappropriate therapy violations
dt_israel_final[
  !is.na(time_to_inap_therapy) &
  !is.na(last_fu_days) &
  time_to_inap_therapy > last_fu_days,
  .(patient_id, time_to_inap_therapy, last_fu_days)
]

# Death violations
dt_israel_final[
  !is.na(days_to_death) &
  !is.na(last_fu_days) &
  days_to_death > last_fu_days,
  .(patient_id, days_to_death, last_fu_days)
]

