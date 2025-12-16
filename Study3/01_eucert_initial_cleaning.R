###################### EU-CERT — Study 3: Cleaning, Events, Exclusions & Descriptives ######################

rm(list = ls())
library(data.table)
library(readxl)
library(lubridate)

###################### 1. PATHS ######################

path_data     <- "eu-cert-icd.csv"
path_smallmap <- "metadata/02_small_map.xlsx"

dt <- fread(path_data, na.strings = c("", NA))

###################### 2. RENAME USING SMALL MAP (harmonised_name) ######################

map_small <- as.data.table(read_excel(path_smallmap, sheet = 1))

map_eu <- map_small[
  !is.na(`EU-CERT`) & `EU-CERT` != "" &
    !is.na(harmonised_name) & harmonised_name != ""
]

rename_map <- map_eu[, .(original = `EU-CERT`, new = harmonised_name)]

for (i in seq_len(nrow(rename_map))) {
  old <- rename_map$original[i]
  new <- rename_map$new[i]
  if (old %in% names(dt)) setnames(dt, old, new)
}

###################### 3. DATE CLEANING ######################

date_vars <- c(
  "icd_implant_date",
  "inapp_shock_date",
  "app_shock_date",
  "death_date",
  "last_fu_date"
)

for (v in date_vars) {
  if (v %in% names(dt)) {
    x <- dt[[v]]
    
    d1 <- suppressWarnings(as.Date(x, format = "%Y-%m-%d"))
    d2 <- suppressWarnings(as.Date(x, format = "%d/%m/%Y"))
    d3 <- suppressWarnings(as.Date("1899-12-30") + as.numeric(x))
    
    out <- d1
    out[is.na(out)] <- d2[is.na(out)]
    out[is.na(out)] <- d3[is.na(out)]
    
    dt[[v]] <- out
  }
}

###################### 4. TIME-TO-EVENT VARIABLES ######################

dt[, days_to_death := fifelse(
  death_flag == "yes",
  as.numeric(death_date - icd_implant_date),
  as.numeric(last_fu_date - icd_implant_date)
)]

dt[days_to_death < 0, days_to_death := NA]

dt[, days_to_app_shock := as.numeric(app_shock_date - icd_implant_date)]
dt[days_to_app_shock < 0, days_to_app_shock := NA]

dt[, status := NA_integer_]
dt[death_flag == "no", status := 0]
dt[death_flag == "yes" & death_type == "non-cardiac death", status := 2]
dt[death_flag == "yes" &
     (death_type %chin% c(
       "sudden cardiac death",
       "non-sudden cardiac death",
       "unknown"
     ) | is.na(death_type)),
   status := 1]

dt[, last_fu_days := as.numeric(last_fu_date - icd_implant_date)]

###################### 5. EU-CERT EXCLUSIONS ######################

valid_icd <- c("VVI", "DDD")
dt <- dt[icd_type %chin% valid_icd]

if ("pat_implant_age" %in% names(dt)) {
  dt <- dt[pat_implant_age >= 18 | is.na(pat_implant_age)]
}

###################### 5b. FOLLOW-UP DEFINITION ######################

dt[, event_or_censor_date := fifelse(
  death_flag == "yes" & !is.na(death_date),
  death_date,
  last_fu_date
)]

dt[, t_followup_days := as.numeric(event_or_censor_date - icd_implant_date)]
dt[t_followup_days < 0, t_followup_days := NA]

###################### 6. VARIABLE SELECTION FOR MERGING ######################

keep_vars <- c(
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
  "sudden_cardiac_death_flag",
  "event_or_censor_date",
  "t_followup_days"
)

keep_vars <- keep_vars[keep_vars %in% names(dt)]
dt2 <- dt[, ..keep_vars]

###################### 7. EXPORT CLEANED EU-CERT DATA ######################

dt_eucert <- copy(dt2)
saveRDS(dt_eucert, "eucert_events_clean.rds")
