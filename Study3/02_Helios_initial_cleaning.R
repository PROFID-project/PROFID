###################### HELIOS — Study 3: Cleaning, Events, Exclusions ######################

rm(list = ls())
library(data.table)
library(readxl)

###################### 1. PATHS ######################

path_helios   <- "Helius.xlsx"
path_smallmap <- "metadata/02_small_map.xlsx"

###################### 2. READ AND MERGE HELIOS SHEETS ######################

dt_target   <- as.data.table(read_excel(path_helios, sheet = "target_pop"))
dt_base     <- as.data.table(read_excel(path_helios, sheet = "Baseline Characteristics"))
dt_pmh_icd  <- as.data.table(read_excel(path_helios, sheet = "Past Medical History (ICD10)"))
dt_pmh_ops  <- as.data.table(read_excel(path_helios, sheet = "Past Medical History (OPS)"))
dt_med      <- as.data.table(read_excel(path_helios, sheet = "Medication (ATC)"))
dt_lab      <- as.data.table(read_excel(path_helios, sheet = "Lab Data"))
dt_imaging  <- as.data.table(read_excel(path_helios, sheet = "Imaging"))
dt_cmr      <- as.data.table(read_excel(path_helios, sheet = "CMR-Scar and GZ"))
dt_ecg      <- as.data.table(read_excel(path_helios, sheet = "ECG"))
dt_icd_q    <- as.data.table(read_excel(path_helios, sheet = "ICD queries"))
dt_outcome  <- as.data.table(read_excel(path_helios, sheet = "Outcome"))

merge_list <- list(
  dt_target, dt_base, dt_pmh_icd, dt_pmh_ops,
  dt_med, dt_lab, dt_imaging, dt_cmr,
  dt_ecg, dt_icd_q, dt_outcome
)

dt_full <- Reduce(function(x, y) merge(x, y, by = "PAT_INDEX", all = TRUE), merge_list)

###################### 3. RENAME USING SMALL MAP ######################

map_small <- as.data.table(read_excel(path_smallmap, sheet = 1))

map_helios <- map_small[
  !is.na(Helios) & Helios != "" &
    !is.na(harmonised_name) & harmonised_name != ""
]

rename_map <- map_helios[, .(original = Helios, new = harmonised_name)]

for (i in seq_len(nrow(rename_map))) {
  old <- rename_map$original[i]
  new <- rename_map$new[i]
  if (old %in% names(dt_full)) setnames(dt_full, old, new)
}

###################### 4. FOLLOW-UP DEFINITION (YEAR → DAYS) ######################

dt_full[, implant_year := as.numeric(icd_implant_date)]
dt_full[, death_year   := as.numeric(death_date)]
dt_full[, fu_year      := as.numeric(last_fu_date)]

dt_full[, followup_years := fifelse(
  death_flag %chin% c("Yes","yes","YES"),
  death_year - implant_year,
  fu_year    - implant_year
)]

dt_full[, t_followup_days := followup_years * 365]
dt_full[t_followup_days < 0, t_followup_days := NA]

str(dt_full)
unique(dt_full$icd_type)
###################### 5. STUDY-3 EXCLUSIONS ######################

dt <- copy(dt_full)

# Device type: KEEP ONLY exact ICD_1 (single) and ICD_2 (dual)
dt <- dt[icd_type %in% c("ICD_1", "ICD_2")]

nrow(dt)

# Age ≥ 18 (if available)
if ("age_icd" %in% names(dt)) {
  dt <- dt[age_icd >= 18 | is.na(age_icd)]
}



###################### 6. SELECT STUDY-3 EVENT VARIABLES ######################

event_vars <- c(
  "patient_id",
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
  "followup_years",
  "t_followup_days",
  "sudden_cardiac_death_flag"
)

event_vars <- event_vars[event_vars %in% names(dt)]
dt_helios_final <- dt[, ..event_vars]
###################### 7. DESCRIPTIVE SNAPSHOT ######################

desc_overview <- data.table(
  n_total        = nrow(dt_helios_final),
  n_deaths       = sum(dt_helios_final$death_flag %chin% c("yes","YES","Yes"), na.rm = TRUE),
  median_fu_days = median(dt_helios_final$t_followup_days, na.rm = TRUE),
  mean_fu_days   = mean(dt_helios_final$t_followup_days, na.rm = TRUE),
  missing_fu     = sum(is.na(dt_helios_final$t_followup_days))
)

desc_overview

###################### 8. EXPORT ######################

saveRDS(dt_helios_final, "helios_events_clean.rds")

