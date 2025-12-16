###################### STUDY 3 — MERGE ALL EVENT DATASETS INTO BASELINE ICD ######################

rm(list = ls())
# Install required packages (run once)
# install.packages(c(
#   "data.table",
#   "readxl",
#   "stringr",
#   "lubridate"
# ))



library(data.table)
library(readxl)
library(stringr)

###################### 1. LOAD CLEANED EVENT DATASETS ######################

dt_helios <- readRDS("helios_events_clean.rds")
dt_eucert <- readRDS("eucert_events_clean.rds")
dt_israel <- readRDS("israel_events_clean.rds")
dt_prose  <- readRDS("prose_events_clean.rds")
dt_lcv    <- readRDS("lcv_events_clean.rds")

###################### 2. LOAD BASELINE STACKED ICD DATASET ######################

dt_icd <- fread("ICD.csv")   

helios_icd  <- sum(startsWith(dt_icd$ID, "HELS_"), na.rm = TRUE)
eucert_icd  <- sum(startsWith(dt_icd$ID, "CERT_"), na.rm = TRUE)
israel_icd  <- sum(startsWith(dt_icd$ID, "ISRL_"), na.rm = TRUE)
prose_icd   <- sum(startsWith(dt_icd$ID, "PRSI_"), na.rm = TRUE)
lcv_icd     <- sum(startsWith(dt_icd$ID, "PRSL_"), na.rm = TRUE)

c(
  HELIOS = helios_icd,
  EUCERT = eucert_icd,
  ISRAEL = israel_icd,
  PROSE  = prose_icd,
  LCV    = lcv_icd
)

cat("Loaded datasets:\n")
cat("  HELIOS:", nrow(dt_helios), "\n")
cat("  EUCERT:", nrow(dt_eucert), "\n")
cat("  ISRAEL:", nrow(dt_israel), "\n")
cat("  PROSE :", nrow(dt_prose),  "\n")
cat("  LCV   :", nrow(dt_lcv),    "\n")
cat("  ICD baseline:", nrow(dt_icd), "\n")


###################### 3. DEFINE FINAL EVENT VARIABLES ######################

vars_events <- c(
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

###################### 4. EXTRACTION FUNCTION ######################

extract_events <- function(dt) {
  cols_present <- intersect(vars_events, names(dt))
  cols_missing <- setdiff(vars_events, names(dt))
  
  for (mc in cols_missing) dt[, (mc) := NA]
  
  dt[, ..vars_events]
}

###################### 5. APPLY EXTRACTION AND LABEL DATASET ######################

dt_helios_small <- extract_events(dt_helios)[, dataset := "HELIOS"]
dt_eucert_small <- extract_events(dt_eucert)[, dataset := "EUCERT"]
dt_israel_small <- extract_events(dt_israel)[, dataset := "ISRAEL"]
dt_prose_small  <- extract_events(dt_prose)[,  dataset := "PROSE"]
dt_lcv_small    <- extract_events(dt_lcv)[,    dataset := "LCV"]


###################### 6. STACK EVENT DATASETS ######################

dt_events_merged <- rbindlist(
  list(
    dt_helios_small,
    dt_eucert_small,
    dt_israel_small,
    dt_prose_small,
    dt_lcv_small
  ),
  use.names = TRUE,
  fill = TRUE
)

str(dt_lcv$patient_id)

dt_events_merged[, patient_id_full :=
                   fifelse(dataset == "EUCERT", paste0("CERT_", patient_id),
                           fifelse(dataset == "HELIOS", paste0("HELS_", patient_id),
                                   fifelse(dataset == "ISRAEL", paste0("ISRL_", patient_id),
                                           fifelse(dataset == "PROSE",  paste0("PRSI_", patient_id),
                                                   fifelse(dataset == "LCV", paste0("PRSL_", patient_id, "_MRI"),
                                                           NA_character_)))))
]

cat("\nMerged event datasets:", nrow(dt_events_merged), "rows\n")

cat("\n===== FOLLOW-UP QC =====\n")
dt_events_merged[
  ,
  .(
    n = .N,
    missing_fu = sum(is.na(t_followup_days)),
    median_fu_years = median(t_followup_days, na.rm = TRUE) / 365.25
  ),
  by = dataset
]

###################### 7. MERGE WITH BASELINE ICD DATASET ######################

dt_final <- merge(
  dt_icd,
  dt_events_merged,
  by.x = "ID",
  by.y = "patient_id_full",
  all.x = TRUE
)

cat("Final merged Study 3 dataset:", nrow(dt_final), "rows\n")
table(dt_final$dataset, useNA = "ifany")
unique(dt_lcv$patient_id)[1:30]
table(dt_final$DB, dt_final$dataset)

###################### 8. QC — VERIFY MERGE SUCCESS ######################


cat("Patients with no event data:", sum(is.na(dt_final$dataset)), "\n")
print(table(dt_final$dataset, useNA = "ifany"))

dt_final[, age_diff := age_icd - Age]
summary(dt_final$age_diff)

print(head(dt_final[!is.na(dataset), .(ID, dataset, patient_id)], 10))


###################### REMOVE DATASETS QWITHOUT EVENT VARIABLES ##################

dt_final_study3 <- dt_final[dataset %in% c("HELIOS", "EUCERT", "ISRAEL", "PROSE", "LCV")]
nrow(dt_final_study3)
table(dt_final_study3$dataset)


###################### REMOVE LCV–PROSE DUPLICATES ######################

dup_map <- fread("PROSE_LCVcommon participant.csv")
setnames(dup_map, c("PROSE_ID","LCV_ID"))

dt_final_study3 <- dt_final_study3[
  !(patient_id %in% dup_map$LCV_ID)
]

cat("Rows after removing LCV duplicates:", nrow(dt_final_study3), "\n")
table(dt_final_study3$dataset)

###################### 9. SAVE FINAL DATASET ######################

saveRDS(dt_final_study3, "study3_final_merged.rds")
fwrite(dt_final_study3, "study3_final_merged.csv")


with(
  dt_final_study3,
  table(is.na(t_followup_days), dataset)
)

