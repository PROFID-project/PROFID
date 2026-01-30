###############################################################################

# EU-CERT-ICD preprocessing (PROFID CDM-ready) 

# 1) Load EU-CERT-ICD registry CSV (patient-level)

# 2) Load dictionary + CDM and (a) convert date columns and (b) rename variables once

# 3) Apply exclusions ( Age < 18; exclude CRT-D; exclude overlap centre)

# 4) QC: check follow-up consistency (length_fu_mortality vs lastfu_date - implant_date)

# 5) Create endpoints in DAYS (NO conversion to months):

#    - Death from any cause (Status_death, Time_death) in DAYS

#    - Inappropriate shock (Status_FIS, Time_FIS) in DAYS

# 6) Add dataset metadata (DB, ID_f, Baseline_type)

# 7) Restrict to CDM + key endpoints and save

# 8) Run HZ basic summary (FINAL STEP)

###############################################################################



# =============================================================================

# 1) Libraries

# =============================================================================

install.packages("tidyverse")
library(tidyverse)
library(dplyr)
library(data.table)
library(tools)
install.packages("openxlsx")
library(openxlsx)



# =============================================================================

# 2) Paths / configuration (EDIT)

# =============================================================================

IN_CSV       <-  "//Charite.de/Centren/AG/f-dhzc-profid/datasets/local/eu-cert-icd/data/original/registry_data_eu-cert-icd_selection_161019-Data-sheet.csv"

IN_DICT_XLSX <- "S:/AG/f-dhzc-profid/datasets/local/eu-cert-icd/data/dictionary/eu-cert-icd-data-dictionary-raw.xlsx"

IN_CDM_CSV   <- "//Charite.de/Centren/AG/f-dhzc-profid/datasets/cdm/profid-common-data-model.csv"



source("S:/AG/f-dhzc-profid/datasets/local/eu-cert-icd/scripts/hz-scripts/hz-basic-summary-statistics.R")



OUT_DIR      <- "T:/EUCID"

OUT_RDS      <- file.path(OUT_DIR, "eu_cert_icd_cdm_ready.rds")

OUT_CSV      <- file.path(OUT_DIR, "eu_cert_icd_cdm_ready.csv")

OUT_QC_CSV   <- file.path(OUT_DIR, "eu_cert_icd_qc_summary.csv")

OUT_BSS_XLSX <- file.path(OUT_DIR, "eu_cert_icd_basic_summary.xlsx")



DB_NAME       <- "CERT"


EXCLUDE_CRT   <- TRUE



dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 3) Load CDM + dictionary (dictionary used ONLY for date conversion)
# =============================================================================
cdm <- fread(IN_CDM_CSV)
setnames(cdm, "Variable name", "Name")

dic <- as.data.table(openxlsx::read.xlsx(IN_DICT_XLSX))
stopifnot(all(c("original_name", "data_type") %in% names(dic)))

# =============================================================================
# 4) Load EU-CERT-ICD CSV
# =============================================================================
dt <- fread(IN_CSV, na.strings = "", stringsAsFactors = FALSE)
setDT(dt)

# =============================================================================
# 5) Convert date columns (dictionary-driven)
# =============================================================================
date_cols <- dic[data_type == "%Y-%m-%d", unique(original_name)]
date_cols <- date_cols[date_cols %in% names(dt)]

for (j in date_cols) {
  na_before <- sum(is.na(dt[[j]]))
  dt[, (j) := as.Date(as.character(get(j)), format = "%Y-%m-%d")]
  stopifnot(na_before == sum(is.na(dt[[j]])))
}

# =============================================================================
# 6) Sample exclusions (RAW EU-CERT names)
# =============================================================================
if ("pat_diag_type" %in% names(dt)) {
  dt <- dt[pat_diag_type == "ischemic"]
}

age_col <- NULL
if ("pat_implant_age" %in% names(dt)) age_col <- "pat_implant_age"
if (is.null(age_col) && "Age" %in% names(dt)) age_col <- "Age"
stopifnot(!is.null(age_col))
dt <- dt[get(age_col) >= 18]

if (EXCLUDE_CRT) {
  # EU-CERT raw often has icd_type; if not, don't crash—print hints
  if ("icd_type" %in% names(dt)) {
    dt <- dt[icd_type != "CRT-D"]
  } else {
    cat("\nWARNING: 'icd_type' not found. Candidates:\n")
    print(grep("icd|crt|device|type", names(dt), ignore.case = TRUE, value = TRUE))
  }
}

if ("ctr_name" %in% names(dt)) {
  dt <- dt[!ctr_name %in% c("Karolinska Institute Stockholm")]
}

# =============================================================================
# 7) QC: follow-up consistency (RAW names) — QC ONLY
# =============================================================================


# Here we verify that Time_death (analysis follow-up time for mortality)
# correctly represents the last follow-up time, by checking that
# length_fu_mortality matches (lastfu_date - icd_implant_date).
#
qc_need <- c("lastfu_date", "icd_implant_date", "length_fu_mortality")
cat("\n--- QC: FU consistency check (mortality) ---\n")
print(data.table(var = qc_need, present = qc_need %in% names(dt)))

if (all(qc_need %in% names(dt))) {
  
  dt[, fu_days_from_dates := as.numeric(lastfu_date - icd_implant_date)]
  dt[, fu_diff_mortality  := fu_days_from_dates - as.numeric(length_fu_mortality)]
  
  cat("\nSummary of fu_diff_mortality:\n")
  print(summary(dt$fu_diff_mortality))
  
  dt[, flag_fu_mortality_mismatch :=
       !is.na(fu_diff_mortality) & abs(fu_diff_mortality) > 7]
  
  cat("\nMismatch counts:\n")
  print(dt[, .N, by = flag_fu_mortality_mismatch])
  
} else {
  dt[, flag_fu_mortality_mismatch := NA]
  cat("\nQC skipped: missing one or more required columns.\n")
}


# =============================================================================

# 9) Endpoints (DAYS ONLY)

# =============================================================================



# --- Inappropriate shock (FIS) ---

dt[, inap_shock := tolower(trimws(as.character(inap_shock)))]



dt[, Status_FIS := NA_integer_]

dt[inap_shock == "no",  Status_FIS := 0L]

dt[inap_shock == "yes", Status_FIS := 1L]



dt[, Time_FIS := as.numeric(length_fu_inap_shock)]



dt[, flag_fis_missing_time := (Status_FIS == 1L & is.na(Time_FIS))]



if (all(c("inap_shock_date", "lastfu_date") %in% names(dt))) {
  
  dt[, flag_fis_after_fu :=
       
       Status_FIS == 1L &
       
       !is.na(inap_shock_date) & !is.na(lastfu_date) &
       
       inap_shock_date > lastfu_date]
  
} else {
  
  dt[, flag_fis_after_fu := NA]
  
}



# --- Death from any cause ---

dt[, death := tolower(trimws(as.character(death)))]



dt[, Status_death := NA_integer_]

dt[death == "no",  Status_death := 0L]

dt[death == "yes", Status_death := 1L]



if ("heart_transplant" %in% names(dt)) {
  
  dt[, heart_transplant := tolower(trimws(as.character(heart_transplant)))]
  
  dt[death == "yes" & heart_transplant == "yes", Status_death := 0L]
  
}



dt[, Time_death := as.numeric(length_fu_mortality)]



dt[, flag_death_missing_time := (Status_death == 1L & is.na(Time_death))]



if (all(c("death_date", "lastfu_date") %in% names(dt))) {
  
  dt[, flag_death_after_fu :=
       
       Status_death == 1L &
       
       !is.na(death_date) & !is.na(lastfu_date) &
       
       death_date > lastfu_date]
  
} else {
  
  dt[, flag_death_after_fu := NA]
  
}



# =============================================================================

# 10) Dataset metadata

# =============================================================================

id_col <- NULL

if ("pat_id" %in% names(dt)) id_col <- "pat_id"

if (is.null(id_col) && "ID" %in% names(dt)) id_col <- "ID"

stopifnot(!is.null(id_col))



dt[, ID := as.character(get(id_col))]

dt[, DB := DB_NAME]

dt[, ID_f := paste0(DB, "-", ID)]




# =============================================================================

# 11) Restrict to CDM + key endpoints

# =============================================================================

key_vars <- c(
  
  "ID", "ID_f", "DB", "Baseline_type",
  
  "Status_death", "Time_death",
  
  "Status_FIS", "Time_FIS","fu_days_from_dates"
  
)



keep_vars <- intersect(names(dt), unique(c(key_vars, cdm$Name)))

dt_cdm <- dt[, ..keep_vars]



setcolorder(
  
  dt_cdm,
  
  c(intersect(key_vars, names(dt_cdm)),
    
    setdiff(names(dt_cdm), key_vars))
  
)



# =============================================================================

# 12) Save processed dataset + QC

# =============================================================================

qc <- dt[, .(
  
  n = .N,
  
  
  
  n_Status_death_NA = sum(is.na(Status_death)),
  
  n_Time_death_NA   = sum(is.na(Time_death)),
  
  n_death_after_fu  = sum(flag_death_after_fu, na.rm = TRUE),
  
  
  
  n_Status_FIS_NA   = sum(is.na(Status_FIS)),
  
  n_Time_FIS_NA     = sum(is.na(Time_FIS)),
  
  n_fis_after_fu    = sum(flag_fis_after_fu, na.rm = TRUE),
  
  
  
  n_fu_mortality_mismatch = sum(flag_fu_mortality_mismatch, na.rm = TRUE)
  
)]

qc

fwrite(qc, OUT_QC_CSV)

saveRDS(dt_cdm, OUT_RDS)

fwrite(dt_cdm, OUT_CSV)



# =============================================================================

# 13) Basic summary statistics (HZ) — FINAL STEP

# =============================================================================

bss <- HZ.eda_describe_data(dt)$summary.non.sdc

HZ.eda_summary_write_excel(bss, "Summary", OUT_BSS_XLSX)

