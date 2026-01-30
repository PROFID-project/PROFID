###############################################################################

# HELIOS preprocessing (PROFID CDM-ready) — EXPOSURE = inappropriate therapy

# PRIMARY OUTCOME = all-cause mortality


# Steps:

# 1) Load libraries

# 2) Load HELIOS Excel sheets and merge to patient-level (PAT_INDEX)

# 3) Load dictionary + CDM and rename variables once

# 4) Sample exclusions (Age < 18; optional CRT exclusion)

# 5) Define follow-up (t_followup_days) + death endpoint (Status_death, Time_death_days)

# 6) Define exposure variables (Status_FIS/Status_ATP + exposure times in days)

# 7) Add metadata (DB, ID_f, Baseline_type)

# 8) Restrict to CDM + key analysis variables and save

# 9) HZ basic summary statistics (FINAL STEP)

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
library(readxl)
library(stringr)



# =============================================================================

# 2) Paths / configuration (EDIT)

# =============================================================================

IN_XLSX     <- "//Charite.de/Centren/AG/f-dhzc-profid/datasets/local/helios-rdb/data/original/Final_delivery.2021-05-20._Ali EDxlsx.xlsx"

IN_DICT_CSV <- "//Charite.de/Centren/AG/f-dhzc-profid/datasets/local/helios-rdb/data/dictionary/helios-data-dictionary-raw.csv"

IN_CDM_CSV  <- "S:/AG/f-dhzc-profid/datasets/local/helios-rdb/scripts/cdm/profid-common-data-model.csv"



source("S:/AG/f-dhzc-profid/datasets/local/helios-rdb/scripts/hz-scripts/hz-basic-summary-statistics.R")



OUT_DIR      <- "T:/HELIOS"

OUT_RDS      <- file.path(OUT_DIR, "helios_cdm_ready.rds")

OUT_CSV      <- file.path(OUT_DIR, "helios_cdm_ready.csv")

OUT_QC_CSV   <- file.path(OUT_DIR, "helios_qc_summary.csv")

OUT_BSS_XLSX <- file.path(OUT_DIR, "helios_basic_summary.xlsx")



DB_NAME       <- "HELS"

EXCLUDE_CRT   <- FALSE



dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)



# =============================================================================

# 3) Load CDM + dictionary

# =============================================================================

cdm <- fread(IN_CDM_CSV)

stopifnot("Variable name" %in% names(cdm))

setnames(cdm, "Variable name", "Name")

stopifnot("Name" %in% names(cdm))



dd <- fread(IN_DICT_CSV)



# Expected dictionary columns: original_name / new_name

if (!all(c("original_name", "new_name") %in% names(dd))) {
  
  if (all(c("raw_name", "cdm_name") %in% names(dd))) {
    
    setnames(dd, c("raw_name", "cdm_name"), c("original_name", "new_name"))
    
  }
  
}

stopifnot(all(c("original_name", "new_name") %in% names(dd)))



# =============================================================================

# 4) Read and merge HELIOS Excel sheets (patient-level; PAT_INDEX)

# =============================================================================

sheets <- openxlsx::getSheetNames(IN_XLSX)



read_sheet_dt <- function(sheet) {
  
  x <- openxlsx::read.xlsx(IN_XLSX, sheet = sheet, detectDates = TRUE)
  
  as.data.table(x)
  
}



req_sheets <- c(
  
  "target_pop",
  
  "Baseline Characteristics",
  
  "Past Medical History (ICD10)",
  
  "Past Medical History (OPS)",
  
  "Medication (ATC)",
  
  "Lab Data",
  
  "ECG",
  
  "ICD queries",
  
  "Outcome",
  
  "CMR-Scar and GZ",
  
  "Imaging"
  
)

stopifnot(length(setdiff(req_sheets, sheets)) == 0)



mi <- read_sheet_dt("target_pop")

stopifnot("PAT_INDEX" %in% names(mi))

stopifnot(!anyDuplicated(mi$PAT_INDEX))



bs <- read_sheet_dt("Baseline Characteristics")

stopifnot("PAT_INDEX" %in% names(bs))

stopifnot(!anyDuplicated(bs$PAT_INDEX))



dt <- merge(mi, bs, by = "PAT_INDEX", all = TRUE)



merge_pat <- function(dt_main, sheet_name) {
  
  x <- read_sheet_dt(sheet_name)
  
  stopifnot("PAT_INDEX" %in% names(x))
  
  stopifnot(!anyDuplicated(x$PAT_INDEX))
  
  merge(dt_main, x, by = "PAT_INDEX", all = TRUE)
  
}



dt <- merge_pat(dt, "Past Medical History (ICD10)")

dt <- merge_pat(dt, "Past Medical History (OPS)")

dt <- merge_pat(dt, "Medication (ATC)")

dt <- merge_pat(dt, "Lab Data")

dt <- merge_pat(dt, "ECG")

dt <- merge_pat(dt, "ICD queries")

dt <- merge_pat(dt, "Outcome")

dt <- merge_pat(dt, "CMR-Scar and GZ")



# Imaging can be multi-row per patient; keep first row per PAT_INDEX (simple)

img <- read_sheet_dt("Imaging")

stopifnot("PAT_INDEX" %in% names(img))

img_first <- img[order(PAT_INDEX)][, lapply(.SD, function(v) v[1]), by = PAT_INDEX]

dt <- merge(dt, img_first, by = "PAT_INDEX", all = TRUE)



# =============================================================================

# 5) Rename variables once (dictionary-driven)

# =============================================================================

dd_use <- dd[original_name %in% names(dt) & !is.na(new_name) & nzchar(new_name)]

if (nrow(dd_use) > 0) {
  
  setnames(dt, dd_use$original_name, dd_use$new_name)
  
}



# =============================================================================

# 6) Sample exclusions

# =============================================================================

age_col <- NULL

if ("Age" %in% names(dt)) age_col <- "Age"

if (is.null(age_col) && "Alter.ICD" %in% names(dt)) age_col <- "Alter.ICD"

stopifnot(!is.null(age_col))



dt <- dt[get(age_col) >= 18]



if (EXCLUDE_CRT) {
  
  stopifnot("ICD_type" %in% names(dt))
  
  dt <- dt[ICD_type %chin% c("SINGLE - Single", "DUAL - Dual")]
  
}



# =============================================================================

# 7) Follow-up and PRIMARY OUTCOME: all-cause mortality (DAYS)

# =============================================================================

# Required:

# - Status_death_cat: "Yes"/"No" (or similar)

# - DAYS2LastFU.ICD, DAYS2DEATH.ICD: time since ICD implant (days)



death_cat_col <- NULL

if ("Status_death_cat" %in% names(dt)) death_cat_col <- "Status_death_cat"

if (is.null(death_cat_col) && "Status_death" %in% names(dt)) death_cat_col <- "Status_death"

stopifnot(!is.null(death_cat_col))



stopifnot(all(c("DAYS2LastFU.ICD", "DAYS2DEATH.ICD") %in% names(dt)))



dt[, DAYS2LastFU.ICD := as.numeric(DAYS2LastFU.ICD)]

dt[, DAYS2DEATH.ICD  := as.numeric(DAYS2DEATH.ICD)]



dt[, tmp_death := tolower(trimws(as.character(get(death_cat_col))))]



dt[, Status_death := NA_integer_]

dt[tmp_death == "no",  Status_death := 0L]

dt[tmp_death == "yes", Status_death := 1L]



dt[, Time_death_days := NA_real_]

dt[Status_death == 0L, Time_death_days := DAYS2LastFU.ICD]

dt[Status_death == 1L, Time_death_days := DAYS2DEATH.ICD]



# Harmonised follow-up in days since ICD implantation (primary time scale)

dt[, t_followup_days := Time_death_days]



# Basic validity

dt[t_followup_days < 0, t_followup_days := NA_real_]

dt[Time_death_days < 0, Time_death_days := NA_real_]



# QC: death after last follow-up (should not happen often)

dt[, flag_death_after_fu :=
     
     Status_death == 1L &
     
     !is.na(DAYS2DEATH.ICD) & !is.na(DAYS2LastFU.ICD) &
     
     DAYS2DEATH.ICD > DAYS2LastFU.ICD]

dt$flag_death_after_fu

dt[, tmp_death := NULL]



# =============================================================================

# 8) EXPOSURE variables (time-from-implant in DAYS)

# =============================================================================

# Exposure 1: Inappropriate shock

stopifnot(all(c("inappropriate_shock", "DAYS2_inappropriate_shock.ICD") %in% names(dt)))



dt[, inappropriate_shock := tolower(trimws(as.character(inappropriate_shock)))]

dt[, DAYS2_inappropriate_shock.ICD := as.numeric(DAYS2_inappropriate_shock.ICD)]



dt[, Status_FIS := NA_integer_]

dt[inappropriate_shock == "no",  Status_FIS := 0L]

dt[inappropriate_shock == "yes", Status_FIS := 1L]



dt[, Time_FIS_days := NA_real_]
dt[Status_FIS == 0L, Time_FIS := DAYS2LastQuery.ICD]
dt[Status_FIS == 1L, Time_FIS_days := DAYS2_inappropriate_shock.ICD]





# QC: exposure time should not exceed follow-up time for mortality

dt[, flag_fis_after_fu :=
     
     Status_FIS == 1L & !is.na(Time_FIS_days) & !is.na(t_followup_days) &
     
     Time_FIS_days > t_followup_days]


dt$flag_fis_after_fu




# =============================================================================

# 9) Metadata

# =============================================================================

dt[, ID := as.character(PAT_INDEX)]

dt[, DB := DB_NAME]

dt[, ID_f := paste0(DB, "-", ID)]



# =============================================================================

# 10) Restrict to CDM + key analysis variables

# =============================================================================

key_vars <- c(
  
  "ID", "ID_f", "DB",
  
  "Status_death", "Time_death_days", "t_followup_days",
  
  "Status_FIS", "Time_FIS_days"
  
  
)



keep_vars <- intersect(names(dt), unique(c(key_vars, cdm$Name)))

dt_cdm <- dt[, ..keep_vars]



setcolorder(
  
  dt_cdm,
  
  c(intersect(key_vars, names(dt_cdm)),
    
    setdiff(names(dt_cdm), key_vars))
  
)



# =============================================================================

# 11) Save processed dataset + QC

# =============================================================================

qc <- dt[, .(
  
  n = .N,
  
  
  
  n_Status_death_NA = sum(is.na(Status_death)),
  
  n_Time_death_NA   = sum(is.na(Time_death_days)),
  
  n_death_after_fu  = sum(flag_death_after_fu, na.rm = TRUE),
  
  n_Status_FIS_NA   = sum(is.na(Status_FIS)),
  
  n_Time_FIS_NA     = sum(is.na(Time_FIS_days)),
  
  n_fis_after_fu    = sum(flag_fis_after_fu, na.rm = TRUE),
  
  n_fu_missing      = sum(is.na(t_followup_days))
  
)]



fwrite(qc, OUT_QC_CSV)

saveRDS(dt_cdm, OUT_RDS)

fwrite(dt_cdm, OUT_CSV)



# =============================================================================

# 12) HZ basic summary statistics (FINAL STEP)

# =============================================================================

bss <- HZ.eda_describe_data(dt)$summary.non.sdc

HZ.eda_summary_write_excel(bss, "Summary", OUT_BSS_XLSX)
