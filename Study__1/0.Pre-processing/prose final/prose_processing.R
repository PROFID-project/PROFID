###############################################################################
# PROSE preprocessing (PROFID CDM-ready)
#
# Workflow:
# 1. Load libraries
# 2. Load raw data + dictionary + CDM
# 3. Rename variables (single source of truth)
# 4. Exclude co-enrolled patients (PROSE-only
# 5. Construct endpoints:
#    - Death from any cause
#    - Inappropriate ICD therapy (FIS)
# 6. Add dataset metadata (DB, ID_f, Baseline_type)
# 7. Restrict to CDM + key endpoint variables
# 8. Save processed dataset
# 9. Run HZ basic summary statistics (FINAL STEP)
#
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
library(haven)


# =============================================================================
# 2) Paths / configuration
# =============================================================================

# Inputs
IN_PROSE_CSV   <- "//Charite.de/Centren/AG/f-dhzc-profid/datasets/local/prose-icd/data/original/FinaltoPROFID_PROSEonlysent_hopkins_prose_study.csv"
IN_JOIN_ID_CSV <- "//Charite.de/Centren/AG/f-dhzc-profid/datasets/local/prose-icd/data/original/FinaltoPROFID_PROSEonlysent_coenrolled.csv"
IN_DICT_XLSX   <- "//Charite.de/Centren/AG/f-dhzc-profid/datasets/local/prose-icd/data/dictionary/prose-only-data-dictionary-raw-v1.xlsx"
IN_CDM_CSV     <- "S:/AG/f-dhzc-profid/datasets/local/prose-lvscd/scripts/cdm/profid-common-data-model.csv"

source("//Charite.de/Centren/AG/f-dhzc-profid/datasets/local/prose-icd/scripts/hz-scripts/hz-basic-summary-statistics.R")

# Outputs
OUT_DIR        <- "T:/PROSE"
OUT_RDS        <- file.path(OUT_DIR, "prose_cdm_ready.rds")
OUT_CSV        <- file.path(OUT_DIR, "prose_cdm_ready.csv")
OUT_QC_CSV     <- file.path(OUT_DIR, "prose_qc_summary.csv")
OUT_BSS_XLSX   <- file.path(OUT_DIR, "prose_basic_summary.xlsx")

# Dataset tags
DB_NAME       <- "PRSE"
#BASELINE_TYPE <- "MI40d"

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
msg <- function(...) cat(sprintf(...), "\n")

# =============================================================================

# 3) Load data

# =============================================================================

dt <- fread(IN_PROSE_CSV, na.strings = c("", "NA", "N/A"))

stopifnot("ID" %in% names(dt))



dd_prose <- as.data.table(read_excel(IN_DICT_XLSX))

stopifnot(all(c("original_name", "new_name") %in% names(dd_prose)))



id_join <- fread(IN_JOIN_ID_CSV)

if (all(c("V1", "V2") %in% names(id_join))) {
  
  setnames(id_join, c("V1", "V2"), c("ID_PROSE", "ID_LVSCD"))
  
}

stopifnot("ID_PROSE" %in% names(id_join))



cdm <- fread(IN_CDM_CSV)

stopifnot("Variable name" %in% names(cdm))

setnames(cdm, "Variable name", "Name")

stopifnot("Name" %in% names(cdm))



# =============================================================================

# 4) Rename variables (apply once)

# =============================================================================

dd_use <- dd_prose[
  
  original_name %in% names(dt) &
    
    !is.na(new_name) & nzchar(new_name)
  
]



if (nrow(dd_use) > 0) {
  
  setnames(dt, dd_use$original_name, dd_use$new_name)
  
}



# =============================================================================

# 5) Sample exclusions

# =============================================================================

stopifnot("Age" %in% names(dt))

dt <- dt[Age >= 18]

# Optional exclusion
EXCLUDE_CRT <- TRUE

if (EXCLUDE_CRT) {
  
  stopifnot("ICD_type" %in% names(dt))
  
  dt <- dt[ICD_type %chin% c("SINGLE - Single", "DUAL - Dual")]
  
}



# =============================================================================

# 6) PROSE-only cohort (exclude co-enrolled)

# =============================================================================

dt <- dt[!(ID %in% id_join$ID_PROSE)]

# =============================================================================
# 6B) PROSE follow-up time (no calendar FU date) — DAYS ONLY
# =============================================================================
# Define follow-up as the latest observed time across available event-time vars.
# This becomes the censoring time for non-events.

fu_vars <- c("Death_days", "days_to_app_shock", "t_inappshock")
fu_vars <- fu_vars[fu_vars %in% names(dt)]   # keep only those that exist

# If you have different names in your raw file, add them here.
stopifnot(length(fu_vars) > 0)

dt[, (fu_vars) := lapply(.SD, as.numeric), .SDcols = fu_vars]

dt[, t_followup_days := do.call(pmax, c(.SD, na.rm = TRUE)), .SDcols = fu_vars]

# If all FU vars are NA, pmax(..., na.rm=TRUE) yields -Inf -> set to NA
dt[is.infinite(t_followup_days), t_followup_days := NA_real_]

# Disallow negative follow-up times
dt[t_followup_days < 0, t_followup_days := NA_real_]


# =============================================================================

# 7) Endpoint: Death from any cause

# =============================================================================

stopifnot(all(c("Death_status", "Death_days") %in% names(dt)))



dt[, Status_death := NA_integer_]

dt[Death_status == "Yes", Status_death := 1L]

dt[Death_status == "No",  Status_death := 0L]



dt[, Time_death := as.numeric(Death_days)]



# =============================================================================

# 8) Exposure variable: Inappropriate ICD therapy (FIS)

# =============================================================================

stopifnot(all(c("inappshock", "t_inappshock", "Death_days") %in% names(dt)))



dt[, inappshock := trimws(as.character(inappshock))]

dt[inappshock %chin% c("", "NA", "N/A"), inappshock := NA_character_]



# Primary assumption:

# Missing shock info + follow-up available => no inappropriate shock

dt[is.na(inappshock) & !is.na(Death_days), inappshock := "No"]



dt[, Status_FIS := NA_integer_]

dt[inappshock == "Yes", Status_FIS := 1L]

dt[inappshock == "No",  Status_FIS := 0L]



dt[, Time_inapp := NA_real_]

dt[Status_FIS == 1L, Time_inapp := as.numeric(t_inappshock)]

dt[Status_FIS == 0L, Time_inapp := as.numeric(Death_days)]



# QC flags

dt[, flag_fis_missing_time :=
     
     Status_FIS == 1L & is.na(Time_inapp)]

dt[, flag_fis_time_after_death :=
     
     !is.na(Time_inapp) & !is.na(Death_days) & Time_inapp > Death_days]

dt$flag_fis_time_after_death

# =============================================================================

# 9) Add dataset metadata

# =============================================================================

dt[, DB := DB_NAME]

dt[, ID_f := paste0(DB, "-", ID)]




# =============================================================================

# 10) Restrict to CDM + key endpoint variables

# =============================================================================

key_vars <- c(
  
  "ID", "ID_f", "DB",
  
  "Status_death", "Time_death",
  
  "Status_FIS", "Time_inapp","t_followup_days"
  
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

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)



qc <- dt[, .(
  
  n = .N,
  
  n_Status_death_NA = sum(is.na(Status_death)),
  
  n_Time_death_NA   = sum(is.na(Time_death)),
  
  n_Status_FIS_NA   = sum(is.na(Status_FIS)),
  
  n_Time_inapp_NA   = sum(is.na(Time_inapp)),
  
  n_fis_missing_time = sum(flag_fis_missing_time, na.rm = TRUE),
  
  n_fis_time_after_death = sum(flag_fis_time_after_death, na.rm = TRUE)
  
)]



fwrite(qc, OUT_QC_CSV)

saveRDS(dt_cdm, OUT_RDS)

fwrite(dt_cdm, OUT_CSV)



# =============================================================================

# 12) Basic summary statistics (HZ) — FINAL STEP

# =============================================================================

bss <- HZ.eda_describe_data(dt)$summary.non.sdc



HZ.eda_summary_write_excel(bss, "Summary", OUT_BSS_XLSX)