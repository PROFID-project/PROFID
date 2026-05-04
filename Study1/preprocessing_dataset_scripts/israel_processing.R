
###############################################################################
# ISRAELI-ICD (or similar) — PREPROCESSING TO PROFID CDM
# Outcome: all-cause mortality
# Exposure: first inappropriate shock (FIS)
# Time scale: DAYS since ICD implantation
###############################################################################

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


# ---------------------------------------------------------------------------
# Inputs (EDIT PATHS)
# ---------------------------------------------------------------------------
IN_RAW_CSV   <- "S:/AG/f-dhzc-profid/datasets/local/israeli-icd/data/original/ICDALL_20170630.csv"
IN_DICT_XLSX <- "S:/AG/f-dhzc-profid/datasets/local/israeli-icd/data/dictionary/israeli-icd-data-dictionary-raw-v3.xlsx"
IN_CDM_CSV   <- "S:/AG/f-dhzc-profid/datasets/local/israeli-icd/scripts/cdm/profid-common-data-model.csv"
IN_HZ_SCRIPT <- "S:/AG/f-dhzc-profid/datasets/local/israeli-icd/scripts/hz-scripts/hz-basic-summary-statistics.R"

OUT_DIR      <- "T:/ISRAEL"
OUT_RDS      <- file.path(OUT_DIR, "processed-isrl-common-data-model.rds")
OUT_CSV      <- file.path(OUT_DIR, "processed-ISRAEL-common-data-model.csv")
OUT_QC_CSV   <- file.path(OUT_DIR, "qc-israel-preprocessing.csv")
OUT_BSS_XLSX <- file.path(OUT_DIR, "basic-data-summaries1-cdm.xlsx")

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

EXCLUDE_CRT      <- FALSE
EXCLUDE_AGE_LT18 <- TRUE



# ---------------------------------------------------------------------------
# Read raw data
# ---------------------------------------------------------------------------
na_vals <- c("NA", "", "N", "Unknown", "UNKNOWN", "MISSING")
dt <- fread(IN_RAW_CSV, na.strings = na_vals, data.table = TRUE, check.names = FALSE)
dt <- as.data.table(dt)
dt_raw <- copy(dt)

# ---------------------------------------------------------------------------
# Rename using data dictionary
# ---------------------------------------------------------------------------
dd <- as.data.table(read_excel(IN_DICT_XLSX))
stopifnot(all(c("original_name", "new_name") %in% names(dd)))
stopifnot(!any(is.na(dd$new_name) | dd$new_name == ""))

stopifnot(identical(names(dt), dd$original_name))
setnames(dt, old = names(dt), new = dd$new_name)

# ---------------------------------------------------------------------------
# Optional date parsing (only if dictionary contains a `data_type` column)
# ---------------------------------------------------------------------------
if ("data_type" %in% names(dd)) {
  
  parse_date_cols <- function(cols, fmt) {
    cols <- intersect(cols, names(dt))
    if (length(cols) == 0L) return(invisible(NULL))
    for (v in cols) dt[, (v) := as.IDate(as.Date(get(v), format = fmt))]
  }
  
  parse_date_cols(dd[data_type == "date_b_d_Y", new_name], "%b %d, %Y")
  parse_date_cols(dd[data_type == "date_dbY",   new_name], "%d/%b/%Y")
}

# ---------------------------------------------------------------------------
# Exclusions (optional)
# ---------------------------------------------------------------------------
if (EXCLUDE_AGE_LT18) {
  stopifnot("Age" %in% names(dt))
  dt <- dt[Age >= 18]
}

if (EXCLUDE_CRT) {
  if ("Device_type" %in% names(dt)) {
    dt <- dt[Device_type %chin% c("SINGLE CHAMBER", "DUAL CHAMBER") | is.na(Device_type)]
  } else if ("ICD_type" %in% names(dt)) {
    dt <- dt[ICD_type %chin% c("SINGLE - Single", "DUAL - Dual")]
  } else {
    stop("EXCLUDE_CRT=TRUE but neither Device_type nor ICD_type exists.")
  }
}

# ---------------------------------------------------------------------------
# Required columns for outcome/exposure logic
# ---------------------------------------------------------------------------
req <- c("Status_last", "Alive_last_FU_days", "Alive_total_days", "Inapp_shock_1st", "Inapp_shock_days")
miss <- setdiff(req, names(dt))
if (length(miss) > 0L) stop("Missing required columns: ", paste(miss, collapse = ", "))

dt[, Alive_last_FU_days := as.numeric(Alive_last_FU_days)]
dt[, Alive_total_days   := as.numeric(Alive_total_days)]
dt[, Inapp_shock_days   := as.numeric(Inapp_shock_days)]

# ---------------------------------------------------------------------------
# Fix: deaths after last FU (align FU to death time)
# ---------------------------------------------------------------------------
dt[, flag_death_after_fu := FALSE]
dt[Status_last == "DIED" &
     !is.na(Alive_total_days) & !is.na(Alive_last_FU_days) &
     Alive_last_FU_days != Alive_total_days,
   flag_death_after_fu := TRUE]

dt[Status_last == "DIED", Alive_last_FU_days := Alive_total_days]

stopifnot(dt[Status_last == "DIED" &
               !is.na(Alive_total_days) & !is.na(Alive_last_FU_days),
             all(Alive_last_FU_days == Alive_total_days)])

# ---------------------------------------------------------------------------
# Outcome: all-cause mortality (Status + time)
# --------------------------------------------------------------------------

## Death from any cause (event = DIED, else censored)
dt[, .N, by = Status_last]

dt[, Status_death := 0L]                 # default: censored
dt[Status_last == "DIED", Status_death := 1L]

# QC: confirm binary and see mapping
dt[, .N, by = .(Status_last, Status_death)]
stopifnot(all(dt$Status_death %in% c(0L, 1L)))


dt[, Time_death_days := Alive_last_FU_days]   # after correction: deaths have FU == death time
dt[, t_followup_days := Alive_last_FU_days]   # harmonised follow-up time (days)


# ------------------------------------------------------------------
# 7) EXPOSURE: FIRST INAPPROPRIATE SHOCK (Status + time)
#   Uses day-offset variable(s) already in DAYS. If naming differs, edit below.
# ------------------------------------------------------------------
stopifnot("Inapp_shock_1st" %in% names(dt))
stopifnot("Inapp_shock_days" %in% names(dt))
stopifnot("Alive_last_FU_days" %in% names(dt))

dt[, Inapp_shock_days := as.numeric(Inapp_shock_days)]

dt[, Status_FIS := NA_integer_]
dt[Inapp_shock_1st == "YES", Status_FIS := 1L]
dt[Inapp_shock_1st == "NO",  Status_FIS := 0L]

# Time to FIS in DAYS: event time if event, else censor at last FU
dt[, Time_FIS_days := NA_real_]
dt[Status_FIS == 1L, Time_FIS_days := Inapp_shock_days]
dt[Status_FIS == 0L, Time_FIS_days := Alive_last_FU_days]
dt[Time_FIS_days < 0, Time_FIS_days := NA_real_]

# Consistency: event must not occur after follow-up
dt[, flag_fis_after_fu := !is.na(Inapp_shock_days) & !is.na(t_followup_days) & Inapp_shock_days > t_followup_days]

dt[flag_fis_after_fu == TRUE,
   `:=`(Status_FIS = 0L,
        Time_FIS_days = t_followup_days)]

# ---------------------------------------------------------------------------
# CDM load (robust Name column)
# ---------------------------------------------------------------------------
cdm <- fread(IN_CDM_CSV)
if (!("Name" %in% names(cdm)) && ("Variable name" %in% names(cdm))) {
  setnames(cdm, "Variable name", "Name")
}
stopifnot("Name" %in% names(cdm))

# ---------------------------------------------------------------------------
# Build CDM dataset
# ---------------------------------------------------------------------------
if (!("DB" %in% names(dt))) dt[, DB := "ISRL"]
if (!("ID" %in% names(dt))) stop("ID column not found. Ensure dictionary maps patient id to 'ID'.")

dt[, ID_f := paste0(DB, "_", ID)]

key_vars <- c(
  "ID", "ID_f", "DB",
  "Status_death", "Time_death_days", "t_followup_days",
  "Status_FIS", "Time_FIS_days"
)

keep_vars <- unique(c(key_vars, cdm$Name))
keep_vars <- intersect(keep_vars, names(dt))

dt_cdm <- dt[, ..keep_vars]

setcolorder(
  dt_cdm,
  c(intersect(key_vars, names(dt_cdm)),
    setdiff(names(dt_cdm), intersect(key_vars, names(dt_cdm))))
)

# ---------------------------------------------------------------------------
# Save outputs + QC
# ---------------------------------------------------------------------------
qc <- rbindlist(list(
  data.table(metric = "n_raw_after_rename", value = nrow(dt_raw)),
  data.table(metric = "n_after_exclusions", value = nrow(dt)),
  data.table(metric = "n_cdm", value = nrow(dt_cdm)),
  data.table(metric = "n_deaths", value = sum(dt$Status_death == 1L, na.rm = TRUE)),
  data.table(metric = "n_fis", value = sum(dt$Status_FIS == 1L, na.rm = TRUE)),
  data.table(metric = "n_death_after_fu_fixed", value = sum(dt$flag_death_after_fu, na.rm = TRUE)),
  data.table(metric = "n_fis_after_fu", value = sum(dt$flag_fis_after_fu, na.rm = TRUE))
))

fwrite(qc, OUT_QC_CSV)
saveRDS(dt_cdm, OUT_RDS)
fwrite(dt_cdm, OUT_CSV)
print(qc)
# ---------------------------------------------------------------------------
# HZ Basic summary statistics (LAST)
# ---------------------------------------------------------------------------
if (file.exists(IN_HZ_SCRIPT)) {
  source(IN_HZ_SCRIPT)
  bss <- HZ.eda_describe_data(dt_cdm)$summary.non.sdc
  HZ.eda_summary_write_excel(bss, "Summary", OUT_BSS_XLSX)
}
