
###############################################################################

# MAIN MERGE SCRIPT — PROFID ICD cohort + harmonised FIS (DAYS ONLY, RDS INPUT)

#

# What this script does:

# 1) Read ICD cohort (RDS) + 4 processed dataset RDS (EU-CERT, HELIOS, ISRAEL, PROSE)

# 2) Harmonise key column names across datasets (rename only; NO unit conversion)

# 3) Stack FIS exposure table (ID_f, DB, Status_death, Time_death_days,

#    t_followup_days, Status_FIS, Time_FIS_days)

# 4) Merge stacked table into ICD cohort by: ICD$ID == stacked$ID_f

# 5) QC: duplicates, match summary, counts before/after, FIS-after-follow-up flag

# 6) Save merged dataset (RDS + CSV) + QC outputs

###############################################################################



install.packages("tidyverse")
library(tidyverse)
library(dplyr)
library(survival)
library(data.table)
install.packages("mice")
library(mice)
install.packages("naniar")
library(naniar)
install.packages("openxlsx")
library(openxlsx)
library(readxl)
install.packages("gt")
library(gt)
install.packages("ggplot2")
library(ggplot2)
library(broom)


# =============================================================================

# 1) PATHS (EDIT)

# =============================================================================

IN_ICD_CSV   <- "//Charite.de/Centren/AG/f-dhzc-profid/Data Transfer to Charite/ICD.csv"



IN_CERT_RDS  <- "T:/EUCID/eu_cert_icd_cdm_ready.rds"

IN_HELS_RDS  <- "T:/HELIOS/helios/helios_cdm_ready.rds"

IN_ISRL_RDS  <- "T:/ISRAEL/israel_final/processed-isrl-common-data-model.rds"

IN_PRSE_RDS  <- "T:/PROSE/prose_cdm_ready.rds"



OUT_DIR      <- "T:/FINAL ICD COHORT"

OUT_RDS      <- file.path(OUT_DIR, "icd_merged1.rds")

OUT_CSV      <- file.path(OUT_DIR, "icd_merged1.csv")

OUT_QC_DIR   <- file.path(OUT_DIR, "qc")



dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

dir.create(OUT_QC_DIR, recursive = TRUE, showWarnings = FALSE)



TARGET_DB <- c("CERT", "HELS", "ISRL", "PRSE")



# =============================================================================

# 2) LOAD DATA

# =============================================================================

cat("\n=== LOAD DATA ===\n")



stopifnot(file.exists(IN_ICD_CSV))

d1 <- fread(IN_ICD_CSV)

stopifnot(all(c("ID","DB") %in% names(d1)))

d1[, ID := as.character(ID)]

d1[, DB := as.character(DB)]

cat("ICD cohort rows:", nrow(d1), "\n")



stopifnot(file.exists(IN_CERT_RDS), file.exists(IN_HELS_RDS),
          
          file.exists(IN_ISRL_RDS), file.exists(IN_PRSE_RDS))



d2 <- readRDS(IN_CERT_RDS); setDT(d2)

d3 <- readRDS(IN_HELS_RDS); setDT(d3)

d4 <- readRDS(IN_ISRL_RDS); setDT(d4)

d5 <- readRDS(IN_PRSE_RDS); setDT(d5)

cat("\nEU-CERT IDs:\n")
head(d2$ID_f, 20)

cat("\nHELIOS IDs:\n")
head(d3$ID_f, 20)

cat("\nISRAEL IDs:\n")
head(d4$ID_f, 20)

cat("\nPROSE IDs:\n")
head(d5$ID_f, 20)


normalize_id_proc <- function(x) {
  x <- as.character(x)
  x <- trimws(x)
  x <- gsub("-", "_", x)
  x
}

# Apply to CERT, HELS, PRSE only
d2[, ID_f := normalize_id_proc(ID_f)]
d3[, ID_f := normalize_id_proc(ID_f)]
d5[, ID_f := normalize_id_proc(ID_f)]
# =============================================================================

# 3) HARMONISE COLUMN NAMES (RENAME ONLY; DAYS ALREADY)

# =============================================================================

# Canonical names we want:

#   ID_f, DB,

#   Status_death, Time_death_days,

#   Status_FIS,   Time_FIS_days,

#   t_followup_days

#

# IMPORTANT: We do NOT change values or convert units here.



cat("\n=== HARMONISE NAMES ===\n")



# --- EU-CERT ---

# from your screenshot: Time_death, fu_days_from_dates, Status_FIS, Time_FIS

setnames(d2, "fu_days_from_dates", "t_followup_days", skip_absent = TRUE)

setnames(d2, "Time_death",         "Time_death_days", skip_absent = TRUE)

setnames(d2, "Time_FIS",           "Time_FIS_days",   skip_absent = TRUE)

if (!("DB" %in% names(d2))) d2[, DB := "CERT"]



# --- HELIOS ---

# already: Time_death_days, t_followup_days, Status_FIS, Time_FIS_days

if (!("DB" %in% names(d3))) d3[, DB := "HELS"]



# --- ISRAEL ---

# already: Time_death_days, t_followup_days, Status_FIS, Time_FIS_days

if (!("DB" %in% names(d4))) d4[, DB := "ISRL"]



# --- PROSE ---

# from your screenshot: Time_death, Time_inapp, t_followup_days

setnames(d5, "Time_death", "Time_death_days", skip_absent = TRUE)

setnames(d5, "Time_inapp", "Time_FIS_days",   skip_absent = TRUE)



# Harmonise PROSE ID prefix to match ICD master (only at start of string)

# Example: "PRSE_123" -> "PRSI_123"

d5[, ID_f := sub("^PRSE_", "PRSI_", ID_f)]



if (!("DB" %in% names(d5))) d5[, DB := "PRSE"]



# =============================================================================

# 4) SANITY CHECK REQUIRED FIELDS PER DATASET

# =============================================================================

req_cols <- c(

  "ID_f","DB",

  "Status_death","Time_death_days",

  "t_followup_days",

  "Status_FIS","Time_FIS_days"

)



check_required <- function(dt, label) {

  miss <- setdiff(req_cols, names(dt))

  if (length(miss) > 0) {

    cat("\nERROR: Missing columns in", label, ":\n")

    print(miss)

    stop("Missing required columns in dataset: ", label)

  }

}



check_required(d2, "EU-CERT")

check_required(d3, "HELIOS")

check_required(d4, "ISRAEL")

check_required(d5, "PROSE")

check_required <- function(dt, label) {
  
  miss <- setdiff(req_cols, names(dt))
  
  
  
  cat("\n------------------------------------------------------------\n")
  
  cat("Checking required columns for:", label, "\n")
  
  cat("------------------------------------------------------------\n")
  
  
  
  if (length(miss) > 0) {
    
    cat("MISSING (", length(miss), "): ", paste(miss, collapse = ", "), "\n", sep = "")
    
    cat("TIP: Here are columns that look related (death/fis/follow):\n")
    
    rel <- grep("death|fis|follow|fu|alive|time", names(dt), ignore.case = TRUE, value = TRUE)
    
    if (length(rel) == 0) rel <- "(none found)"
    
    print(rel)
    
    
    
    flush.console()  # important in RStudio
    
    stop("Missing required columns in dataset: ", label)
    
  } else {
    
    cat("OK: all required columns are present.\n")
    
  }
  
  
  
  invisible(TRUE)
  
}

# =============================================================================

# 5) BUILD STACKED TABLE (d6)

# =============================================================================

cat("\n=== BUILD STACKED TABLE (d6) ===\n")



d6 <- rbindlist(

  list(

    d2[, ..req_cols],

    d3[, ..req_cols],

    d4[, ..req_cols],

    d5[, ..req_cols]

  ),

  use.names = TRUE,

  fill = TRUE

)



cat("Stacked d6 rows:", nrow(d6), "\n")



# Duplicate ID_f check

dup_ids <- d6[, .N, by = ID_f][N > 1]

if (nrow(dup_ids) > 0) {

  cat("\nWARNING: duplicated ID_f found in d6 (showing first 20):\n")

  print(head(dup_ids, 20))

  fwrite(dup_ids, file.path(OUT_QC_DIR, "dup_id_f_in_d6.csv"))

}



# Completeness summary by DB

qc_completeness <- d6[, .(

  n = .N,

  n_Status_FIS_NA = sum(is.na(Status_FIS)),

  n_Time_FIS_NA   = sum(is.na(Time_FIS_days)),

  n_t_followup_NA = sum(is.na(t_followup_days)),

  n_Time_death_NA = sum(is.na(Time_death_days)),
  n_Status_death_NA = sum(is.na(Status_death))
), by = DB][order(DB)]



cat("\nCompleteness by DB:\n")

print(qc_completeness)
gt::gtsave(
  gt::gt(as.data.frame(qc_completeness)),
  file.path(OUT_QC_DIR, "completeness_by_db.html")
)
fwrite(qc_completeness, file.path(OUT_QC_DIR, "completeness_by_db.csv"))



# =============================================================================

# 6) COUNTS BEFORE MERGE (ICD cohort)

# =============================================================================

cat("\n=== COUNTS BEFORE MERGE (ICD cohort) ===\n")



counts_before <- d1[DB %in% TARGET_DB, .(n_patients = uniqueN(ID)), by = DB][order(DB)]

print(counts_before)

fwrite(counts_before, file.path(OUT_QC_DIR, "counts_before.csv"))
gt::gtsave(
  gt::gt(as.data.frame(counts_before)),
  file.path(OUT_QC_DIR, "counts_before.html")
)


# =============================================================================

# 7) MERGE INTO ICD COHORT

# =============================================================================

cat("\n=== MERGE INTO ICD COHORT ===\n")

# Drop DB from d6 so merge does NOT create DB.x / DB.y
d6_merge <- copy(d6)
d6_merge[, DB := NULL]

d1_pre <- copy(d1)



d1_merged <- merge(

  d1_pre,

  d6_merge,

  by.x = "ID",

  by.y = "ID_f",

  all.x = TRUE

)



# Now filter to target DBs
d1_final <- d1_merged[DB %in% TARGET_DB]

cat("Merged dimensions:", paste(dim(d1_merged), collapse = " x "), "\n")
colSums(is.na(d1_final[, c(
     "Status_FIS",
    "Time_FIS_days",
    "t_followup_days",
    "Time_death_days",
    "Status_death"
  ), with = FALSE]))

cat("ICD unique IDs:", uniqueN(d1_pre$ID), "\n")

cat("Stack unique IDs:", uniqueN(d6$ID_f), "\n")

cat("Overlap IDs:", length(intersect(d1_final$ID, d6$ID_f)), "\n")

cat("ICD only:", length(setdiff(d1_pre$ID, d6$ID_f)), "\n")

cat("Stack only:", length(setdiff(d6$ID_f, d1_pre$ID)), "\n")

head(d1_pre$ID, 20)

head(d6$ID_f, 20)

# Filter to target DBs (as per your main script)

d1_final <- d1_merged[DB %in% TARGET_DB]

cat("Matched (Status_death not NA):",
    d1_final[!is.na(Status_death), .N], "\n")

cat("Unmatched (Status_death is NA):",
    d1_final[is.na(Status_death), .N], "\n")
d1_final[is.na(Status_death), .N, by = DB][order(-N)]

# PRSE unmatched examples
d1_final[DB == "PRSE" & is.na(Status_death), head(ID, 50)]

# HELIOS unmatched examples
d1_final[DB == "HELS" & is.na(Status_death), head(ID, 20)]
cat("Filtered dimensions:", paste(dim(d1_final), collapse = " x "), "\n")

d1_analysis <- d1_final[!is.na(Status_death)]
cat("Final ICD analysis cohort:", nrow(d1_analysis), "\n")

# =============================================================================

# 8) QC: ID MATCH SUMMARY

# =============================================================================

cat("\n=== QC: ID MATCH SUMMARY ===\n")



n_overlap_ids  <- length(intersect(d1_pre$ID, d6$ID_f))

ICD_only_ids   <- setdiff(d1_pre$ID, d6$ID_f)

stack_only_ids <- setdiff(d6$ID_f, d1_pre$ID)



id_match_summary <- data.table(

  in_stacked = c(TRUE, TRUE, FALSE),

  in_icd     = c(TRUE, FALSE, TRUE),

  n          = c(n_overlap_ids, length(ICD_only_ids), length(stack_only_ids))

)



print(id_match_summary)

fwrite(id_match_summary, file.path(OUT_QC_DIR, "id_match_summary.csv"))
gt::gtsave(
  
  gt::gt(as.data.frame(id_match_summary)),
  
  file.path(OUT_QC_DIR, "id_match_summary.html")
  
)
d1_merged[!is.na(Status_death), .(
  n_rows = .N,
  n_ids  = uniqueN(ID)
)]

# =============================================================================

# 9) QC: COUNTS AFTER MERGE

# =============================================================================

cat("\n=== QC: COUNTS AFTER MERGE ===\n")



counts_after <- d1_analysis[, .(n_patients = uniqueN(ID)), by = DB][order(DB)]

print(counts_after)

fwrite(counts_after, file.path(OUT_QC_DIR, "counts_after.csv"))
gt::gtsave(
  gt::gt(as.data.frame(counts_after)),
  file.path(OUT_QC_DIR, "counts_after.html")
)


counts_compare <- merge(

  counts_before, counts_after,

  by = "DB", all = TRUE,

  suffixes = c("_before","_after")

)



cat("\nCounts before vs after:\n")

print(counts_compare)

fwrite(counts_compare, file.path(OUT_QC_DIR, "counts_before_after.csv"))
gt::gtsave(
  gt::gt(as.data.frame(counts_compare)),
  file.path(OUT_QC_DIR, "counts_compare.html")
)



# =============================================================================

# 10) QC: FIS AFTER FOLLOW-UP CHECK (DAYS ONLY)

# =============================================================================

cat("\n=== QC: FIS AFTER FOLLOW-UP (DAYS) ===\n")

cat("QC meaning: checks whether recorded FIS event time exceeds the recorded follow-up time.\n")

cat("Rule: if Status_FIS==1 then Time_FIS_days should be <= t_followup_days.\n\n")



d1_qc <- copy(d1_analysis)



d1_qc[, flag_fis_after_fu :=

        Status_FIS == 1L &

        !is.na(Time_FIS_days) &

        !is.na(t_followup_days) &

        Time_FIS_days > t_followup_days

]



cat("flag_fis_after_fu counts:\n")

print(d1_qc[, .N, by = flag_fis_after_fu])



cat("\nflag_fis_after_fu TRUE by DB:\n")

print(d1_qc[flag_fis_after_fu == TRUE, .N, by = DB][order(DB)])



bad <- d1_qc[flag_fis_after_fu == TRUE,

             .(ID, DB, Status_FIS, Time_FIS_days, t_followup_days, Time_death_days)]



if (nrow(bad) > 0) {

  fwrite(bad, file.path(OUT_QC_DIR, "fis_after_fu_rows.csv"))

  cat("\nSaved problematic rows to qc/fis_after_fu_rows.csv\n")

}



# =============================================================================

# 11) SAVE OUTPUTS

# =============================================================================

cat("\n=== SAVE OUTPUTS ===\n")



saveRDS(d1_analysis, OUT_RDS)

fwrite(d1_analysis, OUT_CSV)



cat("Saved merged RDS:", OUT_RDS, "\n")

cat("Saved merged CSV:", OUT_CSV, "\n")



