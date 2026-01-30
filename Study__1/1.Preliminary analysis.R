
###############################################################################

# PRELIMINARY ANALYSIS — MERGED ICD COHORT (CLEAN + REPRODUCIBLE)

#

# What this script does:

# 1) Load merged ICD cohort CSV

# 2) Cohort validity checks (LVEF + duplicate IDs)

# 3) Table 1 inadmissible value exclusions (set values to NA; no row drops)

# 4) ≤40-day baseline rule: set selected baseline measurements to NA + audit

# 5) Variable standardisation:

#    - Convert binary-like fields to ordered factor No < Yes (Table 1)

#    - Create bin_* 0/1 indicators for Cox modelling

#    - Keep NYHA as ordered factor I<II<III<IV

#    - Create bin_sex_male (Male=1, Female=0)

# 6) Save cleaned dataset (csv + rds)

#

# IMPORTANT:

# - Time variables are assumed to already be in DAYS. No unit conversion here.

###############################################################################



# --------------------------- #
## 1.PACKAGE SETUP            
# --------------------------- #
install.packages("tidyverse")
library(tidyverse)
library(dplyr)
library(survival)
library(data.table)
install.packages("openxlsx")
library(openxlsx)
library(readxl)
install.packages("gt")
library(gt)




# =============================================================================

# 1) Paths (EDIT)

# =============================================================================

IN_MERGED_CSV <- "T:/FINAL ICD COHORT/icd_merged1.csv"



OUT_DIR  <- "T:/FINAL ICD COHORT"

OUT_CSV  <- file.path(OUT_DIR, "standardised_data1.csv")

OUT_RDS  <- file.path(OUT_DIR, "standardised_data1.rds")



dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)



# =============================================================================

# 2) Load data

# =============================================================================

stopifnot(file.exists(IN_MERGED_CSV))



dc <- fread(IN_MERGED_CSV)

setDT(dc)
names(dc)
with(dc, table(Status_FIS, Status_death, useNA = "ifany"))
with(dc, table(Status_FIS, Status, useNA = "ifany"))
dc[, ID := as.character(ID)]

dc[, DB := as.character(DB)]



cat("\nLoaded merged ICD cohort\n")

cat("Rows:", nrow(dc), "\n")

cat("DB distribution:\n")

print(dc[, .N, by = DB][order(DB)])



# =============================================================================

# 3) Cohort validity checks (simple)

# =============================================================================

check_cohort_icd <- function(dt) {
  
  
  
  cat("\n========== Cohort Validity Checks ==========\n")
  
  
  
  if ("LVEF" %in% names(dt)) {
    
    lvef_qc <- dt %>%
      
      summarise(
        
        n_total   = n(),
        
        n_include = sum(LVEF <= 35, na.rm = TRUE),
        
        n_exclude = sum(LVEF >  35, na.rm = TRUE),
        
        n_missing = sum(is.na(LVEF))
        
      )
    
    print(lvef_qc)
    
  } else {
    
    cat("LVEF column not present\n")
    
  }
  
  
  
  dup <- dt %>%
    
    group_by(ID) %>%
    
    filter(n() > 1) %>%
    
    ungroup()
  
  
  
  if (nrow(dup) == 0) {
    
    cat("No duplicate IDs found\n")
    
  } else {
    
    cat("WARNING: duplicate IDs detected:", nrow(dup), "rows\n")
    
    print(dup %>% count(ID, sort = TRUE) %>% head(10))
    
  }
  
  
  
  cat("========== Check Complete ==========\n")
  
}



check_cohort_icd(dc)



# =============================================================================

# 4) Audit helper: count newly set to NA

# =============================================================================

check_exclusions <- function(before, after) {
  
  vars <- intersect(names(before), names(after))
  
  sapply(vars, function(v) sum(is.na(after[[v]]) & !is.na(before[[v]])))
  
}



# =============================================================================

# 5) Inadmissible value exclusions (Table 1 rules)

# =============================================================================

apply_inadmissible_exclusion <- function(df) {
  
  
  
  df <- copy(df)
  
  
  
  if ("BMI" %in% names(df))
    
    df[, BMI := fifelse(!is.na(BMI) & (BMI < 12 | BMI > 69), NA_real_, BMI)]
  
  
  
  if ("BUN" %in% names(df))
    
    df[, BUN := fifelse(!is.na(BUN) & BUN > 900, NA_real_, BUN)]
  
  
  
  if ("Haemoglobin" %in% names(df))
    
    df[, Haemoglobin :=
         
         fifelse(!is.na(Haemoglobin) & (Haemoglobin < 2 | Haemoglobin > 110),
                 
                 NA_real_, Haemoglobin)]
  
  
  
  if ("LDL" %in% names(df))
    
    df[, LDL := fifelse(!is.na(LDL) & LDL == 0, NA_real_, LDL)]
  
  
  
  if ("Sodium" %in% names(df))
    
    df[, Sodium := fifelse(!is.na(Sodium) & Sodium < 99, NA_real_, Sodium)]
  
  
  
  if ("Triglycerides" %in% names(df))
    
    df[, Triglycerides :=
         
         fifelse(!is.na(Triglycerides) & Triglycerides < 20,
                 
                 NA_real_, Triglycerides)]
  
  
  
  if ("TSH" %in% names(df))
    
    df[, TSH := fifelse(!is.na(TSH) & TSH == 0, NA_real_, TSH)]
  
  
  
  if ("HR" %in% names(df))
    
    df[, HR := fifelse(!is.na(HR) & (HR < 25 | HR > 140), NA_real_, HR)]
  
  
  
  if ("PR" %in% names(df))
    
    df[, PR := fifelse(!is.na(PR) & (PR <= 50 | PR > 1000), NA_real_, PR)]
  
  
  
  if ("QRS" %in% names(df))
    
    df[, QRS := fifelse(!is.na(QRS) & QRS < 50, NA_real_, QRS)]
  
  
  
  if ("QTc" %in% names(df))
    
    df[, QTc := fifelse(!is.na(QTc) & (QTc <= 250 | QTc > 790),
                        
                        NA_real_, QTc)]
  
  
  
  df
  
}



# =============================================================================

# 6) ≤40-day baseline rule (flag + apply + audit)

# =============================================================================

within_40d_flag <- function(df) {
  
  
  
  n <- nrow(df)
  
  flag <- rep(FALSE, n)
  
  
  
  # A) Baseline_type text labels
  
  if ("Baseline_type" %in% names(df)) {
    
    bt <- tolower(trimws(as.character(df$Baseline_type)))
    
    bt <- gsub("\\s+", " ", bt)
    
    pat <- "(before|pre|early|within\\s*40|<=\\s*40|<\\s*40)"
    
    flag <- flag | grepl(pat, bt)
    
  }
  
  
  
  # B) Years since MI
  
  yr_thr <- 40 / 365.25
  
  
  
  if ("Time_zero_Y" %in% names(df)) {
    
    y <- suppressWarnings(as.numeric(df$Time_zero_Y))
    
    flag <- flag | (!is.na(y) & y <= yr_thr)
    
  }
  
  
  
  if ("Time_zero_Ym" %in% names(df)) {
    
    ym <- suppressWarnings(as.numeric(df$Time_zero_Ym))
    
    flag <- flag | (!is.na(ym) & ym <= yr_thr)
    
  }
  
  
  
  # C) Days since MI
  
  if ("Time_index_MI_CHD" %in% names(df)) {
    
    d <- suppressWarnings(as.numeric(df$Time_index_MI_CHD))
    
    flag <- flag | (!is.na(d) & d <= 40)
    
  }
  
  
  
  flag
  
}



apply_40d_rule <- function(df) {
  
  
  
  vars_40d_all <- c("SBP", "DBP", "CRP", "Troponin_T", "NYHA",
                    
                    "AV_block", "AV_block_II_or_III")
  
  vars_40d <- intersect(vars_40d_all, names(df))
  
  
  
  flag40 <- within_40d_flag(df)
  
  before <- copy(df)
  
  
  
  if (length(vars_40d) > 0 && any(flag40)) {
    
    df[flag40, (vars_40d) := NA]
    
  }
  
  
  
  audit <- check_exclusions(before, df)
  
  
  
  list(
    
    data = df,
    
    within_40d = flag40,
    
    affected_vars = vars_40d,
    
    audit_new_NA = audit
    
  )
  
}



# =============================================================================

# 7) Variable standardisation (binary encoding)

# =============================================================================

to_no_yes <- function(x) {
  
  if (is.factor(x)) x <- as.character(x)
  
  
  
  if (is.numeric(x) || is.integer(x)) {
    
    ux <- unique(stats::na.omit(x))
    
    if (length(ux) > 0 && all(ux %in% c(0, 1))) {
      
      return(factor(ifelse(x == 1, "Yes", "No"),
                    
                    levels = c("No", "Yes"), ordered = TRUE))
      
    }
    
  }
  
  
  
  if (is.logical(x)) {
    
    return(factor(ifelse(x, "Yes", "No"),
                  
                  levels = c("No", "Yes"), ordered = TRUE))
    
  }
  
  
  
  if (is.character(x)) {
    
    z <- tolower(trimws(x))
    
    z[z %in% c("1", "y", "yes", "true", "t")]  <- "yes"
    
    z[z %in% c("0", "n", "no",  "false", "f")] <- "no"
    
    
    
    out <- ifelse(z == "yes", "Yes",
                  
                  ifelse(z == "no", "No", NA_character_))
    
    
    
    return(factor(out, levels = c("No", "Yes"), ordered = TRUE))
    
  }
  
  
  
  x
  
}



to_bin01 <- function(x) {
  
  if (is.factor(x)) x <- as.character(x)
  
  
  
  if (is.character(x)) {
    
    x <- ifelse(x %in% c("Yes", "yes"), 1L,
                
                ifelse(x %in% c("No", "no"), 0L, NA_integer_))
    
  } else if (is.logical(x)) {
    
    x <- ifelse(x, 1L, 0L)
    
  } else if (is.numeric(x) || is.integer(x)) {
    
    x <- ifelse(x %in% c(0, 1), as.integer(x), NA_integer_)
    
  } else {
    
    x <- NA_integer_
    
  }
  
  
  
  as.integer(x)
  
}



standardise_nyha_simple <- function(x) {
  
  factor(as.character(x), levels = c("I", "II", "III", "IV"), ordered = TRUE)
  
}



standardise_for_table1_and_cox <- function(dat, binary_vars, bin_prefix = "bin_") {
  
  
  
  dat <- as.data.frame(dat)
  
  
  
  if ("NYHA" %in% names(dat)) {
    
    dat$NYHA <- standardise_nyha_simple(dat$NYHA)
    
  }
  
  
  
  if ("Sex" %in% names(dat)) {
    
    dat$Sex <- factor(
      
      ifelse(dat$Sex %in% c("Female", "Male"), dat$Sex, NA),
      
      levels = c("Female", "Male")
      
    )
    
    dat[[paste0(bin_prefix, "sex_male")]] <- ifelse(dat$Sex == "Male", 1L,
                                                    
                                                    ifelse(dat$Sex == "Female", 0L, NA_integer_))
    
  }
  
  
  
  bin_present <- intersect(binary_vars, names(dat))
  
  
  
  for (v in bin_present) {
    
    dat[[v]] <- to_no_yes(dat[[v]])
    
  }
  
  
  
  for (v in bin_present) {
    
    out_name <- paste0(bin_prefix, tolower(v))
    
    dat[[out_name]] <- to_bin01(dat[[v]])
    
  }
  
  
  
  dat
  
}



# =============================================================================

# 8) Run the pipeline on dc

# =============================================================================

dc_before <- copy(dc)



dc_step1 <- apply_inadmissible_exclusion(dc)



res_40d <- apply_40d_rule(dc_step1)

dc_clean <- res_40d$data



cat("\n<=40 days flagged:", sum(res_40d$within_40d, na.rm = TRUE),
    
    "out of", nrow(dc_clean), "\n")

cat("Variables affected by 40d rule:", paste(res_40d$affected_vars, collapse = ", "), "\n\n")



audit_40d <- res_40d$audit_new_NA

audit_40d <- audit_40d[audit_40d > 0]

cat("Audit (new NA due to 40d rule):\n")

print(audit_40d)



inadmissible_audit <- check_exclusions(dc, dc_step1)

inadmissible_audit <- inadmissible_audit[inadmissible_audit > 0]

cat("\nAudit (new NA due to inadmissible rules):\n")

print(inadmissible_audit)



# =============================================================================

# 9) Standardise binary vars (Table 1 + Cox)

# =============================================================================

binary_vars <- c(
  
  "Diabetes","Hypertension","HF","COPD","Stroke_TIA","Smoking","FH_CAD","FH_SCD",
  
  "AV_block","AV_block_II_or_III","RBBB","LBBB","Anti_anginal","Cancer","AF_atrial_flutter",
  
  "NSVT","ACE_inhibitor","ARB","ACE_inhibitor_ARB","Beta_blockers","Diuretics",
  
  "Calcium_antagonists","Aldosterone_antagonist","Anti_arrhythmic_III","Anti_coagulant",
  
  "Anti_platelet","Lipid_lowering","Anti_diabetic","Anti_diabetic_oral",
  
  "Anti_diabetic_insulin","Digitalis_glycosides",
  
  "CABG","PCI","Thromolysis_acute","Revascularisation_acute","CABG_acute","PCI_acute",
  
  "MI_location_anterior"
  
)



dat2 <- standardise_for_table1_and_cox(dc_clean, binary_vars, bin_prefix = "bin_")



# Convert back to data.table for saving (optional, but convenient)

dat2 <- as.data.table(dat2)



# =============================================================================

# 10) Save outputs

# =============================================================================

fwrite(dat2, OUT_CSV)

saveRDS(dat2, OUT_RDS)



cat("\nSaved:\n")

cat("CSV:", OUT_CSV, "\n")

cat("RDS:", OUT_RDS, "\n")













