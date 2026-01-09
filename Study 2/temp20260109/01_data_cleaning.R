# ==========================================================
# Project:   Seasonal Patterns in Sudden Cardiac Death (SCD)
#             Incidence Post-Myocardial Infarction (MI)
# Author:    Dr. Alexia Sampri
# Institution: UMBIZO
# Language:  R (tidyverse workflow)
# Version:   0.1.0
# Date:      2025-11-14
# Description:
#   Initial setup of the R project environment and folder structure.
#   This script creates directories, sets key parameters, and logs the setup.
# ==========================================================

# ---------------------
# Load required packages
# ---------------------

install.packages("dplyr")
install.packages("lubridate")
install.packages("timeDate")
install.packages("stringr")
install.packages("survival")
install.packages("tidyverse")
install.packages("naniar")
install.packages("VIM")
install.packages('epitools')
## -----------------------------
## Temporal prep (no custom funcs)
## -----------------------------
## =========================================================
## Packages
## =========================================================
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(timeDate)   # for Easter(); safe + lightweight
  # library(survival)   # for Cox model
  library(tidyverse)
})

# ---------------------
# Define project metadata
# ---------------------
project_name    <- "Seasonal SCD Post-MI"
project_version <- "0.1.0"
project_date    <- Sys.Date()
project_author  <- "Dr. Alexia Sampri"

message("Starting project: ", project_name)
message("Version: ", project_version)
message("Date: ", project_date)

# ---------------------
# Directory structure
# ---------------------
dirs <- c(
  "data/raw",
  "data/processed",
  "scripts",
  "output/tables",
  "output/figures",
  "output/models"
)

for (d in dirs) dir.create(d, showWarnings = FALSE, recursive = TRUE)

# ---------------------
# Save metadata
# ---------------------
project_info <- tibble(
  project_name = project_name,
  version = project_version,
  author = project_author,
  created_on = project_date
)

write_csv(project_info, "project_metadata.csv")

message("✅ Project setup complete.")

setwd("T:/PROFID/data/raw")
## ==========================================================
## Load 4 datasets, check variables, find common vs extra,
## then join into one combined dataset
## ==========================================================

## 1. Load datasets (adjust paths to your files)
df1 <- read.csv("ICD_all.csv")
df2 <- read.csv("ICD.csv")
df3 <- read.csv("NonICD_preserved.csv")
df4 <- read.csv("NonICD_reduced.csv")

## 2. Add dataset identifiers
df1$dataset <- "ICD_all"
df2$dataset <- "ICD"
df3$dataset <- "NonICD_preserved"
df4$dataset <- "NonICD_reduced"

## 3. Check how many columns each dataset has
sapply(list(ICD_all=df1, ICD=df2, NonICD_preserved=df3, NonICD_reduced=df4), ncol)

## 4. Work out common vs all variables
common_vars <- Reduce(intersect, list(names(df1), names(df2), names(df3), names(df4)))
all_vars    <- Reduce(union,     list(names(df1), names(df2), names(df3), names(df4)))

cat("Number of common variables:", length(common_vars), "\n")
cat("Number of total variables (union):", length(all_vars), "\n")

extra_vars <- setdiff(all_vars, common_vars)
cat("Extra variables (only in some datasets):\n")
print(extra_vars)

## 5. Inspect which datasets contain the extra variables
lapply(list(ICD_all=df1, ICD=df2, NonICD_preserved=df3, NonICD_reduced=df4), 
       function(d) intersect(names(d), extra_vars))

## 6. Standardise column sets so they can be bound together
standardise <- function(df, all_cols) {
  for (c in setdiff(all_cols, names(df))) {
    df[[c]] <- NA  # add missing columns as NA
  }
  df <- df[, all_cols]  # reorder columns
  return(df)
}

df1 <- standardise(df1, all_vars)
df2 <- standardise(df2, all_vars)
df3 <- standardise(df3, all_vars)
df4 <- standardise(df4, all_vars)

## 7. Combine all datasets
# df <- bind_rows(df1, df2, df3, df4)
df <- bind_rows(df2, df3, df4)

## 8. Quick check: counts per dataset
table(df$dataset)

## 9. Optional: count NAs in the extra variables (to see which datasets miss them)
colSums(is.na(df[extra_vars]))



# --- Define analytic cohort ---
df <- df %>%
  # Must have a valid baseline date
  filter(!is.na(Time_zero_Ym)) %>%
  # Survival time must be positive
  filter(Survival_time > 0)

# --- Check Status values ---
table(df$Status, useNA = "ifany")

# --- Summary counts ---
total_patients <- nrow(df)
scd_events <- df %>% filter(Status == 1) %>% nrow()
censored <- df %>% filter(Status == 0) %>% nrow()
missing_status <- df %>% filter(is.na(Status)) %>% nrow()

cat("Total included patients: ", total_patients, "\n")
cat("SCD events (Status == 1): ", scd_events, "\n")
cat("Censored (Status == 0): ", censored, "\n")
cat("Missing status (Status == NA): ", missing_status, "\n")


## =========================================================
## Core survival/time fields: coerce & standardise units
## =========================================================
df <- df %>%
  mutate(
    Survival_time = suppressWarnings(as.numeric(Survival_time)),
    Status        = suppressWarnings(as.integer(Status)),
    Time_index_MI_CHD = as.numeric(Time_index_MI_CHD),
    Time_zero_Ym = as.Date(df$Time_zero_Ym),
    Time_zero_Y = as.integer(df$Time_zero_Y)
    
  )
## =========================================================
## Start date: parse Time_zero_Ym (year–month) robustly
##    - Accepts formats: "YYYY-MM", "YYYY/MM", "YYYYMM" (char or numeric)
##    - Builds Time_zero_date as 1st day of that month
##    - Fallback: if parsing fails, use mid-year of Time_zero_Y (YYYY-07-01)
## =========================================================
## ===== REPLACE STEP 2 WITH THIS ROBUST PARSER =====
suppressPackageStartupMessages({
  library(dplyr)
  library(lubridate)
  library(stringr)
})

# Inspect a few raw values (optional)
print(head(unique(df$Time_zero_Y), 20))
print(head(unique(df$Time_zero_Ym), 20))

## 1) Monthly indicators from follow-up start date (date_SCD as Date)
df <- df %>%
  mutate(
    month_num = month(Time_zero_Ym),                               # 1..12
    month_lab = month(Time_zero_Ym, label = TRUE, abbr = TRUE),    # "Jan".."Dec"
    season = case_when(                                            # (recomputed if needed)
      month_num %in% c(12, 1, 2)  ~ "Winter",
      month_num %in% c(3,  4, 5)  ~ "Spring",
      month_num %in% c(6,  7, 8)  ~ "Summer",
      month_num %in% c(9, 10, 11) ~ "Autumn",
      TRUE ~ NA_character_
    ),
    season = factor(season, levels = c("Winter","Spring","Summer","Autumn"))
  )

## (Optional) one-hot month dummies for modeling
# df <- cbind(df, model.matrix(~ month_lab - 1, data = df))
# colnames(df)[(ncol(df)-11):ncol(df)]  # check created dummy names

## 2) Holiday periods
# 2a) Christmas/New Year window: Dec 20 – Jan 5 (inclusive)
df <- df %>%
  mutate(
    is_christmas = ifelse(
      (month(Time_zero_Ym) == 12 & day(Time_zero_Ym) >= 20) |
        (month(Time_zero_Ym) == 1  & day(Time_zero_Ym) <= 5),
      1L, 0L
    )
  )

# 2b) Easter window: Good Friday (Easter-2) – Easter Monday (Easter+1)
easter_date <- as.Date(Easter(year(df$Time_zero_Ym)))
df <- df %>%
  mutate(
    easter_start       = easter_date - days(2),
    easter_end         = easter_date + days(1),
    in_easter_window   = ifelse(Time_zero_Ym >= easter_start &
                                  Time_zero_Ym <= easter_end, 1L, 0L)
  )

## 3) Daylight Saving Time (EU; last Sunday in Mar → last Sunday in Oct)
# last day of March/October for each row's year
ld_march   <- ymd(paste(year(df$Time_zero_Ym), "03", "01", sep = "-")) %m+% months(1) - days(1)
ld_october <- ymd(paste(year(df$Time_zero_Ym), "10", "01", sep = "-")) %m+% months(1) - days(1)

# compute last Sunday (wday: Sun=1,...,Sat=7); offset to last Sunday = (wday(last_day)+6) %% 7
dst_start_date <- ld_march   - days((wday(ld_march)   + 6) %% 7)
dst_end_date   <- ld_october - days((wday(ld_october) + 6) %% 7)

df <- df %>%
  mutate(
    dst_start_date      = dst_start_date,
    dst_end_date        = dst_end_date,
    in_dst              = ifelse(Time_zero_Ym >= dst_start_date &
                                   Time_zero_Ym <  dst_end_date, 1L, 0L),
    is_dst_switch_start = ifelse(Time_zero_Ym == dst_start_date, 1L, 0L),
    is_dst_switch_end   = ifelse(Time_zero_Ym == dst_end_date,   1L, 0L)
  )

## 4) Quick sanity checks (optional)
invisible({
  cat("\n-- Month counts --\n");          print(table(df$month_lab, useNA = "ifany"))
  cat("\n-- Season counts --\n");         print(table(df$season, useNA = "ifany"))
  cat("\n-- Christmas/New Year --\n");    print(table(df$is_christmas_ny, useNA = "ifany"))
  cat("\n-- Easter window --\n");         print(table(df$in_easter_window, useNA = "ifany"))
  cat("\n-- In DST --\n");                print(table(df$in_dst, useNA = "ifany"))
})



## ========================================
## 1.2 Data Quality Checks: Missingness
## ========================================
suppressPackageStartupMessages({
  library(dplyr)
  library(naniar)     # nice tools for missing data exploration
  library(VIM)        # optional: visualisation of missing data
})

## 1) Quick % missing per variable
missing_summary <- df %>%
  summarise(across(everything(), ~ mean(is.na(.))*100)) %>%
  tidyr::pivot_longer(cols = everything(),
                      names_to = "variable", values_to = "pct_missing") %>%
  arrange(desc(pct_missing))

cat("\n--- % missing per variable ---\n")
print(missing_summary)

## 2) Focus on temporal & outcome vars
temporal_vars <- c("Time_zero_Ym", "Time_zero_Y", "Survival_time","Time_index_MI_CHD", "Status")
cat("\n--- Missingness in key temporal/outcome variables ---\n")
print(sapply(df[temporal_vars], function(x) sum(is.na(x))))

## 3) Patterns of missingness (naniar)
cat("\n--- Missingness pattern (naniar) ---\n")
naniar::miss_var_summary(df[, temporal_vars])

# Visual: which variables are missing together
naniar::vis_miss(df[, temporal_vars], sort_miss = TRUE)

## 4) (Optional) Visualise broader missingness structure (VIM)
VIM::aggr(df[temporal_vars], numbers = TRUE, sortVars = TRUE,
          labels = names(df[temporal_vars]), cex.axis = 0.7,
          gap = 3, ylab = c("Missing data","Pattern"))


## ========================================
## Data Quality Checks: Inadmissible values
## Based on Table 1
## ========================================


## Check plausible median/range against expected units
## ========================================
## Apply Table 1 inadmissible value filters
## (skips variables that don't exist)
## ========================================
## list of variables we expect from Table 1
vars_to_check <- c("BMI","SBP","DBP","BUN","Cholesterol","CRP","Haemoglobin",
                   "HbA1c","HDL","LDL","Sodium","Triglycerides","Troponin_T",
                   "TSH","HR","PR","QRS","QTc","NYHA","AV_block","AV_block_II_or_III")


## which are present?
present_vars <- vars_to_check[vars_to_check %in% names(df)]
missing_vars <- setdiff(vars_to_check, names(df))

cat("\n Present in df:\n")
print(present_vars)

cat("\n Missing in df:\n")
print(missing_vars)

## ===========================================
## Unit and plausibility checks
## ===========================================
# function to summarise ranges
qc_summary <- function(x) {
  c(min = suppressWarnings(min(x, na.rm=TRUE)),
    median = suppressWarnings(median(x, na.rm=TRUE)),
    max = suppressWarnings(max(x, na.rm=TRUE)),
    missing = sum(is.na(x)))
}

qc_results <- sapply(vars_to_check[vars_to_check %in% names(df)], 
                     function(v) qc_summary(df[[v]]))
t(qc_results)

df <- df %>%
  mutate(
    BMI = ifelse(BMI < 12 | BMI > 69, NA, BMI),
    SBP = NA, DBP = NA, CRP = NA, Troponin_T = NA,
    NYHA = NA, AV_block = NA, AV_block_II_III = NA,
    BUN = ifelse(BUN > 900, NA, BUN),
    Cholesterol = ifelse(Cholesterol <= 50, NA, Cholesterol),
    Haemoglobin = ifelse(Haemoglobin < 2 | Haemoglobin > 110, NA, Haemoglobin),
    HbA1c = ifelse(HbA1c < 2.5, NA, HbA1c),
    HDL = ifelse(HDL == 0, NA, HDL),
    LDL = ifelse(LDL == 0, NA, LDL),
    Sodium = ifelse(Sodium < 99, NA, Sodium),
    Triglycerides = ifelse(Triglycerides < 20, NA, Triglycerides),
    TSH = ifelse(TSH == 0, NA, TSH),
    HR = ifelse(HR < 25 | HR > 140, NA, HR),
    PR = ifelse(PR <= 50 | PR > 1000, NA, PR),
    QRS = ifelse(QRS < 50, NA, QRS),
    QTc = ifelse(QTc <= 250 | QTc > 790, NA, QTc)
  )


##############################################################
#   TEMPORAL CONSISTENCY CHECKS
#   For: Exploring Seasonal Patterns in SCD Post-MI
#   Goal: Ensure follow-up start → event dates follow logic
##############################################################

library(lubridate)
library(dplyr)

##############################################################
# 2. Create EVENT/CENSORING DATE
#    Survival_time is in MONTHS -> convert to a date
##############################################################

# Round to nearest month (some values include decimals)
df$event_date <- df$Time_zero_Ym %m+% months(round(df$Survival_time))

##############################################################
# 3. TEMPORAL CONSISTENCY CHECKS
##############################################################

### 3.1 Check that event_date occurs AFTER or ON follow-up start
df$check_event_after_start <- df$event_date >= df$Time_zero_Ym

cat("\n--- Check: event_date >= followup_start_date ---\n")
print(table(df$check_event_after_start, useNA = "ifany"))

### 3.2 Check for negative survival times
cat("\n--- Negative Survival Time (should be 0) ---\n")
print(sum(df$Survival_time < 0, na.rm = TRUE))

### 3.3 Check for missing dates
cat("\n--- Missing Dates ---\n")
print(colSums(is.na(df[, c("Time_zero_Ym", "event_date")])))
      
### 3.4 People who DIED must have event_date present
cat("\n--- Deaths with missing event_date ---\n")
missing_event_for_death <- df %>%
  filter(Status %in% c(1,2) & is.na(event_date)) %>%
  nrow()
print(missing_event_for_death)

### 3.5 Same-day events (survival_time = 0)
cat("\n--- Number of same-day events (survival_time = 0) ---\n")
print(df %>% filter(Survival_time == 0) %>% nrow())

##############################################################
# 4. Extract rows with temporal inconsistencies
##############################################################

temporal_issues <- df %>%
  filter(
    event_date < followup_start_date | 
    Survival_time < 0
  )

cat("\n--- Temporal Issues Found ---\n")
print(nrow(temporal_issues))

# Show first few problematic rows
print(head(temporal_issues))

##############################################################
# 5. Summary report text (printed)
##############################################################

cat("\n================================================\n")
cat("TEMPORAL CONSISTENCY SUMMARY\n")
cat("================================================\n")

cat("Total records: ", nrow(df), "\n")
cat("Records with event_date BEFORE followup_start_date: ", 
    sum(df$event_date < df$followup_start_date, na.rm = TRUE), "\n")
cat("Negative survival times: ", sum(df$Survival_time < 0, na.rm = TRUE), "\n")
cat("Deaths with missing event_date: ", missing_event_for_death, "\n")
cat("Same-day events (survival_time = 0): ", 
    df %>% filter(Survival_time == 0) %>% nrow(), "\n")
cat("================================================\n")


## ==========================================================
## Lifestyle Factors Extraction
## ==========================================================

df <- df %>%
  mutate(
    # Smoking status (binary: 1 = Yes, 0 = No)
    lifestyle_smoking = ifelse(Smoking == "Yes", 1,
                               ifelse(Smoking == "No", 0, NA)),
    
    # Alcohol consumption (binary: 1 = Yes, 0 = No)
    lifestyle_alcohol = ifelse(Alcohol == "Yes", 1,
                               ifelse(Alcohol == "No", 0, NA)),
    
    # BMI categories (WHO standard)
    lifestyle_BMIcat = case_when(
      !is.na(BMI) & BMI < 18.5 ~ "Underweight",
      BMI >= 18.5 & BMI < 25   ~ "Normal",
      BMI >= 25   & BMI < 30   ~ "Overweight",
      BMI >= 30   & BMI < 35   ~ "Obese I",
      BMI >= 35   & BMI < 40   ~ "Obese II",
      BMI >= 40               ~ "Obese III",
      TRUE ~ NA_character_
    ),
    
    # Physical activity proxy using NYHA functional class
    lifestyle_physicalactivity = case_when(
      NYHA == "I" ~ "No limitation",
      NYHA == "II" ~ "Mild limitation",
      NYHA %in% c("III", "IV") ~ "Severe limitation",
      TRUE ~ NA_character_
    )
  )

## Quick descriptive check
summary(df$BMI)                       # continuous distribution
table(df$lifestyle_BMIcat, useNA="ifany")
table(df$lifestyle_smoking, useNA="ifany")
table(df$lifestyle_alcohol, useNA="ifany")
table(df$lifestyle_physicalactivity, useNA="ifany")

# ============================================================
# TASK: Clinical Risk Factors Harmonisation
# ============================================================
# ------------------------------
# 1. Standardise LVEF
# ------------------------------
df <- df %>%
mutate(
  # Harmonised LVEF: prefer MRI if available, otherwise ECHO
  LVEF_std = case_when(
    !is.na(MRI_LVEF) ~ MRI_LVEF,
    !is.na(LVEF) ~ LVEF,
    TRUE ~ NA_real_
  ),

  # Flag for source modality
  LVEF_source = case_when(
    !is.na(MRI_LVEF) ~ "MRI",
    !is.na(LVEF) ~ "ECHO",
    TRUE ~ NA_character_
  )
)
#   
#   # ------------------------------
# # 2. Heart Rate measures (ECG)
# # ------------------------------
# ============================================================
# Task: Clinical Risk Factors - Standardise Heart Rate Variability Measures
# ============================================================
# Available variables:
# - HR (bpm), PR (ms), QRS (ms), QTc (ms)
# - AV_block, AV_block_II_or_III, LBBB, RBBB (categorical Yes/No)
#
# Goal:
# 1. Standardise continuous ECG variables (z-scores).
# 2. Convert binary categorical conduction abnormalities to numeric (0/1).
# 3. Create an optional composite indicator for any conduction abnormality.
# ============================================================

df <- df %>%
  mutate(
    # ---- (1) Standardise continuous ECG variables ----
    HR_std   = ifelse(!is.na(HR), scale(HR), NA),      # Heart rate
    PR_std   = ifelse(!is.na(PR), scale(PR), NA),      # PR interval
    QRS_std  = ifelse(!is.na(QRS), scale(QRS), NA),    # QRS duration
    QTc_std  = ifelse(!is.na(QTc), scale(QTc), NA),    # QTc interval
    
    # ---- (2) Convert binary categorical conduction abnormalities to 0/1 ----
    AV_block_bin        = ifelse(AV_block == "Yes", 1,
                                 ifelse(AV_block == "No", 0, NA)),
    AV_block_II_III_bin = ifelse(AV_block_II_or_III == "Yes", 1,
                                 ifelse(AV_block_II_or_III == "No", 0, NA)),
    LBBB_bin            = ifelse(LBBB == "Yes", 1,
                                 ifelse(LBBB == "No", 0, NA)),
    RBBB_bin            = ifelse(RBBB == "Yes", 1,
                                 ifelse(RBBB == "No", 0, NA)),
    
    # ---- (3) Optional: Composite conduction abnormality flag ----
    conduction_abnormality = ifelse(
      AV_block_bin == 1 | AV_block_II_III_bin == 1 | 
        LBBB_bin == 1 | RBBB_bin == 1, 1, 0
    )
  )

# ============================================================
# After running this:
# - HR_std, PR_std, QRS_std, QTc_std = standardised ECG values
# - *_bin = recoded binary flags for conduction abnormalities
# - conduction_abnormality = 1 if any abnormality present, else 0
# ============================================================

  
  # ------------------------------
# 3. Medication Adherence Indicators
# ------------------------------
med_vars <- c("ACE_inhibitor_ARB", "ACE_inhibitor", "ARB", "Aldosterone_antagonist",
              "Anti_anginal", "Anti_arrhythmic_III", "Anti_coagulant", 
              "Anti_diabetic", "Anti_diabetic_insulin", "Anti_diabetic_oral",
              "Anti_platelet", "Beta_blockers", "Calcium_antagonists", 
              "Digitalis_glycosides", "Diuretics", "Lipid_lowering")

# Recode all Yes/No to 1/0
df <- df %>%
  mutate(across(all_of(med_vars), 
                ~ ifelse(. == "Yes", 1, 
                         ifelse(. == "No", 0, NA)), 
                .names = "{.col}_bin"))


# Create adherence score
# Count of medications taken (simple adherence proxy)
# 
# Or percentage adherence (number of Yes / total available)

df <- df %>%
  mutate(
    med_adherence_score = rowSums(across(ends_with("_bin")), na.rm = TRUE),
    med_adherence_pct   = med_adherence_score / length(med_vars) * 100
  )



# ============================================================
# ✅ Quick descriptive checks
# ============================================================

# ============================================================
# Task: Descriptive Checks for Standardised Clinical Risk Factors
# ============================================================
# This step checks the distributions and coding of the newly created
# standardised / harmonised clinical risk factor variables.
# ============================================================


# ---- 1. Check LVEF standardisation (if you created LVEF_final before) ----
df %>%
  summarise(
    n_total    = n(),
    n_missing  = sum(is.na(LVEF_std)),
    mean_LVEF  = mean(LVEF_std, na.rm = TRUE),
    sd_LVEF    = sd(LVEF_std, na.rm = TRUE),
    min_LVEF   = min(LVEF_std, na.rm = TRUE),
    max_LVEF   = max(LVEF_std, na.rm = TRUE)
  )

# ---- 2. Summary of standardised ECG continuous variables ----
df %>%
  summarise(
    mean_HR   = mean(HR, na.rm = TRUE),
    mean_HR_z = mean(HR_std, na.rm = TRUE),  # should be ~0
    sd_HR_z   = sd(HR_std, na.rm = TRUE),    # should be ~1
    
    mean_PR_z = mean(PR_std, na.rm = TRUE),
    sd_PR_z   = sd(PR_std, na.rm = TRUE),
    
    mean_QRS_z = mean(QRS_std, na.rm = TRUE),
    sd_QRS_z   = sd(QRS_std, na.rm = TRUE),
    
    mean_QTc_z = mean(QTc_std, na.rm = TRUE),
    sd_QTc_z   = sd(QTc_std, na.rm = TRUE)
  )

# ---- 3. Frequency tables for binary conduction abnormalities ----
binary_vars <- c("AV_block_bin", "AV_block_II_III_bin", "LBBB_bin", "RBBB_bin", "conduction_abnormality")

for (var in binary_vars) {
  print(var)
  print(table(df[[var]], useNA = "ifany"))
}

# ---- 4. Quick barplot of conduction abnormalities ----
library(ggplot2)
dev.off()

df %>%
  select(all_of(binary_vars)) %>%
  tidyr::pivot_longer(cols = everything(), names_to = "variable", values_to = "value") %>%
  filter(!is.na(value)) %>%
  group_by(variable, value) %>%
  summarise(n = n(), .groups = "drop") %>%
  ggplot(aes(x = variable, y = n, fill = factor(value))) +
  geom_col(position = "dodge") +
  labs(title = "Frequency of Conduction Abnormalities", y = "Count", x = "") +
  scale_fill_manual(values = c("0" = "grey70", "1" = "steelblue"),
                    name = "Value (0=No, 1=Yes)") +
  theme_minimal()


# Save as an RDS file (keeps all data types intact)

saveRDS(df, file = "df_cleaned.rds")

