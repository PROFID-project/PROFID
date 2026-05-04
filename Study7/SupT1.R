###############################################
# Supplementary table 1: Reporting completeness of binary clinical variables by centre (ICD and Non-ICD cohorts)
###############################################

# Packages
req <- c("dplyr","tidyverse","data.table","gt")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Input files
df_ICD <- fread("T:/Data Transfer to Charite/raw/ICD_filtered_with_coords.csv")
df_NR <- fread("T:/Data Transfer to Charite/raw/NonICD_reduced_filtered_with_coords.csv")
df_NP <- fread("T:/Data Transfer to Charite/raw/NonICD_preserved_filtered_with_coords.csv")

# Data preparation
df_ICD <- df_ICD %>%
  filter(CVD_risk_region != 2)

df_ICD$CVD_risk_region <- factor(df_ICD$CVD_risk_region, levels = c(1, 3, 4), labels = c("Low", "High", "Very High")) 

df_NR$CVD_risk_region <- factor(df_NR$CVD_risk_region, levels = c(1, 2, 3, 4), labels = c("Low", "Medium", "High", "Very High"))

df_NP$CVD_risk_region <- factor(df_NP$CVD_risk_region, levels = c(1, 2, 3, 4), labels = c("Low", "Medium", "High", "Very High"))

binary_vars <- c(
  "ACE_inhibitor_ARB",
  "Beta_blockers",
  "Anti_platelet",
  "Lipid_lowering",
  "Anti_diabetic",
  "Anti_coagulant",
  "Diuretics",
  "Smoking",
  "Hypertension"
)

# Reporting completeness of binary clinical variables by centre (ICD cohort)
missing_table_ICD <- df_ICD %>%
  group_by(CVD_risk_region, ctr_name) %>%
  summarise(
    N = n(),
    across(
      all_of(binary_vars),
      ~ N - sum(. %in% c("Yes","No"), na.rm = TRUE)
    ),
    .groups = "drop"
  )
gt_missing_ICD <- missing_table_ICD %>%
  gt() %>%
  tab_header(
    title = md("**Reporting Completeness of Binary Clinical Variables by Centre**"),
    subtitle = md("**ICD Cohort**")
  ) %>%
  cols_label(
    CVD_risk_region = md("**CVD Risk Region**"),
    ctr_name = md("**Centre**"),
    N = md("**n**"),
    ACE_inhibitor_ARB = md("**ACEi/ARB**"),
    Beta_blockers = md("**Beta-blocker**"),
    Anti_platelet = md("**Antiplatelet**"),
    Lipid_lowering = md("**Statin**"),
    Anti_diabetic = md("**Anti-diabetic**"),
    Anti_coagulant = md("**Anti-coagulant**"),
    Diuretics = md("**Diuretics**"),
    Smoking = md("**Smoking Status**"),
    Hypertension = md("**Hypertension Status**")
  )
gtsave(gt_missing_ICD, "med_missingness_ICD.html")

# Reporting completeness of binary clinical variables by centre (Non-ICD reduced cohort - ≤35%)
missing_table_NR <- df_NR %>%
  group_by(CVD_risk_region, DB) %>%
  summarise(
    N = n(),
    across(
      all_of(binary_vars),
      ~ N - sum(. %in% c("Yes","No"), na.rm = TRUE)
    ),
    .groups = "drop"
  )
gt_missing_NR <- missing_table_NR %>%
  gt() %>%
  tab_header(
    title = md("**Reporting Completeness of Binary Clinical Variables by Centre**"),
    subtitle = md("**Non-ICD Cohort (≤35%)**")
  ) %>%
  cols_label(
    CVD_risk_region = md("**CVD Risk Region**"),
    DB = md("**Centre**"),
    N = md("**n**"),
    ACE_inhibitor_ARB = md("**ACEi/ARB**"),
    Beta_blockers = md("**Beta-blocker**"),
    Anti_platelet = md("**Antiplatelet**"),
    Lipid_lowering = md("**Statin**"),
    Anti_diabetic = md("**Anti-diabetic**"),
    Anti_coagulant = md("**Anti-coagulant**"),
    Diuretics = md("**Diuretics**"),
    Smoking = md("**Smoking Status**"),
    Hypertension = md("**Hypertension Status**")
  )
gtsave(gt_missing_NR, "med_missingness_NR.html")

# Reporting completeness of binary clinical variables by centre (Non-ICD preserved cohort - >35%)
missing_table_NP <- df_NP %>%
  group_by(CVD_risk_region, DB) %>%
  summarise(
    N = n(),
    across(
      all_of(binary_vars),
      ~ N - sum(. %in% c("Yes","No"), na.rm = TRUE)
    ),
    .groups = "drop"
  )
gt_missing_NP <- missing_table_NP %>%
  gt() %>%
  tab_header(
    title = md("**Reporting Completeness of Binary Clinical Variables by Centre**"),
    subtitle = md("**Non-ICD Cohort (>35%)**")
  ) %>%
  cols_label(
    CVD_risk_region = md("**CVD Risk Region**"),
    DB = md("**Centre**"),
    N = md("**n**"),
    ACE_inhibitor_ARB = md("**ACEi/ARB**"),
    Beta_blockers = md("**Beta-blocker**"),
    Anti_platelet = md("**Antiplatelet**"),
    Lipid_lowering = md("**Statin**"),
    Anti_diabetic = md("**Anti-diabetic**"),
    Anti_coagulant = md("**Anti-coagulant**"),
    Diuretics = md("**Diuretics**"),
    Smoking = md("**Smoking Status**"),
    Hypertension = md("**Hypertension Status**")
  )
gtsave(gt_missing_NP, "med_missingness_NP.html")