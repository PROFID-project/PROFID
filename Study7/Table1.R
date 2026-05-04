###############################################
# Table 1: Baseline characteristics by cardiovascular risk region (ICD and Non-ICD cohorts)
###############################################

# Packages
req <- c("dplyr","tidyverse","data.table","gt","gtsummary")
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

binary_vars <- c(
  "Hypertension",
  "Diabetes",
  "ACE_inhibitor_ARB",
  "Beta_blockers",
  "Anti_platelet",
  "Lipid_lowering"
)

df_ICD$CVD_risk_region <- factor(df_ICD$CVD_risk_region, levels = c(1, 3, 4), labels = c("Low", "High", "Very High")) 
df_ICD <- df_ICD %>%
  mutate(across(all_of(binary_vars), ~factor(.x, levels = c("No","Yes"))))
df_ICD <- df_ICD %>%
  mutate(Sex = factor(Sex))

df_NR$CVD_risk_region <- factor(df_NR$CVD_risk_region, levels = c(1, 2, 3, 4), labels = c("Low", "Medium", "High", "Very High"))
df_NR <- df_NR %>%
  mutate(across(all_of(binary_vars), ~factor(.x, levels = c("No","Yes"))))
df_NR <- df_NR %>%
  mutate(Sex = factor(Sex))

df_NP$CVD_risk_region <- factor(df_NP$CVD_risk_region, levels = c(1, 2, 3, 4), labels = c("Low", "Medium", "High", "Very High"))
df_NP <- df_NP %>%
  mutate(across(all_of(binary_vars), ~factor(.x, levels = c("No","Yes"))))
df_NP <- df_NP %>%
  mutate(Sex = factor(Sex))

# Table 1: Baseline characteristics by cardiovascular risk region (ICD cohort)
table1_ICD <- df_ICD %>%
  select(
    CVD_risk_region,
    Age,
    Sex,
    LVEF,
    Hypertension,
    Diabetes,
    ACE_inhibitor_ARB,
    Beta_blockers,
    Anti_platelet,
    Lipid_lowering
  ) %>%
  tbl_summary(
    by = CVD_risk_region,
    type = list(
      Age ~ "continuous",
      LVEF ~ "continuous"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      Age ~ "Age (years)",
      Sex ~ "Sex",
      LVEF ~ "LVEF (%)",
      Hypertension ~ "Hypertension",
      Diabetes ~ "Diabetes mellitus",
      ACE_inhibitor_ARB ~ "ACEi / ARB",
      Beta_blockers ~ "Beta-blocker",
      Anti_platelet ~ "Antiplatelet",
      Lipid_lowering ~ "Statin therapy"
    ),
    digits = all_continuous() ~ 1,
    missing = "always",
    missing_text = "Missing"
  ) %>%
  add_n() %>%
  modify_header(label = "**Variable**") %>%
  bold_labels() %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        category = case_when(
          variable %in% c("Age","Sex","LVEF") ~ "Demographics",
          variable %in% c("Hypertension","Diabetes") ~ "Clinical characteristics",
          TRUE ~ "Medications"
        )
      )
  )
table1_ICD_gt <- table1_ICD %>%
  as_gt() %>%
  tab_header(
    title = md("**Baseline Characteristics by Cardiovascular Risk Region**"),
    subtitle = md("ICD Cohort")
  ) 
gtsave(table1_ICD_gt, "table1_ICD.html")

# Table 1: Baseline characteristics by cardiovascular risk region (Non-ICD reduced cohort - ≤35%)
table1_NR <- df_NR %>%
  select(
    CVD_risk_region,
    Age,
    Sex,
    LVEF,
    Hypertension,
    Diabetes,
    ACE_inhibitor_ARB,
    Beta_blockers,
    Anti_platelet,
    Lipid_lowering
  ) %>%
  tbl_summary(
    by = CVD_risk_region,
    type = list(
      Age ~ "continuous",
      LVEF ~ "continuous"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      Age ~ "Age (years)",
      Sex ~ "Sex",
      LVEF ~ "LVEF (%)",
      Hypertension ~ "Hypertension",
      Diabetes ~ "Diabetes mellitus",
      ACE_inhibitor_ARB ~ "ACEi / ARB",
      Beta_blockers ~ "Beta-blocker",
      Anti_platelet ~ "Antiplatelet",
      Lipid_lowering ~ "Statin therapy"
    ),
    digits = all_continuous() ~ 1,
    missing = "always",
    missing_text = "Missing"
  ) %>%
  add_n() %>%
  modify_header(label = "**Variable**") %>%
  bold_labels() %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        category = case_when(
          variable %in% c("Age","Sex","LVEF") ~ "Demographics",
          variable %in% c("Hypertension","Diabetes") ~ "Clinical characteristics",
          TRUE ~ "Medications"
        )
      )
  )
table1_NR_gt <- table1_NR %>%
  as_gt() %>%
  tab_header(
    title = md("**Baseline Characteristics by Cardiovascular Risk Region**"),
    subtitle = md("Non-ICD Cohort (≤35%)")
  ) 
gtsave(table1_NR_gt, "table1_NR.html")

# Table 1: Baseline characteristics by cardiovascular risk region (Non-ICD preserved cohort - >35%)
table1_NP <- df_NP %>%
  select(
    CVD_risk_region,
    Age,
    Sex,
    LVEF,
    Hypertension,
    Diabetes,
    ACE_inhibitor_ARB,
    Beta_blockers,
    Anti_platelet,
    Lipid_lowering
  ) %>%
  tbl_summary(
    by = CVD_risk_region,
    type = list(
      Age ~ "continuous",
      LVEF ~ "continuous"
    ),
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    label = list(
      Age ~ "Age (years)",
      Sex ~ "Sex",
      LVEF ~ "LVEF (%)",
      Hypertension ~ "Hypertension",
      Diabetes ~ "Diabetes mellitus",
      ACE_inhibitor_ARB ~ "ACEi / ARB",
      Beta_blockers ~ "Beta-blocker",
      Anti_platelet ~ "Antiplatelet",
      Lipid_lowering ~ "Statin therapy"
    ),
    digits = all_continuous() ~ 1,
    missing = "always",
    missing_text = "Missing"
  ) %>%
  add_n() %>%
  modify_header(label = "**Variable**") %>%
  bold_labels() %>%
  modify_table_body(
    ~ .x %>%
      mutate(
        category = case_when(
          variable %in% c("Age","Sex","LVEF") ~ "Demographics",
          variable %in% c("Hypertension","Diabetes") ~ "Clinical characteristics",
          TRUE ~ "Medications"
        )
      )
  )
table1_NP_gt <- table1_NP %>%
  as_gt() %>%
  tab_header(
    title = md("**Baseline Characteristics by Cardiovascular Risk Region**"),
    subtitle = md("Non-ICD Cohort (>35%)")
  ) 
gtsave(table1_NP_gt, "table1_NP.html")