###############################################
# Sudden Cardiac Death Event Rates by Cardiovascular Risk Region (ICD and Non-ICD cohorts)
# Supplementary table 2: Outcome reporting completeness by region and centre (ICD and Non-ICD cohorts)
###############################################

# Packages
req <- c("dplyr","tidyverse","data.table","gt","ggplot2","scales")
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

# Event rates by cardiovascular risk region (ICD cohort)
table_events_ICD <- df_ICD %>%
  filter(!is.na(Status)) %>%
  group_by(CVD_risk_region) %>%
  summarise(
    N = n(),
    SCD_events = sum(Status == 1, na.rm = TRUE),
    Other_deaths = sum(Status == 2, na.rm = TRUE),
    SCD_rate = round((SCD_events / N) * 100, 1),
    Other_death_rate = round((Other_deaths / N) * 100, 1),
    .groups = "drop"
  )
events_ICD <- table_events_ICD %>%
  gt() %>%
  tab_header(
    title = md("**Sudden Cardiac Death Events by Cardiovascular Risk Region**"),
    subtitle = md("ICD Cohort")
  ) %>%
  cols_label(
    CVD_risk_region = "Region",
    N = "Patients (N)",
    SCD_events = "SCD Events",
    Other_deaths = "Other Deaths",
    SCD_rate = "SCD Rate (%)",
    Other_death_rate = "Other Death Rate (%)"
  )
gtsave(events_ICD, "events_ICD.html")
# Outcome reporting completeness by region and centre (ICD cohort)
missing_status_multilevel_ICD <- df_ICD %>%
  group_by(CVD_risk_region, ctr_name) %>%
  summarise(
    N = n(),
    missing_status = sum(is.na(Status)),
    missing_pct = round((missing_status / N) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(CVD_risk_region, desc(missing_status))
events_missing_ICD <- missing_status_multilevel_ICD %>%
  gt(groupname_col = "CVD_risk_region") %>%
  cols_label(
    ctr_name = "Centre",
    N = "Patients (N)",
    missing_status = "Missing Outcomes",
    missing_pct = "Missing (%)"
  ) %>%
  tab_header(
    title = md("**Outcome Reporting Completeness by Region and Centre**"),
    subtitle = md("ICD Cohort")
  )
gtsave(events_missing_ICD, "missing_status_multilevel_ICD.html")

# Event rates by cardiovascular risk region (Non-ICD reduced cohort - ≤35%)
table_events_NR <- df_NR %>%
  filter(!is.na(Status)) %>%
  group_by(CVD_risk_region) %>%
  summarise(
    N = n(),
    SCD_events = sum(Status == 1, na.rm = TRUE),
    Other_deaths = sum(Status == 2, na.rm = TRUE),
    SCD_rate = round((SCD_events / N) * 100, 1),
    Other_death_rate = round((Other_deaths / N) * 100, 1),
    .groups = "drop"
  )
events_NR <- table_events_NR %>%
  gt() %>%
  tab_header(
    title = md("**Sudden Cardiac Death Events by Cardiovascular Risk Region**"),
    subtitle = md("Non-ICD Cohort (≤35%)")
  ) %>%
  cols_label(
    CVD_risk_region = "Region",
    N = "Patients (N)",
    SCD_events = "SCD Events",
    Other_deaths = "Other Deaths",
    SCD_rate = "SCD Rate (%)",
    Other_death_rate = "Other Death Rate (%)"
  )
gtsave(events_NR, "events_NR.html")
# Outcome reporting completeness by region and centre (Non-ICD reduced cohort - ≤35%)
missing_status_multilevel_NR <- df_NR %>%
  group_by(CVD_risk_region, DB) %>%
  summarise(
    N = n(),
    missing_status = sum(is.na(Status)),
    missing_pct = round((missing_status / N) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(CVD_risk_region, desc(missing_status))
events_missing_NR <- missing_status_multilevel_NR %>%
  gt(groupname_col = "CVD_risk_region") %>%
  cols_label(
    DB = "Centre",
    N = "Patients (N)",
    missing_status = "Missing Outcomes",
    missing_pct = "Missing (%)"
  ) %>%
  tab_header(
    title = md("**Outcome Reporting Completeness by Region and Centre**"),
    subtitle = md("Non-ICD Cohort (≤35%)")
  )
gtsave(events_missing_NR, "missing_status_multilevel_NR.html")

# Event rates by cardiovascular risk region (Non-ICD preserved cohort - >35%)
table_events_NP <- df_NP %>%
  filter(!is.na(Status)) %>%
  group_by(CVD_risk_region) %>%
  summarise(
    N = n(),
    SCD_events = sum(Status == 1, na.rm = TRUE),
    Other_deaths = sum(Status == 2, na.rm = TRUE),
    SCD_rate = round((SCD_events / N) * 100, 1),
    Other_death_rate = round((Other_deaths / N) * 100, 1),
    .groups = "drop"
  )
events_NP <- table_events_NP %>%
  gt() %>%
  tab_header(
    title = md("**Sudden Cardiac Death Events by Cardiovascular Risk Region**"),
    subtitle = md("Non-ICD Cohort (>35%)")
  ) %>%
  cols_label(
    CVD_risk_region = "Region",
    N = "Patients (N)",
    SCD_events = "SCD Events",
    Other_deaths = "Other Deaths",
    SCD_rate = "SCD Rate (%)",
    Other_death_rate = "Other Death Rate (%)"
  )
gtsave(events_NP, "events_NP.html")
# Outcome reporting completeness by region and centre (Non-ICD preserved cohort - >35%)
missing_status_multilevel_NP <- df_NP %>%
  group_by(CVD_risk_region, DB) %>%
  summarise(
    N = n(),
    missing_status = sum(is.na(Status)),
    missing_pct = round((missing_status / N) * 100, 1),
    .groups = "drop"
  ) %>%
  arrange(CVD_risk_region, desc(missing_status))
events_missing_NP <- missing_status_multilevel_NP %>%
  gt(groupname_col = "CVD_risk_region") %>%
  cols_label(
    DB = "Centre",
    N = "Patients (N)",
    missing_status = "Missing Outcomes",
    missing_pct = "Missing (%)"
  ) %>%
  tab_header(
    title = md("**Outcome Reporting Completeness by Region and Centre**"),
    subtitle = md("Non-ICD Cohort (>35%)")
  )
gtsave(events_missing_NP, "missing_status_multilevel_NP.html")