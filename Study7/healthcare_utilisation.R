###############################################
# Healthcare Utilisation (ICD and Non-ICD cohorts)
# - Table generation
###############################################

# Packages
req <- c("data.table","dplyr","tidyverse","gt")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Input file
df_ICD <- fread("T:/Data Transfer to Charite/imp/ICD_imputed.csv")
df_NR <- fread("T:/Data Transfer to Charite/imp/NonICD_reduced_imputed.csv")
df_NP <- fread("T:/Data Transfer to Charite/imp/NonICD_preserved_imputed.csv")

# Healthcare utilisation patterns
util_ICD <- df_ICD %>%
  group_by(CVD_risk_region) %>%
  summarise(
    n = n(),
    PCI = mean(PCI == 1, na.rm = TRUE) * 100,
    CABG = mean(CABG == 1, na.rm = TRUE) * 100,
    ACE_ARB = mean(ACE_inhibitor_ARB == 1, na.rm = TRUE) * 100,
    Beta_blocker = mean(Beta_blockers == 1, na.rm = TRUE) * 100,
    Diuretics = mean(Diuretics == 1, na.rm = TRUE) * 100,
    Antiplatelet = mean(Anti_platelet == 1, na.rm = TRUE) * 100,
    Lipid_lowering = mean(Lipid_lowering == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(Cohort = "ICD")

util_NR <- df_NR %>%
  group_by(CVD_risk_region) %>%
  summarise(
    n = n(),
    PCI = mean(PCI == 1, na.rm = TRUE) * 100,
    CABG = mean(CABG == 1, na.rm = TRUE) * 100,
    ACE_ARB = mean(ACE_inhibitor_ARB == 1, na.rm = TRUE) * 100,
    Beta_blocker = mean(Beta_blockers == 1, na.rm = TRUE) * 100,
    Diuretics = mean(Diuretics == 1, na.rm = TRUE) * 100,
    Antiplatelet = mean(Anti_platelet == 1, na.rm = TRUE) * 100,
    Lipid_lowering = mean(Lipid_lowering == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(Cohort = "Non-ICD (≤35%)")

util_NP <- df_NP %>%
  group_by(CVD_risk_region) %>%
  summarise(
    n = n(),
    PCI = mean(PCI == 1, na.rm = TRUE) * 100,
    CABG = mean(CABG == 1, na.rm = TRUE) * 100,
    ACE_ARB = mean(ACE_inhibitor_ARB == 1, na.rm = TRUE) * 100,
    Beta_blocker = mean(Beta_blockers == 1, na.rm = TRUE) * 100,
    Diuretics = mean(Diuretics == 1, na.rm = TRUE) * 100,
    Antiplatelet = mean(Anti_platelet == 1, na.rm = TRUE) * 100,
    Lipid_lowering = mean(Lipid_lowering == 1, na.rm = TRUE) * 100,
    .groups = "drop"
  ) %>%
  mutate(Cohort = "Non-ICD (>35%)")

util_stacked <- bind_rows(util_ICD, util_NR, util_NP) %>%
  mutate(across(where(is.numeric), ~ round(.x, 1)))

# GT table
gt_util_stacked <- gt(
  util_stacked,
  groupname_col = "Cohort"
) %>%
  cols_label(
    CVD_risk_region = "Region",
    n = md("**n**"),
    PCI = "PCI (%)",
    CABG = "CABG (%)",
    ACE_ARB = "ACEi/ARB (%)",
    Beta_blocker = "Beta-blocker (%)",
    Diuretics = "Diuretics (%)",
    Antiplatelet = "Antiplatelet (%)",
    Lipid_lowering = "Lipid-lowering (%)"
  ) %>%
  cols_align(
    align = "center",
    columns = c(n, PCI, CABG, ACE_ARB, Beta_blocker,
                Diuretics, Antiplatelet, Lipid_lowering)
  ) %>%
  tab_header(
    title = md("**Healthcare utilisation patterns by region**"),
    subtitle = md("Stacked comparison of ICD and Non-ICD cohorts after imputation")
  ) %>%
  tab_source_note(
    md(
      "Values are percentages of patients within each region. 
      Binary variables coded as 0/1; percentages calculated as mean(x == 1)."
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  tab_options(
    table.font.size = px(13),
    data_row.padding = px(4),
    row_group.padding = px(6)
  )
gtsave(gt_util_stacked, "T:/Data Transfer to Charite/imp/imp_healthcare_utilisation_by_region.html")
