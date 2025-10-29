
install.packages("survival")
install.packages("survminer")
install.packages("dplyr")
install.packages("svglite")
install.packages("rms")

library(survival)
library(survminer)
library(dplyr)
library(rms)

combined <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/combined_dataset.csv")

data_dir <- "S:/AG/f-dhzc-profid/Data Transfer to Charite"
setwd("T:/Dokumente/PROFID/Study6")

missing_summary <- combined %>%
  select(all_of(vars_base)) %>%
  summarise(across(everything(),
                   ~ mean(is.na(.)) * 100))

t(missing_summary)

impute_vars <- c(
  # Outcome
  "Survival_time", "Status",
  
  # Exposure
  "BMI",
  
  # Base confounders
  "Age", "Sex",
  "Diabetes", "Hypertension", "Smoking", "MI_history",
  "LVEF", "NYHA",
  "eGFR", "Haemoglobin",
  "ACE_inhibitor_ARB", "Beta_blockers", "Lipid_lowering",
  "Revascularisation_acute",
  
  #  auxiliary variables
  "NTProBNP", "CRP", "Troponin_T", "HbA1c", "Cholesterol", "LDL", "HDL", "Triglycerides",
  "LV_mass", "MRI_LVEF", "Infarct_size", "Greyzone_size",
  "COPD", "Cancer", "Stroke_TIA", "HF",
  "Anti_platelet", "Anti_coagulant", "Diuretics", "ACE_inhibitor", "ARB",
  "PCI_acute", "CABG_acute", "Thrombolysis_acute",
  "AV_block", "NSVT", "LBBB", "RBBB", "AF_atrial_flutter",
  "Baseline_type", "HasMRI", "Center", "CVD_risk_region", "Baseline_type"
)



imp_data <- combined %>% select(all_of(impute_vars))
imp <- mice(imp_data, m = 20, maxit = 20, seed = 123)

