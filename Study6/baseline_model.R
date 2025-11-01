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


vars_base <- c(
  "Survival_time", "Status", "BMI",
  "Age", "Sex", "Diabetes", "Hypertension", "Smoking", "MI_history",
  "LVEF", "NYHA",
  "eGFR", "Haemoglobin",
  "ACE_inhibitor_ARB", "Beta_blockers", "Lipid_lowering",
  "Revascularisation_acute"
)


