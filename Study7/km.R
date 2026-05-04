###############################################
# Kaplan-Meier curves (ICD and Non-ICD cohorts)
###############################################

# Packages
req <- c("dplyr","tidyverse","data.table","survival","survminer","ggplot2","scales")
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

# Kaplan-Meier curves (ICD cohort)
df_ICD_surv <- df_ICD %>%
  filter(!is.na(Status)) %>%
  mutate(Region = CVD_risk_region) 
km_fit_ICD <- survfit(
  Surv(Survival_time, Status == 1) ~ Region,
  data = df_ICD_surv
)
km_ICD <- ggsurvplot(
  km_fit_ICD,
  data = df_ICD_surv,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#2C7BB6", "#FDAE61", "#D7191C"),
  xlab = "Time (months)",
  ylab = "SCD-free survival probability",
  legend.title = "Region",
  title = "Kaplan-Meier Survival Curves by Cardiovascular Risk Region",
  subtitle = "ICD Cohort",
  ggtheme = theme_minimal()
)
ggsave(filename = "km_plot_ICD.png", plot = km_ICD$plot)

# Kaplan-Meier curves (Non-ICD reduced cohort - ≤35%)
df_NR_surv <- df_NR %>%
  filter(!is.na(Status)) %>%
  mutate(Region = CVD_risk_region) 
km_fit_NR <- survfit(
  Surv(Survival_time, Status == 1) ~ Region,
  data = df_NR_surv
)
km_NR <- ggsurvplot(
  km_fit_NR,
  data = df_NR_surv,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#2C7BB6", "#FFFFBF", "#FDAE61", "#D7191C"),
  xlab = "Time (months)",
  ylab = "SCD-free survival probability",
  legend.title = "Region",
  title = "Kaplan-Meier Survival Curves by Cardiovascular Risk Region",
  subtitle = "Non-ICD Cohort (≤35%)",
  ggtheme = theme_minimal()
)
ggsave(filename = "km_plot_NR.png", plot = km_NR$plot)

# Kaplan-Meier curves (Non-ICD preserved cohort - >35%)
df_NP_surv <- df_NP %>%
  filter(!is.na(Status)) %>%
  mutate(Region = CVD_risk_region) 
km_fit_NP <- survfit(
  Surv(Survival_time, Status == 1) ~ Region,
  data = df_NP_surv
)
km_NP <- ggsurvplot(
  km_fit_NP,
  data = df_NP_surv,
  risk.table = TRUE,
  conf.int = FALSE,
  palette = c("#2C7BB6", "#FFFFBF", "#FDAE61", "#D7191C"),
  xlab = "Time (months)",
  ylab = "SCD-free survival probability",
  legend.title = "Region",
  title = "Kaplan-Meier Survival Curves by Cardiovascular Risk Region",
  subtitle = "Non-ICD Cohort (>35%)",
  ggtheme = theme_minimal()
)
ggsave(filename = "km_plot_NP.png", plot = km_NP$plot)