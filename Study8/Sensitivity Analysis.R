############################################################
# Sensitivity Analysis 1: Complete-case evaluation
############################################################
install.packages(c("survival","riskRegression","gbm","dplyr","broom","MASS","MatchIt","tableone","mice","car","survminer","cobalt"))


library(tidyverse)
library(survival)
library(Hmisc)
library(gbm)

# ----------------------------------------------------------
# PATHS
# ----------------------------------------------------------
BASE <- "T:/PROFID/Study8"
MODELDIR <- file.path(BASE, "Model Validation and Performance/Files")
DATADIR  <- file.path(BASE, "Variable Selection & Model Development/Files")
OUTDIR   <- file.path(BASE, "Sensitivity Analysis/Files")

if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# ----------------------------------------------------------
# LOAD DATA + MODEL
# ----------------------------------------------------------
df_cc <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

final_gbm <- readRDS(file.path(MODELDIR, "Final_GBM_Model.rds"))
tune_res  <- read.csv(file.path(MODELDIR, "GBM_Tuning_Results_Random.csv"))
best_row  <- tune_res[which.max(tune_res$C_index),]

best_params <- list(
  n.trees = best_row$n.trees,
  interaction.depth = best_row$interaction.depth,
  shrinkage = best_row$shrinkage,
  n.minobsinnode = best_row$n.minobsinnode,
  bag.fraction = best_row$bag.fraction
)

# ----------------------------------------------------------
# MODEL PERFORMANCE ON COMPLETE CASES
# ----------------------------------------------------------
cat("\nRunning complete-case sensitivity analysis...\n")

# Predict linear predictor on complete-case dataset
lp <- predict(final_gbm,
              newdata = df_cc,
              n.trees = best_params$n.trees,
              type = "link")

# C-index
cin <- rcorr.cens(-lp, Surv(df_cc$Survival_time, df_cc$Status))["C Index"]

# Calibration slope
cal_fit <- coxph(Surv(Survival_time, Status) ~ lp, data=df_cc)
cal_slope <- as.numeric(cal_fit$coef[1])

# IBS (simple)
risks <- 1 - exp(-exp(lp))
brier <- mean((df_cc$Status - risks)^2)
ibs <- sqrt(brier)

# Store results
out <- data.frame(
  Dataset = "Complete-case only",
  N = nrow(df_cc),
  Events = sum(df_cc$Status),
  C_index = round(cin,4),
  Calibration_Slope = round(cal_slope,4),
  IBS = round(ibs,4)
)

write.csv(out,
          file.path(OUTDIR, "Sensitivity_CompleteCase.csv"),
          row.names = FALSE)

cat("\n✓ Complete-case sensitivity results saved to:\n",
    file.path(OUTDIR, "Sensitivity_CompleteCase.csv"), "\n")
print(out)


###############################################################################
# SENSITIVITY ANALYSIS — ALL OUTPUTS IN ONE FOLDER
###############################################################################

library(tidyverse)
library(survival)
library(Hmisc)
library(riskRegression)
library(ggplot2)
library(survminer)

###############################################################################
# PATHS
###############################################################################

BASE <- "T:/PROFID/Study8"

DATADIR  <- file.path(BASE, "Variable Selection & Model Development/Files")
MODELDIR <- file.path(BASE, "Model Validation and Performance/Files")
OUTDIR   <- file.path(BASE, "Sensitivity Analysis/Files")   # <-- THIS IS WHERE EVERYTHING GOES

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

###############################################################################
# LOAD MAIN DATASET (your modelling dataset)
###############################################################################

df <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.integer(Status),
    event_flag = ifelse(Status == 1, 1, 0)
  )

cat("Loaded dataset:", nrow(df), "rows\n")

###############################################################################
# 1 — COMPLETE-CASE ANALYSIS (Assess missing data impact)
###############################################################################

cat("\n===== 1. COMPLETE CASE ANALYSIS =====\n")

df_cc <- df %>% filter(complete.cases(.))

write.csv(df_cc, file.path(OUTDIR, "CompleteCase_Dataset.csv"), row.names = FALSE)

cat("Complete-case rows:", nrow(df_cc), "\n")

cox_cc <- coxph(Surv(Survival_time, event_flag) ~ ., data=df_cc)

write.csv(
  broom::tidy(cox_cc, exponentiate=TRUE, conf.int=TRUE),
  file.path(OUTDIR, "CompleteCase_Cox_Model.csv"),
  row.names = FALSE
)

lp_cc <- predict(cox_cc, type="lp")
c_cc <- rcorr.cens(-lp_cc, Surv(df_cc$Survival_time, df_cc$event_flag))["C Index"]

write.csv(
  data.frame(C_index=c_cc),
  file.path(OUTDIR, "CompleteCase_CIndex.csv"),
  row.names=FALSE
)

cat("✔ Complete-case C-index:", round(c_cc,3), "\n")

# Plot LP distribution
png(file.path(OUTDIR, "CompleteCase_LP_Distribution.png"), width=900, height=600)
ggplot(data.frame(lp=lp_cc), aes(lp)) +
  geom_density(fill="blue", alpha=.3) +
  theme_minimal() +
  labs(title="Complete Case — Linear Predictor Distribution")
dev.off()


###############################################################################
# 2 — COMPETING-RISK ALTERNATIVE (Cause-Specific Cox)
###############################################################################

cat("\n===== 2. COMPETING RISK ANALYSIS =====\n")

df <- df %>%
  mutate(
    event_SCD = ifelse(Status == 1, 1, 0),
    event_COMP = ifelse(Status == 2, 1, 0)
  )

# Cause-specific Cox for SCD
cs_model <- coxph(Surv(Survival_time, event_SCD) ~ ., data=df)

write.csv(
  broom::tidy(cs_model, exponentiate=TRUE, conf.int=TRUE),
  file.path(OUTDIR, "CauseSpecific_Cox_Model.csv"),
  row.names = FALSE
)

lp_cs <- predict(cs_model)
c_cs <- rcorr.cens(-lp_cs, Surv(df$Survival_time, df$event_SCD))["C Index"]

write.csv(
  data.frame(C_index=c_cs),
  file.path(OUTDIR, "CauseSpecific_CIndex.csv"),
  row.names=FALSE
)

cat("✔ Cause-specific C-index:", round(c_cs,3), "\n")

# CIF Plots
fit_scd <- survfit(Surv(Survival_time, event_SCD) ~ 1, data=df)
p1 <- ggsurvplot(fit_scd, conf.int=TRUE, title="CIF — Sudden Cardiac Death")
ggsave(file.path(OUTDIR, "CIF_SCD.png"), p1$plot, width=7, height=6)

fit_comp <- survfit(Surv(Survival_time, event_COMP) ~ 1, data=df)
p2 <- ggsurvplot(fit_comp, conf.int=TRUE, title="CIF — Competing Event")
ggsave(file.path(OUTDIR, "CIF_COMP.png"), p2$plot, width=7, height=6)

cat("\n===== 3. TEMPORAL ROBUSTNESS =====\n")

# Split datasets
df_old <- df %>% filter(Baseline_type == "MI")
df_new <- df %>% filter(Baseline_type == "MI40d")

# Function to drop all single-level variables
remove_single_level <- function(data) {
  keep <- names(data)[sapply(data, function(x) length(unique(na.omit(x))) > 1)]
  data %>% select(all_of(keep))
}

df_old_clean <- remove_single_level(df_old)
df_new_clean <- remove_single_level(df_new)

# Fit Cox models
cox_old <- coxph(Surv(Survival_time, event_flag) ~ ., data = df_old_clean)
cox_new <- coxph(Surv(Survival_time, event_flag) ~ ., data = df_new_clean)

# Save
write.csv(broom::tidy(cox_old, exponentiate=TRUE, conf.int=TRUE),
          file.path(OUTDIR,"Temporal_MI_Model.csv"), row.names=FALSE)

write.csv(broom::tidy(cox_new, exponentiate=TRUE, conf.int=TRUE),
          file.path(OUTDIR,"Temporal_MI40d_Model.csv"), row.names=FALSE)

cat("✓ Temporal robustness models saved.\n")


###############################################################################
# ==== 4. AF VARIABLE DEFINITIONS SENSITIVITY ====
###############################################################################

cat("\n==== 4. AF VARIABLE DEFINITIONS SENSITIVITY =====\n")

#---------------------------------------
# 1. Broad AF definition (AF alone)
#---------------------------------------
df$AF_broad <- ifelse(df$AF_atrial_flutter == "Yes", 1, 0)

#---------------------------------------
# 2. Strict AF definition (AF + anticoagulant)
#---------------------------------------
df$AF_strict <- ifelse(
  df$AF_atrial_flutter == "Yes" & df$Anti_coagulant == "Yes",
  1, 0
)

# Check levels
cat("\nBroad AF levels:\n")
print(table(df$AF_broad))

cat("\nStrict AF levels:\n")
print(table(df$AF_strict))


###############################################################################
# Function to compute model performance for AF definitions
###############################################################################
run_af_sensitivity <- function(af_var){
  
  cat("\nRunning AF sensitivity for:", af_var, "\n")
  
  df_temp <- df
  
  # Predict linear predictor
  lp <- predict(final_gbm, newdata=df_temp, n.trees=best_params$n.trees)
  
  # C-index
  C_idx <- rcorr.cens(lp, Surv(df_temp$Survival_time, df_temp$Status))["C Index"]
  
  # Calibration slope
  cal_model <- coxph(Surv(Survival_time, Status) ~ lp, data=df_temp)
  cal_slope <- coef(cal_model)[1]
  
  # IBS (simple Brier score approximation)
  risks <- 1 - exp(-exp(lp))
  ibs_val <- mean((df_temp$Status - risks)^2)
  
  out <- data.frame(
    AF_Definition = af_var,
    N = nrow(df_temp),
    Events = sum(df_temp$Status == 1),
    C_index = round(C_idx,4),
    Calibration_Slope = round(cal_slope,4),
    IBS = round(ibs_val,4)
  )
  
  # LP distribution plot
  png(file.path(OUTDIR, paste0("AF_", af_var, "_LP_Distribution.png")),
      width=900, height=600)
  ggplot(data.frame(lp=lp), aes(lp)) +
    geom_density(fill="blue", alpha=0.3) +
    theme_minimal() +
    labs(title=paste("LP Distribution —", af_var))
  dev.off()
  
  return(out)
}


###############################################################################
# Run AF sensitivity analyses
###############################################################################
AF_defs <- c("AF_broad", "AF_strict")

af_results <- lapply(AF_defs, run_af_sensitivity)
af_results <- bind_rows(af_results)

write.csv(af_results,
          file.path(OUTDIR, "Sensitivity_AF_AllResults.csv"),
          row.names = FALSE)

cat("\n✓ AF sensitivity results saved to:",
    file.path(OUTDIR, "Sensitivity_AF_AllResults.csv"), "\n")

print(af_results)

###############################################################################
# 5. VARYING TIME WINDOWS — ACUTE VS CHRONIC RISK FACTORS
###############################################################################

cat("\n===== 5. ACUTE vs CHRONIC RISK FACTORS SENSITIVITY =====\n")

# Load complete-case dataset (already used for other sensitivity analyses)
df <- df_cc # using your imported vs_data_complete.csv
df$event_flag <- ifelse(df$Status == 1, 1, 0)

# Identify acute variables (present in your dataset)
acute_vars <- intersect(
  c("PCI_acute", "Thrombolysis_acute", "Revascularisation_acute"),
  names(df)
)

cat("Acute variables detected:", paste(acute_vars, collapse=", "), "\n")

# Identify all predictors except outcome
all_predictors <- setdiff(names(df), c("Survival_time", "Status", "event_flag"))

# Chronic predictors = all predictors EXCEPT acute
chronic_vars <- setdiff(all_predictors, acute_vars)

# Safety check
cat("Chronic variables count:", length(chronic_vars), "\n")
cat("Total predictors:", length(all_predictors), "\n")

###############################
# 5A — Model with ALL predictors (Acute + Chronic)
###############################
form_all <- as.formula(
  paste("Surv(Survival_time, event_flag) ~", paste(all_predictors, collapse=" + "))
)

cox_all <- coxph(form_all, data = df)

# Compute C-index
lp_all <- predict(cox_all, type="lp")
C_all <- rcorr.cens(-lp_all, Surv(df$Survival_time, df$event_flag))["C Index"]

# Calibration
cal_all <- coxph(Surv(Survival_time, event_flag) ~ lp_all, data=df)$coef[1]

# IBS simple estimate
risk_all <- 1 - exp(-exp(lp_all))
ibs_all <- sqrt(mean((df$event_flag - risk_all)^2))

###############################
# 5B — Model with Chronic-only predictors
###############################
form_chronic <- as.formula(
  paste("Surv(Survival_time, event_flag) ~", paste(chronic_vars, collapse=" + "))
)

cox_chronic <- coxph(form_chronic, data = df)

lp_chronic <- predict(cox_chronic, type="lp")
C_chronic <- rcorr.cens(-lp_chronic, Surv(df$Survival_time, df$event_flag))["C Index"]
cal_chronic <- coxph(Surv(Survival_time, event_flag) ~ lp_chronic, data=df)$coef[1]

risk_chronic <- 1 - exp(-exp(lp_chronic))
ibs_chronic <- sqrt(mean((df$event_flag - risk_chronic)^2))

###############################
# Combine results
###############################
acute_chronic_results <- data.frame(
  Model = c("Full_Acute+Chronic", "Chronic_only"),
  Num_Predictors = c(length(all_predictors), length(chronic_vars)),
  C_index = round(c(C_all, C_chronic), 4),
  Calibration_Slope = round(c(cal_all, cal_chronic), 4),
  IBS = round(c(ibs_all, ibs_chronic), 4)
)

# Save
write.csv(
  acute_chronic_results,
  file.path(OUTDIR, "Sensitivity_AcuteVsChronic.csv"),
  row.names = FALSE
)

cat("\nAcute vs Chronic sensitivity results saved to:\n",
    file.path(OUTDIR, "Sensitivity_AcuteVsChronic.csv"), "\n")

print(acute_chronic_results)

###############################
# Plot LP distributions for comparison
###############################
library(ggplot2)

# LP for full model
png(file.path(OUTDIR, "LP_Full_AcuteChronic.png"), width=900, height=600)
ggplot(data.frame(lp=lp_all), aes(lp)) +
  geom_density(fill="blue", alpha=0.3) +
  theme_minimal() +
  labs(title="LP Distribution — Full Model (Acute + Chronic)")
dev.off()

# LP for chronic-only model
png(file.path(OUTDIR, "LP_ChronicOnly.png"), width=900, height=600)
ggplot(data.frame(lp=lp_chronic), aes(lp)) +
  geom_density(fill="darkgreen", alpha=0.3) +
  theme_minimal() +
  labs(title="LP Distribution — Chronic-Only Model")
dev.off()

