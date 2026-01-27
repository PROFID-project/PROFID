

############################################################
# RISK SCORE DEVELOPMENT SCRIPT
# PROFID – Study 8
############################################################
install.packages(c("timeROC","gbm","survAUC","survival","dplyr","broom","MASS","MatchIt","tableone","mice","car","survminer","cobalt"))

library(dplyr)
library(survival)
library(rms)
library(ggplot2)
library(riskRegression)
library(survAUC)
library("gbm")

############################################################
# PATHS
############################################################

BASE <- "T:/PROFID/Study8"
DATADIR <- file.path(BASE, "Variable Selection & Model Development/Files")
MODELDIR <- file.path(BASE, "Model Validation and Performance/Files")
OUTDIR <- file.path(BASE, "Risk Score Development/Files")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

############################################################
# LOAD DATA + FINAL MODEL
############################################################

df <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

final_gbm <- readRDS(file.path(MODELDIR, "Final_GBM_Model.rds"))
tune_res <- read.csv(file.path(MODELDIR, "GBM_Tuning_Results_Random.csv"))
best_row <- tune_res[which.max(tune_res$C_index),]

best_params <- list(
  n.trees = best_row$n.trees,
  interaction.depth = best_row$interaction.depth,
  shrinkage = best_row$shrinkage,
  n.minobsinnode = best_row$n.minobsinnode,
  bag.fraction = best_row$bag.fraction
)

df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.numeric(Status),
    event_flag = ifelse(Status == 1, 1, 0)
  )

############################################################
# STEP 1: COMPUTE LINEAR PREDICTOR
############################################################

lp <- predict(
  final_gbm,
  newdata = df,
  n.trees = best_params$n.trees,
  type = "link"
)

df$LP <- lp

############################################################
# STEP 2: BUILD INTEGER RISK SCORE
############################################################

# Convert LP to integer points
scale_factor <- 5 # adjust if needed
df$RiskScore <- round((df$LP - min(df$LP)) * scale_factor)

# Save distribution
write.csv(df[,c("LP","RiskScore")],
          file.path(OUTDIR,"RiskScore_Distribution.csv"),
          row.names = FALSE)

############################################################
# STEP 3: CREATE RISK CATEGORIES
############################################################

df <- df %>%
  mutate(
    RiskGroup = case_when(
      RiskScore <= quantile(RiskScore,0.33) ~ "Low",
      RiskScore <= quantile(RiskScore,0.66) ~ "Intermediate",
      TRUE ~ "High"
    )
  )

write.csv(df, file.path(OUTDIR,"RiskScore_with_Groups.csv"), row.names = FALSE)

############################################################
# STEP 4: VALIDATION – TIME-DEPENDENT AUC
############################################################

times <- c(30, 90, 180, 365, 730, 1095)

auc <- AUC.uno(
  Surv(df$Survival_time, df$event_flag),
  Surv(df$Survival_time, df$event_flag),
  df$RiskScore,
  times = times
)

auc_df <- data.frame(
  Time = times,
  AUC = auc$auc
)

write.csv(auc_df, file.path(OUTDIR,"AUC_TimeDependent.csv"), row.names = FALSE)

########################################################
# STEP X: TIME-DEPENDENT AUC (fixed version)
########################################################

# Choose clinically meaningful timepoints
times <- c(30, 60, 90, 120, 150)

# Compute time-dependent AUC (iid=FALSE avoids memory issues)
auc_obj <- timeROC(
  T = df$Survival_time,
  delta = df$event_flag,
  marker = lp,
  cause = 1,
  times = times,
  iid = FALSE
)

# Build results table
auc_df <- data.frame(
  Time = times,
  AUC = auc_obj$AUC
)

# Print to console
cat("\n--- Time-dependent AUC ---\n")
print(auc_df)

# Save AUC results
write.csv(auc_df,
          file.path(OUTDIR, "AUC_TimeDependent.csv"),
          row.names = FALSE)



############################################################
# STEP 5: CALIBRATION SLOPE
############################################################

f_cal <- coxph(Surv(Survival_time, event_flag) ~ RiskScore, data=df)
write.csv(broom::tidy(f_cal, exponentiate=TRUE, conf.int=TRUE),
          file.path(OUTDIR,"CalibrationSlope.csv"), row.names=FALSE)

############################################################
# STEP 6: BRIER SCORE
############################################################

bs <- Score(
  list("RiskScore" = df$RiskScore),
  formula = Surv(Survival_time, event_flag) ~ 1,
  data = df,
  times = times,
  metrics = c("brier"),
  summary = "ibs"
)

write.csv(bs$Brier$score,
          file.path(OUTDIR,"BrierScore.csv"), row.names = FALSE)

############################################################
# STEP 7: KAPLAN–MEIER CURVES OF RISK GROUPS
############################################################
library(survival)
library(survminer)
km_fit <- survfit(Surv(Survival_time, event_flag) ~ RiskGroup, data=df)

png(file.path(OUTDIR,"KM_RiskGroups.png"), width=900, height=700)
ggsurvplot(km_fit, data=df, risk.table=TRUE,
           title="Kaplan–Meier by Risk Group",
           legend.title="Risk Group")
dev.off()


###############################################
# STEP 8 — NOMOGRAM (Final working version)
###############################################

library(survival)
library(rms)

# Determine max time we can safely model
max_follow <- max(df$Survival_time, na.rm = TRUE)
cat("Using follow-up time:", max_follow, "\n")

# Set datadist (required for rms)
dd <- datadist(df)
options(datadist = "dd")

# Fit Cox model
cox_nom <- cph(
  Surv(Survival_time, event_flag) ~ RiskScore,
  data = df,
  x = TRUE,
  y = TRUE,
  surv = TRUE
)

# Define survival function at max follow-up
surv_fun <- function(lp) {
  s <- survest(cox_nom, linear.predictors = lp, times = max_follow)
  return(s$surv)
}

# Build nomogram object
nom_obj <- nomogram(
  cox_nom,
  fun = list(surv_fun),
  funlabel = paste0("Survival at ", round(max_follow,1), " days")
)

# Save nomogram plot
png(file.path(OUTDIR, "Nomogram.png"), width = 900, height = 700)
plot(nom_obj)
dev.off()

cat("Nomogram created using follow-up time =", max_follow, "\n")

############################################################
# STEP 9: PRINT SUMMARY FOR MANUSCRIPT
############################################################

cat("\n\n================ VALIDATION SUMMARY ================\n")

cat("\n--- Time-dependent AUC ---\n")
print(auc_df)

cat("\nAUC at key times:\n")
for(t in times){
  cat(t,"days: ", round(auc_df$AUC[auc_df$Time==t],3),"\n")
}

cat("\n--- Calibration Slope ---\n")
cat("Slope:", round(coef(f_cal)["RiskScore"],3), "\n")

cat("\n--- Brier Score (first rows) ---\n")
print(head(bs$Brier$score))

ibs <- bs$Brier$IBS
cat("\nIntegrated Brier Score:", round(ibs,4), "\n")

cat("\n--- Risk Score Summary ---\n")
print(summary(df$RiskScore))

cat("\nQuantiles:\n")
print(quantile(df$RiskScore, probs=c(.1,.25,.5,.75,.9)))

cat("\n================ END OF SUMMARY ================\n\n")
