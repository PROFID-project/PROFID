

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


############################################################
# FIGURE 6 — AF-ONLY RISK GROUP DISTRIBUTION (FINAL CLEAN)
############################################################

# ----------------------------------------------------------
# LOAD PACKAGES
# ----------------------------------------------------------
library(dplyr)
library(ggplot2)
library(gbm)

# ----------------------------------------------------------
# PATHS
# ----------------------------------------------------------
BASE <- "T:/PROFID/Study8"
DATADIR <- file.path(BASE, "Variable Selection & Model Development/Files")
MODELDIR <- file.path(BASE, "Model Validation and Performance/Files")
OUTDIR <- file.path(BASE, "Risk Score Development/files_2")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# ----------------------------------------------------------
# LOAD DATA + MODEL
# ----------------------------------------------------------
df <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

final_gbm <- readRDS(file.path(MODELDIR, "Final_GBM_Model.rds"))
tune_res <- read.csv(file.path(MODELDIR, "GBM_Tuning_Results_Random.csv"))

best_row <- tune_res[which.max(tune_res$C_index),]

# ----------------------------------------------------------
# PREP DATA
# ----------------------------------------------------------
df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.numeric(Status)
  )

# ----------------------------------------------------------
# COMPUTE LINEAR PREDICTOR
# ----------------------------------------------------------
lp <- predict(
  final_gbm,
  newdata = df,
  n.trees = best_row$n.trees,
  type = "link"
)

df$LP <- lp

# ----------------------------------------------------------
# BUILD RISK SCORE
# ----------------------------------------------------------
scale_factor <- 5
df$RiskScore <- round((df$LP - min(df$LP)) * scale_factor)

# ----------------------------------------------------------
# CREATE RISK GROUPS (TERTILES)
# ----------------------------------------------------------
q1 <- quantile(df$RiskScore, 0.33)
q2 <- quantile(df$RiskScore, 0.66)

df <- df %>%
  mutate(
    RiskGroup = case_when(
      RiskScore <= q1 ~ "Low",
      RiskScore <= q2 ~ "Intermediate",
      TRUE ~ "High"
    )
  )

# ----------------------------------------------------------
# FILTER AF PATIENTS ONLY
# ----------------------------------------------------------
df_AF <- df %>%
  filter(AF_atrial_flutter == "Yes")

cat("AF-only sample size:", nrow(df_AF), "\n")

# ----------------------------------------------------------
# PLOT (WITH COUNTS ON TOP — FIXED)
# ----------------------------------------------------------
p_af <- ggplot(df_AF, aes(x = RiskGroup, fill = Group)) +
  geom_bar(position = position_dodge(width = 0.9)) +
  
  geom_text(
    stat = "count",
    aes(label = after_stat(count)),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 3.5
  ) +
  
  labs(
    title = "Distribution of ICD / EF Category Across Risk Groups (AF Patients Only)",
    x = "Risk Group",
    y = "Number of Patients",
    fill = "Category"
  ) +
  
  scale_fill_manual(
    values = c(
      "ICD" = "#E64B35",
      "NonICD_preserved" = "#4DBBD5",
      "NonICD_reduced" = "#00A087"
    ),
    labels = c(
      "ICD",
      "NonICD_preserved" = "Non-ICD Preserved EF",
      "NonICD_reduced" = "Non-ICD Reduced EF"
    )
  ) +
  
  theme_minimal() +
  theme(
    plot.title = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 11)
  )

# ----------------------------------------------------------
# SAVE FIGURE
# ----------------------------------------------------------
ggsave(
  file.path(OUTDIR, "Figure6_AF_only_with_counts.png"),
  p_af,
  width = 8,
  height = 6,
  dpi = 300
)

# ----------------------------------------------------------
# SAVE DATA (REPRODUCIBILITY)
# ----------------------------------------------------------
write.csv(df_AF,
          file.path(OUTDIR, "Figure6_AF_only_data.csv"),
          row.names = FALSE)

cat("\n✔ DONE: Figure 6 (AF-only) saved in files_2\n")


############################################################
# SCD RATE BY AF STATUS (WITH RISK GROUP CREATION)
############################################################

library(dplyr)
library(ggplot2)
library(gbm)
library(scales)

# ----------------------------------------------------------
# PATHS
# ----------------------------------------------------------
BASE <- "T:/PROFID/Study8"
DATADIR <- file.path(BASE, "Variable Selection & Model Development/Files")
MODELDIR <- file.path(BASE, "Model Validation and Performance/Files")
OUTDIR <- file.path(BASE, "Risk Score Development/files_2")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# ----------------------------------------------------------
# LOAD DATA + MODEL
# ----------------------------------------------------------
df <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

final_gbm <- readRDS(file.path(MODELDIR, "Final_GBM_Model.rds"))
tune_res <- read.csv(file.path(MODELDIR, "GBM_Tuning_Results_Random.csv"))

best_row <- tune_res[which.max(tune_res$C_index),]

# ----------------------------------------------------------
# PREP DATA
# ----------------------------------------------------------
df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.numeric(Status),
    AF_group = ifelse(AF_atrial_flutter == "Yes", "AF Present", "AF Absent"),
    SCD_event = ifelse(Status == 1, 1, 0)
  )

# ----------------------------------------------------------
# COMPUTE RISK SCORE (SAME AS BEFORE)
# ----------------------------------------------------------
lp <- predict(
  final_gbm,
  newdata = df,
  n.trees = best_row$n.trees,
  type = "link"
)

df$RiskScore <- round((lp - min(lp)) * 5)

# ----------------------------------------------------------
# CREATE RISK GROUPS (LOW / INTERMEDIATE / HIGH)
# ----------------------------------------------------------
q1 <- quantile(df$RiskScore, 0.33)
q2 <- quantile(df$RiskScore, 0.66)

df <- df %>%
  mutate(
    RiskGroup = case_when(
      RiskScore <= q1 ~ "Low",
      RiskScore <= q2 ~ "Intermediate",
      TRUE ~ "High"
    )
  )

# ----------------------------------------------------------
# CALCULATE SCD RATE
# ----------------------------------------------------------
summary_df <- df %>%
  group_by(RiskGroup, AF_group) %>%
  summarise(
    N = n(),
    SCD_events = sum(SCD_event),
    SCD_rate = SCD_events / N,
    .groups = "drop"
  )

# ----------------------------------------------------------
# PLOT (SIDE-BY-SIDE, PROPORTIONAL)
# ----------------------------------------------------------
p <- ggplot(summary_df, aes(x = RiskGroup, y = SCD_rate, fill = AF_group)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
  
  geom_text(
    aes(label = paste0(round(SCD_rate * 100, 1), "%")),
    position = position_dodge(width = 0.9),
    vjust = -0.3,
    size = 3.5
  ) +
  
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  
  labs(
    title = "SCD Rate With vs Without AF Across Risk Groups",
    x = "Risk Group",
    y = "SCD Rate (%)",
    fill = "AF Status"
  ) +
  
  scale_fill_manual(
    values = c("AF Absent" = "#2C7FB8", "AF Present" = "#D7301F")
  ) +
  
  theme_minimal()

# ----------------------------------------------------------
# SAVE
# ----------------------------------------------------------
ggsave(
  file.path(OUTDIR, "Figure_SCD_rate_AF.png"),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

cat("✔ DONE — Risk groups + proportional figure created\n")

############################################################
# KAPLAN–MEIER (90-DAY FOLLOW-UP)
# PROFID – Study 8
############################################################

# ----------------------------------------------------------
# LOAD PACKAGES
# ----------------------------------------------------------
library(dplyr)
library(survival)
library(survminer)
library(gbm)
library(scales)

# ----------------------------------------------------------
# PATHS
# ----------------------------------------------------------
BASE <- "T:/PROFID/Study8"
DATADIR <- file.path(BASE, "Variable Selection & Model Development/Files")
MODELDIR <- file.path(BASE, "Model Validation and Performance/Files")
OUTDIR <- file.path(BASE, "Risk Score Development/files_2")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# ----------------------------------------------------------
# LOAD DATA + MODEL
# ----------------------------------------------------------
df <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

final_gbm <- readRDS(file.path(MODELDIR, "Final_GBM_Model.rds"))
tune_res <- read.csv(file.path(MODELDIR, "GBM_Tuning_Results_Random.csv"))

best_row <- tune_res[which.max(tune_res$C_index),]

# ----------------------------------------------------------
# PREP DATA
# ----------------------------------------------------------
df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.numeric(Status),
    event_SCD = ifelse(Status == 1, 1, 0)
  )

# ----------------------------------------------------------
# COMPUTE RISK SCORE (GBM)
# ----------------------------------------------------------
lp <- predict(
  final_gbm,
  newdata = df,
  n.trees = best_row$n.trees,
  type = "link"
)

df$RiskScore <- round((lp - min(lp)) * 5)

# ----------------------------------------------------------
# CREATE RISK GROUPS (TERTILES)
# ----------------------------------------------------------
q1 <- quantile(df$RiskScore, 0.33)
q2 <- quantile(df$RiskScore, 0.66)

df <- df %>%
  mutate(
    RiskGroup = case_when(
      RiskScore <= q1 ~ "Low",
      RiskScore <= q2 ~ "Intermediate",
      TRUE ~ "High"
    )
  )

# ----------------------------------------------------------
# TRUNCATE FOLLOW-UP TO 90 DAYS
# ----------------------------------------------------------
df <- df %>%
  mutate(
    Survival_time_90 = pmin(Survival_time, 90),
    Status_90 = ifelse(Survival_time <= 90, Status, 0),
    event_SCD_90 = ifelse(Status_90 == 1, 1, 0)
  )

# ----------------------------------------------------------
# CREATE KM MODEL
# ----------------------------------------------------------
km_fit <- survfit(
  Surv(Survival_time_90, event_SCD_90) ~ RiskGroup,
  data = df
)

# ----------------------------------------------------------
# 🔥 ADD THIS NEW PLOT CODE HERE
# ----------------------------------------------------------
library(ggplot2)

km_df <- data.frame(
  time = km_fit$time,
  surv = km_fit$surv,
  lower = km_fit$lower,
  upper = km_fit$upper,
  strata = rep(names(km_fit$strata), km_fit$strata)
)

km_df$RiskGroup <- gsub("RiskGroup=", "", km_df$strata)

p <- ggplot(km_df, aes(x = time, y = surv, color = RiskGroup, fill = RiskGroup)) +
  geom_ribbon(aes(ymin = lower, ymax = upper),
              alpha = 0.15,
              color = NA) +
  geom_step(size = 0.8) +
  labs(
    title = "Kaplan–Meier by Risk Group (90-day follow-up)",
    x = "Time (days)",
    y = "Survival probability"
  ) +
  scale_color_manual(values = c(
    "Low" = "#4DBBD5",
    "Intermediate" = "#00A087",
    "High" = "#E64B35"
  )) +
  scale_fill_manual(values = c(
    "Low" = "#4DBBD5",
    "Intermediate" = "#00A087",
    "High" = "#E64B35"
  )) +
  scale_x_continuous(
    breaks = c(0, 10,20,30,40,50,60,70,80,90),
    limits = c(0, 90)
  ) +
  theme_minimal()

# ----------------------------------------------------------
# SAVE
# ----------------------------------------------------------
ggsave(
  file.path(OUTDIR, "KM_RiskGroups_90Days_CLEAN.png"),
  p,
  width = 8,
  height = 6,
  dpi = 300
)

