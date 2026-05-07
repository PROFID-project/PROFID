############################################################
# AF-SPECIFIC GBM MODEL VALIDATION & PERFORMANCE
# PROFID – Study 8
############################################################

rm(list = ls())

##############################
# 1 — Packages
##############################
packages <- c("tidyverse","gbm","survival","Hmisc",
              "riskRegression","pec","rmda","caret",
              "fastshap","shapviz","timeROC")

for(pk in packages){
  if(!requireNamespace(pk, quietly=TRUE)) install.packages(pk)
  library(pk, character.only=TRUE)
}

##############################
# 2 — PATHS
##############################

BASEDIR <- "T:/PROFID/Study8"
INDIR   <- file.path(BASEDIR, "Variable Selection & Model Development/Files")
OUTDIR  <- file.path(BASEDIR, "AF model/files")

if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive=TRUE)

##############################
# 3 — LOAD DATA
##############################

df <- read.csv(file.path(INDIR, "vs_data_complete.csv"))

df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.integer(Status),
    event_flag = ifelse(Status == 0,0,1),
    AF_atrial_flutter = ifelse(AF_atrial_flutter=="Yes",1,0)
  ) %>%
  filter(!is.na(Survival_time), !is.na(Status))

############################################################
# IMPORTANT CHANGE
# KEEP ONLY AF PATIENTS
############################################################

df <- df %>% filter(AF_atrial_flutter == 1)

cat("AF patients:", nrow(df), "\n")
cat("Events:", sum(df$event_flag), "\n")

##############################
# 4 — PREDICTORS
##############################

mandatory <- c("LVEF","Age","BMI","Diabetes","eGFR")

biomarker_list <- c(
  "BUN","Cholesterol","CRP","eGFR","Haemoglobin","HbA1c",
  "HDL","IL6","LDL","NTProBNP","Potassium","Sodium",
  "Triglycerides","Troponin_T","TSH"
)

biomarkers <- intersect(biomarker_list, names(df))

ecg_list <- c("HR","PR","QRS","QTc","AV_block","AV_block_II_or_III","LBBB","RBBB")
ecg_vars <- intersect(ecg_list, names(df))

predictors <- unique(c(mandatory, biomarkers, ecg_vars))

cat("Predictors used:", length(predictors), "\n")

##############################
# 5 — TRAIN / VALIDATION SPLIT
##############################

set.seed(2025)

train_idx <- createDataPartition(df$event_flag, p=0.7, list=FALSE)

train <- df[train_idx,]
valid <- df[-train_idx,]

cat("Train:", nrow(train), " Valid:", nrow(valid), "\n")

##############################
# CONVERT CATEGORICAL → NUMERIC
##############################

train_gbm <- train %>%
  mutate(across(where(is.character), ~as.numeric(as.factor(.))),
         across(where(is.factor),   ~as.numeric(.)))

valid_gbm <- valid %>%
  mutate(across(where(is.character), ~as.numeric(as.factor(.))),
         across(where(is.factor),   ~as.numeric(.)))

##############################
# MODEL FORMULA
##############################

model_formula <- as.formula(
  paste0("Surv(Survival_time, event_flag) ~ ", paste(predictors, collapse=" + "))
)

##############################
# RANDOM HYPERPARAMETER SEARCH
##############################

set.seed(2025)

grid <- expand.grid(
  n.trees = c(500,1000,1500),
  interaction.depth = c(3,4,5),
  shrinkage = c(0.01,0.005),
  n.minobsinnode = c(20,30),
  bag.fraction = c(0.6,0.7)
)

grid <- grid[sample(1:nrow(grid),15),]

tune_results <- data.frame()

for(i in 1:nrow(grid)){
  
  g <- grid[i,]
  
  cat("Tuning", i,"/",nrow(grid),"\n")
  
  fit <- gbm(
    formula=model_formula,
    data=train_gbm,
    distribution="coxph",
    n.trees=g$n.trees,
    interaction.depth=g$interaction.depth,
    shrinkage=g$shrinkage,
    n.minobsinnode=g$n.minobsinnode,
    bag.fraction=g$bag.fraction,
    verbose=FALSE
  )
  
  pred_lp <- predict(fit,newdata=valid_gbm,n.trees=g$n.trees,type="link")
  
  c_index <- rcorr.cens(
    -pred_lp,
    Surv(valid_gbm$Survival_time, valid_gbm$event_flag)
  )["C Index"]
  
  tune_results <- rbind(tune_results,cbind(g,C_index=c_index))
}

write.csv(tune_results,
          file.path(OUTDIR,"AF_GBM_Tuning_Results.csv"),
          row.names=FALSE)

best <- tune_results[which.max(tune_results$C_index),]

cat("Best parameters:\n")
print(best)

##############################
# TRAIN FINAL MODEL
##############################

final_gbm <- gbm(
  formula=model_formula,
  data=train_gbm,
  distribution="coxph",
  n.trees=best$n.trees,
  interaction.depth=best$interaction.depth,
  shrinkage=best$shrinkage,
  n.minobsinnode=best$n.minobsinnode,
  bag.fraction=best$bag.fraction,
  verbose=TRUE
)

saveRDS(final_gbm,file.path(OUTDIR,"AF_GBM_Model.rds"))

##############################
# VALIDATION
##############################

valid_lp <- predict(final_gbm,newdata=valid_gbm,
                    n.trees=best$n.trees,type="link")

c_index <- rcorr.cens(
  -valid_lp,
  Surv(valid_gbm$Survival_time, valid_gbm$event_flag)
)["C Index"]

cat("Validation C-index:", round(c_index,3),"\n")

##############################
# CALIBRATION SLOPE
##############################

cal_model <- coxph(
  Surv(Survival_time,event_flag) ~ lp,
  data=valid_gbm %>% mutate(lp=valid_lp)
)

cal_slope <- coef(cal_model)[1]

##############################
# IBS
##############################

max_follow <- max(valid_gbm$Survival_time)

times <- seq(30,floor(max_follow-1),by=30)

base_fit <- coxph(model_formula,data=train_gbm)

base_haz <- basehaz(base_fit,centered=FALSE)

risk_mat <- matrix(NA, nrow=nrow(valid_gbm), ncol=length(times))

for(j in seq_along(times)){
  
  t <- times[j]
  
  h0 <- max(base_haz$hazard[base_haz$time <= t], na.rm=TRUE)
  
  if(is.infinite(h0) | is.na(h0)) h0 <- 0
  
  S_t <- exp(-h0 * exp(valid_lp))
  
  risk_mat[,j] <- 1 - S_t
}

risk_mat[risk_mat < 0] <- 0
risk_mat[risk_mat > 1] <- 1

score_obj <- Score(
  object = list(GBM = risk_mat),
  formula = Surv(Survival_time,event_flag) ~ 1,
  data = valid_gbm,
  metrics = "Brier",
  summary = "ibs",
  times = times
)############################################################
# EXTRACT IBS
############################################################

ibs_vec <- score_obj$Brier$score$IBS

# If multiple values, take the last one
ibs <- as.numeric(tail(ibs_vec,1))

cat("IBS:", round(ibs,4), "\n")

############################################################
# TIME-DEPENDENT AUC (30–150 MONTHS)
############################################################

library(timeROC)

# Timepoints in MONTHS (since Survival_time is in months)
times <- c(30, 60, 90, 120)

auc_obj <- timeROC(
  T = valid_gbm$Survival_time,
  delta = valid_gbm$event_flag,
  marker = valid_lp,
  cause = 1,
  times = times,
  iid = FALSE
)

auc_table <- data.frame(
  Months = times,
  AUC = round(auc_obj$AUC,3)
)

cat("\n--- Time-dependent AUC ---\n")
print(auc_table)

write.csv(
  auc_table,
  file.path(OUTDIR,"AF_TimeDependent_AUC.csv"),
  row.names = FALSE
)
############################################################
# PLOT TIME-DEPENDENT AUC
############################################################

png(file.path(OUTDIR,"AF_TimeDependent_AUC_Plot.png"),
    width = 900, height = 700)

ggplot(auc_table, aes(x = Months, y = AUC)) +
  geom_point(size = 4) +
  geom_line(linewidth = 1) +
  scale_x_continuous(breaks = auc_table$Months) +
  scale_y_continuous(limits = c(0.5,1)) +
  labs(
    title = "Time-dependent AUC for AF-specific GBM model",
    x = "Time (Months)",
    y = "Area Under the Curve (AUC)"
  ) +
  theme_classic()

dev.off()


##############################
# SAVE PERFORMANCE
##############################

perf <- data.frame(
  C_index=c_index,
  Calibration_Slope=cal_slope,
  IBS=ibs
)

write.csv(perf,
          file.path(OUTDIR,"AF_Model_Performance.csv"),
          row.names=FALSE)

############################################################
# SHAP VALUES
############################################################

df_shap <- train_gbm[,predictors]

pred_fun <- function(model,newdata){
  predict(model,newdata,
          n.trees=best$n.trees,
          type="link")
}

set.seed(2025)

df_sub <- df_shap[sample(1:nrow(df_shap),min(2000,nrow(df_shap))),]

shap_vals <- fastshap::explain(
  final_gbm,
  X=df_sub,
  pred_wrapper=pred_fun,
  nsim=100
)

sv <- shapviz(shap_vals,X=df_sub)

png(file.path(OUTDIR,"AF_SHAP_Summary.png"),width=1000,height=800)

sv_importance(sv)

dev.off()

imp <- data.frame(
  Variable=colnames(shap_vals),
  MeanAbsSHAP=apply(abs(shap_vals),2,mean)
) %>%
  arrange(desc(MeanAbsSHAP))

############################################################
# FOREST PLOT OF TOP PREDICTORS
############################################################

library(ggplot2)

# Select top predictors
top_vars <- imp %>%
  slice(1:10)

png(file.path(OUTDIR,"AF_SHAP_ForestPlot.png"),
    width = 900, height = 700)

ggplot(top_vars,
       aes(x = reorder(Variable, MeanAbsSHAP),
           y = MeanAbsSHAP)) +
  geom_point(size = 4) +
  coord_flip() +
  labs(
    title = "Top Predictors of SCD Risk (AF-specific GBM model)",
    x = "Predictor",
    y = "Mean |SHAP value|"
  ) +
  theme_classic()

dev.off()

write.csv(imp,
          file.path(OUTDIR,"AF_SHAP_Importance.csv"),
          row.names=FALSE)

cat("\nAF MODEL COMPLETE\n")
cat("Results saved to:",OUTDIR,"\n")
