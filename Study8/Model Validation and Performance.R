############################################################
# MODEL VALIDATION & PERFORMANCE — GBM ONLY
# This script performs:
#  1. Load processed dataset (train/valid already created earlier)
#  2. Hyperparameter tuning for GBM (grid search)
#  3. Fit final tuned GBM
#  4. Internal validation: Harrell C-index + bootstrap CI
#  5. Calibration at a fixed time
#  6. Integrated Brier Score
#  7. Decision Curve Analysis
###############################################################

rm(list = ls())

##############################
# 1 — Packages and Paths
##############################
packages <- c("tidyverse","gbm","survival","Hmisc",
              "riskRegression","pec","rmda","caret")

for(pk in packages){
  if(!requireNamespace(pk, quietly=TRUE)) install.packages(pk)
  library(pk, character.only=TRUE)
}

BASEDIR <- "T:/PROFID/Study8"
INDIR   <- file.path(BASEDIR, "Variable Selection & Model Development/Files")
OUTDIR  <- file.path(BASEDIR, "Model Validation and Performance/Files")


df <- read.csv(file.path(INDIR, "vs_data_complete.csv"))


# preprocess survival variables
df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.integer(Status),
    event_flag = ifelse(Status == 0, 0, 1)
  ) %>%
  filter(!is.na(Survival_time), !is.na(Status))

##############################
# 3 — BIOMARKER + ECG SELECTION
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

cat("Biomarkers used:\n"); print(biomarkers)
cat("ECG variables used:\n"); print(ecg_vars)

predictors <- unique(c(mandatory, biomarkers, ecg_vars))
cat("Total predictors used:", length(predictors), "\n")


##############################
# 4 — TRAIN/VALIDATION SPLIT
##############################
set.seed(2025)
train_idx <- createDataPartition(df$event_flag, p = 0.7, list = FALSE)
train <- df[train_idx, ]
valid <- df[-train_idx, ]

cat("Train size:", nrow(train), "| Valid size:", nrow(valid), "\n")

# convert all character/factor variables → numeric for GBM
train_gbm <- train %>% mutate(across(where(is.character), ~as.numeric(as.factor(.))),
                              across(where(is.factor),   ~as.numeric(.)))
valid_gbm <- valid %>% mutate(across(where(is.character), ~as.numeric(as.factor(.))),
                              across(where(is.factor),   ~as.numeric(.)))

##############################
# 5 — CREATE FORMULA
##############################
model_formula <- as.formula(
  paste0("Surv(Survival_time, event_flag) ~ ", paste(predictors, collapse=" + "))
)

cat("Model Formula:\n")
print(model_formula)


###########################################################################
# 6 — RANDOM HYPERPARAMETER SEARCH (15 combinations)
###########################################################################
set.seed(2025)

full_grid <- expand.grid(
  n.trees = c(500, 1000, 1500),
  interaction.depth = c(3, 4, 5),
  shrinkage = c(0.01, 0.005),
  n.minobsinnode = c(20, 30),
  bag.fraction = c(0.6, 0.7)
)

grid <- full_grid[sample(1:nrow(full_grid), 15), ]
cat("Random search grid size:", nrow(grid), "\n")

tune_results <- data.frame()

for(i in 1:nrow(grid)){
  g <- grid[i,]
  cat("Random search:", i, "/", nrow(grid), "\n")
  
  fit <- gbm(
    formula = model_formula,
    data = train_gbm,
    distribution = "coxph",
    n.trees = g$n.trees,
    interaction.depth = g$interaction.depth,
    shrinkage = g$shrinkage,
    n.minobsinnode = g$n.minobsinnode,
    bag.fraction = g$bag.fraction,
    verbose = FALSE
  )
  
  pred_lp <- predict(fit, newdata = valid_gbm, n.trees = g$n.trees, type="link")
  c_index <- rcorr.cens(-pred_lp, Surv(valid_gbm$Survival_time, valid_gbm$event_flag))["C Index"]
  
  tune_results <- rbind(tune_results, cbind(g, C_index = c_index))
}

write.csv(tune_results, file.path(OUTDIR, "GBM_Tuning_Results_Random.csv"), row.names = FALSE)

best <- tune_results[which.max(tune_results$C_index),]
cat("\nBest hyperparameters found:\n")
print(best)


###########################################################################
# 7 — TRAIN FINAL TUNED MODEL
###########################################################################
start_time <- Sys.time()

final_gbm <- gbm(
  formula = model_formula,
  data = train_gbm,
  distribution = "coxph",
  n.trees = best$n.trees,
  interaction.depth = best$interaction.depth,
  shrinkage = best$shrinkage,
  n.minobsinnode = best$n.minobsinnode,
  bag.fraction = best$bag.fraction,
  keep.data = FALSE,
  verbose = TRUE
)

end_time <- Sys.time()
training_time <- end_time - start_time
cat("\nTraining time:", training_time, "\n")

saveRDS(final_gbm, file.path(OUTDIR, "Final_GBM_Model.rds"))

valid_lp <- predict(final_gbm, newdata = valid_gbm,
                    n.trees = best$n.trees, type="link")

c_val <- rcorr.cens(-valid_lp, Surv(valid_gbm$Survival_time, valid_gbm$event_flag))["C Index"]
cat("Validation C-index:", round(c_val, 3), "\n")


# ----------------------------
# Save performance metrics (C-index)
# ----------------------------
perf_df <- data.frame(
  Model = "GBM (tuned)",
  C_index = round(as.numeric(c_val), 4),
  Date_Run = Sys.time(),
  Train_n = nrow(train_gbm),
  Valid_n = nrow(valid_gbm)
)

perf_file <- file.path(OUTDIR, "GBM_Final_Performance.csv")
write.csv(perf_df, perf_file, row.names = FALSE)

cat("📊 Saved final performance →", perf_file, "\n")

best_params
str(best_params)


###########################################################################
# TRUE 10-FOLD CROSS-VALIDATION FOR SURVIVAL GBM
# Using hyperparameters loaded from tuning file
###########################################################################
###########################################################################
# LOAD BEST PARAMETERS FROM FILE (FULL SET)
###########################################################################

tune_file <- file.path(OUTDIR, "GBM_Tuning_Results_Random.csv")

tune_results <- read.csv(tune_file)

# Ensure numeric formatting
tune_results <- tune_results %>%
  mutate(
    n.trees = as.integer(n.trees),
    interaction.depth = as.integer(interaction.depth),
    shrinkage = as.numeric(shrinkage),
    n.minobsinnode = as.integer(n.minobsinnode),
    bag.fraction = as.numeric(bag.fraction),
    C_index = as.numeric(C_index)
  )

# Pick the best row (highest C-index)
best_row <- tune_results[which.max(tune_results$C_index),]

# Extract ALL hyperparameters into a list
best_params <- list(
  n.trees = best_row$n.trees,
  interaction.depth = best_row$interaction.depth,
  shrinkage = best_row$shrinkage,
  n.minobsinnode = best_row$n.minobsinnode,
  bag.fraction = best_row$bag.fraction
)

cat("\n🏆 Best hyperparameters loaded from file:\n")
print(best_params)

###########################################################################
# TRUE 10-FOLD CROSS-VALIDATION FOR SURVIVAL GBM
###########################################################################

library(survival)
library(gbm)
library(Hmisc)
library(caret)

set.seed(2025)

df_cv <- rbind(train_gbm, valid_gbm)

# 1. Create 10 folds
folds <- createFolds(df_cv$event_flag, k = 10, list = TRUE)

# 2. Prepare vector for out-of-fold predictions
lp_cv <- rep(NA, nrow(df_cv))

# 3. Perform 10-fold CV
for (i in 1:10) {
  
  cat("Running fold", i, "of 10...\n")
  
  test_idx  <- folds[[i]]
  train_idx <- setdiff(seq_len(nrow(df_cv)), test_idx)
  
  df_train <- df_cv[train_idx, ]
  df_test  <- df_cv[test_idx, ]
  
  fit <- gbm(
    formula = model_formula,
    data = df_train,
    distribution = "coxph",
    n.trees = best_params$n.trees,
    interaction.depth = best_params$interaction.depth,
    shrinkage = best_params$shrinkage,
    n.minobsinnode = best_params$n.minobsinnode,
    bag.fraction = best_params$bag.fraction,
    keep.data = FALSE,
    verbose = TRUE
  )
  
  lp_fold <- predict(
    fit,
    newdata = df_test,
    n.trees = best_params$n.trees,
    type = "link"
  )
  
  lp_cv[test_idx] <- lp_fold
}

# 4. Compute CV C-index
cv_cindex <- rcorr.cens(
  -lp_cv,
  Surv(df_cv$Survival_time, df_cv$event_flag)
)["C Index"]

cat("\n======================================\n")
cat("⭐ TRUE 10-FOLD CV C-index:", round(cv_cindex, 3), "\n")
cat("======================================\n")

# 5. Save results
write.csv(
  data.frame(
    Model = "GBM (10-fold CV)",
    C_index = round(as.numeric(cv_cindex), 4),
    n_folds = 10,
    Date_Run = Sys.time()
  ),
  file.path(OUTDIR, "GBM_10Fold_CV_Performance.csv"),
  row.names = FALSE
)

cat("📄 Saved → GBM_10Fold_CV_Performance.csv\n")


###############################################################################
# 5.3 PERFORMANCE METRICS
###############################################################################

library(survival)
library(riskRegression)
library(pec)
library(rmda)

cat("\n========================================\n")
cat("📌 5.3 PERFORMANCE METRICS\n")
cat("========================================\n\n")

###############################################################################
# 5.3.1 DISCRIMINATION — HARRELL'S C-INDEX
###############################################################################

final_model_path <- file.path(OUTDIR, "Final_GBM_Model.rds")
final_gbm <- readRDS(final_model_path)

tune_file <- file.path(OUTDIR, "GBM_Tuning_Results_Random.csv")
tune_results <- read.csv(tune_file)

best_row <- tune_results[which.max(tune_results$C_index),]

best_params <- list(
  n.trees = best_row$n.trees,
  interaction.depth = best_row$interaction.depth,
  shrinkage = best_row$shrinkage,
  n.minobsinnode = best_row$n.minobsinnode,
  bag.fraction = best_row$bag.fraction
)



cat("\n---- 5.3.1 Harrell's C-index ----\n")

train_lp <- predict(final_gbm, newdata=train_gbm,
                    n.trees=best_params$n.trees, type="link")

valid_lp <- predict(final_gbm, newdata=valid_gbm,
                    n.trees=best_params$n.trees, type="link")

c_train <- rcorr.cens(-train_lp, Surv(train_gbm$Survival_time, train_gbm$event_flag))["C Index"]
c_valid <- rcorr.cens(-valid_lp, Surv(valid_gbm$Survival_time, valid_gbm$event_flag))["C Index"]

cat("Train C-index:", round(c_train, 3), "\n")
cat("Validation C-index:", round(c_valid, 3), "\n")

write.csv(
  data.frame(
    Train_C = c_train,
    Valid_C = c_valid
  ),
  file.path(OUTDIR, "GBM_Discrimination_CIndex.csv"),
  row.names = FALSE
)

###############################################################################
# 5.3.2 CALIBRATION — Integrated Brier Score (IBS) + Calibration Slope
###############################################################################

cat("\n---- 5.3.2 CALIBRATION (IBS + SLOPE) ----\n")

library(survival)
library(riskRegression)
library(pec)

# 1. Evaluation times
max_follow <- max(valid_gbm$Survival_time, na.rm = TRUE)
eval_times <- seq(30, floor(max_follow - 1), by = 30)
if(length(eval_times) == 0) eval_times <- floor(max_follow/2)

cat("Evaluation times:", paste(eval_times, collapse=", "), "\n")

# 2. Baseline Cox for cumulative hazard
base_fit <- coxph(model_formula, data=train_gbm, x=TRUE)
base_haz <- basehaz(base_fit, centered=FALSE)

# 3. Linear predictor from final GBM
valid_lp <- predict(
  final_gbm,
  newdata = valid_gbm,
  n.trees = best_params$n.trees,
  type = "link"
)

# 4. Convert LP → risk(t)
risk_mat <- matrix(NA_real_, nrow=nrow(valid_gbm), ncol=length(eval_times))

for(j in seq_along(eval_times)){
  t <- eval_times[j]
  h0_t <- max(base_haz$hazard[base_haz$time <= t], na.rm=TRUE)
  if(is.infinite(h0_t) || is.na(h0_t)) h0_t <- 0
  
  S_t <- exp(-h0_t * exp(valid_lp)) # survival probability
  risk_mat[, j] <- 1 - S_t # event risk
}

# Clip to [0,1]
risk_mat[risk_mat < 0] <- 0
risk_mat[risk_mat > 1] <- 1
colnames(risk_mat) <- paste0("t", eval_times)

cat("Risk matrix range:", min(risk_mat), "to", max(risk_mat), "\n")

# 5. Compute IBS with Score()
score_obj <- Score(
  object = list(GBM = risk_mat),
  formula = Surv(Survival_time, event_flag) ~ 1,
  data = valid_gbm,
  metrics = "Brier",
  summary = "ibs",
  times = eval_times,
  conservative = TRUE,
  conf.int = FALSE
)

### 6. Extract IBS SAFELY
ibs_score <- NA

if("Brier" %in% names(score_obj)){
  if("score" %in% names(score_obj$Brier)){
    # Case 1: stored as vector
    if(!is.null(score_obj$Brier$score$ibs)){
      ibs_score <- as.numeric(score_obj$Brier$score$ibs)
    }
    # Case 2: stored as matrix (common!)
    else if(!is.null(score_obj$Brier$score$AppErr$ibs)){
      ibs_score <- as.numeric(score_obj$Brier$score$AppErr$ibs)
    }
  }
}
cat("✅ IBS:", round(ibs_score, 4), "\n")


# ---------------------
# Calibration slope
# ---------------------
cal_slope_fit <- coxph(
  Surv(Survival_time, event_flag) ~ lp,
  data = valid_gbm %>% mutate(lp = valid_lp)
)

cal_slope <- as.numeric(cal_slope_fit$coef[1])
cat("✅ Calibration slope:", round(cal_slope, 3), "\n")
### 7. Save results safely

ibs_vec <- score_obj$Brier$score$IBS
ibs_score <- as.numeric(tail(ibs_vec,1))
out <- data.frame(
  IBS = ibs_score,
  Calibration_Slope = cal_slope,
  Date_Run = as.character(Sys.time())
)

write.csv(
  out,
  file.path(OUTDIR, "GBM_Calibration_IBS_Slope.csv"),
  row.names = FALSE
)

cat("📄 Saved → GBM_Calibration_IBS_Slope.csv\n")






###############################################################################
# 5.3.3 CLINICAL UTILITY — DECISION CURVE ANALYSIS (DCA)
###############################################################################

cat("\n---- 5.3.3 Decision Curve Analysis ----\n")

valid_df_dca <- valid_gbm %>%
  mutate(
    risk = pmin(pmax(1 - exp(-exp(valid_lp)), 0), 1)
  )

dca_res <- decision_curve(
  formula = event_flag ~ risk,
  data = valid_df_dca,
  family = binomial(link="logit"),
  thresholds = seq(0.01, 0.5, by=0.01),
  confidence.intervals = FALSE
)

# Save DCA CSV
write.csv(
  dca_res$derived.data,
  file.path(OUTDIR, "GBM_DCA_Results.csv"),
  row.names = FALSE
)

cat("DCA complete → GBM_DCA_Results.csv\n")

# Plot DCA
png(file.path(OUTDIR, "GBM_DCA_Plot.png"), width=800, height=600)
plot_decision_curve(dca_res, curve.names="GBM", standardize=TRUE)
dev.off()

cat("📉 DCA plot saved → GBM_DCA_Plot.png\n")


# 5. Save results
write.csv(
  as.data.frame(nri_res),
  file.path(OUTDIR, "GBM_NRI_Results.csv"),
  row.names = FALSE
)

cat("📄 NRI saved → GBM_NRI_Results.csv\n")



###############################################################
# SHAP VALUES FOR SAVED GBM SURVIVAL MODEL
# ------------------------------------------------------------
# This script:
#  1. Loads the saved GBM model
#  2. Loads the dataset
#  3. Loads best hyperparameters
#  4. Computes SHAP values (fastshap)
#  5. Saves SHAP summary and dependence plots
#  6. Saves SHAP importance table
###############################################################

rm(list = ls())

##############################
# 1 — PACKAGES
##############################
packages <- c("tidyverse","gbm","survival","fastshap","shapviz")

for(pk in packages){
  if(!requireNamespace(pk, quietly=TRUE)) install.packages(pk)
  library(pk, character.only=TRUE)
}

##############################
# 2 — PATHS
##############################
BASEDIR <- "T:/PROFID/Study8"
INDIR   <- file.path(BASEDIR, "Variable Selection & Model Development/Files")
OUTDIR  <- file.path(BASEDIR, "Model Validation and Performance/Files")

##############################
# 3 — LOAD SAVED MODEL
##############################
final_gbm <- readRDS(file.path(OUTDIR, "Final_GBM_Model.rds"))
cat("Loaded saved GBM model.\n")

##############################
# 4 — LOAD BEST HYPERPARAMETERS
##############################
tune_results <- read.csv(file.path(OUTDIR, "GBM_Tuning_Results_Random.csv"))

best_row <- tune_results[which.max(tune_results$C_index),]

best_params <- list(
  n.trees            = best_row$n.trees,
  interaction.depth  = best_row$interaction.depth,
  shrinkage          = best_row$shrinkage,
  n.minobsinnode     = best_row$n.minobsinnode,
  bag.fraction       = best_row$bag.fraction
)

cat("Best parameters loaded:\n")
print(best_params)

##############################
# 5 — LOAD DATASET
##############################
df <- read.csv(file.path(INDIR, "vs_data_complete.csv"))

# minimal preprocessing (same as during training)
df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.integer(Status),
    event_flag = ifelse(Status == 0, 0, 1)
  )

# Convert categorical to numeric for GBM
df_gbm <- df %>%
  mutate(across(where(is.character), ~as.numeric(as.factor(.))),
         across(where(is.factor),   ~as.numeric(.)))

cat("Dataset loaded and preprocessed.\n")

##############################
# 6 — EXTRACT PREDICTORS USED IN TRAINING
##############################
# (Load from your saved script or recreate)
mandatory <- c("LVEF","Age","BMI","Diabetes","eGFR")

biomarker_list <- c(
  "BUN","Cholesterol","CRP","eGFR","Haemoglobin","HbA1c",
  "HDL","IL6","LDL","NTProBNP","Potassium","Sodium",
  "Triglycerides","Troponin_T","TSH"
)
biomarkers <- intersect(biomarker_list, names(df_gbm))

ecg_list <- c("HR","PR","QRS","QTc","AV_block","AV_block_II_or_III","LBBB","RBBB")
ecg_vars <- intersect(ecg_list, names(df_gbm))

predictors <- unique(c(mandatory, biomarkers, ecg_vars))

df_shap <- df_gbm[, predictors]

##############################
# 7 — DEFINE PREDICTION FUNCTION (required by fastshap)
##############################
pred_fun <- function(model, newdata){
  predict(model, 
          newdata = newdata,
          n.trees = best_params$n.trees,
          type = "link")
}

##############################
# 8 — OPTIONALLY SUBSAMPLE FOR SPEED
##############################
set.seed(2025)
N <- min(2000, nrow(df_shap))
df_sub <- df_shap[sample(1:nrow(df_shap), N), ]

cat("Using", N, "rows for SHAP computation.\n")

##############################
# 9 — COMPUTE SHAP VALUES
##############################
set.seed(2025)
shap_vals <- fastshap::explain(
  object = final_gbm,
  X = df_sub,
  pred_wrapper = pred_fun,
  nsim = 100
)

cat("SHAP computation complete.\n")


##############################
# 10 — WRAP INTO shapviz
##############################
sv <- shapviz(shap_vals, X = df_sub)

##############################
# 11 — SAVE SHAP SUMMARY PLOT
##############################
png(file.path(OUTDIR, "GBM_SHAP_Summary.png"), width=1000, height=800)
sv_importance(sv)      # <-- CORRECT FUNCTION
dev.off()

cat("Saved: GBM_SHAP_Summary.png\n")

##############################
# 12 — SAVE SHAP DEPENDENCE PLOTS
##############################
vars_to_plot <- c("LVEF","Age","NTProBNP","Haemoglobin","Triglycerides")

for(v in vars_to_plot){
  if(v %in% colnames(df_sub)){
    png(file.path(OUTDIR, paste0("GBM_SHAP_", v, ".png")), width=900, height=700)
    sv_dependence(sv, v)
    
    dev.off()
    cat("Saved dependence plot for:", v, "\n")
  }
}

##############################
# 13 — SHAP IMPORTANCE TABLE
##############################
imp <- data.frame(
  Variable = colnames(shap_vals),
  MeanAbsSHAP = apply(abs(shap_vals), 2, mean)
) %>%
  arrange(desc(MeanAbsSHAP))

write.csv(imp, file.path(OUTDIR, "GBM_SHAP_Importance.csv"), row.names = FALSE)

cat("Saved: GBM_SHAP_Importance.csv\n")

###############################################################
cat("\n📌 SHAP analysis complete — all outputs saved in OUTDIR.\n")
###############################################################

