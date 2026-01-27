rm(list = ls())  # clears everything

packages <- c("tidyverse", "survival", "ranger", "gbm",
              "caret","ggfortify" ,"survAUC", "Hmisc", "ggplot2", "viridis", "survminer")
install_if_missing <- function(pk) if(!requireNamespace(pk, quietly=TRUE)) install.packages(pk)
invisible(lapply(packages, install_if_missing))
lapply(packages, library, character.only = TRUE)
install.packages("ranger")
library(ranger)

install.packages("Hmisc")  # run once
library(Hmisc)


# ---- Paths ------------------------------------------------------------


BASEDIR <- "T:/PROFID/Study8"
INFILE  <- file.path(BASEDIR, "Variable Selection & Model Development/Files/vs_data_complete.csv")
OUTDIR  <- file.path(BASEDIR, "AdvancedStats/Files_ML")
MODELFILE <- file.path(OUTDIR, "RSF_model_light.rds")

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
df <- read.csv(INFILE, stringsAsFactors = FALSE)
# ---- Preprocessing ---------------------------------------------------
df <- df %>%
  dplyr::mutate(Survival_time = as.numeric(Survival_time),
                Status = as.integer(Status),
                event_flag = ifelse(Status == 0, 0, 1)) %>%
  dplyr::filter(!is.na(Survival_time), !is.na(Status))


# ---- Predictor Selection ---------------------------------------------
mandatory <- c("LVEF","Age","BMI","Diabetes","eGFR")
biomarker_list <- c("BUN","Cholesterol","CRP","eGFR","Haemoglobin","HbA1c",
                    "HDL","IL6","LDL","NTProBNP","Potassium","Sodium",
                    "Triglycerides","Troponin_T","TSH")

biomarkers <- intersect(biomarker_list, names(df))

cat("Biomarkers present in dataset:\n")
print(biomarkers)

ecg_list <- c("HR","PR","QRS","QTc","AV_block","AV_block_II_or_III",
              "LBBB","RBBB")

ecg_vars <- intersect(ecg_list, names(df))

cat("ECG variables present in dataset:\n")
print(ecg_vars)

write.csv(data.frame(Biomarkers = biomarkers), file.path(OUTDIR, "Biomarkers_present.csv"), row.names = FALSE)

write.csv(data.frame(ECG_vars = ecg_vars), file.path(OUTDIR, "ECG_vars_present.csv"), row.names = FALSE)



predictors <- unique(c(mandatory, biomarkers, ecg_vars))
cat("🧠 Predictors used:", paste(predictors, collapse=", "), "\n")

# ---- 70:30 Split (Development vs Validation) -------------------------
train_idx <- createDataPartition(df$event_flag, p = 0.7, list = FALSE)
train <- df[train_idx, ]
valid <- df[-train_idx, ]
cat("🧩 Split complete: Train =", nrow(train), "| Valid =", nrow(valid), "\n")

# ---- Ensure characters → factors ------------------------------------
train <- train %>% mutate(across(where(is.character), as.factor))
valid <- valid %>% mutate(across(where(is.character), as.factor))

# ---- Create Modeling Formula ----------------------------------------
model_formula <- as.formula(paste0(
  "Surv(Survival_time, event_flag) ~ ", paste(predictors, collapse = " + ")
))
cat("📜 Formula:", deparse(model_formula), "\n")

########################################################################
# 🌲 1️⃣ RANDOM SURVIVAL FOREST (lightweight + memory-optimized)
########################################################################

cat("\n🌲 Fitting Random Survival Forest (lightweight mode)...\n")

set.seed(2025)
rsf_fit <- ranger(
  formula = model_formula,
  data = train,
  num.trees = 10,                # keep very small for testing
  mtry = floor(sqrt(length(predictors))),  # fewer vars per split
  importance = "impurity",
  splitrule = "logrank",         # less RAM than extratrees
  sample.fraction = 0.5,         # use only 50% of data per tree
  num.threads = 2,               # limit threads to control memory
  write.forest = TRUE
)

saveRDS(rsf_fit, file.path(OUTDIR, "RSF_model_light.rds"))
cat("✅ RSF model trained with", rsf_fit$num.trees, "trees.\n")



# ---- Variable Importance ---------------------------------------------
imp <- data.frame(Variable = names(rsf_fit$variable.importance),
                  Importance = rsf_fit$variable.importance) %>%
  arrange(desc(Importance))
write.csv(imp, file.path(OUTDIR, "RSF_Variable_Importance.csv"), row.names = FALSE)

imp %>%
  top_n(10, wt = Importance) %>%
  ggplot(aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_col(fill = "#2E86AB") +
  coord_flip() +
  labs(title = "Top 10 Important Variables — Random Survival Forest",
       x = "Variable", y = "Importance") +
  theme_minimal()
ggsave(file.path(OUTDIR, "RSF_Variable_Importance_Top10.png"), width = 7, height = 5)



########################################################################
# 📈 MODEL EVALUATION (Chunked prediction)
########################################################################
chunk_predict <- function(model, data, chunk_size = 3000) {
  preds <- numeric(nrow(data))
  for (i in seq(1, nrow(data), by = chunk_size)) {
    idx <- i:min(i + chunk_size - 1, nrow(data))
    cat("  Predicting rows", i, "to", max(idx), "...\n")
    
    sub_pred <- tryCatch({
      predict(model, data = data[idx, ])
    }, error = function(e) {
      warning(paste("⚠️ Prediction failed for chunk", i, ":", e$message))
      return(NULL)
    })
    
    if (!is.null(sub_pred$predictions)) {
      preds[idx] <- as.numeric(sub_pred$predictions)
    } else if (!is.null(sub_pred$chf)) {
      preds[idx] <- as.numeric(sub_pred$chf[, ncol(sub_pred$chf)])
    } else if (!is.null(sub_pred$survival)) {
      preds[idx] <- as.numeric(1 - sub_pred$survival[, ncol(sub_pred$survival)])
    }
  }
  return(preds)
}

# ---- Run Prediction ---------------------------------------------------
cat("\n⚡ Running validation predictions (chunked mode)...\n")
valid_pred <- chunk_predict(rsf_fit, valid, chunk_size = 2000)
cat("✅ Predictions complete for", length(valid_pred), "records.\n")

# ---- Evaluation -------------------------------------------------------
eval_df <- valid %>%
  mutate(pred = valid_pred) %>%
  filter(!is.na(Survival_time), !is.na(event_flag), !is.na(pred))

cat("✅ Evaluation dataset ready with", nrow(eval_df), "rows.\n")

# Compute C-index
c_index_rsf <- rcorr.cens(-eval_df$pred,
                          Surv(eval_df$Survival_time, eval_df$event_flag))["C Index"]

cat("\n🔹 Harrell's C-index:", round(c_index_rsf, 3), "\n")

# ---- Save Results -----------------------------------------------------
eval_out <- data.frame(
  Model = "Random Survival Forest",
  C_Index = round(c_index_rsf, 3),
  Date_Run = Sys.time(),
  Valid_Rows = nrow(valid)
)

write.csv(eval_df, file.path(OUTDIR, "RSF_Predictions_Validation.csv"), row.names = FALSE)
write.csv(eval_out, file.path(OUTDIR, "RSF_Model_Performance_OnlyEval.csv"), row.names = FALSE)

cat("\n💾 Results saved to:\n",
    " - RSF_Predictions_Validation.csv\n",
    " - RSF_Model_Performance_OnlyEval.csv\n")



########################################################################
# 🚀 2️⃣ GRADIENT BOOSTING MODEL (GBM) for Survival Analysis
########################################################################

cat("\n🔥 Fitting Gradient Boosting Model (GBM) for Survival data...\n")

set.seed(2025)

# ---- Prepare Data ----------------------------------------------------
# Use the same train and validation sets you created earlier
gbm_formula <- model_formula

# GBM requires numeric predictors — convert factors to numeric
train_gbm <- train %>% mutate(across(where(is.factor), as.numeric))
valid_gbm <- valid %>% mutate(across(where(is.factor), as.numeric))

# ---- Fit Gradient Boosting Model ------------------------------------
gbm_fit <- gbm(
  formula = gbm_formula,
  data = train_gbm,
  distribution = "coxph",   # survival boosting
  n.trees = 100,            # increase to 1000+ for full model later
  interaction.depth = 3,    # allows 3-way interactions
  shrinkage = 0.01,         # learning rate
  n.minobsinnode = 50,      # min obs per terminal node
  bag.fraction = 0.6,       # sub-sampling rate (to reduce overfitting)
  train.fraction = 1.0,     # use all train data
  keep.data = FALSE,
  verbose = TRUE
)

saveRDS(gbm_fit, file.path(OUTDIR, "GBM_model_light.rds"))
cat("✅ GBM model trained with", gbm_fit$n.trees, "trees.\n")

# ---- Variable Importance --------------------------------------------
cat("\n📊 Calculating variable importance (GBM)...\n")
gbm_imp <- summary(gbm_fit, plotit = FALSE)
write.csv(gbm_imp, file.path(OUTDIR, "GBM_Variable_Importance.csv"), row.names = FALSE)

# Plot Top 10 important variables
top10_gbm <- gbm_imp %>% top_n(10, wt = rel.inf)
ggplot(top10_gbm, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_col(fill = "#F39C12") +
  coord_flip() +
  labs(title = "Top 10 Important Variables — Gradient Boosting (GBM)",
       x = "Variable", y = "Relative Influence (%)") +
  theme_minimal()
ggsave(file.path(OUTDIR, "GBM_Variable_Importance_Top10.png"), width = 7, height = 5)

# ---- Model Evaluation -----------------------------------------------
cat("\n📈 Evaluating GBM model on validation set...\n")

# Generate risk predictions
gbm_pred <- predict(gbm_fit, newdata = valid_gbm, n.trees = gbm_fit$n.trees, type = "link")

# Combine for evaluation
eval_gbm <- valid_gbm %>%
  mutate(pred = gbm_pred) %>%
  filter(!is.na(Survival_time), !is.na(event_flag), !is.na(pred))

# Compute Harrell’s C-index
c_index_gbm <- rcorr.cens(-eval_gbm$pred,
                          Surv(eval_gbm$Survival_time, eval_gbm$event_flag))["C Index"]

cat("🔹 GBM C-index (Harrell's):", round(c_index_gbm, 3), "\n")

# ---- Save Results ----------------------------------------------------
eval_out_gbm <- data.frame(
  Model = "Gradient Boosting Machine (GBM)",
  C_Index = round(c_index_gbm, 3),
  Date_Run = Sys.time(),
  Valid_Rows = nrow(valid_gbm)
)

write.csv(eval_gbm, file.path(OUTDIR, "GBM_Predictions_Validation.csv"), row.names = FALSE)
write.csv(eval_out_gbm, file.path(OUTDIR, "GBM_Model_Performance.csv"), row.names = FALSE)

cat("\n💾 Results saved to:\n",
    " - GBM_Predictions_Validation.csv\n",
    " - GBM_Model_Performance.csv\n")

########################################################################
cat("\n✅ Gradient Boosting survival model completed successfully!\n")
########################################################################
########################################################################
# ⚡ Advanced Gradient Boosting Model (GBM) for Survival Analysis
########################################################################

cat("\n🔥 Fitting Advanced Gradient Boosting Model (GBM) for Survival data...\n")

set.seed(2025)

# ---- Prepare Data ----------------------------------------------------
# GBM requires numeric predictors
train_gbm <- train %>% mutate(across(where(is.factor), as.numeric))
valid_gbm <- valid %>% mutate(across(where(is.factor), as.numeric))

# ---- Fit Gradient Boosting Model ------------------------------------
gbm_fit <- gbm(
  formula = model_formula,
  data = train_gbm,
  distribution = "coxph",       # for survival/time-to-event analysis
  n.trees = 1500,               # 🔹 high tree count for stronger learning
  interaction.depth = 4,        # 🔹 allows 4-way interactions
  shrinkage = 0.005,            # 🔹 smaller learning rate = more stability
  n.minobsinnode = 30,          # 🔹 more flexible tree splits
  bag.fraction = 0.7,           # 🔹 stochastic boosting (reduces overfitting)
  train.fraction = 1.0,
  keep.data = FALSE,
  verbose = TRUE
)

saveRDS(gbm_fit, file.path(OUTDIR, "GBM_model_advanced.rds"))
cat("✅ GBM model trained with", gbm_fit$n.trees, "trees.\n")

########################################################################
# 📊 VARIABLE IMPORTANCE
########################################################################

cat("\n📊 Calculating variable importance (GBM)...\n")
gbm_imp <- summary(gbm_fit, plotit = FALSE)
write.csv(gbm_imp, file.path(OUTDIR, "GBM_Variable_Importance_Advanced.csv"), row.names = FALSE)

# Plot top 10 variables
top10_gbm <- gbm_imp %>% top_n(10, wt = rel.inf)
ggplot(top10_gbm, aes(x = reorder(var, rel.inf), y = rel.inf)) +
  geom_col(fill = "#D35400") +
  coord_flip() +
  labs(title = "Top 10 Important Variables — Advanced GBM",
       x = "Variable", y = "Relative Influence (%)") +
  theme_minimal(base_size = 13)
ggsave(file.path(OUTDIR, "GBM_Variable_Importance_Top10_Advanced.png"), width = 7, height = 5)

########################################################################
# 📈 MODEL EVALUATION
########################################################################

cat("\n📈 Evaluating GBM model on validation set...\n")

# Generate predictions
gbm_pred <- predict(gbm_fit, newdata = valid_gbm,
                    n.trees = gbm_fit$n.trees, type = "link")

# Combine results
eval_gbm <- valid_gbm %>%
  mutate(pred = gbm_pred) %>%
  filter(!is.na(Survival_time), !is.na(event_flag), !is.na(pred))

# Compute C-index
c_index_gbm <- rcorr.cens(-eval_gbm$pred,
                          Surv(eval_gbm$Survival_time, eval_gbm$event_flag))["C Index"]

cat("🔹 GBM (advanced) C-index (Harrell's):", round(c_index_gbm, 3), "\n")

########################################################################
# 💾 SAVE RESULTS
########################################################################

eval_out_gbm <- data.frame(
  Model = "Gradient Boosting Machine (Advanced GBM)",
  C_Index = round(c_index_gbm, 3),
  Date_Run = Sys.time(),
  Valid_Rows = nrow(valid_gbm)
)

write.csv(eval_gbm, file.path(OUTDIR, "GBM_Predictions_Validation_Advanced.csv"), row.names = FALSE)
write.csv(eval_out_gbm, file.path(OUTDIR, "GBM_Model_Performance_Advanced.csv"), row.names = FALSE)

cat("\n💾 Results saved to:\n",
    " - GBM_Predictions_Validation_Advanced.csv\n",
    " - GBM_Model_Performance_Advanced.csv\n")

cat("\n✅ Advanced Gradient Boosting Model completed successfully!\n")

# ---- PCA on biomarkers block (if ≥2 biomarkers present) ---------------
cat("Running PCA on biomarker block if possible...\n")
bio_vars <- biomarkers
if (length(bio_vars) >= 2) {
  # take rows with complete biomarker values (PCA needs no NA)
  pca_df <- df %>% dplyr::select(all_of(bio_vars)) %>% na.omit()
  if (nrow(pca_df) < 10) {
    cat("Not enough complete biomarker rows for PCA (need >=10); skipping PCA.\n")
  } else {
    # scale and run PCA
    pca_res <- prcomp(pca_df, center = TRUE, scale. = TRUE)
    # save loadings and scores
    loadings <- as.data.frame(pca_res$rotation)
    write.csv(loadings, file.path(OUTDIR, "Biomarker_PCA_Loadings.csv"), row.names = TRUE)
    scores <- as.data.frame(pca_res$x)
    write.csv(scores, file.path(OUTDIR, "Biomarker_PCA_Scores.csv"), row.names = FALSE)
    
    # Scree plot
    png(file.path(OUTDIR, "Biomarker_PCA_Scree.png"), width = 800, height = 600)
    plot(pca_res, type = "l", main = "PCA Scree Plot — Biomarkers")
    dev.off()
    
    # Save variance explained
    var_expl <- summary(pca_res)$importance
    write.csv(var_expl, file.path(OUTDIR, "Biomarker_PCA_VarianceExplained.csv"), row.names = TRUE)
    
    cat("PCA completed and saved: loadings, scores, scree, variance explained.\n\n")
  }
} else {
  cat("Fewer than 2 biomarkers present; skipping PCA.\n\n")
}

library(ggfortify)

autoplot(pca_res, 
         data = pca_df, 
         loadings = TRUE, 
         loadings.label = TRUE,
         colour = "black",
         loadings.label.size = 5) +
  ggtitle("Biomarker PCA Biplot") +
  theme_minimal(base_size = 14)

loadings <- as.data.frame(pca_res$rotation)
loadings$Variable <- rownames(loadings)

library(reshape2)
load_melt <- melt(loadings, id.vars = "Variable")

ggplot(load_melt, aes(x = Variable, y = value, fill = variable)) +
  geom_bar(stat = "identity", position = "dodge") +
  labs(title = "PCA Loadings — Biomarkers",
       y = "Loading Value",
       x = "Biomarker") +
  theme_minimal(base_size = 14)


# ============================
# 📌 PCA HEATMAP (Biomarkers)
# ============================

library(reshape2)
library(ggplot2)
library(viridis)

if (exists("pca_res")) {
  
  # Extract loadings and reshape
  loadings <- as.data.frame(pca_res$rotation)
  loadings$Biomarker <- rownames(loadings)
  loadings_long <- melt(loadings, id.vars = "Biomarker",
                        variable.name = "PC", value.name = "Loading")
  
  # Create heatmap
  p_heat <- ggplot(loadings_long, aes(x = PC, y = Biomarker, fill = Loading)) +
    geom_tile(color = "white") +
    scale_fill_viridis(option = "plasma", direction = 1) +
    theme_minimal(base_size = 14) +
    labs(
      title = "PCA Heatmap — Biomarker Loadings",
      x = "Principal Component",
      y = "Biomarker",
      fill = "Loading"
    ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      plot.title = element_text(face = "bold", size = 18)
    )
  
  print(p_heat)
  
  # Save file
  ggsave(file.path(OUTDIR, "Biomarker_PCA_Heatmap.png"),
         p_heat, width = 8, height = 5, dpi = 300)
  
  cat("\n✅ PCA heatmap saved to Biomarker_PCA_Heatmap.png\n")
} else {
  cat("⚠️ PCA not run — heatmap skipped.\n")
}
