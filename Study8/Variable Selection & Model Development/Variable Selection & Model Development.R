=

# ---- Setup ----------------------------------------------------------------
# Install packages if needed (uncomment if first time)
install.packages(c("survival","dplyr","broom","MASS","MatchIt","tableone","mice","car","survminer","cobalt"))

library(dplyr)
library(survival)
library(broom)
library(MASS)
library(MatchIt)
library(tableone)
library(mice)
library(car)
library(survminer)

# Output folder (must exist)
OUTDIR <- "T:/PROFID/Study8/Variable Selection & Model Development/Files"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ---- Data load & combine --------------------------------------------------
cat("Loading files...\n")
file1 <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/ICD.csv", stringsAsFactors = FALSE)
file2 <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_preserved.csv", stringsAsFactors = FALSE)
file3 <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_reduced.csv", stringsAsFactors = FALSE)

file1$Group <- "ICD"
file2$Group <- "NonICD_preserved"
file3$Group <- "NonICD_reduced"

combined <- bind_rows(file1, file2, file3)
cat("Combined rows:", nrow(combined), "\n")
cat("Columns:", paste(colnames(combined), collapse = ", "), "\n")

# ---- Safe pre-filter (drop meta columns) ----------------------------------
cols_to_remove <- c("X","DB","ID","Time_zero_Y","Time_zero_Ym","Time_index_MI_CHD","HasMRI")
vs_data <- combined %>%
  dplyr::select(dplyr::all_of(setdiff(names(combined), cols_to_remove))) %>%
  dplyr::filter(!is.na(Status))


cat("Rows after filtering missing Status:", nrow(vs_data), "\n")

# ---- Univariate screening (Cox) -------------------------------------------
cat("Running univariate Cox regressions (SCD event code = 1 assumed)...\n")
surv_obj <- with(vs_data, Surv(Survival_time, Status == 1))
candidate_vars <- setdiff(names(vs_data), c("Survival_time","Status"))

uni_list <- list()
for (var in candidate_vars) {
  f <- as.formula(paste0("surv_obj ~ `", var, "`"))
  res <- tryCatch({
    fit <- coxph(f, data = vs_data)
    tidy(fit) %>% mutate(Variable = var)
  }, error = function(e) {
    NULL
  })
  if (!is.null(res)) uni_list[[var]] <- res
}
uni_results <- bind_rows(uni_list)
write.csv(uni_results, file.path(OUTDIR, "Univariate_All.csv"), row.names = FALSE)
cat("Saved Univariate_All.csv\n")

# select p < 0.20 (per protocol)
uni_selected <- uni_results %>% filter(p.value < 0.20) %>% distinct(Variable) %>% pull(Variable)
write.csv(data.frame(Selected = uni_selected), file.path(OUTDIR, "Univariate_Selected_p20.csv"), row.names = FALSE)
cat("Saved Univariate_Selected_p20.csv (p < 0.20)\n")

# ---- Variable type summary: before conversion ------------------------------
var_types <- sapply(vs_data, class)
var_summary <- data.frame(
  Variable = names(var_types),
  Type = as.character(var_types),
  UniqueValues = sapply(vs_data, function(x) length(unique(na.omit(x)))),
  ExampleValues = sapply(vs_data, function(x) paste0(head(unique(na.omit(x)),3), collapse = ", "))
)
write.csv(var_summary, file.path(OUTDIR, "Variable_Types_Before_Imputation.csv"), row.names = FALSE)
cat("Saved Variable_Types_Before_Imputation.csv\n")

# ---- Convert categorical variables to factors (controlled list) -----------
cat("Converting categorical variables to factors (Yes/No and multi-cat lists)...\n")
yesno_vars <- c("Sex","LBBB","Stroke_TIA","Diabetes","ACE_inhibitor_ARB","Anti_arrhythmic_III","Anti_platelet",
                "Beta_blockers","Diuretics","MI_location_anterior","MI_location_inferior","PCI","CABG","RBBB",
                "NSVT","AF_atrial_flutter","COPD","Hypertension","Smoking","ACE_inhibitor","ARB",
                "Aldosterone_antagonist","Anti_anginal","Anti_diabetic_insulin","Anti_diabetic_oral",
                "Calcium_antagonists","Digitalis_glycosides","FH_CAD","Anti_coagulant","Cancer",
                "AV_block","HF","Anti_diabetic","Lipid_lowering","Alcohol",
                "PCI_acute","Thrombolysis_acute","Revascularisation_acute","MI_history","MI_location_posterior")

multi_cat_vars <- c("NYHA","Diseased_arteries_num_cat","MR","MI_type","Group","Baseline_type")
coded_vars <- c("Status","CVD_risk_region","IsSWHR")

# Convert safely (only if exist)
for (v in intersect(names(vs_data), yesno_vars)) {
  vs_data[[v]] <- trimws(as.character(vs_data[[v]]))
  vs_data[[v]] <- ifelse(vs_data[[v]] %in% c("Yes","No"), vs_data[[v]], vs_data[[v]]) # keep as-is if coded differently
  vs_data[[v]] <- as.factor(vs_data[[v]])
}
for (v in intersect(names(vs_data), multi_cat_vars)) vs_data[[v]] <- as.factor(vs_data[[v]])
for (v in intersect(names(vs_data), coded_vars)) vs_data[[v]] <- as.factor(vs_data[[v]])

# Save post-conversion summary
var_types_after <- sapply(vs_data, class)
var_summary_after <- data.frame(
  Variable = names(var_types_after),
  Type = as.character(var_types_after),
  UniqueValues = sapply(vs_data, function(x) length(unique(na.omit(x)))),
  ExampleValues = sapply(vs_data, function(x) paste0(head(unique(na.omit(x)),3), collapse = ", "))
)
write.csv(var_summary_after, file.path(OUTDIR, "Variable_Types_After_Conversion.csv"), row.names = FALSE)
cat("Saved Variable_Types_After_Conversion.csv\n")

# Optional save of converted dataset
write.csv(vs_data, file.path(OUTDIR, "vs_data_after_conversion.csv"), row.names = FALSE)
cat("Saved vs_data_after_conversion.csv\n")

# ---- Missingness rules & imputation (Heather's rules) ----------------------
# - Variables with > 80% missing: drop (very sparse)
# - Variables with 50-80% missing: consider domain-specific; here we drop >80, keep <=80
# - For history flags (e.g., Cancer, MI_history) default missing -> "No" (0)
# - For remaining variables with missing <=80%: use mice (predictive mean for continuous; logistic for binary)
#
cat("Computing missing proportions...\n")
missing_prop <- sapply(vs_data, function(x) mean(is.na(x)))
missing_df <- data.frame(Variable = names(missing_prop), MissingProp = missing_prop)
write.csv(missing_df, file.path(OUTDIR, "Missingness_Summary.csv"), row.names = FALSE)
cat("Saved Missingness_Summary.csv\n")

# Drop variables with > 80% missing
drop_80 <- names(missing_prop)[missing_prop > 0.80]
if (length(drop_80) > 0) {
  cat("Dropping variables with >80% missingness:\n"); print(drop_80)
  vs_data <- vs_data %>% dplyr::select(!dplyr::all_of(drop_80))
  
}

# For certain clinical-history binary vars, set missing -> "No" (where appropriate)
# Supervisor suggested 'history of cancer' default 0; add other history flags if desired.
history_flags <- intersect(c("Cancer","MI_history","FH_CAD"), names(vs_data))
for (v in history_flags) {
  vs_data[[v]] <- as.character(vs_data[[v]])
  vs_data[[v]][is.na(vs_data[[v]])] <- "No"
  vs_data[[v]] <- as.factor(vs_data[[v]])
}
cat("Applied default 0/No to history flags (if present):", paste(history_flags, collapse=", "), "\n")

# Now: impute remaining variables with missingness <= 80% using mice
# Prepare mice-friendly dataset: only variables we will plausibly impute (exclude IDs/time/outcomes)
impute_vars <- setdiff(names(vs_data), c("Survival_time","Status"))
# Keep small set to avoid huge run times; you can expand as needed
mice_df <- vs_data %>% dplyr::select(dplyr::all_of(impute_vars))

# Convert character factors to factors so mice knows how to impute
mice_df <- mice_df %>% mutate(across(where(is.character), as.factor))

# Quick mice setup: use default methods
cat("Starting mice imputation (this may take a while)...\n")
# set small m to speed up; you can increase m later (e.g., m=5)
imp <- mice(mice_df, m = 1, maxit = 5, printFlag = TRUE)
mice_complete <- complete(imp, action = 1)

# Merge back imputed vars with time/event
vs_data_imputed <- bind_cols(
  vs_data %>% dplyr::select(Survival_time, Status),
  mice_complete
)



write.csv(vs_data_imputed, file.path(OUTDIR, "vs_data_after_imputation.csv"), row.names = FALSE)
save(imp, file = file.path(OUTDIR, "mice_imputation_object.RData"))
cat("Saved imputed dataset and mice object.\n")

# ---- Remove variables with >50% missing (if you still want a stricter filter) -
# (Optional) Here we keep variables with <=50% missing before imputation, already handled
# but keep your earlier rule: drop >50% missing if you prefer.
# vs_data_clean <- vs_data_imputed %>% select(which(colMeans(is.na(.)) <= 0.5))

# ---- Prepare dataset for modeling -----------------------------------------
# Keep mandatory predictors in model regardless
mandatory_vars <- c("LVEF","Age","Diabetes","eGFR")
selected_from_univ <- uni_selected
final_vars <- unique(c(selected_from_univ, mandatory_vars))

cols_to_keep <- intersect(c("Survival_time","Status", final_vars), names(vs_data_imputed))
vs_data_model <- vs_data_imputed %>% dplyr::select(all_of(cols_to_keep))
cat("Model dataset dimensions:", dim(vs_data_model), "\n")
write.csv(vs_data_model, file.path(OUTDIR, "vs_data_model.csv"), row.names = FALSE)

# ---- Drop single-level variables + factor conversion ----------------------
# drop any variable with only one unique non-NA level
single_level_vars <- names(vs_data_model)[sapply(vs_data_model, function(x) length(unique(na.omit(x))) <= 1)]
if (length(single_level_vars) > 0) {
  cat("Dropping single-level variables (non-informative):\n"); print(single_level_vars)
  vs_data_model <- vs_data_model %>% select(-all_of(single_level_vars))
}

# Ensure categorical columns are factors
vs_data_model <- vs_data_model %>% mutate(across(where(is.character), as.factor))

# ---- Create complete-case dataset for stepwise selection -------------------
vars_in_formula <- c("Survival_time","Status", setdiff(names(vs_data_model), c("Survival_time","Status")))
vs_data_complete <- vs_data_model %>% filter(complete.cases(.))
cat("Complete-case rows:", nrow(vs_data_complete), " (out of ", nrow(vs_data_model), ")\n", sep="")

save(vs_data_complete, file = file.path(OUTDIR, "vs_data_complete.RData"))
write.csv(vs_data_complete, file.path(OUTDIR, "vs_data_complete.csv"), row.names = FALSE)

# ---- Fit full Cox model ---------------------------------------------------
formula_str <- paste0("Surv(Survival_time, Status == 1) ~ ",
                      paste(setdiff(names(vs_data_complete), c("Survival_time","Status")), collapse = " + "))
cat("Full model formula:\n", formula_str, "\n")

full_model_complete <- coxph(as.formula(formula_str), data = vs_data_complete)
sink(file.path(OUTDIR, "Cox_Multivariate_Full_Model_Summary.txt"))
print(summary(full_model_complete))
sink()
write.csv(tidy(full_model_complete, exponentiate = TRUE, conf.int = TRUE),
          file.path(OUTDIR, "Cox_Multivariate_Full_Model.csv"), row.names = FALSE)
cat("Saved full model summary and CSV.\n")

# ---- Backward stepwise selection (AIC) ------------------------------------
cat("Running backward stepwise selection (AIC-based). This can take time...\n")
step_model <- step(full_model_complete, direction = "backward", trace = TRUE)
sink(file.path(OUTDIR, "Cox_Multivariate_Stepwise_Model_Summary.txt"))
print(summary(step_model))
sink()
write.csv(tidy(step_model, exponentiate = TRUE, conf.int = TRUE),
          file.path(OUTDIR, "Cox_Multivariate_Stepwise_Model.csv"), row.names = FALSE)
cat("Saved stepwise model outputs.\n")

# ---- ANOVA model comparison (LRT) -----------------------------------------
anova_result <- anova(full_model_complete, step_model, test = "LRT")
# anova_result may be of class "anova" — coerce to data frame for saving
anova_df <- as.data.frame(anova_result)
write.csv(anova_df, file.path(OUTDIR, "Cox_ANOVA_Model_Comparison.csv"), row.names = TRUE)
cat("Saved Cox_ANOVA_Model_Comparison.csv\n")
print(anova_result)

# ---- Multicollinearity / proxy removal -----------------------------------
# numeric correlation matrix
numeric_vars <- vs_data_complete %>% dplyr::select(where(is.numeric)) %>% dplyr::select(-Survival_time, -Status)
cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")
write.csv(cor_matrix, file.path(OUTDIR, "Correlation_Matrix.csv"), row.names = TRUE)
cat("Saved Correlation_Matrix.csv\n")

# find high-correlation pairs (|r| >= 0.8)
high_corr_pairs <- which(abs(cor_matrix) >= 0.8 & abs(cor_matrix) < 1, arr.ind = TRUE)
if (nrow(high_corr_pairs) == 0) {
  cat("No high-correlation pairs (|r|>=0.8) detected among numeric variables.\n")
} else {
  corr_list <- data.frame(
    Var1 = rownames(cor_matrix)[high_corr_pairs[,1]],
    Var2 = colnames(cor_matrix)[high_corr_pairs[,2]],
    Correlation = cor_matrix[high_corr_pairs]
  )
  write.csv(corr_list, file.path(OUTDIR, "Highly_Correlated_Pairs.csv"), row.names = FALSE)
  cat("Saved Highly_Correlated_Pairs.csv\n")
  print(corr_list)
  # Example: choose which to drop: if Cholesterol & LDL correlated, drop LDL (example)
  if ("LDL" %in% names(vs_data_complete) & "Cholesterol" %in% names(vs_data_complete)) {
    cat("Dropping LDL (proxy for Cholesterol) as example decision.\n")
    vs_data_final <- vs_data_complete %>% dplyr::select(-LDL)
    # update formula removing LDL
    formula_final <- update(formula(step_model), . ~ . - LDL)
    final_model <- coxph(formula_final, data = vs_data_final)
  } else {
    vs_data_final <- vs_data_complete
    final_model <- step_model
  }
}

sink(file.path(OUTDIR, "Cox_Final_Model_Summary.txt"))
print(summary(final_model))
sink()
write.csv(tidy(final_model, exponentiate = TRUE, conf.int = TRUE),
          file.path(OUTDIR, "Cox_Final_Model.csv"), row.names = FALSE)
cat("Saved final model summary and CSV.\n")

# ---- Diagnostics (VIF for numeric variables in final model) --------------
# compute VIF where appropriate
cat("Computing VIF (if applicable)...\n")
# Build model matrix from final model to compute VIF (requires no-surv object)
f_terms <- terms(final_model)
mm <- model.matrix(delete.response(f_terms), data = vs_data_final)
# drop intercept column
mm <- mm[, colnames(mm) != "(Intercept)", drop=FALSE]
vif_vals <- tryCatch({
  vif(lm(mm ~ 1))
}, error=function(e) NA)
write.csv(data.frame(VIF = vif_vals), file.path(OUTDIR, "Final_Model_VIF.csv"), row.names = TRUE)
cat("Saved Final_Model_VIF.csv\n")

# ---- Stage 3: Treatment effect modelling (propensity score example) ------
# Example: anti-coagulant as treatment
if ("Anti_coagulant" %in% names(vs_data_final)) {
  cat("Running propensity score example for Anti_coagulant (treated vs untreated)...\n")
  # build PS model (use baseline confounders)
  confounders <- intersect(c("Age","Sex","LVEF","Diabetes","eGFR","BMI","MI_history","HF"), names(vs_data_final))
  ps_formula <- as.formula(paste("Anti_coagulant ~", paste(confounders, collapse = " + ")))
  ps_model <- glm(ps_formula, data = vs_data_final, family = binomial)
  vs_data_final$pscore <- predict(ps_model, type = "response")
  write.csv(data.frame(summary(ps_model)$coefficients), file.path(OUTDIR, "PropensityScore_Model_Coefficients.csv"))
  
  # Match using nearest neighbour
  matchit_obj <- matchit(ps_formula, data = vs_data_final, method = "nearest", distance = vs_data_final$pscore)
  matched <- match.data(matchit_obj)
  write.csv(matched, file.path(OUTDIR, "PropensityScore_Matched_Anti_coagulant.csv"), row.names = FALSE)
  
  # Balance check (SMD)
  tab <- CreateTableOne(vars = confounders, strata = "Anti_coagulant", data = matched, test=FALSE)
  smd_df <- as.data.frame(ExtractSmd(tab))
  write.csv(smd_df, file.path(OUTDIR, "PS_Balance_SMD_Anti_coagulant.csv"), row.names = TRUE)
  cat("Saved propensity score outputs and SMD.\n")
} else {
  cat("Anti_coagulant not in dataset; skipping PS example.\n")
}





# Stage 3: Treatment Effect Modelling and Visualization (Standalone)
# ---------------------------------------------------------------

########################################################################

# ---- Setup ----
library(dplyr)
library(MatchIt)
library(tableone)
library(ggplot2)
library(cobalt)
library(survival)
library(broom)
library(survminer)

OUTDIR <- "T:/PROFID/Study8/Variable Selection & Model Development/Files"
setwd(OUTDIR)

# ---- Load your cleaned dataset ----
cat(" Loading final dataset...\n")
vs_data_final <- read.csv("vs_data_complete.csv")

cat(" Loaded data with", nrow(vs_data_final), "rows and",
    ncol(vs_data_final), "columns.\n")

# Ensure correct types
vs_data_final <- vs_data_final %>%
  mutate(
    across(where(is.character), as.factor),
    Status = as.numeric(Status),
    Survival_time = as.numeric(Survival_time)
  )

# ---- Helper function ----
run_ps_and_visualize <- function(data, treatment_var, confounders, outdir) {
  if (!(treatment_var %in% names(data))) {
    cat("\n Variable", treatment_var, "not found — skipping.\n")
    return(NULL)
  }
  
  cat("\n=====================================================\n")
  cat(" Propensity Score Analysis for:", treatment_var, "\n")
  cat("=====================================================\n")
  
  #  Propensity score model
  ps_formula <- as.formula(paste(treatment_var, "~", paste(confounders, collapse = " + ")))
  ps_model <- glm(ps_formula, data = data, family = binomial)
  data$pscore <- predict(ps_model, type = "response")
  
  #  Visualize PS distribution before matching
  p1 <- ggplot(data, aes(x = pscore, fill = !!sym(treatment_var))) +
    geom_density(alpha = 0.4) +
    labs(title = paste("Propensity Score Distribution —", treatment_var),
         x = "Propensity Score", fill = treatment_var) +
    theme_minimal()
  ggsave(file.path(outdir, paste0("PS_Distribution_", treatment_var, ".png")), p1, width = 7, height = 5)
  
  #  Matching
  matchit_obj <- matchit(ps_formula, data = data, method = "nearest", distance = data$pscore)
  matched <- match.data(matchit_obj)
  cat("Matched sample size:", nrow(matched), "\n")
  
  #  Love plot (SMD balance)
  png(file.path(outdir, paste0("LovePlot_", treatment_var, ".png")), width = 700, height = 500)
  love.plot(matchit_obj, threshold = 0.1, abs = TRUE,
            var.order = "unadjusted",
            title = paste("Love Plot — Covariate Balance for", treatment_var))
  dev.off()
  
  #  Save SMD table
  tab <- CreateTableOne(vars = confounders, strata = treatment_var, data = matched, test = FALSE)
  smd_df <- as.data.frame(ExtractSmd(tab))
  write.csv(smd_df, file.path(outdir, paste0("PS_Balance_SMD_", treatment_var, ".csv")), row.names = TRUE)
  cat("Saved SMD table for", treatment_var, "\n")
  
  #  Cox model on matched data
  cox_model <- coxph(Surv(Survival_time, Status == 1) ~ get(treatment_var), data = matched)
  cox_summary <- summary(cox_model)
  hr <- exp(cox_summary$coefficients[,1])
  ci_low <- cox_summary$conf.int[,3]
  ci_high <- cox_summary$conf.int[,4]
  pval <- cox_summary$coefficients[,5]
  
  cat("\n⚖️  Cox model summary for", treatment_var, "\n")
  cat("Hazard Ratio (HR):", round(hr,3),
      " | 95% CI:", round(ci_low,3), "-", round(ci_high,3),
      " | p =", signif(pval,3), "\n")
  
  # Save matched data & model coefficients
  write.csv(matched, file.path(outdir, paste0("PS_Matched_", treatment_var, ".csv")), row.names = FALSE)
  write.csv(data.frame(summary(ps_model)$coefficients),
            file.path(outdir, paste0("PS_Model_Coefficients_", treatment_var, ".csv")), row.names = TRUE)
  
  invisible(list(matched = matched, cox = cox_summary))
}

# ---- Run for each treatment ----
confounders <- intersect(c("Age","Sex","LVEF","Diabetes","eGFR","BMI","MI_history","HF"), names(vs_data_final))

#  Anticoagulation
res_anticoag <- run_ps_and_visualize(vs_data_final, "Anti_coagulant", confounders, OUTDIR)

#  Rhythm control
res_rhythm <- run_ps_and_visualize(vs_data_final, "Anti_arrhythmic_III", confounders, OUTDIR)

#  Rate control
res_rate <- run_ps_and_visualize(vs_data_final, "Beta_blockers", confounders, OUTDIR)

cat("\nStage 3 Treatment Effect Modelling complete.\n")
