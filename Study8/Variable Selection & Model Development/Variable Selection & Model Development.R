install.packages(c("survival","dplyr","broom","MASS","MatchIt","tableone"))

setwd("T:/PROFID/Study8/Variable Selection & Model Development")
library(survival); library(dplyr); library(broom); library(MASS)
library(MatchIt); library(tableone)

library(survival)
library(survminer)
library(dplyr)
library(broom)
library(cmprsk)

# combining the all of them together

file1 <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/ICD.csv")
file2 <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_preserved.csv")
file3 <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_reduced.csv")


nrow(file1); nrow(file2); nrow(file3)

file1$Group <- "ICD"
file2$Group <- "NonICD_preserved"
file3$Group <- "NonICD_reduced"

combined <- bind_rows(file1, file2, file3)

colnames(combined)

#  Univariate screening


library(dplyr)

# Columns we might want to remove
cols_to_remove <- c("X", "DB", "ID", "Time_zero_Y", "Time_zero_Ym",
                    "Time_index_MI_CHD", "HasMRI")

# Filter dataset safely — only remove columns that exist
vs_data <- combined %>%
  dplyr::select(setdiff(names(.), cols_to_remove)) %>%
  dplyr::filter(!is.na(Status))
surv_obj <- with(vs_data, Surv(Survival_time, Status == 1))

# Identify all candidate predictor variables
candidate_vars <- setdiff(names(vs_data), c("Survival_time", "Status"))

# Initialize a list to store results
uni_list <- list()

# Loop through each variable for univariate Cox regression
for (var in candidate_vars) {
  formula_str <- paste0("surv_obj ~ ", var)
  model_formula <- as.formula(formula_str)
  
  # Try running Cox model; skip variables that fail
  result <- tryCatch({
    fit <- coxph(model_formula, data = vs_data)
    tidy(fit) %>%
      mutate(Variable = var)
  }, error = function(e) NULL)
  
  if (!is.null(result)) uni_list[[var]] <- result
}

# Combine all results into one data frame
uni_results <- bind_rows(uni_list)

# Save full results to CSV
write.csv(uni_results,
          "T:/PROFID/Study8/Variable Selection & Model Development/Files/Univariate_All.csv",
          row.names = FALSE)

# Filter variables with p < 0.20
uni_selected <- uni_results %>%
  filter(p.value < 0.20) %>%
  distinct(Variable) %>%
  pull(Variable)

# Save screened variables (p < 0.20) to CSV
write.csv(data.frame(Selected = uni_selected),
          "T:/PROFID/Study8/Variable Selection & Model Development/Files/Univariate_Selected_p20.csv",
          row.names = FALSE)


# Step 2: Multivariate Model

library(dplyr)
library(survival)
library(broom)
library(tidyr)


selected_vars <- read.csv("T:/PROFID/Study8/Variable Selection & Model Development/Files/Univariate_Selected_p20.csv",
                          stringsAsFactors = FALSE)$Selected

mandatory_vars <- c("LVEF", "Age", "Diabetes", "eGFR")

final_vars <- unique(c(selected_vars, mandatory_vars))
final_vars   # quick print to inspect

# safe subsetting of combined to keep only needed cols + outcome/time
cols_to_keep <- intersect(c("Survival_time", "Status", final_vars), names(combined))

# Use explicit dplyr::select to avoid masking
vs_data <- dplyr::select(combined, dplyr::all_of(cols_to_keep))

# Quick checks
cat("Columns kept:\n"); print(cols_to_keep)
cat("Rows in dataset:", nrow(vs_data), "\n")



# compute missing proportions for each column



missing_prop <- sapply(vs_data, function(x) mean(is.na(x)))
missing_df <- data.frame(Variable = names(missing_prop), MissingProp = missing_prop)
missing_df <- missing_df[order(-missing_df$MissingProp), ]

print(missing_df)   # view all; top rows show worst missingness


# Remove variables with >50% missing values


vs_data_clean <- vs_data[, which(colMeans(is.na(vs_data)) <= 0.5)]

cat("Kept variables with ≤50% missingness. Remaining variables:", ncol(vs_data_clean), "\n")


#  Prepare clean dataset for Cox modeling


vs_data_clean <- vs_data_clean %>%
  filter(!is.na(Status), !is.na(Survival_time)) %>%
  mutate(
    Sex = as.factor(Sex),
    Diabetes = as.factor(Diabetes),
    Group = as.factor(Group)
  )


# Identify and remove variables with only one level
# (prevents contrasts error in Cox model)


single_level_vars <- names(vs_data_clean)[
  sapply(vs_data_clean, function(x) length(unique(na.omit(x))) == 1)
]

if (length(single_level_vars) > 0) {
  cat("⚠️ Removing variables with only one unique value:\n")
  print(single_level_vars)
  vs_data_clean <- vs_data_clean %>% select(-all_of(single_level_vars))
} else {
  cat(" single-level variables found — continuing.\n")
}



#  Build full multivariate Cox model


library(survival)

# Automatically create Cox formula
formula_str <- paste0(
  "Surv(Survival_time, Status == 1) ~ ",
  paste(setdiff(names(vs_data_clean), c("Survival_time", "Status")), collapse = " + ")
)

cat("Model formula:\n", formula_str, "\n")

# --- Diagnostic: find factors with only one level ---
single_level_factors <- names(vs_data_clean)[
  sapply(vs_data_clean, function(x) is.factor(x) && length(unique(na.omit(x))) < 2)
]

if (length(single_level_factors) > 0) {
  cat("⚠️ Dropping factor variables with only one unique level:\n")
  print(single_level_factors)
  vs_data_clean <- vs_data_clean %>% select(-all_of(single_level_factors))
} else {
  cat("✅ No single-level factor variables found.\n")
}

# Verify
str(vs_data_clean)

# --- Step: Convert character columns to factors and remove all-NA variables ---

# Drop columns that are completely NA
vs_data_clean <- vs_data_clean[, colSums(!is.na(vs_data_clean)) > 0]

# Convert all character columns to factors
vs_data_clean <- vs_data_clean %>%
  mutate(across(where(is.character), as.factor))


# Double-check
cat("After cleaning: ", ncol(vs_data_clean), "variables remain.\n")
str(vs_data_clean)



# Fit full model
full_model <- coxph(as.formula(formula_str), data = vs_data_clean)
summary(full_model)


  #  Fix for stepwise selection (ensure complete cases)


# Identify variables used in the full model
vars_in_model <- all.vars(as.formula(formula_str))

# Keep only complete cases for those variables
vs_data_complete <- vs_data_clean %>%
  dplyr::select(all_of(vars_in_model)) %>%
  dplyr::filter(complete.cases(.))

# Save the dataset used in the final model (complete cases)
save(vs_data_complete, file = "T:/PROFID/Study8/Variable Selection & Model Development/Files/vs_data_complete.RData")

# Optional: also save as CSV for inspection
write.csv(vs_data_complete,
          "T:/PROFID/Study8/Variable Selection & Model Development/Files/vs_data_complete.csv",
          row.names = FALSE)

cat("Complete-case dataset created for stepwise selection.\n")
cat("Remaining rows:", nrow(vs_data_complete), "\n")

# Refit full model on this complete dataset
full_model_complete <- coxph(as.formula(formula_str), data = vs_data_complete)


# Backward stepwise selection (AIC-based)


library(MASS)

step_model <- step(full_model_complete, direction = "backward", trace = TRUE)

cat("Stepwise selection complete.\n")
summary(step_model)


#  Save final models (full and reduced)


library(broom)

# Save full model
write.csv(
  tidy(full_model_complete, exponentiate = TRUE, conf.int = TRUE),
  "T:/PROFID/Study8/Variable Selection & Model Development/Files/Cox_Multivariate_Full_Model.csv",
  row.names = FALSE
)

# Save stepwise reduced model
write.csv(
  tidy(step_model, exponentiate = TRUE, conf.int = TRUE),
  "T:/PROFID/Study8/Variable Selection & Model Development/Files/Cox_Multivariate_Stepwise_Model.csv",
  row.names = FALSE
)

cat("\n💾 Both full and stepwise Cox models saved successfully.\n")



# ANOVA - Model Comparison (Full vs Stepwise)


# Compare full model and reduced (stepwise) model
anova_result <- anova(full_model_complete, step_model, test = "LRT")

# Display in console
cat("\n📊 ANOVA model comparison (Likelihood Ratio Test):\n")
print(anova_result)

# Save ANOVA result to file
anova_df <- as.data.frame(anova_result)
write.csv(
  anova_df,
  "T:/PROFID/Study8/Variable Selection & Model Development/Files/Cox_ANOVA_Model_Comparison.csv",
  row.names = TRUE
)

cat("\n💾 ANOVA model comparison results saved successfully.\n")



# Remove proxies for highly correlated variables


load("T:/PROFID/Study8/Variable Selection & Model Development/Files/vs_data_complete.RData")

cat("Dataset loaded successfully. Dimensions:\n")
print(dim(vs_data_complete))

cat("\nChecking for multicollinearity among numeric predictors...\n")

numeric_vars <- vs_data_complete %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-Survival_time, -Status)

cor_matrix <- cor(numeric_vars, use = "pairwise.complete.obs")

# Save full correlation matrix to file
write.csv(cor_matrix,
          "T:/PROFID/Study8/Variable Selection & Model Development/Files/Correlation_Matrix.csv",
          row.names = TRUE)

cat("Correlation matrix saved successfully.\n")

print(cor_matrix)

# Identify pairs with high correlation (|r| ≥ 0.8)
high_corr_pairs <- which(abs(cor_matrix) >= 0.8 & abs(cor_matrix) < 1, arr.ind = TRUE)

if (nrow(high_corr_pairs) == 0) {
  cat("No highly correlated numeric variables detected.\n")
} else {
  cat(" Highly correlated variable pairs detected (|r| ≥ 0.8):\n")
  corr_list <- data.frame(
    Var1 = rownames(cor_matrix)[high_corr_pairs[, 1]],
    Var2 = colnames(cor_matrix)[high_corr_pairs[, 2]],
    Correlation = cor_matrix[high_corr_pairs]
  )
  print(corr_list)
  
  # Save list of highly correlated pairs
  write.csv(corr_list,
            "T:/PROFID/Study8/Variable Selection & Model Development/Files/Highly_Correlated_Pairs.csv",
            row.names = FALSE)
}

