# install.packages("dplyr") 
library(dplyr)

# reading files
file1 <- read.csv("ICD.csv")
file2 <- read.csv("NonICD_preserved.csv")
file3 <- read.csv("NonICD_reduced.csv")

# viewing files
View(file1)
View(file2)
View(file3)

# finding the number of columns
ncol(file1)
ncol(file2)
ncol(file3)

# common columns, Because some variables might exist in one file but not the others.
# For fair comparison, we only keep the shared variables that exist in all three datasets.
# This gives us a consistent base for analysis.

common_cols <- Reduce(intersect, list(colnames(file1), colnames(file2), colnames(file3)))


# Subset to those common columns
file1_common <- file1[, common_cols]
file2_common <- file2[, common_cols]
file3_common <- file3[, common_cols]

# label group: 

file1_common$Group <- "ICD"
file2_common$Group <- "NonICD_preserved"
file3_common$Group <- "NonICD_reduced"

# combining all of them together
combined <- bind_rows(file1_common, file2_common, file3_common)

View(combined)



table(combined$Group)


med_iqr_table <- combined %>%
  group_by(Group) %>%
  summarise(across(
    where(is.numeric),
    list(
      median = ~median(., na.rm = TRUE),
      IQR = ~IQR(., na.rm = TRUE)
    ),
    .names = "{.col}_{.fn}"
  ))

View(med_iqr_table)


write.csv(med_iqr_table, 
          "T:/PROFID/Study8/Descriptive Analysis/files/median_IQR_by_group.csv", 
          row.names = FALSE)


#CATEGORICAL VARIABLES: FREQUENCY + PERCENTAGE


# ---Identify categorical variables ---
categorical_vars <- names(combined)[sapply(combined, function(x) is.character(x) | is.factor(x))]

# --- ⃣ Load tidyr for reshaping (install if missing) ---
# install.packages("tidyr")
library(tidyr)

# ---  Summarize categorical variables by group ---
cat_summary <- combined %>%
  select(Group, all_of(categorical_vars)) %>%
  pivot_longer(-Group, names_to = "Variable", values_to = "Category") %>%
  group_by(Group, Variable, Category) %>%
  summarise(n = n(), .groups = "drop_last") %>%
  mutate(percent = round(100 * n / sum(n), 1)) %>%
  arrange(Variable, Group, desc(n))

# View categorical summary
View(cat_summary)

# ---  Save categorical summary ---
write.csv(cat_summary, 
          "T:/PROFID/Study8/Descriptive Analysis/files/categorical_summary.csv",
          row.names = FALSE)
###############################################################

# BETWEEN-GROUP COMPARISONS

# Mann–Whitney U Test (Continuous) + Chi-Square Test (Categorical)

###############################################################



library(dplyr)





#  Mann–Whitney U Test (Continuous Variables)



# Identify numeric (continuous) variables

numeric_vars <- names(combined)[sapply(combined, is.numeric)]



# Run Mann–Whitney U test for each continuous variable between groups

mann_whitney_results <- lapply(numeric_vars, function(var) {
  data_var <- combined %>%
    select(Group, all_of(var)) %>%
    filter(!is.na(.data[[var]]))  # remove NAs
  
  # Perform Kruskal-Wallis test (non-parametric version of ANOVA for >2 groups)
  test <- kruskal.test(reformulate("Group", var), data = data_var)
  data.frame(
    Variable = var,
    p_value = round(test$p.value, 5)
  )
  
}) %>%
  
  bind_rows()


# Save results

write.csv(mann_whitney_results,
          "T:/PROFID/Study8/Descriptive Analysis/files/mann_whitney_results.csv",
          row.names = FALSE)


#  Chi-Square Test (Categorical Variables)


# Identify categorical variables

categorical_vars <- names(combined)[sapply(combined, function(x) is.character(x) | is.factor(x))]

# Run chi-square test for each categorical variable

chi_square_results <- lapply(categorical_vars, function(var) {
  tbl <- table(combined[[var]], combined$Group)
  
  # Only run if table has >1 level per dimension
  
  if (all(dim(tbl) > 1)) {
    test <- chisq.test(tbl)
    p_val <- round(test$p.value, 5)
  } else {
    p_val <- NA
    
  }
  data.frame(Variable = var, p_value = p_val)
  
}) %>%
  
  bind_rows()
# Save results

write.csv(chi_square_results,
          
          "T:/PROFID/Study8/Descriptive Analysis/files/chi_square_results.csv",
          
          row.names = FALSE)


## after heathers feedback
# ----------------------------------------------------------
# LOAD DATA
# ----------------------------------------------------------
df <- read.csv("T:/PROFID/Study8/Variable Selection & Model Development/Files/vs_data_complete.csv")

# Ensure numeric
df$Survival_time <- as.numeric(df$Survival_time)

# ----------------------------------------------------------
# FUNCTION TO SUMMARISE FOLLOW-UP TIME
# ----------------------------------------------------------
summarise_followup <- function(data, group_name){
  out <- data.frame(
    Group = group_name,
    N = nrow(data),
    Mean_Followup_months   = mean(data$Survival_time, na.rm = TRUE),
    Median_Followup_months = median(data$Survival_time, na.rm = TRUE),
    Min_Followup_months    = min(data$Survival_time, na.rm = TRUE),
    Max_Followup_months    = max(data$Survival_time, na.rm = TRUE),
    IQR_25 = quantile(data$Survival_time, 0.25, na.rm = TRUE),
    IQR_75 = quantile(data$Survival_time, 0.75, na.rm = TRUE)
  )
  return(out)
}

# ----------------------------------------------------------
# SPLIT BY GROUP
# ----------------------------------------------------------
df_ICD        <- df[df$Group == "ICD", ]
df_preserved  <- df[df$Group == "NonICD_preserved", ]
df_reduced    <- df[df$Group == "NonICD_reduced", ]

# ----------------------------------------------------------
# APPLY FUNCTION
# ----------------------------------------------------------
res_all       <- summarise_followup(df, "All_Cohorts")
res_ICD       <- summarise_followup(df_ICD, "ICD")
res_preserved <- summarise_followup(df_preserved, "NonICD_preserved")
res_reduced   <- summarise_followup(df_reduced, "NonICD_reduced")

# ----------------------------------------------------------
# COMBINE RESULTS
# ----------------------------------------------------------
followup_results <- rbind(res_ICD, res_preserved, res_reduced, res_all)

# ----------------------------------------------------------
# SAVE OUTPUT
# ----------------------------------------------------------
outdir <- "T:/PROFID/Study8/Descriptive Analysis/files"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write.csv(followup_results,
          file.path(outdir, "Followup_Time_Summary.csv"),
          row.names = FALSE)

print(followup_results)

## for all variables 

## ----------------------------------------------------------
## LOAD PACKAGES
## ----------------------------------------------------------
library(dplyr)
library(tidyr)

## ----------------------------------------------------------
## LOAD DATA
## ----------------------------------------------------------
df <- read.csv("T:/PROFID/Study8/Variable Selection & Model Development/Files/vs_data_complete.csv")

# Make sure Group exists
stopifnot("Group" %in% names(df))

## ----------------------------------------------------------
## IDENTIFY NUMERIC VARIABLES
## ----------------------------------------------------------
numeric_vars <- names(df)[sapply(df, is.numeric)]

# Optional: exclude variables you do NOT want in the summary
# (e.g. Survival_time, IDs, etc.)
# numeric_vars <- setdiff(numeric_vars, c("Survival_time", "some_ID_var"))

## ----------------------------------------------------------
## 1) SUMMARY FOR ALL COHORTS COMBINED (NO GROUP)
## ----------------------------------------------------------
all_variables_overall <- df %>%
  pivot_longer(cols = all_of(numeric_vars),
               names_to = "Variable",
               values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(
    N      = sum(!is.na(Value)),
    Median = median(Value, na.rm = TRUE),
    Q1     = quantile(Value, 0.25, na.rm = TRUE),
    Q3     = quantile(Value, 0.75, na.rm = TRUE),
    IQR    = IQR(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Group = "All_Cohorts") %>%
  select(Group, everything())

## ----------------------------------------------------------
## 2) SUMMARY BY GROUP (ICD / NonICD_preserved / NonICD_reduced)
## ----------------------------------------------------------
all_variables_by_group <- df %>%
  pivot_longer(cols = all_of(numeric_vars),
               names_to = "Variable",
               values_to = "Value") %>%
  group_by(Group, Variable) %>%
  summarise(
    N      = sum(!is.na(Value)),
    Median = median(Value, na.rm = TRUE),
    Q1     = quantile(Value, 0.25, na.rm = TRUE),
    Q3     = quantile(Value, 0.75, na.rm = TRUE),
    IQR    = IQR(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  select(Group, Variable, N, Median, Q1, Q3, IQR)

## ----------------------------------------------------------
## 3) COMBINE INTO ONE OBJECT: all_variables
## ----------------------------------------------------------
all_variables <- bind_rows(all_variables_by_group, all_variables_overall)

## ----------------------------------------------------------
## 4) SAVE TO CSV (EXPLICIT FILE PATH)
## ----------------------------------------------------------
outdir <- "T:/PROFID/Study8/Descriptive Analysis/files"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

write.csv(
  all_variables,
  file.path(outdir, "all_variables_IQR_median_by_group_and_overall.csv"),
  row.names = FALSE
)

## Quick look in R
print(all_variables)

