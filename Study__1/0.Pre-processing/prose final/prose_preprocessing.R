rm(list = ls())

# Requirements
# ==============================================================================
library(data.table)
library(gt)
library(openxlsx)

source("./utils/hz-basic-summary-statistics.R")

# CONFIGURATION
# ==============================================================================
input_file <- "../../a-datastack/output/datastack-table1.rds"
grouping_variable <- "Status_FIS"

# Define continuous variables (you need to specify these)
continuous_vars <- c(
  "Age", "BMI", "SBP", "DBP", "LVEF", "LVDD",
  "BUN", "Cholesterol", "CRP", "eGFR", "Haemoglobin",
  "HbA1c", "HDL", "IL6", "LDL", "NTProBNP",
  "Potassium", "Sodium", "Triglycerides", "Troponin_T",
  "TSH", "HR", "PR", "QRS", "QTc"
)

# Define variables to EXCLUDE from table
exclude_table_vars <- c(
  "Status",
  "Survival_time",
  "CVD_risk_region",
  "Endpoint_type",
  "Baseline_type",
  "Time_zero_Ym",
  "Time_zero_Y",
  "Time_FIS_days",
  "V1"
)

# Load Data
# ==============================================================================
#ds <- readRDS(input_file)
ds[, group := get(grouping_variable)]

cat("Group distribution:\n")
print(table(ds$group, useNA = "always"))

ds_no <- ds[group == 0]
ds_yes <- ds[group == 1]

# Filter continuous variables to those in dataset
continuous_vars <- continuous_vars[continuous_vars %in% names(ds)]

# Auto-detect binary variables
all_vars <- setdiff(names(ds), 
                    c("ID", "DB", "group", grouping_variable,
                      grep("^bin_", names(ds), value = TRUE)))

binary_vars <- setdiff(all_vars, c("N", continuous_vars, exclude_table_vars))

cat("\nContinuous variables:", length(continuous_vars), "\n")
cat("Binary variables:", length(binary_vars), "\n")
cat("Excluded variables:", length(exclude_table_vars), "\n")

# STEP 1: Assess distributions (skewness-based)
# ==============================================================================
cat("\nAssessing distributions for continuous variables...\n")
cat("Decision rule: |skewness| ≤ 0.5 → MEAN, |skewness| > 0.5 → MEDIAN\n\n")

median_vars <- c()
mean_vars <- c()

for (var in continuous_vars) {
  x <- ds[[var]][!is.na(ds[[var]])]
  
  if (length(x) < 3) {
    median_vars <- c(median_vars, var)
    cat("  ", var, ": MEDIAN (insufficient data)\n", sep = "")
    next
  }
  
  # Check for -Inf or Inf values
  if (any(is.infinite(x))) {
    cat("  ", var, ": WARNING - Contains -Inf or Inf values, removing them\n", sep = "")
    x <- x[is.finite(x)]
    # Update in dataset
    ds[is.infinite(get(var)), (var) := NA]
    ds_no[is.infinite(get(var)), (var) := NA]
    ds_yes[is.infinite(get(var)), (var) := NA]
  }
  
  if (length(x) < 3) {
    median_vars <- c(median_vars, var)
    next
  }
  
  # Calculate skewness: (mean - median) / SD
  skew <- (mean(x) - median(x)) / sd(x)
  
  # Decision: |skewness| ≤ 0.5 → mean, else → median
  if (abs(skew) <= 0.5) {
    mean_vars <- c(mean_vars, var)
    cat("  ", var, ": MEAN (skewness = ", round(skew, 3), ")\n", sep = "")
  } else {
    median_vars <- c(median_vars, var)
    cat("  ", var, ": MEDIAN (skewness = ", round(skew, 3), ")\n", sep = "")
  }
}

cat("\nSummary:\n")
cat("  Variables using MEAN:", length(mean_vars), "\n")
cat("  Variables using MEDIAN:", length(median_vars), "\n")

# STEP 2: Generate summaries using HZ functions
# ==============================================================================
cat("\nGenerating descriptive statistics...\n")

ds_sm <- HZ.eda_describe_data(ds)$summary.sdc
ds_no_sm <- HZ.eda_describe_data(ds_no)$summary.sdc
ds_yes_sm <- HZ.eda_describe_data(ds_yes)$summary.sdc

# Exclude unwanted variables
exclude_vars <- c("ID", "DB", "group", grouping_variable, 
                  exclude_table_vars,
                  grep("^bin_", names(ds), value = TRUE))

ds_sm <- ds_sm[!(Variable %in% exclude_vars)]
ds_no_sm <- ds_no_sm[!(Variable %in% exclude_vars)]
ds_yes_sm <- ds_yes_sm[!(Variable %in% exclude_vars)]

# Merge summaries
table1 <- merge(ds_sm, ds_no_sm, by = c("Variable", "Category"),
                suffixes = c("_total", "_no"), sort = FALSE)
table1 <- merge(table1, ds_yes_sm[, .(Variable, Category, Summary_yes = Summary)],
                by = c("Variable", "Category"), sort = FALSE)

# Calculate missing data CORRECTLY from raw dataset
# ==============================================================================
cat("\nCalculating missing data properly...\n")

all_table_vars <- unique(c("N", continuous_vars, binary_vars))
all_table_vars <- all_table_vars[all_table_vars %in% names(ds)]

missing_data <- data.table(
  Variable = all_table_vars,
  Missing_total = NA_character_
)

for (var in all_table_vars) {
  if (var == "N") {
    missing_data[Variable == var, Missing_total := "0 (0.0%)"]
  } else {
    n_missing <- sum(is.na(ds[[var]]))
    n_total <- nrow(ds)
    pct_missing <- 100 * n_missing / n_total
    missing_data[Variable == var, Missing_total := sprintf("%d (%.1f%%)", n_missing, pct_missing)]
    cat("  ", var, ": ", n_missing, " (", round(pct_missing, 1), "%)\n", sep = "")
  }
}

# STEP 3: Filter based on distribution assessment
# ==============================================================================
cat("\nFiltering summary statistics based on distribution...\n")

# Remove MEDIAN rows for variables that should show MEAN
table1 <- table1[!(Variable %in% mean_vars & grepl("Median", Category))]

# Remove MEAN rows for variables that should show MEDIAN
table1 <- table1[!(Variable %in% median_vars & grepl("Mean", Category))]

# Remove [Min, Max] rows
table1 <- table1[Category != "[Min, Max]"]

# Remove Missing rows (we calculated them properly above)
table1 <- table1[!grepl("Missing", Category)]

# Handle categorical variables properly - keep only ONE category
cat("\nFiltering categorical variables to one category each...\n")

# For Sex: keep only "Male", remove "Female"
if ("Sex" %in% table1$Variable) {
  table1 <- table1[!(Variable == "Sex" & Category == "Female")]
  cat("  Sex: Kept 'Male', removed 'Female'\n")
}

# For all other binary variables
for (var in binary_vars) {
  if (var == "Sex") next
  
  var_categories <- unique(table1[Variable == var, Category])
  
  cat("  ", var, ": categories = ", paste(var_categories, collapse = ", "), sep = "")
  
  if ("Yes" %in% var_categories & "No" %in% var_categories) {
    table1 <- table1[!(Variable == var & Category == "No")]
    cat(" → Kept 'Yes'\n")
  }
  else if ("0" %in% var_categories & "1" %in% var_categories) {
    table1 <- table1[!(Variable == var & Category == "0")]
    cat(" → Kept '1'\n")
  }
  else if (length(var_categories) > 1) {
    keep_category <- var_categories[1]
    remove_categories <- setdiff(var_categories, keep_category)
    table1 <- table1[!(Variable == var & Category %in% remove_categories)]
    cat(" → Kept first category\n")
  }
  else {
    cat(" → Only 1 category\n")
  }
}

cat("  Kept MEAN for", length(mean_vars), "variables\n")
cat("  Kept MEDIAN for", length(median_vars), "variables\n")

# FIX: Statistical tests - DO THIS BEFORE collapsing to one row
# ==============================================================================
cat("\nRunning statistical tests...\n")

# Create a separate data.table to store p-values
pvalue_results <- data.table(
  Variable = character(),
  p_value = numeric(),
  test_method = character()
)

all_test_vars <- unique(table1[Variable != "N", Variable])

for (var in all_test_vars) {
  
  # Determine which version to use for testing
  if (var %in% binary_vars) {
    test_var <- paste0("bin_", var)
    if (!test_var %in% names(ds)) test_var <- var
  } else {
    test_var <- var
  }
  
  x1 <- ds_no[[test_var]][!is.na(ds_no[[test_var]])]
  x2 <- ds_yes[[test_var]][!is.na(ds_yes[[test_var]])]
  
  if (length(x1) < 2 | length(x2) < 2) {
    cat("  ", var, ": Skipped (insufficient data)\n", sep = "")
    next
  }
  
  # Test
  if (var %in% continuous_vars) {
    if (var %in% mean_vars) {
      result <- t.test(x1, x2)
      test_name <- "t-test"
    } else {
      result <- wilcox.test(x1, x2, exact = FALSE)
      test_name <- "Mann-Whitney U"
    }
    p_val <- result$p.value
  } else {
    tbl <- table(c(x1, x2), c(rep(0, length(x1)), rep(1, length(x2))))
    
    # Check if valid table
    if (nrow(tbl) < 2 | ncol(tbl) < 2) {
      cat("  ", var, ": Skipped (invalid table)\n", sep = "")
      next
    }
    
    if (min(tbl) < 5) {
      result <- tryCatch(
        fisher.test(tbl, simulate.p.value = TRUE, B = 10000),
        error = function(e) {
          cat("  ", var, ": Fisher test failed, using Chi-square\n", sep = "")
          chisq.test(tbl)
        }
      )
      test_name <- "Fisher's exact"
    } else {
      result <- chisq.test(tbl)
      test_name <- "Chi-square"
    }
    p_val <- result$p.value
  }
  
  # Store results
  pvalue_results <- rbind(pvalue_results, 
                          data.table(Variable = var, 
                                     p_value = p_val, 
                                     test_method = test_name))
  
  cat("  ", var, ": ", test_name, " (p = ", round(p_val, 4), ")\n", sep = "")
}

# STEP 4: Format and clean table
# ==============================================================================
cat("\nFormatting table...\n")

# Rename columns FIRST
setnames(table1, "Summary_total", "Total")
setnames(table1, "Summary_no", "No Inappropriate Therapy")
setnames(table1, "Summary_yes", "Inappropriate Therapy")

# Now collapse to one row per variable
table1_collapsed <- table1[, .(
  Category = Category[1],
  Total = Total[1],
  `No Inappropriate Therapy` = `No Inappropriate Therapy`[1],
  `Inappropriate Therapy` = `Inappropriate Therapy`[1]
), by = Variable]

# Merge with p-values
table1_collapsed <- merge(table1_collapsed, pvalue_results, 
                          by = "Variable", all.x = TRUE)

# Merge with missing data
table1_collapsed <- merge(table1_collapsed, missing_data, 
                          by = "Variable", all.x = TRUE)

# Create display variable
table1_collapsed[, Variable_display := Variable]

# Clean up category labels
table1_collapsed[Category == "Yes", Category := ""]
table1_collapsed[Category == "Male", Category := ""]
table1_collapsed[Category == "1", Category := ""]
table1_collapsed[grepl("Mean \\(SD\\)", Category), Category := ""]
table1_collapsed[grepl("Median \\(IQR\\)", Category), Category := ""]

# Order variables
var_order <- c("N", continuous_vars, binary_vars)
var_order <- var_order[var_order %in% table1_collapsed$Variable]
table1_collapsed[, Variable := factor(Variable, levels = var_order)]
table1_collapsed <- table1_collapsed[order(Variable)]

# Format p-values
table1_collapsed[, P_value := ifelse(is.na(p_value), "",
                                     ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)))]

# Create final table
final <- table1_collapsed[, .(
  Characteristic = Variable_display,
  Total,
  `No Inappropriate Therapy`,
  `Inappropriate Therapy`,
  `Missing (Total)` = Missing_total,
  `P value` = P_value
)]

# STEP 5: Save outputs
# ==============================================================================
cat("\nSaving outputs...\n")

fwrite(final, "../tables/table-1-baseline.csv", na = "")

# GT table
gt_table <- final %>%
  gt() %>%
  tab_header(title = "Table 1. Baseline Characteristics by Inappropriate ICD Therapy") %>%
  cols_label(
    Characteristic = "Characteristic",
    Total = "Total",
    `No Inappropriate Therapy` = "No Inappropriate Therapy",
    `Inappropriate Therapy` = "Inappropriate Therapy",
    `Missing (Total)` = "Missing (Total)",
    `P value` = md("*P* value")
  ) %>%
  cols_align(align = "left", columns = Characteristic) %>%
  cols_align(align = "center", columns = c(Total, `No Inappropriate Therapy`, 
                                           `Inappropriate Therapy`, 
                                           `Missing (Total)`,
                                           `P value`)) %>%
  cols_width(
    Characteristic ~ px(200),
    Total ~ px(120),
    `No Inappropriate Therapy` ~ px(140),
    `Inappropriate Therapy` ~ px(140),
    `Missing (Total)` ~ px(120),
    `P value` ~ px(100)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Characteristic)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `P value`,
      rows = `P value` != "" & as.numeric(gsub("<", "", `P value`)) < 0.05
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(2)),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_borders(sides = "top", color = "black", weight = px(2)),
    locations = cells_column_labels()
  ) %>%
  tab_footnote(
    footnote = md(paste(
      "Data are n (%) for categorical variables, mean ± SD for approximately normally",
      "distributed continuous variables, and median (IQR) for skewed continuous variables.",
      "Distribution assessed using skewness (|skewness| ≤ 0.5 indicates approximately normal).",
      "Statistical tests: *t*-test for normal continuous variables,",
      "Mann-Whitney *U* test for skewed continuous variables,",
      "χ² test for categorical variables (expected ≥5),",
      "Fisher's exact test for categorical variables (expected <5).",
      "ICD = implantable cardioverter-defibrillator; IQR = interquartile range; SD = standard deviation.",
      sep = " "
    ))
  ) %>%
  tab_options(
    table.font.size = px(12),
    table.width = pct(100),
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold",
    column_labels.font.weight = "bold",
    footnotes.font.size = px(10)
  )
print(gt_table)

gtsave(gt_table, "../tables/table-1-baseline.html")
gtsave(gt_table, "../tables/table-1-baseline.docx")

cat("  Saved: ../tables/table-1-baseline.html\n")
cat("  Saved: ../tables/table-1-baseline.docx\n")

# Excel
wb <- createWorkbook()
addWorksheet(wb, "Table 1")

writeData(wb, 1, "Table 1. Baseline Characteristics by Inappropriate ICD Therapy", 
          startRow = 1, colNames = FALSE)
addStyle(wb, 1, createStyle(fontSize = 14, textDecoration = "bold"), rows = 1, cols = 1)

writeData(wb, 1, final, startRow = 3)

headerStyle <- createStyle(
  fontSize = 11,
  fontName = "Times New Roman",
  textDecoration = "bold",
  halign = "center",
  border = "Bottom",
  borderColour = "black"
)
addStyle(wb, 1, headerStyle, rows = 3, cols = 1:ncol(final), gridExpand = TRUE)

var_rows <- 4:(nrow(final) + 3)
addStyle(wb, 1, createStyle(fontName = "Times New Roman", textDecoration = "bold"),
         rows = var_rows, cols = 1, gridExpand = TRUE, stack = TRUE)

centerStyle <- createStyle(fontName = "Times New Roman", halign = "center")
addStyle(wb, 1, centerStyle, rows = 4:(nrow(final) + 3), cols = 2:ncol(final), 
         gridExpand = TRUE, stack = TRUE)

sig_rows <- which(final$`P value` != "" & 
                    as.numeric(gsub("<", "", final$`P value`)) < 0.05) + 3
if (length(sig_rows) > 0) {
  addStyle(wb, 1, createStyle(fontName = "Times New Roman", textDecoration = "bold"),
           rows = sig_rows, cols = which(names(final) == "P value"),
           gridExpand = TRUE, stack = TRUE)
}

setColWidths(wb, 1, cols = 1, widths = 25)
setColWidths(wb, 1, cols = 2:6, widths = 15)

saveWorkbook(wb, "../tables/table-1-baseline.xlsx", overwrite = TRUE)
cat("  Saved: ../tables/table-1-baseline.xlsx\n")

# Summary
# ==============================================================================
cat("\n========================================\n")
cat("TABLE 1 COMPLETE\n")
cat("========================================\n\n")

cat("SAMPLE SIZE:\n")
cat("  Total: N =", final[Characteristic == "N", Total], "\n")
cat("  No inappropriate therapy:", 
    final[Characteristic == "N", `No Inappropriate Therapy`], "\n")
cat("  Inappropriate therapy:", 
    final[Characteristic == "N", `Inappropriate Therapy`], "\n\n")

cat("VARIABLES:\n")
cat("  Continuous (MEAN):", length(mean_vars), "variables\n")
cat("  Continuous (MEDIAN):", length(median_vars), "variables\n")
cat("  Binary/Categorical:", length(binary_vars), "variables\n")
cat("  Total rows in table:", nrow(final), "\n\n")

n_sig <- sum(!is.na(final$`P value`) & final$`P value` != "" &
               as.numeric(gsub("<", "", final$`P value`)) < 0.05)
cat("STATISTICAL TESTS:\n")
cat("  Significant differences (P < 0.05):", n_sig, "\n\n")

cat("OUTPUT FILES:\n")
cat("  1. ../tables/table-1-baseline.csv\n")
cat("  2. ../tables/table-1-baseline.html\n")
cat("  3. ../tables/table-1-baseline.docx\n")
cat("  4. ../tables/table-1-baseline.xlsx\n\n")

cat("========================================\n\n")

print(final)


rm(list = ls())



# Requirements

# ==============================================================================

library(data.table)

library(gt)

library(openxlsx)



source("./utils/hz-basic-summary-statistics.R")



# CONFIGURATION

# ==============================================================================

input_file <- "../../a-datastack/output/datastack-table1.rds"

grouping_variable <- "status_FIS"



# Define continuous variables (you need to specify these)

continuous_vars <- c(
  
  "Age", "BMI", "SBP", "DBP", "LVEF", "LVDD",
  
  "BUN", "Cholesterol", "CRP", "eGFR", "Haemoglobin",
  
  "HbA1c", "HDL", "IL6", "LDL", "NTProBNP",
  
  "Potassium", "Sodium", "Triglycerides", "Troponin_T",
  
  "TSH", "HR", "PR", "QRS", "QTc"
  
)



# Define variables to EXCLUDE from table
exclude_table_vars <- c(
  "Status",
  "Survival_time",
  "CVD_risk_region",
  "Endpoint_type",
  "Baseline_type",
  "Time_zero_Ym",
  "Time_zero_Y",
  "Time_FIS_days",
  "V1"
)

  
  
  
  
  




# Load Data

# ==============================================================================

ds <- readRDS(input_file)

ds[, group := get(grouping_variable)]



cat("Group distribution:\n")

print(table(ds$group, useNA = "always"))



ds_no <- ds[group == 0]

ds_yes <- ds[group == 1]



# Filter continuous variables to those in dataset

continuous_vars <- continuous_vars[continuous_vars %in% names(ds)]



# Auto-detect binary variables

all_vars <- setdiff(names(ds), 
                    
                    c("ID", "DB", "group", grouping_variable,
                      
                      grep("^bin_", names(ds), value = TRUE)))



binary_vars <- setdiff(all_vars, c("N", continuous_vars, exclude_table_vars))



cat("\nContinuous variables:", length(continuous_vars), "\n")

cat("Binary variables:", length(binary_vars), "\n")

cat("Excluded variables:", length(exclude_table_vars), "\n")



# FIX: Clean data - remove variables with too many missing/infinite values

# ==============================================================================

cat("\nChecking data quality...\n")



vars_to_remove <- c()



# Check continuous variables

for (var in continuous_vars) {
  
  x <- ds[[var]]
  
  n_total <- length(x)
  
  n_missing <- sum(is.na(x))
  
  n_infinite <- sum(is.infinite(x), na.rm = TRUE)
  
  n_valid <- n_total - n_missing - n_infinite
  
  pct_valid <- 100 * n_valid / n_total
  
  
  
  if (pct_valid < 5) {  # Less than 5% valid data
    
    cat("  ", var, ": EXCLUDED - only ", round(pct_valid, 1), "% valid data\n", sep = "")
    
    vars_to_remove <- c(vars_to_remove, var)
    
  } else if (n_infinite > 0) {
    
    cat("  ", var, ": Removing ", n_infinite, " infinite values\n", sep = "")
    
    ds[is.infinite(get(var)), (var) := NA]
    
    ds_no[is.infinite(get(var)), (var) := NA]
    
    ds_yes[is.infinite(get(var)), (var) := NA]
    
  }
  
}



# Check binary variables

for (var in binary_vars) {
  
  x <- ds[[var]]
  
  n_total <- length(x)
  
  n_missing <- sum(is.na(x))
  
  n_valid <- n_total - n_missing
  
  pct_valid <- 100 * n_valid / n_total
  
  
  
  if (pct_valid < 5) {
    
    cat("  ", var, ": EXCLUDED - only ", round(pct_valid, 1), "% valid data\n", sep = "")
    
    vars_to_remove <- c(vars_to_remove, var)
    
  }
  
}



# Remove problematic variables

if (length(vars_to_remove) > 0) {
  
  cat("\nRemoving", length(vars_to_remove), "variables with insufficient data:\n")
  
  cat("  ", paste(vars_to_remove, collapse = ", "), "\n")
  
  
  
  continuous_vars <- setdiff(continuous_vars, vars_to_remove)
  
  binary_vars <- setdiff(binary_vars, vars_to_remove)
  
  exclude_table_vars <- c(exclude_table_vars, vars_to_remove)
  
}



# STEP 1: Assess distributions (skewness-based)

# ==============================================================================

cat("\nAssessing distributions for continuous variables...\n")

cat("Decision rule: |skewness| ≤ 0.5 → MEAN, |skewness| > 0.5 → MEDIAN\n\n")



median_vars <- c()

mean_vars <- c()



for (var in continuous_vars) {
  
  x <- ds[[var]][!is.na(ds[[var]]) & is.finite(ds[[var]])]
  
  
  
  if (length(x) < 3) {
    
    median_vars <- c(median_vars, var)
    
    cat("  ", var, ": MEDIAN (insufficient data)\n", sep = "")
    
    next
    
  }
  
  
  
  # Calculate skewness: (mean - median) / SD
  
  skew <- (mean(x) - median(x)) / sd(x)
  
  
  
  # Decision: |skewness| ≤ 0.5 → mean, else → median
  
  if (abs(skew) <= 0.5) {
    
    mean_vars <- c(mean_vars, var)
    
    cat("  ", var, ": MEAN (skewness = ", round(skew, 3), ")\n", sep = "")
    
  } else {
    
    median_vars <- c(median_vars, var)
    
    cat("  ", var, ": MEDIAN (skewness = ", round(skew, 3), ")\n", sep = "")
    
  }
  
}



cat("\nSummary:\n")

cat("  Variables using MEAN:", length(mean_vars), "\n")

cat("  Variables using MEDIAN:", length(median_vars), "\n")



# STEP 2: Generate summaries using HZ functions

# ==============================================================================

cat("\nGenerating descriptive statistics...\n")



ds_sm <- HZ.eda_describe_data(ds)$summary.sdc

ds_no_sm <- HZ.eda_describe_data(ds_no)$summary.sdc

ds_yes_sm <- HZ.eda_describe_data(ds_yes)$summary.sdc



# Exclude unwanted variables

exclude_vars <- c("ID", "DB", "group", grouping_variable, 
                  
                  exclude_table_vars,
                  
                  grep("^bin_", names(ds), value = TRUE))



ds_sm <- ds_sm[!(Variable %in% exclude_vars)]

ds_no_sm <- ds_no_sm[!(Variable %in% exclude_vars)]

ds_yes_sm <- ds_yes_sm[!(Variable %in% exclude_vars)]



# Merge summaries

table1 <- merge(ds_sm, ds_no_sm, by = c("Variable", "Category"),
                
                suffixes = c("_total", "_no"), sort = FALSE)

table1 <- merge(table1, ds_yes_sm[, .(Variable, Category, Summary_yes = Summary)],
                
                by = c("Variable", "Category"), sort = FALSE)



# Calculate missing data CORRECTLY from raw dataset

# ==============================================================================

cat("\nCalculating missing data properly...\n")



all_table_vars <- unique(c("N", continuous_vars, binary_vars))

all_table_vars <- all_table_vars[all_table_vars %in% names(ds)]



missing_data <- data.table(
  
  Variable = all_table_vars,
  
  Missing_total = NA_character_
  
)



for (var in all_table_vars) {
  
  if (var == "N") {
    
    missing_data[Variable == var, Missing_total := "0 (0.0%)"]
    
  } else {
    
    n_missing <- sum(is.na(ds[[var]]))
    
    n_total <- nrow(ds)
    
    pct_missing <- 100 * n_missing / n_total
    
    missing_data[Variable == var, Missing_total := sprintf("%d (%.1f%%)", n_missing, pct_missing)]
    
  }
  
}



# STEP 3: Filter based on distribution assessment

# ==============================================================================

cat("\nFiltering summary statistics based on distribution...\n")



# Remove MEDIAN rows for variables that should show MEAN

table1 <- table1[!(Variable %in% mean_vars & grepl("Median", Category))]



# Remove MEAN rows for variables that should show MEDIAN

table1 <- table1[!(Variable %in% median_vars & grepl("Mean", Category))]



# Remove [Min, Max] rows

table1 <- table1[Category != "[Min, Max]"]



# Remove Missing rows (we calculated them properly above)

table1 <- table1[!grepl("Missing", Category)]



# FIX: Remove rows with "Inf to -Inf" or similar weird values

table1 <- table1[!grepl("Inf|NaN", Summary_total, ignore.case = TRUE)]

table1 <- table1[!grepl("Inf|NaN", Summary_no, ignore.case = TRUE)]

table1 <- table1[!grepl("Inf|NaN", Summary_yes, ignore.case = TRUE)]



# Handle categorical variables properly - keep only ONE category

cat("\nFiltering categorical variables to one category each...\n")



# For Sex: keep only "Male", remove "Female"

if ("Sex" %in% table1$Variable) {
  
  table1 <- table1[!(Variable == "Sex" & Category == "Female")]
  
  cat("  Sex: Kept 'Male', removed 'Female'\n")
  
}



# For all other binary variables

for (var in binary_vars) {
  
  if (var == "Sex") next
  
  
  
  var_categories <- unique(table1[Variable == var, Category])
  
  
  
  if (length(var_categories) == 0) next
  
  
  
  cat("  ", var, ": categories = ", paste(var_categories, collapse = ", "), sep = "")
  
  
  
  if ("Yes" %in% var_categories & "No" %in% var_categories) {
    
    table1 <- table1[!(Variable == var & Category == "No")]
    
    cat(" → Kept 'Yes'\n")
    
  }
  
  else if ("0" %in% var_categories & "1" %in% var_categories) {
    
    table1 <- table1[!(Variable == var & Category == "0")]
    
    cat(" → Kept '1'\n")
    
  }
  
  else if (length(var_categories) > 1) {
    
    keep_category <- var_categories[1]
    
    remove_categories <- setdiff(var_categories, keep_category)
    
    table1 <- table1[!(Variable == var & Category %in% remove_categories)]
    
    cat(" → Kept first category\n")
    
  }
  
  else {
    
    cat(" → Only 1 category\n")
    
  }
  
}



cat("  Kept MEAN for", length(mean_vars), "variables\n")

cat("  Kept MEDIAN for", length(median_vars), "variables\n")



# Statistical tests - DO THIS BEFORE collapsing to one row

# ==============================================================================

cat("\nRunning statistical tests...\n")



# Create a separate data.table to store p-values

pvalue_results <- data.table(
  
  Variable = character(),
  
  p_value = numeric(),
  
  test_method = character()
  
)



all_test_vars <- unique(table1[Variable != "N", Variable])



for (var in all_test_vars) {
  
  
  
  # Determine which version to use for testing
  
  if (var %in% binary_vars) {
    
    test_var <- paste0("bin_", var)
    
    if (!test_var %in% names(ds)) test_var <- var
    
  } else {
    
    test_var <- var
    
  }
  
  
  
  x1 <- ds_no[[test_var]][!is.na(ds_no[[test_var]]) & is.finite(ds_no[[test_var]])]
  
  x2 <- ds_yes[[test_var]][!is.na(ds_yes[[test_var]]) & is.finite(ds_yes[[test_var]])]
  
  
  
  if (length(x1) < 2 | length(x2) < 2) {
    
    cat("  ", var, ": Skipped (insufficient data)\n", sep = "")
    
    next
    
  }
  
  
  
  # Test
  
  if (var %in% continuous_vars) {
    
    if (var %in% mean_vars) {
      
      result <- t.test(x1, x2)
      
      test_name <- "t-test"
      
    } else {
      
      result <- wilcox.test(x1, x2, exact = FALSE)
      
      test_name <- "Mann-Whitney U"
      
    }
    
    p_val <- result$p.value
    
  } else {
    
    tbl <- table(c(x1, x2), c(rep(0, length(x1)), rep(1, length(x2))))
    
    
    
    if (nrow(tbl) < 2 | ncol(tbl) < 2) {
      
      cat("  ", var, ": Skipped (invalid table)\n", sep = "")
      
      next
      
    }
    
    
    
    if (min(tbl) < 5) {
      
      result <- tryCatch(
        
        fisher.test(tbl, simulate.p.value = TRUE, B = 10000),
        
        error = function(e) {
          
          cat("  ", var, ": Fisher test failed, using Chi-square\n", sep = "")
          
          chisq.test(tbl)
          
        }
        
      )
      
      test_name <- "Fisher's exact"
      
    } else {
      
      result <- chisq.test(tbl)
      
      test_name <- "Chi-square"
      
    }
    
    p_val <- result$p.value
    
  }
  
  
  
  # Store results
  
  pvalue_results <- rbind(pvalue_results, 
                          
                          data.table(Variable = var, 
                                     
                                     p_value = p_val, 
                                     
                                     test_method = test_name))
  
  
  
  cat("  ", var, ": ", test_name, " (p = ", round(p_val, 4), ")\n", sep = "")
  
}



# Format and clean table

# ==============================================================================

cat("\nFormatting table...\n")



# Rename columns FIRST

setnames(table1, "Summary_total", "Total")

setnames(table1, "Summary_no", "No Inappropriate Therapy")

setnames(table1, "Summary_yes", "Inappropriate Therapy")



# Now collapse to one row per variable

table1_collapsed <- table1[, .(
  
  Category = Category[1],
  
  Total = Total[1],
  
  `No Inappropriate Therapy` = `No Inappropriate Therapy`[1],
  
  `Inappropriate Therapy` = `Inappropriate Therapy`[1]
  
), by = Variable]



# Merge with p-values

table1_collapsed <- merge(table1_collapsed, pvalue_results, 
                          
                          by = "Variable", all.x = TRUE)



# Merge with missing data

table1_collapsed <- merge(table1_collapsed, missing_data, 
                          
                          by = "Variable", all.x = TRUE)



# Create display variable

table1_collapsed[, Variable_display := Variable]



# Clean up category labels

table1_collapsed[Category == "Yes", Category := ""]

table1_collapsed[Category == "Male", Category := ""]

table1_collapsed[Category == "1", Category := ""]

table1_collapsed[grepl("Mean \\(SD\\)", Category), Category := ""]

table1_collapsed[grepl("Median \\(IQR\\)", Category), Category := ""]



# Order variables

var_order <- c("N", continuous_vars, binary_vars)

var_order <- var_order[var_order %in% table1_collapsed$Variable]

table1_collapsed[, Variable := factor(Variable, levels = var_order)]

table1_collapsed <- table1_collapsed[order(Variable)]



# Format p-values

table1_collapsed[, P_value := ifelse(is.na(p_value), "",
                                     
                                     ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)))]



# Create final table

final <- table1_collapsed

# Create final table
final <- table1_collapsed[, .(
  Characteristic = Variable_display,
  Total,
  `No Inappropriate Therapy`,
  `Inappropriate Therapy`,
  `Missing (Total)` = Missing_total,
  `P value` = P_value
)]

# Save outputs
# ==============================================================================
cat("\nSaving outputs...\n")

fwrite(final, "../tables/table-1-baseline.csv", na = "")

# GT table
gt_table <- final %>%
  gt() %>%
  tab_header(title = "Table 1. Baseline Characteristics by Inappropriate ICD Therapy") %>%
  cols_label(
    Characteristic = "Characteristic",
    Total = "Total",
    `No Inappropriate Therapy` = "No Inappropriate Therapy",
    `Inappropriate Therapy` = "Inappropriate Therapy",
    `Missing (Total)` = "Missing (Total)",
    `P value` = md("*P* value")
  ) %>%
  cols_align(align = "left", columns = Characteristic) %>%
  cols_align(align = "center", columns = c(Total, `No Inappropriate Therapy`, 
                                           `Inappropriate Therapy`, 
                                           `Missing (Total)`,
                                           `P value`)) %>%
  cols_width(
    Characteristic ~ px(200),
    Total ~ px(120),
    `No Inappropriate Therapy` ~ px(140),
    `Inappropriate Therapy` ~ px(140),
    `Missing (Total)` ~ px(120),
    `P value` ~ px(100)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(columns = Characteristic)
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_body(
      columns = `P value`,
      rows = `P value` != "" & as.numeric(gsub("<", "", `P value`)) < 0.05
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_borders(sides = "bottom", color = "black", weight = px(2)),
    locations = cells_column_labels()
  ) %>%
  tab_style(
    style = cell_borders(sides = "top", color = "black", weight = px(2)),
    locations = cells_column_labels()
  ) %>%
  tab_footnote(
    footnote = md(paste(
      "Data are n (%) for categorical variables, mean ± SD for approximately normally",
      "distributed continuous variables, and median (IQR) for skewed continuous variables.",
      "Distribution assessed using skewness (|skewness| ≤ 0.5 indicates approximately normal).",
      "Statistical tests: *t*-test for normal continuous variables,",
      "Mann-Whitney *U* test for skewed continuous variables,",
      "χ² test for categorical variables (expected ≥5),",
      "Fisher's exact test for categorical variables (expected <5).",
      "Variables with <5% valid data were excluded from analysis.",
      "ICD = implantable cardioverter-defibrillator; IQR = interquartile range; SD = standard deviation.",
      sep = " "
    ))
  ) %>%
  tab_options(
    table.font.size = px(12),
    table.width = pct(100),
    heading.title.font.size = px(14),
    heading.title.font.weight = "bold",
    column_labels.font.weight = "bold",
    footnotes.font.size = px(10)
  )

print(gt_table)
gtsave(gt_table, "../tables/table-1-baseline.html")
gtsave(gt_table, "../tables/table-1-baseline.docx")

cat("  Saved: ../tables/table-1-baseline.html\n")
cat("  Saved: ../tables/table-1-baseline.docx\n")

# Excel
wb <- createWorkbook()
addWorksheet(wb, "Table 1")

writeData(wb, 1, "Table 1. Baseline Characteristics by Inappropriate ICD Therapy", 
          startRow = 1, colNames = FALSE)
addStyle(wb, 1, createStyle(fontSize = 14, textDecoration = "bold"), rows = 1, cols = 1)

writeData(wb, 1, final, startRow = 3)

headerStyle <- createStyle(
  fontSize = 11,
  fontName = "Times New Roman",
  textDecoration = "bold",
  halign = "center",
  border = "Bottom",
  borderColour = "black"
)
addStyle(wb, 1, headerStyle, rows = 3, cols = 1:ncol(final), gridExpand = TRUE)

var_rows <- 4:(nrow(final) + 3)
addStyle(wb, 1, createStyle(fontName = "Times New Roman", textDecoration = "bold"),
         rows = var_rows, cols = 1, gridExpand = TRUE, stack = TRUE)

centerStyle <- createStyle(fontName = "Times New Roman", halign = "center")
addStyle(wb, 1, centerStyle, rows = 4:(nrow(final) + 3), cols = 2:ncol(final), 
         gridExpand = TRUE, stack = TRUE)

sig_rows <- which(final$`P value` != "" & 
                    as.numeric(gsub("<", "", final$`P value`)) < 0.05) + 3
if (length(sig_rows) > 0) {
  addStyle(wb, 1, createStyle(fontName = "Times New Roman", textDecoration = "bold"),
           rows = sig_rows, cols = which(names(final) == "P value"),
           gridExpand = TRUE, stack = TRUE)
}

setColWidths(wb, 1, cols = 1, widths = 25)
setColWidths(wb, 1, cols = 2:6, widths = 15)

saveWorkbook(wb, "../tables/table-1-baseline.xlsx", overwrite = TRUE)
cat("  Saved: ../tables/table-1-baseline.xlsx\n")

# Summary
# ==============================================================================
cat("\n========================================\n")
cat("TABLE 1 COMPLETE\n")
cat("========================================\n\n")

cat("SAMPLE SIZE:\n")
cat("  Total: N =", final[Characteristic == "N", Total], "\n")
cat("  No inappropriate therapy:", 
    final[Characteristic == "N", `No Inappropriate Therapy`], "\n")
cat("  Inappropriate therapy:", 
    final[Characteristic == "N", `Inappropriate Therapy`], "\n\n")

cat("VARIABLES:\n")
cat("  Continuous (MEAN):", length(mean_vars), "variables\n")
cat("  Continuous (MEDIAN):", length(median_vars), "variables\n")
cat("  Binary/Categorical:", length(binary_vars), "variables\n")
cat("  Total rows in table:", nrow(final), "\n\n")

n_sig <- sum(!is.na(final$`P value`) & final$`P value` != "" &
               as.numeric(gsub("<", "", final$`P value`)) < 0.05)
cat("STATISTICAL TESTS:\n")
cat("  Significant differences (P < 0.05):", n_sig, "\n\n")

cat("OUTPUT FILES:\n")
cat("  1. ../tables/table-1-baseline.csv\n")
cat("  2. ../tables/table-1-baseline.html\n")
cat("  3. ../tables/table-1-baseline.docx\n")
cat("  4. ../tables/table-1-baseline.xlsx\n\n")

cat("========================================\n\n")

print(final)


print(table(ds$status_FIS, useNA = "always"))