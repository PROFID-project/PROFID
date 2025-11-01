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


# --- 🔟 Identify categorical variables ---
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

