install.packages("dplyr")  
install.packages("tidyr")  
install.packages("broom")
install.packages("lubridate")
install.packages("mice")
library(dplyr)
library(tidyr)
library(purrr)
library(broom)
library(lubridate)
library

data_dir <- "S:/AG/f-dhzc-profid/Data Transfer to Charite"
setwd("T:/Dokumente/PROFID/Study6")

icd <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/ICD.csv")
nonicd_preserved <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_preserved.csv")
nonicd_reduced <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_reduced.csv")


# Load the imputed object
imp <- readRDS("mice_imputed_data.RDS")

# Extract the first imputed dataset for descriptive analyses
combined <-mice::complete(imp, action = 1)


fulldata <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/combined_dataset.csv")



# Specify the columns you want to summarise
cols_to_check <- c("Cancer", "Stroke_TIA", "Diabetes", "COPD")

# Create counts for selected columns (including NAs)
freq_list <- lapply(fulldata[cols_to_check], function(x) as.data.frame(table(x, useNA = "ifany")))

# Add column identifiers
freq_list <- lapply(names(freq_list), function(nm) {
  tmp <- freq_list[[nm]]
  names(tmp) <- c("value", "count")
  tmp$column <- nm
  tmp
})

# Combine into one dataframe
freq_df <- do.call(rbind, freq_list)

# Reorder columns
freq_df <- freq_df[, c("column", "value", "count")]

# Write to CSV
write.csv(freq_df, "comorbidity_value_counts.csv", row.names = FALSE)


summary(combined$BMI)
summary(combined$Status)
names(combined)

combined <- combined %>%
  mutate(BMI_cat = cut(
    BMI,
    breaks = c(-Inf, 18.5, 25, 30, 35, 40, Inf),
    labels = c("Underweight", "Normal", "Overweight", "Obese I", "Obese II", "Obese III")
  ))

numeric_cols <- sapply(combined, is.numeric)

continuous <- combined[, numeric_cols, drop = FALSE]
names(continuous)

categorical <- combined[, !numeric_cols, drop = FALSE]

move_to_categorical <- c( "Status",  "CVD_risk_region")
categorical[move_to_categorical] <- combined[move_to_categorical]
continuous <- continuous[, !(names(continuous) %in% move_to_categorical), drop = FALSE]


# Assign cohort-level tags
#icd_tag <- icd %>%
#  transmute(ID = ID,
#            ICD_status = 1L,
#            source = "ICD")

#nonicd_pres_tag <- nonicd_preserved %>%
#  transmute(ID = ID,
#            ICD_status = 0L,
#            source = "NonICD_preserved")

#nonicd_red_tag <- nonicd_reduced %>%
#  transmute(ID = ID,
#            ICD_status = 0L,
#            source = "NonICD_reduced")

# Combine and remove duplicates: keep ICD if conflicts
#cohort_tag <- bind_rows(icd_tag, nonicd_pres_tag, nonicd_red_tag) %>%
#  arrange(desc(ICD_status)) %>%      # ICD=1 rows come first
#  distinct(ID, .keep_all = TRUE)

#combined <- combined %>%
#  left_join(cohort_tag, by = "ID")

# Check for duplicates
#sum(duplicated(combined$ID))

# Save cleaned, deduplicated combined dataset
#write.csv(
#  combined,
#  "S:/AG/f-dhzc-profid/Data Transfer to Charite/combined_icd.csv",
#  row.names = FALSE
#)


summary_continuous_tidy <- map_dfr(
  names(continuous),
  ~ tibble(
    Variable = .x,
    Median = median(continuous[[.x]], na.rm = TRUE),
    IQR = IQR(continuous[[.x]], na.rm = TRUE)
  )
)

summary_continuous_tidy


write.csv(summary_continuous_tidy,
          file = "summary_continuous.csv",
          row.names = FALSE)


summary_continuous_by_BMI <- map_dfr(
  names(continuous),
  function(var) {
    combined %>%
      group_by(BMI_cat) %>%
      summarise(
        Median = median(.data[[var]], na.rm = TRUE),
        IQR = IQR(.data[[var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Variable = var)
  }
)

summary_continuous_by_BMI <- summary_continuous_by_BMI %>%
  select(Variable, BMI_cat, Median, IQR)

summary_continuous_by_BMI


outcome_var <- "Status"

safe_kw_subset <- function(df, var, group_col) {
  # Skip if var missing
  if (!(var %in% names(df))) {
    return(tibble(
      Variable = var,
      statistic = NA_real_,
      p.value = NA_real_,
      parameter = NA_real_,
      method = "variable missing"
    ))
  }
  
  x <- df[[var]]
  g <- df[[group_col]]
  
  # convert group to observed character values so factor levels don’t mislead
  g_char <- as.character(g)
  g_unique <- unique(na.omit(g_char))
  
  # check: at least 2 distinct observed Status values
  if (length(g_unique) < 2) {
    return(tibble(
      Variable = var,
      statistic = NA_real_,
      p.value = NA_real_,
      parameter = NA_real_,
      method = "Kruskal-Wallis (not enough groups)"
    ))
  }
  
  # run test safely
  result <- tryCatch(
    {
      test <- kruskal.test(x ~ g_char)
      broom::tidy(test)
    },
    error = function(e) {
      tibble(
        statistic = NA_real_,
        p.value = NA_real_,
        parameter = NA_real_,
        method = paste("Error:", conditionMessage(e))
      )
    }
  )
  
  result %>% mutate(Variable = var)
}


kruskal_by_BMI <- combined %>%
  group_by(BMI_cat) %>%
  group_modify(~ {
    map_dfr(
      names(continuous),
      function(v) safe_kw_subset(.x, v, outcome_var)
    )
  }) %>%
  ungroup() %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

kruskal_by_BMI


summary_continuous_by_BMI <- map_dfr(
  names(continuous),
  function(var) {
    combined %>%
      group_by(BMI_cat) %>%
      summarise(
        Median = median(.data[[var]], na.rm = TRUE),
        IQR = IQR(.data[[var]], na.rm = TRUE),
        .groups = "drop"
      ) %>%
      mutate(Variable = var)
  }
)

summary_continuous_by_BMI <- summary_continuous_by_BMI %>%
  select(Variable, BMI_cat, Median, IQR)

final_summary <- summary_continuous_by_BMI %>%
  # make absolutely sure types match:
  mutate(
    Variable = as.character(Variable),
    BMI_cat  = as.character(BMI_cat)
  ) %>%
  left_join(
    kruskal_by_BMI %>%
      mutate(
        Variable = as.character(Variable),
        BMI_cat  = as.character(BMI_cat)
      ) %>%
      select(Variable, BMI_cat, p.value, p_adj),
    by = c("Variable", "BMI_cat")
  )

write.csv(final_summary,
          file = "continuous_descriptors_imputed.csv",
          row.names = FALSE)

final_summary

length(intersect(icd$ID, nonicd_preserved$ID))
length(intersect(icd$ID, nonicd_reduced$ID))
length(intersect(nonicd_preserved$ID, nonicd_reduced$ID))

library(dplyr)
install.packages("gtsummary")
library(gtsummary)

table1_df <- combined %>%
  mutate(
    BMI_cat = factor(BMI_cat,
                     levels = c("Underweight","Normal","Overweight","Obese I","Obese II","Obese III")
    ),
    # binary outcome for Table 1
    SCD_status = factor(ifelse(Status == 1, 1, 0),
                        levels = c(0, 1),
                        labels = c("No SCD", "SCD")),
    # (optional) 3-level outcome if you also want it displayed somewhere
    Status3 = factor(Status,
                     levels = c(0, 1, 2),
                     labels = c("No event", "SCD", "Other death"))
  ) %>%
  select(
    BMI_cat, SCD_status, # <- include this for the SCD outcome row
    Age, Sex, Diabetes, Hypertension, Smoking, MI_history,
    MI_type, Baseline_type, CVD_risk_region,
    LVEF, eGFR, Haemoglobin,
    Cholesterol, HDL, LDL, Triglycerides,
    ACE_inhibitor_ARB, Beta_blockers, Lipid_lowering,
    Revascularisation_acute
  )

table1 <- table1_df %>%
  tbl_summary(
    by = BMI_cat,
    statistic = list(
      all_continuous() ~ "{median} ({p25}, {p75})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "no"
  ) %>%
  add_overall(last = TRUE, col_label = "Total") %>%  # adds Total column
  add_n(location = "label") %>%                      # adds N in header
  add_p(test = list(
    all_continuous() ~ "kruskal.test",
    all_categorical() ~ "chisq.test"
  )) %>%
  bold_labels()

table1

library(gtsummary)

table1_csv <- as_tibble(table1, col_labels = TRUE)

write.csv(
  table1_csv,
  file = "Table1_baseline_characteristics.csv",
  row.names = FALSE
)

