install.packages("dplyr")  
install.packages("tidyr")  
install.packages("broom")
library(dplyr)
library(tidyr)
library(purrr)
library(broom)

icd <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/ICD.csv")
nonicd_preserved <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_preserved.csv")
nonicd_reduced <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_reduced.csv")
combined <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/combined_dataset.csv")

data_dir <- "S:/AG/f-dhzc-profid/Data Transfer to Charite"
setwd("//charite.de/homes/h02/clco10/Dokumente/PROFID/Study6s")

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

move_to_categorical <- c("X", "Status", "Time_zero_Y", "CVD_risk_region", "IsSWHR")
categorical[move_to_categorical] <- combined[move_to_categorical]
continuous <- continuous[, !(names(continuous) %in% move_to_categorical), drop = FALSE]


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
          file = "continuous_descriptors.csv",
          row.names = FALSE)

final_summary



