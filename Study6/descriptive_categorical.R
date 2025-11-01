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
setwd("//charite.de/homes/h02/clco10/Dokumente/PROFID/Study6")

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
categorical[move_to_categorical] <- continuous[move_to_categorical]
continuous <- continuous[, !(names(continuous) %in% move_to_categorical), drop = FALSE]


categorical_vars <- setdiff(names(categorical), c("X","BMI_cat","X","DB","Status","ID"))

cat_summary <- combined %>%
  select(all_of(c("BMI_cat", categorical_vars))) %>%
  mutate(across(all_of(categorical_vars), as.character)) %>%
  pivot_longer(
    -BMI_cat,
    names_to = "Variable",
    values_to = "Level"
  ) %>%
  group_by(BMI_cat, Variable, Level) %>%
  summarise(Count = n(), .groups = "drop_last") %>%
  mutate(Percent = 100 * Count / sum(Count)) %>%
  ungroup()

chi_sq_results <- map_dfr(
  categorical_vars,
  function(var) {
    tbl <- table(combined[[var]], combined$BMI_cat)
    if (min(dim(tbl)) < 2) {
      return(tibble(
        Variable = var,
        statistic = NA_real_,
        p.value = NA_real_,
        method = "Chi-squared (not enough categories)"
      ))
    }
    test <- tryCatch(chisq.test(tbl), error = function(e) NULL)
    if (is.null(test)) {
      tibble(
        Variable = var,
        statistic = NA_real_,
        p.value = NA_real_,
        method = "Chi-squared (error)"
      )
    } else {
      tidy(test) %>%
        mutate(Variable = var)
    }
  }
) %>%
  mutate(p_adj = p.adjust(p.value, method = "fdr"))

cat_final <- cat_summary %>%
  left_join(chi_sq_results %>% select(Variable, p.value, p_adj), by = "Variable")

cat_final_pretty <- cat_final %>%
  mutate(
    Percent = sprintf("%.1f%%", Percent),
    p.value_rounded = signif(p.value, 3),
    p_adj_rounded = signif(p_adj, 3)
  ) %>%
  select(Variable, Level, BMI_cat, Count, Percent, p.value_rounded, p_adj_rounded)

write.csv(cat_final_pretty,
          file = "final_categorical_summary.csv",
          row.names = FALSE)

