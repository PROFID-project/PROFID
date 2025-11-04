install.packages("survival")
install.packages("survminer")
install.packages("dplyr")
install.packages("svglite")

library(survival)
library(survminer)
library(dplyr)

combined <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/combined_dataset.csv")

data_dir <- "S:/AG/f-dhzc-profid/Data Transfer to Charite"
setwd("T:/Dokumente/PROFID/Study6")


combined_KM <- combined %>%
  mutate(
    Status_KM = ifelse(Status == 1, 1, 0),
    BMI_cat = as.factor(BMI_cat)
  )

fit_KM <- survfit(Surv(Survival_time, Status_KM) ~ BMI_cat, data = combined_KM)


km_plot <- ggsurvplot(
  fit_KM,
  data = combined_KM,
  risk.table = TRUE,
  pval = TRUE,
  conf.int = TRUE,
  ylim = c(0.5,1.00),
  break.time.by = 20,
  xlab = "Time (days)",
  ylab = "SCD-free survival probability",
  legend.title = "BMI Category",
  surv.median.line = "hv",
  palette = "Dark2",
  ggtheme = theme_minimal()
)

ggsave(
  filename = "KM_Survival_by_BMI.svg",
  plot = km_plot$plot,
  width = 7,
  height = 5,
  device = "svg"
)

km_plot

summary_fit <- summary(fit_KM)
write.csv(as.data.frame(summary_fit), "KM_Survival_Summary.csv", row.names = FALSE)


logrank_test <- survdiff(Surv(Survival_time, Status_KM) ~ BMI_cat, data = combined_KM)
logrank_test


# Extract results
logrank_results <- data.frame(
  chi_squared = logrank_test$chisq,
  degrees_of_freedom = length(logrank_test$n) - 1,
  p_value = 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
)
p_value

# Write to CSV
write.csv(logrank_results, "KM_LogRank_Test.csv", row.names = FALSE)


