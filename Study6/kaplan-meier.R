library(survival)
library(survminer)
library(dplyr)

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

logrank_test <- survdiff(Surv(Survival_time, Status_KM) ~ BMI_cat, data = combined_KM)
logrank_test

p_value <- 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
p_value

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

