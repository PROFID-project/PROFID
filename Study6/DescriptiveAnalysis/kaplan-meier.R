library(survival)
library(survminer)
library(dplyr)

combined <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/combined_dataset.csv")
setwd("T:/Dokumente/PROFID/Study6")

run_km_at_horizon <- function(df, horizon_months, out_prefix = "KM") {
  
  df_h <- df %>%
    mutate(
      BMI_cat = as.factor(BMI_cat),
      time_h  = pmin(Survival_time, horizon_months),
      event_h = ifelse(Survival_time <= horizon_months & Status == 1, 1, 0)
    )
  
  # sanity checks
  stopifnot(max(df_h$time_h, na.rm = TRUE) <= horizon_months)
  stopifnot(sum(df_h$event_h == 1 & df_h$Survival_time > horizon_months, na.rm = TRUE) == 0)
  
  fit <- survfit(Surv(time_h, event_h) ~ BMI_cat, data = df_h)
  
  km_plot <- ggsurvplot(
    fit,
    data = df_h,
    risk.table = TRUE,
    pval = TRUE,
    conf.int = FALSE,
    censor = FALSE,
    ylim = c(0.8, 1.00),
    break.time.by = 20,
    xlim = c(0, horizon_months),
    xlab = "Time (months)",
    ylab = "SCD-free survival probability",
    legend.title = "BMI Category",
    palette = "Dark2",
    ggtheme = theme_minimal()
  )
  
  # Save plot
  ggsave(
    filename = sprintf("%s_Survival_by_BMI_%dmo.pdf", out_prefix, horizon_months),
    plot = km_plot$plot,
    width = 7,
    height = 5,
    device = "svg"
  )
  
  # PNG (high-resolution for slides / review)
  ggsave(
    filename = sprintf("%s_Survival_by_BMI_%dmo.png", out_prefix, horizon_months),
    plot = km_plot$plot,
    width = 7,
    height = 5,
    dpi = 300
  )
  
  km_summary_df <- survminer::surv_summary(fit, data = df_h)
  write.csv(km_summary_df,
            sprintf("%s_Survival_Summary_%dmo.csv", out_prefix, horizon_months),
            row.names = FALSE)
  
  
  
  # Log-rank test
  logrank_test <- survdiff(Surv(time_h, event_h) ~ BMI_cat, data = df_h)
  
  logrank_results <- data.frame(
    horizon_months = horizon_months,
    chi_squared = unname(logrank_test$chisq),
    degrees_of_freedom = length(logrank_test$n) - 1,
    p_value = 1 - pchisq(logrank_test$chisq, length(logrank_test$n) - 1)
  )
  
  write.csv(logrank_results,
            sprintf("%s_LogRank_Test_%dmo.csv", out_prefix, horizon_months),
            row.names = FALSE)
  
  return(list(
    horizon = horizon_months,
    n = nrow(df_h),
    events_within_horizon = sum(df_h$event_h == 1, na.rm = TRUE),
    logrank = logrank_results
  ))
}

res_60  <- run_km_at_horizon(combined, 60)
res_90  <- run_km_at_horizon(combined, 90)
res_120 <- run_km_at_horizon(combined, 120)

# Optional: combine log-rank results into one table
all_logrank <- bind_rows(res_60$logrank, res_90$logrank, res_120$logrank)
write.csv(all_logrank, "KM_LogRank_AllHorizons.csv", row.names = FALSE)

all_logrank

