###############################################
# Model validation (ICD and Non-ICD cohorts)
###############################################

# Packages
req <- c("data.table","dplyr","tidyverse","lme4","purrr","ggplot2","gt")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Data preparation
data_prep <- function(df_train) {
  df_train <- df_train %>%
    mutate(
      SCD_event = ifelse(Status == 1, 1, 0)
    ) %>%
    filter(
      !is.na(ctr_id),
      !is.na(SCD_event)
    )
}

# Model fitting
fit_model4 <- function(df_train) {
  glmer(
    SCD_event ~ Age + Sex + LVEF + Diabetes + Hypertension +
      (1 | ctr_id),
    data = df_train,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa")
  )
}

# Prediction function
predict_risk <- function(model, df_valid) {
  df_valid <- df_valid %>%
    mutate(
      pred_risk = predict(
        model,
        newdata = df_valid,
        type = "response",
        allow.new.levels = TRUE
      )
    ) %>%
    filter(!is.na(pred_risk))
}

# Calibration dataframe
make_calibration_df <- function(df, bins = 10) {
  df <- df %>%
    mutate(risk_bin = ntile(pred_risk, bins)) %>%
    group_by(risk_bin) %>%
    summarise(
      mean_pred = mean(pred_risk),
      obs_rate  = mean(SCD_event),
      n = n(),
      .groups = "drop"
    )
}

# Validation metrics
calc_validation_metrics <- function(df) {
  # Calibration intercept & slope
  cal_model <- glm(
    SCD_event ~ qlogis(pred_risk),
    data = df,
    family = binomial
  )
  tibble(
    intercept = coef(cal_model)[1],
    slope     = coef(cal_model)[2],
    brier     = mean((df$SCD_event - df$pred_risk)^2)
  )
}

# Run validation
run_validation <- function(df_train, df_valid, cohort, data_type) {
  df_train <- data_prep(df_train)
  df_valid <- data_prep(df_valid)
  model <- fit_model4(df_train)
  df_p  <- predict_risk(model, df_valid)
  list(
    calibration = make_calibration_df(df_p) %>%
      mutate(cohort = cohort, data_type = data_type),
    metrics = calc_validation_metrics(df_p) %>%
      mutate(cohort = cohort, data_type = data_type)
  )
}

# Results
results <- list(
  run_validation(fread("T:/Data Transfer to Charite/model_splits/ICD_raw_train.csv"), fread("T:/Data Transfer to Charite/model_splits/ICD_raw_valid.csv"), "ICD", "Raw"),
  run_validation(fread("T:/Data Transfer to Charite/model_splits/ICD_imp_train.csv"), fread("T:/Data Transfer to Charite/model_splits/ICD_imp_valid.csv"), "ICD", "Imputed"),
  run_validation(fread("T:/Data Transfer to Charite/model_splits/NR_raw_train.csv"), fread("T:/Data Transfer to Charite/model_splits/NR_raw_valid.csv"), "Non-ICD ≤35%", "Raw"),
  run_validation(fread("T:/Data Transfer to Charite/model_splits/NR_imp_train.csv"), fread("T:/Data Transfer to Charite/model_splits/NR_imp_valid.csv"), "Non-ICD ≤35%", "Imputed"),
  run_validation(fread("T:/Data Transfer to Charite/model_splits/NP_raw_train.csv"), fread("T:/Data Transfer to Charite/model_splits/NP_raw_valid.csv"), "Non-ICD >35%", "Raw"),
  run_validation(fread("T:/Data Transfer to Charite/model_splits/NP_imp_train.csv"), fread("T:/Data Transfer to Charite/model_splits/NP_imp_valid.csv"), "Non-ICD >35%", "Imputed")
)

# Calibration and metrics
cal_all <- map(results, "calibration") %>% bind_rows()
metrics_all <- map(results, "metrics") %>% bind_rows()

# Calibration plots
p <- ggplot(cal_all, aes(mean_pred, obs_rate)) +
  geom_point(size = 2.5, colour = "#0072B2") +
  geom_line(colour = "#0072B2") +
  geom_abline(
    intercept = 0, slope = 1,
    linetype = "dashed", colour = "grey40"
  ) +
  facet_grid(cohort ~ data_type) +
  labs(
    title = "Internal Validation: Calibration of Predicted SCD Risk",
    x = "Mean predicted risk",
    y = "Observed event rate"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "bold"),
    strip.text = element_text(face = "bold")
  )
ggsave(p, filename = "T:/Data Transfer to Charite/calibration_plot.png")

# GT table
validation_gt <- metrics_all %>%
  gt() %>%
  tab_header(
    title = "Internal Validation Metrics",
    subtitle = "Calibration performance using 90:10 centre-stratified split"
  ) %>%
  fmt_number(
    columns = c(intercept, slope, brier),
    decimals = 3
  ) %>%
  cols_label(
    intercept = "Calibration intercept",
    slope = "Calibration slope",
    brier = "Brier score"
  ) %>%
  tab_spanner(
    label = "Calibration",
    columns = c(intercept, slope)
  )
gtsave(validation_gt, "T:/Data Transfer to Charite/validation_metrics.html")

