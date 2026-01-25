###############################################
# Model fitting (ICD and Non-ICD cohorts)
# - Model 4: Cross-level interactions
###############################################

# Packages
req <- c("data.table","dplyr","tidyverse","lme4","broom.mixed","gt")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Model 4 function
run_model_4 <- function(df, dataset_name) {
  df <- df %>%
    mutate(
      SCD_event = ifelse(Status == 1, 1, 0),
      sex = factor(Sex),
      CVD_risk_region = factor(CVD_risk_region)
    ) %>%
    filter(
      !is.na(ctr_id),
      !is.na(CVD_risk_region),
      !is.na(Status)
    )
  # Safety check
  if (sum(df$SCD_event) == 0) {
    message("Skipping ", dataset_name, ": no SCD events")
    return(NULL)
  }
  model <- glmer(
    SCD_event ~ 
      Age * CVD_risk_region +
      LVEF * CVD_risk_region +
      Sex + Diabetes + Beta_blockers +
      (1 | ctr_id),
    data = df,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa")
  )
  vc <- as.data.frame(VarCorr(model))
  centre_var <- vc$vcov[vc$grp == "ctr_id"]
  icc <- centre_var / (centre_var + (pi^2 / 3))
  tidy(model, effects = "fixed", conf.int = TRUE, exponentiate = TRUE) %>%
    filter(grepl(":", term)) %>%   # keep only interactions
    mutate(
      dataset = dataset_name,
      centre_variance = round(centre_var, 3),
      ICC = round(icc, 3),
      AIC = round(AIC(model), 1),
      BIC = round(BIC(model), 1)
    )
}

# Run the function
model4_results <- bind_rows(
  run_model_4(fread("T:/Data Transfer to Charite/model_splits/ICD_raw_train.csv"),"Raw ICD cohort"),
  run_model_4(fread("T:/Data Transfer to Charite/model_splits/ICD_imp_train.csv"),"Imputed ICD cohort"),
  run_model_4(fread("T:/Data Transfer to Charite/model_splits/NR_raw_train.csv"),"Raw Non-ICD cohort (≤35%)"),
  run_model_4(fread("T:/Data Transfer to Charite/model_splits/NR_imp_train.csv"),"Imputed Non-ICD cohort (≤35%)"),
  run_model_4(fread("T:/Data Transfer to Charite/model_splits/NP_raw_train.csv"),"Raw Non-ICD cohort (>35%)"),
  run_model_4(fread("T:/Data Transfer to Charite/model_splits/NP_imp_train.csv"),"Imputed Non-ICD cohort (>35%)")
)

# GT table
gt_model4 <- model4_results %>%
  select(
    dataset,
    term,
    estimate,
    conf.low,
    conf.high,
    p.value,
    centre_variance,
    ICC,
    AIC,
    BIC
  ) %>%
  mutate(
    OR_CI = sprintf("%.2f (%.2f–%.2f)", estimate, conf.low, conf.high),
    p.value = sprintf("%.3f", p.value)
  ) %>%
  select(-estimate, -conf.low, -conf.high) %>%
  gt(groupname_col = "dataset") %>%
  tab_header(
    title = "Model 4: Cross-level Interaction Multilevel Model",
    subtitle = "Interactions between patient-level risk factors and geographic regions"
  ) %>%
  cols_label(
    dataset = "Dataset",
    term = "Interaction Term",
    OR_CI = "OR (95% CI)",
    p.value = "P-value",
    centre_variance = "Centre-level variance",
    ICC = "ICC",
    AIC = "AIC",
    BIC = "BIC"
  ) %>%
  tab_source_note(
    source_note = "Logistic mixed-effects models with random intercepts for centre. ICC calculated assuming logistic variance π²/3."
  )
gtsave(gt_model4, "model4.html")
