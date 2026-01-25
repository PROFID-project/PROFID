###############################################
# Model fitting (ICD and Non-ICD cohorts)
# - Model 1: Null model
###############################################

# Packages
req <- c("data.table","dplyr","tidyverse","coxme","broom.mixed","lme4","gt")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Model 1 function
run_model_1 <- function(df, dataset_name) {
  df <- df %>%
    mutate(
      SCD_event = ifelse(Status == 1, 1, 0)
    ) %>%
    filter(!is.na(Status))
  # Safety checks
  if (n_distinct(df$ctr_id) < 2) return(NULL)
  if (sum(df$SCD_event, na.rm = TRUE) == 0) return(NULL)
  model <- glmer(
    SCD_event ~ 1 + (1 | ctr_id),
    data = df,
    family = binomial,
    control = glmerControl(optimizer = "bobyqa")
  )
  # Extract random effects variance
  re_var <- as.data.frame(VarCorr(model)) %>%
    filter(grp == "ctr_id") %>%
    pull(vcov)
  tibble(
    Dataset = dataset_name,
    N_patients = nrow(df),
    N_centres = n_distinct(df$ctr_id),
    N_events = sum(df$SCD_event),
    Centre_variance = round(re_var, 3),
    ICC = round(re_var / (re_var + (pi^2 / 3)), 3),
    AIC = round(AIC(model), 1),
    BIC = round(BIC(model), 1)
  )
}

# Run the function
model1_results <- bind_rows(
  run_model_1(fread("T:/Data Transfer to Charite/model_splits/ICD_raw_train.csv"),"Raw ICD cohort"),
  run_model_1(fread("T:/Data Transfer to Charite/model_splits/ICD_cmr_train.csv"),"CMR ICD cohort"),
  run_model_1(fread("T:/Data Transfer to Charite/model_splits/ICD_imp_train.csv"),"Imputed ICD cohort"),
  run_model_1(fread("T:/Data Transfer to Charite/model_splits/NR_raw_train.csv"),"Raw Non-ICD cohort (≤35%)"),
  run_model_1(fread("T:/Data Transfer to Charite/model_splits/NR_cmr_train.csv"),"CMR Non-ICD cohort (≤35%)"),
  run_model_1(fread("T:/Data Transfer to Charite/model_splits/NR_imp_train.csv"),"Imputed Non-ICD cohort (≤35%)"),
  run_model_1(fread("T:/Data Transfer to Charite/model_splits/NP_raw_train.csv"),"Raw Non-ICD cohort (>35%)"),
  run_model_1(fread("T:/Data Transfer to Charite/model_splits/NP_cmr_train.csv"),"CMR Non-ICD cohort (>35%)"),
  run_model_1(fread("T:/Data Transfer to Charite/model_splits/NP_imp_train.csv"),"Imputed Non-ICD cohort (>35%)")
)

# GT table
gt_model1 <- model1_results %>%
  gt() %>%
  tab_header(
    title = "Model 1: Null Multilevel Model",
    subtitle = "Random intercept for centre (baseline geographic variation in SCD risk)"
  ) %>%
  cols_label(
    Dataset = "Dataset",
    N_patients = "Patients (N)",
    N_centres = "Centres (N)",
    N_events = "SCD events (N)",
    Centre_variance = "Centre-level variance",
    ICC = "ICC",
    AIC = "AIC",
    BIC = "BIC"
  ) %>%
  fmt_number(
    columns = c(Centre_variance, ICC),
    decimals = 3
  ) %>%
  tab_spanner(
    label = "Model diagnostics",
    columns = c(AIC, BIC)
  ) %>%
  opt_all_caps() %>%
  opt_table_outline() %>%
  tab_source_note(
    source_note = "ICC calculated assuming logistic mixed model variance (π² / 3)."
  )
gtsave(gt_model1, "model1.html")
