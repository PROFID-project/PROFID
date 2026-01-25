###############################################
# Model fitting (ICD and Non-ICD cohorts)
# - Model 2: Patient-level model
###############################################

# Packages
req <- c("data.table","dplyr","tidyverse","lme4","gt")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Model 2 function
run_model_2 <- function(df, dataset_name) {
  df <- df %>%
    mutate(
      SCD_event = ifelse(Status == 1, 1, 0)
    ) %>%
    filter(!is.na(Status)) %>%
    drop_na(
      SCD_event, Age, Sex, LVEF,
      Diabetes, Hypertension, Smoking,
      ACE_inhibitor_ARB, ACE_inhibitor, ARB, Aldosterone_antagonist, 
      Anti_anginal, Anti_arrhythmic_III, Anti_coagulant, Anti_diabetic, 
      Anti_diabetic_insulin, Anti_diabetic_oral, Anti_platelet, Beta_blockers, 
      Calcium_antagonists, Digitalis_glycosides, Diuretics, Lipid_lowering,
      ctr_id
    )
  # Safety check
  if (sum(df$SCD_event) == 0) {
    message("Skipping ", dataset_name, ": no SCD events")
    return(NULL)
  }
  model <- glmer(
    SCD_event ~
      Age + Sex + LVEF +
      Diabetes + Hypertension + Smoking +
      ACE_inhibitor_ARB + ACE_inhibitor + ARB + Aldosterone_antagonist + 
      Anti_anginal + Anti_arrhythmic_III + Anti_coagulant + Anti_diabetic + 
      Anti_diabetic_insulin + Anti_diabetic_oral + Anti_platelet + Beta_blockers + 
      Calcium_antagonists + Digitalis_glycosides + Diuretics + Lipid_lowering +
      (1 | ctr_id),
    data = df,
    family = binomial,
    control = glmerControl(
      optimizer = "bobyqa",
      optCtrl = list(maxfun = 2e5)
    )
  )
  # Extract centre-level variance
  var_centre <- as.data.frame(VarCorr(model))$vcov[1]
  # ICC for logistic mixed model
  ICC <- var_centre / (var_centre + (pi^2 / 3))
  tibble(
    dataset = dataset_name,
    patients_n = nrow(df),
    centres_n = n_distinct(df$ctr_id),
    SCD_events = sum(df$SCD_event),
    centre_variance = round(var_centre, 3),
    ICC = round(ICC, 3),
    AIC = round(AIC(model), 1),
    BIC = round(BIC(model), 1)
  )
}

# Run the function
model2_results <- bind_rows(
  run_model_2(fread("T:/Data Transfer to Charite/model_splits/ICD_raw_train.csv"),"Raw ICD cohort"),
  run_model_2(fread("T:/Data Transfer to Charite/model_splits/ICD_cmr_train.csv"),"CMR ICD cohort"),
  run_model_2(fread("T:/Data Transfer to Charite/model_splits/ICD_imp_train.csv"),"Imputed ICD cohort"),
  run_model_2(fread("T:/Data Transfer to Charite/model_splits/NR_raw_train.csv"),"Raw Non-ICD cohort (≤35%)"),
  run_model_2(fread("T:/Data Transfer to Charite/model_splits/NR_cmr_train.csv"),"CMR Non-ICD cohort (≤35%)"),
  run_model_2(fread("T:/Data Transfer to Charite/model_splits/NR_imp_train.csv"),"Imputed Non-ICD cohort (≤35%)"),
  run_model_2(fread("T:/Data Transfer to Charite/model_splits/NP_raw_train.csv"),"Raw Non-ICD cohort (>35%)"),
  run_model_2(fread("T:/Data Transfer to Charite/model_splits/NP_cmr_train.csv"),"CMR Non-ICD cohort (>35%)"),
  run_model_2(fread("T:/Data Transfer to Charite/model_splits/NP_imp_train.csv"),"Imputed Non-ICD cohort (>35%)")
)

# GT table
gt_model2 <- model2_results %>%
  gt() %>%
  tab_header(
    title = "Model 2: Patient-Level Multilevel Model",
    subtitle = "Adjusted for age, sex, LVEF, comorbidities, and medications"
  ) %>%
  cols_label(
    dataset = "Dataset",
    patients_n = "Patients (N)",
    centres_n = "Centres (N)",
    SCD_events = "SCD events (N)",
    centre_variance = "Centre-level variance",
    ICC = "ICC",
    AIC = "AIC",
    BIC = "BIC"
  ) %>%
  fmt_number(
    columns = c(centre_variance, ICC),
    decimals = 3
  ) %>%
  tab_source_note(
    source_note = "ICC calculated assuming logistic mixed model variance (π²/3)."
  )
gtsave(gt_model2, "model2.html")