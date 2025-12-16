library(mice)
library(survival)
library(rms)      # for rcs()

## Load the imputed mids object and saved fits
data          <- readRDS("mice_imputed_data_extended.RDS")
fit_list_cs1  <- readRDS("fit_list_cs1_extended.RDS")     # per-imputation fits
pool_fit_cs1  <- readRDS("pool_fit_cs1_extended.RDS")     # pooled object
curve_mice    <- read.csv("cox_RCS_BMI_curve_cs1_extended.csv")  # pooled BMI–HR curve
head(curve_mice)


make_subgroup_vars <- function(data) {
  data$age65 <- factor(ifelse(data$Age >= 65, ">=65", "<65"),
                       levels = c("<65", ">=65"))
  
  data$mi40 <- factor(ifelse(data$Survival_time >= 40, ">=40d", "<40d"),
                      levels = c("<40d", ">=40d"))
  
  data$LVEF_num <- as.numeric(as.character(data$LVEF))
  
  data$lvef_cat <- cut(data$LVEF,
                       breaks = c(-Inf, 30, 40, Inf),
                       right = FALSE,
                       labels = c("<30", "30-39", ">=40"))
  
  data$lvef_cat <- droplevels(data$lvef_cat)
  
  
  # Ensure binary subgroup vars are factors too 
  data$icd <- as.factor(data$ICD_status)
  data$diabetes <- as.factor(data$Diabetes)
  
  data
}

mi_lrt_interaction <- function(imp, time_var, event_var,
                               bmi_var,
                               covars,          # character vector of extended covariates (excluding subgroup var)
                               subgroup_var,    # e.g. "age65"
                               ties = "efron") {
  
  # BMI spline term: rcs(bmi, 4) = restricted cubic spline with 4 knots
  bmi_term <- sprintf("rcs(%s, 4)", bmi_var)
  
  # base includes subgroup main effect
  rhs_base <- paste(c(bmi_term, covars, subgroup_var), collapse = " + ")
  
  # interaction model: BMI spline x subgroup
  rhs_int  <- paste(c(sprintf("(%s) * %s", bmi_term, subgroup_var), covars), collapse = " + ")
  
  f_base <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, event_var, rhs_base))
  f_int  <- as.formula(sprintf("Surv(%s, %s) ~ %s", time_var, event_var, rhs_int))
  
  fit0 <- with(imp, {
    d <- make_subgroup_vars(data)
    coxph(f_base, data = d, ties = ties, x = TRUE)
  })
  
  fit1 <- with(imp, {
    d <- make_subgroup_vars(data)
    coxph(f_int, data = d, ties = ties, x = TRUE)
  })
  
  # Pooled likelihood ratio test for nested models under MI:
  test <- D2(fit1, fit0)
  
  # return a tidy row
  data.frame(
    subgroup = subgroup_var,
    statistic = unname(test$statistic),
    df = unname(test$df),
    p_value = unname(test$p.value),
    stringsAsFactors = FALSE
  )
}

extended_covars <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute", "Cholesterol", "HDL", "LDL", "Triglycerides", "Stroke_TIA", "ICD_status", "Cancer", "COPD_cat"
)

time_var  <- "Survival_time"     
event_var <- "SCD_Status"   
bmi_var   <- "BMI"        
s


results <- do.call(rbind, list(
  mi_lrt_interaction(imp, time_var, event_var, bmi_var, extended_covars, "icd"),
  mi_lrt_interaction(imp, time_var, event_var, bmi_var, extended_covars, "diabetes"),
  mi_lrt_interaction(imp, time_var, event_var, bmi_var, extended_covars, "age65"),
  mi_lrt_interaction(imp, time_var, event_var, bmi_var, extended_covars, "mi40"),
))
results
write.csv(results, "subgroup_interactions_mi_LRT_D2.csv", row.names = FALSE)

covars_no_lvef <- setdiff(extended_covars, "LVEF")

lvef_result <- mi_lrt_interaction(
  imp, time_var, event_var, bmi_var,
  covars_no_lvef,
  subgroup_var = "lvef_cat"
)

results2 <- rbind(results, lvef_result)
results2
write.csv(results2, "subgroup_interactions_mi_LRT_D2.csv", row.names = FALSE)

