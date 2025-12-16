library(mice)
library(survival)
library(rms)


# ---- Load imputed data ----
imp <- readRDS("mice_imputed_data_extended.RDS")
stopifnot(is.mids(imp))

# For rcs(): set datadist once from any completed dataset
dd <- datadist(complete(imp, 1))
options(datadist = "dd")

completed <- complete(imp, action = "all")  # list of completed datasets


# ---- model variables ----
time_var  <- "Survival_time"
event_var <- "Status_cs1"

extended_covars <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute","Cholesterol","HDL","LDL","Triglycerides",
  "Stroke_TIA","ICD_status","Cancer","COPD_cat"
)

# remove ICD_status from covars because subgroup captures it

covars_no_icd <- extended_covars


run_subgroup_LRT_linear <- function(completed, time_var, event_var, covars,
                                          subgroup_name, setup_fun) {
  
  # Keep spline main effect, but only linear interaction term
  f_base <- as.formula(paste(
    sprintf("Surv(%s, %s)", time_var, event_var),
    "~",
    paste(c("rcs(BMI, 4)", covars, subgroup_name), collapse = " + ")
  ))
  
  f_int <- as.formula(paste(
    sprintf("Surv(%s, %s)", time_var, event_var),
    "~",
    paste(c("rcs(BMI, 4)", covars, subgroup_name, sprintf("BMI:%s", subgroup_name)),
          collapse = " + ")
  ))
  
  fits0 <- lapply(completed, function(d) {
    d <- setup_fun(d)
    coxph(f_base, data = d, ties = "efron", x = TRUE)
  })
  
  fits1 <- lapply(completed, function(d) {
    d <- setup_fun(d)
    coxph(f_int, data = d, ties = "efron", x = TRUE)
  })
  
  fit0 <- list(analyses = fits0); class(fit0) <- "mira"
  fit1 <- list(analyses = fits1); class(fit1) <- "mira"
  
  test <- D2(fit1, fit0)
  
  data.frame(
    subgroup = subgroup_name,
    statistic = unname(test$statistic),
    df = unname(test$df),
    p_value = unname(test$p.value),
    stringsAsFactors = FALSE
  )
}



# ---- 1) ICD vs non-ICD ----
covars_no_icd <- setdiff(extended_covars, "ICD_status")

setup_icd <- function(d) {
  d$icd <- factor(d$ICD_status)
  d
}

res_icd <- run_subgroup_LRT_linear(
  completed, time_var, event_var,
  covars = covars_no_icd,
  subgroup_name = "icd",
  setup_fun = setup_icd
)


# ---- 2) Diabetes yes/no ----
covars_no_diabetes <- setdiff(extended_covars, "Diabetes")

setup_diabetes <- function(d) {
  d$diab <- factor(d$Diabetes)
  d
}

res_diab <- run_subgroup_LRT_linear(
  completed, time_var, event_var,
  covars = covars_no_diabetes,
  subgroup_name = "diab",
  setup_fun = setup_diabetes)


# ---- 3) Age <65 vs >=65 ----
covars_no_age <- setdiff(extended_covars, "Age")

setup_age <- function(d) {
  d$age65 <- factor(ifelse(d$Age >= 65, ">=65", "<65"),
                    levels = c("<65", ">=65"))
  d
}

res_age <- run_subgroup_LRT_linear(
  completed, time_var, event_var,
  covars = covars_no_age,
  subgroup_name = "age65",
  setup_fun = setup_age
)


# ---- 4) LVEF categories: <30, 30-39, >=40 ----
# Remove continuous LVEF from covars for this test (avoid duplication/collinearity)
covars_no_lvef <- setdiff(extended_covars, "LVEF")

setup_lvef <- function(d) {
  d$LVEF_num <- suppressWarnings(as.numeric(as.character(d$LVEF)))
  d$lvef_cat <- cut(d$LVEF_num,
                    breaks = c(-Inf, 30, 40, Inf),
                    right = FALSE,
                    labels = c("<30", "30-39", ">=40"))
  d$lvef_cat <- droplevels(d$lvef_cat)
  d
}

res_lvef <- run_subgroup_LRT_linear(
  completed, time_var, event_var,
  covars = covars_no_lvef,
  subgroup_name = "lvef_cat",
  setup_fun = setup_lvef
)


# ---- Combine + save ----
results <- rbind(res_icd, res_diab, res_age, res_lvef)
print(results)

write.csv(results, "subgroup_interactions_mi_LRT_D2.csv", row.names = FALSE)
