
install.packages("survival")
install.packages("survminer")
install.packages("dplyr")
install.packages("svglite")
install.packages("rms")
install.packages("mice")


library(survival)
library(survminer)
library(dplyr)
library(rms)
library(mice)

combined <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/combined_dataset.csv")

data_dir <- "S:/AG/f-dhzc-profid/Data Transfer to Charite"
setwd("T:/Dokumente/PROFID/Study6")

vars_base <- c(
  "Survival_time", "Status", "BMI",
  "Age", "Sex", "Diabetes", "Hypertension", "Smoking", "MI_history",
  "LVEF", "NYHA",
  "eGFR", "Haemoglobin",
  "ACE_inhibitor_ARB", "Beta_blockers", "Lipid_lowering",
  "Revascularisation_acute"
)


missing_summary <- combined %>%
  select(all_of(vars_base)) %>%
  summarise(across(everything(),
                   ~ mean(is.na(.)) * 100))

t(missing_summary)




# Remove individuals with missing or implausible BMI and missing outcome data 
combined_clean <- combined %>%
  mutate(
    BMI    = suppressWarnings(as.numeric(BMI)),
    Status = suppressWarnings(as.numeric(Status))
  ) %>%
  filter(
    !is.na(BMI) &                
      !is.na(Status) &      
      !is.na(Survival_time)      
  )



write.csv(combined_clean, "combined_BMI_outcomefiltered.csv", row.names = FALSE)


impute_vars <- c(
  # Outcome
  "Survival_time", "Status",
  
  # Exposure
  "BMI",
  
  # Base confounders
  "Age", "Sex",
  "Diabetes", "Hypertension", "Smoking", "MI_history",
  "LVEF",
  "eGFR", "Haemoglobin",
  "ACE_inhibitor_ARB", "Beta_blockers", "Lipid_lowering",
  "Revascularisation_acute",
  
  #  auxiliary variables
  "Cholesterol", "LDL", "HDL", "Triglycerides",
  "Stroke_TIA", "HF",
  "Anti_platelet", "Anti_coagulant", "Diuretics", "Anti_anginal", "Calcium_antagonists",
  "Aldosterone_antagonist", "Digitalis_glycosides",
  "PCI_acute", "PCI", "CABG_acute", "CABG", "Thrombolysis_acute",
  "LBBB", "RBBB", "AF_atrial_flutter",
  "Baseline_type", "HasMRI", "CVD_risk_region", 
  "Anti_diabetic", "Anti_diabetic_oral", "Anti_diabetic_insulin",
  "MI_type"
)

imp_vars <- impute_vars 

imp_data <- combined_clean %>% select(all_of(impute_vars))

raw_subset <- combined_clean %>% select(all_of(imp_vars))

imo <- raw_subset %>%
  mutate(BMI = suppressWarnings(as.numeric(BMI))) %>%
  filter(!is.na(BMI))

pred <- make.predictorMatrix(imp_data)

# outcomes as predictors-only
pred[c("Survival_time","Status"), ] <- 0

# don't impute BMI, but allow it to predict others
pred["BMI", ] <- 0         # BMI will not be imputed
pred[, "BMI"] <- 1         # BMI can be used as a predictor


complete_cols <- names(imp_data)[colMeans(is.na(imp_data)) == 0]
if (length(complete_cols)){
  pred[complete_cols, ] <- 0   #  complete rows wont be imputed 
}

# let mice autoselect methods
meth <- make.method(imp_data)
meth["Survival_time"] <- ""
meth["Status"] <- ""


imp <- mice(imp_data, m = 20, maxit = 50, seed = 123, predictorMatrix = pred, 
            method = meth, printFlag = TRUE)

# check some key variables 
plot(imp, c("BMI", "eGFR", "Haemoglobin", "LDL"))

densityplot(imp, ~ BMI)
densityplot(imp, ~ eGFR)
stripplot(imp, BMI ~ .imp)

saveRDS(imp, "mice_imputed_data.RDS")

completed_imp1 <- complete(imp, action = 1)
write.csv(completed_imp1, "imputed_data_1_for_plots.csv", row.names = FALSE)

