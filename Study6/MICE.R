
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

categorical_vars <- c(
  # Base confounders
  "Sex",
  "Diabetes", "Hypertension", "Smoking", "MI_history",
  "ACE_inhibitor_ARB", "Beta_blockers", "Lipid_lowering",
  "Revascularisation_acute",
  
  #  auxiliary variables
  "Stroke_TIA", "HF",
  "Anti_platelet", "Anti_coagulant", "Diuretics", "Anti_anginal", "Calcium_antagonists",
  "Aldosterone_antagonist", "Digitalis_glycosides",
  "PCI_acute", "PCI", "CABG_acute", "CABG", "Thrombolysis_acute",
  "LBBB", "RBBB", "AF_atrial_flutter",
  "Baseline_type", "HasMRI", "CVD_risk_region", 
  "Anti_diabetic", "Anti_diabetic_oral", "Anti_diabetic_insulin",
  "MI_type"
)

for (v in categorical_vars) {
  if (v %in% names(combined_clean)) {
    combined_clean[[v]] <- as.factor(combined_clean[[v]])
  }
}

sapply(combined_clean[categorical_vars], class)


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

# --- Lipid stabilisation ---
lipids <- intersect(c("Cholesterol", "LDL", "HDL", "Triglycerides"), names(imp_data))


# let mice autoselect methods
meth <- make.method(imp_data)
meth["Survival_time"] <- ""
meth["Status"] <- ""

# Use PMM for the lipids
meth[lipids] <- "pmm"

# Restrict lipid predictors (no cross-prediction among lipids)
pred[lipids, ] <- 0
lipid_base_preds <- intersect(
  c("Age","Sex","BMI","Diabetes","Hypertension","Smoking",
    "Lipid_lowering","ACE_inhibitor","ARB","Beta_blockers",
    "eGFR","Haemoglobin","MI_history"),
  names(imp_data)
)
pred[lipids, lipid_base_preds] <- 1

set.seed(123)
imp <- mice(
  imp_data,
  m = 20,
  maxit = 20,       
  method = meth,
  predictorMatrix = pred,
  pmm.k = 10,         # more donor candidates
  printFlag = TRUE
)

# Diagnostics
plot(imp, c(intersect(c("HDL","LDL","Cholesterol"), names(imp$imp)),
            intersect("logTrig", names(imp$imp))))
densityplot(imp, ~ HDL + LDL + Cholesterol)
if ("logTrig" %in% names(imp$data)) densityplot(imp, ~ logTrig)


# check some key variables 
plot(imp, c("Age", "eGFR", "Haemoglobin", "LDL", "HDL", "Triglycerides", "LVEF", "Cholesterol"))

stripplot(imp, Diabetes ~ .imp, pch = 20, cex = 1.2)
stripplot(imp, Sex ~ .imp, pch = 20, cex = 1.2)

densityplot(imp, ~ as.numeric(NYHA))

densityplot(imp, ~ BMI)
densityplot(imp, ~ eGFR)
stripplot(imp, BMI ~ .imp)

saveRDS(imp, "mice_imputed_data.RDS")

completed_imp1 <- complete(imp, action = 1)
write.csv(completed_imp1, "imputed_data_1_for_plots.csv", row.names = FALSE)

