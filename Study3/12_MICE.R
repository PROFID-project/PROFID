##########################SAP sensitivity: MICE (m=20) for covariates <30% missing#####################################################

# install.packages("data.table")
#install.packages("cmprsk")
# install.packages("survival")
# install.packages("readxl")
# install.packages ("mice")

library(data.table)
library(cmprsk)
library(survival)
library(readxl)
library (mice)


dt <- as.data.table(readRDS("study3_analysis_final.rds"))


# Exposure
dt[, device_group := factor(device_group, levels = c("Dual","Single"))]

# Dataset (strata)
dt[, dataset := factor(dataset)]

#EVENT
dt[, event_inapp_shock := as.integer(event_inapp_shock)]

#sex
dt[, Sex := factor(Sex)]
if (all(c("Female","Male") %in% levels(dt$Sex))) dt[, Sex := relevel(Sex, ref = "Female")]

# Diabetes 
dt[, Diabetes := factor(Diabetes, levels = c("No","Yes"))]

# FOLLOW UP
dt[, t_followup_days_final := as.numeric(t_followup_days_final)]

############# SAP eligibility: keep only vars with <30% missingness for imputation #

# Candidate predictors 
cand <- c("device_group","age_icd","Sex","LVEF","Diabetes","dataset",
          "t_followup_days_final","event_inapp_shock")

# Compute missingness 
miss <- sapply(dt[, ..cand], function(x) mean(is.na(x)))
print(round(miss, 3))

# SAP threshold
eligible_impute <- names(miss)[miss > 0 & miss < 0.30]
eligible_impute
#  Diabetes NOT IMP
eligible_impute <- setdiff(eligible_impute, "Diabetes")

# Variables with >=30% missing (excluded from MICE per SAP)
excluded_high_miss <- names(miss)[miss >= 0.30]
excluded_high_miss

# ######## Build MICE data frame ##########
imp_vars <- cand  
imp_df <- as.data.frame(dt[, ..imp_vars])

############# MICE setup ###################
ini <- mice(imp_df, maxit = 0, printFlag = FALSE)
meth <- ini$method
pred <- ini$predictorMatrix

no_impute <- c("t_followup_days_final","event_inapp_shock","device_group","dataset")
meth[no_impute] <- ""

# Methods: pmm for numeric; factors left unimputed here
if ("LVEF" %in% names(meth))     meth["LVEF"]     <- "pmm"
if ("age_icd" %in% names(meth))  meth["age_icd"]  <- ""  
if ("Diabetes" %in% names(meth)) meth["Diabetes"] <- ""
if ("Sex" %in% names(meth))      meth["Sex"]      <- ""


for (v in setdiff(names(meth), eligible_impute)) {
  if (!(v %in% no_impute)) meth[v] <- ""
}


diag(pred) <- 0


##########Run MICE (m = 20)############
set.seed(20260113)
imp <- mice(
  imp_df,
  m = 20,
  method = meth,
  predictorMatrix = pred,
  maxit = 20,
  printFlag = TRUE
)

# Qc
print(imp)
stripplot(imp, LVEF ~ .imp, pch = 20, cex = 0.6)

imp$method

########## Cox model (same covariate set as main MV model)

fit_imp <- with(
  imp,
  coxph(Surv(t_followup_days_final, event_inapp_shock) ~
          device_group + age_icd + Sex + LVEF  + strata(dataset))
)

pooled <- pool(fit_imp)
summary_pooled <- summary(pooled, conf.int = TRUE, exponentiate = TRUE)
print(summary_pooled)

summary_pooled[summary_pooled$term == "device_groupSingle", ]

