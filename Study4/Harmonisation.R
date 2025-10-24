########################################################################################################################
########################################################################################################################
# projet: UmBIZO- PROFID (SCD post-MI, age- stratified)
# Script: 01_pipeline_professionnel.R
# Author: Amina Boudamaana
# Date: 2025-10-17
# 
#
# Purpose: 
#       1)  Harmonisation of ICD& non-ICD cohorts
#       2)  Derived variables (binary standardisation, *_log1p companions)
#    

# SAP references: Variable selection & Harmonisation
######################################################"

## ========================================================================================================
## 0. Reproducibility & Session Info
## ========================================================================================================
#file.edit('~/.Renviron')

RStudio.Version()$version ##  �e2024.9.0.375�f
R.version.string #  " version 4.4.1 (2024-06-14 ucrt)"

set.seed(20251017)

## ========================================================================================================
## 1. Packages
## ========================================================================================================
# install. packages
install.packages(c("tidyverse","janitor","lubridate","survival","cmprsk","mice","survminer", "simkid", "stringr", "readr"))
library(tidyverse)
library(dplyr)
library(janitor)
library(stringr)
library(readr)

## ========================================================================================================
## 2. Data import
## ========================================================================================================
getwd()

ICD_all          <- read.csv("ICD_all.csv")  #(Reference)
ICD              <- read.csv("ICD.csv")
NonICD_preserved <-read.csv("NonICD_preserved.csv")
NonICD_reduced   <-read.csv("NonICD_reduced.csv")


cat("ICD:",              nrow(ICD),             "rows | ", ncol(ICD), "cols\n",              sep= "")  ## ICD:7543 rows | 75 cols
cat("NonICD_preserved:", nrow(NonICD_preserved),"rows | ", ncol(NonICD_preserved), "cols\n", sep= "")  ## NonICD_preserved:107603 rows | 88cols
cat("NonICD_reduced:",   nrow(NonICD_reduced),  "rows | ", ncol(NonICD_reduced), "cols\n",   sep= "")  ## NonICD_reduced:25058rows | 88cols


str(ICD) ; str(NonICD_preserved); str(NonICD_reduced)
## ========================================================================================================
## 3. Core derived variables (cohort flags)
## ========================================================================================================
# ICD_status :   identifies implanted-ICD registry vs non-ICD cohorts (selection bias later)
# LVEF_category: cohort-defined preserved ">35%" vs reduced "<=35%" ( used for descriptive stratification)
ICD$ICD_status               <- 1
NonICD_preserved$ICD_status  <- 0
NonICD_reduced$ICD_status    <- 0


ICD$LVEF_category              <- "Reduced"     # ICD cohort (post-IM, reduced EF)

NonICD_preserved$LVEF_category <- "Preserved"   # Non-ICD preserved EF

NonICD_reduced$LVEF_category   <- "Reduced"     # Non-ICD reduced EF


# Clean factors with stable ordering for reporting 
ICD$LVEF_category               <- as.factor(ICD$LVEF_category)
NonICD_preserved$LVEF_category  <- as.factor(NonICD_preserved$LVEF_category)
NonICD_reduced$LVEF_category    <- as.factor(NonICD_reduced$LVEF_category)

## ========================================================================================================
## 4. Harmonisation 
## ========================================================================================================
as_num<- function(x) suppressWarnings(as.numeric(x))

stopifnot(exists("ICD"), exists("NonICD_preserved"), exists("NonICD_reduced"))


# Numeric variables we wante to coerce when present (union of my columns)
numeric_vars<- c("Age","Survival_time", "Time_zero_Y", "Time_zero_Ym", "LVEF", "MRI_LVEF", "LVEDV", 
                 "LVESV","LV_mass", "SBP", "DBP", "HR", "PR", "QRS", "QTc", "BMI", "eGFR", "BUN", 
                 "Sodium", "Potassium", "Haemoglobin", "CRP", "Troponin_T", "NTProBNP", "Cholesterol", 
                 "HDL", "LDL", "Triglycerides", "Infarct_size", "Total_scar","Greyzone_size"
)

force_numeric<- function(df, var){
  for(v in intersect(var, names(df))) df[[v]] <- as_num(df[[v]])
  df
}

#  Align columns across datasets:
# - Add any missing columns (as NA)
# - Reorder columns identifically (so rbin is safe)

align_to <- function(df, cols){
  miss<- setdiff(cols, names(df))
  for(m in miss) df[[m]]<- NA
  df[, cols, drop=FALSE]
}

#  Build the "union" of existing columns + must have (ICD_status, Lvef_category)
present_all  <- sort(Reduce(union, list(names(ICD), names(NonICD_preserved), names(NonICD_reduced))))
must_have    <- c ("ICD_status", "LVEF_category")
final_col    <- sort(unique(c(present_all, must_have)))

## ========================================================================================================
## 5. Coercions +column alignment + Merge
## ========================================================================================================

ICD               <- force_numeric(ICD,numeric_vars)
NonICD_preserved  <- force_numeric(NonICD_preserved, numeric_vars)
NonICD_reduced    <- force_numeric(NonICD_reduced, numeric_vars)

ICD_align  <- align_to (ICD, final_col)
NP_align   <- align_to(NonICD_preserved, final_col)
NR_align   <- align_to(NonICD_reduced,final_col)

# 4/ Safe row-bind after identical columns order
dat <- rbind(ICD_align,NP_align, NR_align)

# Sanity  checks 
stopifnot(identical(names(ICD_align) , names(NP_align)))
stopifnot(identical(names(ICD_align) , names(NR_align)))

identical(names(ICD_align), names(NP_align))
identical(names(ICD_align), names(NR_align))
cat("Merged rows", nrow(dat),  "| variables", ncol(dat) ,"\n")  # dat==> 140204 obs. 93 variables

## ========================================================================================================
## 6.  Column order 
## ========================================================================================================

move_to  <-  function(df, cols_front) {
  front    <-  unique(intersect(cols_front, names(df)))
  rest     <-  setdiff(names(df), front)
  df[, c(front, rest), drop=FALSE]
}

#  Define the order
preferred_front  <- c(
  "ID", "DB", "Age", "Sex", "ICD_status", "LVEF", "LVEF_category","Status","Survival_time", "Time_zero_Y", "Time_zeo_Ym", "baseline_type")

dat   <- move_to(dat,preferred_front)
table(ICD$LVEF_category,                 useNA = "ifany")
table(NonICD_preserved$LVEF_category,    useNA = "ifany")
table(NonICD_reduced$LVEF_category,      useNA = "ifany")
table(dat$ICD_status, dat$LVEF_category, useNA = "ifany")
unique(dat$LVEF_category)

print(table(dat$LVEF_category,useNA = "ifany"))
print(table(dat$ICD_status, useNA = "ifany"))
by(dat$LVEF_category, dat$ICD_status, function(x) head(x,3))
set.seed(1); print(sample(dat$LVEF_category, 20))



## ========================================================================================================
## 7.  Choose which variables to log1p
## ========================================================================================================
# Rationale for log1p candidate selection:
## We do NOT log-tranform every numeric field, Instead, we add *_log1p only for continuous variables
#  that are strictly non-negative and show clear right-skew,where a multiplicative scale is clinically sensible and improves model fit!


#  1) Safely coerce to numeric:
## We parse each column to numeric with suppressWarnings(); character columns
#  are treated as "numeric-like" only if >=70% of values can be parded as finite numbers.
#  This avoids spuriously analysing  codes/labels as numbers.
as_num<-function(x) suppressWarnings(as.numeric(x))


# 2) Quartile (Bowly) skewness: robuste to outliers:
##We compute Bowley (quartile) skewness:
# (Q3+Q1-2*Q2)/(Q3-Q1), which depends on quartiles (Q1,Q2,Q3) and is  robust to outliers.
bowley_skew <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) <20) return(NA_real_)
  qs <- quantile(x, c(.25, .5, .75),  na.rm = TRUE, names = FALSE)
  if (qs[3] ==qs[1]) return(0)
  (qs[3]+ qs[1] - 2*qs[2])  / (qs[3] -qs[1])
}

# 3) Build a metrics table for numeric-like variables:

numeric_like <- names(dat)[sapply(dat, function(z){
  # Consider numeric /integer or character that can be mostly parsed as numeric
  is.numeric(z) || is.integer(z) || (is.character(z) && {
    tmp  <- suppressWarnings(as.numeric(z)); mean(is.finite(tmp), na.rm=TRUE) >0.7
  })
})]

metrics  <- do.call(rbind, lapply(numeric_like, function(v){
  x <- as_num(dat[[v]])
  finite  <- is.finite(x)
  x2 <- x[finite]
  if (length(x2) ==0) return(NULL)
  q <- quantile(x2, c(.25,.5,.75,.95,.99), names=FALSE)
  data.frame(
    variable= v,
    n_nonmiss       = sum(!is.na(x)),                      # count of no-missing values 
    prop_finite     = mean(finite),                  
    prop_negative   = mean(x2 <0),                        # proportion<0 (log not defined -> should  be 0 )
    prop_zero       = mean(x2 ==0),                       # proportion =0 (log1p handles 0)
    min             = min(x2),
    median          = q[2],                               # scale/upper-tail summary
    max             = max(x2),                            # scale/upper-tail summary
    q1              = q[1],
    q3              = q[3],
    q95             = q[4],                               # scale/upper-tail summary
    q99             = q[5],
    bowley_skew     = bowley_skew(x2),                    # robust asymmetry
    r_max_med       = ifelse(q[2] >0, max(x2)/ q[2], Inf), # max/median (extreme tail vs center)
    r_q95_q50       = ifelse(q[2] >0, q[4]/q[2], Inf),    # upper-tail vs central tendency
    stringsAsFactors= FALSE 
  )
}))

cand_report <- metrics
print(head(cand_report, 50), digits = 3)  # Only 38 variables metrics
# 4) Heuristic decision rule: Hard-exclude variables where log is rarely useful 
## A variable is flagged for log1p if: 
# prop_negative ==0 AND 
# at least two of {bowley_skew> 0.6, r_q95_q50 >3, r_max_med >10}
#This captures clearly right-skewed, positive distribution (NTProBNP, Troponin_T, CRP, Triglycerides, 
# often TSH/BUN,  and some sacr/volume metrics) while leaving physiologically bounded or clinically  linear variables
# Age, LVEF, SBP/DBP, HR, PR, QRS, QTc, electrolytes, lipides, eGFR) on the original scale.

hard_exclude  <- c( "Age", "LVEF", "SBP", "DBP", "HR", "PR", "QRS", "QTc",
                    "Sodium", "Haemoglobin", "Cholesterol","HDL", "LDL",
                    "ICD_status", "Status", "Time_index_MI_CHD", "Triglycerides", "TSH")


# 5) Heuristic to flag "very right-skewed & positive" candidates
cand <- subset(metrics,
               prop_negative <= 0.01 &  # should be >= 0 for log1p
                 (bowley_skew   > 0.6 |      # strongly right-skewed
                    r_q95_q50    > 3   |      # 95th >> median
                    r_max_med    > 10))       # extreme range
cand_vars <- intersect(trimws(cand$variable), names(dat))

# 6) Remove hard-excluded variables and keep only columns that exist in dat
cand_vars <- setdiff(cand_vars, hard_exclude)


# 7) Report top candidates sorted by skewness

cand_report <- cand[order(-cand$bowley_skew), c("variable", "n_nonmiss", "prop_zero",
                                                "bowley_skew", "r_q95_q50", "r_max_med", "median", "q95", "max")]
print(head(cand_report, 20),  digits=3)
## Based on the screening metrics, we will log-trabsform NTProBNP, Troponin_T, CRP, and Greyzone_size, 
# as each shows marked right-skew (meetss >2 of: r_q95_q50 >3, r_max_med >10, Bowley_skew >0.6) and contains no negatives

## The SAP requires the inclusion of BMI and eGFR log-transformed versions for model stability, 
#   regardless of the statistical skewness test results
vars_sap <- c("BMI", "eGFR")

# Merge the automatic selection (cand_vars) with the required variables (SA)
cand_vars_final <- unique(c(cand_vars, vars_sap))

# verification of the final list
cat("Final_list of log-transformed candidates (SAP+skewness): \n")
print(cand_vars_final)

# 8) create log1p companions ONLY for the chosen candidates (the final shortlist):
## Apply log1p only to variables that (i) existe in 'dat' and (ii) passed the skewness screen('cand_vars').
# Values <0 are set to NA because log1p is not appropriate dor negatives on these clinical scales
# The function returns 'dat' with new columns named <var>_log1p


add_log1p <- function (df, vars){
  for(v in intersect(vars, names(df))){
    x<- as_num(df[[v]]); x[x<  0]<- NA     ## log1p undefined for negatives
    df[[paste0(v, "_log1p")]]<- log1p(x)
  }
  df
}

dat  <- add_log1p(dat, cand_vars_final)
# 9) Sanity check: list of log columns created:

created_logs  <- paste0(cand_vars_final, "_log1p")
created_logs  <- created_logs[created_logs %in% names(dat)]
cat ("Created log1p variables:\n"); print(created_logs) 
# Created log1p variables:
#  [1] "CRP_log1p"   "Greyzone_size_log1p" "NTProBNP_log1p"   "Troponin_T_log1p, "BMI_log1p","eGFR_log1p" "Infract_size_log1p" 

## ========================================================================================================
## 8.  Standardise binary indicators to an ordered factor {No, Yes}
## ========================================================================================================

# Different cohorts encode binaries in many ways (0/1, TRUE/FALSE),
# "yes"/"no", "Y/N") For modelling and reporting we harmonise them
# to a single, ordered factor with levels c("No", "Yes")

## - Preserves NA values- Converts safely when values are clearly binary; otherwise leaves the vector as a factor (fallback) 
# so we do not silently recode non-binary values

to_no_yes <- function (x) {
  ## Standarise to factor ("No", "Yes")
  if (is.factor(x)) x <- as.character(x)
  
  ## 1/ Numeric/integer 0/1 already
  if (is.numeric(x)|| is.integer(x)) {
    if (all(na.omit(x)%in% c(0,1))) {  
      return(factor(x, levels = c(0,1), labels = c("No", "Yes")))
    }
  }
  ## 2/ Character "0"/"1"
  if (is.character(x) && all(na.omit(x)%in% c("0","1"))) {
    y <- as.numeric(x)
    return(factory(y, levels=c(0,1), labels=c("No", "Yes")))
  }
  ## 3/ Logical TRUE/FALSE
  if (is.logical(x)) {
    return(factor(ifelse(x, "Yes", "No"), levels=c("No", "Yes")))   
  }
  ## 4/ Character "yes/no" variants 
  if (is.character(x)) {
    z <- tolower(trimws(x))
    z[z %in% c("y", "yes", "true")]<- "yes"
    z[z %in% c("n", "no", "false")]<- "no"
    out <- ifelse(z== "yes", "Yes", ifelse(z=="no", "No", NA))
    return(factor(out, levels = c("No", "Yes")))
  }
  ##fallback:  we return as factor without forcing a binary mapping 
  return(as.factor(x))
}
# Apply the standardization ONLY to variables that are present in 'dat' 
for (v in intersect(priority_bin, names(dat))) {
  dat[[v]] <- to_no_yes(dat[[v]])
}

## List binary variables now standardized and show quick tables

yesno_vars <- names(dat)[sapply(dat, function(x) is.factor(x) && all(levels(x) %in%c("No", "Yes")))]
print(yesno_vars)
for (v in head(yesno_vars, 5)){cat ("\n==", v, "==\n"); print(table(dat[[v]], useNA = "ifany"))}

## Priority binaires to standardise if they exist in data##
priority_bin <-c(
  "Diabetes","Hypertension","HF", "COPD","Stroke_TIA","Smoking","FH_CAD","FH_SCD",
  "AV_block","AV_blokc_II_or_III", "RBBB","LBBB",
  ## Medications
  "ACE_inhibitor","ARB","ACE_inhibitor_ARB","Beta_blockers","Diuretics","Calcium_antagonistes",
  "Aldosterone_antagonist","Anti_arrhythmic_III","Anti_coagulant","Anti_platelet","Lipid_lowering",
  "Anti_diabtetic","Anti_diabetic_oral","Anti_diabetic_insulin",
  
  "CABG", "PCI", "Thromolysis_acute","Revascularisation_acute","CABG_acute","PCI_acute","HaMRI"
)
for (v in intersect(priority_bin, names(dat))) {
  dat[[v]] <- to_no_yes(dat[[v]])
}

## Expected log1p list per SAP+ Swekness
expected <- c("BMI", "eGFR", "NTProBNP", "CRP", "Troponin_T", "Greyzone_size","Infarct_size")

## Check source &log1p presence
chk<- data.frame(
  var = expected,
  source= expected %in% names(dat), 
  log1p=paste0(expected, "_log1p")%in% names(dat),
  stringsAsFactors = FALSE
)
print(chk)

