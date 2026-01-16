########################################################################################################################
########################################################################################################################
# projet: UmBIZO- PROFID (SCD post-MI, age- stratified)
# Script: 01_pipeline_professionnel.R
# Author: Amina Boudamaana
# Date: 2025-10-24
# 
#
# Purpose: 
#       1)  Harmonization of ICD& non-ICD cohorts
#       2)  Data Quality Control (Table 1 bounds + <= 40day rule)
#       3)  Derived variables (binary standardization, *_log1p companions)
#       4)  Age stratification (<= 50, 51-65, 66-75, >75)

# SAP references: Variable selection & Harmonization, Data Quality Control
######################################################"

## ========================================================================================================
## 0. Reproducibility & Session Info
## ========================================================================================================
#file.edit('~/.Renviron')
set.seed(20251224)
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


# Numeric variables we want to coerce when present (union of my columns)
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
# - Reorder columns identically (so rbind is safe)

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
  "ID", "DB", "Age", "Sex", "ICD_status", "LVEF", "LVEF_category","Status","Survival_time", "Time_zero_Y", "Time_zero_Ym", "Baseline_type")

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
    #  are treated as "numeric-like" only if >=70% of values can be parsed as finite numbers.
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
    #This captures clearly right-skewed, positive distribution (NTProBNP, Troponin_T, CRP, Greyzone_size,
   # often TSH/BUN,  and some sacr/volume metrics) while leaving physiologically bounded or clinically  linear variables
   # Age, LVEF, SBP/DBP, HR, PR, QRS, QTc, electrolytes, lipides, eGFR) on the original scale.

hard_exclude  <- c( "Age", "LVEF", "SBP", "DBP", "HR", "PR", "QRS", "QTc",
                    "Sodium", "Haemoglobin", "Cholesterol","HDL", "LDL",
                    "ICD_status", "Status", "Time_index_MI_CHD", "Triglycerides", "TSH", "Infarct_size")


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
       ## Based on the screening metrics, we will log-transform NTProBNP, Troponin_T, CRP, and Greyzone_size, 
       # as each shows marked right-skew (meets >2 of: r_q95_q50 >3, r_max_med >10, Bowley_skew >0.6) and contains no negatives
 
  ##   The SAP requires the inclusion of BMI and eGFR log-transformed versions for model stability, 
   #   regardless of the statistical skewness test results
vars_sap <- c("BMI", "eGFR")

# Merge the automatic selection (cand_vars) with the required variables (SA)
cand_vars_final <- unique(c(cand_vars, vars_sap))

# verification of the final list
cat("Final_list of log-transformed candidates (SAP+skewness): \n")
print(cand_vars_final)
           
# 8) create log1p companions ONLY for the chosen candidates (the final shortlist):
  ## Apply log1p only to variables that (i) exist in 'dat' and (ii) passed the skewness screen('cand_vars').
  # Values <0 are set to NA because log1p is not appropriate for negatives on these clinical scales
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
                  #  [1] "CRP_log1p"   "Greyzone_size_log1p" "NTProBNP_log1p"   "Troponin_T_log1p, "BMI_log1p","eGFR_log1p" 

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
    return(factor(y, levels=c(0,1), labels=c("No", "Yes")))
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
# Apply the standardization ONLY to variables that are present in 'dat' 
for (v in intersect(priority_bin, names(dat))) {
  dat[[v]] <- to_no_yes(dat[[v]])
}
## Check source bin_var presence :  List binary variables now standardized and show quick tables
bin_chk<- sapply(intersect(priority_bin, names(dat)), function(v)is.factor(dat[[v]]) && identical(levels(dat[[v]]), c("No", "Yes")))
print(bin_chk)
## List binary variables now standardized and show quick tables

yesno_vars <- names(dat)[sapply(dat, function(x) is.factor(x) && all(levels(x) %in%c("No", "Yes")))]
print(yesno_vars)

for (v in head(yesno_vars, 5)){cat ("\n==", v, "==\n"); print(table(dat[[v]], useNA = "ifany"))}

## Expected log1p list per SAP+ Swekness
expected <- c("BMI", "eGFR", "NTProBNP", "CRP", "Troponin_T", "Greyzone_size")

## Check source &log1p presence
log_chk<- data.frame(
  var = expected,
  source= expected %in% names(dat), 
  log1p=paste0(expected, "_log1p")%in% names(dat),
  stringsAsFactors = FALSE
)
print(log_chk)
## =================================================================================================================================

##  DATA QUALITY- Table 1 invalide-value bounds 

## =================================================================================================================================
   # for each variable listed in Table 1 , we set clearly inadmissible values to NA on the original clinical scale,
   # Boundry inclusions follow Table 1 exactly 
   # Important: the "measured <=40 days  after MI- rule is applied elsewhere and is intentionally Not handled in this function
## -----------------------------------------------------------------------------------------------------------------------------------
## numeric coercion without stopping on non-numeric  strings
as_num  <- function(x) suppressWarnings(as.numeric(x)) 
# lower/upper cutoffs with control of inclusiveness
## 1A. Inadmissible bounds from Table 1   ##
set_bounds <- function(x, lower=-Inf, upper= Inf,
                       include_lower=TRUE, include_upper=TRUE) {
                         
  v<- as_num(x)
  # lower bound
  if (is.finite(lower)) {
    if (include_lower) v[v <= lower] <- NA else v[v < lower] <- NA
  }
  # upper bound
  if (is.finite(upper)){
    if (include_upper) v[v >= upper] <- NA else v[v> upper] <- NA
  }
  v
}
## Apply Table 1 bounds only if columns  exist
apply_table1_bounds <- function(df) {
 
      ## BMI< 12 or >69 -> NA 
  if ("BMI" %in% names(df)) {
    v<- as_num(df$BMI); v[v <12 | v> 69] <- NA; df$BMI <-v
  }
     ## BUN >900 mmo/L --> NA
  if ("BUN" %in% names(df)){  
    df$BUN <- set_bounds(df$BUN, upper = 900, include_upper= FALSE)
    } 
    ## Cholesterol <= 50 mg/dL -> NA
  if ("Cholesterol"  %in% names(df)) {
    df$Cholesterol <- set_bounds(df$Cholesterol, lower = 50, include_lower= TRUE)
  } 
   ## Haemoglobin <2 or >110 g/dL -> NA
  if ("Haemoglobin" %in% names(df)) 
    v<- as_num(df$Haemoglobin); v[v < 2 | v> 110] <- NA; df$Haemoglobin <-v
  ## HbA1c < 2.5% -> NA 
  if ("HbA1c"  %in% names(df)) {
    df$HbA1c <- set_bounds(df$HbA1c, lower = 2.5, include_lower= FALSE)
  }
  ## HDL 0 mg/dL -> NA (<=0)
  if ("HDL"  %in% names(df)) {
    df$HDL <- set_bounds(df$HDL, lower = 0, include_lower= TRUE)
  } 
  ## LDL 0 mg/dL --> NA (<= 0)
  if ("LDL" %in% names(df)) {
    df$LDL <- set_bounds(df$LDL, lower = 0, include_lower= TRUE)
  }
  ## Sodium < 99 mmol/L ->NA 
  if ("Sodium" %in% names(df)) {
    df$Sodium <- set_bounds(df$Sodium, lower = 99, include_lower= FALSE)
  } 
  ## Triglycerides <20 mg/dL -> NA (20 is OK)
  if ("Triglycerides" %in% names(df)) {
    df$Triglycerides <- set_bounds(df$Triglycerides, lower = 20, include_lower= FALSE)
  }
  ## TSH 0 mU/L -> NA (<=0)
  if ("TSH" %in%  names(df)){
    df$TSH <- set_bounds(df$TSH, lower = 0, include_lower= TRUE)
  } 
  ## HR < 25 or >140 bpm -> NA  
  if ("HR" %in% names(df)) {
    v<- as_num(df$HR); v[v < 25 | v> 140] <- NA; df$HR<-v
 }
  ## PR < 50 or >1000 ms -> NA
  if ("PR" %in% names(df)) {
    v<- as_num(df$PR); v[v < 50 | v> 1000] <- NA; df$PR<-v
 }
  ## QRS <50 ms-> NA (50 Okay)
  if ("QRS" %in% names(df)) {
    df$QRS <- set_bounds(df$QRS, lower = 50, include_lower= FALSE)
  }
  ## QTc <= 250 or >790 ms -> NA
  if ("QTc" %in% names(df)) {
    v<- as_num(df$QTc); v[v <= 250 | v> 790] <- NA; df$QTC<-v
  }
  df
}
# Save pre-Quality Controle (pre-QC) copy for audit
preQC <- dat
# Apply table 1 bounds
dat <- apply_table1_bounds(dat)
# Recompute log1p companions after QC
relog1p <- function(df, vars){
  for (v in intersect (vars, names(df))) {
    x<- suppressWarnings(as.numeric(df[[v]]))
    x[x <0] <- NA # log undefined for negatives
    df[[paste0(v, "_log1p")]]<- log1p(x)
  }  
  df
}
dat<- relog1p(dat, c("BMI", "eGFR", "NTProBNP", "CRP", "Troponin_T", "Greyzone_size"))
 
# Audit: How many values were set to NA by QC??
na_delta <- function (before, after)sum(is.na(after) & !is.na(before))
vars_chek <- c("BMI", "BUN","Cholesterol", "Haemoglobin", "HbA1c", "HDL","LDL","Sodium", "Triglycerides","TSH","HR", "PR", "QRS", "QTc" )
audit <- do.call(rbind, lapply(vars_chek, function(v){
  if (!v %in% names(dat)) return(NULL)
  b<- suppressWarnings(as.numeric(preQC[[v]]))
  a<- suppressWarnings(as.numeric(dat[[v]]))
  data.frame(variable=v,
             n_set_to_NA= na_delta(b,a),
             n_nonmiss_post=sum(is.finite(a)))
}))
cat("Table-1 audit (values set to NA by QC):\n"); print(audit, row.names =FALSE)


## =================================================================================================================================
## Rule: Set measurements taken <=40 days after MI to NA
## Applies to: SBP, DBP, CRP,Troponin_T, NYHA, AV_blok, AV_blok_II_III
## Time source supported ( OR logic):
#                                  1/ "Time_index_MI_CHD" (days since MI to Baseline)
#                                  2/ Time_zero_Y" (years since MI) -< convert to days
#                                  3/ "Time_zero_Ym " (Years since MI -< convert to day)
#                                  4/ " Baseline_type" ("pre-40", "within 40", "early)
## We audit how many values are set to NA and recompute log1p only for CRP/Troponin_T 
## =================================================================================================================================
# 1/ Robust flag: measurement taken <=40 days from MI
within_40d_flag <- function(df){
  n<- nrow(df)
  flag<- rep(FALSE, n)
  
# A/ Use Baseline_type : match 'before', 'pre-40', '<=40', <40', 'within 40', 'early' 
  if("Baseline_type" %in% names(df)){
    bt<- tolower(trimws(as.character(df$Baseline_type)))
    bt<- gsub("\\s+", " ", bt) ## normalize spaces##
   ## 'pre' only when followed by -40 or 40 (avoid matching 'preserved')
    pat<- "(before|pre[ -]?40|<=\\s*40|<\\s*40|within\\s*40|within40|early|<=\\s*40)"
        flag<- flag | grepl(pat,bt)
  }
# B/ Use numeric time sine MI in years(Time_zero_Y)(covert threshold into years =40/365.25 years)
  year_threshold <- 40/365.25
  
  if("Time_zero_Y" %in% names(df)){
    y<- suppressWarnings(as.numeric(df$Time_Y))
    flag<- flag | (is.na(y) & y< year_threshold)
  }

# C/ If we also have another "years since MI" fiels  => same threshold
  if("Time_zero_Ym" %in% names(df)){
    ym<- suppressWarnings(as.numeric(df$Time_zero_Ym))
    flag<- flag | (is.na(ym) & ym <year_threshold)
  }
# D/ Explicit days since MI(Time_index_MI_CHD) last  
  if("Time_index_MI_CHD" %in% names(df)) {
    ty<- suppressWarnings(as.numeric(df$Time_index_MI_CHD))
    flag<- flag | (is.na(ty) & ty<year_threshold)
  }
  flag
}
flag40 <- within_40d_flag(dat)
cat("Rows flagged as <=40 days:", sum(flag40), "out of", nrow(dat), "\n")

## 2/ Variables affected by the <=40 days rule (Table1) 
vars_40d_all <- c("SBP", "DBP", "CRP", "Troponin_T", "NYHA", "AV_block", "AV_block_II_or_III")
vars_40d<-intersect(vars_40d_all, names(dat)) ### keep only those that exist


##NYHA harmonisation (SAP) 
## Convert various encodings to an ordered factor: I < II < III < IV
standardise_nyha <- function(x){
  z <- tolower(trimws(as.character(x)))
  # remove words/separators and collapse spaces
  z <- gsub("class|nyha|:|;|,|-|_", "", z)
  z <- gsub("\\s+", "", z)
  # keep only roman/digit patterns
  z <- gsub("[^ivx\\d]", "", z)
  
  # map numeric & roman variants to canonical levels
  z[z %in% c("1","i")] <- "I"
  z[z %in% c("2","ii")] <- "II"
  z[z %in% c("3","iii")] <- "III"
  z[z %in% c("4","iv")] <- "IV"
  
  # anything else -> NA
  z[!z %in% c("I","II","III","IV")] <- NA
  
  factor(z, levels = c("I","II","III","IV"), ordered = TRUE)
}

if ("NYHA" %in% names(dat)) {
  dat$NYHA <- standardise_nyha(dat$NYHA)
  cat("NYHA levels after standardisation:\n")
  print(table(dat$NYHA, useNA = "ifany"))
}
## ----------------------------------------------------

## 3/ Count how many values will be set to NA (pre-application) 
count_40d_na <- function(df, var, flag){
  x<- df[[var]]
  data.frame(
    variable        =var,
    nonNA_before    =sum(!is.na(x)),
    stringsAsFactors = FALSE
  )
}
report_40d <- do.call(rbind, lapply(vars_40d, function(v) count_40d_na(dat, v, flag40)))
print(report_40d)     # variable    nonNA_before
     #                 SBP          8906
    #                  DBP         8904
    #                  CRP         6416
    #              Troponin_T      6216
    #                  NYHA        15008
    #             AV_block         8010
    #   AV_block_II_or_III         8254

## 4/ Apply the <= 40 days rule: set values to NA where flag40 is TRUE
for(v in vars_40d){
  dat[[v]][flag40] <-NA
}
## 5/ Sanity check after application 
report_40d_after<- data.frame(
  variable  = vars_40d,
  nonNA_after =sapply(vars_40d, function(v) sum(!is.na(dat[[v]]))),
  stringsAsFactors = FALSE
)
print(report_40d_after)      # variable    nonNA_after
#                 SBP          8906
#                  DBP         8904
#                  CRP         6416
#              Troponin_T      6216
#                  NYHA        15008
#             AV_block         8010
#   AV_block_II_or_III         8254

## Recompute *_log1p companions  
relog1p<- function(df, vars){
  for(v in intersect( vars, names(df))){
    x<- suppressWarnings(as.numeric(df[[v]]))
    x[x<0] <- NA
    df[[paste0(v, "_log1p")]] <- log1p(x)
  }
  df
}
dat<- relog1p(dat, c("CRP", "Troponin_T")) # We only recomputed CRP_log1p and Troponin_T_log1p because 
                                           # The <= 40 days rule affects only these two variables among those with log transformation. 
                                           # Recomputing the other _log1p variables isn't necessary

summary(dat$Time_index_MI_CHD)    # Min=1314 (years), so the shortest delay is =480 days (>40 days)
sum(flag40, na.rm = TRUE)         # ==0 is consistent: no measurements <= 40 days, so nothing needs to be set to NA under this rule
## =================================================================================================================================
##   7 x IQR  outlier screening (report only, do deletion here)
## APPLY ONLY TO CONTINUOUS VARIABLES WITHOUT EXPLICIT CLINICAL CUT-OFFS
## WE compute Q1/Q3/IQR on numeric , non-missing values and report how many observations fall below/ above the 7*IQR fences. 

## =================================================================================================================================
iqr_report_one<- function(x) {
  x<- suppressWarnings(as.numeric(x))
  x<- x[!is.finite(x)]     # we keep numeric, non-missing
  # skip very small samples (robustness)!
  if (length(x) <20) 
    return(c(n=length(x), 
             q1=NA, q3=NA, iqr=NA,
             low=NA, up=NA, n_low=NA, n_up=NA))
  qs<- quantile(x, c(.25, .75), na.rm = TRUE, names = FALSE)
  q1<- qs[1]; q3<- qs[2]
  iqr<- q3-q1
  low<- q1-7*iqr
  up<- q3+7*iqr
  c(n=length(x), 
    q1=q1, q3=q3, iqr=iqr,
    low=low, up=up,
    n_low=sum(x <low),
    n_up=sum(x>up))
}

## A/ 7*IQR outlier report( no deletion yet) 
# Choose continuous vars likely to benefit from 7*IQR rule (no clinical cutoffs) 
iqr_vars <- intersect(c(
  "Total_scar", "NTProBNP"  # labs with huge skew but no explicit clinical bound in Table  1 NTProBNP
), names(dat))

#  Compute 7*IQR fences on numeric, non-missing values  
iqr_report<- do.call(rbind, lapply(iqr_vars,function(v){
  out<- iqr_report_one (dat[[v]])
  data.frame(variable=v, t(out), row.names = NULL, check.names = FALSE)
}))
print(iqr_report)
# Recompute log1p only for biomarkers affected by a change (if any) 
if("NTProBNP" %in% names(dat)){
  x<- suppressWarnings(as.numeric(dat$NTProBNP))
  x[x<0] <- NA
  dat$NTProBNP_log1p<- log1p(x)
}
##Quick sanity check ##
summary(dat$NTProBNP_log1p)
 #  NTProBNP_log1p summary (after recomputation):
 # Min= 0.40, Q1= 2.80, Median=3.75, Q3=4.84, Max=9.26
 # NA'S=133.220 -this reflects limited availability for NTProBNP in the cohort and is NOT due to the 7*IQR screening (no values were set to NA)
## =================================================================================================================================
##        AGE STRATIFICATION
# <= 50 = young adults
# 51-65 = middle-aged
# 66-75 = older adults
# >75   = elderly
## =================================================================================================================================
# Create the 4 age cohorts per SAP: 
dat$age_group<- cut(as.numeric(dat$Age),
                    breaks=c(-Inf,50, 65, 75, Inf),
                      labels=c("<=50", "51-65", "66-75", ">75"),
                    right=TRUE)

# Ensure a consistent ordered factor( same order everywhere)
dat$age_group <- factor(dat$age_group,
                        levels = c("<=50", "51-65", "66-75", ">75"),
                        ordered = TRUE)

# A descriptive label version 
dat$age_group_desc <- factor(dat$age_group,
                        levels = c("<=50", "51-65", "66-75", ">75"),
                        labels = c("<=50  (young adults)",
                                   "51-65 (middle-aged)",
                                   "66-75 (older adults)", 
                                   ">75   (elderly)"),
                        ordered = TRUE)

                        
## Sanity check 
stopifnot(sum(!is.na(dat$Age)) == sum(!is.na(dat$age_group)))
print(table(dat$age_group, dat$ICD_status,     useNA = "ifany")) #  ICD_status             0      1
                                                                  # young-adult "<=50"     9507   845
                                                                  # middle-aged "51-65"    39728  3061
                                                                  # Older adults "66-75"   40063  2558
                                                                  # elderly      ">75"     43363  1077
                                                                  # <NA>                   0      2
print(table(dat$age_group, dat$LVEF_category,  useNA = "ifany"))
  ##  LVEF_category        Reduced Preserved
 #     young-adult "<=50"     1922      8430
 #     middle-aged  "51-65"   8532      34257
 #     Older-adults "66-75"   9705      32916
 #    elderly       ">75"     12440     32000
 #    <NA>                    2         0
