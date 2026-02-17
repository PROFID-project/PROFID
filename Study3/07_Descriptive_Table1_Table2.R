
#######  BASELINE DESCRIPTIVES & DEVICE COMPARISON #######

#rm(list = ls())
# #
# install.packages("readxl")
# install.packages("data.table")
# install.packages("tableone")
# install.packages("psych")
# install.packages("stats")


library(readxl)
library(data.table)
library(tableone)
library(psych)
library(stats)

####### LOAD FINAL MERGED DATASET #######
 
dt <- readRDS("study3_final_merged.rds")

####### COHORT 1: BASELINE DESCRIPTIVES (ALL ELIGIBLE DATASETS)

dt_desc <- dt[
  dataset %in% c("EUCERT", "HELIOS", "PROSE", "ISRAEL")  
]
names(dt_desc)

####### DEFINE DEVICE GROUP (SINGLE VS DUAL) #######

dt_desc[, device_group := NA_character_]

dt_desc[grepl("\\b(single|vvi|icd_1)\\b", tolower(icd_type)),
         device_group := "Single"]

dt_desc[grepl("\\b(dual|ddd|icd_2)\\b", tolower(icd_type)),
  device_group := "Dual"]

# QC: unclassified
dt_desc[is.na(device_group), .N, by = dataset]

# Restrict to single vs dual only
dt_desc <- dt_desc[!is.na(device_group)]
table(dt_desc$device_group)
table(dt_desc$Hypertension, dt_desc$device_group)
table(dt_desc$dataset)
####### VARIABLES FOR TABLE 1 #######

vars_continuous <- c("age_icd","BMI","LVEF","eGFR")

vars_categorical <- c("Sex","NYHA","AF_atrial_flutter", "Diabetes",
  "Hypertension", "COPD", "Stroke_TIA", "Smoking",  "ACE_inhibitor_ARB", "Beta_blockers",
  "Anti_arrhythmic_III","Anti_platelet","Anti_coagulant")

vars_table1  <- c(vars_continuous, vars_categorical)

####### NORMALISE YES/NO FIELDS #######
yn_keep_na <- function(x){
  x <- tolower(trimws(as.character(x)))
  fifelse(
    x %in% c("yes","y","1","true","t"), "Yes",
    fifelse(x %in% c("no","n","0","false","f"), "No", NA_character_)
  )
}

# Apply to all med/comorbidity yes/no fields  (keeps NA as NA)
yn_vars <- intersect(
  c("AF_atrial_flutter","Diabetes","Hypertension","Stroke_TIA","Smoking",
    "ACE_inhibitor_ARB","Beta_blockers","Anti_arrhythmic_III","Anti_platelet","Anti_coagulant"),
  names(dt_desc)
)
dt_desc[, (yn_vars) := lapply(.SD, yn_keep_na), .SDcols = yn_vars]

####### SPECIAL RULE: recode NA -> "No" ONLY for meds with partial missingness #######
med_recode_vars <- intersect(
  c("ACE_inhibitor_ARB","Beta_blockers","Anti_platelet",
    "Anti_arrhythmic_III","Anti_coagulant"),
  names(dt_desc)
)

for (v in med_recode_vars) {
  
  # registries where the variable is available (i.e., not 100% missing)
  avail_ds <- dt_desc[, .(all_na = all(is.na(get(v)))), by = dataset][all_na == FALSE, dataset]
  
  # only within those registries: patient-level NA => "No"
  dt_desc[dataset %in% avail_ds & is.na(get(v)), (v) := "No"]
}
dt_desc[, .(n_yes=sum(get("Anti_coagulant")=="Yes", na.rm=TRUE),
            n_no =sum(get("Anti_coagulant")=="No",  na.rm=TRUE),
            n_na =sum(is.na(get("Anti_coagulant")))),
        by=dataset]


####### TABLE 1 — BASELINE DESCRIPTIVES (OVERALL + BY DEVICE) #######

table1 <- CreateTableOne(
  vars        = vars_table1,
  strata      = "device_group",
  data        = dt_desc,
  factorVars  = vars_categorical,
  includeNA  = TRUE,
  addOverall = TRUE
)

print(
  table1,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE
)


####### COHORT 2: DEVICE COMPARISON / INFERENTIAL ANALYSES

dt_desc[, t_followup_days_final := NA_real_]

# Israel-specific follow-up (if present)
if ("t_followup_days_israel" %in% names(dt_desc)) {
  dt_desc[dataset == "ISRAEL",
          t_followup_days_final := t_followup_days_israel]
}

# All other datasets
dt_desc[dataset != "ISRAEL",
        t_followup_days_final := t_followup_days]

dt_analysis <- dt_desc   


#######  NORMALITY (PER DEVICE GROUP) #######

normality_results <- rbindlist(
  lapply(vars_continuous, function(v) {
    dt_analysis[!is.na(get(v)),
                .(
                  variable = v,
                  p_shapiro = if (.N >= 3 & .N <= 5000)
                    shapiro.test(get(v))$p.value
                  else NA_real_
                ),
                by = device_group
    ]
  })
)

normality_results

non_normal_vars <- unique(
  normality_results[p_shapiro < 0.05, variable]
)

non_normal_vars

####### CONTINUOUS VARIABLE TESTS (T-TEST VS WILCOXON) #######

continuous_tests <- rbindlist(
  lapply(vars_continuous, function(v) {
    
    x <- dt_analysis[device_group == "Single", get(v)]
    y <- dt_analysis[device_group == "Dual",   get(v)]
    
    if (v %in% non_normal_vars) {
      test <- wilcox.test(x, y)
      test_type <- "Wilcoxon rank-sum"
      stat <- unname(test$statistic)
    } else {
      test <- t.test(x, y)
      test_type <- "Student t-test"
      stat <- unname(test$statistic)
    }
    
    data.table(
      variable = v,
      test = test_type,
      statistic = stat,
      p_value = test$p.value
    )
  })
)

continuous_tests

####### CATEGORICAL VARIABLE TESTS (CHI-SQUARE / FISHER) #######

categorical_tests <- rbindlist(
  lapply(vars_categorical, function(v) {
    
    sub <- dt_analysis[!is.na(get(v)) & !is.na(device_group)]
    
    # If after filtering there is no data, skip
    if (nrow(sub) == 0) {
      return(data.table(
        variable = v,
        test = NA_character_,
        chi_square = NA_real_,
        p_value = NA_real_
      ))
    }
    
    tab <- table(sub[[v]], sub$device_group)
    
    # If table has any zero margins, skip safely
    if (any(dim(tab) < 2) || all(tab == 0)) {
      return(data.table(
        variable = v,
        test = NA_character_,
        chi_square = NA_real_,
        p_value = NA_real_
      ))
    }
    
    if (any(tab < 5)) {
      test <- fisher.test(tab)
      data.table(
        variable = v,
        test = "Fisher exact",
        chi_square = NA_real_,
        p_value = test$p.value
      )
    } else {
      test <- chisq.test(tab)
      data.table(
        variable = v,
        test = "Chi-square",
        chi_square = unname(test$statistic),
        p_value = test$p.value
      )
    }
  })
)



categorical_tests

####### FOLLOW-UP SUMMARY (YEARS) #######

followup_summary <- dt_analysis[
  ,
  .(
    n = .N,
    median_years = median(t_followup_days_final, na.rm = TRUE) / 365.25,
    iqr_years = IQR(t_followup_days_final, na.rm = TRUE) / 365.25
  ),
  by = device_group
]

followup_summary

####### MISSINGNESS BY DEVICE GROUP (%) #######

missing_by_group <- dt_analysis[
  ,
  lapply(.SD, function(x) round(mean(is.na(x)) * 100, 2)),
  by = device_group,
  .SDcols = vars_table1
]

missing_by_group


####### SAVE ANALYSIS COHORT #######

saveRDS(dt_analysis, "study3_device_analysis_cohort.rds")
fwrite(missing_by_group, "supplement_missingness_by_device.csv")
fwrite(followup_summary, "followup_summary_by_device.csv")


######################### Table 2 #######################


dt_desc[, t_followup_days_final := NA_real_]

# Israel-specific follow-up (if present)
if ("t_followup_days_israel" %in% names(dt_desc)) {
  dt_desc[dataset == "ISRAEL",
          t_followup_days_final := t_followup_days_israel]
}

# All other datasets
dt_desc[dataset != "ISRAEL",
        t_followup_days_final := t_followup_days]

dt_analysis <- dt_desc


table(dt_analysis$Inapp_ATP)

########## DEFINE TABLE 2 DATASET (NO FOLLOW-UP FILTER) ##########

dt_table2 <- dt_analysis[
  dataset %in% c("EUCERT", "HELIOS", "PROSE", "ISRAEL")
]

########## CREATE new_SCD FROM Status (0=alive/censored, 1=SCD, 2=non-SCD) ##########
dt_table2[, new_Status3 := NA_character_]

if ("Status" %in% names(dt_analysis)) {
  dt_analysis[, new_Status3 := ifelse(
    Status == 0, "Alive/censored",
    ifelse(Status == 1, "Sudden cardiac death",
           ifelse(Status == 2, "Non-sudden cardiac death", NA_character_)
    )
  )]
}


########## NORMALISE ENDPOINT FLAGS (Yes / No ) ##########

endpoint_flags <- c(
  "inapp_shock_flag",
  "app_shock_flag",
  "death_flag"
)

endpoint_flags <- intersect(endpoint_flags, names(dt_table2))

dt_table2[, (endpoint_flags) := lapply(.SD, function(x) {
  x <- tolower(trimws(as.character(x)))
  ifelse(
    x %in% c("yes", "y", "1"), "Yes",
    ifelse(x %in% c("no", "n", "0"), "No", NA)
  )
}), .SDcols = endpoint_flags]

########## APPLY "MISSING INDICATOR = NO EVENT" FOR INAPP SHOCK (MAIN) ##########
if ("inapp_shock_flag" %in% names(dt_table2)) {
  dt_table2[, inapp_shock_flag_raw := inapp_shock_flag]
  dt_table2[is.na(inapp_shock_flag), inapp_shock_flag := "No"]
}

########## DEFINE TIME AND EVENT VARIABLES ##########

time_vars <- c(
  "t_followup_days_final",
  "days_to_inapp_shock",
  "days_to_app_shock",
  "days_to_death"
)
event_vars <- intersect(
  c(endpoint_flags, "new_Status3", "Innap_ATP_clean"),
  names(dt_table2)
)

vars_table2 <- c(time_vars, event_vars)
########## TABLE 2 — DESCRIPTIVES BY DEVICE GROUP ##########

table2 <- CreateTableOne(
  vars        = vars_table2,
  strata      = "device_group",
  data        = dt_table2,
  factorVars  = event_vars,
  includeNA  = TRUE,
  addOverall = TRUE
)

table2_print <- print(
  table2,
  nonnormal = time_vars,
  showAllLevels = TRUE,
  quote = FALSE,
  noSpaces = TRUE
)

table2_print


################INAP ATP ############
str(dt_analysis$Inapp_ATP)
setDT(dt_analysis)  # make sure it's a data.table

setDT(dt_analysis)

dt_analysis[, Innap_ATP_clean := {
  x <- tolower(trimws(as.character(.SD[[1]])))
  ifelse(
    x %in% c("yes","y","1"), "Yes",
    ifelse(x %in% c("no","n","0"), "No", NA_character_)
  )
}, .SDcols = "Inapp_ATP"]
x <- table(dt_analysis$Innap_ATP_clean)
prop.table(x)

y <- table(dt_analysis$Innap_ATP_clean, dt_analysis$device_group,useNA="ifany" )
y
prop.table(y,2)

chisq.test(dt_analysis$new_Status3, dt_analysis$device_group)


table(dt_analysis$new_Status3, dt_analysis$device_group)


########## FIX 3: MAKE “NA = NO EVENT” FOR INAP SHOCK ##########

# create binary versions JUST for incidence rates (0/1)
ir_endpoints <- c("inapp_shock_flag","app_shock_flag","death_flag")
ir_endpoints <- intersect(ir_endpoints, names(dt_table2))

for (v in ir_endpoints) {
  newv <- paste0(v, "_bin")
  dt_table2[, (newv) := fifelse(get(v) == "Yes", 1L, 0L)]  # NA -> 0
}

incidence_rates <- rbindlist(
  lapply(ir_endpoints, function(v) {
    vbin <- paste0(v, "_bin")
    dt_table2[
      ,
      {
        events <- sum(get(vbin), na.rm = TRUE)
        py     <- sum(t_followup_days_final, na.rm = TRUE) / 365.25
        rate   <- events / py * 100
        ci     <- poisson.test(events, py)$conf.int * 100
        
        .(
          events = events,
          person_years = round(py, 1),
          rate_per_100py = round(rate, 2),
          lci_95 = round(ci[1], 2),
          uci_95 = round(ci[2], 2)
        )
      },
      by = device_group
    ][, endpoint := v]
  })
)

########## INCIDENCE RATES PER 100 PERSON-YEARS ##########

incidence_rates <- rbindlist(
  lapply(
    c("inapp_shock_flag", "app_shock_flag", "death_flag"),
    function(v) {
      dt_table2[
        ,
        {
          events <- sum(get(v) == "Yes", na.rm = TRUE)
          py     <- sum(t_followup_days_final, na.rm = TRUE) / 365.25
          rate   <- events / py * 100
          ci     <- poisson.test(events, py)$conf.int * 100
          
          .(
            events = events,
            person_years = round(py, 1),
            rate_per_100py = round(rate, 2),
            lci_95 = round(ci[1], 2),
            uci_95 = round(ci[2], 2)
          )
        },
        by = device_group
      ][, endpoint := v]
    }
  )
)

fwrite(
  incidence_rates,
  "incidence_rates_per_100py_by_device_primary_all_datasets.csv"
)

########## QC 1: RAW ENDPOINT VALUES BY DATASET ##########

dt_table2[
  ,
  .N,
  by = .(dataset, raw_value = as.character(inapp_shock_flag))
][order(dataset, raw_value)]

########## QC 2: FOLLOW-UP SUMMARY BY DATASET ##########

dt_table2[
  ,
  .(
    n = .N,
    median_fu_days = median(t_followup_days_final, na.rm = TRUE),
    min_fu_days = min(t_followup_days_final, na.rm = TRUE)
  ),
  by = dataset
]

########## QC 3: SHOCK DATE PRESENT BUT FLAG MISSING ##########

dt_table2[
  is.na(inapp_shock_flag) & !is.na(inapp_shock_date),
  .N,
  by = dataset
]

########## QC 4: FOLLOW-UP BY INAPPROPRIATE SHOCK STATUS ##########

dt_table2[
  ,
  .(
    n = .N,
    median_fu_days = median(t_followup_days_final, na.rm = TRUE),
    iqr_fu_days = IQR(t_followup_days_final, na.rm = TRUE)
  ),
  by = .(dataset, inapp_shock_flag)
]

########## END OF TABLE 2 — PRIMARY ANALYSIS ##########


# Check what datasets exist BEFORE the table2 filter
dt_desc[, .N, by = dataset][order(dataset)]

# Check how many Israel rows you have + whether follow-up is missing
dt_desc[dataset == "ISRAEL",
        .(n = .N,
          n_fu_final_nonmissing = sum(!is.na(t_followup_days_final)),
          n_fu_israel_nonmissing = sum(!is.na(t_followup_days_israel)),
          n_fu_nonisrael_nonmissing = sum(!is.na(t_followup_days))
        )
]

# What datasets survive the dt_table2 filter?
dt_table2[, .N, by = dataset][order(dataset)]

##################### deaths after inappropriate shock #####################

# 1)  QC: mortality availability in dt_table2 
dt_table2[, .(
  n_total = .N,
  deaths_yes = sum(death_flag == "Yes", na.rm = TRUE),
  deaths_no  = sum(death_flag == "No",  na.rm = TRUE),
  deaths_na  = sum(is.na(death_flag))
)]

# 2) Ensure time-to-death is numeric 
if ("days_to_death" %in% names(dt_table2)) {
  dt_table2[, days_to_death_num := suppressWarnings(as.numeric(days_to_death))]
}

# 3) Descriptive: deaths among those with ≥1 inappropriate shock

inapp_subset <- dt_table2[inapp_shock_flag == "Yes"]

impact_inapp_on_death <- inapp_subset[, .(
  n_inapp = .N,
  deaths_inapp = sum(death_flag == "Yes", na.rm = TRUE),
  deaths_inapp_pct = round(100 * mean(death_flag == "Yes", na.rm = TRUE), 1),
  deaths_missing = sum(is.na(death_flag))
)]

impact_inapp_on_death

#check time-ordering (death after inapp shock) if i have days_to_inapp_shock

if (all(c("days_to_inapp_shock","days_to_death") %in% names(dt_table2))) {
  dt_table2[, days_to_inapp_num := suppressWarnings(as.numeric(days_to_inapp_shock))]
  dt_table2[, days_to_death_num := suppressWarnings(as.numeric(days_to_death))]
  
  death_after_inapp <- dt_table2[
    inapp_shock_flag == "Yes" & death_flag == "Yes" &
      !is.na(days_to_inapp_num) & !is.na(days_to_death_num),
    .(
      n_deaths_with_times = .N,
      n_death_after_inapp = sum(days_to_death_num >= days_to_inapp_num),
      n_death_before_inapp = sum(days_to_death_num <  days_to_inapp_num)
    )
  ]
  
  death_after_inapp
}

# 5) by-device breakdown 
impact_inapp_on_death_by_device <- dt_table2[
  inapp_shock_flag == "Yes",
  .(
    n_inapp = .N,
    deaths_inapp = sum(death_flag == "Yes", na.rm = TRUE)
  ),
  by = device_group
]

impact_inapp_on_death_by_device

############ Burden of inappropriate shocks per patient (recurrent)##############################################################

# Ensure totals are numeric
for (v in c("total_inapp_shock","total_app_shock","total_inapp_ATP","total_ATP")) {
  if (v %in% names(dt_table2)) dt_table2[, (v) := suppressWarnings(as.numeric(get(v)))]
}

# Burden summaries for inappropriate shocks
burden_inapp <- dt_table2[
  ,
  .(
    n = .N,
    n_any_inapp = sum(inapp_shock_flag == "Yes", na.rm = TRUE),
    prop_any_inapp = round(100 * mean(inapp_shock_flag == "Yes", na.rm = TRUE), 1),
    
    # Among those with inapp shock recorded as Yes:
    median_inapp_shocks = median(total_inapp_shock[inapp_shock_flag == "Yes"], na.rm = TRUE),
    iqr_inapp_shocks    = IQR(total_inapp_shock[inapp_shock_flag == "Yes"], na.rm = TRUE),
    prop_ge2_inapp      = round(100 * mean(total_inapp_shock[inapp_shock_flag == "Yes"] >= 2, na.rm = TRUE), 1),
    
    missing_total_inapp = sum(is.na(total_inapp_shock))
  ),
  by = device_group
]

burden_inapp

table(dt_desc$dataset, dt_desc$total_inapp_shock)


saveRDS(
  dt_analysis,
  "study3_analysis_cohort_v2_final_followup.rds"
)


