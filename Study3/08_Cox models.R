# STUDY 3 — SURVIVAL QC & COX MODELS (INAPPROPRIATE SHOCK)

# install.packages("survival")

library(data.table)
library(survival)


#################### 1. LOAD DATA ####################
dt <- readRDS("study3_analysis_cohort_v2_final_followup.rds")

#################### 2. HARMONISE ICD TYPE → DEVICE GROUP ####################
dt[, .N, by = .(dataset, device_group)][order(dataset, device_group)]

# Single-chamber ICDs
dt[
  tolower(icd_type) %in% c("vvi", "icd_1", "single - single", "single chamber"),
  device_group := "Single"
]

# Dual-chamber ICDs
dt[
  tolower(icd_type) %in% c("ddd", "icd_2", "dual - dual", "dual chamber"),
  device_group := "Dual"
]

# QC
dt[, .N, by = .(dataset, device_group)]
dt[
  is.na(device_group),
  .N,
  by = .(dataset, icd_type)
]

#################### 3. HARMONISE INAPPROPRIATE SHOCK FLAG ####################
# Raw values:
# EUCERT: yes / no / NA
# HELIOS, PROSE: Yes / No
# ISRAEL: YES / NO / UNKNOWN

dt[, inapp_shock_flag_std :=
     tolower(trimws(as.character(inapp_shock_flag)))]

#################### 4. DEFINE EVENT INDICATOR ####################
# YES = 1
# NO / UNKNOWN / NA = 0

dt[, event_inapp_shock :=
     fifelse(inapp_shock_flag_std == "yes", 1L, 0L, na = 0L)
]

# QC
dt[, .N, by = .(dataset, event_inapp_shock)]


dt[
  event_inapp_shock == 1 &
    !is.na(days_to_inapp_shock) &
    days_to_inapp_shock > t_followup_days_final,
  days_to_inapp_shock := t_followup_days_final
]

########## REMOVE INVALID FOLLOW-UP ##########
# (No time at risk)

dt <- dt[t_followup_days_final > 0]


########## PRIMARY ANALYSIS COHORT (WITH ISRAEL) ##########

dt_primary <- dt[
  dataset %in% c("EUCERT", "HELIOS", "PROSE", "ISRAEL") &
    !is.na(device_group) &
    !is.na(t_followup_days_final)
]

prop.table(table(dt_primary$device_group, dt_primary$event_inapp_shock),1)
table(dt_primary$device_group, dt_primary$event_inapp_shock)
prop.table(table(dt_primary$device_group, dt_primary$inapp_shock_flag),1)

dt_primary[, .(
  n = .N,
  events = sum(event_inapp_shock)
), by = .(dataset, device_group)][order(dataset, device_group)]

########## PRIMARY COX MODEL ##########
# Stratified by dataset

cox_primary <- coxph(
  Surv(t_followup_days_final, event_inapp_shock) ~
    device_group + strata(dataset),
  data = dt_primary
)

summary(cox_primary)

########## SENSITIVITY ANALYSIS (EXCLUDING ISRAEL) ##########

dt_sens_no_israel <- dt[
  dataset %in% c("EUCERT", "HELIOS", "PROSE") &
    !is.na(device_group) &
    !is.na(t_followup_days_final)
]

cox_sens_no_israel <- coxph(
  Surv(t_followup_days_final, event_inapp_shock) ~
    device_group + strata(dataset),
  data = dt_sens_no_israel
)

summary(cox_sens_no_israel)


################# SENSITIVITY ANALYSIS (EXCLUDING EUCERT (TOO MANY NA'S IN INAPP SHOCK VARIABLE) ######
dt_inapp_valid <- dt[
  dataset %in% c("HELIOS","ISRAEL","PROSE") &
    !is.na(device_group) &
    !is.na(t_followup_days_final)
]

dt_inapp_valid[, .(n=.N, events=sum(event_inapp_shock)), by=.(dataset, device_group)]

cox_valid <- coxph(
  Surv(t_followup_days_final, event_inapp_shock) ~
    device_group +  strata(dataset),
  data = dt_inapp_valid
)

summary(cox_valid)


############ SENSITIVITY ANALYSIS (EXCLUDING ALL NA INAPP SHOCKS)######

dt_noNA <- dt[!is.na(inapp_shock_flag)]

cox_noNA <- coxph(
  Surv(t_followup_days_final, event_inapp_shock) ~
    device_group + strata(dataset),
  data = dt_noNA
)
summary(cox_noNA)


########## ADJUSTED MODEL (SUPPLEMENTARY) ##########

cox_adj <- coxph(
  Surv(t_followup_days_final, event_inapp_shock) ~
    device_group + age_icd + Sex + strata(dataset),
  data = dt_primary
)

summary(cox_adj)


########## QC 0: REQUIRED VARIABLES PRESENT ##########

req <- c("dataset","device_group","t_followup_days_final","event_inapp_shock")
missing_req <- setdiff(req, names(dt_primary))
stopifnot(length(missing_req) == 0)

########## QC 1: NO NEGATIVE / ZERO FOLLOW-UP ##########

dt_primary[t_followup_days_final <= 0, .N]   # should be 0
stopifnot(dt_primary[t_followup_days_final <= 0, .N] == 0)

########## QC 2: EVENT CODING IS BINARY ##########

dt_primary[, .N, by = event_inapp_shock][order(event_inapp_shock)]
stopifnot(all(na.omit(unique(dt_primary$event_inapp_shock)) %in% c(0L,1L)))

########## QC 3: EVENT COUNTS BY DATASET × DEVICE ##########

dt_primary[
  ,
  .(n = .N, events = sum(event_inapp_shock, na.rm = TRUE)),
  by = .(dataset, device_group)
][order(dataset, device_group)]


########## QC 5: EVENT TIME <= FOLLOW-UP (ONLY IF days_to_* EXISTS) ##########

dt_primary[
  event_inapp_shock == 1 &
    !is.na(days_to_inapp_shock) &
    days_to_inapp_shock > t_followup_days_final,
  .N
]  

########## QC 6: ENOUGH EVENTS PER STRATUM ##########

dt_primary[
  ,
  .(events = sum(event_inapp_shock)),
  by = dataset
][order(events)]

########## QC 7: PROPORTIONAL HAZARDS ASSUMPTION ##########


fit <- coxph(
  Surv(t_followup_days_final, event_inapp_shock) ~ device_group + strata(dataset),
  data = dt_primary
)

ph <- cox.zph(fit)
print(ph)     
 plot(ph)       


# --- CHECK 1: registry-level Cox (within each dataset) ---
by_dataset_cox <- dt_primary[
  ,
  {
    fit <- coxph(Surv(t_followup_days_final, event_inapp_shock) ~ device_group, data = .SD)
    s <- summary(fit)
    list(
      HR = s$conf.int["device_groupSingle","exp(coef)"],
      lci = s$conf.int["device_groupSingle","lower .95"],
      uci = s$conf.int["device_groupSingle","upper .95"],
      p   = s$coefficients["device_groupSingle","Pr(>|z|)"],
      N = .N,
      events = sum(event_inapp_shock)
    )
  },
  by = dataset
]

by_dataset_cox

# --- CHECK 2: leave-one-dataset-out sensitivity ---

datasets <- unique(dt_primary$dataset)

loo_results <- lapply(datasets, function(d){
  tmp <- dt_primary[dataset != d]
  
  fit <- coxph(
    Surv(t_followup_days_final, event_inapp_shock) ~
      device_group + strata(dataset),
    data = tmp
  )
  
  s <- summary(fit)
  
  data.table(
    excluded = d,
    HR = s$conf.int["device_groupSingle","exp(coef)"],
    lci = s$conf.int["device_groupSingle","lower .95"],
    uci = s$conf.int["device_groupSingle","upper .95"],
    p   = s$coefficients["device_groupSingle","Pr(>|z|)"],
    N = nrow(tmp),
    events = sum(tmp$event_inapp_shock)
  )
})

loo_table <- rbindlist(loo_results)
loo_table


