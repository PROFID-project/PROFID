
###########################SECONDARY ANALYSIS ####################################################

library(data.table)
library(survival)
#install.packages("mice")

dt_final <- as.data.table(readRDS("study3_analysis_final.rds"))

#  Quick inventory / missingness / coding 
# -----------------------------
vars_check <- c(
  "age_icd","Sex","BMI","eGFR","LVEF","Diabetes","AF_atrial_flutter","NYHA","Smoking",
  "Hypertension","Stroke_TIA","MR","device_group","dataset",
  "t_followup_days_final","event_inapp_shock","death_flag"
)
vars_check <- vars_check[vars_check %in% names(dt_final)]

print(round(colMeans(is.na(dt_final[, ..vars_check])), 3))

vars_code <- c("Sex","Diabetes","Hypertension","Stroke_TIA","AF_atrial_flutter",
               "NYHA","Smoking","MR","device_group")
vars_code <- vars_code[vars_code %in% names(dt_final)]

print(lapply(dt_final[, ..vars_code], function(x) table(x, useNA = "ifany")))

######### Factor refs + derived subgroup variables ##########
dt_final[, device_group := factor(device_group, levels = c("Dual","Single"))]

dt_final[, Sex := factor(Sex)]
if (all(c("Female","Male") %in% levels(dt_final$Sex))) {
  dt_final[, Sex := relevel(Sex, ref = "Female")]
}

# Age subgroup <65 vs >=65
dt_final[, age65 := fifelse(age_icd >= 65, 1L, 0L)]

# LVEF subgroup <30 vs >=30
dt_final[, lvef30 := fifelse(LVEF < 30, 1L, 0L)]

# NYHA collapse I–II vs III–IV (keep NA)
if ("NYHA" %in% names(dt_final)) {
  dt_final[, nyha_bin := fifelse(
    NYHA %in% c("I","II",1,2), "I-II",
    fifelse(NYHA %in% c("III","IV",3,4), "III-IV", NA_character_)
  )]
  dt_final[, nyha_bin := factor(nyha_bin, levels = c("I-II","III-IV"))]
}

# AF as factor No/Yes
if ("AF_atrial_flutter" %in% names(dt_final)) {
  dt_final[, AF_atrial_flutter := factor(AF_atrial_flutter, levels = c("No","Yes"))]
}

# Diabetes as factor No/Yes
if ("Diabetes" %in% names(dt_final)) {
  dt_final[, Diabetes := factor(Diabetes, levels = c("No","Yes"))]
}

# Stroke/TIA as factor No/Yes (for QC only; not used in modelling)
if ("Stroke_TIA" %in% names(dt_final)) {
  dt_final[, Stroke_TIA := factor(Stroke_TIA, levels = c("No","Yes"))]
}

######## Univariable Cox screening (stratified by dataset)  ###############

run_uv <- function(var) {
  d <- dt_final[complete.cases(dt_final[, .(t_followup_days_final, event_inapp_shock, dataset, get(var))])]
  
  if (nrow(d) == 0) return(NULL)
  if (is.factor(d[[var]]) && nlevels(d[[var]]) < 2) return(NULL)
  if (!is.factor(d[[var]]) && length(unique(d[[var]])) < 2) return(NULL)
  
  f <- as.formula(paste0("Surv(t_followup_days_final, event_inapp_shock) ~ ", var, " + strata(dataset)"))
  fit <- coxph(f, data = d)
  s <- summary(fit)
  
  data.table(
    var = var,
    term = rownames(s$coef),
    n = s$n,
    events = s$nevent,
    HR = s$coef[, "exp(coef)"],
    lci = s$conf.int[, "lower .95"],
    uci = s$conf.int[, "upper .95"],
    p = s$coef[, "Pr(>|z|)"],
    c_index = s$concordance[1]
  )
}

uv_vars <- c(
  "device_group",
  "age_icd",
  "Sex",
  "LVEF",
  "Diabetes",
  "BMI",
  "eGFR",
  "AF_atrial_flutter",
  "nyha_bin"
)
uv_vars <- uv_vars[uv_vars %in% names(dt_final)]

uv_res <- rbindlist(lapply(uv_vars, run_uv), fill = TRUE)
uv_res <- uv_res[order(p)]

cat("\n--- Univariable Cox (sorted by p) ---\n")
print(uv_res)

cat("\n--- Univariable Cox (p < .10) ---\n")
print(uv_res[p < 0.10])

############## Multivariable Cox complete-case: 
## device_group + age_icd + Sex + LVEF + Diabetes; stratified by dataset ##########
mv_vars <- c("device_group","age_icd","Sex","LVEF","Diabetes")
mv_vars <- mv_vars[mv_vars %in% names(dt_final)]

d_mv <- dt_final[complete.cases(dt_final[, c("t_followup_days_final","event_inapp_shock","dataset", mv_vars), with = FALSE])]

f_mv <- as.formula(
  paste0("Surv(t_followup_days_final, event_inapp_shock) ~ ",
         paste(mv_vars, collapse = " + "),
         " + strata(dataset)")
)

fit_mv <- coxph(f_mv, data = d_mv)
s_mv <- summary(fit_mv)

cat("\n--- Multivariable Cox (complete-case) ---\n")
print(s_mv)

cat("\n--- PH check ---\n")
ph_mv <- cox.zph(fit_mv)
print(ph_mv)

cat("\nHarrell's C (MV Cox): ", round(s_mv$concordance[1], 3),
    " SE: ", round(s_mv$concordance[2], 3), "\n", sep = "")

############### Subgroup interaction tests ###############
run_interaction <- function(subvar) {
  if (!(subvar %in% names(dt_final))) return(NULL)
  
  d <- dt_final[complete.cases(dt_final[, .(t_followup_days_final, event_inapp_shock, dataset, device_group, get(subvar))])]
  if (nrow(d) == 0) return(NULL)
  
  if (!is.factor(d[[subvar]]) && !is.numeric(d[[subvar]]) && !is.integer(d[[subvar]])) return(NULL)
  
  f <- as.formula(paste0("Surv(t_followup_days_final, event_inapp_shock) ~ device_group * ", subvar, " + strata(dataset)"))
  fit <- coxph(f, data = d)
  s <- summary(fit)
  
  int_row <- grep("device_groupSingle:", rownames(s$coef), value = TRUE)
  if (length(int_row) != 1) return(NULL)
  
  data.table(
    subgroup = subvar,
    n = s$n,
    events = s$nevent,
    interaction_term = int_row,
    HR_int = s$coef[int_row, "exp(coef)"],
    lci = s$conf.int[int_row, "lower .95"],
    uci = s$conf.int[int_row, "upper .95"],
    p_int = s$coef[int_row, "Pr(>|z|)"]
  )
}

int_vars <- c("age65","lvef30","AF_atrial_flutter")
int_vars <- int_vars[int_vars %in% names(dt_final)]

int_res <- rbindlist(lapply(int_vars, run_interaction), fill = TRUE)

cat("\n--- Interaction tests ---\n")
print(int_res)

######### Stroke_TIA QC: separation check (not modelled) ######
if ("Stroke_TIA" %in% names(dt_final)) {
  cat("\n--- Stroke_TIA QC: events by level (for separation check) ---\n")
  print(
    dt_final[!is.na(Stroke_TIA),
             .(n = .N, events = sum(event_inapp_shock), event_rate = mean(event_inapp_shock)),
             by = Stroke_TIA]
  )
}

