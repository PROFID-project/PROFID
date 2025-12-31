############################################################
# MICE imputation per cohort (Phase 1 rules)
# - Impute only variables with ≤80% missingness
# - Don't impute outcomes; don't impute comorbidities/meds (missing=0)
# - Use PMM for numeric, logreg for binary, polyreg for multiclass
############################################################

# Packages 
req <- c("data.table","mice","dplyr","stringr")
inst <- setdiff(req, rownames(installed.packages()))
if (length(inst)) install.packages(inst, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Input files 
paths <- list(
  ICD              = "ICD_filtered_with_coords.csv",
  NonICD_reduced   = "NonICD_reduced_filtered_with_coords.csv",
  NonICD_preserved = "NonICD_preserved_filtered_with_coords.csv"
)
stopifnot(all(file.exists(unlist(paths))))

# Outcomes (NEVER impute) 
outcome_vars <- c("Status","Survival_time")

# Technical/ID columns (exclude from model) 
id_vars <- c("X","DB","ID","Time_zero_Y","Time_zero_Ym","Time_zero_y","Time_zero_Ym_") # safe extras; keep those that exist

# “Missing = absence” binary indicators (DON’T impute) 
# based on your screenshots (ICD + NonICD columns)
bin_no_impute <- c(
  # comorbidities / history
  "Diabetes","Hypertension","Smoking","COPD","FH_SCD","FH_CAD","Stroke_TIA",
  "AF_atrial_flutter","NSVT","RBBB","LBBB","AV_block","AV_block_II_or_III","Alcohol","Cancer",
  # procedures / acute care (present in NonICD set too)
  "PCI","CABG","revascularisation_acute","thrombolysis_acute","MI_history",
  # meds
  "ACE_inhibitor_ARB","ACE_inhibitor","ARB","Beta_blockers","Diuretics",
  "Anti_platelet","Anticoagulant","Anti_anginal","Calcium_antagonists",
  "Digitalis_glycosides","Anti_diabetic_insulin","Anti_diabetic_oral","Lipid_lowering",
  # MRI / echo presence flags
  "HasMRI"
)

# Helper: cap to existing columns
cap_to_existing <- function(df, vars) intersect(vars, names(df))

# Read data 
read_dt <- function(p) as.data.frame(data.table::fread(p, data.table = FALSE))
dfs <- lapply(paths, read_dt)

# Main function to impute one cohort 
impute_cohort <- function(df, cohort_name, m = 5, maxit = 10, seed = 123) {

  # 1) Set missing -> 0 for “no-impute” binaries that exist
  no_imp <- cap_to_existing(df, bin_no_impute)
  for (v in no_imp) {
    # coerce to numeric 0/1 if not already; treat "Yes"/"No" etc.
    if (is.character(df[[v]]) || is.factor(df[[v]])) {
      df[[v]] <- ifelse(tolower(as.character(df[[v]])) %in% c("1","yes","y","true","t"), 1,
                        ifelse(tolower(as.character(df[[v]])) %in% c("0","no","n","false","f"), 0, NA))
    }
    df[[v]][is.na(df[[v]])] <- 0
  }

  # 2) Compute missingness and select variables with ≤80% missing
  miss_pct <- sapply(df, function(x) mean(is.na(x)))
  vars_le80 <- names(miss_pct)[miss_pct <= 0.80]

  # 3) Build the set to be **imputed** (remove outcomes, no-impute binaries, IDs)
  do_not_model <- union(outcome_vars, union(bin_no_impute, id_vars))
  vars_for_mice <- setdiff(vars_le80, cap_to_existing(df, do_not_model))

  # 4) Cast obvious binary (0/1) to factor for logreg/polyreg decisions later
  is_binary01 <- function(x) {
    if (is.numeric(x)) {
      u <- unique(x[!is.na(x)])
      all(u %in% c(0,1)) && length(u) <= 2
    } else FALSE
  }
  for (nm in vars_for_mice) {
    if (is_binary01(df[[nm]])) df[[nm]] <- factor(df[[nm]], levels = c(0,1))
  }

  # 5) Prepare predictor matrix & methods
  #   - variables to impute: those in vars_for_mice that still have missingness
  imp_vars <- vars_for_mice[sapply(df[vars_for_mice], function(x) any(is.na(x)))]
  if (length(imp_vars) == 0) {
    message("[", cohort_name, "] No variables require imputation under the rules. Returning original df.")
    return(list(
      completed = df,
      log = list(
        cohort = cohort_name,
        dropped_gt80 = names(miss_pct)[miss_pct > 0.80],
        not_imputed_missing_eq_no = cap_to_existing(df, bin_no_impute),
        outcomes_excluded = cap_to_existing(df, outcome_vars),
        ids_excluded = cap_to_existing(df, id_vars),
        imputed_vars = character(0)
      )
    ))
  }

  meth <- make.method(df)
  pred <- make.predictorMatrix(df)

  # Exclude IDs & outcomes & no-impute vars from being imputed and from predictors
  excl <- cap_to_existing(df, union(outcome_vars, union(bin_no_impute, id_vars)))
  meth[excl] <- ""                       # never impute these
  pred[, excl] <- 0                      # don't use them as predictors either (conservative)

  # Set method by variable type
  for (v in imp_vars) {
    x <- df[[v]]
    if (is.numeric(x)) {
      meth[v] <- "pmm"
    } else if (is.factor(x)) {
      if (nlevels(x) == 2) meth[v] <- "logreg" else meth[v] <- "polyreg"
    } else if (is.character(x)) {
      # coerce to factor and use polyreg
      df[[v]] <- factor(x)
      meth[v] <- if (nlevels(df[[v]]) == 2) "logreg" else "polyreg"
    } else {
      # fallback
      meth[v] <- "pmm"
    }
  }

  # Variables with no missing get empty method so mice skips imputing them
  no_miss <- setdiff(vars_for_mice, imp_vars)
  meth[no_miss] <- ""

  # Keep predictors only among vars_for_mice (plus a few strong predictors if you wish)
  keep_cols <- union(vars_for_mice, setdiff(names(df), vars_for_mice))
  pred[ , setdiff(colnames(pred), vars_for_mice)] <- 0

  # 6) Run MICE
  set.seed(seed)
  imp <- mice(df, m = m, maxit = maxit, method = meth, predictorMatrix = pred, printFlag = TRUE)

  # 7) Take first completed dataset (you can pool later in models if desired)
  completed <- complete(imp, 1)

  # 8) Return completed data + a log
  log_list <- list(
    cohort = cohort_name,
    dropped_gt80 = names(miss_pct)[miss_pct > 0.80],
    not_imputed_missing_eq_no = cap_to_existing(df, bin_no_impute),
    outcomes_excluded = cap_to_existing(df, outcome_vars),
    ids_excluded = cap_to_existing(df, id_vars),
    imputed_vars = imp_vars
  )

  list(completed = completed, imp = imp, log = log_list)
}

# Run for the three cohorts 
res_ICD <- impute_cohort(dfs$ICD, "ICD")
res_NR  <- impute_cohort(dfs$NonICD_reduced, "NonICD_reduced")
res_NP  <- impute_cohort(dfs$NonICD_preserved, "NonICD_preserved")

# Save completed datasets 
data.table::fwrite(res_ICD$completed, "ICD_imputed.csv")
data.table::fwrite(res_NR$completed,  "NonICD_reduced_imputed.csv")
data.table::fwrite(res_NP$completed,  "NonICD_preserved_imputed.csv")

# Save quick logs (what we imputed / dropped) 
log_to_dt <- function(L) {
  data.table::data.table(
    cohort = L$cohort,
    dropped_gt80 = paste(L$dropped_gt80, collapse = "; "),
    not_imputed_missing_eq_no = paste(L$not_imputed_missing_eq_no, collapse = "; "),
    outcomes_excluded = paste(L$outcomes_excluded, collapse = "; "),
    ids_excluded = paste(L$ids_excluded, collapse = "; "),
    imputed_vars = paste(L$imputed_vars, collapse = "; ")
  )
}
logs <- rbindlist(list(log_to_dt(res_ICD$log),
                       log_to_dt(res_NR$log),
                       log_to_dt(res_NP$log)),
                  fill = TRUE)
data.table::fwrite(logs, "MICE_log_summary.csv")

message("Done. Wrote: ICD_imputed.csv, NonICD_reduced_imputed.csv, NonICD_preserved_imputed.csv, and MICE_log_summary.csv")
