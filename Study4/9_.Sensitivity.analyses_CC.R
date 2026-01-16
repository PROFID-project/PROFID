##################################################################################################################
##################################################################################################################
# projet: UmBIZO- PROFID_Study4 (SCD post-MI, age- stratified)
# Script: 09_Sensitivity_analyses.R
# Author: Amina Boudamaana
# =========================================================================================
## STUDY 4 — Sensitivity analyses (COX) on complete cases (df_handled_cc)
## Scenarios:
## 1) Alternative age cut-offs
## 2) Exclude patients missing key variables
## 3) Restrict to patients with complete echocardiography
## 4) Separate analyses: ICD vs. non-ICD
## =========================================================================================
suppressPackageStartupMessages({ library(survival) })
# ---------- 0) Input dataset ----------
stopifnot(exists("df_handled_cc"), is.data.frame(df_handled_cc))
df <- df_handled_cc

# Harmonize survival variables (years + binary SCD event for Cox)
if (!"time_years" %in% names(df)) {
  cand_t <- c("Survival_time","survival_time","ftime_mo_int","ftime","time_months")
  nm_t <- cand_t[cand_t %in% names(df)][1]
  stopifnot(!is.na(nm_t))
  years <- if (grepl("year", nm_t, ignore.case = TRUE)) as.numeric(df[[nm_t]])
  else as.numeric(df[[nm_t]])/12
  df$time_years <- years
}
if (!"event_scd" %in% names(df)) {
  cand_s <- c("Status","status","fstatus","status_cr")
  nm_s <- cand_s[cand_s %in% names(df)][1]
  stopifnot(!is.na(nm_s))
  st <- as.integer(df[[nm_s]])
  df$event_scd <- as.integer(st == 1L)
}
# ===============================================================================================================
# ---------- 1) Output directory & helpers ----------
# ===============================================================================================================
ensure_writable <- function(path){
  ok <- tryCatch({
    if (!dir.exists(path)) dir.create(path, recursive = TRUE, showWarnings = FALSE)
    tf <- tempfile(tmpdir = path); writeLines("ok", tf); unlink(tf); TRUE
  }, error = function(e) FALSE)
  ok
}
TAB_DIR <- "T:/study_4/Results_tables"
if (!ensure_writable(TAB_DIR)) TAB_DIR <- file.path(getwd(), "Results_tables")
invisible(ensure_writable(TAB_DIR))
message("[Export] tables -> ", normalizePath(TAB_DIR, winslash = "/", mustWork = FALSE))

safe_write_csv <- function(x, file){
  if (is.null(x) || !NROW(x)) return(invisible(FALSE))
  tryCatch(utils::write.csv(x, file, row.names = FALSE),
           error = function(e) message("Could not write ", basename(file), ": ", e$message))
}

pick_present <- function(x, nm = names(df)) intersect(x, nm)

# Recode logical/character to No/Yes factor
recode_yesno <- function(z){
  if (is.factor(z)) z <- as.character(z)
  yes <- c("1","y","Y","yes","Yes","YES","true","True","TRUE")
  no <- c("0","n","N","no","No","NO","false","False","FALSE")
  out <- ifelse(z %in% yes, "Yes", ifelse(z %in% no, "No", NA))
  factor(out, levels = c("No","Yes"))
}
harmonise_binary <- function(d, cols){
  for (v in cols) if (v %in% names(d)){
    z <- d[[v]]
    if (is.logical(z) || is.character(z) || is.factor(z)) d[[v]] <- recode_yesno(z)
  }
  d
}
# Canonicalize medication synonyms to prevent duplicates in models/tables
canon_med_cols <- function(d, prefer = "Anti_coagulant"){
  syn <- intersect(c("Anti_coagulant","Ant_coag","anticoag","anticoagulant",
                     "Anticoagulant_any","anticoagulant_any"), names(d))
  if (length(syn)) {
    zlist <- lapply(d[syn], function(x) { # tout en texte
      x <- as.character(x)
      x <- trimws(tolower(x))
      ifelse(x %in% c("1","y","yes","true"), "Yes",
             ifelse(x %in% c("0","n","no","false"), "No", NA))
    })
    z <- Reduce(function(a,b) ifelse(is.na(a), b, a), zlist) # priorité 1er non-NA
    d[[prefer]] <- factor(z, levels = c("No","Yes"))
    # on SUPPRIME seulement les colonnes synonymes autres que 'prefer'
    to_drop <- setdiff(syn, prefer)
    if (length(to_drop)) d[to_drop] <- NULL
  }
  d
}
df_handled_cc <- canon_med_cols(df_handled_cc, prefer = "Anti_coagulant")
cc_mask <- function(d, vars){
  need <- c("time_years","event_scd", vars)
  complete.cases(d[, need, drop = FALSE]) &
    is.finite(d$time_years) & d$time_years > 0 & d$event_scd %in% c(0,1)
}
hr_table_from_fit <- function(fit){
  s <- summary(fit); co <- s$coefficients
  data.frame(
    term = rownames(co),
    HR = exp(co[,"coef"]),
    LCL_95 = exp(co[,"coef"] - 1.96 * co[,"se(coef)"]),
    UCL_95 = exp(co[,"coef"] + 1.96 * co[,"se(coef)"]),
    p_value= co[,"Pr(>|z|)"],
    row.names = NULL, check.names = FALSE
  )
}
# ====================================================================================================
# ---------- 2) Age helpers ----------
# ====================================================================================================
norm_age <- function(x){
  x <- trimws(as.character(x))
  x <- gsub("\u2264","<=", x); x <- gsub("\u2013|\u2014|–|—","-", x); x <- gsub("\\s+","", x)
  x[x %in% c("<=50","<=50y","<=50yrs","<=50years")] <- "<=50"
  x[x %in% c("51-65","51–65")] <- "51-65"
  x[x %in% c("66-75","66–75")] <- "66-75"
  x[grepl("^> ?75|^>=?76|\\b76\\+\\b", x)] <- ">75"
  factor(x, levels = c("<=50","51-65","66-75",">75"), ordered = TRUE)
}
get_age_years <- function(d){
  for (nm in c("Age","age","age_years","Age_years")) {
    if (nm %in% names(d) && is.numeric(d[[nm]])) return(as.numeric(d[[nm]]))
  }
  rep(NA_real_, nrow(d))
}
make_ageband_from_cutoffs <- function(age_years, cuts){
  stopifnot(length(cuts) == 3L)
  labs <- c(paste0("<=", cuts[1]),
            paste0(cuts[1]+1,"-", cuts[2]),
            paste0(cuts[2]+1,"-", cuts[3]),
            paste0(">", cuts[3]))
  cutp <- c(-Inf, cuts, Inf)
  out <- cut(age_years, breaks = cutp, labels = labs,
             right = TRUE, include.lowest = TRUE, ordered_result = TRUE)
  factor(as.character(out), levels = unique(as.character(out)), ordered = TRUE)
}
# ===================================================================================================
# ---------- 3) Variables used by models ----------
# ===================================================================================================
CONTINUOUS_CANDIDATES <- c("LVEF","eGFR","BMI","Haemoglobin","Cholesterol")
get_final_vars <- function(band){
  if (exists("results") && is.list(results) && !is.null(results[[band]]$final_vars)) {
    v <- results[[band]]$final_vars
    return(v[v %in% names(df)])
  }
  fs_csv <- file.path(TAB_DIR, "final_summary.csv")
  if (file.exists(fs_csv)) {
    fs <- try(read.csv(fs_csv, stringsAsFactors = FALSE), silent = TRUE)
    if (!inherits(fs,"try-error") && all(c("age_band","final_vars") %in% names(fs))) {
      row <- fs[match(band, fs$age_band), , drop = FALSE]
      if (NROW(row)) {
        v <- trimws(unlist(strsplit(row$final_vars, ",")))
        return(v[nzchar(v) & v %in% names(df)])
      }
    }
  }
  unique(pick_present(c(
    "LVEF","eGFR","BMI","Haemoglobin",
    "Diabetes","Hypertension","ACE_inhibitor","ARB","Anti_coagulant"
  )))
}
# =============================================================================================================
# ---------- 4) Scenario-specific dataset builders ----------
# =============================================================================================================
build_cc_keyvars <- function(d){
  key <- pick_present(c("LVEF","eGFR","Diabetes"))
                        
  d <- canon_med_cols(d)
  d <- harmonise_binary(d, key)
  keep <- if (length(key)) complete.cases(d[, key, drop = FALSE]) else TRUE
  d[keep, , drop = FALSE]
}
build_cc_echo <- function(d){
  echo_candidates <- pick_present(c("LVEF","LVEDD","LAVI","LVESD","MR_severity","Mitral_regurgitation"))
  if (!length(echo_candidates)) echo_candidates <- "LVEF"
  keep <- complete.cases(d[, echo_candidates, drop = FALSE])
  d[keep, , drop = FALSE]
}
detect_icd_col <- function(d){
  cands <- pick_present(c("ICD","ICD_implanted","ICD_at_discharge","ICD_status","ICD_present"), names(d))
  if (!length(cands)) return(NA_character_)
  cands[1]
}
subset_icd <- function(d, want_icd = TRUE){
  icd_col <- detect_icd_col(d)
  if (is.na(icd_col)) return(NULL)
  x <- d[[icd_col]]
  if (is.character(x) || is.logical(x) || is.factor(x)) x <- recode_yesno(x)
  keep <- if (is.factor(x)) x == "Yes" else as.numeric(x) == 1
  if (!want_icd) keep <- !keep
  d[keep & !is.na(keep), , drop = FALSE]
}
## base = df_handled_cc -> on ajoute time_years / event_scd
stopifnot(exists("df_handled_cc"), is.data.frame(df_handled_cc))
df <- within(df_handled_cc, {
  time_years <- as.numeric(ftime_mo_int)/12
  event_scd <- as.integer(fstatus == 1L) # 1 = SCD, 0/2 = autre/censuré
})
# ==============================================================================================
# ---------- 5) Analysis engine ----------
# ==============================================================================================
analyze_dataset <- function(d, scenario_label, alt_age = NULL){
  stopifnot(is.data.frame(d),
            all(c("time_years","event_scd") %in% names(d)))
  ## build age_group & normalize ---
  if (!is.null(alt_age)) {
    agey <- get_age_years(d)
    if (!any(is.finite(agey)))
      stop("Alternative cut-offs requested but numeric 'Age' not found.")
    d$age_group <- make_ageband_from_cutoffs(agey, alt_age)
  } else {
    if ("age_group" %in% names(d)) {
      d$age_group <- norm_age(d$age_group)
      if (all(is.na(d$age_group))) {
        agey <- get_age_years(d)
        d$age_group <- make_ageband_from_cutoffs(agey, c(50,65,75))
      }
    } else {
      agey <- get_age_years(d)
      d$age_group <- make_ageband_from_cutoffs(agey, c(50,65,75))
    }
  }
  
  ## Keep rows where age_group is non-NA
  d <- d[!is.na(d$age_group), , drop = FALSE]
  bands <- levels(d$age_group)
  
  out_HR <- list()
  out_SUM <- list()
  
  for (b in bands){
    dd <- d[d$age_group == b, , drop = FALSE]
    
    ## Final variables for this strip/section
    vars <- get_final_vars(b)
    
    ## Harmonize binary and filter in CC on these variables
    dd <- harmonise_binary(dd, vars)
    ok <- cc_mask(dd, vars)
    dd <- dd[ok, , drop = FALSE]
    
    ev <- sum(dd$event_scd == 1L, na.rm = TRUE)
    if (NROW(dd) < 25 || ev < 5 || length(vars) == 0){
      out_SUM[[b]] <- data.frame(
        scenario = scenario_label, age_band = b,
        N_CC = NROW(dd), events = ev,
        EPV = ifelse(length(vars)>0, round(ev/length(vars),1), NA_real_),
        n_vars = length(vars),
        final_vars = paste(vars, collapse = ", "),
        note = "insufficient data"
      )
      next
    }
    
    fm <- as.formula(paste("Surv(time_years, event_scd) ~", paste(vars, collapse = " + ")))
    fit <- try(coxph(fm, data = dd, ties = "efron", na.action = na.omit), silent = TRUE)
    if (inherits(fit, "try-error")){
      out_SUM[[b]] <- data.frame(
        scenario = scenario_label, age_band = b,
        N_CC = NROW(dd), events = ev,
        EPV = ifelse(length(vars)>0, round(ev/length(vars),1), NA_real_),
        n_vars = length(vars),
        final_vars = paste(vars, collapse = ", "),
        note = "fit error"
      )
      next
    }
    hr <- hr_table_from_fit(fit)
    hr$scenario <- scenario_label
    hr$age_band <- b
    out_HR[[b]] <- hr
    
    out_SUM[[b]] <- data.frame(
      scenario = scenario_label, age_band = b,
      N_CC = NROW(dd), events = ev,
      EPV = round(ev / max(1, length(vars)), 1),
      n_vars = length(vars),
      final_vars = paste(vars, collapse = ", "),
      note = "ok"
    )
  }
  
  list(
    hr = if (length(out_HR)) do.call(rbind, out_HR) else NULL,
    summary = if (length(out_SUM)) do.call(rbind, out_SUM) else NULL
  )
}
# ============================================================================================
    # ---------- 6) Run all scenarios & export ----------
# ============================================================================================
    ALL_HR <- list(); ALL_SUM <- list()
    
    # (1) Alternative age cut-offs
    AGE_SCHEMES <- list(
      ALT_55_64_74 = c(55, 64, 74), # <=55, 56–64, 65–74, >74
      ALT_60_70_75 = c(60, 70, 75) # <=60, 61–70, 71–75, >75
    )
    if (length(AGE_SCHEMES)){
      for (nm in names(AGE_SCHEMES)){
        res <- analyze_dataset(df, scenario_label = paste0("AGE_ALT_", nm), alt_age = AGE_SCHEMES[[nm]])
        if (!is.null(res$hr)) { fn <- file.path(TAB_DIR, paste0("sens_COX_", "AGE_ALT_", nm, ".csv")); safe_write_csv(res$hr, fn); ALL_HR[[paste0("AGE_ALT_", nm)]] <- res$hr }
        if (!is.null(res$summary)) { fn <- file.path(TAB_DIR, paste0("sens_COX_SUM_", "AGE_ALT_", nm, ".csv")); safe_write_csv(res$summary, fn); ALL_SUM[[paste0("AGE_ALT_", nm)]] <- res$summary }
      }
    }
    # (2) Exclude patients with missing key variables
    df_cc <- build_cc_keyvars(df)
    res_cc <- analyze_dataset(df_cc, scenario_label = "CC_KEYVARS")
    if (!is.null(res_cc$hr)) { safe_write_csv(res_cc$hr, file.path(TAB_DIR, "sens_COX_CC_KEYVARS.csv")); ALL_HR[["CC_KEYVARS"]] <- res_cc$hr }
    if (!is.null(res_cc$summary)) { safe_write_csv(res_cc$summary, file.path(TAB_DIR, "sens_COX_SUM_CC_KEYVARS.csv")); ALL_SUM[["CC_KEYVARS"]] <- res_cc$summary }
    
    # (3) Restrict to complete echocardiography
    df_echo <- build_cc_echo(df)
    res_echo <- analyze_dataset(df_echo, scenario_label = "COMPLETE_ECHO")
    if (!is.null(res_echo$hr)) { safe_write_csv(res_echo$hr, file.path(TAB_DIR, "sens_COX_COMPLETE_ECHO.csv")); ALL_HR[["COMPLETE_ECHO"]] <- res_echo$hr }
    if (!is.null(res_echo$summary)) { safe_write_csv(res_echo$summary, file.path(TAB_DIR, "sens_COX_SUM_COMPLETE_ECHO.csv")); ALL_SUM[["COMPLETE_ECHO"]] <- res_echo$summary }
    
    # (4) ICD vs non-ICD
    df_icd <- subset_icd(df, want_icd = TRUE)
    if (!is.null(df_icd)){
      res_icd <- analyze_dataset(df_icd, scenario_label = "ICD_ONLY")
      if (!is.null(res_icd$hr)) { safe_write_csv(res_icd$hr, file.path(TAB_DIR, "sens_COX_ICD_ONLY.csv")); ALL_HR[["ICD_ONLY"]] <- res_icd$hr }
      if (!is.null(res_icd$summary)) { safe_write_csv(res_icd$summary, file.path(TAB_DIR, "sens_COX_SUM_ICD_ONLY.csv")); ALL_SUM[["ICD_ONLY"]] <- res_icd$summary }
    }
    df_noicd <- subset_icd(df, want_icd = FALSE)
    if (!is.null(df_noicd)){
      res_noicd <- analyze_dataset(df_noicd, scenario_label = "NO_ICD")
      if (!is.null(res_noicd$hr)) { safe_write_csv(res_noicd$hr, file.path(TAB_DIR, "sens_COX_NO_ICD.csv")); ALL_HR[["NO_ICD"]] <- res_noicd$hr }
      if (!is.null(res_noicd$summary)) { safe_write_csv(res_noicd$summary, file.path(TAB_DIR, "sens_COX_SUM_NO_ICD.csv")); ALL_SUM[["NO_ICD"]] <- res_noicd$summary }
    }
    
    # Stacked comparison tables
    HR_ALL <- if (length(ALL_HR)) do.call(rbind, ALL_HR) else NULL
    SUM_ALL <- if (length(ALL_SUM)) do.call(rbind, ALL_SUM) else NULL
    safe_write_csv(HR_ALL, file.path(TAB_DIR, "sens_COX_ALL.csv"))
    safe_write_csv(SUM_ALL, file.path(TAB_DIR, "sens_COX_SUM_ALL.csv"))
    
    message("\n[Done] Sensitivity analysis CSVs written in: ",
            normalizePath(TAB_DIR, winslash = "/", mustWork = FALSE))

    sum_all <- read.csv(file.path(TAB_DIR, "sens_COX_SUM_ALL.csv"))
    aggregate(cbind(N_CC, events, EPV) ~ scenario + age_band, sum_all, function(x) c(min=min(x), max=max(x)))
    
   
