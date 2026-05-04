###############################################################################

# TABLE 1 — Baseline characteristics by inappropriate ICD shock status

#

# Input:  master_clean_dataset.rds 
# Output: Table1_baseline.html/.csv/.rds

#

# Structure:

#   Characteristics | Overall (N=3084) | Shocked (n=67) | No shock (n=3017) |

#   Missing n (%)

#

# Rules:

#   - NO p-values (observational study, severely unbalanced groups 67 vs 3017)

#   - Missing column shows per-variable missingness within analysis cohort

#   - Exposure (Status_FIS) missing NOT shown — already 0 by construction

#   - Binary:     n (%) — denominator is non-missing observations in that group

#   - Continuous: mean +/- SD throughout (consistent with cardiology convention)

#   - NYHA shown as full 4-level categorical variable with sub-rows

#   - Cancer excluded (>80% missing — documented in supplementary missing table)

#   - Variables >80% missing excluded per pre-specified threshold

###############################################################################



# =============================================================================

# STEP 0: Install and load required packages

# =============================================================================

packages <- c("data.table", "gt")



for (pkg in packages) {
  
  if (!requireNamespace(pkg, quietly = TRUE)) {
    
    install.packages(pkg, dependencies = TRUE)
    
  }
  
  library(pkg, character.only = TRUE)
  
}



if (!exists("as.data.table")) stop("data.table failed to load — restart R and try again")

cat("All packages loaded successfully\n")



# =============================================================================

# STEP 1: File paths

# =============================================================================

IN_RDS <- "T:/Study_1/master_clean_dataset1.rds"

OUTDIR <- "T:/Study_1/Supplementary_data"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)



# =============================================================================

# STEP 2: Settings

# =============================================================================

GRP_VAR          <- "Status_FIS"

DIGITS_CONT      <- 1

SPECIAL_MISS     <- c(9, 99, 999, -1, -9)

MIN_NONMISS_PROP <- 0.01

MISS_THRESHOLD   <- 80   # variables >=80% missing excluded from Table 1



# =============================================================================

# STEP 3: Load master clean dataset and verify exposure has no missing values

# =============================================================================

d <- as.data.table(readRDS(IN_RDS))

stopifnot(GRP_VAR %in% names(d))



n_exp_miss <- sum(is.na(d[[GRP_VAR]]))

cat(sprintf("Exposure (Status_FIS) missing: %d (expect 0)\n", n_exp_miss))

if (n_exp_miss > 0) stop("Exposure variable has missing values — check master_clean_dataset.rds")



# Ensure strict 0/1

d[, (GRP_VAR) := as.integer(get(GRP_VAR) > 0)]



N_ANALYSIS <- nrow(d)

n1 <- d[get(GRP_VAR) == 1, .N]   # shocked

n0 <- d[get(GRP_VAR) == 0, .N]   # non-shocked



cat(sprintf("Analysis cohort: N=%d (shocked n=%d, no shock n=%d)\n",
            
            N_ANALYSIS, n1, n0))



# Column header labels

COL_OVERALL <- sprintf("Overall (N=%s)",  format(N_ANALYSIS, big.mark=","))

COL_SHOCK   <- sprintf("Patients with inappropriate shock (n=%d)",  n1)

COL_NOSHOCK <- sprintf("Patients without inappropriate shock (n=%s)", format(n0, big.mark=","))

COL_MISS    <- "Missing, n (%)"



# =============================================================================

# STEP 4: Clean binary variables to strict 0/1

# =============================================================================

BIN_VARS <- grep("^bin_", names(d), value = TRUE)



to_bin01 <- function(x) {
  
  if (is.factor(x))  x <- as.character(x)
  
  if (is.character(x)) {
    
    xx <- toupper(trimws(x))
    
    xx[xx %in% c("", "NA", "N/A", "NULL", "UNKNOWN", "UNK", "MISSING")] <- NA_character_
    
    xx[xx %in% c("YES","Y","TRUE","T","PRESENT","POSITIVE",
                 
                 "HX","HISTORY","HAS","WITH")]                           <- "1"
    
    xx[xx %in% c("NO","N","FALSE","F","ABSENT","NEGATIVE",
                 
                 "NONE","WITHOUT")]                                       <- "0"
    
    suppressWarnings(xn <- as.numeric(xx))
    
    x <- xn
    
  }
  
  if (is.logical(x)) x <- as.integer(x)
  
  if (is.numeric(x) || is.integer(x)) {
    
    x[x %in% SPECIAL_MISS] <- NA
    
    ux <- sort(unique(x[!is.na(x)]))
    
    # Handle 1/2 coding (1=yes, 2=no)
    
    if (length(ux) > 0 && all(ux %in% c(1, 2)))
      
      return(ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L)))
    
    # Strict 0/1
    
    return(ifelse(is.na(x), NA_integer_,
                  
                  ifelse(x %in% c(0, 1), as.integer(x), NA_integer_)))
    
  }
  
  rep(NA_integer_, length(x))
  
}



for (v in BIN_VARS) d[, (v) := to_bin01(get(v))]



# =============================================================================

# STEP 5: NYHA — standardise to ordered factor I/II/III/IV

# =============================================================================

if ("NYHA" %in% names(d)) {
  
  ny <- toupper(trimws(as.character(d$NYHA)))
  
  ny <- fifelse(ny %in% c("1","I"),   "I",
                
                fifelse(ny %in% c("2","II"),  "II",
                        
                        fifelse(ny %in% c("3","III"), "III",
                                
                                fifelse(ny %in% c("4","IV"),  "IV", NA_character_))))
  
  d[, NYHA := factor(ny, levels = c("I","II","III","IV"))]
  
}



# =============================================================================

# STEP 6: Exclude variables exceeding 80% missingness threshold

# =============================================================================

pct_missing <- function(x) 100 * sum(is.na(x)) / length(x)



vars_excluded_80 <- names(d)[sapply(names(d), function(v) pct_missing(d[[v]]) >= MISS_THRESHOLD)]

cat(sprintf("\nVariables excluded (>=80%% missing): %s\n",
            
            paste(vars_excluded_80, collapse = ", ")))



# =============================================================================

# STEP 7: Helper functions

# =============================================================================



# Continuous: mean +/- SD (used throughout — consistent with cardiology convention)

fmt_mean_sd <- function(x, digits = DIGITS_CONT) {
  
  x <- x[!is.na(x)]
  
  if (!length(x)) return("\u2014")
  
  sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"), mean(x), sd(x))
  
}



# Note: median/IQR and Bowley skewness not used in Table 1

# Skewness checks and log transformations are modelling decisions only

# Table 1 reports raw values as mean +/- SD throughout (cardiology convention)



# Is variable continuous (10+ unique non-missing values)

is_continuous <- function(x) {
  
  if (!(is.numeric(x) || is.integer(x))) return(FALSE)
  
  length(unique(x[!is.na(x)])) >= 10
  
}



# Enough non-missing data to include row

has_enough_data <- function(x, min_prop = MIN_NONMISS_PROP) {
  
  mean(!is.na(x)) >= min_prop
  
}



# Binary: n (%) — denominator is non-missing in that group

fmt_yes_n_pct <- function(x01) {
  
  denom <- sum(!is.na(x01))
  
  if (denom == 0) return("\u2014")
  
  k <- sum(x01 == 1, na.rm = TRUE)
  
  sprintf("%d (%.1f%%)", k, 100 * k / denom)
  
}



# Missing: n (%) out of total analysis cohort

miss_str <- function(x) {
  
  m <- sum(is.na(x))
  
  p <- 100 * m / length(x)
  
  sprintf("%d (%.1f%%)", m, p)
  
}



# Build a single table row

make_row <- function(label, overall, shock, noshock, miss,
                     
                     is_section = FALSE, is_sublevel = FALSE) {
  
  data.table(
    
    Characteristics = label,
    
    Overall         = overall,
    
    Shock           = shock,
    
    `No shock`      = noshock,
    
    Missing         = miss,
    
    is_section      = is_section,
    
    is_sublevel     = is_sublevel
    
  )
  
}



# NYHA: build multi-row categorical block (one row per level)

make_nyha_rows <- function(x, g, miss_val) {
  
  rows <- list()
  
  # Header row — no values, just label
  
  rows[[1]] <- make_row("NYHA class", "", "", "", miss_val,
                        
                        is_section = FALSE, is_sublevel = FALSE)
  
  # One sub-row per level
  
  for (lv in c("I","II","III","IV")) {
    
    x_bin    <- as.integer(!is.na(x) & x == lv)
    
    x_bin[is.na(x)] <- NA_integer_
    
    
    
    x0 <- x_bin[g == 0]
    
    x1 <- x_bin[g == 1]
    
    
    
    rows[[length(rows) + 1]] <- make_row(
      
      label   = lv,
      
      overall = fmt_yes_n_pct(x_bin),
      
      shock   = fmt_yes_n_pct(x1),
      
      noshock = fmt_yes_n_pct(x0),
      
      miss    = "",           # missing shown once on header row
      
      is_section  = FALSE,
      
      is_sublevel = TRUE
      
    )
    
  }
  
  rows
  
}



# =============================================================================

# STEP 8: Variable labels

# =============================================================================

pretty_bin <- c(
  
  bin_sex_male               = "Male, n (%)",
  
  bin_diabetes               = "Diabetes, n (%)",
  
  bin_hypertension           = "Hypertension, n (%)",
  
  bin_stroke_tia             = "Prior stroke/TIA, n (%)",
  
  bin_copd                   = "COPD, n (%)",
  
  bin_af_atrial_flutter      = "Atrial fibrillation/flutter, n (%)",
  
  bin_nsvt                   = "NSVT, n (%)",
  
  bin_lbbb                   = "LBBB, n (%)",
  
  bin_rbbb                   = "RBBB, n (%)",
  
  bin_av_block               = "AV block, n (%)",
  
  bin_av_block_ii_or_iii     = "AV block II-III, n (%)",
  
  bin_pci                    = "Prior PCI, n (%)",
  
  bin_cabg                   = "Prior CABG, n (%)",
  
  bin_mi_location_anterior   = "Anterior MI, n (%)",
  
  bin_smoking                = "Smoking, n (%)",
  
  bin_beta_blockers          = "Beta-blocker, n (%)",
  
  bin_ace_inhibitor          = "ACE inhibitor, n (%)",
  
  bin_ace_inhibitor_arb      = "ACEi/ARB, n (%)",
  
  bin_arb                    = "ARB, n (%)",
  
  bin_diuretics              = "Diuretics, n (%)",
  
  bin_aldosterone_antagonist = "Aldosterone antagonist, n (%)",
  
  bin_calcium_antagonists    = "Calcium antagonist, n (%)",
  
  bin_anti_platelet          = "Antiplatelet therapy, n (%)",
  
  bin_anti_coagulant         = "Anticoagulant therapy, n (%)",
  
  bin_lipid_lowering         = "Lipid-lowering therapy, n (%)",
  
  bin_anti_arrhythmic_iii    = "Class III antiarrhythmic, n (%)",
  
  bin_digitalis_glycosides   = "Digitalis glycosides, n (%)"
  
)



pretty_nonbin <- c(
  
  Age         = "Age (years)",
  
  LVEF        = "LVEF",
  
  eGFR        = "eGFR (mL/min/1.73m\u00b2)",
  
  HR          = "Heart rate (bpm)",
  
  SBP         = "Systolic BP (mmHg)",
  
  DBP         = "Diastolic BP (mmHg)",
  
  BMI         = "BMI (kg/m\u00b2)",
  
  QRS         = "QRS duration (ms)",
  
  QTc         = "QTc (ms)",
  
  BUN         = "BUN",
  
  Sodium      = "Sodium (mmol/L)",
  
  Potassium   = "Potassium (mmol/L)",
  
  Calcium     = "Calcium (mmol/L)",
  
  Haemoglobin = "Haemoglobin (g/dL)",
  
  CRP         = "CRP",
  
  TSH         = "TSH",
  
  Troponin_T  = "Troponin T"
  
)



pretty_label <- function(v) {
  
  if (v %in% names(pretty_bin))    return(pretty_bin[[v]])
  
  if (v %in% names(pretty_nonbin)) return(pretty_nonbin[[v]])
  
  v
  
}



# =============================================================================

# STEP 9: Define sections

# NOTE: Cancer excluded (>80% missing — documented in supplementary missing table)

#       NYHA handled separately as multi-level categorical

# =============================================================================

sections <- list(
  
  "Demographic" = c(
    
    "Age", "bin_sex_male"
    
  ),
  
  "Clinical parameters" = c(
    
    "LVEF", "eGFR", "HR", "SBP", "DBP", "BMI"
    
  ),
  
  "Comorbidities" = c(
    
    "bin_diabetes", "bin_hypertension", "bin_stroke_tia", "bin_copd"
    
  ),
  
  "Arrhythmia / Conduction / ECG" = c(
    
    "bin_af_atrial_flutter", "bin_nsvt", "bin_lbbb", "bin_rbbb",
    
    "bin_av_block", "bin_av_block_ii_or_iii", "QRS", "QTc"
    
  ),
  
  "Medical history" = c(
    
    "bin_pci", "bin_cabg", "bin_mi_location_anterior", "bin_smoking"
    
  ),
  
  "Medications" = c(
    
    "bin_beta_blockers", "bin_ace_inhibitor", "bin_ace_inhibitor_arb",
    
    "bin_arb", "bin_diuretics", "bin_aldosterone_antagonist",
    
    "bin_calcium_antagonists", "bin_anti_platelet", "bin_anti_coagulant",
    
    "bin_lipid_lowering", "bin_anti_arrhythmic_iii", "bin_digitalis_glycosides"
    
  ),
  
  "Laboratory" = c(
    
    "BUN", "Sodium", "Potassium", "Calcium",
    
    "Haemoglobin", "CRP", "TSH", "Troponin_T"
    
  )
  
)



# Remove variables that:

# (a) are not present in the dataset

# (b) exceed 80% missingness threshold

sections <- lapply(sections, function(vs) {
  
  vs[vs %in% names(d) & !vs %in% vars_excluded_80]
  
})



# =============================================================================

# STEP 10: Build Table 1

# =============================================================================

cat("\nBuilding Table 1...\n")



table1_rows <- list()

g <- d[[GRP_VAR]]



for (sec in names(sections)) {
  
  
  
  # Section header
  
  table1_rows[[length(table1_rows) + 1]] <-
    
    make_row(sec, "", "", "", "", is_section = TRUE)
  
  
  
  # NYHA inserted first in Clinical parameters section
  
  if (sec == "Clinical parameters" && "NYHA" %in% names(d)) {
    
    nyha_miss <- miss_str(d$NYHA)
    
    nyha_rows <- make_nyha_rows(d$NYHA, g, nyha_miss)
    
    table1_rows <- c(table1_rows, nyha_rows)
    
  }
  
  
  
  for (v in sections[[sec]]) {
    
    
    
    x <- d[[v]]
    
    if (!has_enough_data(x)) next
    
    
    
    x0    <- x[g == 0]
    
    x1    <- x[g == 1]
    
    label <- pretty_label(v)
    
    missv <- miss_str(x)
    
    
    
    # --- Binary ---
    
    if (grepl("^bin_", v)) {
      
      table1_rows[[length(table1_rows) + 1]] <- make_row(
        
        label,
        
        fmt_yes_n_pct(x),
        
        fmt_yes_n_pct(x1),
        
        fmt_yes_n_pct(x0),
        
        missv
        
      )
      
      next
      
    }
    
    
    
    # --- Continuous: mean +/- SD throughout ---
    
    if (is_continuous(x)) {
      
      table1_rows[[length(table1_rows) + 1]] <- make_row(
        
        label,
        
        fmt_mean_sd(x),
        
        fmt_mean_sd(x1),
        
        fmt_mean_sd(x0),
        
        missv
        
      )
      
      next
      
    }
    
    
    
    # --- Other categorical (fallback) ---
    
    # Not expected for current variable set but kept for safety
    
    xf <- droplevels(as.factor(x))
    
    lvls <- levels(xf)
    
    # Header row
    
    table1_rows[[length(table1_rows) + 1]] <- make_row(
      
      label, "", "", "", missv,
      
      is_section = FALSE, is_sublevel = FALSE
      
    )
    
    for (lv in lvls) {
      
      x_bin    <- as.integer(!is.na(x) & as.character(x) == lv)
      
      x_bin[is.na(x)] <- NA_integer_
      
      table1_rows[[length(table1_rows) + 1]] <- make_row(
        
        lv,
        
        fmt_yes_n_pct(x_bin),
        
        fmt_yes_n_pct(x_bin[g == 1]),
        
        fmt_yes_n_pct(x_bin[g == 0]),
        
        "",
        
        is_sublevel = TRUE
        
      )
      
    }
    
  }
  
}



table1 <- rbindlist(table1_rows, fill = TRUE)

table1[is.na(is_section),  is_section  := FALSE]

table1[is.na(is_sublevel), is_sublevel := FALSE]



# Display labels: section headers bold (handled by gt styling),

# sub-levels indented, regular variables indented once

table1[, Characteristics_disp := fcase(
  
  is_section  == TRUE,  Characteristics,
  
  is_sublevel == TRUE,  paste0("      ", Characteristics),  # deeper indent for NYHA levels
  
  default               = paste0("   ", Characteristics)
  
)]



gt_df <- table1[, .(
  
  Characteristics = Characteristics_disp,
  
  Overall,
  
  Shock,
  
  `No shock`,
  
  Missing,
  
  is_section,
  
  is_sublevel
  
)]



# =============================================================================

# STEP 11: Render table using gt

# =============================================================================

footnote_txt <- paste(
  
  "Categorical variables reported as n (%);",
  
  "denominator is the number of non-missing observations within each group.",
  
  "Continuous variables reported as mean \u00b1 SD.",
  
  "Missing values reflect per-variable missingness within the analysis cohort",
  
  sprintf("(N=%s).", format(N_ANALYSIS, big.mark=",")),
  
  "Variables with \u226580% missing data were excluded from this table;",
  
  "a full missing data summary is provided in the supplementary material.",
  
  "Missing values in covariates were addressed using multiple imputation",
  
  "by chained equations (MICE) in all regression analyses.",
  
  "Inappropriate ICD shock status was resolved using pre-specified harmonisation rules;",
  
  "no missing values remained for the exposure variable in the analysis cohort.",
  
  "p-values are not reported: this is an observational study with severely unbalanced",
  
  sprintf("group sizes (n=%d vs n=%s);", n1, format(n0, big.mark=",")),
  
  "in this setting p-values conflate statistical significance with sample size",
  
  "and do not assist in assessing confounding (Senn, 1994; Austin and Brunner, 2004)."
  
)



gt_tbl <- gt(gt_df) |>
  
  cols_label(
    
    Characteristics = "Characteristics",
    
    Overall         = COL_OVERALL,
    
    Shock           = COL_SHOCK,
    
    `No shock`      = COL_NOSHOCK,
    
    Missing         = COL_MISS
    
  ) |>
  
  cols_hide(columns = c(is_section, is_sublevel)) |>
  
  tab_header(
    
    title    = "Table 1. Baseline characteristics by inappropriate ICD shock status",
    
    subtitle = md("Values are n (%) for categorical variables and mean \u00b1 SD for continuous variables.")
    
  ) |>
  
  tab_source_note(source_note = footnote_txt) |>
  
  cols_align(align = "left",   columns = Characteristics) |>
  
  cols_align(align = "center", columns = c(Overall, Shock, `No shock`, Missing)) |>
  
  tab_options(
    
    table.font.size            = 10,
    
    heading.title.font.size    = 12,
    
    heading.subtitle.font.size = 10,
    
    table.width                = pct(100),
    
    data_row.padding           = px(3)
    
  ) |>
  
  # Section header rows — grey background, bold
  
  tab_style(
    
    style     = list(cell_fill(color = "#F2F2F2"),
                     
                     cell_text(weight = "bold")),
    
    locations = cells_body(rows = gt_df$is_section == TRUE)
    
  ) |>
  
  # NYHA header row (label only, no values) — light italic
  
  tab_style(
    
    style     = cell_text(style = "italic"),
    
    locations = cells_body(
      
      rows    = gt_df$is_section == FALSE &
        
        gt_df$is_sublevel == FALSE &
        
        gt_df$Overall == "" &
        
        gt_df$Missing != ""
      
    )
    
  ) |>
  
  # NYHA sub-level rows — slightly lighter text
  
  tab_style(
    
    style     = cell_text(color = "#333333"),
    
    locations = cells_body(rows = gt_df$is_sublevel == TRUE)
    
  )



# =============================================================================

# STEP 12: Save outputs

# =============================================================================

html_out <- file.path(OUTDIR, "Table1_baseline.html")

csv_out  <- file.path(OUTDIR, "Table1_baseline.csv")

rds_out  <- file.path(OUTDIR, "Table1_baseline.rds")



gtsave(gt_tbl, html_out)

fwrite(table1[, .(Characteristics, Overall, Shock, `No shock`, Missing)], csv_out)

saveRDS(table1, rds_out)



cat("\nDone. Files saved:\n")

cat("  Table 1 HTML: ", html_out, "\n", sep = "")

cat("  Table 1 CSV:  ", csv_out,  "\n", sep = "")

cat("  Table 1 RDS:  ", rds_out,  "\n", sep = "")



# =============================================================================

# STEP 13: Supplementary Missing Data Table


# Simple descriptive table: variable, N missing, % missing

# Sorted by % missing descending

# Outcome/exposure variables in separate section at top

# No handling categories — those are documented in the MICE script

# =============================================================================

cat("\nBuilding Supplementary Missing Data Table...\n")



# Outcome / exposure variables — listed separately at top

OUTCOME_EXPOSURE_VARS <- intersect(
  
  c("Status_death", "Time_death_days", "Status_FIS",
    
    "Time_FIS_days", "Survival_time", "t_followup_days"),
  
  names(d)
  
)



# All other variables

OTHER_VARS <- setdiff(names(d), OUTCOME_EXPOSURE_VARS)



# Pretty label for missingness table — strips bin_ prefix and cleans up name

pretty_miss_label <- function(v) {
  
  lbl <- pretty_label(v)               # use existing pretty_bin / pretty_nonbin lookup
  
  # strip ", n (%)" suffix if present — not needed in a missingness table
  
  lbl <- gsub(",\\s*n\\s*\\(%\\)", "", lbl)
  
  trimws(lbl)
  
}



# Build missingness row for a single variable

build_miss_row <- function(v, section_label) {
  
  x <- d[[v]]
  
  data.table(
    
    Section     = section_label,
    
    Variable    = pretty_miss_label(v),   # human-readable label, no bin_ prefix
    
    N_missing   = sum(is.na(x)),
    
    Pct_missing = round(100 * sum(is.na(x)) / length(x), 1)
    
  )
  
}



# Build clean variable list for missingness table:

# 1. Use pretty_bin and pretty_nonbin keys — human readable, no bin_ prefix

# 2. Exclude log1p companions — derived modelling variables, not clinical vars

# 3. Exclude NYHA_34 — derived binary, NYHA already included

# 4. No other bin_ variables or raw dataset internals



# Variables to exclude explicitly

EXCLUDE_FROM_MISS <- c(
  
  "NYHA_34",                                    # derived from NYHA
  
  GRP_VAR,                                      # exposure — shown in outcome section
  
  grep("_log1p$", names(d), value = TRUE)       # all log1p companions
  
)



# Covariate variables: only those with pretty labels (known clinical variables)

cov_vars_bin    <- names(pretty_bin)[names(pretty_bin) %in% names(d)]

cov_vars_nonbin <- names(pretty_nonbin)[names(pretty_nonbin) %in% names(d)]



# Add NYHA once using raw variable

if ("NYHA" %in% names(d)) cov_vars_nonbin <- c("NYHA", cov_vars_nonbin)



# Final covariate list — no duplicates, no log1p, no derived variables

ALL_COV_VARS <- unique(c(cov_vars_bin, cov_vars_nonbin))

ALL_COV_VARS <- ALL_COV_VARS[
  
  ALL_COV_VARS %in% names(d) &
    
    !ALL_COV_VARS %in% EXCLUDE_FROM_MISS
  
]



# Build and sort each section by % missing descending

miss_oe    <- rbindlist(lapply(OUTCOME_EXPOSURE_VARS, build_miss_row, "Outcome / Exposure"))

miss_other <- rbindlist(lapply(ALL_COV_VARS,          build_miss_row, "Covariate"))



miss_oe    <- miss_oe[order(-Pct_missing)]

miss_other <- miss_other[order(-Pct_missing)]



miss_full  <- rbind(miss_oe, miss_other)



# Console summary

cat(sprintf("Total variables:              %d\n", nrow(miss_full)))

cat(sprintf("Variables with any missing:   %d\n", miss_full[N_missing > 0, .N]))

cat(sprintf("Variables with 0%% missing:   %d\n", miss_full[N_missing == 0, .N]))



# -----------------------------------------------------------------------------

# Build gt table

# -----------------------------------------------------------------------------

miss_gt <- gt(miss_full, groupname_col = "Section") |>
  
  cols_label(
    
    Variable    = "Variable",
    
    N_missing   = "N Missing",
    
    Pct_missing = "% Missing"
    
  ) |>
  
  tab_header(
    
    title    = md("**Supplementary Table. Missing data summary**"),
    
    subtitle = md(sprintf(
      
      "Per-variable missingness in the analysis cohort (N = %s), sorted by %% missing (descending).",
      
      format(N_ANALYSIS, big.mark = ",")
      
    ))
    
  ) |>
  
  # Highlight variables >= 80% missing
  
  tab_style(
    
    style     = list(cell_fill(color = "#fff3cd"),
                     
                     cell_text(weight = "bold")),
    
    locations = cells_body(rows = Pct_missing >= MISS_THRESHOLD)
    
  ) |>
  
  # Section group headers
  
  tab_style(
    
    style     = list(cell_fill(color = "#f8f9fa"),
                     
                     cell_text(weight = "bold")),
    
    locations = cells_row_groups()
    
  ) |>
  
  tab_footnote(
    
    footnote  = md("Variables highlighted in **yellow** had \u226580% missing data and were excluded from Table 1 and from multiple imputation. Full details of missing data handling are provided in the statistical analysis plan."),
    
    locations = cells_column_labels(columns = Pct_missing)
    
  ) |>
  
  cols_align(align = "center", columns = c(N_missing, Pct_missing)) |>
  
  cols_align(align = "left",   columns = Variable) |>
  
  opt_row_striping() |>
  
  tab_options(
    
    table.font.size              = 10,
    
    heading.title.font.size      = 12,
    
    heading.subtitle.font.size   = 10,
    
    table.width                  = pct(50),
    
    data_row.padding             = px(3),
    
    row_group.font.weight        = "bold",
    
    column_labels.font.weight    = "bold"
    
  )



# -----------------------------------------------------------------------------

# Save supplementary missing data table

# -----------------------------------------------------------------------------

miss_html <- file.path(OUTDIR, "Supplementary_Missing_Data.html")

miss_csv  <- file.path(OUTDIR, "Supplementary_Missing_Data.csv")



gtsave(miss_gt, miss_html)

fwrite(miss_full[, .(Section, Variable, N_missing, Pct_missing)], miss_csv)



cat("\nSupplementary missing data table saved:\n")

cat("  HTML: ", miss_html, "\n", sep = "")

cat("  CSV:  ", miss_csv,  "\n", sep = "")





















