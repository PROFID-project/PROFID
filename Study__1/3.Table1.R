
###############################################################################
# TABLE 1 — Publication-ready baseline characteristics by inappropriate shocks
# UPDATED (robust bin recode + analysis-set N + clear footnote)
#
# Output columns:
#   Characteristics | Overall | Shock | No shock | Missing, n (%)[Inapp. shock missing: ...] | p value*
#
# Key rules (reviewer-proof):
#   - Overall/Group columns are computed on ANALYSIS SET (non-missing inapp. shock)
#   - Missing column is per-variable missingness within the ANALYSIS SET
#   - Exposure missingness is shown ONLY in the Missing column header
#   - p values are COMPLETE-CASE per variable (computed on rows with non-missing
#     exposure AND non-missing variable)
#   - Display formats:
#       * Categorical (incl. binary): n (%)
#       * Continuous: mean ± SD if approx symmetric, else median [IQR]
###############################################################################


# -------------------------------------------------------------------------
# 0) CONFIG
# -------------------------------------------------------------------------
IN_RDS  <- "T:/FINAL ICD COHORT/standardised_data1.rds"
OUTDIR  <- "T:/FINAL ICD COHORT/table1_outputs"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

GRP_VAR     <- "Status_FIS"      # exposure (0/1/NA)
BOWLEY_THR  <- 0.2               # symmetry threshold (Bowley skewness)
DIGITS_CONT <- 1                 # digits for continuous summaries

# Treat these as missing codes (esp. for bin_* vars)
SPECIAL_MISS <- c(9, 99, 999, -1, -9)

# Option: hide variables with extremely sparse availability
MIN_NONMISS_PROP <- 0.01  # 1% available minimum to show row

# -------------------------------------------------------------------------
# 1) LOAD RAW + EXPOSURE MISSINGNESS (HEADER ONLY)
# -------------------------------------------------------------------------
d_raw <- readRDS(IN_RDS)
setDT(d_raw)
stopifnot(GRP_VAR %in% names(d_raw))
N_TOTAL    <- nrow(d_raw)
N_EXP_MISS <- d_raw[is.na(get(GRP_VAR)), .N]
P_EXP_MISS <- 100 * N_EXP_MISS / N_TOTAL

# Analysis set = exposure observed
d <- d_raw[!is.na(get(GRP_VAR))]
d[, (GRP_VAR) := as.integer(get(GRP_VAR) > 0)]  # strict 0/1

n0 <- d[get(GRP_VAR) == 0, .N]
n1 <- d[get(GRP_VAR) == 1, .N]
N_ANALYSIS <- n0 + n1

# Header labels (Overall N reflects analysis set)
COL_OVERALL <- sprintf("Overall (N=%d)", N_ANALYSIS)
COL_1       <- sprintf("Patients with inappropriate shock (N=%d)", n1)
COL_0       <- sprintf("Patients without inappropriate shock (N=%d)", n0)

# Missing column header includes ONLY exposure missingness, renamed
COL_MISS <- sprintf("Missing, n (%%) [Inapp. shock missing: %d (%.1f%%)]",
                    N_EXP_MISS, P_EXP_MISS)

cat(sprintf("\nRAW cohort N=%d; Inapp. shock missing=%d (%.1f%%)\n", N_TOTAL, N_EXP_MISS, P_EXP_MISS))
cat(sprintf("ANALYSIS set: shock=%d, no shock=%d (N=%d)\n\n", n1, n0, N_ANALYSIS))

# -------------------------------------------------------------------------
# 2) ROBUST BINARY CLEANING (fixes cancer/diabetes weirdness)
# -------------------------------------------------------------------------
BIN_VARS <- grep("^bin_", names(d), value = TRUE)

to_bin01 <- function(x) {
  # factor -> character
  if (is.factor(x)) x <- as.character(x)
  
  # character handling
  if (is.character(x)) {
    xx <- toupper(trimws(x))
    
    # obvious missing strings
    xx[xx %in% c("", "NA", "N/A", "NULL", "UNKNOWN", "UNK", "MISSING")] <- NA_character_
    
    # yes-ish / no-ish
    xx[xx %in% c("YES","Y","TRUE","T","PRESENT","POSITIVE","HX","HISTORY","HAS","WITH")] <- "1"
    xx[xx %in% c("NO","N","FALSE","F","ABSENT","NEGATIVE","NONE","WITHOUT")] <- "0"
    
    suppressWarnings(xn <- as.numeric(xx))
    # if numeric conversion works, use numeric; otherwise keep as-is (will become NA below)
    x <- xn
  }
  
  if (is.logical(x)) x <- as.integer(x)
  
  if (is.numeric(x) || is.integer(x)) {
    # apply special missing codes
    x[x %in% SPECIAL_MISS] <- NA
    
    # handle common 1/2 coding: 1=yes, 2=no
    ux <- sort(unique(x[!is.na(x)]))
    if (length(ux) > 0 && all(ux %in% c(1,2))) {
      return(ifelse(is.na(x), NA_integer_, ifelse(x == 1, 1L, 0L)))
    }
    
    # strict 0/1
    return(ifelse(is.na(x), NA_integer_,
                  ifelse(x %in% c(0,1), as.integer(x), NA_integer_)))
  }
  
  rep(NA_integer_, length(x))
}

for (v in BIN_VARS) d[, (v) := to_bin01(get(v))]

# QC for cancer (prints to console; does not change table)
if ("bin_cancer" %in% names(d)) {
  cat("\n=== QC: bin_cancer (analysis set) ===\n")
  print(table(d$bin_cancer, useNA = "always"))
  cat("\nBy exposure group:\n")
  print(d[, .(
    n = .N,
    nonmiss = sum(!is.na(bin_cancer)),
    miss = sum(is.na(bin_cancer)),
    yes = sum(bin_cancer == 1, na.rm = TRUE),
    no  = sum(bin_cancer == 0, na.rm = TRUE)
  ), by = get(GRP_VAR)])
  cat("=== END QC ===\n\n")
}

# -------------------------------------------------------------------------
# 3) NYHA: I/II/III/IV -> NYHA_34 (III–IV vs I–II)
# -------------------------------------------------------------------------
if ("NYHA" %in% names(d)) {
  ny <- toupper(trimws(as.character(d$NYHA)))
  ny <- fifelse(ny %in% c("1","I"), "I",
                fifelse(ny %in% c("2","II"), "II",
                        fifelse(ny %in% c("3","III"), "III",
                                fifelse(ny %in% c("4","IV"), "IV", NA_character_))))
  d[, NYHA := factor(ny, levels = c("I","II","III","IV"))]
  d[, NYHA_34 := fifelse(is.na(NYHA), NA_integer_, as.integer(NYHA %in% c("III","IV")))]
}

FORCE_BINARY_YES <- c("NYHA_34")

# -------------------------------------------------------------------------
# 4) HELPERS
# -------------------------------------------------------------------------
p_fmt <- function(p) {
  if (is.na(p)) return("")
  if (p < 0.001) return("<0.001")
  sprintf("%.3f", p)
}

fmt_mean_sd <- function(x, digits = DIGITS_CONT) {
  x <- x[!is.na(x)]
  if (!length(x)) return("\u2014")
  sprintf(paste0("%.", digits, "f \u00b1 %.", digits, "f"), mean(x), sd(x))
}

fmt_median_iqr <- function(x, digits = DIGITS_CONT) {
  x <- x[!is.na(x)]
  if (!length(x)) return("\u2014")
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE, type = 2, na.rm = TRUE)
  sprintf(paste0("%.", digits, "f [%.", digits, "f, %.", digits, "f]"), q[2], q[1], q[3])
}

bowley_skew <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 20) return(NA_real_)
  q <- quantile(x, c(0.25, 0.5, 0.75), names = FALSE, type = 2, na.rm = TRUE)
  iqr <- q[3] - q[1]
  if (iqr == 0) return(NA_real_)
  (q[3] + q[1] - 2*q[2]) / iqr
}

is_symmetric_enough <- function(x, thr = BOWLEY_THR) {
  b <- bowley_skew(x)
  if (is.na(b)) return(FALSE)
  abs(b) < thr
}

is_continuous <- function(x) {
  if (!(is.numeric(x) || is.integer(x))) return(FALSE)
  ux <- length(unique(x[!is.na(x)]))
  ux >= 10
}

has_enough_data <- function(x, min_prop = MIN_NONMISS_PROP) {
  mean(!is.na(x)) >= min_prop
}

fmt_yes_n_pct <- function(x01) {
  denom <- sum(!is.na(x01))  # variable-specific complete-case denom
  if (denom == 0) return("\u2014")
  k <- sum(x01 == 1, na.rm = TRUE)
  sprintf("%d (%.1f%%)", k, 100 * k / denom)
}

miss_str <- function(x) {
  m <- sum(is.na(x))
  p <- 100 * m / length(x)
  sprintf("%d (%.1f%%)", m, p)
}

p_binary <- function(g01, x01) {
  keep <- !is.na(g01) & !is.na(x01)
  tab <- table(g01[keep], x01[keep])
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  chi <- suppressWarnings(try(chisq.test(tab, correct = FALSE), silent = TRUE))
  if (inherits(chi, "try-error") || any(chi$expected < 5)) return(fisher.test(tab)$p.value)
  chi$p.value
}

p_categorical <- function(g01, xfac) {
  keep <- !is.na(g01) & !is.na(xfac)
  xfac <- droplevels(factor(xfac[keep]))
  gg <- factor(g01[keep], levels = c(0, 1))
  if (nlevels(xfac) < 2) return(NA_real_)
  tab <- table(gg, xfac)
  if (nrow(tab) < 2 || ncol(tab) < 2) return(NA_real_)
  chi <- suppressWarnings(try(chisq.test(tab, correct = FALSE), silent = TRUE))
  if (inherits(chi, "try-error") || any(chi$expected < 5)) return(fisher.test(tab)$p.value)
  chi$p.value
}

make_row <- function(characteristics, overall, shock, noshock, miss, pval) {
  data.table(
    Characteristics = characteristics,
    Overall = overall,
    Shock = shock,
    `No shock` = noshock,
    Missing = miss,
    `p value*` = pval
  )
}

# -------------------------------------------------------------------------
# 5) LABELS
# -------------------------------------------------------------------------
pretty_bin <- c(
  bin_sex_male              = "Male, n (%)",
  bin_diabetes              = "Diabetes, n (%)",
  bin_hypertension          = "Hypertension, n (%)",
  bin_stroke_tia            = "Prior stroke/TIA, n (%)",
  bin_copd                  = "COPD, n (%)",
  bin_cancer                = "Cancer, n (%)",
  bin_af_atrial_flutter     = "Atrial flutter/fibrillation, n (%)",
  bin_nsvt                  = "NSVT, n (%)",
  bin_lbbb                  = "LBBB, n (%)",
  bin_rbbb                  = "RBBB, n (%)",
  bin_av_block              = "AV block, n (%)",
  bin_av_block_ii_or_iii    = "AV block II–III, n (%)",
  bin_pci                   = "Prior PCI, n (%)",
  bin_cabg                  = "Prior CABG, n (%)",
  bin_mi_location_anterior  = "Anterior MI, n (%)",
  bin_smoking               = "Smoking, n (%)",
  bin_beta_blockers          = "β-blocker, n (%)",
  bin_ace_inhibitor          = "ACE inhibitor, n (%)",
  bin_ace_inhibitor_arb      = "ACEi/ARB, n (%)",
  bin_arb                    = "ARB, n (%)",
  bin_diuretics              = "Diuretics, n (%)",
  bin_aldosterone_antagonist = "Aldosterone antagonist, n (%)",
  bin_calcium_antagonists    = "Calcium antagonist, n (%)",
  bin_anti_platelet          = "Antiplatelet therapy, n (%)",
  bin_anti_coagulant         = "Anticoagulant therapy, n (%)",
  bin_lipid_lowering         = "Lipid-lowering therapy, n (%)",
  bin_anti_arrhythmic_iii     = "Class III antiarrhythmic, n (%)",
  bin_digitalis_glycosides    = "Digitalis glycosides, n (%)"
)

pretty_nonbin <- c(
  Age = "Age (years)",
  NYHA_34 = "NYHA class III–IV, n (%)",
  LVEF = "LVEF (%)",
  eGFR = "eGFR",
  HR = "Heart rate",
  SBP = "Systolic BP",
  DBP = "Diastolic BP",
  BMI = "BMI",
  QRS = "QRS duration",
  QTc = "QTc",
  BUN = "BUN",
  Sodium = "Sodium",
  Potassium = "Potassium",
  Calcium = "Calcium",
  Haemoglobin = "Haemoglobin",
  CRP = "CRP",
  TSH = "TSH",
  Troponin_T = "Troponin T"
)

pretty_label <- function(v) {
  if (v %in% names(pretty_bin)) return(pretty_bin[[v]])
  if (v %in% names(pretty_nonbin)) return(pretty_nonbin[[v]])
  v
}

# -------------------------------------------------------------------------
# 6) SECTIONS (edit as needed)
# -------------------------------------------------------------------------
sections <- list(
  "Demographic" = c("Age", "bin_sex_male"),
  "Clinical parameters" = c("NYHA_34", "LVEF", "eGFR", "HR", "SBP", "DBP", "BMI"),
  "Comorbidities" = c("bin_diabetes", "bin_hypertension", "bin_stroke_tia", "bin_copd", "bin_cancer"),
  "Arrhythmia / conduction / ECG" = c("bin_af_atrial_flutter", "bin_nsvt", "bin_lbbb", "bin_rbbb",
                                      "bin_av_block", "bin_av_block_ii_or_iii", "QRS", "QTc"),
  "Medical History" = c("bin_pci", "bin_cabg", "bin_mi_location_anterior", "bin_smoking"),
  "Medications" = c("bin_beta_blockers", "bin_ace_inhibitor", "bin_ace_inhibitor_arb", "bin_arb",
                    "bin_diuretics", "bin_aldosterone_antagonist", "bin_calcium_antagonists",
                    "bin_anti_platelet", "bin_anti_coagulant", "bin_lipid_lowering",
                    "bin_anti_arrhythmic_iii", "bin_digitalis_glycosides"),
  "Laboratory" = c("BUN", "Sodium", "Potassium", "Calcium", "Haemoglobin", "CRP", "TSH", "Troponin_T")
)
sections <- lapply(sections, function(vs) vs[vs %in% names(d)])

# -------------------------------------------------------------------------
# 7) BUILD TABLE 1
# -------------------------------------------------------------------------
table1_rows <- list()
g <- d[[GRP_VAR]]

for (sec in names(sections)) {
  
  # Section header row
  table1_rows[[length(table1_rows) + 1]] <- make_row(sec, "", "", "", "", "")
  
  for (v in sections[[sec]]) {
    
    x <- d[[v]]
    if (!has_enough_data(x)) next
    
    x0 <- x[g == 0]
    x1 <- x[g == 1]
    
    label <- pretty_label(v)
    missv <- miss_str(x)  # missing within analysis set
    
    # Binary variables
    if (grepl("^bin_", v) || v %in% FORCE_BINARY_YES) {
      overall <- fmt_yes_n_pct(x)
      shock   <- fmt_yes_n_pct(x1)
      noshock <- fmt_yes_n_pct(x0)
      p <- p_binary(g, x)
      
      table1_rows[[length(table1_rows) + 1]] <- make_row(
        label, overall, shock, noshock, missv, p_fmt(p)
      )
      next
    }
    
    # Continuous variables
    if (is_continuous(x)) {
      if (is_symmetric_enough(x, thr = BOWLEY_THR)) {
        overall <- fmt_mean_sd(x)
        shock   <- fmt_mean_sd(x1)
        noshock <- fmt_mean_sd(x0)
        tt <- try(t.test(x ~ g), silent = TRUE)
        p <- if (inherits(tt, "try-error")) NA_real_ else tt$p.value
      } else {
        overall <- fmt_median_iqr(x)
        shock   <- fmt_median_iqr(x1)
        noshock <- fmt_median_iqr(x0)
        wt <- try(wilcox.test(x ~ g, exact = FALSE), silent = TRUE)
        p <- if (inherits(wt, "try-error")) NA_real_ else wt$p.value
      }
      
      table1_rows[[length(table1_rows) + 1]] <- make_row(
        label, overall, shock, noshock, missv, p_fmt(p)
      )
      next
    }
    
    # Multi-level categorical (if any)
    xf <- droplevels(as.factor(x))
    
    make_level_summary <- function(xx) {
      xx <- xx[!is.na(xx)]
      n <- length(xx)
      if (n == 0) return("\u2014")
      tt2 <- table(xx)
      paste(sprintf("%s: %d (%.1f%%)", names(tt2), as.integer(tt2),
                    100 * as.integer(tt2) / n), collapse = "; ")
    }
    
    overall <- make_level_summary(xf)
    shock   <- make_level_summary(xf[g == 1])
    noshock <- make_level_summary(xf[g == 0])
    p <- p_categorical(g, xf)
    
    table1_rows[[length(table1_rows) + 1]] <- make_row(
      label, overall, shock, noshock, missv, p_fmt(p)
    )
  }
}

table1 <- rbindlist(table1_rows, fill = TRUE)

# Indent non-section rows + mark section rows
table1[, is_section := (Overall == "" & Shock == "" & `No shock` == "" & Missing == "" & `p value*` == "")]
table1[, Characteristics_disp := ifelse(is_section, Characteristics, paste0("  ", Characteristics))]

gt_df <- table1[, .(
  Characteristics = Characteristics_disp,
  Overall, Shock, `No shock`, Missing, `p value*`,
  is_section
)]

# -------------------------------------------------------------------------
# 8) RENDER (gt) + FOOTNOTE (your requested explanations)
# -------------------------------------------------------------------------
footnote_txt <- paste(
  "*Categorical variables are reported as n (%), where % is calculated using the number of non-missing observations for that variable within each group.",
  "Continuous variables are reported as mean \u00b1 SD when approximately symmetric and as median [IQR] when skewed.",
  "Missing values are shown as n (%) within the analysis set.",
  "p values are based on complete-case analysis for each variable (i.e., using only patients with non-missing exposure and non-missing values for that variable).*"
)

gt_tbl <- gt(gt_df) |>
  cols_label(
    Characteristics = "Characteristics",
    Overall = COL_OVERALL,
    Shock = COL_1,
    `No shock` = COL_0,
    Missing = COL_MISS,
    `p value*` = "p value*"
  ) |>
  cols_hide(columns = is_section) |>
  tab_header(
    title = "Table 1. Baseline characteristics of patients by inappropriate shocks"
  ) |>
  tab_source_note(source_note = footnote_txt) |>
  cols_align(align = "left", columns = Characteristics) |>
  cols_align(align = "center", columns = c(Overall, Shock, `No shock`, Missing, `p value*`)) |>
  tab_options(
    table.font.size = 10,
    heading.title.font.size = 12,
    table.width = pct(100),
    data_row.padding = px(3)
  ) |>
  tab_style(
    style = list(cell_fill(color = "#F2F2F2"), cell_text(weight = "bold")),
    locations = cells_body(rows = gt_df$is_section %in% TRUE)
  )

gt_tbl

# -------------------------------------------------------------------------
# 9) SAVE OUTPUTS
# -------------------------------------------------------------------------
html_out <- file.path(OUTDIR, "Table1_baseline_inapp_shocks1.html")
csv_out  <- file.path(OUTDIR, "Table1_baseline_inapp_shocks1.csv")
rds_out  <- file.path(OUTDIR, "Table1_baseline_inapp_shocks1.rds")

gtsave(gt_tbl, html_out)
fwrite(table1[, -c("is_section","Characteristics_disp"), with = FALSE], csv_out)
saveRDS(table1, rds_out)

cat("\nSaved:\n")
cat(" - ", html_out, "\n", sep = "")
cat(" - ", csv_out, "\n", sep = "")
cat(" - ", rds_out, "\n", sep = "")
###############################################################################

