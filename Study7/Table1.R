###############################################
# Table 1 (ICD vs Non-ICD cohorts)
# - Uses HZ.eda_describe_data() per-cohort
# - Robust union merge of SDC summaries
# - Outputs CSV + nicely formatted Word/HTML with gt
###############################################

# Packages ---------------------------------------------------------------
req <- c("data.table","dplyr","openxlsx","stringr","gt")
inst <- setdiff(req, rownames(installed.packages()))
if (length(inst)) install.packages(inst, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# HZ functions -----------------------------------------------------------
# If you have the HZ functions in a file, source it; otherwise assume
# they're already in the session.
if (!exists("HZ.eda_describe_data")) {
  hz_path <- "hz-basic-summary-statistics.R"
  if (file.exists(hz_path)) {
    source(hz_path)
  } else {
    stop("HZ.eda_describe_data() not found. Source 'hz-basic-summary-statistics.R' or paste the HZ functions into the session.")
  }
}

# Input files ------------------------------------------------------------
files <- list(
  ICD              = "ICD.csv",
  NonICD_reduced   = "NonICD_reduced.csv",   # Non-ICD, LVEF <= 35%
  NonICD_preserved = "NonICD_preserved.csv"  # Non-ICD, LVEF > 35%
)
missing_files <- names(files)[!file.exists(unlist(files))]
if (length(missing_files)) stop("Missing input file(s): ", paste(missing_files, collapse = ", "))

# Read and summarise each dataset with HZ (SDC) --------------------------
read_df <- function(p) as.data.frame(data.table::fread(p, data.table = FALSE))
dfs <- lapply(files, read_df)
hz_summ <- lapply(dfs, function(df) HZ.eda_describe_data(df)$summary.sdc)
# keep as data.table
hz_summ <- lapply(hz_summ, as.data.table)
names(hz_summ) <- names(files)

# Pull N (sample size) for table headers ---------------------------------
get_N_from_summary <- function(dt) {
  as.integer(gsub("^[<≥]*([0-9]+).*$", "\\1", dt[Variable == "N" & Category == "N", Summary][1]))
}
Ns <- sapply(hz_summ, get_N_from_summary)

# Build union table of (Variable, Category) across datasets --------------
all_keys <- unique(rbindlist(lapply(hz_summ, function(x) x[, .(Variable, Category)]), use.names = TRUE))
union_tab <- copy(all_keys)
for (nm in names(hz_summ)) {
  tmp <- hz_summ[[nm]][, .(Variable, Category, Summary)]
  setnames(tmp, "Summary", nm)
  union_tab <- merge(union_tab, tmp, by = c("Variable","Category"), all.x = TRUE, sort = FALSE)
}

# Helper: choose a variable by aliases (first that exists) ---------------
pick_var <- function(aliases, available_vars) {
  hit <- aliases[aliases %in% available_vars]
  if (length(hit)) hit[1] else NA_character_
}

# Aliases to tolerate naming differences ---------------------------------
aliases <- list(
  Age   = c("Age","age"),
  Sex   = c("Sex","sex"),
  BMI   = c("BMI","bmi","Body mass index"),
  Diabetes = c("Diabetes","diabetes"),
  Hypertension = c("Hypertension","hypertension"),
  Smoking = c("Smoking","smoking"),
  eGFR  = c("eGFR","eGFR_(mL/min/1.73_m2)","Estimated_glomerular_filtration_rate"),
  LBBB  = c("LBBB","Left bundle branch block","Left_bundle_branch_block"),
  ACE_inhibitor_ARB = c("ACE_inhibitor_ARB","ACE_inhibitor","ARB","ACE_inhibitors","Angiotensin-converting enzyme inhibitors and angiotensin receptor blockers"),
  Beta_blockers = c("Beta_blockers","Beta blockers"),
  Diuretics     = c("Diuretics"),
  Anti_platelet = c("Anti_platelet","Antiplatelet drugs","Antiplatelet"),
  Anticoagulant = c("Anticoagulant","Oral anticoagulants","Anticoagulants"),
  Lipid_lowering= c("Lipid_lowering","Lipid-lowering medication","Statins"),
  PCI   = c("PCI","Percutaneous coronary intervention","Prior percutaneous coronary intervention"),
  CABG  = c("CABG","Prior coronary bypass graft surgery","Coronary artery bypass graft surgery"),
  LVEF  = c("LVEF","LV_EF","Left ventricular ejection fraction (%)","Left ventricular ejection fraction"),
  .dummy = character(0)
)

# Find canonical variable names actually present in your union table ------
available_vars <- unique(union_tab$Variable)
canon <- lapply(aliases, pick_var, available_vars = available_vars)
# Remove missing alias hits
canon <- canon[!vapply(canon, is.na, logical(1))]

# Prefer categories: continuous -> Mean (SD), fallback Median (IQR) ------
get_cat_pref <- function(var) {
  if (var %in% c(canon$Age, canon$BMI, canon$LVEF)) return(c("Mean (SD)","Median (IQR)"))
  if (var %in% c(canon$eGFR)) return(c("Median (IQR)","Mean (SD)"))  # like reference
  # categorical positives
  if (var %in% c(canon$Sex)) return("Male")
  if (var %in% c(canon$Diabetes, canon$Hypertension, canon$Smoking,
                 canon$LBBB, canon$ACE_inhibitor_ARB, canon$Beta_blockers,
                 canon$Diuretics, canon$Anti_platelet, canon$Anticoagulant,
                 canon$Lipid_lowering, canon$PCI, canon$CABG)) return("Yes")
  return(c("Mean (SD)","Median (IQR)"))
}

# Fetch a single display row for a variable (best available category) -----
fetch_row <- function(var, cat_pref) {
  # try preferred categories in order
  cats <- if (is.character(cat_pref)) cat_pref else unlist(cat_pref)
  hit <- NULL
  for (c in cats) {
    hit <- union_tab[Variable == var & Category == c]
    if (nrow(hit) > 0) break
  }
  if (is.null(hit) || nrow(hit) == 0) {
    # if nothing, return NA row (will show "-")
    return(data.frame(Variable = var, Category = cats[1], ICD = NA, NonICD_reduced = NA, NonICD_preserved = NA, stringsAsFactors = FALSE))
  }
  as.data.frame(hit[, .(Variable, Category, ICD, NonICD_reduced, NonICD_preserved)])
}

# Compute overall Missing n (%) for a variable across the three datasets --
parse_count <- function(x) {
  # Extract leading number from strings like "287 (3.8%)", "<5 (...)", "≥10 (...)"
  if (is.na(x) || !nzchar(x)) return(NA_real_)
  as.numeric(str_replace(x, "^[<≥]*([0-9]+).*$", "\\1"))
}
overall_missing <- function(var) {
  miss_row <- union_tab[Variable == var & Category == "Missing"]
  if (nrow(miss_row) == 0) return("-")
  n1 <- parse_count(miss_row$ICD)
  n2 <- parse_count(miss_row$NonICD_reduced)
  n3 <- parse_count(miss_row$NonICD_preserved)
  total_n <- sum(Ns, na.rm = TRUE)
  missing_n <- sum(c(n1, n2, n3), na.rm = TRUE)
  if (is.na(missing_n)) return("-")
  pct <- if (total_n > 0) round(100 * missing_n / total_n, 1) else NA_real_
  paste0(missing_n, " (", pct, "%)")
}

# Build the table blocks to match the reference groups -------------------
# (Section order & variables mirror the screenshot)
sections <- list(
  "Demographics" = c(canon$Age, canon$Sex),
  "Medical history" = c(canon$PCI, canon$CABG, canon$Smoking, canon$Diabetes),
  "Clinical characteristics" = c(canon$BMI, canon$Hypertension, canon$eGFR),
  "Electrocardiography" = c(canon$LBBB),
  "Medication" = c(canon$ACE_inhibitor_ARB, canon$Beta_blockers, canon$Diuretics,
                   canon$Anti_platelet, canon$Anticoagulant, canon$Lipid_lowering),
  "Echocardiography" = c(canon$LVEF)
)

# Assemble rows
rows <- list()
for (grp in names(sections)) {
  for (v in sections[[grp]]) {
    if (is.null(v) || is.na(v)) next
    row <- fetch_row(v, get_cat_pref(v))
    row$Group <- grp
    row$`Missing data n (%)` <- overall_missing(v)
    rows[[length(rows)+1]] <- row
  }
}
tbl <- bind_rows(rows) |>
  mutate(across(c(ICD, NonICD_reduced, NonICD_preserved), ~ ifelse(is.na(.x) | .x == "", "-", .x))) |>
  select(Group, Variable, Category, ICD, NonICD_reduced, NonICD_preserved, `Missing data n (%)`)

# Write CSV (machine-readable) ------------------------------------------
out_csv <- "Table1.csv"
fwrite(tbl, out_csv)
message("Saved CSV: ", out_csv)

# Pretty table with gt (Word + HTML) ------------------------------------
hdr_icd  <- paste0("ICD patients<br/>(n = ", Ns["ICD"], ")")
hdr_nr   <- paste0("Non-ICD patients<br/>≤35% (n = ", Ns["NonICD_reduced"], ")")
hdr_np   <- paste0("Non-ICD patients<br/>>35% (n = ", Ns["NonICD_preserved"], ")")

gt_tbl <- tbl |>
  gt(groupname_col = "Group") |>
  tab_header(title = md("**Table 1: Baseline characteristics**"),
             subtitle = md("Excluding cardiac magnetic resonance imaging parameters")) |>
  cols_label(
    Variable = md("**Variable**"),
    Category = md("**Statistic / Level**"),
    ICD = md(hdr_icd),
    NonICD_reduced = md(hdr_nr),
    NonICD_preserved = md(hdr_np),
    `Missing data n (%)` = md("**Missing data n (%)**")
  ) |>
  tab_style(style = cell_text(weight = "bold"),
            locations = cells_row_groups()) |>
  cols_align(align = "center", columns = c(ICD, NonICD_reduced, NonICD_preserved, `Missing data n (%)`)) |>
  tab_options(table.font.size = px(13),
              data_row.padding = px(4))

gtsave(gt_tbl, "Table1.html")

# Console preview --------------------------------------------------------
print(gt_tbl)
