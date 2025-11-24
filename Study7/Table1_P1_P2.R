###############################################################################
# Produces publication-style Table 1 (Phase1 + Phase2)
# - Uses HZ.eda_describe_data() for SDC-aware summaries
# - Compact (publication) tables for Phase1 and Phase2
# - CSV + HTML outputs
###############################################################################

# Packages
req <- c("data.table","dplyr","gt","htmltools","DT","htmlwidgets","glue","stringr","openxlsx")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Source HZ functions (edit if path differs)
hz_path <- "hz-basic-summary-statistics.R"
if (file.exists(hz_path)) {
  source(hz_path)
} else {
  if (!exists("HZ.eda_describe_data")) stop("HZ.eda_describe_data() not found. Put hz-basic-summary-statistics.R in working dir or source it.")
}

# Input files - edit paths if stored elsewhere
file_ICD  <- "ICD.csv"
file_NR   <- "NonICD_reduced.csv"
file_NP   <- "NonICD_preserved.csv"

stopifnot(file.exists(file_ICD), file.exists(file_NR), file.exists(file_NP))

# Read raw datasets
df_ICD  <- data.table::fread(file_ICD, data.table = FALSE)
df_NR   <- data.table::fread(file_NR,  data.table = FALSE)
df_NP   <- data.table::fread(file_NP,  data.table = FALSE)

# Helper: detect HasMRI variants
find_mri_col <- function(df) {
  possible <- c("HasMRI","Has_MRI","hasMRI","HAS_MRI","has_mri","Has_Mri")
  present  <- intersect(possible, colnames(df))
  if (length(present) == 0) return(NA_character_)
  present[1]
}
col_icd_mri <- find_mri_col(df_ICD)
col_nr_mri  <- find_mri_col(df_NR)
col_np_mri  <- find_mri_col(df_NP)

# Compute HZ summaries (SDC)
s_icd <- HZ.eda_describe_data(df_ICD)$summary.sdc
s_nr  <- HZ.eda_describe_data(df_NR)$summary.sdc
s_np  <- HZ.eda_describe_data(df_NP)$summary.sdc

# Helper to build merged union table from HZ outputs
build_union <- function(s_icd, s_nr, s_np) {
  s_icd <- as.data.table(s_icd); s_nr <- as.data.table(s_nr); s_np <- as.data.table(s_np)
  all_keys <- unique(rbindlist(list(s_icd[,.(Variable,Category)], s_nr[,.(Variable,Category)], s_np[,.(Variable,Category)]), use.names = TRUE))
  master <- copy(all_keys)
  master <- merge(master, s_icd[, .(Variable, Category, Summary)], by = c("Variable","Category"), all.x = TRUE, sort = FALSE)
  setnames(master, "Summary", "ICD")
  master <- merge(master, s_nr[, .(Variable, Category, Summary)], by = c("Variable","Category"), all.x = TRUE, sort = FALSE)
  setnames(master, "Summary", "NonICD_reduced")
  master <- merge(master, s_np[, .(Variable, Category, Summary)], by = c("Variable","Category"), all.x = TRUE, sort = FALSE)
  setnames(master, "Summary", "NonICD_preserved")
  setcolorder(master, c("Variable","Category","ICD","NonICD_reduced","NonICD_preserved"))
  return(master)
}

table1_phase1_union <- build_union(s_icd, s_nr, s_np)

# Pick preferred rows for compact table
# Preference: Mean (SD) for continuous, Median (IQR) when preferred, for categorical keep positive (Male/Yes)
pick_row <- function(merged, var, prefer = c("Mean (SD)","Median (IQR)","Male","Yes")) {
  subs <- merged[Variable == var]
  if (nrow(subs) == 0) return(NULL)
  for (p in prefer) {
    if (any(subs$Category == p)) return(subs[Category == p,])
  }
  # fallback: first non-NA summary row
  nm <- subs[!(is.na(ICD) & is.na(NonICD_reduced) & is.na(NonICD_preserved))]
  if (nrow(nm) > 0) return(nm[1,])
  return(subs[1,])
}

# Define the groups & variable 
core_groups <- list(
  "Demographics" = c("Age","Sex"),
  "Medical history" = c("PCI","CABG","Smoking","Diabetes"),
  "Clinical characteristics" = c("BMI","Hypertension","eGFR"),
  "Electrocardiography" = c("LBBB"),
  "Medication" = c("ACE_inhibitor_ARB","Beta_blockers","Diuretics","Anti_platelet","Anticoagulant","Lipid_lowering"),
  "Echocardiography" = c("LVEF")
)

# Build compact rows for Phase 1
compact_rows_p1 <- list()
for (grp in names(core_groups)) {
  for (v in core_groups[[grp]]) {
    r <- pick_row(table1_phase1_union, v, prefer = c("Mean (SD)","Median (IQR)","Male","Yes"))
    if (!is.null(r)) {
      r$Group <- grp
      compact_rows_p1[[length(compact_rows_p1)+1]] <- r
    }
  }
}
table1_phase1_compact <- if (length(compact_rows_p1)>0) rbindlist(compact_rows_p1, fill = TRUE) else data.table::data.table()

# Add Missing data n (%) column: find row where Category == "Missing" for each variable and sum across cohorts
compute_missing_col <- function(merged, var) {
  mrow <- merged[Variable == var & Category == "Missing"]
  if (nrow(mrow)==0) return("-")
  parse_leading_num <- function(x) {
    if (is.na(x) || !nzchar(x)) return(NA_real_)
    as.numeric(gsub("^[<≥]*([0-9]+).*$", "\\1", x))
  }
  n1 <- parse_leading_num(mrow$ICD[1]); n2 <- parse_leading_num(mrow$NonICD_reduced[1]); n3 <- parse_leading_num(mrow$NonICD_preserved[1])
  total_missing_n <- sum(c(n1,n2,n3), na.rm = TRUE)
  total_N <- as.numeric(gsub("^[^0-9]*([0-9]+).*$","\\1", table1_phase1_union[Variable=="N" & Category=="N", ICD][1])) # ICD N
  # safer: fetch Ns from HZ summaries for each cohort
  fetch_N <- function(s) {
    if (nrow(s[s$Variable=="N" & s$Category=="N",])>0) {
      as.numeric(gsub("^[^0-9]*([0-9]+).*$","\\1", s[s$Variable=="N" & s$Category=="N", Summary][1]))
    } else NA_real_
  }
  N_icd <- fetch_N(as.data.table(s_icd)); N_nr <- fetch_N(as.data.table(s_nr)); N_np <- fetch_N(as.data.table(s_np))
  totalN <- sum(c(N_icd, N_nr, N_np), na.rm = TRUE)
  if (is.na(totalN) || totalN == 0) return(as.character(total_missing_n))
  pct <- round(100 * total_missing_n / totalN, 1)
  paste0(total_missing_n, " (", pct, "%)")
}

# Attach missing column for each compact row
if (nrow(table1_phase1_compact) > 0) {
  table1_phase1_compact[, `Missing data n (%)` := sapply(Variable, function(v) compute_missing_col(table1_phase1_union, v))]
  setcolorder(table1_phase1_compact, c("Group","Variable","Category","ICD","NonICD_reduced","NonICD_preserved","Missing data n (%)"))
  data.table::fwrite(table1_phase1_compact, "Table1_Phase1.csv")
} else {
  message("No Phase1 compact rows built; check variable names in core_groups.")
}

# -------------------------
# Phase 2 (CMR subset) - filter HasMRI == 1
# -------------------------
filter_has_mri <- function(df, colname) {
  if (is.na(colname)) return(df[0,,drop=FALSE])
  col <- df[[colname]]
  ok <- (!is.na(col)) & ((col == 1) | (tolower(as.character(col)) %in% c("1","true","t","yes","y")))
  df[which(ok), , drop = FALSE]
}
df_ICD_cmr <- filter_has_mri(df_ICD, col_icd_mri)
df_NR_cmr  <- filter_has_mri(df_NR,  col_nr_mri)
df_NP_cmr  <- filter_has_mri(df_NP,  col_np_mri)
message("CMR subset sizes — ICD: ", nrow(df_ICD_cmr), "  NR: ", nrow(df_NR_cmr), "  NP: ", nrow(df_NP_cmr))

# HZ summaries for CMR subset (if empty, create empty dt)
s_icd_cmr <- if (nrow(df_ICD_cmr)>0) HZ.eda_describe_data(df_ICD_cmr)$summary.sdc else data.table::data.table()
s_nr_cmr  <- if (nrow(df_NR_cmr)>0)  HZ.eda_describe_data(df_NR_cmr)$summary.sdc  else data.table::data.table()
s_np_cmr  <- if (nrow(df_NP_cmr)>0)  HZ.eda_describe_data(df_NP_cmr)$summary.sdc  else data.table::data.table()

table1_phase2_union <- build_union(s_icd_cmr, s_nr_cmr, s_np_cmr)

# Build compact rows for Phase2
compact_rows_p2 <- list()
for (grp in names(core_groups)) {
  for (v in core_groups[[grp]]) {
    r <- pick_row(table1_phase2_union, v, prefer = c("Mean (SD)","Median (IQR)","Male","Yes"))
    if (!is.null(r)) {
      r$Group <- grp
      compact_rows_p2[[length(compact_rows_p2)+1]] <- r
    }
  }
}
table1_phase2_compact <- if(length(compact_rows_p2)>0) rbindlist(compact_rows_p2, fill=TRUE) else data.table::data.table()
if (nrow(table1_phase2_compact)>0) {
  table1_phase2_compact[, `Missing data n (%)` := sapply(Variable, function(v) compute_missing_col(table1_phase2_union, v))]
  setcolorder(table1_phase2_compact, c("Group","Variable","Category","ICD","NonICD_reduced","NonICD_preserved","Missing data n (%)"))
  data.table::fwrite(table1_phase2_compact, "Table1_Phase2.csv")
} else {
  message("Phase2 compact empty (no CMR subset or core variables missing).")
}

# -------------------------
# Create readable GT tables
# -------------------------
make_gt_html <- function(dt_compact, title, subtitle, out_html) {
  if (nrow(dt_compact) == 0) {
    message("No rows to render for ", title)
    return(NULL)
  }
  # clean Group -> header rows
  gt_tbl <- gt(dt_compact, groupname_col = "Group") %>%
    tab_header(title = md(glue::glue("**{title}**")), subtitle = md(subtitle)) %>%
    cols_label(Variable = md("**Variable**"), Category = md("**Statistic / Level**"),
               ICD = md(glue::glue("ICD patients<br/>(n = {as.integer(ifelse(!is.na(nrow(df_ICD)), nrow(df_ICD), NA))})")),
               NonICD_reduced = md(glue::glue("Non-ICD patients<br/>≤35% (n = {as.integer(ifelse(!is.na(nrow(df_NR)), nrow(df_NR), NA))})")),
               NonICD_preserved = md(glue::glue("Non-ICD patients<br/>>35% (n = {as.integer(ifelse(!is.na(nrow(df_NP)), nrow(df_NP), NA))})")),
               `Missing data n (%)` = md("**Missing data n (%)**")) %>%
    tab_style(style = cell_text(weight = "bold"), locations = cells_row_groups()) %>%
    tab_options(table.font.size = px(13), data_row.padding = px(4)) %>%
    cols_align(align = "center", columns = vars(ICD, NonICD_reduced, NonICD_preserved, `Missing data n (%)`))
  # save HTML
  gtsave(gt_tbl, out_html)
  message("Wrote HTML:", out_html)
  return(out_html)
}

# Production of phase HTMLs
html1 <- make_gt_html(as.data.frame(table1_phase1_compact), "Table 1: Phase 1: Baseline characteristics",
                      "Excluding cardiac magnetic resonance imaging parameters", "Table1_Phase1.html")
html2 <- make_gt_html(as.data.frame(table1_phase2_compact), "Table 1: Phase 2: Baseline characteristics with CMR (HasMRI == 1)",
                      "Including cardiac magnetic resonance imaging parameters", "Table1_Phase2.html")

message("All done. Files written: Table1_Phase1.csv/html, Table1_Phase2.csv/html")
