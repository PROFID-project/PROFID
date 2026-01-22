###############################################
# Temporal Consistency (ICD and Non-ICD cohorts)
# - Table generation
###############################################

# Packages
req <- c("data.table","dplyr","gt")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Input file
df <- fread("NR_cmr.csv")

# Temporal consistency summary
temporal_mi_summary <- df %>%
  filter(!is.na(Time_index_MI_CHD)) %>%
  group_by(DB) %>%
  summarise(
    n = n(),
    median_months = median(Time_index_MI_CHD),
    Q1_months = quantile(Time_index_MI_CHD, 0.25),
    Q3_months = quantile(Time_index_MI_CHD, 0.75),
    min_months = min(Time_index_MI_CHD),
    max_months = max(Time_index_MI_CHD),
    .groups = "drop"
  ) %>%
  mutate(
    median_months = round(median_months, 1),
    Q1_months = round(Q1_months, 1),
    Q3_months = round(Q3_months, 1),
    min_months = round(min_months, 1),
    max_months = round(max_months, 1),
    IQR_months = paste0(Q1_months, "–", Q3_months),
    flag_remote_MI = max_months > 12
  ) %>%
  select(
    DB,
    n,
    median_months,
    IQR_months,
    min_months,
    max_months,
    flag_remote_MI
  )

# GT table
gt_temporal_mi <- gt(temporal_mi_summary) %>%
  cols_label(
    DB = "Database",
    n = md("**n**"),
    median_months = md("**Median (months)**"),
    IQR_months = md("**IQR (months)**"),
    min_months = md("**Min**"),
    max_months = md("**Max**"),
    flag_remote_MI = md("**Max >12 months**")
  ) %>%
  cols_align(
    align = "center",
    columns = c(
      n,
      median_months,
      IQR_months,
      min_months,
      max_months,
      flag_remote_MI
    )
  ) %>%
  tab_header(
    title = md("**Time since most recent MI/CHD at baseline — Non-ICD cohort (≤35%) with CMR**"),
    subtitle = md(
      "Centre-level distributions used as a proxy temporal consistency check"
    )
  ) %>%
  tab_source_note(
    md(
      "`Time_index_MI_CHD` represents time (in months) from most recent MI/CHD to baseline assessment. 
      This table is descriptive and does not represent exact post-MI baseline timing."
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_column_labels(everything())
  ) %>%
  tab_options(
    table.font.size = px(13),
    data_row.padding = px(4)
  )
gtsave(gt_temporal_mi, "NR_cmr_time_since_MI_baseline_by_centre.html")