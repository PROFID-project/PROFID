# =============================================================================
# Table X: Construction of mortality, exposure, and follow-up variables by dataset
# (gt HTML output)
# =============================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(gt)
})

OUT_QC_DIR <- "T:/Study_1/Pre-processing"
dir.create(OUT_QC_DIR, showWarnings = FALSE, recursive = TRUE)

map_dt <- data.table(
  Dataset = c("PROSE", "HELIOS", "EU-CERT", "ISRAEL"),
  `Raw mortality variable` = c(
    "Death_status",
    "Status_death_cat",
    "death",
    "Status_last"
  ),
  `Status_death definition` = c(
    "Yes→1; No→0",
    "Death→1; Alive→0",
    "Yes→1; No→0",
    "DIED→1; else→0"
  ),
  `Time_death_days source` = c(
    "Death_days",
    "DAYS2DEATH.ICD (death); DAYS2LastFU.ICD (alive)",
    "length_fu_mortality",
    "Alive_last_FU_days (aligned for deceased)"
  ),
  `Raw exposure variable` = c(
    "inappshock",
    "inappropriate_shock",
    "inap_shock",
    "Inapp_shock_1st"
  ),
  `Status_FIS definition` = c(
    "Yes→1; No→0; if missing but FU available → assume No",
    "Yes→1; No→0",
    "Yes→1; No→0",
    "YES→1; NO→0"
  ),
  `Time_FIS_days source` = c(
    "t_inappshock (exposed); Death_days (unexposed)",
    "DAYS2_inappropriate_shock.ICD (exposed); DAYS2LastQuery.ICD (unexposed)",
    "length_fu_inap_shock",
    "Inapp_shock_days (exposed); Alive_last_FU_days (unexposed)"
  ),
  `Follow-up definition` = c(
    "t_followup_days = max(Death_days, days_to_app_shock, t_inappshock)",
    "t_followup_days = Time_death_days (time since ICD implant)",
    "Registry follow-up duration",
    "t_followup_days = Alive_last_FU_days"
  )
)

gt_tbl <- gt(map_dt) |>
  tab_header(
    title = "Table X. Construction of mortality, exposure, and follow-up variables by dataset"
  ) |>
  tab_options(
    table.font.size = px(11),
    data_row.padding = px(3)
  )

gtsave(
  gt_tbl,
  filename = file.path(OUT_QC_DIR, "table_variable_construction_by_dataset.html")
)
