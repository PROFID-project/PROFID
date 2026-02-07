###############################################################################
# CRUDE INCIDENCE (per 100 PY) + POWER (approx) — PRIMARY Cox (all-cause death)
# - Reads CDM-ready analysis RDS
# - QC: fixes FIS events occurring after follow-up (optional)
# - Computes crude incidence for death + FIS (Poisson exact 95% CI)
# - Computes post-hoc power for HR grid + MDHR (Schoenfeld approximation)
# - Saves CSV + HTML (gt) outputs
###############################################################################


# Load required packages
suppressPackageStartupMessages({
  library(data.table)
  install.packages("powerSurvEpi")
  library(powerSurvEpi)
  install.packages("survminer")
  library(survminer)
  
})


# Load dat2 (output from standardise_for_table1_and_cox)
ds1 <-readRDS("T:/FINAL ICD COHORT/standardised_data1.rds")

summary(ds1$Time_death_days)
# Work on analysis copy
dst_qc <- copy(ds1)

if (!is.data.table(dst_qc)) {
  setDT(dst_qc)
}

# Flag cases where FIS time exceeds follow-up time
dst_qc[, flag_FIS_after_fu :=
         Status_FIS == 1 &
         !is.na(Time_FIS_days) &
         !is.na(Survival_time) &
         Time_FIS_days > (Survival_time*30.44)
]
dst_qc[, .N, by = flag_FIS_after_fu]

# Fix ONLY those 12 cases

dst_qc[flag_FIS_after_fu == TRUE, `:=`(
  
  Status_FIS = 0,
  
  Time_FIS_days = NA_real_
  
)]



# (Optional) drop helper flag

dst_qc[, flag_FIS_after_fu := NULL]

cat(sprintf("  ✓ Dataset loaded: %s observations, %s variables\n", 
            format(nrow(dst_qc), big.mark = ","), 
            ncol(ds1)))


#--------------------------------------------------

# Work on a copy (do NOT touch original dst)

#--------------------------------------------------

ds <- copy(dst_qc)



# --- (1) Function: Poisson exact 95% CI for incidence rates ---
calculate_crude_rate <- function(
    dt,
    time_var,
    event_var,
    event_label,
    time_unit = 365.25,   # DAYS -> YEARS
    multiplier = 100      # per 100 person-years
) {
  ds <- as.data.table(dt)
  
  total_time_yrs   <- sum(ds[[time_var]],  na.rm = TRUE) / time_unit
  total_events     <- sum(ds[[event_var]], na.rm = TRUE)
  rate_per_unit    <- (total_events / total_time_yrs) * multiplier
  
  # Poisson exact 95% CI
  ci_lower <- (qchisq(0.025, 2 * total_events) / 2) / total_time_yrs * multiplier
  ci_upper <- (qchisq(0.975, 2 * (total_events + 1)) / 2) / total_time_yrs * multiplier
  
  data.table(
    outcome      = event_label,
    n_total      = nrow(ds),
    n_events     = total_events,
    percent      = round(total_events / nrow(ds) * 100, 1),
    person_years = round(total_time_yrs, 1),
    rate_100py   = round(rate_per_unit, 2),
    ci_lower     = round(ci_lower, 2),
    ci_upper     = round(ci_upper, 2),
    ci_method    = "Poisson exact"
  )
}

# --- (2) Calculate crude death incidence ---
death_crude <- calculate_crude_rate(
  dt          = ds,
  time_var    = "Time_death_days",
  event_var   = "Status_death",
  event_label = "All-cause mortality",
  time_unit   = 365.25,
  multiplier  = 100
)

# --- (3) Save CSV (optional, good practice) ---
data.table::fwrite(death_crude, file.path(OUTDIR, "results_crude_death.csv"))

# --- (4) Save HTML with gt (temporary conversion) ---
gtsave(
  gt::gt(death_crude) |>
    gt::cols_label(
      outcome      = "Outcome",
      n_total      = "N",
      n_events     = "Deaths",
      percent      = "% deaths",
      person_years = "Person-years",
      rate_100py   = "Crude rate (per 100 PY)",
      ci_lower     = "95% CI (Lower)",
      ci_upper     = "95% CI (Upper)",
      ci_method    = "CI method"
    ),
  filename = file.path(OUTDIR, "results_crude_death.html"),
  inline_css = TRUE
)

#------------------------------------------------------------------------------
# Inappropriate ICD Therapy distribution and incidence
#------------------------------------------------------------------------------

cat("  --- Inappropriate ICD Therapy Distribution ---\n")

# Distribution
exposure_summary <- ds[, .(N = .N, 
                           Percent = round(.N / nrow(ds) * 100, 1)), 
                       by = Status_FIS]
exposure_summary[, Label := fifelse(Status_FIS == 0, 
                                    "No Inappropriate ICD Therapy",
                                    "Inappropriate ICD Therapy")]
setcolorder(exposure_summary, c("Status_FIS", "Label", "N", "Percent"))

print(exposure_summary, row.names = FALSE)

# Calculate incidence rate for inappropriate therapy
cat("\n  --- Inappropriate ICD Therapy Incidence Rate ---\n")

# Keep only rows with valid time-at-risk
ds <- ds[!is.na(Time_FIS_days) & Time_FIS_days > 0]

# Total person-time (years)
total_time_fis <- sum(ds$Time_FIS_days) / 365.25

# Total FIS events
total_fis_events <- sum(ds$Status_FIS == 1, na.rm = TRUE)

# Crude incidence rate per 100 PY
fis_rate <- (total_fis_events / total_time_fis) * 100

# Poisson exact 95% CI
fis_ci_lower <- (qchisq(0.025, 2 * total_fis_events) / 2) / total_time_fis * 100
fis_ci_upper <- (qchisq(0.975, 2 * (total_fis_events + 1)) / 2) / total_time_fis * 100

# Create comprehensive summary for export
inappropriate_therapy_incidence <- data.table(
  outcome = "Inappropriate ICD Therapy",
  n_total = nrow(ds),
  n_events = total_fis_events,
  percent_events = round(total_fis_events / nrow(ds) * 100, 1),
  person_years = round(total_time_fis, 1),
  crude_rate = round(fis_rate, 1),
  ci_lower = round(fis_ci_lower, 1),
  ci_upper = round(fis_ci_upper, 1),
  ci_method = "Poisson exact"
)

inappropriate_therapy_incidence[, rate_95ci := sprintf("%.1f (%.1f-%.1f)", 
                                                       crude_rate, ci_lower, ci_upper)]
inappropriate_therapy_incidence[, summary := sprintf("N = %d (%.1f%%), incidence rate %.1f per 100 person-years",
                                                     n_events, percent_events, crude_rate)]
print(inappropriate_therapy_incidence)
# Export
fwrite(exposure_summary, "results/exposure_summary.csv")
fwrite(inappropriate_therapy_incidence, "results/inappropriate_therapy_incidence.csv")


gtsave(
  gt::gt(exposure_summary),
  filename = "results/exposure_summary.html",
  inline_css = TRUE
)

gtsave(
  gt::gt(inappropriate_therapy_incidence),
  filename = "results/inappropriate_therapy_incidence.html",
  inline_css = TRUE
)


# ============================================================
# POWER ANALYSIS for Cox Regression (Manual Calculation)
# Does NOT rely on powerSurvEpi package functions
# Based on Schoenfeld (1983) and Hsieh & Lavori (2000)
# ============================================================

library(data.table)
library(gt)

cat("Power Analysis (manual calculation for Cox regression):\n\n")

# Parameters
alpha <- 0.05
target_power <- 0.80
HR_grid <- c(1.3, 1.4, 1.5)

# Data from your dataset
D_events <- sum(ds$Status_death == 1, na.rm = TRUE)
p_exposed <- mean(ds$Status_FIS == 1, na.rm = TRUE)

cat(sprintf("  Outcome events (D) = %d\n", D_events))
cat(sprintf("  Exposure prevalence (p) = %.4f (%.2f%%)\n\n", 
            p_exposed, p_exposed * 100))

# ============================================================
# FUNCTION: Calculate power for Cox regression
# Based on Schoenfeld formula
# ============================================================
calc_cox_power <- function(n_events, prop_exposed, HR, alpha = 0.05, two_sided = TRUE) {
  # Critical value
  z_alpha <- ifelse(two_sided, qnorm(1 - alpha/2), qnorm(1 - alpha))
  
  # Non-centrality parameter
  # delta = log(HR) * sqrt(n_events * p * (1-p))
  delta <- log(HR) * sqrt(n_events * prop_exposed * (1 - prop_exposed))
  
  # Power calculation
  if(two_sided) {
    power <- pnorm(delta - z_alpha) + pnorm(-delta - z_alpha)
  } else {
    power <- pnorm(delta - z_alpha)
  }
  
  return(power)
}

# ============================================================
# FUNCTION: Calculate minimum detectable HR
# ============================================================
calc_mdhr <- function(n_events, prop_exposed, target_power = 0.80, 
                      alpha = 0.05, two_sided = TRUE) {
  z_alpha <- ifelse(two_sided, qnorm(1 - alpha/2), qnorm(1 - alpha))
  z_beta <- qnorm(target_power)
  
  # Solve for log(HR)
  # log(HR) = (z_alpha + z_beta) / sqrt(n_events * p * (1-p))
  log_hr <- (z_alpha + z_beta) / sqrt(n_events * prop_exposed * (1 - prop_exposed))
  
  mdhr <- exp(log_hr)
  return(mdhr)
}

# ============================================================
# POST-HOC POWER for hypothesized HRs
# ============================================================
cat("POST-HOC POWER ANALYSIS\n")
cat("=======================\n")

power_vec <- sapply(HR_grid, function(hr) {
  calc_cox_power(
    n_events = D_events,
    prop_exposed = p_exposed,
    HR = hr,
    alpha = alpha,
    two_sided = TRUE
  )
})

power_tab <- data.table(
  HR = HR_grid,
  posthoc_power = round(power_vec, 3),
  interpretation = ifelse(power_vec >= 0.80, "Adequate", "Inadequate")
)

print(power_tab)
cat("\n")

# ============================================================
# MINIMUM DETECTABLE HR
# ============================================================
mdhr <- calc_mdhr(
  n_events = D_events,
  prop_exposed = p_exposed,
  target_power = target_power,
  alpha = alpha,
  two_sided = TRUE
)

cat(sprintf("MINIMUM DETECTABLE HR (MDHR)\n"))
cat(sprintf("=============================\n"))
cat(sprintf("At %.0f%% power and alpha=%.2f: HR = %.2f\n\n", 
            target_power * 100, alpha, mdhr))

# Alternative: also show power needed to detect HR=1.3
hr_1_3_power <- calc_cox_power(D_events, p_exposed, 1.3, alpha)
cat(sprintf("Note: Current study has %.1f%% power to detect HR=1.3\n", 
            hr_1_3_power * 100))

# ============================================================
# ADDITIONAL POWER CURVE (optional)
# ============================================================
hr_curve <- seq(1.1, 2.0, by = 0.1)
power_curve <- sapply(hr_curve, function(hr) {
  calc_cox_power(D_events, p_exposed, hr, alpha)
})

power_curve_tab <- data.table(
  HR = hr_curve,
  Power = round(power_curve, 3)
)

cat("\nFULL POWER CURVE\n")
cat("================\n")
print(power_curve_tab)

# ============================================================
# SAVE OUTPUTS
# ============================================================
# Create output directory if needed
if(!exists("OUTDIR")) {
  OUTDIR <- "output"
}
if(!dir.exists(OUTDIR)) {
  dir.create(OUTDIR, recursive = TRUE)
}

# Save CSVs
fwrite(power_tab, file.path(OUTDIR, "power_posthoc_table.csv"))
fwrite(power_curve_tab, file.path(OUTDIR, "power_curve_table.csv"))

power_summary <- data.table(
  outcome_events_D = D_events,
  exposure_prevalence_p = round(p_exposed, 4),
  alpha = alpha,
  target_power = target_power,
  MDHR = round(mdhr, 2),
  power_for_HR_1.3 = round(hr_1_3_power, 3)
)
fwrite(power_summary, file.path(OUTDIR, "power_mdhr_summary.csv"))

# Save HTML tables
gt::gtsave(
  gt::gt(power_tab) |>
    gt::tab_header(
      title = "Post-hoc Power Analysis",
      subtitle = sprintf("Based on %d events and %.1f%% exposure prevalence", 
                         D_events, p_exposed * 100)
    ) |>
    gt::cols_label(
      HR = "Hazard Ratio",
      posthoc_power = "Power",
      interpretation = "Interpretation"
    ) |>
    gt::tab_source_note(
      source_note = "Power calculated using Schoenfeld (1983) formula for Cox regression"
    ),
  filename = file.path(OUTDIR, "power_posthoc_table.html"),
  inline_css = TRUE
)

gt::gtsave(
  gt::gt(power_summary) |>
    gt::tab_header(
      title = "Minimum Detectable Hazard Ratio (MDHR)",
      subtitle = sprintf("Two-sided test, alpha=%.2f", alpha)
    ) |>
    gt::cols_label(
      outcome_events_D = "Events (D)",
      exposure_prevalence_p = "Exposure prev. (p)",
      alpha = "Alpha",
      target_power = "Target power",
      MDHR = "MDHR",
      power_for_HR_1.3 = "Power for HR=1.3"
    ) |>
    gt::fmt_number(
      columns = c(exposure_prevalence_p, power_for_HR_1.3),
      decimals = 3
    ),
  filename = file.path(OUTDIR, "power_mdhr_summary.html"),
  inline_css = TRUE
)

gt::gtsave(
  gt::gt(power_curve_tab) |>
    gt::tab_header(
      title = "Power Curve",
      subtitle = "Power to detect different hazard ratios"
    ) |>
    gt::cols_label(
      HR = "Hazard Ratio",
      Power = "Statistical Power"
    ) |>
    gt::data_color(
      columns = Power,
      colors = scales::col_numeric(
        palette = c("red", "orange", "lightgreen", "darkgreen"),
        domain = c(0, 1)
      )
    ),
  filename = file.path(OUTDIR, "power_curve_table.html"),
  inline_css = TRUE
)

cat("\n✓ Output files saved to:", OUTDIR, "\n")
cat("  - power_posthoc_table.csv/.html\n")
cat("  - power_curve_table.csv/.html\n")
cat("  - power_mdhr_summary.csv/.html\n")






