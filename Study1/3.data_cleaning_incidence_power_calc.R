###############################################################################
# CRUDE INCIDENCE (per 100 PY) + POWER (approx) — PRIMARY Cox (all-cause death)
#   - Data cleaning and cohort derivation (Table S1)
#   - Crude incidence rates for mortality and inappropriate ICD shock
#   - Post-hoc power analysis and MDHR (Schoenfeld approximation)
###############################################################################

# Required packages ---------------------------------------------------------
pkgs <- c(
  "tidyverse","dplyr" ,"survival", "data.table", "grid","gridExtra",
  "openxlsx", "readxl"
)

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))
install.packages("gt")
library(gt)


# -------------------------------------------------------------------------
# 1) Set file paths
# -------------------------------------------------------------------------
INFILE  <- "T:/FINAL ICD COHORT/Transformed_data1.rds"
OUTFILE <- "T:/Study_1/master_clean_dataset1.rds"
OUTDIR  <- "T:/Study_1/Supplementary_data"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# -------------------------------------------------------------------------
# 2) Load data
# -------------------------------------------------------------------------
ds0 <- as.data.table(readRDS(INFILE))
cat(sprintf("✓ Dataset loaded: %s observations, %s variables\n",
            format(nrow(ds0), big.mark=","), ncol(ds0)))
summary(ds0$Time_FIS_days)

cat("\n--- Quick check: Time_death_days summary (days) ---\n")
print(summary(ds0$Time_death_days))
print(summary)(ds0$t_followup_days)
table(ds0$Status_FIS)
#Time_death_days and follow up time is same.

# -------------------------------------------------------------------------
# 3) Data quality checks (reporting only — nothing is changed here)

# We look for common data problems and report how many patients
# are affected by each issue. No data is modified in this step.
# Problems found here are properly handled in Step 4.
# Note :*FIS = first inappropriate shock
# -------------------------------------------------------------------------

ds_qc <- as.data.table(copy(ds0))

# Flag 1: FIS occurring after follow-up (days-based)
ds_qc[, flag_FIS_after_fu :=
        Status_FIS == 1L &
        !is.na(Time_FIS_days) &
        !is.na(Time_death_days) &
        Time_FIS_days > Time_death_days]

# Flag 2: FIS time exceeds Survival time(months-based check)
ds_qc[, flag_FIS_exceeds_fu :=
        Status_FIS == 1L &
        !is.na(Time_FIS_days) &
        !is.na(Survival_time) &
        Time_FIS_days > (Survival_time * 30.44)]

# Flag 3: Exposed but FIS time missing
ds_qc[, flag_FIS_ambiguous :=
        Status_FIS == 1L &
        is.na(Time_FIS_days)]

# Flag 4: Both Status_FIS and Time_FIS_days missing
ds_qc[, flag_both_missing :=
        is.na(Status_FIS) &
        is.na(Time_FIS_days)]

# Flag 5: Dead but missing follow-up time
#         Critical: these deaths cannot be counted — if > 0 then
#         mortality rate is slightly underestimated
ds_qc[, flag_dead_missing_fu :=
        is.na(Time_death_days) &
        !is.na(Status_death) &
        Status_death == 1L]

# Flag 6: Valid follow-up time but missing Status_death
#         Standard practice: treat as censored (assumed alive at last contact)
ds_qc[, flag_status_missing_time_valid :=
        is.na(Status_death) &
        !is.na(Time_death_days) &
        Time_death_days > 0]

# Report all flags
cat("\nFlag 1 — FIS after follow-up (days-based):\n")
print(ds_qc[, .N, by = flag_FIS_after_fu])

cat("\nFlag 2 — FIS time exceeds follow-up (months-based, x30.44):\n")
print(ds_qc[, .N, by = flag_FIS_exceeds_fu])

cat("\nFlag 3 — Exposed but FIS time missing:\n")
print(ds_qc[, .N, by = flag_FIS_ambiguous])

cat("\nFlag 4 — Both Status_FIS and Time_FIS_days missing:\n")
print(ds_qc[, .N, by = flag_both_missing])

cat("\nFlag 5 — Dead but missing follow-up time:\n")
print(ds_qc[, .N, by = flag_dead_missing_fu])
cat("  NOTE: Cannot contribute person-time without follow-up duration.\n")
cat("        If deaths > 0, mortality rate is slightly underestimated.\n")

cat("\nFlag 6 — Missing Status_death but valid follow-up time:\n")
print(ds_qc[, .N, by = flag_status_missing_time_valid])
cat("  NOTE: Treated as censored (Status_death = 0) — standard practice.\n")
cat("        Unknown vital status = assumed alive at last known contact.\n")

# Cleanup QC copy
ds_qc[, c("flag_FIS_after_fu",  "flag_FIS_exceeds_fu",
          "flag_FIS_ambiguous", "flag_both_missing",
          "flag_dead_missing_fu", "flag_status_missing_time_valid") := NULL]

# -------------------------------------------------------------------------
# 4)Clean the dataset

# We apply a fixed set of pre-specified rules to handle missing values 
# in the outcome (death) and exposure (inappropriate shock) variables. 
# These rules are applied ONCE here and the cleaned dataset is saved 
# for use in ALL downstream analyses (Table 1, Cox model, etc.) 

# # OUTCOME (death) RULES: 

# O1: Remove patients with no follow-up time recorded 
# We cannot calculate person-time without knowing how long a patient was followed.
# Therefore, these patients must be excluded 
# O2: If vital status is missing but follow-up time exists → mark as alive (censored) 
# We assume the patient was alive at their last known contact 


# # EXPOSURE (inappropriate shock) RULES: 

# E1: If both shock status and shock timing are missing → classify as no shock 
# No evidence of a shock anywhere in the data; conservative assumption 
# E2: If shock timing is recorded but status is missing → classify as shocked 
# A recorded time implies an event occurred 
# E3: Patients with no shock should have no shock timing either → set to missing 
# E4: If marked as shocked but timing is unknown → REMOVE from dataset
# We cannot determine when the shock happened, so we cannot correctly classify their 
# follow-up time as exposed or unexposed.# We apply a fixed set of pre-specified rules to handle missing values # in the outcome (death) and exposure (inappropriate shock) variables. # These rules are applied ONCE here and the cleaned dataset is saved # for use in ALL downstream analyses (Table 1, Cox model, etc.) # # OUTCOME (death) RULES: # O1: Remove patients with no follow-up time recorded # → We cannot calculate person-time without knowing how long # a patient was followed these patients must be excluded # O2: If vital status is missing but follow-up time exists → mark as alive (censored) # → We assume the patient was alive at their last known contact # # EXPOSURE (inappropriate shock) RULES: # E1: If both shock status and shock timing are missing → classify as no shock # → No evidence of a shock anywhere in the data; conservative assumption # E2: If shock timing is recorded but status is missing → classify as shocked # → A recorded time implies an event occurred # E3: Patients with no shock should have no shock timing either → set to missing # E4: If marked as shocked but timing is unknown → REMOVE from dataset # → We cannot determine when the shock happened, so we cannot correctly # classify their follow-up time as exposed or unexposed # → Keeping them would make our results inconsistent # E5: If shock is recorded at or after the end of follow-up → reclassify as no shock # → A shock outside the observation window is not a valid study event # =============================================================================
# → Keeping them would make our results inconsistent 
# E5: If shock is recorded at or after the end of follow-up → reclassify as no shock 
# → A shock outside the observation window is not a valid study event 
# =============================================================================

ds_clean <- as.data.table(copy(ds0))

n_o1 <- sum(is.na(ds_clean$Time_death_days) | ds_clean$Time_death_days <= 0,
            na.rm = TRUE)
ds_clean <- ds_clean[!is.na(Time_death_days) & Time_death_days > 0]
cat(sprintf("O1 - Excluded (missing/zero follow-up time):        %d\n", n_o1))

# O2: missing vital status with valid follow-up => censored
n_o2 <- sum(is.na(ds_clean$Status_death) & !is.na(ds_clean$Time_death_days) &
              ds_clean$Time_death_days > 0)
ds_clean[is.na(Status_death) & !is.na(Time_death_days) & Time_death_days > 0,
         Status_death := 0L]
cat(sprintf("O2 - Recoded missing vital status => censored:      %d\n", n_o2))

# E1: both FIS variables missing => unexposed
# Censoring time assigned as Time_death_days (last known follow-up)
n_e1 <- sum(is.na(ds_clean$Status_FIS) & is.na(ds_clean$Time_FIS_days))
ds_clean[is.na(Status_FIS) & is.na(Time_FIS_days),
         `:=`(Status_FIS    = 0L,
              Time_FIS_days = Time_death_days)]
cat(sprintf("E1 - Both FIS missing => unexposed:                 %d\n", n_e1))

# E2: FIS time present but status missing => exposed
n_e2 <- sum(is.na(ds_clean$Status_FIS) & !is.na(ds_clean$Time_FIS_days))
ds_clean[is.na(Status_FIS) & !is.na(Time_FIS_days), Status_FIS := 1L]
cat(sprintf("E2 - FIS time present, status missing => exposed:   %d\n", n_e2))

# E4: exposed but timing unknown => HARD EXCLUSION
# Cannot place shock in time - exposure status unclassifiable
# Excluded from all analyses to ensure consistent population
n_e4 <- sum(ds_clean$Status_FIS == 1L & is.na(ds_clean$Time_FIS_days),
            na.rm = TRUE)
ds_clean <- ds_clean[!(Status_FIS == 1L & is.na(Time_FIS_days))]
cat(sprintf("E4 - Exposed, timing unknown => excluded:           %d\n", n_e4))

# E5: FIS at/after end of follow-up => recode as unexposed
# Restricted to exposed patients only
# Censoring time assigned as Time_death_days
n_e5 <- sum(ds_clean$Status_FIS == 1L & !is.na(ds_clean$Time_FIS_days) &
              ds_clean$Time_FIS_days >= ds_clean$Time_death_days, na.rm = TRUE)
ds_clean[Status_FIS == 1L & !is.na(Time_FIS_days) &
           Time_FIS_days >= Time_death_days,
         `:=`(Status_FIS    = 0L,
              Time_FIS_days = Time_death_days)]
cat(sprintf("E5 - FIS at/after follow-up => unexposed:           %d\n", n_e5))

# E6: unexposed patients with missing FIS time => assign censoring time
# Status explicitly recorded as unexposed but timing not recorded in source data
# Censoring time assigned as Time_death_days (last known follow-up)
n_e6 <- sum(ds_clean$Status_FIS == 0L & is.na(ds_clean$Time_FIS_days),
            na.rm = TRUE)
ds_clean[Status_FIS == 0L & is.na(Time_FIS_days),
         Time_FIS_days := Time_death_days]
cat(sprintf("E6 - Unexposed, FIS time missing => assign censoring time: %d\n", n_e6))
# Verification
cat(sprintf("\nRaw cohort N:                          %d\n", nrow(ds0)))
cat(sprintf("Excluded (Rule O1):                    %d\n", n_o1))
cat(sprintf("Excluded (Rule E4):                    %d\n", n_e4))
cat(sprintf("Total excluded:                        %d\n", n_o1 + n_e4))
cat(sprintf("Final clean dataset N:                 %d\n", nrow(ds_clean)))
cat("NOTE: Single consistent population used for ALL analyses\n")
cat(sprintf("Time_death_days NAs remaining:         %d (expect 0)\n",
            sum(is.na(ds_clean$Time_death_days))))
cat(sprintf("Status_death NAs remaining:            %d (expect 0)\n",
            sum(is.na(ds_clean$Status_death))))
cat(sprintf("Status_FIS NAs remaining:              %d (expect 0)\n",
            sum(is.na(ds_clean$Status_FIS))))
cat(sprintf("Time_FIS_days NAs remaining:           %d (expect 0)\n",
            sum(is.na(ds_clean$Time_FIS_days))))
cat(sprintf("Exposed (Status_FIS == 1):             %d (%.1f%%)\n",
            sum(ds_clean$Status_FIS == 1L),
            100 * mean(ds_clean$Status_FIS == 1L)))

# Save master clean dataset - input for ALL downstream analyses
saveRDS(ds_clean, OUTFILE)
cat(sprintf("\n✓ Master clean dataset saved: %s\n", OUTFILE))

# -------------------------------------------------------------------------

# 4b) Table S1 — Cohort derivation and exposure classification (gt)

# -------------------------------------------------------------------------



n_raw_exp <- sum(ds_clean$Status_FIS == 1L) + n_e4 + n_e5

n_final_exp <- sum(ds_clean$Status_FIS == 1L)



tbl_s1_data <- data.frame(
  
  Step = c("1", "2", "3", "", "1", "2", "3", ""),
  
  Description = c(
    
    "Raw pooled cohort after harmonisation",
    
    "Excluded: missing or zero follow-up time",
    
    "Missing vital status with valid follow-up time \u2014 treated as censored",
    
    "Final analytic cohort",
    
    "Patients with inappropriate shock recorded in raw data",
    
    "Excluded: confirmed inappropriate shock but exact timing unknown",
    
    "Reclassified as unexposed: shock recorded at or after end of follow-up",
    
    "Final exposed patients in analytic cohort"
    
  ),
  
  `N patients` = c(
    
    format(nrow(ds0),   big.mark = ","),
    
    paste0("\u2212", n_o1),
    
    as.character(n_o2),
    
    format(nrow(ds_clean), big.mark = ","),
    
    as.character(n_raw_exp),
    
    paste0("\u2212", n_e4),
    
    paste0("\u2212", n_e5),
    
    as.character(n_final_exp)
    
  ),
  
  check.names = FALSE,
  
  stringsAsFactors = FALSE
  
)



table_S1 <- tbl_s1_data |>
  
  gt() |>
  
  tab_header(title = "Table S1. Cohort derivation and exposure classification") |>
  
  tab_row_group(label = "OUTCOME VARIABLE (ALL-CAUSE MORTALITY)",          rows = 1:4) |>
  
  tab_row_group(label = "EXPOSURE VARIABLE (FIRST INAPPROPRIATE ICD SHOCK)", rows = 5:8) |>
  
  tab_style(
    
    style     = list(cell_text(weight = "bold"), cell_fill(color = "#F2F2F2")),
    
    locations = cells_body(rows = c(4, 8))
    
  ) |>
  
  tab_style(
    
    style     = list(cell_fill(color = "#D9D9D9"), cell_text(weight = "bold")),
    
    locations = cells_row_groups()
    
  ) |>
  
  cols_align(align = "center", columns = c(Step, `N patients`)) |>
  
  tab_footnote(
    
    footnote = paste0(
      
      "Raw cohort N = ", format(nrow(ds0), big.mark = ","),
      
      " from four pooled registries (EU-CERT-ICD, HELIOS, PROSE-ICD, ISRAEL-ICD). ",
      
      "Total excluded = ", n_o1 + n_e4,
      
      " (missing follow-up time: n = ", n_o1,
      
      "; confirmed shock with unknown timing: n = ", n_e4, "). ",
      
      "Final analytic cohort N = ", format(nrow(ds_clean), big.mark = ","), ". ",
      
      "Of ", n_raw_exp, " patients with an inappropriate shock recorded, ",
      
      n_final_exp, " (", round(n_final_exp / nrow(ds_clean) * 100, 1),
      
      "%) were classified as exposed in all analyses. FIS, first inappropriate shock."
      
    ),
    
    locations = cells_column_labels(columns = Description)
    
  )



gtsave(table_S1, file.path(OUTDIR, "TableS1_CohortDerivation.html"), inline_css = TRUE)

cat(sprintf("\u2713 Table S1 saved: %s\n", file.path(OUTDIR, "TableS1_CohortDerivation.html")))

cat("\u2713 Table S1 gt object available as: table_S1\n")


# -------------------------------------------------------------------------
# 5) Helper: Poisson exact 95% CI for crude incidence rates
# -------------------------------------------------------------------------
poisson_ci_rate <- function(n_events, person_years, multiplier = 100) {
  list(
    rate = n_events / person_years * multiplier,
    lo   = qchisq(0.025, 2 * n_events)       / 2 / person_years * multiplier,
    hi   = qchisq(0.975, 2 * (n_events + 1)) / 2 / person_years * multiplier
  )
}

save_gt <- function(dt, stem) {
  fwrite(dt, file.path(OUTDIR, paste0(stem, ".csv")))
  gtsave(gt(dt), file.path(OUTDIR, paste0(stem, ".html")), inline_css = TRUE)
}

# -------------------------------------------------------------------------
# 6) Create the analysis dataset 
# We use the clean dataset directly since all problematic patients 
# have already been removed or corrected, one single dataset is used 
# for both mortality and inappropriate shock calculations. 
# This ensures both rates come from exactly the same group of patients.
# -------------------------------------------------------------------------
cat("\n--- Defining base dataset for all incidence calculations ---\n")

ds_base <- as.data.table(copy(ds_clean))

cat(sprintf("Base dataset N (mortality + FIS):   %d\n", nrow(ds_base)))
cat(sprintf("Status_death NA remaining:          %d (expect 0)\n",
            ds_base[is.na(Status_death), .N]))
cat(sprintf("Status_FIS NA remaining:            %d (expect 0)\n",
            ds_base[is.na(Status_FIS), .N]))
cat(sprintf("Exposed (Status_FIS == 1):          %d (%.1f%%)\n",
            ds_base[Status_FIS == 1L, .N],
            ds_base[Status_FIS == 1L, .N] / nrow(ds_base) * 100))

# -------------------------------------------------------------------------
# 7) Calculate all-cause mortality incidence rate
# We count how many patients died and divide by total follow-up time
# across all patients (in person-years).
# Result: number of deaths per 100 person-years, with 95% CI.
# -------------------------------------------------------------------------
cat("\n--- All-cause Mortality Incidence Rate ---\n")

D    <- sum(ds_base$Status_death)
PY_d <- sum(ds_base$Time_death_days) / 365.25
ci_d <- poisson_ci_rate(D, PY_d)

death_crude <- data.table(
  outcome      = "All-cause mortality",
  n_total      = nrow(ds_base),
  n_events     = D,
  percent      = round(D / nrow(ds_base) * 100, 1),
  person_years = round(PY_d, 1),
  rate_100py   = round(ci_d$rate, 2),
  ci_lower     = round(ci_d$lo, 2),
  ci_upper     = round(ci_d$hi, 2),
  ci_method    = "Poisson exact"
)

print(death_crude)
save_gt(death_crude, "results_crude_death")

# -------------------------------------------------------------------------
# 8) Calculate inappropriate ICD shock incidence rate
# We count how many patients had a first inappropriate shock and
# divide by the same total follow-up time used for mortality above.
# Using the same denominator for both rates makes them directly
# comparable — both expressed over the same person-time.
# Result: number of first shocks per 100 person-years, with 95% CI.
# -------------------------------------------------------------------------
cat("\n--- Inappropriate ICD Therapy Incidence Rate ---\n")

FIS  <- sum(ds_base$Status_FIS == 1L)
PY_f <- sum(ds_base$Time_death_days) / 365.25   # identical to PY_d

ci_f <- poisson_ci_rate(FIS, PY_f)

inappropriate_therapy_incidence <- data.table(
  outcome        = "Inappropriate ICD therapy (first shock)",
  n_total        = nrow(ds_base),
  n_events       = FIS,
  percent_events = round(FIS / nrow(ds_base) * 100, 1),
  person_years   = round(PY_f, 1),
  crude_rate     = round(ci_f$rate, 2),
  ci_lower       = round(ci_f$lo, 2),
  ci_upper       = round(ci_f$hi, 2),
  ci_method      = "Poisson exact",
  rate_95ci      = sprintf("%.2f (%.2f–%.2f)", ci_f$rate, ci_f$lo, ci_f$hi)
)

print(inappropriate_therapy_incidence)
save_gt(inappropriate_therapy_incidence, "inappropriate_therapy_incidence")

cat(sprintf("\nDenominator check — PY_d == PY_f: %.1f == %.1f (expect identical)\n",
            PY_d, PY_f))

# -------------------------------------------------------------------------
# 9) deaths among exposed vs unexposed patients.
# -------------------------------------------------------------------------
cat("\n--- Deaths by inappropriate shock status (reviewer comment AS4) ---\n")

ctx_tab <- ds_base[!is.na(Status_death),
                   .(n_patients    = .N,
                     deaths        = sum(Status_death == 1L),
                     death_percent = round(mean(Status_death == 1L) * 100, 1)),
                   by = Status_FIS
][, exposure_label := fifelse(Status_FIS == 1L,
                              "Inappropriate shock (ever)",
                              "No inappropriate shock")
][, .(Status_FIS, exposure_label, n_patients, deaths, death_percent)]

print(ctx_tab)
cat(sprintf(
  "\nManuscript line: Of the %d patients experiencing inappropriate shock,\n%d deaths occurred during follow-up (%.1f%%).\n",
  ctx_tab[Status_FIS == 1L, n_patients],
  ctx_tab[Status_FIS == 1L, deaths],
  ctx_tab[Status_FIS == 1L, death_percent]
))

save_gt(ctx_tab, "deaths_by_exposure")

# -------------------------------------------------------------------------
# 10) Inappropriate shock distribution (cohort-based)
# -------------------------------------------------------------------------

exposure_summary <- ds_base[,
                            .(N = .N, Percent = round(.N / nrow(ds_base) * 100, 1)),
                            by = Status_FIS
][, Label := fifelse(Status_FIS == 0L,
                     "No Inappropriate shock",
                     "Inappropriate shock")
][, .(Status_FIS, Label, N, Percent)]

print(exposure_summary)
save_gt(exposure_summary, "exposure_summary")

# -------------------------------------------------------------------------
# 11) POWER ANALYSIS — Schoenfeld approximation
# D_events and p_exposed derived from base dataset
# for consistency with incidence calculations above
# This tells us whether the study had enough statistical power to detect a  
# meaningful effect of inappropriate shock on mortality.
#
# We calculate:
# - Post-hoc power: given our observed data, how likely were we to
#   detect hazard ratios of 1.3, 1.4, and 1.5?
# - MDHR (Minimum Detectable Hazard Ratio): what is the smallest
#   hazard ratio we could reliably detect at 80% power?
# -------------------------------------------------------------------------
cat("\n--- Power Analysis (Schoenfeld approximation for Cox regression) ---\n")

alpha        <- 0.05
target_power <- 0.80
HR_grid      <- c(1.3, 1.4, 1.5)

D_events  <- D                        # total deaths from section 6
p_exposed <- FIS / nrow(ds_base)      # exposure prevalence from section 7

cat(sprintf("Outcome events (D) = %d\n", D_events))
cat(sprintf("Exposure prevalence (p) = %.4f (%.2f%%)\n", p_exposed, p_exposed * 100))

calc_cox_power <- function(n_events, prop_exposed, HR, alpha = 0.05) {
  z   <- qnorm(1 - alpha / 2)
  lam <- log(HR) * sqrt(n_events * prop_exposed * (1 - prop_exposed))
  pnorm(lam - z) + pnorm(-lam - z)
}

calc_mdhr <- function(n_events, prop_exposed, target_power = 0.80, alpha = 0.05) {
  exp((qnorm(1 - alpha / 2) + qnorm(target_power)) /
        sqrt(n_events * prop_exposed * (1 - prop_exposed)))
}

power_tab <- data.table(
  HR            = HR_grid,
  posthoc_power = round(sapply(HR_grid, \(h) calc_cox_power(D_events, p_exposed, h, alpha)), 3)
)[, interpretation := fifelse(posthoc_power >= 0.80, "Adequate", "Inadequate")]

print(power_tab)

mdhr <- calc_mdhr(D_events, p_exposed, target_power, alpha)
cat(sprintf("\nMinimum Detectable HR (MDHR) at %.0f%% power: %.2f\n",
            target_power * 100, mdhr))
cat(sprintf("Power to detect HR = 1.3: %.1f%%\n",
            calc_cox_power(D_events, p_exposed, 1.3, alpha) * 100))

power_summary <- data.table(
  outcome_events_D      = D_events,
  exposure_prevalence_p = round(p_exposed, 4),
  alpha                 = alpha,
  target_power          = target_power,
  MDHR                  = round(mdhr, 2),
  power_for_HR_1.3      = round(calc_cox_power(D_events, p_exposed, 1.3, alpha), 3)
)

power_curve_tab <- data.table(HR = seq(1.1, 2.0, by = 0.1))[
  , Power := round(sapply(HR, \(h) calc_cox_power(D_events, p_exposed, h, alpha)), 3)]

save_gt(power_tab,       "power_posthoc_table")
save_gt(power_summary,   "power_mdhr_summary")
save_gt(power_curve_tab, "power_curve_table")

cat("\n✓ Output files saved to:", OUTDIR, "\n")
cat("  - master_clean_dataset.rds\n")
cat("  - results_crude_death.csv/.html\n")
cat("  - inappropriate_therapy_incidence.csv/.html\n")
cat("  - deaths_by_exposure.csv/.html\n")
cat("  - exposure_summary.csv/.html\n")
cat("  - power_posthoc_table.csv/.html\n")
cat("  - power_curve_table.csv/.html\n")
cat("  - power_mdhr_summary.csv/.html\n")














