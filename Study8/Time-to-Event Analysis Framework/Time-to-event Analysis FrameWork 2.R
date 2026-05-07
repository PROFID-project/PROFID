###############################################################
# 3.2 Time-to-Event Analysis Framework
###############################################################

install.packages("fastcmprsk")
install.packages("caret")
library(survival)
library(cmprsk)
library(ggplot2)

# =============================================================
# Load data
# =============================================================
datadir <- "T:/PROFID/Study8/Variable Selection & Model Development/Files"
df <- read.csv(file.path(datadir, "vs_data_complete.csv"))

# Convert variables
df$Survival_time <- as.numeric(df$Survival_time)
df$Status <- as.numeric(df$Status)

# Create event indicators
df$event_SCD <- ifelse(df$Status == 1, 1, 0)
df$event_COMP <- ifelse(df$Status == 2, 1, 0)

# Output folder
outdir <- "T:/PROFID/Study8/Time-to-event Analysis Framework 2/files"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

###############################################################
# ===================== COX PH MODEL ==========================
###############################################################

cox_model <- coxph(Surv(Survival_time, event_SCD) ~ ., data = df)

cat("\n================ Cox Model Summary ================\n")
print(summary(cox_model))

write.csv(broom::tidy(cox_model, conf.int = TRUE, exponentiate = TRUE),
          file.path(outdir,"cox_model_results.csv"),
          row.names = FALSE)

###############################################################
# ============ CHECK PROPORTIONAL HAZARDS =====================
###############################################################

cox_zph <- cox.zph(cox_model)

cat("\n================ Schoenfeld PH Test ================\n")
print(cox_zph)

# Save Schoenfeld residual test table
write.csv(as.data.frame(cox_zph$table),
          file.path(outdir, "schoenfeld_test_results.csv"),
          row.names = TRUE)

###############################################################
# ============ SCHOENFELD RESIDUAL PLOTS =====================
###############################################################

# Plot ALL Schoenfeld residual curves in one PDF
pdf(file.path(outdir, "Schoenfeld_Residuals_All.pdf"), width = 10, height = 10)
plot(cox_zph)
dev.off()

# Also save individual PNG files for each covariate
cov_names <- rownames(cox_zph$table)[-nrow(cox_zph$table)] # remove GLOBAL row

for (cov in cov_names) {
  png(file.path(outdir, paste0("Schoenfeld_", cov, ".png")),
      width = 900, height = 700)
  
  plot(cox_zph, var = which(cov_names == cov),
       main = paste("Schoenfeld Residuals:", cov))
  
  dev.off()
}


###############################################################################
# FINE–GRAY FULL DATA (WITH TIME TRACKER)
###############################################################################

library(dplyr)
library(cmprsk)
library(survival)

# ---------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------
BASEDIR <- "T:/PROFID/Study8"
SENS_DIR <- file.path(BASEDIR, "Sensitivity Analysis/Files")
if(!dir.exists(SENS_DIR)) dir.create(SENS_DIR, recursive = TRUE)

# ---------------------------------------------------------------
# LOAD FULL DATA
# ---------------------------------------------------------------
df <- read.csv(file.path(
  "T:/PROFID/Study8/Variable Selection & Model Development/Files",
  "vs_data_complete.csv"
))

cat("Full dataset size:", nrow(df), "rows\n")

# ---------------------------------------------------------------
# KEEP ONLY BIOMARKERS + MANDATORY + ECG + OUTCOME
# ---------------------------------------------------------------
mandatory_vars <- c("Age", "LVEF", "Diabetes", "eGFR")
biomarkers <- c("Haemoglobin", "Cholesterol", "HDL", "LDL", "Triglycerides")
ecg_vars <- c("LBBB", "RBBB", "AF_atrial_flutter")

keep_vars <- unique(c(
  "Survival_time", "Status",
  mandatory_vars, biomarkers, ecg_vars
))

df <- df %>% dplyr::select(all_of(keep_vars))

# ---------------------------------------------------------------
# FACTOR CONVERSION
# ---------------------------------------------------------------
df <- df %>%
  mutate(
    Status = as.numeric(Status),
    across(where(is.character), as.factor),
    across(all_of(c(ecg_vars, "Diabetes")), as.factor)
  )

# ---------------------------------------------------------------
# REMOVE HIGH CORRELATION > 0.80
# ---------------------------------------------------------------
num_df <- df %>% select(where(is.numeric)) %>% select(-Survival_time, -Status)
cor_mat <- cor(num_df, use="pairwise.complete.obs")

high_corr <- which(abs(cor_mat) > 0.80 & abs(cor_mat) < 1, arr.ind = TRUE)
drop_vars <- c()

if (nrow(high_corr) > 0) {
  for (i in seq_len(nrow(high_corr))) {
    v1 <- rownames(cor_mat)[high_corr[i,1]]
    v2 <- colnames(cor_mat)[high_corr[i,2]]
    drop_vars <- c(drop_vars, v2)
  }
  drop_vars <- unique(drop_vars)
  df <- df %>% select(!all_of(drop_vars))
}

cat("Dropped due to high correlation:\n")
print(drop_vars)

write.csv(data.frame(Dropped = drop_vars),
          file.path(SENS_DIR, "Full_Dropped_HighCorr.csv"),
          row.names = FALSE)

# ---------------------------------------------------------------
# PREPARE MATRIX
# ---------------------------------------------------------------
ftime <- df$Survival_time
fstatus <- df$Status
X <- df %>% select(-Survival_time, -Status)

model_matrix <- model.matrix(~ ., data=X)[,-1]

# ---------------------------------------------------------------
# TIME TRACKER FUNCTION
# ---------------------------------------------------------------
progress_timer <- function(step, total_steps, start_time) {
  now <- Sys.time()
  elapsed <- as.numeric(difftime(now, start_time, units="secs"))
  avg_time_per_step <- elapsed / step
  remaining_steps <- total_steps - step
  est_remaining <- remaining_steps * avg_time_per_step
  eta <- now + est_remaining
  
  cat(sprintf("⏳ Step %d/%d | Elapsed: %.1f sec | Remaining: %.1f sec (~%.1f min) | ETA: %s\n",
              step, total_steps, elapsed, est_remaining, est_remaining/60, format(eta, "%H:%M:%S")))
}

# ---------------------------------------------------------------
# RUN FINE–GRAY MODELS
# ---------------------------------------------------------------
total_steps <- 2
start_time <- Sys.time()

cat("\n========================\n")
cat(" Running Fine–Gray FULL DATA\n")
cat("========================\n\n")

### STEP 1 — SCD (cause = 1)
cat("Step 1: SCD model...\n")

fg_scd <- tryCatch(
  {
    crr(
      ftime = ftime,
      fstatus = fstatus,
      cov1 = model_matrix,
      failcode = 1,
      cencode = 0
    )
  },
  error = function(e) e
)

progress_timer(step = 1, total_steps = total_steps, start_time)

### STEP 2 — non-SCD (cause = 2)
cat("\nStep 2: non-SCD model...\n")

fg_nonscd <- tryCatch(
  {
    crr(
      ftime = ftime,
      fstatus = fstatus,
      cov1 = model_matrix,
      failcode = 2,
      cencode = 0
    )
  },
  error = function(e) e
)

progress_timer(step = 2, total_steps = total_steps, start_time)

end_time <- Sys.time()
cat("\nTotal runtime:", difftime(end_time, start_time, units="mins"), "minutes\n")

# ---------------------------------------------------------------
# SAVE RESULTS
# ---------------------------------------------------------------
capture.output(fg_scd, file=file.path(SENS_DIR,"Full_FG_SCD.txt"))
capture.output(fg_nonscd, file=file.path(SENS_DIR,"Full_FG_nonSCD.txt"))

cat("\nResults saved to:", SENS_DIR, "\n")


###############################################################
# 3.2 TIME-TO-EVENT ANALYSIS (PH-Corrected to 120 Days)
###############################################################

# Required packages
packages <- c("survival", "cmprsk", "ggplot2", "broom")
lapply(packages, require, character.only = TRUE)

###############################################################
# ✔ PATHS
###############################################################

BASE <- "T:/PROFID/Study8"
DATADIR <- file.path(BASE, "Variable Selection & Model Development/Files")
OUTDIR <- file.path(BASE, "Time-to-event Analysis Framework 2/files")

if (!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

###############################################################
# ✔ LOAD DATA
###############################################################

df <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

df$Survival_time <- as.numeric(df$Survival_time)
df$Status <- as.numeric(df$Status)

# Create event indicators
df$event_SCD <- ifelse(df$Status == 1, 1, 0)
df$event_COMP <- ifelse(df$Status == 2, 1, 0)

###############################################################
# ✔ STEP 1: Restrict Follow-up to 120 Days (PH-Corrected)
###############################################################

df_120 <- df %>%
  mutate(
    Survival_time = pmin(Survival_time, 120),
    Status_120 = ifelse(Survival_time >= 120, 0, Status),
    event_SCD_120 = ifelse(Status_120 == 1, 1, 0),
    event_COMP_120 = ifelse(Status_120 == 2, 1, 0)
  )

write.csv(df_120, file.path(OUTDIR, "Dataset_Truncated120Days.csv"), row.names = FALSE)

cat("✔ Follow-up truncated to 120 days and dataset saved.\n")

###############################################################
# ✔ STEP 2: Cox Proportional Hazards Model (PH-corrected)
###############################################################

cox_model_120 <- coxph(Surv(Survival_time, event_SCD_120) ~ ., data = df_120)

summary(cox_model_120)

# Save summary table
write.csv(
  broom::tidy(cox_model_120, exponentiate = TRUE, conf.int = TRUE),
  file.path(OUTDIR, "CoxModel_120Days.csv"),
  row.names = FALSE
)

cat("✔ Cox model fitted and saved.\n")

###############################################################
# ✔ STEP 3: Proportional Hazards Test (Schoenfeld)
###############################################################

cox_zph_120 <- cox.zph(cox_model_120)

# Save PH test table
write.csv(
  as.data.frame(cox_zph_120$table),
  file.path(OUTDIR, "Schoenfeld_PH_Test_120Days.csv")
)

# Save ALL Schoenfeld plots in one PDF
pdf(file.path(OUTDIR, "Schoenfeld_All_120Days.pdf"), width = 9, height = 9)
plot(cox_zph_120)
dev.off()

cat("✔ Schoenfeld global & individual plots saved.\n")

# Save individual PNG plots
covs <- rownames(cox_zph_120$table)[-nrow(cox_zph_120$table)]  # remove GLOBAL row

for (v in covs) {
  png(file.path(OUTDIR, paste0("Schoenfeld_", v, "_120Days.png")),
      width = 1000, height = 800)
  plot(cox_zph_120, var = which(covs == v),
       main = paste("Schoenfeld Residuals (120 days):", v))
  dev.off()
}

cat("✔ Individual Schoenfeld plots saved.\n")

###############################################################
# ✔ IDENTIFY WHICH COVARIATES VIOLATE THE PH ASSUMPTION
###############################################################

library(dplyr)
library(splines)
library(broom)

# Convert zph table to dataframe
ph_df <- as.data.frame(cox_zph_120$table)

# Add variable names
ph_df$Variable <- rownames(ph_df)

# Remove GLOBAL row
ph_df <- ph_df[ph_df$Variable != "GLOBAL", ]

# Keep only Variable and p columns (using base R to avoid dplyr conflicts)
ph_df <- ph_df[, c("Variable", "p")]

# Mark violations
ph_df$PH_Violation <- ifelse(ph_df$p < 0.05, "Yes", "No")

# Save report
write.csv(ph_df,
          file.path(OUTDIR, "PH_Violation_Report_120Days.csv"),
          row.names = FALSE)

cat("\nPH violations identified:\n")
print(ph_df)



###############################################################
# ✔ SPLINE-ADJUST ONLY NUMERIC VIOLATING VARIABLES
###############################################################

# Identify numeric variables
numeric_covariates <- names(df_120)[sapply(df_120, is.numeric)]

# Numeric variables violating PH
violating_numeric <- ph_df %>%
  filter(PH_Violation == "Yes", Variable %in% numeric_covariates) %>%
  pull(Variable)

cat("\nNumeric PH violations:\n")
print(violating_numeric)

# If none → exit
if (length(violating_numeric) == 0) {
  
  cat("\nNo numeric PH-violating variables → no spline correction required.\n")
  
} else {
  
  cat("\nFitting spline-adjusted Cox model for:\n")
  print(violating_numeric)
  
  # Only apply spline to continuous numeric variables
  # Remove any 0/1 binary numeric variables
  violating_numeric <- violating_numeric[
    sapply(df_120[violating_numeric], function(x) length(unique(x)) > 5)
  ]
  
  cat("\nSpline will be applied to:\n")
  print(violating_numeric)
  
  spline_terms <- paste0("ns(", violating_numeric, ", df = 3)", collapse = " + ")
  
  # Other (non-violating) covariates
  other_vars <- ph_df %>%
    filter(PH_Violation == "No") %>%
    pull(Variable)
  
  other_terms <- paste(other_vars, collapse = " + ")
  
  # Final formula
  final_formula <- as.formula(
    paste0("Surv(Survival_time, event_SCD_120) ~ ",
           other_terms, " + ", spline_terms)
  )
  
  # Fit model
  spline_model <- coxph(final_formula, data = df_120)
  
  # Save output
  write.csv(
    broom::tidy(spline_model, exponentiate = TRUE, conf.int = TRUE),
    file.path(OUTDIR, "CoxModel_120Days_SplineAdjusted.csv"),
    row.names = FALSE
  )
  
  cat("\n✔ Spline-adjusted model saved: CoxModel_120Days_SplineAdjusted.csv\n")
}


###############################################################
# ✔ STEP 4: Competing-Risk Curves (CIF)
###############################################################

# CIF for SCD & competing death
fit_scd <- survfit(Surv(Survival_time, event_SCD_120) ~ 1, data = df_120)
fit_comp <- survfit(Surv(Survival_time, event_COMP_120) ~ 1, data = df_120)

p1 <- ggsurvplot(fit_scd, conf.int = TRUE, ggtheme = theme_minimal(),
                 title = "CIF – Sudden Cardiac Death (120-Day Truncation)")
ggsave(file.path(OUTDIR, "CIF_SCD_120Days.png"), p1$plot, width = 7, height = 6)

p2 <- ggsurvplot(fit_comp, conf.int = TRUE, ggtheme = theme_minimal(),
                 title = "CIF – Competing Mortality (120-Day Truncation)")
ggsave(file.path(OUTDIR, "CIF_CompetingDeath_120Days.png"), p2$plot, width = 7, height = 6)

cat("✔ CIF plots saved.\n")

###############################################################
# ✔ STEP 5: Save model diagnostics
###############################################################

model_diag <- data.frame(
  Global_PH_pvalue = cox_zph_120$table["GLOBAL", "p"],
  Num_Covariates = length(cox_model_120$coefficients),
  Num_Observations = nrow(df_120)
)

write.csv(model_diag, file.path(OUTDIR, "ModelDiagnostics_120Days.csv"), row.names = FALSE)

cat("\n================ DONE: All files saved ================\n")


###############################################################
# ✔ IDENTIFY WHICH COVARIATES VIOLATED PH ASSUMPTION
###############################################################

# Convert Schoenfeld table to data frame
ph_df <- as.data.frame(cox_zph_120$table)

# Remove the GLOBAL row
ph_df <- ph_df[rownames(ph_df) != "GLOBAL", ]

# Extract p-values
ph_df$Variable <- rownames(ph_df)
ph_df <- ph_df %>% select(Variable, p)

# Mark violations (p < 0.05)
ph_df$PH_Violation <- ifelse(ph_df$p < 0.05, "Yes", "No")

# Count total violations
num_violations <- sum(ph_df$PH_Violation == "Yes")

cat("\n============================\n")
cat("PH Assumption Violations\n")
cat("============================\n")
cat("Total variables tested: ", nrow(ph_df), "\n")
cat("Number violating PH: ", num_violations, "\n\n")

# Print list of violating variables
if(num_violations > 0){
  cat("Variables violating PH assumption:\n")
  print(ph_df %>% filter(PH_Violation == "Yes"))
} else {
  cat("No variables violated PH.\n")
}

# Save results
write.csv(ph_df,
          file.path(OUTDIR, "PH_Violation_Report_120Days.csv"),
          row.names = FALSE)

cat("\nPH violation report saved → PH_Violation_Report_120Days.csv\n")


###############################################################
# ✔ IDENTIFY COVARIATES THAT VIOLATE PH (Schoenfeld test)
###############################################################

library(dplyr)
library(splines)   # needed for spline adjustment

# Convert Schoenfeld object to clean df
ph_df <- as.data.frame(cox_zph_120$table)

# Add rownames as variable names
ph_df$Variable <- rownames(ph_df)

# Remove global test
ph_df <- ph_df %>% filter(Variable != "GLOBAL")

# Extract only variable name + p-value
ph_df_clean <- ph_df %>% select(Variable, p)

# Mark PH violations (p < 0.05)
ph_df_clean$PH_violation <- ifelse(ph_df_clean$p < 0.05, "Yes", "No")

# Count violations
num_violations <- sum(ph_df_clean$PH_violation == "Yes")

cat("\n============================\n")
cat(" PH ASSUMPTION VIOLATIONS\n")
cat("============================\n")
cat("Total variables tested:", nrow(ph_df_clean), "\n")
cat("Variables violating PH:", num_violations, "\n\n")

print(ph_df_clean %>% filter(PH_violation == "Yes"))

# Save PH violation table
write.csv(ph_df_clean,
          file.path(OUTDIR, "PH_Violation_Report_120Days.csv"),
          row.names = FALSE)

cat("\n✔ PH violation report saved.\n")


###############################################################
# ✔ REFIT COX MODEL WITH SPLINES FOR VIOLATING VARIABLES
###############################################################

# Identify violating variables
violators <- ph_df_clean %>% filter(PH_violation == "Yes") %>% pull(Variable)

cat("\nVariables requiring spline correction:\n")
print(violators)

# Example: if LVEF violates PH, include spline term
# We apply **3-knot restricted cubic splines** (rcs)

if ("LVEF" %in% violators) {
    cox_spline <- coxph(Surv(Survival_time, event_SCD_120) ~ 
                          ns(LVEF, df = 3) + .,
                        data = df_120)

    # Save spline-adjusted model
    write.csv(
      broom::tidy(cox_spline, exponentiate = TRUE, conf.int = TRUE),
      file.path(OUTDIR, "CoxModel_120Days_SPLINE_Adjusted.csv"),
      row.names = FALSE
    )

    cat("\n✔ Cox model refitted with spline adjustment for LVEF\n")
} else {
    cat("\nNo spline correction needed — no continuous variables violated PH.\n")
}


###############################################################
# 1. SELECT ONLY VALID BASELINE VARIABLES
###############################################################

baseline_vars <- c(
  "LVEF",
  "Sex",
  "Age",
  "BMI",
  "eGFR",
  "LBBB",
  "RBBB",
  "HF",
  "Diabetes",
  "Hypertension",
  "Smoking",
  "Cholesterol",
  "HDL",
  "LDL",
  "Triglycerides",
  "ACE_inhibitor_ARB",
  "Beta_blockers",
  "Diuretics",
  "Haemoglobin",
  "Group",
  "MI_history",
  "MI_type",
  "PCI",
  "CABG"
)
# Adjust to your dataset columns

df_tv <- df_120[, c("Survival_time", "event_SCD_120", baseline_vars)]

# Interaction term
df_tv$log_time <- log(df_tv$Survival_time + 1)

###############################################################
# 2. FIT CLEAN COX MODEL
###############################################################

cox_tv <- coxph(
  as.formula(
    paste(
      "Surv(Survival_time, event_SCD_120) ~ LVEF * log_time +",
      paste(setdiff(baseline_vars, "LVEF"), collapse = " + ")
    )
  ),
  data = df_tv
)

summary(cox_tv)

###############################################################
# 3. EXTRACT COEFFICIENTS
###############################################################
beta1 <- coef(cox_tv)["LVEF"]
beta2 <- coef(cox_tv)["LVEF:log_time"]

if (is.na(beta1) | is.na(beta2)) {
  stop("❗ LVEF coefficients are NA. The model still cannot estimate them.")
}

###############################################################
# 4. COMPUTE TIME-VARYING EFFECT
###############################################################

times <- 1:120
tv_df <- data.frame(
  Time = times,
  LVEF_Effect = beta1 + beta2 * log(times)
)

###############################################################
# 5. PLOT
###############################################################
library(ggplot2)

p <- ggplot(tv_df, aes(x = Time, y = LVEF_Effect)) +
  geom_line(linewidth = 1.1, color = "blue") +
  geom_hline(yintercept = beta1, linetype = "dashed", color = "red") +
  labs(
    title = "Time-Varying Effect of LVEF",
    x = "Follow-up Time (days)",
    y = "Coefficient for LVEF"
  ) +
  theme_minimal()

ggsave(
  filename = file.path(OUTDIR, "LVEF_TimeVarying_Coefficient.png"),
  plot = p,
  width = 8, height = 6
)

cat("\n✔ Non-empty PNG saved.\n")




