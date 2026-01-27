# Install once
install.packages(c("survival", "survminer", "dplyr", "broom"))
if(!"cmprsk" %in% installed.packages()[,"Package"]) install.packages("cmprsk")

library(survival)
library(survminer)
library(dplyr)
library(broom)
library(cmprsk)
# combining the all of them together

file1 <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/ICD.csv")
file2 <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_preserved.csv")
file3 <- read.csv("S:/AG/f-dhzc-profid/Data Transfer to Charite/NonICD_reduced.csv")


nrow(file1); nrow(file2); nrow(file3)

file1$Group <- "ICD"
file2$Group <- "NonICD_preserved"
file3$Group <- "NonICD_reduced"

combined <- bind_rows(file1, file2, file3)

# removing missing event data
combined <- combined %>%
  filter(!is.na(Status)) %>%       # remove missing event data
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.integer(Status),
    Sex = factor(Sex),
    Group = factor(Group)
  )

# Check that grouping variable has 3 levels
table(combined$Group)


# Create the survival object for SCD

surv_obj <- with(combined, Surv(time = Survival_time, event = (Status == 1)))
head(surv_obj)
# Cox Model

cox_model <- coxph(Surv(Survival_time, Status == 1) ~ Age + Sex + LVEF + Group, data = combined)
summary(cox_model)

# Save Cox Model Results

# Convert the Cox model summary into a tidy data frame
cox_summary <- broom::tidy(cox_model, exponentiate = TRUE, conf.int = TRUE)

# Save it to your Study 8 folder
write.csv(cox_summary,
          "T:/PROFID/Study8/Time-to-Event Analysis Framework/cox_model_results.csv",
          row.names = FALSE)


View(cox_summary)


# 1Check proportional hazards assumption using Schoenfeld residuals
cox_zph <- cox.zph(cox_model)
print(cox_zph)

# 2️ Plot Schoenfeld residuals for visual inspection
# (each variable should have a roughly horizontal line if PH assumption holds)
plot(cox_zph)


# Save proportional hazards assumption test results correctly
ph_results <- as.data.frame(cox_zph$table)  # Extracts the table part

# Save to your study folder
write.csv(ph_results,
          "T:/PROFID/Study8/Time-to-Event Analysis Framework/Files/PH_Assumption_Results.csv",
          row.names = TRUE)

# Save Schoenfeld residuals plot
# Create separate PNGs for each variable's Schoenfeld residual plot

var_names <- rownames(cox_zph$table)  # Extract variable names (Age, Sex, etc.)

for (i in seq_along(var_names)) {
  var <- var_names[i]
  
  # Skip the GLOBAL test (it has no plot)
  if (var == "GLOBAL") next
  
  print(paste("Saving plot for:", var))
  
  # Define the output file path
  file_path <- paste0("T:/PROFID/Study8/Time-to-Event Analysis Framework/Files/Schoenfeld_", var, ".png")
  
  # Save the plot
  png(file_path, width = 800, height = 600)
  plot(cox_zph[i], main = paste("Schoenfeld Residuals for", var))
  dev.off()
}

# Save log–log survival plot

png("T:/PROFID/Study8/Time-to-Event Analysis Framework/loglog_survival_plot.png", width = 800, height = 600)
ggsurvplot(
  
  survfit(Surv(Survival_time, Status == 1) ~ Group, data = combined),
  fun = "cloglog",
  palette = "Dark2",
  title = "Log–Log Survival Curves by Group",
  xlab = "Log(Time)",
  ylab = "Log(-Log(Survival Probability))"
)

dev.off()

# Time-varying coefficient model

cox_model_tv <- coxph(Surv(Survival_time, Status == 1) ~ Age * Survival_time + Sex * Survival_time + LVEF * Survival_time + Group * Survival_time, data = combined)



# Summarize the model

summary(cox_model_tv)


# Stratified Cox model

cox_model_strat <- coxph(Surv(Survival_time, Status == 1) ~ Age + Sex + LVEF + strata(Group), data = combined)



# Summarize the model

summary(cox_model_strat)


# --- SAVE TIME-VARYING COEFFICIENT MODEL RESULTS ---

# Convert the summary to a tidy table
tv_summary <- broom::tidy(cox_model_tv, exponentiate = TRUE, conf.int = TRUE)

# Save to CSV
write.csv(tv_summary,
          "T:/PROFID/Study8/Time-to-Event Analysis Framework/Files/Cox_TimeVarying_Model_Results.csv",
          row.names = FALSE)


# --- SAVE STRATIFIED COX MODEL RESULTS ---

# Convert to tidy data frame
strat_summary <- broom::tidy(cox_model_strat, exponentiate = TRUE, conf.int = TRUE)

# Save to CSV
write.csv(strat_summary,
          "T:/PROFID/Study8/Time-to-Event Analysis Framework/Files/Cox_Stratified_Model_Results.csv",
          row.names = FALSE)


# Competing risks analysis (Fine-Gray subdistribution hazards for non-SCD mortality)

table(combined$Status, useNA = "ifany")


# Clean data for Fine–Gray model
fg_data <- combined %>%
  select(Survival_time, Status, Age, Sex, LVEF, Group) %>%
  filter(complete.cases(.))

# Check how many rows remain
nrow(fg_data)

# Fine–Gray Competing Risks Model
library(cmprsk)



# Step 2: Create a 10% random sample
set.seed(42)  # for reproducibility
fg_data_10 <- fg_data %>% sample_frac(0.10)
nrow(fg_data_10)


# Step 3: Define a reusable model-running function
run_fg_models <- function(data, label) {
  message("Running Fine–Gray models for: ", label)
  
  # --- Fine–Gray for SCD (failcode = 1)
  fg_scd <- crr(
    ftime = data$Survival_time,
    fstatus = data$Status,
    cov1 = model.matrix(~ Age + Sex + LVEF + Group, data = data)[, -1],
    failcode = 1,   # SCD
    cencode = 0,
    variance = TRUE
  )
  
  # --- Fine–Gray for non-SCD (failcode = 2)
  fg_nonscd <- crr(
    ftime = data$Survival_time,
    fstatus = data$Status,
    cov1 = model.matrix(~ Age + Sex + LVEF + Group, data = data)[, -1],
    failcode = 2,   # Non-SCD
    cencode = 0,
    variance = TRUE
  )
  
  # Step 4: Save model summaries as text files
  sink(paste0("T:/PROFID/Study8/Time-to-Event Analysis Framework/Files/FineGray_", label, "_SCD_Summary.txt"))
  cat("Fine–Gray Model for Sudden Cardiac Death (SCD, failcode = 1)\n\n")
  print(summary(fg_scd))
  sink()
  
  sink(paste0("T:/PROFID/Study8/Time-to-Event Analysis Framework/Files/FineGray_", label, "_NonSCD_Summary.txt"))
  cat("Fine–Gray Model for Non-Sudden Cardiac Death (Non-SCD, failcode = 2)\n\n")
  print(summary(fg_nonscd))
  sink()
  
  # Step 5: Create tidy data frames for easy comparison
  tidy_fg <- function(model) {
    data.frame(
      Variable = names(model$coef),
      Estimate = model$coef,
      sHR = exp(model$coef),
      SE = sqrt(diag(model$var)),
      Z = model$coef / sqrt(diag(model$var)),
      P_value = 1 - pchisq((model$coef / sqrt(diag(model$var)))^2, df = 1)
    )
  }
  
  fg_scd_tbl <- tidy_fg(fg_scd)
  fg_nonscd_tbl <- tidy_fg(fg_nonscd)
  
  # Step 6: Save results to CSV
  write.csv(fg_scd_tbl,
            paste0("T:/PROFID/Study8/Time-to-Event Analysis Framework/Files/FineGray_", label, "_SCD.csv"),
            row.names = FALSE)
  
  write.csv(fg_nonscd_tbl,
            paste0("T:/PROFID/Study8/Time-to-Event Analysis Framework/Files/FineGray_", label, "_NonSCD.csv"),
            row.names = FALSE)
  
  message("Models completed and saved for: ", label)
}


# Step 7: Run on 10% subset (fast check)
run_fg_models(fg_data_10, "10pct")

# step 8: Run on full dataset (final model — may take longer)
run_fg_models(fg_data, "100pct")