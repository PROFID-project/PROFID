# Install once
# install.packages(c("survival", "survminer", "dplyr", "broom"))
# if(!"cmprsk" %in% installed.packages()[,"Package"]) install.packages("cmprsk")

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




