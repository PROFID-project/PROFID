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




