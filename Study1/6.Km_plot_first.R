###############################################################################

# KM 1: Time to first inappropriate ICD intervention (TRAIN imputed data)

###############################################################################



library(mice)

library(data.table)

library(survival)

library(ggplot2)

OUTDIR <- "T:/Study_1/Descriptive_analysis"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

t1_train <- readRDS(TRAIN_RDS)
t1_test  <- readRDS(TEST_RDS)

# Provide your mids objects here (recommended) ------------------------------
TRAIN_RDS <- "T:/Imputed_data/mice_train_object3.rds"  # <-- CHANGE

t2_train <- readRDS(TRAIN_RDS)
t1_test  <- readRDS(TEST_RDS)

# Use one imputation for plotting (standard practice)

d <- as.data.table(complete(t2_train, 1))



# Required variables

# Status_FIS      : 1 = inappropriate shock, 0 = none

# Time_FIS_days   : time to inappropriate shock or censoring



d <- d[
  
  !is.na(Status_FIS) &
    
    !is.na(Time_FIS_days) &
    
    Time_FIS_days >= 0
  
]



# KM fit

fit <- survfit(
  
  Surv(Time_FIS_days, Status_FIS) ~ 1,
  
  data = d
  
)



# Convert to data frame for ggplot

s <- summary(fit)

km_df <- data.table(
  
  time  = s$time,
  
  surv  = s$surv,
  
  lower = s$lower,
  
  upper = s$upper
  
)



# Plot

ggplot(km_df, aes(x = time, y = surv)) +
  
  geom_step(linewidth = 1) +
  
  geom_step(aes(y = lower), linetype = 2) +
  
  geom_step(aes(y = upper), linetype = 2) +
  
  labs(
    
    title = "Time to first inappropriate ICD intervention",
    
    x = "Days since baseline",
    
    y = "Probability of remaining free of inappropriate intervention"
    
  ) +
  
  theme_minimal()

