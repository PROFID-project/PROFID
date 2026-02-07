###############################################################################
# KM 2: Overall survival by EVER inappropriate ICD therapy (TRAIN imputed data)
###############################################################################

library(mice)
library(data.table)
library(survival)
library(ggplot2)


# Provide your mids objects here (recommended) ------------------------------
TRAIN_RDS <- "T:/Imputed_data/mice_train_object3.rds"  # <-- CHANGE

t3_train <- readRDS(TRAIN_RDS)

# Use one imputation for plotting
d <- as.data.table(complete(t3_train, 1))

# Required variables
# Status_death  : 1 = death, 0 = alive
# Survival_time: time to death or censoring
# Status_FIS   : 1 = ever inappropriate shock, 0 = never

d <- d[
  !is.na(Status_death) &
    !is.na(Survival_time) &
    !is.na(Status_FIS) &
    Survival_time >= 0
]

# Group definition
d[, shock_group := ifelse(
  Status_FIS == 1,
  "Received inappropriate shock",
  "No inappropriate shock"
)]

# KM fit
fit <- survfit(
  Surv(Survival_time, Status_death) ~ shock_group,
  data = d
)

# Plot
plot(
  fit,
  col = c("red", "blue"),
  lwd = 2,
  xlab = "Days since baseline",
  ylab = "Overall survival probability",
  main = "Overall survival by inappropriate shock status"
)

legend(
  "bottomleft",
  legend = levels(factor(d$shock_group)),
  col = c("red", "blue"),
  lwd = 2,
  bty = "n"
)
