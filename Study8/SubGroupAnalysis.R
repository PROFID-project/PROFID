###############################################################################
# SUBGROUP ANALYSIS — GBM SURVIVAL MODEL
###############################################################################

rm(list = ls())
library(tidyverse)
library(survival)
library(Hmisc)

# ---------------------------------------------------------------------------
# PATHS
# ---------------------------------------------------------------------------
BASEDIR <- "T:/PROFID/Study8"
DATADIR <- file.path(BASEDIR, "Variable Selection & Model Development/Files")
MODELDIR <- file.path(BASEDIR, "Model Validation and Performance/Files")
OUTDIR <- file.path(BASEDIR, "SubGroup Analysis/Files")

if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# ---------------------------------------------------------------------------
# LOAD DATA + FINAL MODEL
# ---------------------------------------------------------------------------
df <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    event_flag = ifelse(Status == 0, 0, 1),
    AF_atrial_flutter = ifelse(AF_atrial_flutter == "Yes", 1, 0)
  )

final_gbm <- readRDS(file.path(MODELDIR, "Final_GBM_Model.rds"))

tune_res <- read.csv(file.path(MODELDIR, "GBM_Tuning_Results_Random.csv"))
best_row <- tune_res[which.max(tune_res$C_index),]

best_params <- list(
  n.trees = best_row$n.trees,
  interaction.depth = best_row$interaction.depth,
  shrinkage = best_row$shrinkage,
  n.minobsinnode = best_row$n.minobsinnode,
  bag.fraction = best_row$bag.fraction
)

###############################################################################
# FUNCTION TO RUN SUBGROUP
###############################################################################

run_subgroup <- function(data, subgroup_name, subgroup_value){
  
  cat("\n=============================\n")
  cat(" Subgroup:", subgroup_name, "=", subgroup_value, "\n")
  cat("=============================\n")
  
  dat <- data %>% filter(.data[[subgroup_name]] == subgroup_value)
  
  if(sum(dat$event_flag) < 5){
    cat(" ⚠ Too few events. Included but marked as UNRELIABLE.\n")
  }
  
  # Predict linear predictor
  lp <- predict(final_gbm,
                newdata = dat,
                n.trees = best_params$n.trees,
                type = "link")
  
  # C-index
  cin <- rcorr.cens(-lp, Surv(dat$Survival_time, dat$event_flag))["C Index"]
  
  # Calibration slope
  cal_fit <- coxph(Surv(Survival_time, event_flag) ~ lp, data = dat)
  cal_slope <- as.numeric(cal_fit$coef[1])
  
  # IBS quick approximation
  risks <- 1 - exp(-exp(lp))
  brier <- mean((dat$event_flag - risks)^2)
  ibs <- sqrt(brier)
  
  out <- data.frame(
    Subgroup = subgroup_name,
    Level = as.character(subgroup_value),
    N = nrow(dat),
    Events = sum(dat$event_flag),
    C_index = round(cin,4),
    Calibration_Slope = round(cal_slope,4),
    IBS = round(ibs,4)
  )
  
  write.csv(out,
            file.path(OUTDIR, paste0("Subgroup_", subgroup_name, "_", subgroup_value, ".csv")),
            row.names = FALSE)
  
  return(out)
}

###############################################################################
# RUN SUBGROUPS
###############################################################################

results <- list()

# 1. ICD subgroup (3 groups)
group_levels <- c("ICD", "NonICD_preserved", "NonICD_reduced")

results[["Group"]] <- bind_rows(
  lapply(group_levels, function(g) run_subgroup(df, "Group", g))
)

# 2. AF vs non-AF
results[["AF"]] <- bind_rows(
  run_subgroup(df, "AF_atrial_flutter", 1),
  run_subgroup(df, "AF_atrial_flutter", 0)
)

# 3. Regions
regions <- sort(unique(df$CVD_risk_region))

results[["Region"]] <- bind_rows(
  lapply(regions, function(r) run_subgroup(df, "CVD_risk_region", r))
)

# 4. Baseline type (MI vs MI40d)
periods <- unique(df$Baseline_type)

results[["Baseline_type"]] <- bind_rows(
  lapply(periods, function(t) run_subgroup(df, "Baseline_type", t))
)

###############################################################################
# SAVE COMBINED RESULTS
###############################################################################

summary_df <- bind_rows(results)

write.csv(summary_df,
          file.path(OUTDIR, "Subgroup_Analysis_AllResults.csv"),
          row.names = FALSE)

cat("\n\n ALL SUBGROUP RESULTS SAVED TO: ", OUTDIR, "\n")