
###############################################################################
# AUC TIME-DEPENDENT DISCRIMINATION PLOT
###############################################################################

library(ggplot2)
library(readr)
library(dplyr)

AUC_FILE <- "T:/PROFID/Study8/Risk Score Development/Files/AUC_TimeDependent.csv"
OUTDIR <- "T:/PROFID/Study8/Plots"

# Load AUC results
auc_df <- read_csv(AUC_FILE)

# Ensure time is numeric
auc_df <- auc_df %>% mutate(Time = as.numeric(Time))

p_auc <- ggplot(auc_df, aes(x = Time, y = AUC)) +
  geom_line(size = 1.3, colour = "#0072B2") +
  geom_point(size = 2.5, colour = "#D55E00") +
  theme_minimal(base_size = 16) +
  labs(
    title = "Time-Dependent AUC for SCD Prediction",
    x = "Time (days)",
    y = "AUC"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(file.path(OUTDIR, "AUC_TimeDependent.png"), p_auc, width = 8, height = 5, dpi = 300)


###############################################################################
# SUBGROUP ANALYSIS BAR PLOT (C-INDEX)
###############################################################################

library(ggplot2)
library(readr)
library(dplyr)

SUB_FILE <- "T:/PROFID/Study8/SubGroup Analysis/Files/Subgroup_Analysis_AllResults.csv"

dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)

sub <- read_csv(SUB_FILE)

# Create a label combining subgroup + level
sub <- sub %>% mutate(Label = paste(Subgroup, Level, sep = " - "))

p_sub <- ggplot(sub, aes(x = reorder(Label, C_index), y = C_index)) +
  geom_col(fill = "#0072B2") +
  coord_flip() +
  theme_minimal(base_size = 15) +
  labs(
    title = "Subgroup C-index Comparison",
    x = "Subgroup",
    y = "C-index"
  ) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title = element_text(face = "bold")
  )

ggsave(file.path(OUTDIR, "Subgroup_Cindex.png"),
       p_sub, width = 8, height = 6, dpi = 300)


############################################################
# RISK SCORE PLOTTING SCRIPT
# PROFID – Study 8
# Creates all manuscript-quality plots for the RiskScore
############################################################

library(dplyr)
library(ggplot2)
library(survival)
library(survminer)
library(rms)

############################################################
# PATHS
############################################################

BASE <- "T:/PROFID/Study8"
DATADIR <- file.path(BASE, "Variable Selection & Model Development/Files")
MODELDIR <- file.path(BASE, "Model Validation and Performance/Files")
OUTDIR <- file.path(BASE, "Risk Score Development/Files")
PLOTDIR <- file.path(OUTDIR, "Plots")

if (!dir.exists(PLOTDIR)) dir.create(PLOTDIR, recursive = TRUE)

############################################################
# LOAD DATA + FINAL MODEL
############################################################

df <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

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

df <- df %>%
  mutate(
    Survival_time = as.numeric(Survival_time),
    Status = as.numeric(Status),
    event_flag = ifelse(Status == 1, 1, 0)
  )

############################################################
# COMPUTE LP + RISKSCORE
############################################################

df$LP <- predict(
  final_gbm,
  newdata = df,
  n.trees = best_params$n.trees,
  type = "link"
)

scale_factor <- 5
df$RiskScore <- round((df$LP - min(df$LP)) * scale_factor)

df <- df %>%
  mutate(
    RiskGroup = case_when(
      RiskScore <= quantile(RiskScore, 0.33) ~ "Low",
      RiskScore <= quantile(RiskScore, 0.66) ~ "Intermediate",
      TRUE ~ "High"
    )
  )

############################################################
# 1. HISTOGRAM OF RISKSCORE
############################################################

png(file.path(PLOTDIR, "RiskScore_Histogram.png"), width=1000, height=800)
ggplot(df, aes(x = RiskScore)) +
  geom_histogram(binwidth = 1, fill = "#2E86AB", color = "white") +
  theme_minimal(base_size = 18) +
  labs(title = "Distribution of RiskScore",
       x = "RiskScore", y = "Count")
dev.off()

############################################################
# 2. BOXPLOT OF RISKSCORE BY EVENT STATUS
############################################################

png(file.path(PLOTDIR, "RiskScore_by_Event.png"), width=1000, height=800)
ggplot(df, aes(x = as.factor(event_flag), y = RiskScore, fill = as.factor(event_flag))) +
  geom_boxplot() +
  scale_fill_manual(values = c("#6ab04c", "#eb4d4b")) +
  theme_minimal(base_size = 18) +
  labs(title = "RiskScore by Event Status",
       x = "Event (0 = No, 1 = Yes)", y = "RiskScore")
dev.off()

############################################################
# 3. KAPLAN–MEIER CURVES BY RISKGROUP
############################################################

km_fit <- survfit(Surv(Survival_time, event_flag) ~ RiskGroup, data=df)

png(file.path(PLOTDIR, "KM_RiskGroups.png"), width=1200, height=900)
ggsurvplot(
  km_fit, data=df, risk.table=TRUE,
  palette = c("#2ecc71", "#f1c40f", "#e74c3c"),
  title = "Kaplan–Meier Curves by Risk Group",
  legend.title = "Risk Group"
)
dev.off()

############################################################
# 4. RISKSCORE vs PREDICTED SURVIVAL (150 DAYS)
############################################################

dd <- datadist(df); options(datadist="dd")

cox_nom <- cph(
  Surv(Survival_time, event_flag) ~ LP,
  data = df,
  x = TRUE, y = TRUE, surv = TRUE
)

pred_time <- 150
lp_grid <- seq(min(df$LP), max(df$LP), length.out = 200)

surv_values <- sapply(lp_grid, function(lp){
  s <- survest(cox_nom, newdata = data.frame(LP = lp), times = pred_time)
  s$surv
})

plot_df <- data.frame(LP = lp_grid, Survival150 = surv_values)

png(file.path(PLOTDIR, "RiskScore_vs_Survival150.png"), width=1000, height=800)
ggplot(plot_df, aes(x = LP, y = Survival150)) +
  geom_line(size=1.4, color="#1f78b4") +
  theme_minimal(base_size = 18) +
  labs(title = "Predicted Survival Probability at 150 Days",
       x = "Linear Predictor (LP)",
       y = "Predicted Survival at 150 Days")
dev.off()

############################################################
# 5. EVENT RATE BY RISKSCORE
############################################################

event_df <- df %>%
  group_by(RiskScore) %>%
  summarise(EventRate = mean(event_flag))

png(file.path(PLOTDIR, "EventRate_by_RiskScore.png"), width=1000, height=800)
ggplot(event_df, aes(x = RiskScore, y = EventRate)) +
  geom_line(size=1.4, color="#d35400") +
  geom_point(size=2) +
  theme_minimal(base_size = 18) +
  labs(title = "Observed Event Rate by RiskScore",
       x = "RiskScore", y = "Event Rate")
dev.off()

############################################################
# 8. BAR PLOT: NUMBER OF PATIENTS IN EACH RISK GROUP
############################################################

risk_counts <- df %>%
  group_by(RiskGroup) %>%
  summarise(N = n()) %>%
  mutate(RiskGroup = factor(RiskGroup, levels = c("Low","Intermediate","High")))

png(file.path(PLOTDIR, "RiskGroup_Counts.png"), width=1000, height=800)
ggplot(risk_counts, aes(x = RiskGroup, y = N, fill = RiskGroup)) +
  geom_col(width = 0.7) +
  scale_fill_manual(values = c("Low" = "#2ecc71",
                               "Intermediate" = "#f1c40f",
                               "High" = "#e74c3c")) +
  theme_minimal(base_size = 18) +
  labs(title = "Number of Patients in Each Risk Group",
       x = "Risk Group",
       y = "Number of Patients") +
  geom_text(aes(label = N), vjust = -0.8, size = 6)
dev.off()

############################################################
# GROUPED BAR PLOT: Risk Group × ICD / EF Category
############################################################

library(dplyr)
library(ggplot2)

# Paths
BASE <- "T:/PROFID/Study8"
OUTDIR <- file.path(BASE, "Risk Score Development/Files")
PLOTDIR <- file.path(OUTDIR, "Plots")

if (!dir.exists(PLOTDIR)) dir.create(PLOTDIR, recursive = TRUE)

# Load dataset
df <- read.csv(file.path(OUTDIR, "RiskScore_with_Groups.csv"))

# Create summary table (counts)
count_df <- df %>%
  group_by(RiskGroup, Group) %>%
  summarise(N = n(), .groups = "drop")

# Order factors for nicer plotting
count_df$RiskGroup <- factor(count_df$RiskGroup, 
                             levels = c("Low", "Intermediate", "High"))

count_df$Group <- factor(count_df$Group,
                         levels = c("ICD", "NonICD_preserved", "NonICD_reduced"),
                         labels = c("ICD", "Non-ICD Preserved EF", "Non-ICD Reduced EF"))

# Plot
p <- ggplot(count_df, aes(x = RiskGroup, y = N, fill = Group)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.75) +
  theme_minimal(base_size = 18) +
  scale_fill_manual(values = c("#e74c3c", "#3498db", "#2ecc71")) +
  labs(
    title = "Distribution of ICD / EF Category Across Risk Groups",
    x = "Risk Group",
    y = "Number of Patients",
    fill = "Category"
  )

# Save plot
ggsave(
  filename = file.path(PLOTDIR, "RiskGroup_by_Category.png"),
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)
############################################################
# DESCRIPTIVE ANALYSIS TABLE
############################################################

library(dplyr)
library(readr)

# Paths
BASE <- "T:/PROFID/Study8"
DATADIR <- file.path(BASE, "Variable Selection & Model Development/Files")
OUTDIR  <- file.path(BASE, "T:\PROFID\Study8\Plots")

if(!dir.exists(OUTDIR)) dir.create(OUTDIR, recursive = TRUE)

# Load dataset
df <- read.csv(file.path(DATADIR, "vs_data_complete.csv"))

# Prepare variables
df <- df %>%
  mutate(
    event_flag = ifelse(Status == 1, 1, 0),
    AF_flag    = ifelse(AF_atrial_flutter == "Yes", 1, 0)
  )

############################################################
# TABLE 1: GROUP-LEVEL SUMMARY
############################################################

table1 <- df %>%
  group_by(Group) %>%
  summarise(
    N = n(),
    SCD_events = sum(event_flag),
    AF_cases = sum(AF_flag),
    SCD_in_AF = sum(event_flag == 1 & AF_flag == 1),
    
    # Percentages
    Pct_SCD = round(100 * SCD_events / N, 1),
    Pct_AF = round(100 * AF_cases / N, 1),
    Pct_SCD_in_AF = round(100 * SCD_in_AF / AF_cases, 1)
  )

# Save table
write.csv(table1,
          file.path(OUTDIR, "Descriptive_Table_GroupSummary.csv"),
          row.names = FALSE)

cat("Saved: Descriptive_Table_GroupSummary.csv\n")

############################################################
# TABLE 2: OVERALL COUNTS
############################################################

overall <- data.frame(
  Total_N = nrow(df),
  Total_SCD = sum(df$event_flag),
  Total_AF = sum(df$AF_flag),
  SCD_in_AF = sum(df$event_flag == 1 & df$AF_flag == 1)
)

write.csv(overall,
          file.path(OUTDIR, "Descriptive_Table_OverallCounts.csv"),
          row.names = FALSE)

cat("Saved: Descriptive_Table_OverallCounts.csv\n")

############################################################
# OPTIONAL: Baseline numeric summaries (age, LVEF, BMI)
############################################################

baseline_table <- df %>%
  group_by(Group) %>%
  summarise(
    Mean_Age = round(mean(Age, na.rm=TRUE),1),
    Mean_LVEF = round(mean(LVEF, na.rm=TRUE),1),
    Mean_BMI = round(mean(BMI, na.rm=TRUE),1)
  )

write.csv(baseline_table,
          file.path(OUTDIR, "Descriptive_Table_BaselineFeatures.csv"),
          row.names = FALSE)

cat("Saved: Descriptive_Table_BaselineFeatures.csv\n")


############################################################
# FIXED BAR CHART:
# SCD events WITH vs WITHOUT AF across Risk Groups
############################################################

library(dplyr)
library(ggplot2)

# Load data
df <- read.csv("T:/PROFID/Study8/Risk Score Development/Files/RiskScore_with_Groups.csv")

# AF variable cleanup
df$AF_flag <- ifelse(df$AF_atrial_flutter == "Yes", "AF Present", "AF Absent")

# Summarize SCD counts by RiskGroup × AF status
plot_df <- df %>%
  filter(event_flag == 1) %>%
  group_by(RiskGroup, AF_flag) %>%
  summarise(SCD_events = n(), .groups = "drop")

# Make sure order is correct
plot_df$RiskGroup <- factor(plot_df$RiskGroup,
                            levels = c("Low", "Intermediate", "High"))
plot_df$AF_flag <- factor(plot_df$AF_flag,
                          levels = c("AF Absent", "AF Present"))

# Plot
p <- ggplot(plot_df, aes(x = RiskGroup, y = SCD_events, fill = AF_flag)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = SCD_events),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 5) +
  scale_fill_manual(values = c("#2980B9", "#C0392B")) +
  theme_minimal(base_size = 18) +
  labs(
    title = "SCD Events With vs Without AF Across Risk Groups",
    x = "Risk Group",
    y = "Number of SCD Events",
    fill = "AF Status"
  )

# Save plot
ggsave(
  filename = "T:/PROFID/Study8/Risk Score Development/Files/Plots/SCD_by_AF_and_RiskGroup_FIXED.png",
  plot = p,
  width = 12,
  height = 8,
  dpi = 300
)
