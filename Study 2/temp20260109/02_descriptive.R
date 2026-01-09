library(dplyr)
library(tidyr)
library(tibble)
library(lubridate)
library(ggplot2)

setwd("T:/PROFID/data/processed")
# Load the dataset back into R
df <- readRDS("df_cleaned.rds")

# ------------------------------
# Define variables
# ------------------------------
cont_vars <- c("Age", "BMI", "SBP", "DBP", "LVEF_std",
               "LVDD", "Cholesterol", "HDL", "LDL", "Triglycerides",
               "Sodium", "Potassium", "BUN", "HbA1c", "Haemoglobin", "CRP")

cat_vars <- c("Sex", "NYHA", "MI_type", "MI_history", "HF", "Stroke_TIA",
              "Diabetes", "Hypertension", "Smoking", "Alcohol",
              "FH_CAD", "FH_SCD")

# ------------------------------
# Continuous summary
# ------------------------------
cont_summary <- df %>%
  select(all_of(c("season", cont_vars))) %>%
  pivot_longer(-season, names_to = "Variable", values_to = "Value") %>%
  group_by(season, Variable) %>%
  summarise(
    n    = sum(!is.na(Value)),                  # number of non-missing values
    mean = mean(Value, na.rm = TRUE),
    sd   = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Value = as.character(sprintf("%.1f ± %.1f (n=%d)", mean, sd, n))
  ) %>%
  select(Variable, season, Value) %>%
  pivot_wider(names_from = season, values_from = Value)

# Ensure character type across the board
cont_summary[] <- lapply(cont_summary, as.character)

# ------------------------------
# Categorical summary (expand levels)
# ------------------------------
cat_summary <- df %>%
  select(all_of(c("season", cat_vars))) %>%
  pivot_longer(-season, names_to = "Variable", values_to = "Level") %>%
  filter(!is.na(Level)) %>%
  group_by(season, Variable, Level) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(season, Variable) %>%
  mutate(
    pct   = round(100 * N / sum(N), 1),
    Value = paste0(N, " (", pct, "%)")
  ) %>%
  ungroup() %>%
  mutate(Variable = paste0(Variable, "=", Level)) %>%
  select(-Level, -N, -pct) %>%
  pivot_wider(
    names_from  = season,
    values_from = Value,
    values_fill = ""        # make sure empty cells are ""
  )

# force all to character to avoid list/char mismatch
cat_summary <- mutate_all(cat_summary, as.character)



# ------------------------------
# Combine continuous + categorical
# ------------------------------
baseline_table <- bind_rows(cont_summary, cat_summary)

# Save
write.csv(baseline_table, "~/output/tables/baseline_characteristics_by_season.csv", row.names = FALSE)


cont_summary <- df %>%
  select(all_of(c("season", cont_vars))) %>%
  pivot_longer(-season, names_to = "Variable", values_to = "Value") %>%
  group_by(season, Variable) %>%
  summarise(
    Value = sprintf("%.1f ± %.1f", mean(Value, na.rm = TRUE), sd(Value, na.rm = TRUE)),
    .groups = "drop"
  )

# Add Overall
cont_overall <- df %>%
  select(all_of(cont_vars)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(
    Value = sprintf("%.1f ± %.1f", mean(Value, na.rm = TRUE), sd(Value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(season = "Overall")

cont_summary <- bind_rows(cont_summary, cont_overall) %>%
  pivot_wider(names_from = season, values_from = Value)

cat_summary <- df %>%
  select(all_of(c("season", cat_vars))) %>%
  pivot_longer(-season, names_to = "Variable", values_to = "Level") %>%
  filter(!is.na(Level)) %>%
  group_by(season, Variable, Level) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(season, Variable) %>%
  mutate(
    pct = round(100 * N / sum(N), 1),
    Value = sprintf("%d (%.1f%%)", N, pct)
  ) %>%
  ungroup()

# Add Overall
cat_overall <- df %>%
  pivot_longer(all_of(cat_vars), names_to = "Variable", values_to = "Level") %>%
  filter(!is.na(Level)) %>%
  group_by(Variable, Level) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Variable) %>%
  mutate(
    pct = round(100 * N / sum(N), 1),
    Value = sprintf("%d (%.1f%%)", N, pct),
    season = "Overall"
  )

cat_summary <- bind_rows(cat_summary, cat_overall) %>%
  mutate(Variable = paste0(Variable, "=", Level)) %>%
  select(-Level, -N, -pct) %>%
  pivot_wider(names_from = season, values_from = Value, values_fill = "")

# Make sure everything is character
cat_summary[] <- lapply(cat_summary, as.character)


baseline_table <- bind_rows(cont_summary, cat_summary)

# Save
write.csv(baseline_table, "~/output/tables/baseline_characteristics_nocounts.csv", row.names = FALSE)


# ------------------------------
# Continuous summary
# ------------------------------
cont_summary <- df %>%
  select(all_of(c("season", cont_vars))) %>%
  pivot_longer(-season_mi, names_to = "Variable", values_to = "Value") %>%
  group_by(season, Variable) %>%
  summarise(
    n    = sum(!is.na(Value)),                  # number of non-missing values
    mean = mean(Value, na.rm = TRUE),
    sd   = sd(Value, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    Value = as.character(sprintf("%.1f ± %.1f (n=%d)", mean, sd, n))
  ) %>%
  select(Variable, season_mi, Value) %>%
  pivot_wider(names_from = season_mi, values_from = Value)

# Ensure character type across the board
cont_summary[] <- lapply(cont_summary, as.character)

# ------------------------------
# Categorical summary (expand levels)
# ------------------------------
cat_summary <- df %>%
  select(all_of(c("season_mi", cat_vars))) %>%
  pivot_longer(-season_mi, names_to = "Variable", values_to = "Level") %>%
  filter(!is.na(Level)) %>%
  group_by(season_mi, Variable, Level) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(season_mi, Variable) %>%
  mutate(
    pct   = round(100 * N / sum(N), 1),
    Value = paste0(N, " (", pct, "%)")
  ) %>%
  ungroup() %>%
  mutate(Variable = paste0(Variable, "=", Level)) %>%
  select(-Level, -N, -pct) %>%
  pivot_wider(
    names_from  = season_mi,
    values_from = Value,
    values_fill = ""        # make sure empty cells are ""
  )

# force all to character to avoid list/char mismatch
cat_summary <- mutate_all(cat_summary, as.character)



# ------------------------------
# Combine continuous + categorical
# ------------------------------
baseline_table <- bind_rows(cont_summary, cat_summary)

# Save
write.csv(baseline_table, "baseline_characteristics_by_season_MI.csv", row.names = FALSE)


cont_summary <- df %>%
  select(all_of(c("season", cont_vars))) %>%
  pivot_longer(-season, names_to = "Variable", values_to = "Value") %>%
  group_by(season, Variable) %>%
  summarise(
    Value = sprintf("%.1f ± %.1f", mean(Value, na.rm = TRUE), sd(Value, na.rm = TRUE)),
    .groups = "drop"
  )

# Add Overall
cont_overall <- df %>%
  select(all_of(cont_vars)) %>%
  pivot_longer(everything(), names_to = "Variable", values_to = "Value") %>%
  group_by(Variable) %>%
  summarise(
    Value = sprintf("%.1f ± %.1f", mean(Value, na.rm = TRUE), sd(Value, na.rm = TRUE)),
    .groups = "drop"
  ) %>%
  mutate(season = "Overall")

cont_summary <- bind_rows(cont_summary, cont_overall) %>%
  pivot_wider(names_from = season, values_from = Value)

cat_summary <- df %>%
  select(all_of(c("season", cat_vars))) %>%
  pivot_longer(-season, names_to = "Variable", values_to = "Level") %>%
  filter(!is.na(Level)) %>%
  group_by(season, Variable, Level) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(season, Variable) %>%
  mutate(
    pct = round(100 * N / sum(N), 1),
    Value = sprintf("%d (%.1f%%)", N, pct)
  ) %>%
  ungroup()

# Add Overall
cat_overall <- df %>%
  pivot_longer(all_of(cat_vars), names_to = "Variable", values_to = "Level") %>%
  filter(!is.na(Level)) %>%
  group_by(Variable, Level) %>%
  summarise(N = n(), .groups = "drop") %>%
  group_by(Variable) %>%
  mutate(
    pct = round(100 * N / sum(N), 1),
    Value = sprintf("%d (%.1f%%)", N, pct),
    season = "Overall"
  )

cat_summary <- bind_rows(cat_summary, cat_overall) %>%
  mutate(Variable = paste0(Variable, "=", Level)) %>%
  select(-Level, -N, -pct) %>%
  pivot_wider(names_from = season, values_from = Value, values_fill = "")

# Make sure everything is character
cat_summary[] <- lapply(cat_summary, as.character)


baseline_table <- bind_rows(cont_summary, cat_summary)

# Save
write.csv(baseline_table, "baseline_characteristics_by_season_andoverall.csv", row.names = FALSE)


# ================================================================
# Task 2.2: Event Distribution
# Calculate SCD incidence rates by month and season 
# with 95% Poisson confidence intervals
# ================================================================

# Load packages
library(dplyr)
library(epitools)   # for poisson.exact CI
library(lubridate)
library(purrr)
############################################################
# EVENT DISTRIBUTION
# Calculate SCD incidence rates by month and season
# with 95% Poisson confidence intervals
############################################################
df_surv <- df   # change to df if needed

#-----------------------------------------------------------
# 1. Basic preparation: event indicator + person-time (years)
#-----------------------------------------------------------

df_surv <- df_surv %>%
  mutate(
    # SCD event indicator: 1 = SCD, 0 = no SCD
    event_scd = ifelse(Status == 1, 1, 0),
    
    # Person-time in years (Survival_time is in months)
    py = Survival_time / 12
  )

# Optional quick checks
summary(df_surv$Survival_time)
table(df_surv$Status, useNA = "ifany")
table(df_surv$event_date, useNA = "ifany")

#-----------------------------------------------------------
# 2. Make sure Season and Month variables exist
#    Season: based on follow-up start date
#    Month: calendar month of follow-up start date (1–12)
#-----------------------------------------------------------

# If you already have followup_start_date and Season, skip this block.
if (!"followup_start_date" %in% names(df_surv)) {
  df_surv$Time_zero_Ym <- as.character(df_surv$Time_zero_Ym)
  df_surv$followup_start_date <- as.Date(paste0(df_surv$Time_zero_Ym, "-01"))
}

if (!"Season" %in% names(df_surv)) {
  df_surv$month_start <- month(df_surv$followup_start_date)
  
  df_surv$Season <- case_when(
    df_surv$month_start %in% c(12, 1, 2)  ~ "Winter",
    df_surv$month_start %in% c(3, 4, 5)   ~ "Spring",
    df_surv$month_start %in% c(6, 7, 8)   ~ "Summer",
    df_surv$month_start %in% c(9, 10, 11) ~ "Autumn",
    TRUE ~ NA_character_
  )
  
  df_surv$Season <- factor(df_surv$Season,
                           levels = c("Winter", "Spring", "Summer", "Autumn"))
}

# Month number (1–12) and label (Jan, Feb, ...)
df_surv$month_start <- month(df_surv$followup_start_date)
df_surv$month_label <- factor(
  df_surv$month_start,
  levels = 1:12,
  labels = c("Jan","Feb","Mar","Apr","May","Jun",
             "Jul","Aug","Sep","Oct","Nov","Dec")
)

#-----------------------------------------------------------
# Helper function: get Poisson 95% CI per 1000 person-years
#-----------------------------------------------------------

get_poisson_ci <- function(events, py) {
  # If no person-time, return NA
  if (is.na(py) || py <= 0) return(c(NA_real_, NA_real_))
  
  pt <- poisson.test(events, T = py)  # base R
  c(pt$conf.int[1] * 1000, pt$conf.int[2] * 1000)
}

#-----------------------------------------------------------
# 3. Incidence rates by SEASON
#-----------------------------------------------------------

inc_season <- df_surv %>%
  group_by(Season) %>%
  summarise(
    events = sum(event_scd, na.rm = TRUE),
    person_time = sum(py, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    rate_per_1000 = (events / person_time) * 1000
  )

# Apply Poisson CI for each season
ci_season <- mapply(
  get_poisson_ci,
  events = inc_season$events,
  py     = inc_season$person_time
)

inc_season$lower_95 <- ci_season[1, ]
inc_season$upper_95 <- ci_season[2, ]

# Optional: overall rate (all seasons combined)
overall <- df_surv %>%
  summarise(
    Season = "Overall",
    events = sum(event_scd, na.rm = TRUE),
    person_time = sum(py, na.rm = TRUE)
  ) %>%
  mutate(
    rate_per_1000 = (events / person_time) * 1000
  )

ci_overall <- get_poisson_ci(overall$events, overall$person_time)
overall$lower_95 <- ci_overall[1]
overall$upper_95 <- ci_overall[2]

inc_season <- bind_rows(inc_season, overall)

# View season-level incidence
inc_season

write.csv(inc_season, "T:/PROFID/output/incidence_by_season.csv", row.names=FALSE)

#-----------------------------------------------------------
# 4. Incidence rates by CALENDAR MONTH of follow-up start
#-----------------------------------------------------------

inc_month <- df_surv %>%
  group_by(month_start, month_label) %>%
  summarise(
    events = sum(event_scd, na.rm = TRUE),
    person_time = sum(py, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(month_start) %>%
  mutate(
    rate_per_1000 = (events / person_time) * 1000
  )

ci_month <- mapply(
  get_poisson_ci,
  events = inc_month$events,
  py     = inc_month$person_time
)

inc_month$lower_95 <- ci_month[1, ]
inc_month$upper_95 <- ci_month[2, ]

# View month-level incidence
inc_month

write.csv(inc_month, "T:/PROFID/output/incidence_by_month.csv", row.names = FALSE)


##############################################
# 3. PLOT — FOREST PLOT FOR SEASONAL INCIDENCE
##############################################

p_season <- ggplot(inc_season, aes(x=Season, y=rate_per_1000)) +
  geom_point(size=3, color="black") +
  geom_errorbar(aes(ymin=lower_95, ymax=upper_95), width=0.2) +
  labs(
    title="SCD Incidence Rate by Season",
    x="Season", y="Incidence rate per 1000 person-years"
  ) +
  theme_minimal(base_size=14)

ggsave("T:/PROFID/output/season_incidence_forest.png", p_season, width=7, height=5)

##############################################
# 4. PLOT — MONTHLY LINE PLOT WITH CI RIBBON
##############################################

p_month <- ggplot(inc_month, aes(x=month_start, y=rate_per_1000)) +
  geom_line(size=1) +
  geom_point(size=2) +
  geom_ribbon(aes(ymin=lower_95, ymax=upper_95), alpha=0.2) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  labs(
    title="SCD Incidence Rate by Month",
    x="Month", y="Incidence rate per 1000 person-years"
  ) +
  theme_minimal(base_size=14)

ggsave("T:/PROFID/output/monthly_incidence_plot.png", p_month, width=8, height=5)


# --- Combine all ---
inc_all <- bind_rows(overall, inc_season, inc_month)

print(inc_all)

# Save
write.csv(inc_all, "T:/PROFID/output/inc_all.csv", row.names = FALSE)


# Generate time series plots showing SCD events over calendar time
# Construct Event Date from Survival Time



df_events <- df %>%
  mutate(
    # convert follow-up start to date
    start_date = as.Date(paste0(Time_zero_Ym, "-01")),
    
    # approximate event date
    event_date = start_date + Survival_time * 30.44,
    
    # extract year-month for aggregation
    event_month = floor_date(event_date, "month"),
    
    # define SCD event
    event_scd = ifelse(Status == 1, 1, 0)
  )

# Aggregate Monthly SCD Events

scd_monthly <- df_events %>%
  group_by(event_month) %>%
  summarise(
    scd_events = sum(event_scd, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  arrange(event_month)

# Check if dates are correct
summary(df_events$event_date)
sum(is.na(df_events$event_date))

# Time Series Plot

p_trend <- ggplot(scd_monthly, aes(x = event_month, y = scd_events)) +
  geom_line(size = 1, color = "steelblue") +
  geom_point(size = 2, color = "steelblue") +
  labs(
    title = "Temporal Trends in SCD Events Over Calendar Time",
    x = "Calendar Time",
    y = "Number of SCD Events"
  ) +
  theme_minimal(base_size = 14)

library(ggplot2)
library(scales)

p_trend <- ggplot(scd_monthly, aes(x = event_month, y = scd_events)) +
  
  # Line
  geom_line(linewidth = 0.8, colour = "#2C3E50") +
  
  # # Optional trend smoother
  # geom_smooth(method = "loess", se = FALSE, colour = "#E74C3C", linewidth = 0.9) +
  # 
  # Axes formatting
  scale_x_date(date_breaks = "1 year",
               date_labels = "%Y",
               expand = expansion(mult = c(0.01, 0.01))) +
  scale_y_continuous(expand = expansion(mult = c(0.02, 0.05))) +
  
  labs(
    title = "Monthly Sudden Cardiac Death (SCD) Events Over Time",
    x = "Calendar Year",
    y = "Number of SCD Events"
  ) +
  
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank()
  )

p_trend


print(p_trend)

ggsave("T:/PROFID/output/temporal_trends_scd.png", p_trend, width = 10, height = 5)

# # EXAMPLE: If we want to visually highlight seasons:
# seasons_df <- data.frame(
#   start = as.Date(c("2000-12-01","2001-12-01","2002-12-01")),
#   end   = as.Date(c("2001-02-28","2002-02-28","2003-02-28"))
# )
# 
# p_trend +
#   geom_rect(data = seasons_df,
#             aes(xmin = start, xmax = end, ymin = -Inf, ymax = Inf),
#             fill = "#BDC3C7", alpha = 0.2, inherit.aes = FALSE)


# Task: Covariate Balance — Assess distribution of key risk factors (LVEF, age, sex, comorbidities) across seasonal group

library(dplyr)
library(tableone)
library(ggplot2)

df_cb <- df

# Continuous variables
cont_vars <- c("Age", "LVEF_std")

# Categorical variables
cat_vars <- c(
  "Sex",
  "Diabetes",
  "Hypertension",
  "HF",
  "MI_history",
  "AF_atrial_flutter",
  "COPD",
  "Stroke_TIA"
)

vars <- c(cont_vars, cat_vars)

table_cb <- CreateTableOne(
  vars = vars,
  strata = "season",
  data = df_cb,
  factorVars = cat_vars
)

print(table_cb, showAllLevels = TRUE, quote = FALSE, noSpaces = TRUE)

# Boxplots for continuous variables
# Age
ggplot(df_cb, aes(x = season, y = Age, fill = season)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(title = "Age Distribution by Season",
       x = "Season", y = "Age (years)") +
  guides(fill = FALSE)

# LVEF
ggplot(df_cb, aes(x = season, y = LVEF_std, fill = season)) +
  geom_boxplot(alpha = 0.7) +
  theme_minimal(base_size = 14) +
  labs(title = "LVEF Distribution by Season",
       x = "Season", y = "LVEF (%)") +
  guides(fill = FALSE)

# Bar plots for categorical variables

ggplot(df_cb, aes(x = season, fill = Sex)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  labs(title = "Sex Distribution by Season",
       y = "Percentage",
       x = "Season")

ggplot(df_cb, aes(x = season, fill = Sex)) +
  geom_bar(position = "fill") +
  scale_y_continuous(labels = scales::percent) +
  theme_minimal(base_size = 14) +
  labs(title = "Sex Distribution by Season",
       y = "Percentage",
       x = "Season")

library(cobalt)

bal.tab(season ~ Age + LVEF_std + Sex + Diabetes + Hypertension +
          HF + MI_history + AF_atrial_flutter + COPD + Stroke_TIA,
        data = df_cb)



# Create Table 1 stratified by season
table1 <- CreateTableOne(vars = covariates, strata = "season", data = df, test = TRUE)

# Print nicely
print(table1, showAllLevels = TRUE, quote = TRUE, nospaces = TRUE)

# CODE TO: Automatically generate all covariate balance plots
# Saves every plot as a high-quality PNG
# Handles continuous + categorical variables
# Organises files into folders so everything is tidy

library(dplyr)
library(ggplot2)
library(forcats)

# -------------------------------
# Prep: Create output folders
# -------------------------------
dir.create("plots", showWarnings = FALSE)
dir.create("plots/covariate_balance", showWarnings = FALSE)
dir.create("plots/covariate_balance/continuous", showWarnings = FALSE)
dir.create("plots/covariate_balance/categorical", showWarnings = FALSE)

# -------------------------------
# Data used for plotting
# -------------------------------
df_cb <- df %>%
  mutate(
    season = factor(
      season,
      levels = c("Winter", "Spring", "Summer", "Autumn")
    )
  )

# -------------------------------
# Define variables
# -------------------------------

# Continuous
continuous_vars <- c("Age", "LVEF_std")

# Categorical
categorical_vars <- c(
  "Sex",
  "Diabetes",
  "Hypertension",
  "HF",
  "MI_history",
  "AF_atrial_flutter",
  "COPD",
  "Stroke_TIA"
)

# -------------------------------
# LOOP FOR CONTINUOUS VARIABLES
# -------------------------------

for (var in continuous_vars) {
  
  p <- ggplot(df_cb, aes_string(x = "season", y = var, fill = "season")) +
    geom_boxplot(alpha = 0.7, outlier.alpha = 0.2) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(var, "Distribution by Season"),
      x = "Season",
      y = var
    ) +
    guides(fill = "none")
  
  # Save plot
  ggsave(
    filename = paste0("plots/covariate_balance/continuous/", var, "_by_season.png"),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}

# -------------------------------
# LOOP FOR CATEGORICAL VARIABLES
# -------------------------------

for (var in categorical_vars) {
  
  p <- df_cb %>%
    filter(!is.na(.data[[var]])) %>%   # Ensure no NA in variable
    ggplot(aes(x = season, fill = .data[[var]])) +
    geom_bar(position = "fill") +
    scale_y_continuous(labels = scales::percent_format()) +
    theme_minimal(base_size = 14) +
    labs(
      title = paste(var, "Distribution by Season"),
      x = "Season",
      y = "Percentage",
      fill = var
    )
  
  # Save plot
  ggsave(
    filename = paste0("plots/covariate_balance/categorical/", var, "_by_season.png"),
    plot = p,
    width = 8,
    height = 5,
    dpi = 300
  )
}
