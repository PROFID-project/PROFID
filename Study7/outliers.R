###############################################
# Outlier Detection (ICD and Non-ICD cohorts)
# - Visual representations of distribution of biomarkers and physiological variables
# - Geographic-specific outliers detection algorithms 
# - Outliers documentation and heatmaps of outliers behaviour
###############################################

# Packages
req <- c("data.table","dplyr","tidyverse","ggplot2")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Input files
df <- fread("NonICD_reduced_imputed.csv")

# Distribution plots
# re-run for all the biomarkers and physiological variables of all the three datasets
p1 <- ggplot(df, aes(x = HR)) +
  geom_histogram(
    bins = 50,
    fill = "#55A868",
    color = "black",
    linewidth = 0.2
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "Distribution of HR",
    subtitle = "Non-ICD patients (>35%)",
    x = "HR (bpm)",
    y = "Count"
  )
ggsave(p1, filename = "NP_HR_distribution.png")

# Geographic-specific outliers detection algorithms
# Median Absolute Deviation (MAD) for skewed distributions
flag_outliers_mad_geo <- function(df, var, geo, k = 3.5, log_transform = TRUE) {
  var_sym <- rlang::ensym(var)
  geo_sym <- rlang::ensym(geo)
  df %>%
    mutate(
      !!paste0("log_", var_sym) :=
        if (log_transform) log1p(!!var_sym) else !!var_sym
    ) %>%
    group_by(!!geo_sym) %>%
    mutate(
      .median = median(!!sym(paste0("log_", var_sym)), na.rm = TRUE),
      .mad = mad(!!sym(paste0("log_", var_sym)), constant = 1, na.rm = TRUE),
      !!paste0("outlier_", var_sym) :=
        ifelse(
          .mad == 0 | is.na(!!sym(paste0("log_", var_sym))),
          FALSE,
          abs(!!sym(paste0("log_", var_sym)) - .median) / .mad > k
        )
    ) %>%
    ungroup() %>%
    select(-.median, -.mad)
}
mad_vars <- c("NTProBNP", "Troponin_T", "CRP")
for (v in mad_vars) {
  df <- flag_outliers_mad_geo(
    df,
    var = !!sym(v),
    geo = DB,
    k = 3.5,
    log_transform = TRUE
  )
}
# Interquartile Range (IQR) for normal distributions
flag_outliers_iqr_geo <- function(df, var, geo) {
  var_sym <- rlang::ensym(var)
  geo_sym <- rlang::ensym(geo)
  df %>%
    group_by(!!geo_sym) %>%
    mutate(
      Q1 = quantile(!!var_sym, 0.25, na.rm = TRUE),
      Q3 = quantile(!!var_sym, 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      !!paste0("outlier_", var_sym) :=
        !!var_sym < (Q1 - 1.5 * IQR) |
        !!var_sym > (Q3 + 1.5 * IQR)
    ) %>%
    ungroup() %>%
    select(-Q1, -Q3, -IQR)
}
iqr_vars <- c("LVEF", "BMI", "SBP", "DBP", "HR")
for (v in iqr_vars) {
  df <- flag_outliers_iqr_geo(
    df,
    var = !!sym(v),
    geo = DB
  )
}
all_vars <- c(mad_vars, iqr_vars)

# Outliers plots
# re-run for all the biomarkers and physiological variables of all the three datasets
p2 <- ggplot(df, aes(x = DB, y = NTProBNP)) +
  geom_boxplot(outlier.shape = NA, fill = "white", color = "black") +
  geom_point(
    data = df %>% filter(outlier_NTProBNP),
    aes(x = DB, y = NTProBNP),
    color = "red",
    size = 1,
    alpha = 0.6
  ) +
  coord_flip() +
  theme_classic(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  ) +
  labs(
    title = "NTProBNP outliers by centre",
    subtitle = "Non-ICD patients (≤35%)",
    x = "Centre",
    y = "NTProBNP (pmol/L)"
  )
ggsave(p2, filename = "NR_NTProBNP_outliers.png")

# Outliers documentation
outlier_doc <- lapply(all_vars, function(v) {
  outlier_col <- paste0("outlier_", v)
  df %>%
    group_by(DB) %>%
    summarise(
      variable = v,
      n = n(),
      n_outliers = sum(.data[[outlier_col]]),
      outlier_rate = n_outliers / n,
      max_value = max(.data[[v]], na.rm = TRUE),
      p99 = quantile(.data[[v]], 0.99, na.rm = TRUE),
      .groups = "drop"
    )
}) %>%
  bind_rows()
fwrite(outlier_doc, "NR_Centre_level_outlier_documentation.csv")
# order centres by overall outlier burden
outlier_doc <- outlier_doc %>%
  group_by(DB) %>%
  mutate(mean_rate = mean(outlier_rate, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(DB = reorder(DB, mean_rate))
# outliers rate heatmap
p3 <- ggplot(outlier_doc,
       aes(x = variable, y = DB, fill = outlier_rate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white",
    high = "#B2182B",
    name = "Outlier rate"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "Centre-level outlier rates by variable",
    subtitle = "Non-ICD patients (≤35%)",
    x = "Variable",
    y = "Centre"
  )
ggsave(p3, filename = "NR_Centre_level_outlier_heatmap.png")
# order centres by overall p99 burden
outlier_doc <- outlier_doc %>%
  group_by(DB) %>%
  mutate(mean_p99 = mean(p99, na.rm = TRUE)) %>%
  ungroup() %>%
  mutate(DB = reorder(DB, mean_p99))
# outliers 99th percentile heatmap
p4 <- ggplot(outlier_doc,
       aes(x = variable, y = DB, fill = p99)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradient(
    low = "white",
    high = "#2166AC",
    name = "99th percentile"
  ) +
  theme_classic(base_size = 13) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    legend.position = "right"
  ) +
  labs(
    title = "Centre-level 99th percentile values by variable",
    subtitle = "Non-ICD patients (≤35%)",
    x = "Variable",
    y = "Centre"
  )
ggsave(p4, filename = "NR_Centre_level_99th_percentile_heatmap.png")