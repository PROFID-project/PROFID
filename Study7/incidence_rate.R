###############################################
# Incidence rates (ICD and Non-ICD cohorts)
# - Table generation
# - Heatmap
###############################################

# Packages
req <- c("data.table","dplyr","tidyverse","gt","ggplot2","scales")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Input file
df_ICD <- fread("T:/Data Transfer to Charite/imp/ICD_imputed.csv")
df_NR <- fread("T:/Data Transfer to Charite/imp/NonICD_reduced_imputed.csv")
df_NP <- fread("T:/Data Transfer to Charite/imp/NonICD_preserved_imputed.csv")

# Incidence rates
build_incidence_all <- function(df, cohort_label) {

  # --- Age groups ---
  df <- df %>%
    mutate(
      age_group = cut(
        Age,
        breaks = c(0, 50, 60, 70, 80, Inf),
        right = FALSE
      )
    )

  # --- Crude incidence ---
  crude <- df %>%
    group_by(CVD_risk_region) %>%
    summarise(
      events = sum(Status == 1, na.rm = TRUE),
      py = sum(Survival_time, na.rm = TRUE),
      crude_rate = events / py * 100000,
      crude_lci = (qchisq(0.025, 2 * events) / 2) / py * 100000,
      crude_uci = (qchisq(0.975, 2 * (events + 1)) / 2) / py * 100000,
      .groups = "drop"
    )

  # --- Age-specific rates ---
  age_specific <- df %>%
    group_by(CVD_risk_region, age_group) %>%
    summarise(
      events = sum(Status == 1, na.rm = TRUE),
      py = sum(Survival_time, na.rm = TRUE),
      rate = events / py,
      .groups = "drop"
    )

  # --- Standard population (study population) ---
  std_pop <- df %>%
    group_by(age_group) %>%
    summarise(std_py = sum(Survival_time, na.rm = TRUE), .groups = "drop")

  # --- Age-standardised rate + variance ---
  age_std <- age_specific %>%
    left_join(std_pop, by = "age_group") %>%
    mutate(
      weight = std_py / sum(std_py),
      var_comp = (events / (py^2)) * (weight^2)
    ) %>%
    group_by(CVD_risk_region) %>%
    summarise(
      age_std_rate = sum(rate * weight) * 100000,
      variance = sum(var_comp) * (100000^2),
      .groups = "drop"
    ) %>%
    mutate(
      se = sqrt(variance),
      age_lci = age_std_rate - 1.96 * se,
      age_uci = age_std_rate + 1.96 * se
    )

  # --- Combine ---
  out <- crude %>%
    left_join(age_std, by = "CVD_risk_region") %>%
    mutate(
      Cohort = cohort_label,
      crude_rate = round(crude_rate, 1),
      crude_CI = paste0(round(crude_lci, 1), "–", round(crude_uci, 1)),
      age_std_rate = round(age_std_rate, 1),
      age_CI = paste0(round(age_lci, 1), "–", round(age_uci, 1))
    ) %>%
    select(
      Cohort,
      CVD_risk_region,
      events,
      py,
      crude_rate,
      crude_CI,
      age_std_rate,
      age_CI
    )

  out
}

inc_ICD <- build_incidence_all(df_ICD, "ICD")
inc_NR <- build_incidence_all(df_NR, "Non-ICD (≤35%)")
inc_NP <- build_incidence_all(df_NP, "Non-ICD (>35%)")

incidence_stacked <- bind_rows(inc_ICD, inc_NR, inc_NP)

gt_incidence <- gt(
  incidence_stacked,
  groupname_col = "Cohort"
) %>%
  cols_label(
    CVD_risk_region = "Region",
    events = "SCD events",
    py = "Person-years",
    crude_rate = md("**Crude rate**"),
    crude_CI = "Crude 95% CI",
    age_std_rate = md("**Age-standardised rate**"),
    age_CI = "Age-std 95% CI"
  ) %>%
  cols_align(
    align = "center",
    columns = c(events, py, crude_rate, crude_CI, age_std_rate, age_CI)
  ) %>%
  tab_header(
    title = md("**Sudden cardiac death incidence by region**"),
    subtitle = md("Stacked ICD vs Non-ICD cohorts after imputation (rates per 100,000 person-years)")
  ) %>%
  tab_source_note(
    md(
      "Crude and age-standardised incidence rates with 95% confidence intervals. 
       Age standardisation used direct standardisation with the study population as the standard."
    )
  ) %>%
  tab_style(
    style = cell_text(weight = "bold"),
    locations = cells_row_groups()
  ) %>%
  data_color(
    columns = crude_rate,
    colors = col_numeric(
      palette = c("#f7fbff", "#c6dbef", "#6baed6", "#2171b5"),
      domain = range(incidence_stacked$crude_rate, na.rm = TRUE)
    )
  ) %>%
  data_color(
    columns = age_std_rate,
    colors = col_numeric(
      palette = c("#fff5f0", "#fcbba1", "#fb6a4a", "#cb181d"),
      domain = range(incidence_stacked$age_std_rate, na.rm = TRUE)
    )
  ) %>%
  tab_options(
    table.font.size = px(13),
    data_row.padding = px(4),
    row_group.padding = px(6)
  )
gtsave(gt_incidence, "T:/Data Transfer to Charite/imp/imp_SCD_incidence_stacked_crude_and_age_std.html")
