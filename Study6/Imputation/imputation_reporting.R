library(mice)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

imp <- readRDS("mice_imputed_data.RDS")
pre <- imp$data
post_long <- complete(imp, action = "long", include = TRUE) %>%
  filter(.imp != 0)  # exclude the original incomplete data

# ---------- Continuous summary (mean + SD) ----------
num_vars <- names(pre)[sapply(pre, is.numeric)]
num_vars <- setdiff(num_vars, c(".imp", ".id"))

cont_summary <- lapply(num_vars, function(v){
  pre_x <- pre[[v]]
  pre_missing_n <- sum(is.na(pre_x))
  pre_missing_pct <- 100 * pre_missing_n / length(pre_x)
  
  pre_mean <- mean(pre_x, na.rm = TRUE)
  pre_sd   <- sd(pre_x,   na.rm = TRUE)
  
  # post: compute mean/sd within each imputation, then pool by averaging
  post_stats <- post_long %>%
    group_by(.imp) %>%
    summarise(
      mean_i = mean(.data[[v]], na.rm = TRUE),
      sd_i   = sd(.data[[v]],   na.rm = TRUE),
      .groups = "drop"
    )
  
  post_mean <- mean(post_stats$mean_i, na.rm = TRUE)
  post_sd   <- mean(post_stats$sd_i,   na.rm = TRUE)  # descriptive pooling of SD
  
  tibble(
    variable = v,
    missing_n = pre_missing_n,
    missing_pct = pre_missing_pct,
    pre_mean = pre_mean,
    pre_sd = pre_sd,
    post_mean = post_mean,
    post_sd = post_sd
  )
}) %>% bind_rows()

# round for reporting (tweak decimals as you like)
cont_summary_out <- cont_summary %>%
  mutate(
    missing_pct = round(missing_pct, 2),
    pre_mean = round(pre_mean, 3),
    pre_sd   = round(pre_sd, 3),
    post_mean = round(post_mean, 3),
    post_sd   = round(post_sd, 3),
    mean_diff = round(post_mean - pre_mean, 3),
    sd_diff   = round(post_sd - pre_sd, 3)
  ) %>%
  arrange(desc(missing_pct))

write.csv(cont_summary_out, "mice_pre_post_continuous_mean_sd.csv", row.names = FALSE)

# ---------- Plot 1: Missingness (%), numeric variables ----------
p_miss <- cont_summary_out %>%
  ggplot(aes(x = reorder(variable, missing_pct), y = missing_pct)) +
  geom_col() +
  coord_flip() +
  labs(x = NULL, y = "% missing (pre-imputation)",
       title = "Missingness by variable (numeric)") +
  theme_bw()

ggsave("plot_missingness_numeric.png", p_miss, width = 8, height = 10, dpi = 300)

# ---------- Plot 2: Pre vs post mean (scatter) ----------
p_mean <- cont_summary_out %>%
  ggplot(aes(x = pre_mean, y = post_mean)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Pre-imputation mean", y = "Post-imputation mean (pooled)",
       title = "Pre vs post means (numeric)") +
  theme_bw()

ggsave("plot_pre_post_means.png", p_mean, width = 6, height = 5, dpi = 300)

# ---------- Plot 3: Pre vs post SD (scatter) ----------
p_sd <- cont_summary_out %>%
  ggplot(aes(x = pre_sd, y = post_sd)) +
  geom_point() +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
  labs(x = "Pre-imputation SD", y = "Post-imputation SD (pooled)",
       title = "Pre vs post SDs (numeric)") +
  theme_bw()

ggsave("plot_pre_post_sds.png", p_sd, width = 6, height = 5, dpi = 300)

library(ggplot2)
library(dplyr)

imp <- readRDS("mice_imputed_data.RDS")
pre <- imp$data
post_long <- complete(imp, action = "long", include = TRUE) %>% filter(.imp != 0)

num_vars <- names(pre)[sapply(pre, is.numeric)]
num_vars <- setdiff(num_vars, c(".imp", ".id"))

# pick top K by missingness
missing_rank <- tibble(
  variable = num_vars,
  missing_pct = sapply(num_vars, function(v) mean(is.na(pre[[v]])) * 100)
) %>% arrange(desc(missing_pct))

top_k <- 12  # change as needed
vars_to_plot <- missing_rank$variable[1:min(top_k, nrow(missing_rank))]

for(v in vars_to_plot){
  pre_df <- tibble(value = pre[[v]], dataset = "Pre (observed)") %>% filter(!is.na(value))
  
  # pool post across imputations by sampling equal numbers from each imputation (keeps it light)
  post_df <- post_long %>%
    group_by(.imp) %>%
    slice_sample(prop = 0.1, replace = FALSE) %>%  # 10% per imputation
    ungroup() %>%
    transmute(value = .data[[v]], dataset = "Post (imputed)")
  
  
  plot_df <- bind_rows(pre_df, post_df)
  
  p <- ggplot(plot_df, aes(x = value, linetype = dataset)) +
    geom_density() +
    labs(title = paste("Distribution check:", v), x = v, y = "Density") +
    theme_bw()
  
  ggsave(paste0("density_pre_post_", v, ".png"), p, width = 6, height = 4, dpi = 300)
}

library(dplyr)
library(tidyr)
library(ggplot2)

imp <- readRDS("mice_imputed_data.RDS")
pre <- imp$data
post_long <- complete(imp, action = "long", include = TRUE) %>% filter(.imp != 0)

cat_vars <- names(pre)[sapply(pre, function(x) is.factor(x) || is.character(x) || is.logical(x))]

for(v in cat_vars){
  pre_tab <- prop.table(table(pre[[v]], useNA = "no"))
  pre_df <- tibble(level = names(pre_tab), prop = as.numeric(pre_tab), dataset = "Pre (observed)")
  
  post_df <- post_long %>%
    group_by(.imp, .data[[v]]) %>%
    summarise(n = n(), .groups = "drop_last") %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    rename(level = .data[[v]]) %>%
    group_by(level) %>%
    summarise(prop = mean(prop, na.rm = TRUE), .groups = "drop") %>%
    mutate(dataset = "Post (imputed)")
  
  plot_df <- bind_rows(pre_df, post_df) %>%
    mutate(level = as.character(level))
  
  p <- ggplot(plot_df, aes(x = level, y = prop, fill = dataset)) +
    geom_col(position = "dodge") +
    labs(title = paste("Proportions:", v), x = NULL, y = "Proportion") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  ggsave(paste0("prop_pre_post_", v, ".png"), p, width = 7, height = 4, dpi = 300)
}



# Choose 3 continuous + 1 categorical (edit as needed)
cont_vars <- c("eGFR", "Triglycerides", "Haemoglobin")   # or swap in Triglycerides, etc.
cat_var   <- "Smoking"                         # or "Smoking_status" etc.

# ---- Helper: long format with imputation indicator ----
# .imp = 0 is the original data; .imp = 1..m are completed datasets
long <- complete(imp, action = "long", include = TRUE) %>%
  as_tibble()

# Identify which cells were imputed in the original data
# imp$where is a logical matrix: rows = cases, cols = variables
where <- as_tibble(imp$where, .name_repair = "minimal") %>%
  mutate(.id = row_number())

# Add row id to long so we can match back to 'where'
if (!(".id" %in% names(long))) {
  # mice long usually has .id already; if not, create it from rownames
  long <- long %>% mutate(.id = as.integer(.id))
}

# ---- Helper function: density overlay for a continuous variable ----
plot_density_imp <- function(var) {
  # observed = values that were originally observed (not missing) in dat
  # imputed  = values filled in for originally-missing cells, pooled across imputations
  
  # Safety: skip if not in where/long
  if (!(var %in% names(where)) || !(var %in% names(long))) {
    return(ggplot() + theme_void() + labs(title = paste(var, "(not found)")))
  }
  
  df_obs <- long %>%
    filter(.imp == 0) %>%                      # original data
    select(.id, value = all_of(var)) %>%
    left_join(where %>% select(.id, was_imputed = all_of(var)), by = ".id") %>%
    filter(!is.na(value), was_imputed == FALSE) %>%
    mutate(source = "Observed")
  
  df_imp <- long %>%
    filter(.imp != 0) %>%                      # completed datasets
    select(.id, value = all_of(var)) %>%
    left_join(where %>% select(.id, was_imputed = all_of(var)), by = ".id") %>%
    filter(!is.na(value), was_imputed == TRUE) %>%
    mutate(source = "Imputed")
  
  df <- bind_rows(df_obs, df_imp)
  
  ggplot(df, aes(x = value, linetype = source)) +
    geom_density(linewidth = 0.8, na.rm = TRUE) +
    labs(
      title = var,
      x = NULL,
      y = "Density",
      linetype = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

# ---- Helper function: bar plot for categorical variable ----
plot_bar_imp <- function(var) {
  if (!(var %in% names(where)) || !(var %in% names(long))) {
    return(ggplot() + theme_void() + labs(title = paste(var, "(not found)")))
  }
  
  df_obs <- long %>%
    filter(.imp == 0) %>%
    select(.id, value = all_of(var)) %>%
    left_join(where %>% select(.id, was_imputed = all_of(var)), by = ".id") %>%
    filter(!is.na(value), was_imputed == FALSE) %>%
    mutate(source = "Observed")
  
  df_imp <- long %>%
    filter(.imp != 0) %>%
    select(.id, value = all_of(var)) %>%
    left_join(where %>% select(.id, was_imputed = all_of(var)), by = ".id") %>%
    filter(!is.na(value), was_imputed == TRUE) %>%
    mutate(source = "Imputed")
  
  df <- bind_rows(df_obs, df_imp) %>%
    mutate(value = as.factor(value)) %>%
    count(source, value) %>%
    group_by(source) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup()
  
  ggplot(df, aes(x = value, y = prop, fill = source)) +
    geom_col(position = "dodge", width = 0.8) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
    labs(
      title = var,
      x = NULL,
      y = "Proportion",
      fill = NULL
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "bottom")
}

# ---- Build panels ----
p1 <- plot_density_imp(cont_vars[1])
p2 <- plot_density_imp(cont_vars[2])
p3 <- plot_density_imp(cont_vars[3])
p4 <- plot_bar_imp(cat_var)

# Combined figure (2x2)
figSX <- (p1 | p2) / (p3 | p4) +
  plot_annotation(
    title = "Imputation diagnostics: observed vs imputed distributions",
    subtitle = "Observed = originally non-missing values; Imputed = values filled in for originally missing cells (pooled across imputations)."
  )

# Print to viewer
figSX

# Save
ggsave("Figure_SX_imputation_diagnostics.png", figSX, width = 12, height = 8, dpi = 300)
ggsave("Figure_SX_imputation_diagnostics.pdf", figSX, width = 12, height = 8)

