###############################################
# Moran's I statistic (ICD and Non-ICD cohorts)
# - Crude rate computation
# - Combine GT table
###############################################

# Packages
req <- c("data.table","dplyr","tidyverse","sf","spdep","tibble","gt")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# ICD crude rate
ICD_ctr <- fread("T:/Data Transfer to Charite/raw/ICD_filtered_with_coords.csv") %>%
  select(ctr_name, DB, latitude, longitude) %>%
  distinct()
ICD_age_std_rate <- fread("T:/Data Transfer to Charite/raw/ICD_filtered_with_coords.csv") %>%
  group_by(ctr_name) %>%
  summarise(
    events = sum(Status == 1, na.rm = TRUE),
    person_years = sum(Survival_time, na.rm = TRUE),
    crude_rate = events / person_years * 100000,
    .groups = "drop"
  )
ICD_geo_df <- ICD_age_std_rate %>%
  left_join(
    ICD_ctr %>% select(ctr_name, DB, latitude, longitude),
    by = "ctr_name"
  )
ICD_geo_df <- ICD_geo_df %>%
  mutate(
    crude_rate_plot = ifelse(is.na(crude_rate), 0, crude_rate),
    rate_missing = is.na(crude_rate)
  )
fwrite(ICD_geo_df, "T:/Data Transfer to Charite/raw/ICD_raw_geo_df.csv")
ICD_geo_sf <- ICD_geo_df %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Non-ICD reduced crude rate
NR_ctr <- fread("T:/Data Transfer to Charite/raw/NonICD_reduced_filtered_with_coords.csv") %>%
  select(DB, latitude, longitude) %>%
  distinct()
NR_age_std_rate <- fread("T:/Data Transfer to Charite/raw/NonICD_reduced_filtered_with_coords.csv") %>%
  group_by(DB) %>%
  summarise(
    events = sum(Status == 1, na.rm = TRUE),
    person_years = sum(Survival_time, na.rm = TRUE),
    crude_rate = events / person_years * 100000,
    .groups = "drop"
  )
NR_geo_df <- NR_age_std_rate %>%
  left_join(
    NR_ctr %>% select(DB, latitude, longitude),
    by = "DB"
  )
NR_geo_df <- NR_geo_df %>%
  mutate(
    crude_rate_plot = ifelse(is.na(crude_rate), 0, crude_rate),
    rate_missing = is.na(crude_rate)
  )
fwrite(NR_geo_df, "T:/Data Transfer to Charite/raw/NR_raw_geo_df.csv")
NR_geo_sf <- NR_geo_df %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Non-ICD preserved crude rate
NP_ctr <- fread("T:/Data Transfer to Charite/raw/NonICD_preserved_filtered_with_coords.csv") %>%
  select(DB, latitude, longitude) %>%
  distinct()
NP_age_std_rate <- fread("T:/Data Transfer to Charite/raw/NonICD_preserved_filtered_with_coords.csv") %>%
  group_by(DB) %>%
  summarise(
    events = sum(Status == 1, na.rm = TRUE),
    person_years = sum(Survival_time, na.rm = TRUE),
    crude_rate = events / person_years * 100000,
    .groups = "drop"
  )
NP_geo_df <- NP_age_std_rate %>%
  left_join(
    NP_ctr %>% select(DB, latitude, longitude),
    by = "DB"
  )
NP_geo_df <- NP_geo_df %>%
  mutate(
    crude_rate_plot = ifelse(is.na(crude_rate), 0, crude_rate),
    rate_missing = is.na(crude_rate)
  )
fwrite(NP_geo_df, "T:/Data Transfer to Charite/raw/NP_raw_geo_df.csv")
NP_geo_sf <- NP_geo_df %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

# Moran's I function
run_moran_one <- function(geo_sf, cohort_label, k = 4) {
  # Drop centres with missing incidence
  geo_sf_clean <- geo_sf %>% filter(!is.na(crude_rate))
  n_centres <- nrow(geo_sf_clean)
  # If too few centres, return NA row (important!)
  if (n_centres < (k + 1)) {
    return(
      tibble(
        Cohort = cohort_label,
        n_centres = n_centres,
        Moran_I = NA_real_,
        Expected_I = NA_real_,
        Variance = NA_real_,
        z_value = NA_real_,
        p_value = NA_real_,
        Interpretation = "Too few centres for spatial analysis"
      )
    )
  }
  coords <- st_coordinates(geo_sf_clean)
  nb <- knn2nb(knearneigh(coords, k = k))
  lw <- nb2listw(nb, style = "W", zero.policy = TRUE)
  moran_res <- moran.test(
    geo_sf_clean$crude_rate,
    lw,
    zero.policy = TRUE
  )
  tibble(
    Cohort = cohort_label,
    n_centres = n_centres,
    Moran_I = as.numeric(moran_res$estimate["Moran I statistic"]),
    Expected_I = as.numeric(moran_res$estimate["Expectation"]),
    Variance = as.numeric(moran_res$estimate["Variance"]),
    z_value = as.numeric(moran_res$statistic),
    p_value = moran_res$p.value,
    Interpretation = ifelse(
      moran_res$p.value < 0.05,
      "Evidence of global spatial autocorrelation",
      "No global spatial autocorrelation"
    )
  )
}
moran_raw <- bind_rows(
  run_moran_one(ICD_geo_sf,  "ICD"),
  run_moran_one(NR_geo_sf,   "Non-ICD ≤35%"),
  run_moran_one(NP_geo_sf,   "Non-ICD >35%")
)

# GT table
make_moran_gt <- function(moran_df, table_title) {
  moran_df %>%
    mutate(
      Moran_I   = round(Moran_I, 3),
      Expected_I = round(Expected_I, 3),
      z_value   = round(z_value, 2),
      p_value   = round(p_value, 3)
    ) %>%
    gt() %>%
    tab_header(
      title = table_title,
      subtitle = "Global spatial autocorrelation (Moran’s I)"
    ) %>%
    cols_label(
      Cohort = "Cohort",
      n_centres = "Centres (n)",
      Moran_I = "Moran’s I",
      Expected_I = "Expected I",
      Variance = "Variance",
      z_value = "Z statistic",
      p_value = "P-value",
      Interpretation = "Interpretation"
    ) %>%
    fmt_number(
      columns = c(Moran_I, Expected_I, Variance, z_value, p_value),
      decimals = 3
    ) %>%
    tab_source_note(
      source_note = md(
        "*Moran’s I computed using k-nearest neighbour spatial weights (k = 4).  
        Results for Non-ICD cohorts should be interpreted cautiously due to small numbers of contributing centres.*"
      )
    )
}
gt_moran_raw <- make_moran_gt(
  moran_raw,
  "Moran’s I – Raw data"
)
gtsave(gt_moran_raw, "T:/Data Transfer to Charite/raw/raw_morans.html")

