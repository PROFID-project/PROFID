##########################################################################
# Function: Split datasets into train and validation sets based on centre
##########################################################################

# Packages
req <- c("data.table","dplyr", "readr")
inst <- setdiff(req, rownames(installed.packages()))
if (length(inst)) install.packages(inst, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Function
split_by_centre_and_save <- function(
  df,
  dataset_name,
  id_col = "ID",
  centre_col = "DB",
  train_prop = 0.9,
  seed = 2026,
  out_dir = "data_splits"
) {
  
  # Create output directory if it doesn't exist
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
  }
  
  set.seed(seed)
  
  # Train set (centre-stratified)
  df_train <- df %>%
    group_by(.data[[centre_col]]) %>%
    slice_sample(prop = train_prop) %>%
    ungroup()
  
  # Validation set
  df_valid <- df %>%
    anti_join(df_train, by = id_col)
  
  # File paths
  train_path <- file.path(out_dir, paste0(dataset_name, "_train.csv"))
  valid_path <- file.path(out_dir, paste0(dataset_name, "_valid.csv"))
  
  # Save CSVs
  write_csv(df_train, train_path)
  write_csv(df_valid, valid_path)
  
  # Diagnostics
  diagnostics <- tibble::tibble(
    dataset = dataset_name,
    N_total = nrow(df),
    N_train = nrow(df_train),
    N_valid = nrow(df_valid),
    train_fraction = round(nrow(df_train) / nrow(df), 3),
    valid_fraction = round(nrow(df_valid) / nrow(df), 3),
    centres_total = n_distinct(df[[centre_col]]),
    centres_train = n_distinct(df_train[[centre_col]]),
    centres_valid = n_distinct(df_valid[[centre_col]])
  )
  
  return(list(
    train = df_train,
    valid = df_valid,
    diagnostics = diagnostics,
    files = c(train = train_path, valid = valid_path)
  ))
}

split_result <- split_by_centre_and_save(
  fread("T:/Data Transfer to Charite/imp/NonICD_preserved_imputed.csv"),
  dataset_name = "NP_imp",
  out_dir = "T:/Data Transfer to Charite/model_splits"
)

split_result$diagnostics

harmonise_centre_column <- function(df) {
  if ("ctr_name" %in% names(df)) {
    df <- df %>% rename(centre_id = ctr_name)
  } else if ("DB" %in% names(df)) {
    df <- df %>% rename(centre_id = DB)
  } else {
    stop("No centre column found (ctr_name or DB)")
  }
  df
}

fread("T:/Data Transfer to Charite/model_splits/NP_imp_valid.csv") <- harmonise_centre_column(fread("T:/Data Transfer to Charite/model_splits/NP_imp_valid.csv"))
