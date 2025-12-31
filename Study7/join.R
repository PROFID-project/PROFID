############################################################
# Function: Join patient-level data with centre coordinates
############################################################

# Packages
req <- c("data.table","dplyr")
inst <- setdiff(req, rownames(installed.packages()))
if (length(inst)) install.packages(inst, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Function
join_with_centres <- function(
  patient_file,
  centre_file = "centreNonICD.csv",    # for ICD, centreIDSEUCERT.csv and for Non-ICD, centreNonICD.csv
  output_file = NULL
) {
  message("Processing: ", patient_file)
  
  # Read data
  df_pat <- fread(patient_file)
  df_ctr <- fread(centre_file)
  
  # Safety checks
  stopifnot("DB" %in% names(df_pat))
  stopifnot("ctr_name" %in% names(df_ctr))
  
  # Select only needed columns from centre file
  df_ctr <- df_ctr %>%
    select(ctr_name, latitude, longitude)
  
  # Perform LEFT JOIN (VLOOKUP equivalent)
  df_joined <- df_pat %>%
    left_join(df_ctr, by = c("DB" = "ctr_name"))
  
  # Optional: warn if unmatched IDs exist
  n_unmatched <- sum(is.na(df_joined$ctr_name))
  if (n_unmatched > 0) {
    warning(n_unmatched, " patients had no matching centre information.")
  }
  
  # Output filename
  if (is.null(output_file)) {
    output_file <- sub("\\.csv$", "_with_coords.csv", patient_file)
  }
  
  # Save result
  fwrite(df_joined, output_file)
  message("Saved: ", output_file)
  
  invisible(df_joined)
}

join_with_centres("NonICD_reduced_filtered.csv")
join_with_centres("NonICD_preserved_filtered.csv")