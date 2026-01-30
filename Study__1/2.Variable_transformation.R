###############################################################################
# Log Transformation of Right-Skewed Variables
# Complete visual check: Boxplots + Histograms
###############################################################################
library(ggplot2)
install.packages("gridExtra")
library(gridExtra) 
library(data.table)

cat("\n=== LOG TRANSFORMATION ===\n\n")

###############################################################################
# 1. LOAD DATA
###############################################################################

dlt <- readRDS("T:/FINAL ICD COHORT/standardised_data1.rds")
setDT(dlt)
names(dlt)
cat(sprintf("Loaded: %d observations, %d variables\n\n", nrow(dlt), ncol(dlt)))

###############################################################################
# 2. IDENTIFY CANDIDATES
###############################################################################

all_numeric <- names(dlt)[sapply(dlt, function(x) is.numeric(x) | is.integer(x))]

exclude <- c(
  "ID", "DB","V1",
  "Survival_time", "Time_FIS_days", "Time_index_MI_CHD","Time_death_days","Status_death","t_followup_days",
  "Status", "Status_FIS", "Time_zero_Ym","Time_zero_Y","CVD_risk_region",
  grep("^bin_", names(dlt), value = TRUE),
  grep("_log1p$", names(dlt), value = TRUE)
)

candidates <- setdiff(all_numeric, exclude)
candidates


cat(sprintf("Candidates: %d variables\n\n", length(candidates)))


###############################################################################
# 2. VISUAL CHECK - HISTOGRAMS
###############################################################################

cat("Creating histograms for visual inspection...\n")
cat("Look at Plots pane - use arrows to navigate\n\n")

for (v in candidates) {
  
  p <- ggplot(dlt, aes(x = .data[[v]])) +
    geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
    labs(
      title = paste0("Distribution: ", v),
      subtitle = paste0("Min: ", round(min(dlt[[v]], na.rm = TRUE), 2),
                        " | Max: ", round(max(dlt[[v]], na.rm = TRUE), 2)),
      x = v,
      y = "Frequency"
    ) +
    theme_classic(base_size = 14)
  
  print(p)
}

cat("Review the histograms in Plots pane\n")
cat("Note which variables look RIGHT-SKEWED (long tail to the right)\n\n")

###############################################################################
# 3. CALCULATE BOWLEY SKEWNESS
###############################################################################

bowley_skew <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 20) return(NA)
  q <- quantile(x, c(0.25, 0.5, 0.75))
  iqr <- q[3] - q[1]
  if (iqr == 0) return(NA)
  (q[3] + q[1] - 2*q[2]) / iqr
}

results <- rbindlist(lapply(candidates, function(v) {
  x <- dlt[[v]]
  data.table(
    variable = v,
    n_unique = length(unique(x[!is.na(x)])),
    min = round(min(x, na.rm = TRUE), 2),
    bowley = round(bowley_skew(x), 3)
  )
}))

results <- results[order(-bowley)]

cat("Bowley Skewness (sorted highest to lowest):\n\n")
print(results[, .(variable, bowley)])

cat("\nInterpretation:\n")
cat("  > 0.4  : Highly right-skewed\n")
cat("  0.2-0.4: Moderately right-skewed\n")
cat("  < 0.2  : Not very skewed\n")
cat("  < 0    : LEFT-skewed (don't transform)\n\n")

###############################################################################
# 4. DECIDE & TRANSFORM
###############################################################################

# Automatic: Bowley >= 0.2, positive values
results[, transform := !is.na(bowley) & bowley >= 0.2 & min >= 0 & n_unique >= 10]

vars_to_transform <- results[transform == TRUE, variable]

cat(sprintf("Transforming %d variables (Bowley >= 0.2):\n\n", length(vars_to_transform)))

for (v in vars_to_transform) {
  new_var <- paste0(v, "_log1p")
  dlt[, (new_var) := log1p(get(v))]
  
  bowley_after <- bowley_skew(dlt[[new_var]])
  
  # NEW CONDITION: drop transform if skew flips negative
  
  if (!is.na(bowley_after) && bowley_after < 0) {
    
    dlt[, (new_var) := NULL]                 # remove transformed column
    
    results[variable == v, `:=`(
      
      transform = FALSE,
      
      bowley_after = round(bowley_after, 3),
      
      reason = "log over-corrects (Bowley < 0)"
      
    )]
    
    
    
    cat(sprintf("  ✗ %s NOT transformed (Bowley: %.3f → %.3f)\n",
                
                v,
                
                results[variable == v, bowley],
                
                bowley_after))
    
  } else {
    
    results[variable == v, bowley_after := round(bowley_after, 3)]
    
    
    
    cat(sprintf("  ✓ %s → %s (Bowley: %.3f → %.3f)\n",
                
                v, new_var,
                
                results[variable == v, bowley],
                
                bowley_after))
    
  }
  
}
 

###############################################################################
# 5. VERIFY - BEFORE/AFTER PLOTS
###############################################################################

if (length(vars_to_transform) > 0) {
  
  cat("\nCreating before/after comparison plots...\n")
  cat("Look at Plots pane - use arrows to navigate\n\n")
  
  for (v in vars_to_transform) {
    
    log_var <- paste0(v, "_log1p")
    # If transformed column was not created (or was removed), skip safely
    if (!log_var %in% names(dlt)) {
      cat(sprintf("Skipping %s (missing %s)\n", v, log_var))
      next
    }
    bowley_orig <- results[variable == v, bowley]
    bowley_log <- results[variable == v, bowley_after]
    
    plot_data <- data.table(
      Original = dlt[[v]],
      Transformed = dlt[[log_var]]
    )
    
    # Boxplots
    p1 <- ggplot(plot_data, aes(y = Original)) +
      geom_boxplot(fill = "steelblue", alpha = 0.7) +
      labs(title = paste0(v, " (Original)"),
           subtitle = sprintf("Bowley: %.3f", bowley_orig)) +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    p2 <- ggplot(plot_data, aes(y = Transformed)) +
      geom_boxplot(fill = "forestgreen", alpha = 0.7) +
      labs(title = paste0(log_var, " (Transformed)"),
           subtitle = sprintf("Bowley: %.3f", bowley_log)) +
      theme_classic(base_size = 12) +
      theme(axis.text.x = element_blank(), axis.ticks.x = element_blank())
    
    # Histograms
    p3 <- ggplot(plot_data, aes(x = Original)) +
      geom_histogram(bins = 50, fill = "steelblue", color = "white", alpha = 0.7) +
      labs(x = v) +
      theme_classic(base_size = 12)
    
    p4 <- ggplot(plot_data, aes(x = Transformed)) +
      geom_histogram(bins = 50, fill = "forestgreen", color = "white", alpha = 0.7) +
      labs(x = log_var) +
      theme_classic(base_size = 12)
    
    # Display
    print(grid.arrange(p1, p2, p3, p4, ncol = 2,
                       top = paste0("Transformation: ", v)))
  }
}

###############################################################################
# Age bins (SAP) - 4 cohorts: <=50, 51-65, 66-75, >75
###############################################################################

# Ensure Age is numeric
dlt[, Age := as.numeric(Age)]

# 1) Create age_group per SAP
# right = TRUE means intervals are: (-Inf,50], (50,65], (65,75], (75,Inf]
dlt[, age_group := cut(
  Age,
  breaks = c(-Inf, 50, 65, 75, Inf),
  labels = c("<=50", "51-65", "66-75", ">75"),
  right  = TRUE
)]

# 2) Ensure consistent ordered factor (same order everywhere)
dlt[, age_group := factor(
  age_group,
  levels  = c("<=50", "51-65", "66-75", ">75"),
  ordered = TRUE
)]

# 3) Descriptive label version (nice for Table 1)
dlt[, age_group_desc := factor(
  age_group,
  levels  = c("<=50", "51-65", "66-75", ">75"),
  labels  = c("<=50 (young adults)",
              "51-65 (middle-aged)",
              "66-75 (older adults)",
              ">75 (elderly)"),
  ordered = TRUE
)]

# 4) Sanity check: missing Age => missing age_group (should match counts)
stopifnot(sum(is.na(dlt$Age)) == sum(is.na(dlt$age_group)))

# Quick QC
print(table(dlt$age_group, useNA = "ifany"))

# Optional: cross-tab with a variable if present (edit name to your column)
# Example: ICD_status or Status_FIS etc.
if ("Status_FIS" %in% names(dlt)) {
  print(table(dlt$age_group, dlt$Status_FIS, useNA = "ifany"))
}

saveRDS(dlt, "T:/FINAL ICD COHORT/Transformed_data1.rds")

