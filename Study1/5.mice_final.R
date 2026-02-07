install.packages("tidyverse")
library(tidyverse)
library(dplyr)
library(survival)
library(data.table)
install.packages("mice")
library(mice)
install.packages("VIM")
library(VIM)
install.packages("naniar")
library(naniar)
install.packages("openxlsx")
library(openxlsx)
library(readxl)
install.packages("gt")
library(gt)
install.packages("ggplot2")
library(ggplot2)


set.seed(123)

d <- readRDS("T:/FINAL ICD COHORT/Transformed_data1.rds")
setDT(d)

cat(sprintf("Original dataset: %s observations, %s variables\n\n", 
            format(nrow(d), big.mark = ","), ncol(d)))

#==============================================================================
# STEP 1) TRAIN/TEST SPLIT BY UNIQUE ID (NO LEAKAGE)
#==============================================================================
cat("STEP 1) Train/test split...\n")

set.seed(123)
ids <- unique(d$ID)
train_ids <- sample(ids, size = floor(0.7 * length(ids)), replace = FALSE)
test_ids  <- setdiff(ids, train_ids)

train_data <- d[ID %in% train_ids]
test_data  <- d[ID %in% test_ids]
# --- Explicit copies (to avoid accidental by-reference issues) ---
train_data_copy <- copy(train_data)
test_data_copy  <- copy(test_data)

saveRDS(train_data_copy, "T:/FINAL ICD COHORT/train_raw_id_split1.rds")

saveRDS(test_data_copy,  "T:/FINAL ICD COHORT/test_raw_id_split1.rds")

cat(sprintf("  Train: %d rows (%d IDs)\n", nrow(train_data), length(train_ids)))
cat(sprintf("  Test : %d rows (%d IDs)\n\n", nrow(test_data), length(test_ids)))


#==============================================================================
# STEP 2) DEFINE VARIABLE GROUPS (EDIT HERE ONLY IF NEEDED)
#==============================================================================
cat("STEP 2) Defining variable groups...\n")

id_vars       <- c("ID", "DB")
outcome_vars  <- c("Survival_time", "Status","Status_death","Time_death_days","t_followup_days")          # kept as predictors; NOT imputed
exposure_vars <- c("Status_FIS", "Time_FIS_days")      # NOT imputed; NOT predictors

time_zero_vars <- c("time_zero", "time_zero_ym", "time_index_MI")
time_zero_vars <- time_zero_vars[time_zero_vars %in% names(d)]

log_vars <- grep("_log1p$", names(d), value = TRUE)
original_of_log <- sub("_log1p$", "", log_vars)

bin_vars <- grep("^bin_", names(d), value = TRUE)

all_numeric <- names(d)[sapply(d, function(x) is.numeric(x) || is.integer(x))]
continuous_vars <- setdiff(
  all_numeric,
  c(id_vars, outcome_vars, exposure_vars, time_zero_vars, bin_vars, original_of_log)
)

cat_vars <- names(d)[sapply(d, function(x) is.factor(x) || is.character(x))]

#==============================================================================
# STEP 3) BUILD IMPUTATION DATASETS (SAME COLUMNS IN TRAIN + TEST)
#==============================================================================
cat("STEP 3) Building train_imp_data/test_imp_data...\n")

vars_for_mice <- unique(c(
  id_vars, outcome_vars, exposure_vars, time_zero_vars,
  continuous_vars, cat_vars,
  original_of_log, log_vars,
  bin_vars
))

vars_for_mice <- intersect(vars_for_mice, names(train_data))
vars_for_mice <- intersect(vars_for_mice, names(test_data))

train_imp_data <- as.data.frame(train_data[, ..vars_for_mice])
test_imp_data  <- as.data.frame(test_data[,  ..vars_for_mice])

cat(sprintf("  Variables used: %d\n\n", length(vars_for_mice)))

#==============================================================================

# STEP 3.5) ASSESS MISSINGNESS (TRAIN + TEST) + SAVE GGPLOT FIGURES

#==============================================================================

cat("STEP 3.5) Assess missingness (train + test)...\n")



DIR_FIG <- "T:/Imputed_data/Figures"

dir.create(DIR_FIG, recursive = TRUE, showWarnings = FALSE)

miss_table <- function(df, label){
  data.table(
    cohort      = label,
    variable    = names(df),
    n_missing   = sapply(df, function(x) sum(is.na(x))),
    pct_missing = round(sapply(df, function(x) mean(is.na(x))) * 100, 2)
  )[order(-pct_missing)]
}





plot_missingness <- function(miss_dt,
                             
                             label,
                             
                             top_n = 40,
                             
                             out_dir = DIR_FIG,
                             
                             width = 10,
                             
                             height = 9,
                             
                             dpi = 300,
                             
                             use_percent = TRUE,       # TRUE -> % missing, FALSE -> count missing
                             
                             wrap_width = 25){
  
  
  
  stopifnot(all(c("variable", "n_missing", "pct_missing") %in% names(miss_dt)))
  
  
  
  dd <- data.table::copy(miss_dt)
  
  
  
  # Keep only variables with missingness
  
  dd <- dd[pct_missing > 0]
  
  
  
  # Choose y-axis metric
  
  if (use_percent) {
    
    dd[, value := pct_missing]
    
    xlab <- "% missing"
    
    xlim <- c(0, 100)
    
    xbreaks <- seq(0, 100, 20)
    
  } else {
    
    dd[, value := n_missing]
    
    xlab <- "Number of missing values"
    
    xlim <- c(0, max(dd$value, na.rm = TRUE) * 1.05)
    
    xbreaks <- waiver()
    
  }
  
  
  
  # Take top N by missingness (value)
  
  dd <- dd[order(-value)]
  
  dd <- head(dd, top_n)
  
  
  
  # Pretty labels (replace underscores + wrap)
  
  lab <- gsub("_", " ", dd$variable)
  
  if (requireNamespace("stringr", quietly = TRUE)) {
    
    lab <- stringr::str_wrap(lab, width = wrap_width)
    
  }
  
  dd[, variable_lab := lab]
  
  
  
  # Lollipop plot (layout like your screenshot)
  
  p <- ggplot2::ggplot(dd, ggplot2::aes(y = stats::reorder(variable_lab, value), x = value)) +
    
    ggplot2::geom_segment(
      
      ggplot2::aes(x = 0, xend = value, yend = variable_lab),
      
      linewidth = 0.35,
      
      colour = "#1F4EAA"
      
    ) +
    
    ggplot2::geom_point(
      
      size = 1.8,
      
      colour = "#1F4EAA"
      
    ) +
    
    ggplot2::labs(
      
      title = paste0("Missingness by variable (", label, ")"),
      
      x = xlab,
      
      y = NULL
      
    ) +
    
    ggplot2::theme_minimal(base_size = 12) +
    
    ggplot2::theme(
      
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      
      axis.text.y = ggplot2::element_text(size = 9),
      
      panel.grid.minor = ggplot2::element_blank(),
      
      panel.grid.major.y = ggplot2::element_blank(),
      
      panel.grid.major.x = ggplot2::element_line(linewidth = 0.25, colour = "grey88"),
      
      plot.margin = ggplot2::margin(10, 25, 10, 10)
      
    ) +
    
    ggplot2::scale_x_continuous(limits = xlim, breaks = xbreaks, expand = ggplot2::expansion(mult = c(0, 0.02)))
  
  
  
  # Save
  
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
  
  ggplot2::ggsave(
    
    filename = file.path(out_dir, paste0("missingness_lollipop_top", top_n, "_", label, ".png")),
    
    plot = p,
    
    width = width,
    
    height = height,
    
    dpi = dpi
    
  )
  outfile <- file.path(out_dir,
                       paste0("missingness_lollipop_top", top_n, "_", label, ".png"))
  
  ggplot2::ggsave(outfile, p, width = width, height = height, dpi = dpi)
  
  message("Saved missingness plot to: ", outfile)
  
  
  return(p)
  
}





miss_train <- miss_table(train_imp_data, "TRAIN")

miss_test  <- miss_table(test_imp_data,  "TEST")



# Print top missing vars (quick console view)

cat("\nTop missing variables (TRAIN):\n")

print(miss_train[1:min(.N, 20)], row.names = FALSE)

cat("\nTop missing variables (TEST):\n")

print(miss_test[1:min(.N, 20)], row.names = FALSE)



# Save plots

plot_missingness(miss_train, "TRAIN")

plot_missingness(miss_test,  "TEST")



# Optional: save missingness tables

DIR_LOG <- "T:/Imputed_data/Logs"

dir.create(DIR_LOG, recursive = TRUE, showWarnings = FALSE)

fwrite(miss_train, file.path(DIR_LOG, "missingness_train1.csv"))

fwrite(miss_test,  file.path(DIR_LOG, "missingness_test1.csv"))



cat(sprintf("\n  Saved plots to: %s\n", DIR_FIG))

cat(sprintf("  Saved tables to: %s\n\n", DIR_LOG))



#==============================================================================
# STEP 4) DROP VARIABLES WITH >80% MISSING (TRAIN-BASED; APPLY TO BOTH)
#==============================================================================
cat("STEP 4) Dropping vars with >80% missing (train-based)...\n")
# Variables that must NEVER be dropped
protected_vars <- c(
  "Status_death",
  "Time_death_days",
  "t_followup_years"
)

# Calculate missingness on training data
miss_prop <- colMeans(is.na(train_imp_data))
miss_prop
drop_80 <- names(which(colMeans(is.na(train_imp_data)) > 0.80))

if (length(drop_80)) {
  cat(sprintf("  Dropping %d vars: %s\n\n", length(drop_80), paste(drop_80, collapse = ", ")))
  train_imp_data <- train_imp_data[, !names(train_imp_data) %in% drop_80, drop = FALSE]
  test_imp_data  <- test_imp_data[,  !names(test_imp_data)  %in% drop_80, drop = FALSE]
} else {
  cat("  None dropped\n\n")
}


#==============================================================================
# STEP 5) REMOVE DUPLICATE BINARY VARIABLES (KEEP bin_* DROP BASE VAR)
#==============================================================================
cat("STEP 5) Removing duplicate base vars when bin_* exists...\n")

names(train_imp_data) <- trimws(names(train_imp_data))
names(test_imp_data)  <- trimws(names(test_imp_data))

bin_vars_dup  <- grep("\\bbin_", names(train_imp_data), value = TRUE, ignore.case = TRUE)
base_vars_dup <- sub("\\bbin_", "", bin_vars_dup, ignore.case = TRUE)

drop_dups <- names(train_imp_data)[tolower(names(train_imp_data)) %in% tolower(base_vars_dup)]

if (length(drop_dups)) {
  cat(sprintf("  Dropping %d base vars: %s\n\n", length(drop_dups), paste(drop_dups, collapse = ", ")))
  train_imp_data <- train_imp_data[, !names(train_imp_data) %in% drop_dups, drop = FALSE]
  test_imp_data  <- test_imp_data[,  !names(test_imp_data)  %in% drop_dups, drop = FALSE]
} else {
  cat("  No duplicates found\n\n")
}

# Refresh lists after column drops
log_vars <- grep("_log1p$", names(train_imp_data), value = TRUE)
original_of_log <- sub("_log1p$", "", log_vars)
time_zero_vars <- time_zero_vars[time_zero_vars %in% names(train_imp_data)]

#==============================================================================
# STEP 6) OPTIONAL DIAGNOSTICS: Little's MCAR + MAR quick screen (TRAIN + TEST)
#==============================================================================
has_missing <- function(x) anyNA(x)
missing_rate <- function(x) mean(is.na(x))

safe_glm_p <- function(y, x){
  dat <- data.frame(y = y, x = x)
  dat <- dat[complete.cases(dat), , drop = FALSE]
  if (nrow(dat) < 50) return(NA_real_)
  if (length(unique(dat$y)) < 2) return(NA_real_)
  if (is.factor(dat$x) && nlevels(droplevels(dat$x)) < 2) return(NA_real_)
  fit <- try(glm(y ~ x, data = dat, family = binomial()), silent = TRUE)
  if (inherits(fit, "try-error")) return(NA_real_)
  sm <- summary(fit)$coefficients
  if (nrow(sm) < 2) return(NA_real_)
  sm[2, 4]
}

run_mcar_mar <- function(df, label){
  
  cat(sprintf("=== Diagnostics: %s ===\n", label))
  
  # Little's MCAR test (numeric only)
  if (requireNamespace("naniar", quietly = TRUE)) {
    
    numeric_for_test <- unique(c(
      intersect(continuous_vars, names(df)),
      intersect(log_vars,        names(df))
    ))
    
    mcar_data <- df[, numeric_for_test, drop = FALSE]
    
    cat("Little's MCAR test (numeric only):\n")
    if (ncol(mcar_data) >= 2) {
      mcar <- try(naniar::mcar_test(mcar_data), silent = TRUE)
      if (!inherits(mcar, "try-error")) print(mcar) else cat("  Test could not be performed\n")
    } else {
      cat("  Not enough numeric variables for Little's test\n")
    }
    
  } else {
    cat("Little's MCAR test skipped (install 'naniar' to run it)\n")
  }
  
  # MAR quick screen
  cat("\nMAR quick screen (missingness ~ predictors):\n")
  
  vars_with_missing <- names(df)[vapply(df, has_missing, logical(1))]
  predictors <- names(df)[vapply(df, function(x) missing_rate(x) < 0.5, logical(1))]
  
  excluded_preds <- c(id_vars, exposure_vars, outcome_vars, time_zero_vars)
  predictors <- setdiff(predictors, excluded_preds)
  
  if (length(vars_with_missing) == 0 || length(predictors) == 0) {
    cat("  No variables to test\n\n")
    return(invisible(NULL))
  }
  
  pvals <- c()
  for (v in vars_with_missing) {
    y <- as.integer(is.na(df[[v]]))
    if (length(unique(y)) < 2) next
    for (p in predictors) {
      if (p == v) next
      pv <- safe_glm_p(y, df[[p]])
      if (!is.na(pv)) pvals <- c(pvals, pv)
    }
  }
  
  if (length(pvals) == 0) {
    cat("  No regressions could be performed\n\n")
  } else {
    n_tests <- length(pvals)
    n_sig   <- sum(pvals < 0.05)
    pct_sig <- round(100 * n_sig / n_tests, 1)
    
    cat(sprintf("  Variables with missing: %d\n", length(vars_with_missing)))
    cat(sprintf("  Predictors tested     : %d\n", length(predictors)))
    cat(sprintf("  Regressions performed : %d\n", n_tests))
    cat(sprintf("  Significant assoc.    : %d (%.1f%%)\n", n_sig, pct_sig))
    cat(if (n_sig > 0)
      "  → Missingness associated with observed variables (MAR plausible)\n\n"
      else
        "  → No associations found (MCAR possible, but not proven)\n\n"
    )
  }
  
  invisible(NULL)
}

cat("STEP 6) Diagnostics (optional)...\n")
run_mcar_mar(train_imp_data, "TRAIN")
run_mcar_mar(test_imp_data,  "TEST")
data.frame(
  var = c("bin_anti_diabetic","bin_anti_diabetic_insulin"),
  has_na = sapply(train_imp_data[, c("bin_anti_diabetic","bin_anti_diabetic_insulin")],
                  function(x) anyNA(x)),
  n_levels = sapply(train_imp_data[, c("bin_anti_diabetic","bin_anti_diabetic_insulin")],
                    function(x) length(unique(na.omit(x)))),
  class = sapply(train_imp_data[, c("bin_anti_diabetic","bin_anti_diabetic_insulin")],
                 class)
)
#==============================================================================
# STEP 7) MICE HELPERS + RUN TRAIN/TEST
#   (kept after the data steps, as you prefer)
#==============================================================================
char_to_factor <- function(d){
  for (v in names(d)) if (is.character(d[[v]])) d[[v]] <- factor(d[[v]])
  d
}
to_num <- function(x) suppressWarnings(as.numeric(as.character(x)))

build_mice_spec <- function(dat,
                            id_vars,
                            outcome_vars,
                            exposure_vars,
                            time_zero_vars,
                            log_vars,
                            original_of_log){
  
  dat <- char_to_factor(as.data.frame(dat))
  
  ini  <- mice(dat, maxit = 0, printFlag = FALSE)
  meth <- ini$method
  pred <- ini$predictorMatrix
  diag(pred) <- 0
  
  no_impute <- intersect(c(id_vars, outcome_vars, exposure_vars, time_zero_vars), names(dat))
  meth[no_impute] <- ""
  
  pred[, intersect(id_vars, colnames(pred))] <- 0
  pred[, intersect(time_zero_vars, colnames(pred))] <- 0
  pred[, intersect(outcome_vars, colnames(pred))] <- 1
  pred[, intersect(exposure_vars, colnames(pred))] <- 0
  
  log_vars  <- intersect(log_vars, names(dat))
  base_vars <- intersect(original_of_log, names(dat))
  passive_skipped <- character(0)
  
  for (lv in log_vars) {
    base <- base_vars[paste0(base_vars, "_log1p") == lv]
    if (length(base) == 1) {
      dat[[base]] <- to_num(dat[[base]])
      dat[[lv]]   <- to_num(dat[[lv]])
      meth[lv] <- sprintf("~I(log1p(%s))", base)
      pred[lv, ] <- 0
      pred[lv, base] <- 1
    } else {
      meth[lv] <- ""
      pred[, lv] <- 0
      pred[lv, ] <- 0
      passive_skipped <- c(passive_skipped, lv)
    }
  }
  
  for (v in names(dat)) {
    if (meth[v] == "" || grepl("^~I\\(", meth[v])) next
    x <- dat[[v]]
    
    if (is.factor(x)) {
      nl <- nlevels(droplevels(x))
      if (nl < 2) {
        meth[v] <- ""
        pred[, v] <- 0; pred[v, ] <- 0
      } else if (nl == 2) {
        meth[v] <- "logreg"
      } else {
        meth[v] <- "polyreg"
      }
      next
    }
    
    if (is.numeric(x) || is.integer(x)) {
      ux <- unique(x[!is.na(x)])
      if (length(ux) == 2 && all(ux %in% c(0, 1))) {
        dat[[v]] <- factor(x, levels = c(0, 1))
        meth[v] <- "logreg"
      } else {
        meth[v] <- "pmm"
      }
      next
    }
    
    meth[v] <- "pmm"
  }
  
  meth[names(dat)[!sapply(dat, anyNA)]] <- ""
  targets <- names(meth)[meth != "" & !grepl("^~I\\(", meth)]
  
  list(dat = dat, meth = meth, pred = pred, targets = targets,
       passive_skipped = passive_skipped, no_impute = no_impute)
}

mice_guardrail <- function(spec, seed = 123, tries_max = 3){
  
  dat <- spec$dat
  meth <- spec$meth
  pred <- spec$pred
  targets <- spec$targets
  
  tries <- 0
  removed_total <- 0
  
  repeat {
    set.seed(seed)
    imp_test <- mice(dat, m = 2, maxit = 5,
                     method = meth, predictorMatrix = pred,
                     printFlag = FALSE)
    cmpl <- complete(imp_test, 1)
    leftover <- if (length(targets)) colSums(is.na(cmpl[, targets, drop = FALSE])) else integer(0)
    
    if (!length(leftover) || all(leftover == 0) || tries >= tries_max) break
    
    tries <- tries + 1
    idx_bad <- which(rowSums(is.na(cmpl[, targets, drop = FALSE])) > 0)
    
    cat(sprintf("  Guardrail: removing %d non-imputable rows (try %d)\n",
                length(idx_bad), tries))
    
    removed_total <- removed_total + length(idx_bad)
    dat <- dat[-idx_bad, , drop = FALSE]
    
    spec2 <- build_mice_spec(dat, id_vars, outcome_vars, exposure_vars, time_zero_vars, log_vars, original_of_log)
    dat <- spec2$dat; meth <- spec2$meth; pred <- spec2$pred; targets <- spec2$targets
  }
  
  list(dat = dat, meth = meth, pred = pred, targets = targets, removed_total = removed_total)
}

run_mice <- function(df, label, seed = 123){
  
  cat(sprintf("\nSTEP 7) Running MICE for %s...\n", label))
  
  spec <- build_mice_spec(df, id_vars, outcome_vars, exposure_vars, time_zero_vars, log_vars, original_of_log)
  
  if (length(spec$passive_skipped)) {
    cat("  Passive *_log1p skipped (missing base): ",
        paste(spec$passive_skipped, collapse = ", "), "\n", sep = "")
  }
  
  cat(sprintf("  Targets (non-passive imputed vars): %d\n", length(spec$targets)))
  
  #grd <- mice_guardrail(spec, seed = seed)
  
  # ========== FORCE IMPUTATION FOR SPECIFIC BINARY VARIABLES ==========
  # Force bin_anti_diabetic
  # Fix: bin_anti_diabetic is perfectly collinear with bin_diabetes (r=1.0 in observed data)
  # Cross-tab shows 1728 cases where both=No and 54 cases where both=Yes, zero disagreements
  # MICE cannot impute one variable using itself, so we directly copy bin_diabetes values
  if ("bin_anti_diabetic" %in% names(spec$dat) && "bin_diabetes" %in% names(spec$dat)) {
    missing_idx <- is.na(spec$dat$bin_anti_diabetic)
    spec$dat$bin_anti_diabetic[missing_idx] <- spec$dat$bin_diabetes[missing_idx]
    spec$meth["bin_anti_diabetic"] <- ""
    cat(sprintf("  [%s] ✓ Filled bin_anti_diabetic from bin_diabetes (%d values)\n", 
                label, sum(missing_idx)))
  }
  
  
  
  
  # Force bin_anti_diabetic_insulin (if it has 2+ levels, otherwise skip)
  if ("bin_anti_diabetic_insulin" %in% names(spec$dat) && anyNA(spec$dat$bin_anti_diabetic_insulin)) {
    unique_vals <- unique(spec$dat$bin_anti_diabetic_insulin[!is.na(spec$dat$bin_anti_diabetic_insulin)])
    
    if (length(unique_vals) >= 2) {
      spec$meth["bin_anti_diabetic_insulin"] <- "logreg"
      spec$pred["bin_anti_diabetic_insulin", ] <- 1
      spec$pred["bin_anti_diabetic_insulin", c(id_vars, outcome_vars, exposure_vars, time_zero_vars, "bin_anti_diabetic_insulin")] <- 0
      cat(sprintf("  [%s] Forcing bin_anti_diabetic_insulin imputation with logreg\n", label))
    } else {
      # Single level - just fill with the constant
      spec$dat$bin_anti_diabetic_insulin[is.na(spec$dat$bin_anti_diabetic_insulin)] <- unique_vals[1]
      spec$meth["bin_anti_diabetic_insulin"] <- ""
      cat(sprintf("  [%s] Pre-filled bin_anti_diabetic_insulin with constant %s\n", label, unique_vals[1]))
    }
  }
  # ====================================================================
  grd <- mice_guardrail(spec, seed = seed)
  
  set.seed(seed)
  imp <- mice(grd$dat, m = 20, maxit = 10,
              method = grd$meth, predictorMatrix = grd$pred,
              printFlag = TRUE)
  
  d1 <- complete(imp, 1)
  if (length(grd$targets)) {
    stopifnot(all(colSums(is.na(d1[, grd$targets, drop = FALSE])) == 0))
  }
  
  cat(sprintf("  ✓ %s done. Rows removed by guardrail: %d\n", label, grd$removed_total))
  list(imp = imp, removed_total = grd$removed_total, targets = grd$targets)
}

# Run MICE (TRAIN + TEST)
res_train <- run_mice(train_imp_data, "TRAIN", seed = 123)
res_test  <- run_mice(test_imp_data,  "TEST",  seed = 123)

imp_train <- res_train$imp
imp_test  <- res_test$imp
imp_train$loggedEvents
colSums(is.na(complete(imp_train)))
colSums(is.na(complete(imp_test)))
# TRAIN
train_before <- train_imp_data
train_after  <- complete(imp_train, 1)

# TEST
test_before  <- test_imp_data
test_after   <- complete(imp_test, 1)  

cat("\nTEST – remaining NAs:\n")
print(colSums(is.na(test_after))[colSums(is.na(test_after)) > 0])

cat("\n=== NA CHECK AFTER IMPUTATION ===\n")

cat("TRAIN – remaining NAs:\n")
print(colSums(is.na(train_after))[colSums(is.na(train_after)) > 0])

# Save outputs
DIR_MI <- "T:/Imputed_data"
dir.create(DIR_MI, recursive = TRUE, showWarnings = FALSE)

saveRDS(imp_train, file.path(DIR_MI, "mice_train_object3.rds"))
saveRDS(imp_test,  file.path(DIR_MI, "mice_test_object3.rds"))

fwrite(as.data.table(complete(imp_train, "long", include = TRUE)),
       file.path(DIR_MI, "train_imputed_long1.csv"))
fwrite(as.data.table(complete(imp_test,  "long", include = TRUE)),
       file.path(DIR_MI, "test_imputed_long1.csv"))

cat("\nSaved:\n",
    "  - mice_train_object.rds\n  - mice_test_object.rds\n",
    "  - train_imputed_long.csv\n  - test_imputed_long.csv\n", sep = "")


#==============================================================================
# STEP 11) MICE AUDIT LOG
#==============================================================================
cat("\n=== MICE AUDIT LOG ===\n")

create_audit_log <- function(imp_obj, original_data, label, seed) {
  
  # Variables not imputed (structural)
  not_imputed_vars <- names(imp_obj$method)[imp_obj$method == ""]
  
  # Variables imputed
  imputed_vars <- names(imp_obj$method)[imp_obj$method != "" & !grepl("^~I\\(", imp_obj$method)]
  
  # Missingness BEFORE imputation
  miss_before <- colSums(is.na(original_data))
  
  # Missingness AFTER imputation (use first completed dataset)
  comp1 <- complete(imp_obj, 1)
  miss_after <- colSums(is.na(comp1))
  
  # Create audit log
  audit_log <- data.table(
    variable        = names(miss_before),
    missing_before  = miss_before,
    missing_after   = miss_after,
    imputed         = names(miss_before) %in% imputed_vars,
    method          = imp_obj$method[names(miss_before)]
  )
  
  # Add global metadata
  audit_log[, `:=`(
    cohort          = label,
    n_individuals   = nrow(original_data),
    n_imputations   = imp_obj$m,
    max_iterations  = imp_obj$iteration,
    seed            = seed
  )]
  
  # Order by imputed status and missing count
  setorder(audit_log, -imputed, -missing_before)
  
  return(audit_log)
}

# Create audit logs for TRAIN and TEST
audit_train <- create_audit_log(imp_train, train_imp_data, "TRAIN", seed = 123)
audit_test  <- create_audit_log(imp_test,  test_imp_data,  "TEST",  seed = 123)

# Print summaries
for (audit in list(list(log = audit_train, label = "TRAIN"), 
                   list(log = audit_test,  label = "TEST"))) {
  
  cat(sprintf("\n--- %s ---\n", audit$label))
  cat(sprintf("Variables imputed      : %d\n", sum(audit$log$imputed)))
  cat(sprintf("Variables NOT imputed  : %d\n", sum(!audit$log$imputed)))
  cat(sprintf("Total missing BEFORE   : %s\n", format(sum(audit$log$missing_before), big.mark=",")))
  cat(sprintf("Total missing AFTER    : %s\n", format(sum(audit$log$missing_after), big.mark=",")))
  cat(sprintf("Reduction              : %s (%.1f%%)\n", 
              format(sum(audit$log$missing_before) - sum(audit$log$missing_after), big.mark=","),
              100 * (sum(audit$log$missing_before) - sum(audit$log$missing_after)) / sum(audit$log$missing_before)))
}

# Save audit logs
fwrite(audit_train, file.path(DIR_MI, "mice_audit_log_TRAIN1.csv"))
fwrite(audit_test,  file.path(DIR_MI, "mice_audit_log_TEST2.csv"))

# Combine for a full audit
audit_combined <- rbindlist(list(audit_train, audit_test))
fwrite(audit_combined, file.path(DIR_MI, "mice_audit_log_COMBINED1.csv"))

cat("\nMICE audit logs saved:\n")
cat("  - mice_audit_log_TRAIN.csv\n")
cat("  - mice_audit_log_TEST.csv\n")
cat("  - mice_audit_log_COMBINED.csv\n")
