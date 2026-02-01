library(mice)
library(survival)
library(splines)
library(dplyr)

setwd("T:/Dokumente/PROFID/Study6")

# Load imputed object
imp2 <- readRDS("mice_imputed_data_extended.RDS")

# Recompute knots (same as your main extended model code)
q   <- quantile(imp2$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])
K4  <- as.numeric(q[c(2,3,4,5)])

# 90-month horizon variables (match your updated extended model)
horizon <- 90

make_horizon <- function(d) {
  tt <- d$Survival_time
  ss <- d$Status_cs1
  ss <- suppressWarnings(as.integer(as.character(ss)))
  
  d$Survival_time_h <- pmin(tt, horizon)
  d$Status_cs1_h    <- ifelse(!is.na(ss) & ss == 1L & tt <= horizon, 1L, 0L)
  d
}

# UPDATED extended model formula (remove COPD/Cancer; match your current extended model)
form_ext_h <- Surv(Survival_time_h, Status_cs1_h) ~
  ns(BMI, knots = K4, Boundary.knots = BND) +
  Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
  LVEF + eGFR + Haemoglobin +
  ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
  Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
  Stroke_TIA + ICD_status

# Fit across imputations (use imputation 1 for residual plots, like you did before)
fit_ext_h <- with(imp2, {
  d <- make_horizon(data)
  coxph(form_ext_h, data = d, ties = "efron", x = TRUE)
})

fm1 <- fit_ext_h$analyses[[1]]

ph_tbl_list <- lapply(fit_ext_h$analyses, function(fm) as.data.frame(cox.zph(fm, transform="km")$table))

ph_pvals <- Reduce(function(a,b) cbind(a["p"], b["p"]), ph_tbl_list)
colnames(ph_pvals) <- paste0("imp", seq_len(ncol(ph_pvals)))
ph_pvals$term <- rownames(ph_tbl_list[[1]])
ph_pvals <- ph_pvals[, c("term", grep("^imp", names(ph_pvals), value = TRUE))]

viol_rate <- ph_pvals %>%
  filter(term != "GLOBAL") %>%
  rowwise() %>%
  mutate(prop_p_lt_0.05 = mean(c_across(starts_with("imp")) < 0.05, na.rm = TRUE)) %>%
  ungroup() %>%
  arrange(desc(prop_p_lt_0.05))

write.csv(ph_pvals,  "PH_pvalues_per_imputation_extended_90mo.csv", row.names = FALSE)
write.csv(viol_rate, "PH_violation_rate_extended_90mo.csv", row.names = FALSE)


# Schoenfeld test
czph <- cox.zph(fm1, transform = "log")  # keep "km" to match your previous scripts

# Helper: plot all rows matching a pattern (CRUCIAL for spline terms)
plot_terms <- function(zph_obj, pattern, xlim = c(0, 90)) {
  idx <- grep(pattern, rownames(zph_obj$table))
  if (length(idx) == 0) {
    message("No coefficient rows match: ", pattern)
  } else {
    for (i in idx) {
      term_name <- rownames(zph_obj$table)[i]
      plot(zph_obj[i], resid = TRUE, se = TRUE, xlim = xlim, main = term_name,
           xlab = "Time (months)", ylab = "Scaled Schoenfeld residuals")
    }
  }
}

# --- A) Selected plots (BMI spline + key covariates) ---
pdf("PH_selected_covariates_imp1_extended_90mo.pdf", width = 8.5, height = 10)
par(mfrow = c(3,2), mar = c(4,4,2.5,1))

plot_terms(czph, "^ns\\(BMI", xlim = c(0, 90))  # plots all BMI spline basis rows
plot_terms(czph, "^Age",      xlim = c(0, 90))
plot_terms(czph, "^eGFR",     xlim = c(0, 90))
plot_terms(czph, "^Haemoglobin", xlim = c(0, 90))
plot_terms(czph, "^LVEF",     xlim = c(0, 90))  # optional

dev.off()
par(mfrow = c(1,1))

# --- B) All terms panel (multi-page, 6 per page) ---
terms_idx <- which(rownames(czph$table) != "GLOBAL")

pdf("PH_all_terms_imp1_extended_90mo.pdf", width = 8.5, height = 10)
par(mfrow = c(3,2), mar = c(4,4,2.5,1))

for (k in seq_along(terms_idx)) {
  i <- terms_idx[k]
  plot(czph[i], resid = TRUE, se = TRUE, xlim = c(0, 90),
       main = rownames(czph$table)[i],
       xlab = "Time (months)", ylab = "Scaled Schoenfeld residuals")
  
  # start a new page every 6 plots
  if (k %% 6 == 0) {
    par(mfrow = c(3,2))
  }
}
dev.off()
par(mfrow = c(1,1))
