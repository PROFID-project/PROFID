
# PH assumptions – extended imputed model


# Clear environment (optional)
rm(list = ls())

# 1. Setup ----
library(survival)
library(dplyr)
library(ggplot2)
library(survminer)
library(patchwork)

# Set working directory
setwd("T:/Dokumente/PROFID/Study6")

# 2. Load objects ----
# Imputed mids object with extended variables
imp2         <- readRDS("mice_imputed_data_extended.RDS")

# List of Cox models across imputations (extended BMI spline model)
fit_list_cs1 <- readRDS("fit_list_cs1_extended.RDS")


# Recompute knots
q   <- quantile(imp2$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])              # boundary knots 5th & 95th
K4  <- as.numeric(q[c(2,3,4,5)])          # inner knots 10,35,65,90


# Quick sanity check
length(fit_list_cs1$analyses)
fit_list_cs1$analyses[[1]]
names(imp2$data)[1:20]


# 3. Schoenfeld residual tests per imputation ----

# Run cox.zph on each imputed model
ph_tbl_list <- lapply(fit_list_cs1$analyses, function(fm) {
  as.data.frame(cox.zph(fm)$table)
})

# Combine p-values across imputations
# Each table has rows = terms (including "GLOBAL"), columns = rho, chisq, p
ph_pvals <- Reduce(
  function(a, b) cbind(a["p"], b["p"]),
  ph_tbl_list
)

colnames(ph_pvals) <- paste0("imp", seq_len(ncol(ph_pvals)))
ph_pvals$term <- rownames(ph_tbl_list[[1]])
ph_pvals <- ph_pvals[, c("term", grep("^imp", names(ph_pvals), value = TRUE))]

# Proportion of imputations with p < 0.05 (evidence of PH violation)
viol_rate <- ph_pvals %>%
  filter(term != "GLOBAL") %>%  # drop global if you want per-covariate view
  rowwise() %>%
  mutate(
    prop_p_lt_0.05 = mean(c_across(starts_with("imp")) < 0.05, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  arrange(desc(prop_p_lt_0.05))

viol_rate

# Save PH test outputs
write.csv(ph_pvals,
          "PH_pvalues_per_imputation_extended_new.csv",
          row.names = FALSE)

write.csv(viol_rate,
          "PH_violation_rate_extended_new.csv",
          row.names = FALSE)


fm1  <- fit_list_cs1$analyses[[1]]
czph <- cox.zph(fm1)

# Helper to plot all rows that match a variable name
plot_terms <- function(zph_obj, pattern, title = NULL) {
  idx <- grep(pattern, rownames(zph_obj$table))
  if (length(idx) == 0) {
    message("No coefficient rows match: ", pattern)
  } else {
    for (i in idx) {
      plot(zph_obj[i], main = paste0(pattern, " — ", rownames(zph_obj$table)[i]))
    }
  }
}

pdf("PH_selected_covariates_imp1_extended.pdf", width = 8, height = 10)

# BMI (all spline basis terms)
plot_terms(czph, "^ns\\(BMI")

# Continuous variables
plot_terms(czph, "^Age")
plot_terms(czph, "^LVEF")
plot_terms(czph, "^eGFR")

# Binary/factors
plot_terms(czph, "^Sex")
plot_terms(czph, "^ICD_status")

dev.off()


# 5. Log(-log) survival plots for key groups ----
# Use the stacked imputed data (imp2$data) as a single "pseudo-dataset".
# These are descriptive PH checks, not formal tests.

dat_all <- imp2$data

# Example 1: Sex
sf_sex <- survfit(Surv(Survival_time, Status_cs1) ~ Sex, data = dat_all)

p_sex <- ggsurvplot(
  sf_sex,
  fun = "cloglog",
  xlab = "Time",
  ylab = "log(-log(Survival))",
  legend.title = "Sex",
  legend.labs = levels(dat_all$Sex),
  palette = c("#1F77B4", "#FF7F0E")
)$plot +
  ggtitle("log(-log(Survival)) by Sex")

# Example 2: ICD status
sf_icd <- survfit(Surv(Survival_time, Status_cs1) ~ ICD_status, data = dat_all)

p_icd <- ggsurvplot(
  sf_icd,
  fun = "cloglog",
  xlab = "Time",
  ylab = "log(-log(Survival))",
  legend.title = "ICD status",
  legend.labs = levels(dat_all$ICD_status),
  palette = c("#2CA02C", "#D62728")
)$plot +
  ggtitle("log(-log(Survival)) by ICD status")

# Example 3: LVEF grouped (e.g. <35, 35–50, >50)
dat_all <- dat_all %>%
  mutate(
    LVEF_group = cut(
      LVEF,
      breaks = c(-Inf, 35, 50, Inf),
      labels = c("<35%", "35–50%", ">50%")
    )
  )

sf_lvef <- survfit(Surv(Survival_time, Status_cs1) ~ LVEF_group, data = dat_all)

p_lvef <- ggsurvplot(
  sf_lvef,
  fun = "cloglog",
  xlab = "Time",
  ylab = "log(-log(Survival))",
  legend.title = "LVEF group",
  legend.labs = levels(dat_all$LVEF_group),
  palette = c("#9467BD", "#8C564B", "#17BECF")
)$plot +
  ggtitle("log(-log(Survival)) by LVEF group")

# Combine plots into one panel and save
panel_loglog <- (p_sex | p_icd) / p_lvef

ggsave(
  "loglog_PH_checks_extended_imputed.pdf",
  panel_loglog,
  width = 10, height = 8
)




czph <- cox.zph(fm1)

# Helper to plot all rows that match a variable name
plot_terms <- function(zph_obj, pattern, title = NULL) {
  idx <- grep(pattern, rownames(zph_obj$table))
  if (length(idx) == 0) {
    message("No coefficient rows match: ", pattern)
  } else {
    for (i in idx) {
      plot(zph_obj[i], main = paste0(pattern, " — ", rownames(zph_obj$table)[i]))
    }
  }
}

pdf("PH_selected_covariates_imp1_extended.pdf", width = 8, height = 10)

terms_all <- rownames(czph$table)
terms_all <- terms_all[terms_all != "GLOBAL"]  # drop global row

pdf("PH_Schoenfeld_all_covariates_panels_imp1_extended.pdf", width = 8, height = 10)

# choose how many per page
n_per_page <- 6
n_rows     <- 3
n_cols     <- 2

par(mfrow = c(n_rows, n_cols))

for (i in seq_along(terms_all)) {
  term_i <- terms_all[i]
  plot(czph[term_i], main = term_i)
  
  # when we’ve filled a page, start a new one
  if (i %% n_per_page == 0 && i != length(terms_all)) {
    par(mfrow = c(n_rows, n_cols))  # new page, reset layout
  }
}

dev.off()
par(mfrow = c(1, 1))


# BMI (all spline basis terms)
plot_terms(czph, "^ns\\(BMI")

# Continuous variables
plot_terms(czph, "^Age")
plot_terms(czph, "^LVEF")
plot_terms(czph, "^eGFR")

# Binary/factors
plot_terms(czph, "^Sex")
plot_terms(czph, "^ICD_status")

dev.off()

# Variables with high violation rate
head(viol_rate, 10)

# Time-varying coefficient model: Age, LVEF, eGFR
fit_tvc_small <- with(
  imp2,
  coxph(
    Surv(Survival_time, Status_cs1) ~
      ns(BMI, knots = K4, Boundary.knots = BND) +
      Sex + ICD_status +
      Age + LVEF + eGFR +
      tt(Age) + tt(LVEF) + tt(eGFR),
    x = FALSE, y = FALSE,
    tt = function(x,t,...) x*log(t)
  )
)

# Pooled TVC model
pool_fit_tvc <- pool(fit_tvc_small)
summary(pool_fit_tvc)

saveRDS(fit_tvc,      "fit_list_cs1_extended_TVC.RDS")
saveRDS(pool_fit_tvc, "pool_fit_cs1_extended_TVC.RDS")

# C-index per imputation for TVC model
cindex_tvc <- sapply(fit_tvc_small$analyses, function(fm) {
  y  <- fm$y
  lp <- predict(fm, type = "lp")
  sc <- survConcordance(Surv(y[,1], y[,2]) ~ lp)
  c(C = sc$concordance, Var = sc$var)
})

cindex_tvc <- t(cindex_tvc)
C_tvc  <- mean(cindex_tvc[,"C"])
Var_tvc <- mean(cindex_tvc[,"Var"])
SE_tvc  <- sqrt(Var_tvc)
CI_tvc  <- C_tvc + c(-1.96, 1.96)*SE_tvc

cindex_tvc_summary <- data.frame(
  Model  = "Extended + TVC(Age,LVEF,eGFR)",
  C_index = C_tvc,
  SE      = SE_tvc,
  LCL     = CI_tvc[1],
  UCL     = CI_tvc[2]
)

cindex_tvc_summary
write.csv(cindex_tvc_summary,
          "C_index_extended_TVC_cs1.csv",
          row.names = FALSE)

