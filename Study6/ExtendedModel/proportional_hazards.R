
library(survival)
library(dplyr)
library(ggplot2)
library(survminer)
library(patchwork)

# Set working directory
setwd("T:/Dokumente/PROFID/Study6")

# Load objects ----
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


# Schoenfeld residual tests per imputation ----

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


# Log(-log) survival plots for key groups 


dat_all <- imp2$data

# Sex
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

# ICD status
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

#  LVEF grouped (e.g. <35, 35–50, >50)
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

# plot all rows that match a variable name
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


terms_all <- rownames(czph$table)
terms_all <- terms_all[terms_all != "GLOBAL"]  # drop global row

fm1  <- fit_list_cs1$analyses[[1]]
czph <- cox.zph(fm1)

# all terms except GLOBAL
idx <- which(rownames(czph$table) != "GLOBAL")

pdf("PH_Schoenfeld_all_covariates_panels_imp1_extended.pdf",
    width = 8, height = 10)

n_rows <- 3
n_cols <- 2
par(mfrow = c(n_rows, n_cols))

for (i in idx) {
  term_name <- rownames(czph$table)[i]
  plot(czph[i], main = term_name)
 
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



czph <- cox.zph(fm1, transform = "km")

pdf("PH_Schoenfeld_panels_2x2.pdf", width = 8, height = 10)

par(
  mfrow = c(2,2),
  mar = c(4.5, 4.5, 3, 1),
  oma = c(0, 0, 2, 0)
)

## ---- a) BMI spline (custom y-axis) ----
plot(
  czph,
  var = "ns(BMI, knots = K4, Boundary.knots = BND)",
  resid = TRUE,
  se = TRUE,
  ylab = "",
  main = "a) BMI spline"
)

mtext(
  "Scaled Schoenfeld residuals",
  side = 2,
  line = 3
)

## ---- b) Age ----
plot(
  czph,
  var = "Age",
  resid = TRUE,
  se = TRUE,
  ylab = "",
  main = "b) Age"
)

mtext(
  "Scaled Schoenfeld residuals",
  side = 2,
  line = 3
)

## ---- c) LVEF ----
plot(
  czph,
  var = "LVEF",
  resid = TRUE,
  se = TRUE,
  ylab = "",
  main = "c) LVEF"
)

mtext(
  "Scaled Schoenfeld residuals",
  side = 2,
  line = 3
)

## ---- d) eGFR ----
plot(
  czph,
  var = "eGFR",
  resid = TRUE,
  se = TRUE,
  ylab = "",
  main = "d) eGFR"
)

mtext(
  "Scaled Schoenfeld residuals",
  side = 2,
  line = 3
)

## ---- Global PH test ----
mtext(
  paste0(
    "Global proportional hazards test: p = ",
    format.pval(czph$table["GLOBAL","p"], digits = 3, eps = 1e-3)
  ),
  outer = TRUE,
  cex = 0.9
)

dev.off()
