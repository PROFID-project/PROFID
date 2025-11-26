library(survival)
library(splines)
library(dplyr)

setwd("T:/Dokumente/PROFID/Study6")

## Load the imputed mids object and saved fits
imp2          <- readRDS("mice_imputed_data_extended.RDS")
fit_list_cs1  <- readRDS("fit_list_cs1_extended.RDS")     # per-imputation fits
pool_fit_cs1  <- readRDS("pool_fit_cs1_extended.RDS")     # pooled object
curve_mice    <- read.csv("cox_RCS_BMI_curve_cs1_extended.csv")  # pooled BMI–HR curve
head(curve_mice)

q   <- quantile(imp2$data$BMI, probs = c(.05,.10,.35,.65,.90,.95), na.rm = TRUE)
BND <- as.numeric(q[c(1,6)])              # boundary knots 5th & 95th
K4  <- as.numeric(q[c(2,3,4,5)])          # inner knots 10,35,65,90

# Original non-imputed data
preimp <- read.csv("combined_BMI_outcomefiltered.csv")

# Recode Cancer and COPD_cat as before
preimp_cc <- preimp %>%
  mutate(
    # Cancer: treat missing as "no" (0)
    Cancer = dplyr::case_when(
      Cancer %in% c("Yes", 1L) ~ 1L,
      is.na(Cancer)            ~ 0L,
      TRUE                     ~ 0L
    ),
    # COPD: make 3-level factor (No / Yes / Missing)
    COPD_cat = dplyr::case_when(
      is.na(COPD)             ~ "Missing",
      COPD %in% c("Yes", 1L)  ~ "Yes",
      TRUE                    ~ "No"
    ),
    COPD_cat = factor(COPD_cat, levels = c("No","Yes","Missing")),
    # Cause-specific event indicator: 1 = Status==1, else 0
    Status_cs1 = ifelse(Status == 1, 1L, 0L)
  )

# Variables used in the extended spline model
vars_model <- c(
  "Age","Sex",
  "Diabetes","Hypertension","Smoking","MI_history",
  "LVEF",
  "eGFR","Haemoglobin",
  "ACE_inhibitor_ARB","Beta_blockers","Lipid_lowering",
  "Revascularisation_acute", "Cholesterol", "HDL", "LDL", "Triglycerides",
  "COPD_cat", "Cancer", "Stroke_TIA", "ICD_status",
  "BMI", "Survival_time", "Status_cs1"
)

# Complete-case dataset: only rows with no missing values in model vars
cca_data <- preimp_cc %>%
  dplyr::filter(stats::complete.cases(dplyr::across(all_of(vars_model))))

nrow(preimp_cc); nrow(cca_data)   # check how many you lose

cox_cca <- coxph(
  Surv(Survival_time, Status_cs1) ~
    ns(BMI, knots = K4, Boundary.knots = BND) +
    Age + Sex + Diabetes + Hypertension + Smoking + MI_history +
    LVEF  + eGFR + Haemoglobin +
    ACE_inhibitor_ARB + Beta_blockers + Lipid_lowering +
    Revascularisation_acute + Cholesterol + HDL + LDL + Triglycerides +
    COPD_cat + Cancer + Stroke_TIA + ICD_status,
  data = cca_data,
  x = TRUE, y = TRUE
)

summary(cox_cca)
saveRDS(cox_cca, "cox_RCS_cs1_extended_CCA.RDS")

# BMI grid within same boundaries as MICE analysis
bmi_min <- max(min(cca_data$BMI, na.rm = TRUE), BND[1])
bmi_max <- min(max(cca_data$BMI, na.rm = TRUE), BND[2])
bmi_grid <- seq(bmi_min, bmi_max, length.out = 400)

# BMI grid within same boundaries as MICE analysis
bmi_min <- max(min(cca_data$BMI, na.rm = TRUE), BND[1])
bmi_max <- min(max(cca_data$BMI, na.rm = TRUE), BND[2])
bmi_grid <- seq(bmi_min, bmi_max, length.out = 400)

# Covariate variables (exclude BMI + time + event)
covariate_vars <- setdiff(vars_model, c("BMI", "Survival_time", "Status_cs1"))

# 1) Take the first complete-case row as a template profile
template_row <- cca_data[1, covariate_vars, drop = FALSE]

# 2) Repeat that row for every BMI grid point
template_covs <- template_row[rep(1, length(bmi_grid)), , drop = FALSE]

# 3) Combine BMI grid with the covariate template
newdata_cca <- cbind(BMI = bmi_grid, template_covs)


pred_cca <- predict(cox_cca, newdata = newdata_cca, type = "lp", se.fit = TRUE)
lp       <- pred_cca$fit
se_lp    <- pred_cca$se.fit

# Centre at BMI = 25 (same as your MICE curve)
ref_idx  <- which.min(abs(bmi_grid - 25))
lp_ref   <- lp[ref_idx]

logHR    <- lp - lp_ref
HR       <- exp(logHR)
LCL      <- exp(logHR - 1.96 * se_lp)
UCL      <- exp(logHR + 1.96 * se_lp)

curve_cca <- data.frame(BMI = bmi_grid, HR = HR, LCL = LCL, UCL = UCL)
write.csv(curve_cca, "cox_RCS_BMI_curve_cs1_extended_CCA.csv", row.names = FALSE)

curve_cca  <- read.csv("cox_RCS_BMI_curve_cs1_extended_CCA.csv")
curve_mice <- read.csv("cox_RCS_BMI_curve_cs1_extended.csv")

pdf("BMI_curve_MICE_vs_CCA_extended.pdf", width = 7, height = 5)

plot(curve_mice$BMI, curve_mice$HR, type = "n",
     ylim = range(c(curve_mice$LCL, curve_mice$UCL,
                    curve_cca$LCL, curve_cca$UCL)),
     xlab = "BMI (kg/m²)",
     ylab = "Hazard ratio (ref = 25)",
     main = "Restricted cubic spline: MICE vs Complete-case")

# MICE shaded CI
polygon(c(curve_mice$BMI, rev(curve_mice$BMI)),
        c(curve_mice$LCL, rev(curve_mice$UCL)),
        col = adjustcolor("darkblue", alpha.f = 0.15), border = NA)
lines(curve_mice$BMI, curve_mice$HR, col = "darkblue", lwd = 2)

# CCA shaded CI
polygon(c(curve_cca$BMI, rev(curve_cca$BMI)),
        c(curve_cca$LCL, rev(curve_cca$UCL)),
        col = adjustcolor("firebrick", alpha.f = 0.15), border = NA)
lines(curve_cca$BMI, curve_cca$HR, col = "firebrick", lwd = 2)

abline(h = 1,  lty = 3, col = "gray50")
abline(v = 25, lty = 3, col = "gray80")

legend("topleft",
       legend = c("MICE (main analysis)", "Complete case"),
       col    = c("darkblue", "firebrick"),
       lwd    = 2,
       bty    = "n")

dev.off()

## Full cohort size (pre-imputation)
n_full   <- nrow(preimp_cc)
events_full <- sum(preimp_cc$Status_cs1 == 1, na.rm = TRUE)

## Complete-case cohort
n_cca    <- nrow(cca_data)
events_cca <- sum(cca_data$Status_cs1 == 1, na.rm = TRUE)

## Retention percentages
pct_n_retained     <- 100 * n_cca / n_full
pct_events_retained <- 100 * events_cca / events_full

sample_summary <- data.frame(
  cohort          = c("Full (pre-imputation)", "Complete case"),
  N               = c(n_full, n_cca),
  Events          = c(events_full, events_cca),
  Percent_of_full = c(100, pct_n_retained),
  Events_percent_of_full = c(100, pct_events_retained)
)

sample_summary

library(survival)

# Linear predictor for each subject
lp_cca <- predict(cox_cca, type = "lp")

# Concordance calculation
sc_cca <- survConcordance(
  Surv(Survival_time, Status_cs1) ~ lp_cca,
  data = cca_data
)

# C-index
c_cca <- sc_cca$concordance

# Standard error: prefer std.err if present, otherwise derive from var
if (!is.null(sc_cca$std.err)) {
  se_cca <- sc_cca$std.err
} else {
  se_cca <- sqrt(as.numeric(sc_cca$var))
}

ci_cca <- c_cca + c(-1.96, 1.96) * se_cca

cindex_cca_summary <- data.frame(
  Method = "Complete case",
  C_index = c_cca,
  SE      = se_cca,
  LCL     = ci_cca[1],
  UCL     = ci_cca[2]
)

cindex_cca_summary

# Save to disk
write.csv(cindex_cca_summary,
          "C_index_complete_case_cs1_extended.csv",
          row.names = FALSE)
saveRDS(cindex_cca_summary,
        "C_index_complete_case_cs1_extended.RDS")

preimp_cc$complete_case <- complete.cases(preimp_cc[vars_model])

library(dplyr)

baseline_compare <- preimp_cc %>%
  group_by(complete_case) %>%
  summarise(
    N = n(),
    Events = sum(Status_cs1 == 1, na.rm = TRUE),
    Event_rate = mean(Status_cs1 == 1, na.rm = TRUE),
    Age_mean = mean(Age, na.rm = TRUE),
    BMI_mean = mean(BMI, na.rm = TRUE),
    Diabetes_pct = mean(Diabetes == 1, na.rm = TRUE) * 100,
    Hypertension_pct = mean(Hypertension == 1, na.rm = TRUE) * 100,
    LVEF_mean = mean(LVEF, na.rm = TRUE),
    eGFR_mean = mean(eGFR, na.rm = TRUE),
    Stroke_pct = mean(Stroke_TIA == 1, na.rm = TRUE) * 100,
    MI_history_pct = mean(MI_history == 1, na.rm = TRUE) * 100
  )
baseline_compare

table(preimp_cc$Status_cs1, preimp_cc$complete_case)
prop.table(table(preimp_cc$Status_cs1, preimp_cc$complete_case), 2)

library(ggplot2)
library(patchwork)



continuous_vars <- c(
  "BMI","Age","LVEF","eGFR","Haemoglobin",
  "Cholesterol","HDL","LDL","Triglycerides"
)

plot_list <- lapply(continuous_vars, function(v) {
  ggplot(preimp_cc, aes(x = .data[[v]], fill = complete_case)) +
    geom_density(alpha = 0.35) +
    labs(
      title = v,
      x = v, y = "Density"
    ) +
    scale_fill_manual(
      values = c("#E55372","#1F9E89"),
      labels = c("Excluded","Complete")
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
})

# Combine into a grid
panel_cont <- wrap_plots(plot_list, ncol = 3)
panel_cont

cat_vars <- c("Sex","Diabetes","Hypertension","Smoking","COPD_cat","Stroke_TIA")

plot_list_cat <- lapply(cat_vars, function(v) {
  
  df_v <- preimp_cc %>%
    group_by(complete_case, !!sym(v)) %>%
    summarise(n = n(), .groups="drop") %>%
    group_by(complete_case) %>%
    mutate(p = n / sum(n))
  
  ggplot(df_v, aes(x = !!sym(v), y = p, fill = complete_case)) +
    geom_bar(stat = "identity", position = "dodge") +
    labs(title = v, x = v, y = "Proportion") +
    scale_fill_manual(
      values = c("#E55372","#1F9E89"),
      labels = c("Excluded","Complete")
    ) +
    theme_minimal(base_size = 11) +
    theme(legend.position = "none")
})

panel_cat <- wrap_plots(plot_list_cat, ncol = 3)
panel_cat

legend_plot <- ggplot(preimp_cc, aes(x=Age, fill=complete_case)) +
  geom_density(alpha=0.35) +
  scale_fill_manual(values=c("#E55372","#1F9E89"), labels=c("Excluded","Complete")) +
  theme_minimal() +
  theme(legend.position="bottom", legend.title=element_blank()) +
  guides(fill=guide_legend(nrow=1)) +
  labs(x="", y="", title=NULL)

panel_cont_final <- panel_cont / legend_plot + plot_layout(heights=c(1,0.2))
panel_cont_final

ggsave(
  filename = "panel_density_continuous_CCA_vs_excluded.pdf",
  plot = panel_cont,
  width = 13, height = 8,
  dpi = 300
)

ggsave(
  filename = "panel_bar_categorical_CCA_vs_excluded.pdf",
  plot = panel_cat,
  width = 13, height = 7,
  dpi = 300

)
