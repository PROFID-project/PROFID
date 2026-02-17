
################ SAP sensitivity #########

library(data.table)
library(survival)
library(cmprsk)

dt_final <- as.data.table(readRDS("study3_analysis_final.rds"))

#  6-month filter (SAP)
min_fu_days <- 180

# QC:  excluded COUNTS
dt_final[, fu_lt_6mo := t_followup_days_final < min_fu_days]
cat("\n--- SAP 6-month filter QC ---\n")
print(dt_final[, .(
  n_total = .N,
  n_lt_6mo = sum(fu_lt_6mo, na.rm = TRUE),
  pct_lt_6mo = round(100 * mean(fu_lt_6mo, na.rm = TRUE), 2),
  events_total = sum(event_inapp_shock, na.rm = TRUE),
  events_lt_6mo = sum(event_inapp_shock == 1 & fu_lt_6mo, na.rm = TRUE)
)])

# Sensitivity cohort: >= 180 days
dt_sens_6mo <- dt_final[t_followup_days_final >= min_fu_days]

cat("\n--- Sensitivity cohort size (>=180d) ---\n")
print(dt_sens_6mo[, .(
  n = .N,
  events = sum(event_inapp_shock, na.rm = TRUE)
)])

#  dataset x device breakdown
cat("\n--- Sens cohort: events by dataset x device ---\n")
print(dt_sens_6mo[, .(
  n = .N,
  events = sum(event_inapp_shock, na.rm = TRUE)
), by = .(dataset, device_group)][order(dataset, device_group)])

### Incidence rates per 100 person-years ###########
incidence_rates_6mo <- rbindlist(
  lapply(c("event_inapp_shock", "app_shock_flag", "death_flag"), function(v) {
    
    # events definition:
    # - event_inapp_shock is already 0/1
    # - app_shock_flag / death_flag are "Yes"/"No"/NA in your pipeline
    if (v == "event_inapp_shock") {
      dt_sens_6mo[, .(
        events = sum(event_inapp_shock == 1, na.rm = TRUE),
        person_years = round(sum(t_followup_days_final, na.rm = TRUE) / 365.25, 1)
      ), by = device_group][
        , endpoint := "inapp_shock_flag"
      ][
        , {
          py <- person_years
          ev <- events
          rate <- ev / py * 100
          ci <- poisson.test(ev, py)$conf.int * 100
          .(events = ev,
            person_years = py,
            rate_per_100py = round(rate, 2),
            lci_95 = round(ci[1], 2),
            uci_95 = round(ci[2], 2),
            endpoint = endpoint)
        }, by = device_group
      ]
    } else {
      
      # Recode flags to binary for rate calculation
      tmp <- copy(dt_sens_6mo)
      tmp[, flag_bin := fifelse(tolower(as.character(get(v))) == "yes", 1L, 0L)]
      
      tmp[, .(
        events = sum(flag_bin, na.rm = TRUE),
        person_years = round(sum(t_followup_days_final, na.rm = TRUE) / 365.25, 1)
      ), by = device_group][
        , endpoint := v
      ][
        , {
          py <- person_years
          ev <- events
          rate <- ev / py * 100
          ci <- poisson.test(ev, py)$conf.int * 100
          .(events = ev,
            person_years = py,
            rate_per_100py = round(rate, 2),
            lci_95 = round(ci[1], 2),
            uci_95 = round(ci[2], 2),
            endpoint = endpoint)
        }, by = device_group
      ]
    }
  })
)

# Tidy endpoint names 
incidence_rates_6mo[endpoint == "death_flag", endpoint := "death_flag"]
incidence_rates_6mo[endpoint == "app_shock_flag", endpoint := "app_shock_flag"]

print(incidence_rates_6mo)

fwrite(incidence_rates_6mo, "incidence_rates_per_100py_by_device_sens_min180d.csv")

############# Cox model (stratified by dataset) — sensitivity###################
cox_sens_6mo <- coxph(
  Surv(t_followup_days_final, event_inapp_shock) ~ device_group + strata(dataset),
  data = dt_sens_6mo
)
cat("\n--- Cox sensitivity (>=180d) ---\n")
print(summary(cox_sens_6mo))

# PH test
ph_sens_6mo <- cox.zph(cox_sens_6mo)
cat("\n--- PH test (>=180d) ---\n")
print(ph_sens_6mo)

####################################### Kaplan–Meier + log-rank — sensitivity ####################################
km_sens_6mo <- survfit(
  Surv(t_followup_days_final, event_inapp_shock) ~ device_group,
  data = dt_sens_6mo
)

logrank_sens_6mo <- survdiff(
  Surv(t_followup_days_final, event_inapp_shock) ~ device_group,
  data = dt_sens_6mo
)

cat("\n--- Log-rank (>=180d) ---\n")
print(logrank_sens_6mo)

png("figure_km_inapp_shock_by_device_sens_min180d.png", width = 1400, height = 900, res = 160)
plot(
  km_sens_6mo,
  col = c("red", "blue"),
  lwd = 2,
  xlab = "Days since ICD implantation",
  ylab = "Event-free survival (inappropriate shock)",
  mark.time = TRUE
)
legend("bottomleft",
       legend = levels(factor(dt_sens_6mo$device_group)),
       col = c("red", "blue"), lwd = 2, bty = "n")
dev.off()

################################Fine–Gray competing risks — sensitivity (>=180d ##########################################
dt_sens_6mo[, fg_event := 0L]
dt_sens_6mo[event_inapp_shock == 1, fg_event := 1L]
dt_sens_6mo[event_inapp_shock == 0 & tolower(as.character(death_flag)) == "yes", fg_event := 2L]

cat("\n--- Fine–Gray event table (>=180d) ---\n")
print(table(dt_sens_6mo$fg_event))

X <- model.matrix(~ device_group, data = dt_sens_6mo)[, -1, drop = FALSE]

fg_sens_6mo <- crr(
  ftime   = dt_sens_6mo$t_followup_days_final,
  fstatus = dt_sens_6mo$fg_event,
  cov1    = X
)

cat("\n--- Fine–Gray sensitivity (>=180d) ---\n")
print(summary(fg_sens_6mo))

# Extract sHR, CI, p
beta <- fg_sens_6mo$coef[1]
se   <- sqrt(fg_sens_6mo$var[1, 1])
sHR  <- exp(beta)
lci  <- exp(beta - 1.96 * se)
uci  <- exp(beta + 1.96 * se)
pval <- 2 * pnorm(-abs(beta / se))

cat(
  "\nFine–Gray sensitivity (>=180d): ",
  "sHR = ", round(sHR, 2),
  ", 95% CI [", round(lci, 2), ", ", round(uci, 2), "]",
  ", p = ", signif(pval, 3),
  "\n", sep = ""
)

# CIF plot 
cif_sens_6mo <- with(
  dt_sens_6mo,
  cuminc(t_followup_days_final, fg_event, group = device_group)
)

png("figure_cif_inapp_shock_by_device_sens_min180d.png", width = 1400, height = 900, res = 160)
plot(
  cif_sens_6mo,
  lwd = 2,
  xlab = "Days since ICD implantation",
  ylab = "Cumulative incidence (inappropriate shock)"
)
dev.off()



############ outocmes

dt_sens_6mo[, .(
  n = .N,
  events = sum(event_inapp_shock)
)]

dt_sens_6mo[, .(
  n = .N,
  events = sum(event_inapp_shock)
), by = device_group]

summary(cox_sens_6mo)
cox.zph(cox_sens_6mo)
summary(fg_sens_6mo)


saveRDS(dt_sens_6mo, "study3_analysis_sensitivity_min180d.rds")
