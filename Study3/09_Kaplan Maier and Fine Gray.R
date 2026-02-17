#
# install.packages("data.table")
# install.packages("cmprsk")
# install.packages("survival")
# install.packages("readxl")
# install.packages ("mice")

library(data.table)
library(cmprsk)
library(survival)
library(readxl)
library (mice)

dt_final <- readRDS("study3_analysis_final.rds")

km_fit <- survfit(
  Surv(t_followup_days_final, event_inapp_shock) ~ device_group,
  data = dt_final
)

survdiff(
  Surv(t_followup_days_final, event_inapp_shock) ~ device_group,
  data = dt_final
)


plot(
  km_fit,
  col = c("red", "blue"),
  lwd = 2,
  xlab = "Days since ICD implantation",
  ylab = "Event-free survival (inappropriate shock)",
  xlim = c(0, 4000),
  mark.time = TRUE
)

legend("bottomleft", legend = levels(factor(dt_final$device_group)),
       col = c("red", "blue"), lwd = 2)



##### quality check 

table(
  dt_final$device_group,
  is.na(dt_final$inapp_shock_flag)
)
summary(dt_final[event_inapp_shock == 1, t_followup_days_final])

names(dt_final)


dt_final <- readRDS("study3_analysis_final.rds")
dt_final <- as.data.table(dt_final)

# Fine–Gray event variable
# 0 = censored
# 1 = inappropriate shock 
# 2 = death before inapproiate shock - competing
dt_final[, fg_event := 0L]
dt_final[event_inapp_shock == 1, fg_event := 1L]
dt_final[event_inapp_shock == 0 & death_flag == "Yes", fg_event := 2L]

# QC  
cat("N =", nrow(dt_final), "\n")
print(table(dt_final$fg_event)) 

stopifnot(all(!is.na(dt_final$t_followup_days_final)))
stopifnot(all(dt_final$t_followup_days_final > 0))
stopifnot(all(dt_final$fg_event %in% c(0L, 1L, 2L)))

# Fine–Gray: sHR for Single vs Dual (Dual reference)
X <- model.matrix(~ device_group, data = dt_final)[, -1, drop = FALSE]

fg_fit <- crr(
  ftime   = dt_final$t_followup_days_final,
  fstatus = dt_final$fg_event,
  cov1    = X
)

fg_sum <- summary(fg_fit)
print(fg_sum)

# Extract subdistribution HR (sHR), 95% CI, and p-value
beta <- fg_fit$coef[1]
se   <- sqrt(fg_fit$var[1, 1])

sHR  <- exp(beta)
lci  <- exp(beta - 1.96 * se)
uci  <- exp(beta + 1.96 * se)
pval <- 2 * pnorm(-abs(beta / se))

cat(
  "\nFine–Gray (Single vs Dual): ",
  "sHR = ", round(sHR, 2),
  ", 95% CI [", round(lci, 2), ", ", round(uci, 2), "]",
  ", p = ", signif(pval, 3),
  "\n", sep = ""
)

cif <- with(
  dt_final,
  cuminc(
    t_followup_days_final,
    fg_event,
    group = device_group
  )
)

cif_inapp <- cif[grep(" 1$", names(cif))]

plot(
  cif_inapp,
  lwd = 2,
  lty = c(1, 2),
  curvlab = c("Dual-chamber ICD", "Single-chamber ICD"),
  xlab = "Days since ICD implantation",
  ylab = "Cumulative incidence of inappropriate shock",
  xlim = c(0, 5500),
  ylim = c(0, 1)
)

