###############################################
# TASK 2.4 – ENVIRONMENTAL & LIFESTYLE FACTORS
# WEATHER PATTERN INTEGRATION – FULL PIPELINE
###############################################

library(dplyr)
library(lubridate)
library(splines)
library(survival)

########################################################
# 0. INPUTS REQUIRED:
#
# (A) df_cleaned: patient-level dataset, including:
#     ID, DB, center_id, event_date, time, event
#
# (B) climate_all: stacked daily climate dataset for centres you have:
#     center_id, date, temp_mean, temp_min, temp_max,
#     pressure, daylight, sunshine
#
# (C) centres: table linking centre_id to DB, Country, City, etc.
#
########################################################
df_cleaned <- readRDS("T:/PROFID/data/processed/df_cleaned.rds")




# Read climate data
climate_all <- readRDS("T:/PROFID/data/raw/era5_city_raw_until2022/climate_all_centers.rds")

# Read the metadata for research centres
research_centres <- read.csv("T:/PROFID/data/raw/research_centres_UPDATED.txt") %>%
  rename(center_id = ID)

########################################################
# 1. DEFINE REGISTRY TYPES: SINGLE vs MULTI-CENTRE
########################################################

single_regs <- c("ASTN","ATMS","HELS","ISAR","NANC",
                 "OLMC","SLSN")

multi_regs  <- c("DOIT","FREN","ISRL",
                 "PRDT","PRSE",
                 "SWHR")

exclude_regs <- c("MDRT", "MDII", "SHFT", "DRVT")


# Attach CERT center IDs

cert <- read.csv(
  "T:/PROFID/data/raw/centreIDSEUCERT.csv",
  stringsAsFactors = FALSE
)



df_env_base <- df_cleaned %>%
  filter(!DB %in% exclude_regs)


df_env_base <- df_env_base %>%
  mutate(
    CERT_pat_id = trimws(as.character(CERT_pat_id))
  )

cert <- cert %>%
  mutate(
    pat_id = trimws(as.character(pat_id))
  )

df_env_base <- df_env_base %>%
  select(-ctr_name) %>%   # remove old join
  left_join(cert, by = c("CERT_pat_id" = "pat_id"))


# check
df_env_base %>%
  filter(DB == "CERT") %>%
  summarise(
    total = n(),
    missing = sum(is.na(ctr_name))
  )


########################################################
# 2. PREPARE CLIMATE DATA
########################################################
df_env <- df_env_base %>%
  filter(!DB %in% exclude_regs) %>%   # KEEP ONLY needed registries
  filter(!is.na(event_date),
         !is.na(Time_zero_Ym),
         !is.na(Status))

df_env <- df_env %>%
  select(DB, event_date, Survival_time, Time_zero_Ym, Status,
         ID, Age, Sex, BMI, LVEF_std, ctr_name)

df_env %>%
  filter(DB == "CERT") %>%
  summarise(
    n_CERT_included = n(),
    n_unique_centres = n_distinct(ctr_name)
  )

library(stringi)

df_env <- df_env %>%
  mutate(
    ctr_name_clean = stri_trans_general(ctr_name, "Latin-ASCII")
  )

df_env <- df_env %>%
  left_join(
    centers,
    by = c("ctr_name_clean" = "Center")
  )


min_date <- min(df_env$event_date, na.rm = TRUE)  # "1995-03-01"
max_date <- max(df_env$event_date, na.rm = TRUE)  # "2021-06-17"


climate_env <- climate_all %>%
  filter(date >= min_date, date <= max_date)



# 2.1 Add DB and Country to climate dataset
climate_env <- climate_env %>%
  left_join(research_centres %>% select(center_id, DB, Country),
            by = "center_id")

# 2.2 Convert dates and add seasons
climate_env <- climate_env %>%
  mutate(
    date   = as.Date(date),
    season = factor(quarters(date),
                    labels = c("Winter","Spring","Summer","Autumn"))
  )



########################################################
# 3. CREATE COUNTRY-LEVEL SEASONAL CLIMATE (MULTI-CENTRE)
########################################################

climate_country_season <- climate_env %>%
  filter(DB %in% multi_regs) %>%
  group_by(DB, Country, season) %>%
  summarise(
    temp_mean = mean(temp_mean, na.rm = TRUE),
    temp_min  = mean(temp_min,  na.rm = TRUE),
    temp_max  = mean(temp_max,  na.rm = TRUE),
    pressure  = mean(pressure,  na.rm = TRUE),
    daylight  = mean(daylight,  na.rm = TRUE),
    sunshine  = mean(sunshine,  na.rm = TRUE),
    .groups   = "drop"
  )



########################################################
# 4. CREATE CITY-LEVEL DAILY CLIMATE (SINGLE-CENTRE)
########################################################

climate_city_daily <- climate_env %>%
  filter(DB %in% single_regs)



########################################################
# 5. PREPARE PATIENT DATA FOR MERGING
########################################################

df_env <- df_env %>%
  mutate(
    event_date = as.Date(event_date),
    season = factor(quarters(event_date),
                    labels=c("Winter","Spring","Summer","Autumn"))
  )

# Add country to patient data (needed for multi-centre merge)
df_env <- df_env %>%
  left_join(research_centres %>% select(DB, Country), 
            by = "DB")




########################################################
# 6. LINK CLIMATE TO SINGLE-CENTRE PATIENTS (DAILY)
########################################################

single_env <- df_env %>%
  filter(DB %in% single_regs) %>%
  left_join(climate_city_daily,
            by = c("event_date" = "date",
                   "DB" = "DB"))



########################################################
# 7. LINK CLIMATE TO MULTI-CENTRE PATIENTS (COUNTRY-SEASON)
########################################################

multi_env <- df_env %>%
  filter(DB %in% multi_regs) %>%
  left_join(climate_country_season,
            by = c("DB", "Country", "season"))



########################################################
# 8. COMBINE FINAL ENVIRONMENTAL DATASET
########################################################

final_env <- bind_rows(single_env, multi_env)



########################################################
# 9. MAIN ANALYSES
########################################################
final_env <- df_env
########################################################
# 9A. Temperature Extremes – ns() spline
########################################################
final_env <- final_env %>%
  mutate(event = ifelse(Status == 1, 1, 0))


model_temp_extremes <- coxph(
  Surv(Survival_time, event) ~ ns(temp_mean, df=4) + strata(DB),
  data = final_env
)

summary(model_temp_extremes)

library(splines)

# Create a sequence of temperatures to predict on
temp_range <- seq(min(final_env$temp_mean, na.rm=TRUE),
                  max(final_env$temp_mean, na.rm=TRUE),
                  length = 300)

pred_temp <- data.frame(
  temp_mean = temp_range,
  DB = final_env$DB[1]  # needed for strata; value doesn't matter
)

# Predict linear predictor
lp <- predict(model_temp_extremes, newdata = pred_temp, type = "lp")
pred_temp$HR <- exp(lp)

plot(pred_temp$temp_mean, pred_temp$HR, type="l",
     xlab="Temperature (°C)",
     ylab="Relative Hazard of SCD",
     main="Temperature–SCD Risk Curve (Spline)")

summary(final_env$temp_mean)
quantile(final_env$temp_mean, probs = c(0.01, 0.05, 0.95, 0.99), na.rm = TRUE)

predict(model_temp_extremes, 
        newdata = data.frame(temp_mean = c(-10, 0, 20), DB = final_env$DB[1]),
        type = "risk")
abline(h=1, col="red", lty=2)


########################################################
# 9B. Barometric Pressure
########################################################

model_pressure <- coxph(
  Surv(time=Survival_time, event) ~ scale(pressure) + strata(DB),
  data = final_env
)

summary(model_pressure)

library(ggplot2)
library(dplyr)

# Pressure range
p1  <- quantile(final_env$pressure, 0.01, na.rm=TRUE)
p99 <- quantile(final_env$pressure, 0.99, na.rm=TRUE)
pressure_seq <- seq(p1, p99, length.out = 300)

pred_p <- data.frame(
  pressure = pressure_seq,
  DB = final_env$DB[1]
)

pp <- predict(model_pressure, newdata = pred_p, type = "lp", se.fit = TRUE)

pred_p <- pred_p %>%
  mutate(
    HR      = exp(pp$fit),
    HR_low  = exp(pp$fit - 1.96 * pp$se.fit),
    HR_high = exp(pp$fit + 1.96 * pp$se.fit)
  )

ggplot(pred_p, aes(x=pressure, y=HR)) +
  geom_ribbon(aes(ymin=HR_low, ymax=HR_high), fill="#C7E9C0", alpha=0.4) +
  geom_line(size=1.2, color="#238B45") +
  geom_hline(yintercept=1, linetype="dashed", color="red") +
  labs(
    title = "Barometric Pressure–SCD Risk Curve",
    subtitle = "Cox model with registry stratification",
    x = "Atmospheric Pressure (hPa)",
    y = "Relative Hazard of SCD"
  ) +
  theme_minimal(base_size=15) +
  theme(plot.title = element_text(face="bold", size=18))



########################################################
# 9C. Daylight Hours (Photoperiod)
########################################################

model_daylight <- coxph(
  Surv(time=Survival_time, event) ~ daylight + strata(DB),
  data = final_env
)

summary(model_daylight)

# Daylight range
d1  <- quantile(final_env$daylight, 0.01, na.rm=TRUE)
d99 <- quantile(final_env$daylight, 0.99, na.rm=TRUE)
day_seq <- seq(d1, d99, length.out = 300)

pred_d <- data.frame(
  daylight = day_seq,
  DB = final_env$DB[1]
)

pd <- predict(model_daylight, newdata = pred_d, type="lp", se.fit=TRUE)

pred_d <- pred_d %>%
  mutate(
    HR = exp(pd$fit),
    HR_low = exp(pd$fit - 1.96 * pd$se.fit),
    HR_high = exp(pd$fit + 1.96 * pd$se.fit)
  )

ggplot(pred_d, aes(x = daylight, y = HR)) +
  geom_ribbon(aes(ymin=HR_low, ymax=HR_high), fill="#9ECAE1", alpha=0.4) +
  geom_line(size=1.2, color="#08519C") +
  geom_hline(yintercept=1, linetype="dashed", color="red") +
  labs(
    title = "Daylight Duration–SCD Risk Curve",
    subtitle = "Cox model with registry stratification",
    x = "Daylight Hours",
    y = "Relative Hazard of SCD"
  ) +
  theme_minimal(base_size=15) +
  theme(plot.title = element_text(face="bold", size=18))


########################################################
# 9D. Sunshine Duration
########################################################

model_sunshine <- coxph(
  Surv(time=Survival_time, event) ~ sunshine + strata(DB),
  data = final_env
)

summary(model_sunshine)

# Sunshine range
s1  <- quantile(final_env$sunshine, 0.01, na.rm=TRUE)
s99 <- quantile(final_env$sunshine, 0.99, na.rm=TRUE)
sun_seq <- seq(s1, s99, length.out = 300)

pred_s <- data.frame(
  sunshine = sun_seq,
  DB = final_env$DB[1]
)

ps <- predict(model_sunshine, newdata = pred_s, type="lp", se.fit=TRUE)

pred_s <- pred_s %>%
  mutate(
    HR = exp(ps$fit),
    HR_low = exp(ps$fit - 1.96 * ps$se.fit),
    HR_high = exp(ps$fit + 1.96 * ps$se.fit)
  )

ggplot(pred_s, aes(x = sunshine, y = HR)) +
  geom_ribbon(aes(ymin = HR_low, ymax = HR_high), fill="#FDD0A2", alpha=0.4) +
  geom_line(size=1.2, color="#D94801") +
  geom_hline(yintercept = 1, linetype = "dashed", color="red") +
  labs(
    title = "Sunshine Duration–SCD Risk Curve",
    subtitle = "Cox model with registry stratification",
    x = "Sunshine (hours/day)",
    y = "Relative Hazard of SCD"
  ) +
  theme_minimal(base_size=15) +
  theme(plot.title = element_text(face="bold", size=18))


########################################################
# 9E. Air Quality (variable does not exist)
########################################################

# if ("air_quality" %in% colnames(final_env)) {
#   model_air <- coxph(
#     Surv(time, event) ~ air_quality + strata(DB),
#     data = final_env
#   )
#   print(summary(model_air))
# }



########################################################
# 10. SENSITIVITY ANALYSES
########################################################

# 10A. Single-centre only
model_single_only <- coxph(
  Surv(time=Survival_time,event) ~ ns(temp_mean, df=4),
  data = final_env %>% filter(DB %in% single_regs)
)

summary(model_single_only)

# 10B. Multi-centre only
model_multi_only <- coxph(
  Surv(time=Survival_time,event) ~ temp_mean + pressure + daylight,
  data = final_env %>% filter(DB %in% multi_regs)
)

summary(model_multi_only)


# 10C. Winter-only (cold effects)
model_winter <- coxph(
  Surv(time=Survival_time,event) ~ temp_mean + pressure + daylight + strata(DB),
  data = final_env %>% filter(season == "Winter")
)

summary(model_winter)



########################################################
# 11. SAVE OUTPUT DATASET & MODELS
########################################################

saveRDS(final_env, "final_environmental_dataset.rds")
saveRDS(model_temp_extremes, "model_temperature.rds")
saveRDS(model_pressure, "model_pressure.rds")
saveRDS(model_daylight, "model_daylight.rds")
saveRDS(model_sunshine, "model_sunshine.rds")

########################################################
# END OF PIPELINE
########################################################
# 
# 4-Panel Combined Figure (2×2 Layout)

library(ggplot2)
library(patchwork)

# Temperature range (only where data exist)
t1  <- quantile(final_env$temp_mean, 0.01, na.rm=TRUE)
t99 <- quantile(final_env$temp_mean, 0.99, na.rm=TRUE)
temp_seq <- seq(t1, t99, length.out = 300)

pred_temp <- data.frame(
  temp_mean = temp_seq,
  DB = final_env$DB[1]
)

pt <- predict(model_temp_extremes, newdata = pred_temp, type = "lp", se.fit = TRUE)

pred_temp <- pred_temp %>%
  mutate(
    HR      = exp(pt$fit),
    HR_low  = exp(pt$fit - 1.96 * pt$se.fit),
    HR_high = exp(pt$fit + 1.96 * pt$se.fit)
  )


# Panel A – Temperature
p_temp <- ggplot(pred_temp, aes(x=temp_mean, y=HR)) +
  geom_ribbon(aes(ymin=HR_low, ymax=HR_high), fill="#BDD7E7", alpha=0.4) +
  geom_line(size=1.2, colour="#08519C") +
  geom_hline(yintercept=1, linetype="dashed", colour="red") +
  labs(
    title="A) Temperature–SCD Risk (Spline)",
    x="Temperature (°C)",
    y="HR"
  ) +
  theme_minimal(base_size=14)

# Panel B – Pressure
p_press <- ggplot(pred_p, aes(x=pressure, y=HR)) +
  geom_ribbon(aes(ymin=HR_low, ymax=HR_high), fill="#C7E9C0", alpha=0.4) +
  geom_line(size=1.2, colour="#238B45") +
  geom_hline(yintercept=1, linetype="dashed", colour="red") +
  labs(
    title="B) Barometric Pressure–SCD Risk",
    x="Pressure (hPa)",
    y="HR"
  ) +
  theme_minimal(base_size=14)

# Panel C – Daylight
p_day <- ggplot(pred_d, aes(x=daylight, y=HR)) +
  geom_ribbon(aes(ymin=HR_low, ymax=HR_high), fill="#9ECAE1", alpha=0.4) +
  geom_line(size=1.2, colour="#08519C") +
  geom_hline(yintercept=1, linetype="dashed", colour="red") +
  labs(
    title="C) Daylight–SCD Risk",
    x="Daylight Hours",
    y="HR"
  ) +
  theme_minimal(base_size=14)

# Panel D – Sunshine
p_sun <- ggplot(pred_s, aes(x=sunshine, y=HR)) +
  geom_ribbon(aes(ymin=HR_low, ymax=HR_high), fill="#FDD0A2", alpha=0.4) +
  geom_line(size=1.2, colour="#D94801") +
  geom_hline(yintercept=1, linetype="dashed", colour="red") +
  labs(
    title="D) Sunshine–SCD Risk",
    x="Sunshine (hrs/day)",
    y="HR"
  ) +
  theme_minimal(base_size=14)

# Combined figure (2x2)
combined_plot <- (p_temp | p_press) / (p_day | p_sun)

combined_plot


ggplot() +
  geom_ribbon(data=pred_single, aes(x=temp_mean, ymin=HR_low, ymax=HR_high, fill="Single-centre"), alpha=0.25) +
  geom_line(data=pred_single, aes(x=temp_mean, y=HR, colour="Single-centre"), size=1.2) +
  geom_ribbon(data=pred_multi, aes(x=temp_mean, ymin=HR_low, ymax=HR_high, fill="Multi-centre"), alpha=0.25) +
  geom_line(data=pred_multi, aes(x=temp_mean, y=HR, colour="Multi-centre"), size=1.2) +
  geom_hline(yintercept=1, linetype="dashed", colour="red") +
  labs(
    title="Temperature–SCD Risk by Registry Type",
    subtitle="Daily-resolution climate data vs seasonal multi-centre averages",
    x="Temperature (°C)",
    y="HR",
    colour="Registry type",
    fill="Registry type"
  ) +
  scale_colour_manual(values=c("Single-centre"="#08519C","Multi-centre"="#CB181D")) +
  scale_fill_manual(values=c("Single-centre"="#9ECAE1","Multi-centre"="#FCBBA1")) +
  theme_minimal(base_size=15)


ggplot() +
  geom_ribbon(data=pred_winter, aes(x=temp_mean, ymin=HR_low, ymax=HR_high, fill="Winter"), alpha=0.25) +
  geom_line(data=pred_winter, aes(x=temp_mean, y=HR, colour="Winter"), size=1.2) +
  geom_ribbon(data=pred_summer, aes(x=temp_mean, ymin=HR_low, ymax=HR_high, fill="Summer"), alpha=0.25) +
  geom_line(data=pred_summer, aes(x=temp_mean, y=HR, colour="Summer"), size=1.2) +
  geom_hline(yintercept=1, linetype="dashed", colour="red") +
  labs(
    title="Temperature–SCD Risk by Season",
    subtitle="Colder periods demonstrate stronger risk elevation",
    x="Temperature (°C)",
    y="HR",
    colour="Season",
    fill="Season"
  ) +
  scale_colour_manual(values=c("Winter"="#08519C","Summer"="#D94801")) +
  scale_fill_manual(values=c("Winter"="#9ECAE1","Summer"="#FDD0A2")) +
  theme_minimal(base_size=15)




valid_single <- single_regs

season_env <- final_env %>% 
  filter(DB %in% valid_single)

season_env <- final_env %>% 
  filter(DB %in% single_regs)

winter_env <- season_env %>% filter(season == "Winter")
summer_env <- season_env %>% filter(season == "Summer")

model_temp_winter <- coxph(
  Surv(Survival_time, event) ~ ns(temp_mean, df=2) + strata(DB),
  data = winter_env
)

model_temp_summer <- coxph(
  Surv(Survival_time, event) ~ ns(temp_mean, df=2) + strata(DB),
  data = summer_env
)

