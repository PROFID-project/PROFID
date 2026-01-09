df_env_linked <- df_astn %>%
  filter(Status == 1 & !is.na(event_date)) %>%
  left_join(climate_merged, by = c("event_date" = "date", "DB" = "DB"))

library(slider)

climate_lagged <- climate_merged %>%
  group_by(center_id) %>%
  arrange(center_id, date) %>%
  mutate(temp_mean_lag0_3 = slide_dbl(temp_mean, mean, .before = 3, .complete = TRUE),
         pressure_lag0_3 = slide_dbl(pressure, mean, .before = 3, .complete = TRUE),
         daylight_lag0_3 = slide_dbl(daylight, mean, .before = 3, .complete = TRUE))

df_env_linked <- df_astn %>%
  filter(Status == 1 & !is.na(event_date)) %>%
  left_join(climate_lagged, by = c("event_date" = "date", "DB" = "DB"))

ggplot(df_env_linked, aes(x = temp_mean)) +
  geom_histogram(bins = 30, fill = "skyblue") +
  labs(title = "Distribution of Temperature on SCD Event Dates")

df_env_linked$scd_event <- ifelse(df_astn$Status == 1, 1, 0)

log_model <- glm(scd_event ~ temp_mean + pressure + daylight + Age + Sex + LVEF_std,
                 data = df_env_linked, family = "binomial")
summary(log_model)

library(survival)
library(dlnm)

cox_env <- coxph(Surv(Survival_time, Status == 1) ~ temp_mean + pressure + daylight + Age + Sex + LVEF_std,
                 data = df_env_linked)
summary(cox_env)
