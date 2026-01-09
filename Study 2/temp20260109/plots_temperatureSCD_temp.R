library(ggplot2)
library(dplyr)
library(splines)
library(broom)

#############################################
# 1. Identify the realistic temperature range
#############################################

temp_p1  <- quantile(final_env$temp_mean, 0.01, na.rm = TRUE)
temp_p99 <- quantile(final_env$temp_mean, 0.99, na.rm = TRUE)

temp_seq <- seq(temp_p1, temp_p99, length.out = 300)

pred_data <- data.frame(
  temp_mean = temp_seq,
  DB = final_env$DB[1]   # dummy value for strata
)

#############################################
# 2. Predict log-HR + 95% CI
#############################################

pred <- cbind(
  pred_data,
  predict(model_temp_extremes, newdata = pred_data,
          type = "lp", se.fit = TRUE)
)

pred <- pred %>%
  mutate(
    HR      = exp(fit),
    HR_low  = exp(fit - 1.96 * se.fit),
    HR_high = exp(fit + 1.96 * se.fit)
  )

#############################################
# 3. Find the nadir (minimum risk temperature)
#############################################

nadir_temp <- pred$temp_mean[which.min(pred$HR)]

#############################################
# 4. Publication-ready plot
#############################################

ggplot(pred, aes(x = temp_mean, y = HR)) +
  geom_ribbon(aes(ymin = HR_low, ymax = HR_high),
              fill = "#A6CEE3", alpha = 0.3) +
  geom_line(size = 1.2, colour = "#1F78B4") +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red", size = 0.8) +
  geom_vline(xintercept = nadir_temp, linetype = "dotted",
             colour = "darkgreen", size = 0.8) +
  annotate("text", x = nadir_temp, y = max(pred$HR),
           label = paste0("Nadir ≈ ", round(nadir_temp, 1), "°C"),
           colour = "darkgreen", angle = 90, vjust = -0.5, size = 4.5) +
  labs(
    x = "Temperature (°C)",
    y = "Relative Hazard of SCD",
    title = "Temperature–SCD Risk Curve",
    subtitle = "Spline-based Cox model with registry stratification",
    caption = "Shaded area = 95% CI"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold", size = 18),
    plot.subtitle = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14)
  )

library(dplyr)
library(splines)
library(ggplot2)

## Split datasets
single_env <- final_env %>% filter(DB %in% single_regs)
multi_env  <- final_env %>% filter(DB %in% multi_regs)

## Fit models separately
model_temp_single <- coxph(
  Surv(Survival_time, event) ~ ns(temp_mean, df = 4) + strata(DB),
  data = single_env
)

model_temp_multi <- coxph(
  Surv(Survival_time, event) ~ ns(temp_mean, df = 4) + strata(DB),
  data = multi_env
)

## Helper to make prediction grid from a model + data subset
make_pred_df <- function(model, dat, group_label) {
  t1  <- quantile(dat$temp_mean, 0.01, na.rm = TRUE)
  t99 <- quantile(dat$temp_mean, 0.99, na.rm = TRUE)
  temp_seq <- seq(t1, t99, length.out = 300)
  
  newdat <- data.frame(
    temp_mean = temp_seq,
    DB = dat$DB[1]   # dummy for strata term
  )
  
  p <- predict(model, newdata = newdat, type = "lp", se.fit = TRUE)
  
  newdat %>%
    mutate(
      HR      = exp(p$fit),
      HR_low  = exp(p$fit - 1.96 * p$se.fit),
      HR_high = exp(p$fit + 1.96 * p$se.fit),
      group   = group_label
    )
}

pred_single <- make_pred_df(model_temp_single, single_env, "Single-centre")
pred_multi  <- make_pred_df(model_temp_multi,  multi_env,  "Multi-centre")

pred_both <- bind_rows(pred_single, pred_multi)

## Plot
ggplot(pred_both, aes(x = temp_mean, y = HR, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = HR_low, ymax = HR_high), alpha = 0.20, colour = NA) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  labs(
    x = "Temperature (°C)",
    y = "Relative Hazard of SCD",
    colour = "Registry type",
    fill   = "Registry type",
    title = "Temperature–SCD Risk by Registry Type",
    subtitle = "Single-centre vs multi-centre registries"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold")
  )     


## Keep only Winter & Summer
ws_env <- final_env %>% filter(season %in% c("Winter","Summer"))

winter_env <- ws_env %>% filter(season == "Winter")
summer_env <- ws_env %>% filter(season == "Summer")

model_temp_winter <- coxph(
  Surv(Survival_time, event) ~ temp_mean  + strata(DB),
  data = winter_env
)

model_temp_summer <- coxph(
  Surv(Survival_time, event) ~ temp_mean  + strata(DB),
  data = summer_env
)

## Re-use helper function
pred_winter <- make_pred_df(model_temp_winter, winter_env, "Winter")
pred_summer <- make_pred_df(model_temp_summer, summer_env, "Summer")

pred_ws <- bind_rows(pred_winter, pred_summer)

ggplot(pred_ws, aes(x = temp_mean, y = HR, colour = group, fill = group)) +
  geom_ribbon(aes(ymin = HR_low, ymax = HR_high), alpha = 0.20, colour = NA) +
  geom_line(size = 1.1) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "red") +
  labs(
    x = "Temperature (°C)",
    y = "Relative Hazard of SCD",
    colour = "Season",
    fill   = "Season",
    title = "Temperature–SCD Risk by Season",
    subtitle = "Overlay of Winter and Summer curves"
  ) +
  theme_minimal(base_size = 15) +
  theme(
    plot.title = element_text(face = "bold")
  )


library(purrr)
library(broom)

## Optional: scale temp per 5°C for interpretability
final_env <- final_env %>%
  mutate(temp_mean_5 = temp_mean / 5)

## Fit a separate Cox model per registry: HR per 5°C increase
registry_fits <- final_env %>%
  group_by(DB) %>%
  group_split() %>%
  map_df(function(dat) {
    # require at least, say, 20 events per DB
    if (sum(dat$event) < 20) {
      return(NULL)
    }
    
    fit <- coxph(
      Surv(Survival_time, event) ~ temp_mean_5,
      data = dat
    )
    
    broom::tidy(fit) %>%
      filter(term == "temp_mean_5") %>%
      mutate(
        DB   = unique(dat$DB),
        HR   = exp(estimate),
        HR_lo = exp(estimate - 1.96 * std.error),
        HR_hi = exp(estimate + 1.96 * std.error)
      )
  })

## Order registries by HR
registry_fits <- registry_fits %>%
  arrange(HR) %>%
  mutate(DB = factor(DB, levels = DB))

## Forest plot
ggplot(registry_fits,
       aes(x = DB, y = HR, ymin = HR_lo, ymax = HR_hi)) +
  geom_hline(yintercept = 1, linetype = "dashed", colour = "grey50") +
  geom_pointrange(size = 0.6) +
  coord_flip() +
  labs(
    x = "Registry",
    y = "HR per 5°C increase in temperature",
    title = "Temperature–SCD Association by Registry"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    axis.title.y = element_text(margin = margin(r = 10))
  )

