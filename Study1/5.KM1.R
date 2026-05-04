###############################################################################

# KM 1: Time to first inappropriate ICD shock

###############################################################################
pkgs <- c(
  "tidyverse","dplyr" ,"survival", "data.table", "mice", "naniar", "grid","gridExtra",
  "openxlsx", "readxl", "gt", "ggplot2"
)

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))
install.packages("ggsurvplot")
library(ggsurvplot)
install.packages("survminer")
library(survminer)
install.packages("patchwork")
library(patchwork)
library(scales)

# ── Paths ──────────────────────────────────────────────────────────────────
INFILE <- "T:/Study_1/master_clean_dataset1.rds"
OUTDIR <- "T:/Study_1/Supplementary_data"

dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ── Load full cohort ───────────────────────────────────────────────────────
d <- as.data.table(readRDS(INFILE))
cat(sprintf("Full cohort loaded: N = %d\n", nrow(d)))

# ── Missing Status_FIS treated as unexposed (consistent with main analysis) ─
#d[is.na(Status_FIS), Status_FIS := 0L]

# ── Inclusion: valid follow-up time only ───────────────────────────────────
d_km <- d[!is.na(Time_FIS_days) & Time_FIS_days > 0]
cat(sprintf("Patients included in KM: N = %d\n", nrow(d_km)))
cat(sprintf("Exposed (Status_FIS = 1): N = %d (%.1f%%)\n",
            sum(d_km$Status_FIS),
            mean(d_km$Status_FIS) * 100))

# ── Convert days to years ──────────────────────────────────────────────────
d_km[, Time_FIS_years := Time_FIS_days / 365.25]

# ── KM fit ─────────────────────────────────────────────────────────────────
fit <- survfit(Surv(Time_FIS_years, Status_FIS) ~ 1, data = d_km)

# ── CHECK 1: Median estimability ───────────────────────────────────────────
# Median = time at which survival drops to 0.50 (i.e. cumulative incidence reaches 50%)
# If the curve never reaches 50% cumulative incidence, quantile() returns NA
med <- quantile(fit, probs = 0.5)$quantile

cat("\n--- Median estimability check ---\n")
cat(sprintf("Maximum cumulative incidence reached: %.1f%%\n",
            (1 - min(fit$surv)) * 100))
cat(sprintf("Required for median estimability:     50.0%%\n"))

if (is.na(med)) {
  cat("RESULT: Median NOT estimable — curve never reached 50% cumulative incidence\n")
  cat("        This is expected given the low event rate (~2.2%)\n")
  median_text <- "not estimable (cumulative incidence did not reach 50%)"
} else {
  cat(sprintf("RESULT: Median estimable = %.2f years\n", med))
  median_text <- sprintf("%.2f years", med)
}

# ── CHECK 2: Cumulative incidence at 1, 3, 5 years ────────────────────────
cat("\n--- Cumulative incidence at 1, 3, 5 years ---\n")
sm <- summary(fit, times = c(1, 3, 5), extend = TRUE)

ci_tbl <- data.table(
  Year              = c(1, 3, 5),
  N_at_risk         = sm$n.risk,
  Cum_incidence_pct = round((1 - sm$surv)  * 100, 1),
  CI_lower_pct      = round((1 - sm$upper) * 100, 1),
  CI_upper_pct      = round((1 - sm$lower) * 100, 1)
)
print(ci_tbl)

# Build caption text from results
ci_caption <- paste(
  sprintf("%s%% (95%% CI %s\u2013%s%%) at %d year%s",
          ci_tbl$Cum_incidence_pct,
          ci_tbl$CI_lower_pct,
          ci_tbl$CI_upper_pct,
          ci_tbl$Year,
          ifelse(ci_tbl$Year == 1, "", "s")),
  collapse = ", "
)

caption_text <- paste0(
  "Full analysis cohort (N = ", nrow(d_km), "). ",
  "Median time to first inappropriate ICD shock: ", median_text, ". ",
  "Cumulative incidence: ", ci_caption, ". ",
  "These curves are descriptive and do not account for the time-dependent ",
  "nature of exposure."
)

cat(sprintf("\nCaption text:\n%s\n", caption_text))
# ── Build KM data frame for ggplot ─────────────────────────────────────────
km_df <- data.table(
  time  = c(0, fit$time),
  surv  = c(1, fit$surv),
  lower = c(1, fit$lower),
  upper = c(1, fit$upper)
)
km_df[, `:=`(
  cuminc = 1 - surv,
  ci_lo  = 1 - upper,
  ci_hi  = 1 - lower
)]

# ── Risk table data ────────────────────────────────────────────────────────
risk_times <- 0:5
sm_risk    <- summary(fit, times = risk_times, extend = TRUE)

risk_df <- data.table(
  time   = sm_risk$time,
  n_risk = sm_risk$n.risk
)

# ── Main KM plot ───────────────────────────────────────────────────────────
p_main <- ggplot(km_df, aes(x = time, y = cuminc)) +
  geom_step(aes(y = ci_lo), colour = "black",
            linewidth = 0.5, linetype = "dashed") +
  geom_step(aes(y = ci_hi), colour = "black",
            linewidth = 0.5, linetype = "dashed") +
  geom_step(colour = "black", linewidth = 0.9) +
  scale_x_continuous(
    breaks = 0:5,
    limits = c(0, 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, NA),
    expand = c(0.02, 0)
  ) +
  coord_cartesian(xlim = c(0, 5), clip = "off") +  # clip = "off" key fix
  labs(
    x = NULL,
    y = "Cumulative incidence of first\ninappropriate ICD shock"
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text        = element_text(size = 11, colour = "black"),
    axis.title.y     = element_text(size = 11, margin = margin(r = 8)),
    axis.text.x      = element_blank(),
    axis.ticks.x     = element_blank(),
    axis.line        = element_line(colour = "black")
  )

# ── Risk table panel ───────────────────────────────────────────────────────
p_risk <- ggplot(risk_df, aes(x = time, y = 0.5, label = n_risk)) +
  geom_text(size = 3.8, fontface = "plain", colour = "black") +
  annotate("text", x = -0.35, y = 0.5,
           label = "Number\nat risk",
           size = 3.5, hjust = 1, fontface = "plain",
           colour = "grey30") +
  scale_x_continuous(
    breaks = 0:5,
    limits = c(0, 5),
    labels = c("0", "1", "2", "3", "4", "5"),
    expand = c(0, 0)
  ) +
  scale_y_continuous(limits = c(0, 1)) +
  coord_cartesian(xlim = c(0, 5), clip = "off") +  # clip = "off" key fix
  labs(x = "Years", y = NULL) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x      = element_text(size = 11, colour = "black"),
    axis.title.x     = element_text(size = 11, margin = margin(t = 6)),
    axis.text.y      = element_blank(),
    axis.ticks.y     = element_blank(),
    axis.line.y      = element_blank(),
    axis.line.x      = element_line(colour = "black"),
    panel.grid       = element_blank()
  )

# ── Combine with forced alignment via & operator ───────────────────────────
p_final <- p_main / p_risk +
  plot_layout(heights = c(5, 1)) &
  theme(plot.margin = margin(5, 15, 5, 100))

# ── Save ───────────────────────────────────────────────────────────────────
out_png <- file.path(OUTDIR, "KM_FigA_time_to_FIS_fullcohort.png")
out_pdf <- file.path(OUTDIR, "KM_FigA_time_to_FIS_fullcohort.pdf")

suppressWarnings({
  ggsave(out_png, plot = p_final, width = 10, height = 5.5, dpi = 300)
  ggsave(out_pdf, plot = p_final, width = 10, height = 5.5)
})

cat(sprintf("\n✓ Saved:\n  %s\n  %s\n", out_png, out_pdf))
# ── Figure legend text for manuscript (NOT embedded in figure) ─────────────
cat("\n--- Figure legend (paste into manuscript) ---\n")
cat(sprintf(
  "Figure A. Kaplan-Meier curve depicting the cumulative incidence of first
inappropriate ICD shock over time. Full analysis cohort (N = %d).
Median time to first inappropriate ICD shock was not estimable, as the
cumulative incidence did not reach 50%% during follow-up, consistent with
the low overall event rate. Cumulative incidence at 1, 3, and 5 years
was %.1f%% (95%% CI %.1f%%\u2013%.1f%%), %.1f%% (95%% CI %.1f%%\u2013%.1f%%),
and %.1f%% (95%% CI %.1f%%\u2013%.1f%%), respectively. Dashed lines represent
95%% confidence intervals. These curves are descriptive and do not account
for the time-dependent nature of the exposure.\n",
  nrow(d_km),
  ci_tbl[Year==1, Cum_incidence_pct],
  ci_tbl[Year==1, CI_lower_pct],
  ci_tbl[Year==1, CI_upper_pct],
  ci_tbl[Year==3, Cum_incidence_pct],
  ci_tbl[Year==3, CI_lower_pct],
  ci_tbl[Year==3, CI_upper_pct],
  ci_tbl[Year==5, Cum_incidence_pct],
  ci_tbl[Year==5, CI_lower_pct],
  ci_tbl[Year==5, CI_upper_pct]
))