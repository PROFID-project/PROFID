###############################################################################

# KM 2: Comparison overall survival in patients who experienced inappropriate ICD shock 
#       and who did not.

###############################################################################
pkgs <- c(
  "tidyverse","dplyr" ,"survival", "data.table", "mice", "naniar", "grid","gridExtra",
  "openxlsx", "readxl", "gt", "ggplot2"
)

invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))


# ── Paths ──────────────────────────────────────────────────────────────────
INFILE <- "T:/Study_1/master_clean_dataset1.rds"
OUTDIR <- "T:/Study_1/Supplementary_data"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)

# ── Load full cohort ───────────────────────────────────────────────────────
d <- as.data.table(readRDS(INFILE))
cat(sprintf("Full cohort loaded: N = %d\n", nrow(d)))

# ── Missing Status_FIS treated as unexposed ────────────────────────────────
d[is.na(Status_FIS), Status_FIS := 0L]

# ── Inclusion: valid Time_death_days and death status ─────────────────────
d_km <- d[
  !is.na(Status_death)    &
    !is.na(Time_death_days) &
    !is.na(Status_FIS)      &
    Time_death_days > 0
]
cat(sprintf("Patients included in KM: N = %d\n", nrow(d_km)))
cat(sprintf("Shocked (Status_FIS = 1):     N = %d (%.1f%%)\n",
            d_km[Status_FIS == 1L, .N],
            d_km[Status_FIS == 1L, .N] / nrow(d_km) * 100))
cat(sprintf("Not shocked (Status_FIS = 0): N = %d (%.1f%%)\n",
            d_km[Status_FIS == 0L, .N],
            d_km[Status_FIS == 0L, .N] / nrow(d_km) * 100))

# ── Convert Time_death_days to years ──────────────────────────────────────
d_km[, Time_death_years := Time_death_days / 365.25]

# ── Group label ───────────────────────────────────────────────────────────
d_km[, shock_group := fifelse(
  Status_FIS == 1L,
  "Inappropriate shock",
  "No inappropriate shock"
)]
d_km[, shock_group := factor(shock_group,
                             levels = c("No inappropriate shock",
                                        "Inappropriate shock"))]







# ── KM fit ─────────────────────────────────────────────────────────────────
fit <- survfit(Surv(Time_death_years, Status_death) ~ shock_group, data = d_km)

# ── Deaths by group — print for manuscript ────────────────────────────────
cat("\n--- Deaths by shock group ---\n")
death_tab <- d_km[, .(
  N      = .N,
  Deaths = sum(Status_death),
  Pct    = round(mean(Status_death) * 100, 1)
), by = shock_group]
print(death_tab)

# ── Build KM data frame for ggplot ─────────────────────────────────────────
s     <- summary(fit)
km_df <- data.table(
  time  = s$time,
  surv  = s$surv,
  lower = s$lower,
  upper = s$upper,
  group = s$strata
)

# Clean group labels
km_df[, group := gsub("shock_group=", "", group)]
km_df[, group := factor(group,
                        levels = c("No inappropriate shock",
                                   "Inappropriate shock"))]

# Add time = 0 rows for clean start
t0 <- data.table(
  time  = c(0, 0),
  surv  = c(1, 1),
  lower = c(1, 1),
  upper = c(1, 1),
  group = factor(c("No inappropriate shock", "Inappropriate shock"),
                 levels = c("No inappropriate shock", "Inappropriate shock"))
)
km_df <- rbind(t0, km_df)

# ── Colour scheme ──────────────────────────────────────────────────────────
col_map <- c("No inappropriate shock" = "#2c6e96",   # blue
             "Inappropriate shock"    = "#c0392b")   # red

# ── Plot ───────────────────────────────────────────────────────────────────
p_final <- ggplot(km_df, aes(x = time, y = surv,
                             group = group, colour = group)) +
  geom_step(linewidth = 0.9) +
  scale_colour_manual(values = col_map) +
  scale_x_continuous(
    breaks = 0:5,
    limits = c(0, 5),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    labels = percent_format(accuracy = 1),
    limits = c(0, 1),
    breaks = seq(0, 1, 0.25),    # 0, 25, 50, 75, 100%
    expand = c(0.02, 0)
  ) +
  coord_cartesian(xlim = c(0, 5), clip = "off") +
  labs(
    x      = "Years",
    y      = "Overall survival",
    colour = NULL
  ) +
  theme_classic(base_size = 12) +
  theme(
    axis.text         = element_text(size = 11, colour = "black"),
    axis.title        = element_text(size = 11),
    axis.title.x      = element_text(margin = margin(t = 6)),
    axis.title.y      = element_text(margin = margin(r = 8)),
    axis.line         = element_line(colour = "black"),
    legend.position   = c(0.75, 0.15),
    legend.text       = element_text(size = 10),
    legend.key.width  = unit(1.5, "cm"),
    legend.background = element_blank(),
    plot.margin       = margin(10, 15, 10, 15)
  )

# ── Save ───────────────────────────────────────────────────────────────────
out_png <- file.path(OUTDIR, "KM_FigB_survival_by_shock_fullcohort.png")
out_pdf <- file.path(OUTDIR, "KM_FigB_survival_by_shock_fullcohort.pdf")

suppressWarnings({
  ggsave(out_png, plot = p_final, width = 7.5, height = 5.5, dpi = 300)
  ggsave(out_pdf, plot = p_final, width = 7.5, height = 5.5)
})

cat(sprintf("\n✓ Saved:\n  %s\n  %s\n", out_png, out_pdf))

# ── Figure legend text for manuscript (NOT embedded in figure) ─────────────
cat("\n--- Figure legend (paste into manuscript) ---\n")
cat(sprintf(
  "Figure S1. Kaplan-Meier curves comparing overall survival between patients
who experienced inappropriate ICD shock (N = %d) and those who did not
(N = %d). Full analysis cohort (N = %d). Follow-up time measured from
ICD implantation. Blue line = no inappropriate shock; red line =
inappropriate shock. These curves do not account for the time-dependent
nature of the exposure and should be interpreted as descriptive
summaries only.\n",
  d_km[Status_FIS == 1L, .N],
  d_km[Status_FIS == 0L, .N],
  nrow(d_km)
))