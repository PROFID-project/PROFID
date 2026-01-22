###############################################
# Graduated symbol map and Getis–Ord Gi* (ICD and Non-ICD cohorts)
# - Graduated symbol map
# - Getis–Ord Gi* analysis
###############################################

# Packages
req <- c("data.table","dplyr","tidyverse","ggplot2","sf","spdep","rnaturalearth","rnaturalearthdata","ggrepel")
new_pkgs <- setdiff(req, installed.packages()[,"Package"])
if(length(new_pkgs)) install.packages(new_pkgs, quiet = TRUE)
invisible(lapply(req, library, character.only = TRUE))

# Input files
df <- fread("T:/Data Transfer to Charite/imp/ICD_imp_geo_df.csv")

# Graduated symbol map
geo_sf <- df %>%
  st_as_sf(coords = c("longitude", "latitude"), crs = 4326, remove = FALSE)

world_sf <- ne_countries(
  scale = "medium",
  returnclass = "sf"
)

p1 <- ggplot() +
  geom_sf(
    data = world_sf,
    fill = "grey95",
    colour = "grey80",
    linewidth = 0.2
  ) +
  geom_sf(
    data = geo_sf,
    aes(
      size = crude_rate_plot,
      colour = crude_rate_plot,
      shape = rate_missing
    ),
    alpha = 0.85
  ) +
  geom_text_repel(
    data = geo_sf,
    aes(
      geometry = geometry,
      label = ctr_name
    ),
    stat = "sf_coordinates",
    size = 2,
    min.segment.length = 0,
    box.padding = 0.3,
    seed = 123,
    max.overlaps = 20
  ) +
  scale_colour_viridis_c(
    option = "plasma",
    name = "Incidence rate\n(per 100,000 PY)"
  ) +
  scale_size_continuous(
    range = c(3, 10),
    name = "Relative incidence\n(per 100,000 PY)"
  ) +
  scale_shape_manual(
    values = c("FALSE" = 16, "TRUE" = 1),
    name = "Missing rate"
  ) +
  guides(
    size = guide_legend(
        override.aes = list(alpha = 0.6)
      )
  ) +
  coord_sf(
    xlim = c(-15, 35),
    ylim = c(35, 70)
  ) +
  labs(
    title = "Centre-level crude SCD incidence rates",
    subtitle = "ICD cohort after imputation",
    caption = "Each point represents a recruiting centre\nCentres without abbreviations belongs to CERT database"
    # caption = "Each point represents a recruiting centre"    #for Non-ICD cohorts
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 16),
    plot.subtitle = element_text(size = 12),
    panel.grid.major = element_line(colour = "transparent"),
    axis.title = element_blank(),
    legend.position = "right",
    legend.box = "vertical",
    legend.title = element_text(size = 11),
    legend.text = element_text(size = 10)
  )
ggsave(p1, filename = "T:/Data Transfer to Charite/imp/ICD_imp_gsm.png")

# Getis–Ord Gi* analysis
geo_sf_fil <- geo_sf %>%
  filter(!is.na(crude_rate))

coords <- st_coordinates(geo_sf_fil)
lw <- nb2listw(knn2nb(knearneigh(coords, k = 4)), style = "W", zero.policy = TRUE)

geo_sf_fil <- geo_sf_fil %>%
  mutate(
    Gi_star = as.numeric(localG(geo_sf_fil$crude_rate, lw, zero.policy = TRUE)),
    Gi_category = case_when(
      Gi_star >= 1.96  ~ "Hotspot",
      Gi_star <= -1.96 ~ "Cold spot",
      TRUE             ~ "Not significant"
    )
  )

p2 <- ggplot() +
  geom_sf(
    data = world_sf, 
    fill = "grey95", 
    color = "grey80", 
    linewidth = 0.2
  ) +
  geom_sf(
    data = geo_sf_fil,
    aes(color = Gi_category),
    size = 2
  ) +
  geom_text_repel(
    data = geo_sf_fil,
    aes(
      x = longitude,
      y = latitude,
      label = ctr_name
    ),
    size = 2,
    max.overlaps = 20
  ) +
  scale_color_manual(
    name = "Getis–Ord Gi*",
    values = c(
      "Hotspot" = "#D73027",
      "Cold spot" = "#4575B4",
      "Not significant" = "grey50"
    )
  ) +
  labs(
    title = "Centre-level clustering of crude SCD incidence",
    subtitle = "ICD cohort after imputation",
    caption = "Red = high-incidence clusters; Blue = low-incidence clusters\nCentres without abbreviations belongs to CERT database"
  ) +
  coord_sf(
    xlim = c(-15, 35), 
    ylim = c(35, 70)
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "right",
    plot.title = element_text(face = "bold"),
    panel.grid = element_blank()
  )
ggsave(p2, filename = "T:/Data Transfer to Charite/imp/ICD_imp_hotspot.png")