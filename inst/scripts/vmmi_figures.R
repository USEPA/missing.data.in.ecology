# vmmi figures
library(tidyverse)
library(here)
library(sf)
library(cowplot)

dat_full <- read_csv(here("inst", "data", "nwca_2016.csv"), guess_max = Inf)

dat <- dat_full %>%
  mutate(SAMPLEABLE_WATER = if_else(SAMPLEABLE_H2O == "Yes", 1, 0)) %>%
  mutate(IS_WOODY = if_else(WETCLS_GRP %in% c("PRLW", "EW"), 1, 0)) %>%
  mutate(NTL_RESULT = if_else(is.na(NTL_RESULT), 0, log(NTL_RESULT))) %>%
  select(VMMI_2016, IS_WOODY, SAMPLEABLE_WATER, NTL_RESULT,
         STRESS_VEGRMV, PALT_SOHARD, PALT_SOMODF, LON_DD83, LAT_DD83, XCOORD, YCOORD) %>%
  mutate(PALT_SOMODF = as.numeric(PALT_SOMODF), PALT_SOHARD = as.numeric(PALT_SOHARD),
         STRESS_VEGRMV = factor(STRESS_VEGRMV, levels = c("Low", "Moderate", "High"))) %>%
  filter(!is.na(PALT_SOMODF), !is.na(VMMI_2016))

# map of vmmi
dat_sf <- dat %>%
  st_as_sf(coords = c("LON_DD83", "LAT_DD83"), crs = 4269)

usa <- map_data("usa")
usa_sf <- usa %>%
  st_as_sf(crs = 4326, coords = c("long", "lat")) %>%
  st_transform(crs = 4269)
usa_sf_coords <- st_coordinates(usa_sf)
usa_sf$long <- usa_sf_coords[, 1]
usa_sf$lat <- usa_sf_coords[, 2]



vmmi_sf_plot <- ggplot(dat_sf, aes(color = VMMI_2016)) +
  # group aesthetic needed so points are connected
  geom_path(data = usa_sf, aes(x = long, y = lat, group = group), color = "black", inherit.aes = FALSE, alpha = 0.4) +
  geom_sf(size = 1) +
  scale_color_viridis_c(name = "VMMI", option = "H", begin = 0, end = 1) +
  theme_bw(base_size = 12) +
  labs(x = "", y = "") +
  theme(
    legend.position = "right",
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  )

# default 6.03 x 4.97 in
ggsave(
  filename = here("inst", "figures", "vmmi-0.jpeg"),
  plot = vmmi_sf_plot,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)

dat_sf <- dat_sf %>%
  mutate(SAMPLEABLE = if_else(SAMPLEABLE_WATER == 1, "Yes", "No"))
sampleable_sf_plot <- ggplot(dat_sf, aes(color = SAMPLEABLE)) +
  # group aesthetic needed so points are connected
  geom_path(data = usa_sf, aes(x = long, y = lat, group = group), color = "black", inherit.aes = FALSE, alpha = 0.4) +
  geom_sf(size = 1) +
  scale_color_viridis_d(name = "Surface Water Present", option = "B", begin = 0.2, end = 0.8) +
  theme_bw(base_size = 16) +
  labs(x = "", y = "") +
  theme(
    legend.position = "top",
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

# default 6.03 x 4.97 in
ggsave(
  filename = here("inst", "figures", "sampleable-0.jpeg"),
  plot = sampleable_sf_plot,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)


dat_sf <- dat_sf %>%
  mutate(NTL_PLOT = case_when(
    SAMPLEABLE_WATER == 0 ~ "No Water",
    NTL_RESULT < 0 ~ "[-3, 0)",
    between(NTL_RESULT, 0, 2) ~ "[0, 2)",
    between(NTL_RESULT, 2, 5) ~ "[2, 5)"
  )) %>%
  mutate(NTL_PLOT = factor(NTL_PLOT, levels = c("No Water", "[-3, 0)", "[0, 2)", "[2, 5)")))
sampleable_sf_plot2 <- ggplot(dat_sf, aes(color = NTL_PLOT)) +
  # group aesthetic needed so points are connected
  geom_path(data = usa_sf, aes(x = long, y = lat, group = group), color = "black", inherit.aes = FALSE, alpha = 0.4) +
  geom_sf(size = 1) +
  scale_color_viridis_d(name = "Log Nitrogen", option = "H", begin = 0.9, end = 0) +
  theme_bw(base_size = 12) +
  labs(x = "", y = "") +
  theme(
    legend.position = "right",
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    panel.grid = element_blank(),
    axis.ticks = element_blank()
  ) +
  guides(color = guide_legend(override.aes = list(size = 3)))

ggsave(
  filename = here("inst", "figures", "sampleable-1.jpeg"),
  plot = sampleable_sf_plot2,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)

#
nitrogen_plot <- ggplot(dat_sf %>% filter(SAMPLEABLE_WATER == 1), aes(x = NTL_RESULT, y = VMMI_2016)) +
  geom_point() +
  theme_bw(base_size = 16) +
  labs(x = "Log Nitrogen (Surface Water Present)", y = "VMMI") +
  theme(
  )

ggsave(
  filename = here("inst", "figures", "nitrogen-0.jpeg"),
  plot = nitrogen_plot,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)

nitrogen_plot2 <- ggplot(dat_sf, aes(x = SAMPLEABLE, y = VMMI_2016)) +
  geom_boxplot() +
  theme_bw(base_size = 16) +
  labs(x = "Surface Water Presence", y = "VMMI") +
  theme(
  )

ggsave(
  filename = here("inst", "figures", "nitrogen-1.jpeg"),
  plot = nitrogen_plot2,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)


vmmi_nitr <- plot_grid(vmmi_sf_plot, sampleable_sf_plot2, ncol = 1, align = "v")

ggsave(
  filename = here("inst", "figures", "vmmi-nitr-0.jpeg"),
  plot = vmmi_nitr,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)

loocv_miss <- read_csv(here("inst", "output", "vmmi", "loocv_miss.csv"))
loocv_cc <- read_csv(here("inst", "output", "vmmi", "loocv_cc.csv"))

predr2 <- left_join(loocv_miss, loocv_cc, by = "response") %>%
  pivot_longer(`sp-missing`:`none-cc`)

predr2_plot <- ggplot(predr2, aes(x = value, y = response)) +
  geom_point() +
  geom_abline(color = "red", linewidth = 1) +
  lims(x = c(0, 100), y = c(0, 100)) +
  facet_wrap(~ name)

