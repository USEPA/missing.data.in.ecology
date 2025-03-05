library(tidyverse)
library(spsurvey)
library(here)
library(sf)
library(maps)
library(mapdata)

# read in NWCA data
nwca <- read_csv(here("inst", "data", "nwca_2016.csv"))

# # add a national estimate variable
nwca <- nwca |>
  mutate(NATIONAL = "National")

# change name of NA
nwca <- nwca %>%
  mutate(NTL_COND = if_else(NTL_COND == "Not Assessed", "No Water", NTL_COND))

nwca_sf <- nwca %>%
  st_as_sf(crs = 4269, coords = c("LON_DD83", "LAT_DD83")) %>%
  # st_transform(crs = 5070) %>%
  mutate(NTL_COND = factor(NTL_COND, levels = c("Good", "Fair", "Poor", "No Water")))

# display tables
table(nwca_sf$NTL_COND)

usa <- map_data("usa")
usa_sf <- usa %>%
  st_as_sf(crs = 4326, coords = c("long", "lat")) %>%
  st_transform(crs = 4269)
usa_sf_coords <- st_coordinates(usa_sf)
usa_sf$long <- usa_sf_coords[, 1]
usa_sf$lat <- usa_sf_coords[, 2]

nwca_sf_plot <- ggplot(nwca_sf, aes(color = NTL_COND)) +
  # group aesthetic needed so points are connected
  geom_path(data = usa_sf, aes(x = long, y = lat, group = group), color = "black", inherit.aes = FALSE, alpha = 0.4) +
  geom_sf(size = 0.8) +
  scale_color_viridis_d(name = "Category", option = "A", begin = 0.3, end = 0.9) +
  facet_wrap(~ NTL_COND) +
  labs(x = "Nitrogen Condition", y = "") +
  theme_bw(base_size = 16) +
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
  filename = here("inst", "figures", "water-present-0.jpeg"),
  plot = nwca_sf_plot,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)

# estimate nitrogen condition
ntl_cond_miss <- cat_analysis(
  dframe = nwca,
  vars = "NTL_COND",
  subpops = c("NATIONAL"),
  weight = "WGT_TP",
  xcoord = "XCOORD",
  ycoord = "YCOORD"
)

ntl_cond_miss <- ntl_cond_miss |>
  filter(Category != "Total") |>
  mutate(Category = factor(Category, levels = c("Good", "Fair", "Poor", "Missing")))

# table of estimates
ntl_cond_miss |>
  select(Subpopulation, Category, nResp, Estimate.P, MarginofError.P)

ntl_plot_miss <- ggplot(ntl_cond_miss, aes(x = Category, y = Estimate.P, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = LCB95Pct.P, ymax = UCB95Pct.P), color = "black", linewidth = 0.75, width = 0.5) +
  labs(y = "Percent", x = "Nitrogen Condition Class") +
  scale_fill_viridis_d(option = "A", begin = 0.3, end = 0.9) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank()
  )

# default 6.03 x 4.97 in
ggsave(
  filename = here("inst", "figures", "water-present-1.jpeg"),
  plot = ntl_plot_miss,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)

nwca <- nwca |>
  mutate(SW_SAMPLEABLE = replace_na(SW_SAMPLEABLE, "M")) |>
  mutate(CHEM_NOT_COLLECTED = replace_na(CHEM_NOT_COLLECTED, "M")) |>
  mutate(MICX_NOT_COLLECTED = replace_na(MICX_NOT_COLLECTED, "M")) |>
  mutate(WCHL_NOT_COLLECTED = replace_na(WCHL_NOT_COLLECTED, "M")) |>
  mutate(CMW_NOT_COLLECTED = str_c(CHEM_NOT_COLLECTED, MICX_NOT_COLLECTED,
                                    WCHL_NOT_COLLECTED)) |>
  mutate(WATER_PRESENT = case_when(
    CMW_NOT_COLLECTED == "YYY" ~ "No",
    CMW_NOT_COLLECTED == "MMY" ~ "Yes",
    CMW_NOT_COLLECTED == "MMM" ~ "Yes")
  )

# add missing no water
nwca <- nwca |>
  mutate(NTL_MISSING_REASON = case_when(
    NTL_COND != "Missing" ~ "Not Missing",
    NTL_COND == "Missing" & WATER_PRESENT == "No" ~ "No Water",
    NTL_COND == "Missing" & WATER_PRESENT == "Yes" ~ "Other"
  ))

# reason for missing
nwca_missing <- nwca |>
  filter(NTL_MISSING_REASON != "Not Missing")
# only missing because of no water (rest of commented code is unnecessary)
table(nwca_missing$NTL_MISSING_REASON)

h20_cond_miss <- cat_analysis(
  dframe = nwca,
  vars = "WATER_PRESENT",
  subpops = c("NATIONAL"),
  weight = "WGT_TP",
  xcoord = "XCOORD",
  ycoord = "YCOORD"
)

h20_cond_miss <- h20_cond_miss |>
  filter(Category != "Total") |>
  mutate(Category = factor(Category, levels = c("Yes", "No")))

# library(viridis)
# magma(12)

h20_plot <- ggplot(h20_cond_miss, aes(x = Category, y = Estimate.P, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = LCB95Pct.P, ymax = UCB95Pct.P), color = "black", linewidth = 0.75, width = 0.5) +
  labs(y = "Percent", x = "Water Present") +
  scale_fill_manual(values = c("Yes" = "#A5CFE3", "No" = "#FEC98DFF")) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank()
  )

ggsave(
  filename = here("inst", "figures", "water-present-2.jpeg"),
  plot = h20_plot,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)

ntl_cond_nomiss <- cat_analysis(
  dframe = nwca |> filter(WATER_PRESENT == "Yes"),
  vars = "NTL_COND",
  subpops = c("NATIONAL"),
  weight = "WGT_TP",
  xcoord = "XCOORD",
  ycoord = "YCOORD"
)

ntl_cond_nomiss <- ntl_cond_nomiss |>
  filter(Category != "Total") |>
  mutate(Category = factor(Category, levels = c("Good", "Fair", "Poor", "Missing")))

ntl_plot_nomiss <- ggplot(ntl_cond_nomiss, aes(x = Category, y = Estimate.P, fill = Category)) +
  geom_bar(stat = "identity", width = 0.5) +
  geom_errorbar(aes(ymin = LCB95Pct.P, ymax = UCB95Pct.P), color = "black", linewidth = 0.75, width = 0.5) +
  labs(y = "Percent (Water Present)", x = "Nitrogen Condition Class") +
  scale_fill_viridis_d(option = "A", begin = 0.3, end = 0.7) +
  theme_bw(base_size = 20) +
  theme(
    legend.position = "top",
    axis.text.x = element_blank()
  )

ggsave(
  filename = here("inst", "figures", "water-present-3.jpeg"),
  plot = ntl_plot_nomiss,
  dpi = 300,
  units = "in",
  width = 6.03,
  height = 4.97
)

write_csv(ntl_cond_miss, here("inst", "output", "water-present", "ntl_cond_with_missing.csv"))
write_csv(h20_cond_miss, here("inst", "output", "water-present", "water_present.csv"))
write_csv(ntl_cond_nomiss, here("inst", "output", "water-present", "ntl_cond_water_present.csv"))
