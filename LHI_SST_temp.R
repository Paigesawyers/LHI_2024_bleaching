##===============================
## 2024 MARINE HEATWAVE ANALYSIS – Lord Howe Island
## SST / In Situ & Satellite Data
## Consistent NOAA-style DHW calculations
##===============================

# -----------------------------
# 0. Clean environment
# -----------------------------
rm(list = ls())

setwd("/users/paigesawyers/Desktop/LHI/LHI_sawyers_et_al_2026/paige_all_files/")

# -----------------------------
# 1. Load required packages
# -----------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  tidyverse, lubridate, janitor, zoo, data.table,
  ggplot2, patchwork, viridis, scales, cowplot, paletteer, dplyr
)

# -----------------------------
# 2. Define paths & global parameters
# -----------------------------
insitu_start <- ymd("2023-12-01")
insitu_end   <- ymd("2024-05-31")
sat_years    <- c(1998, 2010, 2019, 2024)
dhw_window_days <- 84

# -----------------------------
# 3. Load site MMMs
# -----------------------------
mmms <- read_csv("LHI_mmms.csv") %>%
  clean_names() %>%
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")))

# -----------------------------
# NOAA-style DHW helper
# HotSpot = temperature - MMM, only where temperature >= MMM + 1
# DHW = sum HotSpots over previous 84 days / 7
# -----------------------------
calc_dhw <- function(hotspot, width = 84, partial = FALSE) {
  zoo::rollapply(
    hotspot,
    width = width,
    FUN = function(x) sum(x, na.rm = TRUE) / 7,
    fill = NA,
    align = "right",
    partial = partial
  )
}

# -----------------------------
# Theme
# -----------------------------
theme_noaa <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(size = 12),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 12),
      axis.title = element_text(face = "bold", size = 12),
      axis.text = element_text(size = 12),
      plot.title = element_text(face = "bold", size = 12),
      plot.subtitle = element_text(face = "bold", size = 12),
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12)
    )
}

scale_dhw <- function(df, dhw_col = "dhw", sst_max = 30) {
  df %>%
    group_by(site) %>%
    mutate(
      scale_factor = sst_max / 20,
      dhw_scaled   = .data[[dhw_col]] * scale_factor,
      dhw4_scaled  = 4 * scale_factor,
      dhw8_scaled  = 8 * scale_factor
    ) %>%
    ungroup()
}

# -----------------------------
# 4. Load raw in situ logger data
# -----------------------------
insitu_raw <- read_csv("LHIMP_DATA_in_situ_loggers.csv") %>%
  clean_names() %>%
  filter(date_time != "") %>%
  mutate(
    date_time = dmy_hms(paste0(date_time, ":00")),
    date_only = as_date(date_time)
  ) %>%
  mutate(across(where(is.numeric), ~na_if(., 0)))

# -----------------------------
# 5. Convert in situ data to long format
# -----------------------------
insitu_long <- insitu_raw %>%
  pivot_longer(
    observatory_rock:coral_gardens,
    names_to = "site",
    values_to = "temperature"
  ) %>%
  filter(!is.na(temperature)) %>%
  mutate(site = str_to_lower(str_replace_all(site, " ", "_"))) %>%
  left_join(mmms, by = "site") %>%
  filter(!is.na(mmm)) %>%
  filter(date_only >= insitu_start, date_only <= insitu_end) %>%
  filter(site != "comets_hole")

# -----------------------------
# 6. Daily in situ mean, max, min temperatures
# -----------------------------
insitu_daily <- insitu_long %>%
  group_by(site, date_only, mmm) %>%
  summarise(
    mean_temp = mean(temperature, na.rm = TRUE),
    max_temp  = max(temperature, na.rm = TRUE),
    min_temp  = min(temperature, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    anomaly_mean = mean_temp - mmm,
    anomaly_max  = max_temp  - mmm,
    hotspot_mean = if_else(mean_temp >= mmm + 1, mean_temp - mmm, 0),
    hotspot_max  = if_else(max_temp  >= mmm + 1, max_temp  - mmm, 0)
  )

# -----------------------------
# 7. In situ 2024 max/min summaries
# -----------------------------
insitu_2024 <- insitu_daily %>%
  filter(year(date_only) == 2024)

insitu_sites <- c(
  "coral_gardens", "sylphs_hole", "north_bay",
  "horseshoe_reef", "potholes_reef"
)

insitu_max_temp_2024 <- insitu_2024 %>%
  filter(site %in% insitu_sites) %>%
  group_by(site) %>%
  slice_max(max_temp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(site, max_temp, date_of_max = date_only)

insitu_min_temp_2024 <- insitu_2024 %>%
  filter(site %in% insitu_sites) %>%
  group_by(site) %>%
  slice_min(min_temp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(site, min_temp, date_of_min = date_only)

# Highest daily mean temperature across all sites

highest_mean_temp_by_site <- insitu_daily %>%
  group_by(site) %>%
  slice_max(mean_temp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(site, date_only, mean_temp)

highest_mean_temp_by_site

# -----------------------------
# 8. In situ DHW calculations
# Mean-DHW object uses daily mean temperature HotSpots
# Max-DHW object uses daily maximum temperature HotSpots
# -----------------------------
insitu_dhw_mean <- insitu_daily %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(
    dhw = calc_dhw(hotspot_mean, width = dhw_window_days, partial = TRUE),
    dhw = replace_na(dhw, 0)
  ) %>%
  ungroup()

insitu_dhw_max <- insitu_daily %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(
    dhw = calc_dhw(hotspot_max, width = dhw_window_days, partial = TRUE),
    dhw = replace_na(dhw, 0)
  ) %>%
  ungroup()

# Default in situ DHW for mean-condition analyses
insitu_dhw <- insitu_dhw_mean

# -----------------------------
# 9. In situ peak DHW summaries
# -----------------------------
peak_insitu_mean_dhw <- insitu_dhw_mean %>%
  group_by(site) %>%
  summarise(
    peak_mean_dhw = max(dhw, na.rm = TRUE),
    date_peak_mean_dhw = date_only[which.max(dhw)],
    .groups = "drop"
  )

peak_insitu_max_dhw <- insitu_dhw_max %>%
  group_by(site) %>%
  summarise(
    peak_max_dhw = max(dhw, na.rm = TRUE),
    date_peak_max_dhw = date_only[which.max(dhw)],
    .groups = "drop"
  )

peak_insitu_mean_dhw
peak_insitu_max_dhw

# -----------------------------
# 10. Satellite full time series preprocessing
# Do NOT filter to 2024 before calculating DHW.
# DHW needs the previous 84 days.
# -----------------------------
satellite_all <- read_csv("LHI_sst_locations_timeseries.csv") %>%
  mutate(date_only = dmy(date)) %>%
  pivot_longer(
    cols = where(is.numeric),
    names_to = "site",
    values_to = "sst"
  ) %>%
  mutate(
    site = str_to_lower(str_replace_all(site, " ", "_")),
    site = case_when(site == "stevens_hole" ~ "coral_gardens", TRUE ~ site)
  ) %>%
  left_join(mmms, by = "site") %>%
  filter(!is.na(mmm), site != "le_meurthe") %>%
  mutate(
    year = year(date_only),
    anomaly = sst - mmm,
    hotspot = if_else(sst >= mmm + 1, sst - mmm, 0),
    above_mmm = anomaly > 0
  )

# -----------------------------
# 11. Satellite DHW
# NOAA-style DHW calculated separately for each site
# -----------------------------
satellite_dhw_all <- satellite_all %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(
    dhw = calc_dhw(hotspot, width = dhw_window_days, partial = FALSE)
  ) %>%
  ungroup()

# 2024 satellite data after DHW is calculated
satellite <- satellite_all %>%
  filter(year == 2024, date_only >= ymd("2024-01-01"), date_only <= ymd("2024-05-31"))

satellite_dhw <- satellite_dhw_all %>%
  filter(year == 2024, date_only >= ymd("2024-01-01"), date_only <= ymd("2024-05-31"))

# -----------------------------
# 12. Satellite 2024 max/min SST summaries
# -----------------------------
satellite_max_temp_2024 <- satellite %>%
  group_by(site) %>%
  summarise(
    max_temp = max(sst, na.rm = TRUE),
    date_of_max = date_only[which.max(sst)],
    .groups = "drop"
  )

satellite_min_temp_2024 <- satellite %>%
  filter(date_only <= ymd("2024-04-30")) %>%
  group_by(site) %>%
  summarise(
    min_temp = min(sst, na.rm = TRUE),
    date_of_min = date_only[which.min(sst)],
    .groups = "drop"
  )

# -----------------------------
# 13. Satellite peak DHW summaries
# Peak satellite DHW = maximum site/pixel DHW
# Mean satellite DHW = regional average across sites
# -----------------------------
satellite_peak_dhw_duration <- satellite_dhw %>%
  group_by(site) %>%
  filter(dhw == max(dhw, na.rm = TRUE)) %>%
  summarise(
    peak_dhw = unique(dhw),
    start_date = min(date_only),
    end_date = max(date_only),
    duration_days = n(),
    .groups = "drop"
  )

years_of_interest <- c(1998, 2010, 2019, 2024)

peak_satellite_dhw_by_year <- satellite_dhw_all %>%
  filter(year %in% years_of_interest) %>%
  group_by(year) %>%
  slice_max(dhw, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(year, site, date_only, peak_satellite_dhw = dhw) %>%
  arrange(year)

peak_mean_satellite_dhw_by_year <- satellite_dhw_all %>%
  filter(year %in% years_of_interest) %>%
  group_by(year, date_only) %>%
  summarise(
    mean_sst = mean(sst, na.rm = TRUE),
    mean_dhw = mean(dhw, na.rm = TRUE),
    max_dhw = max(dhw, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  group_by(year) %>%
  slice_max(mean_dhw, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  arrange(year)

peak_satellite_dhw_by_year
peak_mean_satellite_dhw_by_year

##===============================
## FIGURE 1: MEAN In Situ & Satellite SST + DHW SUP FIGURES
##===============================

insitu_plot_df <- insitu_dhw_mean %>%
  mutate(site = if_else(site == "coral_gardens", "stephens_hole", site))

satellite_plot_df <- satellite_dhw %>%
  rename(mean_temp = sst) %>%
  scale_dhw(dhw_col = "dhw", sst_max = 30) %>%
  filter(date_only <= as.Date("2024-04-30")) %>%
  mutate(site = if_else(site == "coral_gardens", "stephens_hole", site))

max_insitu_dhw_mean <- max(insitu_plot_df$dhw, na.rm = TRUE)

insitu_plot_df <- insitu_plot_df %>%
  mutate(
    dhw_scaled  = dhw / max_insitu_dhw_mean * 30,
    dhw4_scaled = 4 / max_insitu_dhw_mean * 30,
    dhw8_scaled = 8 / max_insitu_dhw_mean * 30
  )

make_labels_fixed <- function(df) {
  df %>%
    group_by(site) %>%
    summarise(
      x_label = as.Date("2024-04-30") - 10,
      mmm_y   = max(mmm, na.rm = TRUE),
      mmm1_y  = max(mmm + 1, na.rm = TRUE),
      dhw4_y  = 4,
      dhw8_y  = 8,
      .groups = "drop"
    )
}

insitu_labels_fixed <- make_labels_fixed(insitu_plot_df) %>%
  mutate(
    mmm_y  = mmm_y  - 0.5,
    mmm1_y = mmm1_y + 0.5,
    dhw4_y = dhw4_y - 0.3,
    dhw8_y = dhw8_y + 0.3
  )

satellite_labels_fixed <- make_labels_fixed(satellite_plot_df) %>%
  mutate(
    mmm_y  = mmm_y  - 0.5,
    mmm1_y = mmm1_y + 0.5,
    dhw4_y = dhw4_y - 0.3,
    dhw8_y = dhw8_y + 0.3
  )

fig_insitu_mean <- ggplot(insitu_plot_df, aes(x = date_only)) +
  geom_line(aes(y = mean_temp), colour = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm), colour = "pink", linetype = "dashed") +
  geom_line(aes(y = mmm + 1), colour = "red", linetype = "dashed") +
  geom_rect(
    data = insitu_plot_df %>% filter(hotspot_mean >= 1),
    aes(xmin = date_only - 0.5, xmax = date_only + 0.5, ymin = -Inf, ymax = Inf),
    fill = "orange",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_line(aes(y = dhw_scaled), colour = "blue", linewidth = 1) +
  geom_hline(aes(yintercept = dhw4_scaled), colour = "lightblue", linetype = "dashed") +
  geom_hline(aes(yintercept = dhw8_scaled), colour = "blue", linetype = "dashed") +
  geom_text(data = insitu_labels_fixed, aes(x = x_label, y = mmm_y, label = "MMM"),
            colour = "pink", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels_fixed, aes(x = x_label, y = mmm1_y, label = "MMM + 1°C"),
            colour = "red", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels_fixed, aes(x = x_label, y = dhw4_y, label = "eDHW 4"),
            colour = "lightblue", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels_fixed, aes(x = x_label, y = dhw8_y, label = "eDHW 8"),
            colour = "blue", hjust = 1.1, size = 3, fontface = "bold") +
  facet_wrap(~ site, scales = "free_y") +
  scale_x_date(
    limits = c(as.Date("2024-01-01"), max(insitu_plot_df$date_only)),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    name = "IST (°C)",
    sec.axis = sec_axis(
      trans = ~ . * (max_insitu_dhw_mean / 30),
      name = "eDHW (°C-weeks)"
    )
  ) +
  coord_cartesian(ylim = c(NA, 30)) +
  labs(
    x = "Date",
    title = expression(italic("In situ") ~ "daily mean IST and experimental Degree Heating Weeks")
  ) +
  theme_noaa()

fig_satellite_mean <- ggplot(satellite_plot_df, aes(x = date_only)) +
  geom_line(aes(y = mean_temp), colour = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm), colour = "pink", linetype = "dashed") +
  geom_line(aes(y = mmm + 1), colour = "red", linetype = "dashed") +
  geom_rect(
    data = satellite_plot_df %>% filter(hotspot >= 1),
    aes(xmin = date_only - 0.5, xmax = date_only + 0.5, ymin = -Inf, ymax = Inf),
    fill = "orange",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_line(aes(y = dhw_scaled), colour = "blue", linewidth = 1, na.rm = TRUE) +
  geom_hline(aes(yintercept = dhw4_scaled), colour = "lightblue", linetype = "dashed") +
  geom_hline(aes(yintercept = dhw8_scaled), colour = "blue", linetype = "dashed") +
  geom_text(data = satellite_labels_fixed, aes(x = x_label, y = mmm_y, label = "MMM"),
            colour = "pink", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels_fixed, aes(x = x_label, y = mmm1_y, label = "MMM + 1°C"),
            colour = "red", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels_fixed, aes(x = x_label, y = dhw4_y, label = "eDHW 4"),
            colour = "lightblue", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels_fixed, aes(x = x_label, y = dhw8_y, label = "eDHW 8"),
            colour = "blue", hjust = 1.1, size = 3, fontface = "bold") +
  facet_wrap(~ site, scales = "free_y") +
  scale_x_date(
    limits = c(as.Date("2024-01-01"), as.Date("2024-04-30")),
    breaks = as.Date(c("2024-01-01", "2024-02-01", "2024-03-01", "2024-04-01")),
    labels = c("Jan", "Feb", "Mar", "Apr"),
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    name = "SST (°C)",
    sec.axis = sec_axis(
      transform = ~ . * (20 / 30),
      name = "DHWs (°C-weeks)"
    )
  ) +
  coord_cartesian(ylim = c(NA, 30)) +
  labs(
    x = "Date",
    title = "Satellite SST and DHWs"
  ) +
  theme_noaa() +
  theme(plot.title = element_text(face = "plain", size = 12))

fig_combined_mean <- fig_insitu_mean / fig_satellite_mean + plot_layout(heights = c(1, 1))
fig_combined_mean

##===============================
## FIGURE 2: MAX In Situ & Satellite SST + DHW SUP FIGURES
##===============================

insitu_plot_df <- insitu_dhw_max %>%
  mutate(site = if_else(site == "coral_gardens", "stephens_hole", site))

satellite_plot_df <- satellite_dhw %>%
  rename(max_temp = sst) %>%
  scale_dhw(dhw_col = "dhw", sst_max = 30) %>%
  filter(date_only <= as.Date("2024-04-30")) %>%
  mutate(site = if_else(site == "coral_gardens", "stephens_hole", site))

max_insitu_dhw <- max(insitu_plot_df$dhw, na.rm = TRUE)

insitu_plot_df <- insitu_plot_df %>%
  mutate(
    dhw_scaled  = dhw / max_insitu_dhw * 30,
    dhw4_scaled = 4 / max_insitu_dhw * 30,
    dhw8_scaled = 8 / max_insitu_dhw * 30
  )

insitu_labels_fixed <- make_labels_fixed(insitu_plot_df) %>%
  mutate(
    mmm_y  = mmm_y  - 0.5,
    mmm1_y = mmm1_y + 0.5,
    dhw4_y = dhw4_y - 0.3,
    dhw8_y = dhw8_y + 0.3
  )

satellite_labels_fixed <- make_labels_fixed(satellite_plot_df) %>%
  mutate(
    mmm_y  = mmm_y  - 0.5,
    mmm1_y = mmm1_y + 0.5,
    dhw4_y = dhw4_y - 0.3,
    dhw8_y = dhw8_y + 0.3
  )

fig_insitu_max <- ggplot(insitu_plot_df, aes(x = date_only)) +
  geom_line(aes(y = max_temp), colour = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm), colour = "pink", linetype = "dashed") +
  geom_line(aes(y = mmm + 1), colour = "red", linetype = "dashed") +
  geom_rect(
    data = insitu_plot_df %>% filter(hotspot_max >= 1),
    aes(xmin = date_only - 0.5, xmax = date_only + 0.5, ymin = -Inf, ymax = Inf),
    fill = "orange",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_line(aes(y = dhw_scaled), colour = "blue", linewidth = 1) +
  geom_hline(aes(yintercept = dhw4_scaled), colour = "lightblue", linetype = "dashed") +
  geom_hline(aes(yintercept = dhw8_scaled), colour = "blue", linetype = "dashed") +
  geom_text(data = insitu_labels_fixed, aes(x = x_label, y = mmm_y, label = "MMM"),
            colour = "pink", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels_fixed, aes(x = x_label, y = mmm1_y, label = "MMM + 1°C"),
            colour = "red", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels_fixed, aes(x = x_label, y = dhw4_y, label = "eDHW 4"),
            colour = "lightblue", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels_fixed, aes(x = x_label, y = dhw8_y, label = "eDHW 8"),
            colour = "blue", hjust = 1.1, size = 3, fontface = "bold") +
  facet_wrap(~ site, scales = "free_y") +
  scale_x_date(
    limits = c(as.Date("2024-01-01"), max(insitu_plot_df$date_only)),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  scale_y_continuous(
    name = "IST (°C)",
    sec.axis = sec_axis(
      trans = ~ . * (max_insitu_dhw / 30),
      name = "eDHW (°C-weeks)"
    )
  ) +
  coord_cartesian(ylim = c(NA, 30)) +
  labs(
    x = "Date",
    title = expression(italic("In situ") ~ "daily maximum IST and experimental Degree Heating Weeks")
  ) +
  theme_noaa()

fig_satellite_max <- ggplot(satellite_plot_df, aes(x = date_only)) +
  geom_line(aes(y = max_temp), colour = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm), colour = "pink", linetype = "dashed") +
  geom_line(aes(y = mmm + 1), colour = "red", linetype = "dashed") +
  geom_rect(
    data = satellite_plot_df %>% filter(hotspot >= 1),
    aes(xmin = date_only - 0.5, xmax = date_only + 0.5, ymin = -Inf, ymax = Inf),
    fill = "orange",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  geom_line(aes(y = dhw_scaled), colour = "blue", linewidth = 1, na.rm = TRUE) +
  geom_hline(aes(yintercept = dhw4_scaled), colour = "lightblue", linetype = "dashed") +
  geom_hline(aes(yintercept = dhw8_scaled), colour = "blue", linetype = "dashed") +
  geom_text(data = satellite_labels_fixed, aes(x = x_label, y = mmm_y, label = "MMM"),
            colour = "pink", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels_fixed, aes(x = x_label, y = mmm1_y, label = "MMM + 1°C"),
            colour = "red", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels_fixed, aes(x = x_label, y = dhw4_y, label = "eDHW 4"),
            colour = "lightblue", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels_fixed, aes(x = x_label, y = dhw8_y, label = "eDHW 8"),
            colour = "blue", hjust = 1.1, size = 3, fontface = "bold") +
  facet_wrap(~ site, scales = "free_y") +
  scale_x_date(
    limits = c(as.Date("2024-01-01"), as.Date("2024-04-30")),
    breaks = as.Date(c("2024-01-01", "2024-02-01", "2024-03-01", "2024-04-01")),
    labels = c("Jan", "Feb", "Mar", "Apr"),
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    name = "SST (°C)",
    sec.axis = sec_axis(
      transform = ~ . * (20 / 30),
      name = "DHWs (°C-weeks)"
    )
  ) +
  coord_cartesian(ylim = c(NA, 30)) +
  labs(
    x = "Date",
    title = "Satellite SST and DHWs"
  ) +
  theme_noaa() +
  theme(plot.title = element_text(face = "plain", size = 12))

fig_combined_max <- 
  (fig_insitu_max +
     labs(tag = "a") +
     theme(
       plot.tag = element_text(size = 14, face = "bold"),
       plot.tag.position = c(0.01, 0.98)
     )) /
  (fig_satellite_max +
     labs(tag = "b") +
     theme(
       plot.tag = element_text(size = 14, face = "bold"),
       plot.tag.position = c(0.01, 0.98)
     )) +
  plot_layout(heights = c(1, 1))

fig_combined_max

# ggsave(
#   "fig_combined_max.png",
#   fig_combined_max,
#   width = 12,
#   height = 15,
#   dpi = 300,
#   bg = "white"
# )

##===============================
## FIGURE 3: Combined satellite mean SST + mean DHW, 2024
##===============================

combined_df <- satellite_dhw %>%
  group_by(date_only) %>%
  summarise(
    mean_temp  = mean(sst, na.rm = TRUE),
    mmm        = mean(mmm, na.rm = TRUE),
    mean_dhw   = mean(dhw, na.rm = TRUE),
    hotspot    = max(hotspot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dhw_scaled = mean_dhw,
    dhw4_scaled = 4,
    dhw8_scaled = 8
  )

combined_df_full <- combined_df %>%
  complete(date_only = seq(as.Date("2024-01-01"), as.Date("2024-05-31"), by = "day")) %>%
  fill(mean_temp, mmm, mean_dhw, dhw_scaled, dhw4_scaled, dhw8_scaled, hotspot, .direction = "down")

fig_combined_sites <- ggplot(combined_df_full, aes(x = date_only)) +
  geom_line(aes(y = mean_temp), colour = "grey50", linewidth = 0.8) +
  geom_line(aes(y = mmm), colour = "#D61C59", linetype = "dashed") +
  geom_line(aes(y = mmm + 1), colour = "red", linetype = "dashed") +
  geom_line(aes(y = dhw_scaled), colour = "black", linewidth = 1) +
  geom_hline(aes(yintercept = dhw4_scaled), colour = "#5DADE2", linetype = "dashed") +
  geom_hline(aes(yintercept = dhw8_scaled), colour = "#2166AC", linetype = "dashed") +
  annotate("text", x = as.Date("2024-05-25"), y = max(combined_df$mmm) - 0.5,
           label = "MMM", colour = "#D61C59", hjust = 1.1, size = 4, fontface = "bold") +
  annotate("text", x = as.Date("2024-05-25"), y = max(combined_df$mmm) + 1.5,
           label = "MMM + 1°C", colour = "red", hjust = 1.1, size = 4, fontface = "bold") +
  annotate("text", x = as.Date("2024-05-25"), y = 3.5,
           label = "DHW 4°C", colour = "#5DADE2", hjust = 1.1, size = 4, fontface = "bold") +
  annotate("text", x = as.Date("2024-05-25"), y = 7.5,
           label = "DHW 8°C", colour = "#2166AC", hjust = 1.1, size = 4, fontface = "bold") +
  scale_x_date(
    limits = c(as.Date("2024-01-01"), as.Date("2024-05-31")),
    breaks = as.Date(c("2024-01-01","2024-02-01","2024-03-01","2024-04-01","2024-05-01")),
    labels = c("Jan 2024","Feb 2024","Mar 2024","Apr 2024","May 2024"),
    expand = expansion(add = 0)
  ) +
  scale_y_continuous(
    name = "SST (°C)",
    sec.axis = sec_axis(
      trans = ~ .,
      name = "DHWs (°C-weeks)"
    )
  ) +
  coord_cartesian(ylim = c(NA, 30)) +
  labs(x = "Month") +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1, margin = margin(t = 5)),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 12, face = "plain")
  )

season_df <- data.frame(
  date_only = seq(as.Date("2024-01-01"), as.Date("2024-05-31"), by = "day")
) %>%
  mutate(
    doy = yday(date_only),
    season_value = 1 - (doy - 1) / (151 - 1) * 0.5
  )

fig_combined_sites_season <- fig_combined_sites +
  geom_tile(
    data = season_df,
    aes(x = date_only, y = -2, fill = season_value),
    height = 1.5,
    inherit.aes = FALSE
  ) +
  scale_fill_gradientn(
    colours = c("#053061", "#2166AC", "#E34A33", "#67001F"),
    values = scales::rescale(c(0, 0.5, 0.8, 1)),
    limits = c(0, 1),
    guide = "none"
  ) +
  annotate("text", x = as.Date("2024-01-15"), y = -0.8,
           label = "Summer", colour = "#67001F", size = 4.2, fontface = "bold") +
  annotate("text", x = as.Date("2024-03-15"), y = -0.8,
           label = "Autumn", colour = "#E34A33", size = 4.2, fontface = "bold") +
  annotate("text", x = as.Date("2024-05-15"), y = -0.8,
           label = "Winter", colour = "#053061", size = 4.2, fontface = "bold") +
  coord_cartesian(ylim = c(-2, 30)) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    plot.title = element_text(size = 12)
  )

fig_combined_sites_season

# ggsave(
#   "fig_combined_sites_season.png",
#   fig_combined_sites_season,
#   width = 10,
#   height = 6,
#   dpi = 300,
#   bg = "white"
# )

##===============================
## RESULTS FIGURE: Mean in situ sites + mean satellite conditions
##===============================

plot_start_sst <- insitu_start
plot_end_sst   <- insitu_end
plot_start_dhw <- insitu_start
plot_end_dhw   <- insitu_end

insitu_site_map <- c(
  "coral_gardens"  = "Stephens Hole",
  "sylphs_hole"    = "Sylphs Hole",
  "north_bay"      = "North Bay",
  "horseshoe_reef" = "Horseshoe Reef",
  "potholes_reef"  = "Potholes"
)

site_levels <- c(
  "Satellite",
  "Stephens Hole",
  "Sylphs Hole",
  "North Bay",
  "Horseshoe Reef",
  "Potholes"
)

site_colors <- c(
  "Satellite"      = "black",
  "Stephens Hole"  = "#B54BFC",
  "Sylphs Hole"    = "#4378E9",
  "North Bay"      = "#FA5454",
  "Horseshoe Reef" = "#29DABA",
  "Potholes"       = "#EC9B31"
)

site_linewidths <- c(
  "Satellite"      = 1.2,
  "Stephens Hole"  = 0.8,
  "Sylphs Hole"    = 0.8,
  "North Bay"      = 0.8,
  "Horseshoe Reef" = 0.8,
  "Potholes"       = 0.8
)

insitu_sst_plot_data <- insitu_daily %>%
  filter(
    site %in% names(insitu_site_map),
    date_only >= plot_start_sst,
    date_only <= plot_end_sst
  ) %>%
  transmute(site, date_only, temp = mean_temp)

satellite_sst_plot_data <- satellite %>%
  filter(date_only >= plot_start_sst, date_only <= plot_end_sst) %>%
  group_by(date_only) %>%
  summarise(temp = mean(sst, na.rm = TRUE), .groups = "drop") %>%
  mutate(site = "Satellite")

sst_combined <- bind_rows(
  insitu_sst_plot_data,
  satellite_sst_plot_data
) %>%
  mutate(
    site = dplyr::recode(site, !!!insitu_site_map, "Satellite" = "Satellite"),
    site = factor(site, levels = site_levels)
  )

insitu_dhw_plot_data <- insitu_dhw_mean %>%
  filter(
    site %in% names(insitu_site_map),
    date_only >= plot_start_dhw,
    date_only <= plot_end_dhw
  )

satellite_dhw_mean <- satellite_dhw %>%
  filter(date_only >= plot_start_dhw, date_only <= plot_end_dhw) %>%
  group_by(date_only) %>%
  summarise(dhw = mean(dhw, na.rm = TRUE), .groups = "drop") %>%
  mutate(site = "Satellite")

dhw_combined <- bind_rows(
  insitu_dhw_plot_data %>% select(site, date_only, dhw),
  satellite_dhw_mean
) %>%
  mutate(
    site = dplyr::recode(site, !!!insitu_site_map, "Satellite" = "Satellite"),
    site = factor(site, levels = site_levels)
  )

fig_sst_filtered <- ggplot(sst_combined, aes(x = date_only, y = temp, colour = site)) +
  geom_line(aes(size = site), na.rm = TRUE) +
  geom_hline(yintercept = 24.06, colour = "black", linetype = "dashed", linewidth = 0.8) +
  annotate("text", x = as.Date("2024-04-28"), y = 24.06 + 0.15,
           label = "MMM", colour = "black", hjust = 1, size = 4.2) +
  scale_color_manual(values = site_colors, breaks = site_levels) +
  scale_size_manual(values = site_linewidths, guide = "none") +
  scale_x_date(
    breaks = seq(as.Date("2023-12-01"), as.Date("2024-05-01"), by = "1 month"),
    labels = label_date("%b %Y"),
    limits = c(as.Date("2023-12-01"), as.Date("2024-05-01"))
  ) +
  labs(x = "Month", y = "Daily mean SST (°C)", colour = "Site / Source") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "none"
  )

fig_dhw_filtered <- ggplot(dhw_combined, aes(x = date_only, y = dhw, colour = site)) +
  geom_line(aes(size = site), na.rm = TRUE) +
  scale_color_manual(values = site_colors, breaks = site_levels) +
  scale_size_manual(values = site_linewidths, guide = "none") +
  scale_x_date(
    breaks = seq(as.Date("2023-12-01"), as.Date("2024-05-01"), by = "1 month"),
    labels = label_date("%b %Y"),
    limits = c(as.Date("2023-12-01"), as.Date("2024-05-01"))
  ) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "Month", y = "Degree Heating Weeks (°C-weeks)") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.position = "none"
  )

shared_legend <- get_legend(
  ggplot(sst_combined, aes(x = date_only, y = temp, colour = site)) +
    geom_line(size = 1, na.rm = TRUE) +
    scale_color_manual(values = site_colors, breaks = site_levels) +
    labs(colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.2))) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12),
      legend.title = element_text(size = 12),
      legend.key = element_blank()
    )
)

fig_sst_combined <- plot_grid(
  shared_legend,
  plot_grid(fig_sst_filtered, fig_dhw_filtered, ncol = 2, align = "hv"),
  ncol = 1,
  rel_heights = c(0.1, 1)
)

fig_sst_combined

##===============================
## PAST BLEACHING EVENTS RESULTS FIGURE
## Mean SST and mean DHW across satellite sites
##===============================

df_sst <- satellite_dhw_all %>%
  select(date_only, site, sst, dhw, mmm) %>%
  rename(sdhw = dhw)

years_of_interest <- c(1998, 2010, 2019, 2024)

df_sst_filtered <- df_sst %>%
  filter(year(date_only) %in% years_of_interest, site != "le_meurthe") %>%
  mutate(
    year = year(date_only),
    day_of_year = yday(date_only)
  ) %>%
  filter(month(date_only) <= 7)

month_starts <- c(1, 32, 60, 91, 121, 152, 182)
month_labels <- month.abb[1:7]
MMM <- 24.06

plot_sst_year <- df_sst_filtered %>%
  group_by(year, day_of_year) %>%
  summarise(sst = mean(sst, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_year, y = sst, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = MMM, linetype = "dotted", color = "black", size = 0.8) +
  annotate("text", x = 5, y = MMM + 0.15, label = "MMM", color = "black", hjust = 0) +
  labs(title = "Satellite SST", x = "Month", y = "SST (°C)", color = "Year") +
  scale_x_continuous(breaks = month_starts, labels = month_labels) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5)
  )

plot_dhw_mean <- df_sst_filtered %>%
  group_by(year, day_of_year) %>%
  summarise(dhw = mean(sdhw, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_year, y = dhw, color = factor(year))) +
  geom_line(size = 1) +
  labs(title = "Satellite DHW", x = "Month", y = "DHW (°C-weeks)", color = "Year") +
  scale_x_continuous(breaks = month_starts, labels = month_labels) +
  theme_classic(base_size = 12) +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12)
  )

combined_mean <- plot_sst_year + plot_dhw_mean +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

combined_mean

mean_dhw_summary <- df_sst_filtered %>%
  group_by(year, day_of_year) %>%
  summarise(mean_sdhw = mean(sdhw, na.rm = TRUE), .groups = "drop_last") %>%
  slice_max(mean_sdhw, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(year, peak_day_of_year = day_of_year, peak_mean_sdhw = mean_sdhw)

mean_dhw_summary

##===============================
## Final combined figure
##===============================

year_colors <- c(
  "1998" = "#A6CEE3",
  "2010" = "#B2DF8A",
  "2019" = "#FB9A99",
  "2024" = "#FDBF6F"
)

month_starts_sst <- c(1,32,60,91,121,152,182,213,244,274,305,335)

plot_sst_year <- df_sst %>%
  filter(year(date_only) %in% years_of_interest, site != "le_meurthe") %>%
  mutate(
    year = year(date_only),
    day_of_year = yday(date_only)
  ) %>%
  group_by(year, day_of_year) %>%
  summarise(sst = mean(sst, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_year, y = sst, color = factor(year))) +
  geom_line(size = 1.2) +
  geom_hline(
    yintercept = 24.06,
    colour = "black",
    linetype = "dashed",
    linewidth = 0.8
  ) +
  annotate(
    "text",
    x = 355,
    y = MMM + 0.15,
    label = "MMM",
    color = "black",
    hjust = 1,
    size = 4.2
  ) +
  scale_x_continuous(
    breaks = month_starts_sst,
    labels = paste(month.abb, "2024")
  ) +
  scale_color_manual(values = year_colors) +
  labs(title = "Satellite SST", x = "Month", y = "SST (°C)", color = "Year") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5, vjust = -14)
  )

plot_dhw_mean_dec_jul <- df_sst %>%
  filter(year(date_only) %in% years_of_interest, site != "le_meurthe") %>%
  mutate(
    year = year(date_only),
    day_of_year_adj = yday(date_only) + ifelse(month(date_only) == 12, -365, 0)
  ) %>%
  filter(month(date_only) == 12 | month(date_only) <= 7) %>%
  group_by(year, day_of_year_adj) %>%
  summarise(dhw = max(sdhw, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_year_adj, y = dhw, color = factor(year))) +
  geom_line(size = 1.2) +
  scale_x_continuous(
    breaks = c(-30, 1, 32, 60, 91, 121, 152, 182),
    labels = c(
      "Dec 2023", "Jan 2024", "Feb 2024", "Mar 2024",
      "Apr 2024", "May 2024", "Jun 2024", "Jul 2024"
    )
  ) +
  scale_color_manual(values = year_colors) +
  labs(title = "Satellite DHWs", x = "Month", y = "DHW (°C-weeks)", color = "Year") +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(size = 12, angle = 45, hjust = 1),
    axis.text.y = element_text(size = 12),
    axis.title.x = element_text(size = 12),
    axis.title.y = element_text(size = 12),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 12),
    plot.title = element_text(size = 12, hjust = 0.5, vjust = -14)
  )

fig_updated_combined_mean <-
  (plot_sst_year + labs(x = NULL)) +
  (plot_dhw_mean_dec_jul + labs(x = NULL)) +
  plot_layout(ncol = 2, guides = "collect") &
  theme(legend.position = "top")

fig_sst_combined <- plot_grid(
  shared_legend,
  plot_grid(
    fig_sst_filtered + labs(x = NULL),
    fig_dhw_filtered + labs(x = NULL),
    ncol = 2,
    align = "hv"
  ),
  ncol = 1,
  rel_heights = c(0.1, 1)
)

top_fig <-
  ggdraw(fig_sst_combined) +
  draw_plot_label(
    label = c("a", "b"),
    x = c(0.02, 0.52),
    y = c(0.92, 0.92),
    size = 16,
    fontface = "bold"
  ) +
  draw_label(
    "Month",
    x = 0.5,
    y = 0.015,
    size = 12
  )

bottom_fig <-
  ggdraw(fig_updated_combined_mean) +
  draw_plot_label(
    label = c("c", "d"),
    x = c(0.02, 0.52),
    y = c(0.88, 0.88),
    size = 16,
    fontface = "bold"
  ) +
  draw_label(
    "Month",
    x = 0.5,
    y = 0.015,
    size = 12
  )

final_figure <-
  plot_grid(
    top_fig,
    bottom_fig,
    ncol = 1,
    rel_heights = c(1, 1),
    align = "v"
  )

final_figure

# ggsave(
#   "final_figure.png",
#   final_figure,
#   width = 12,
#   height = 15,
#   dpi = 300,
#   bg = "white"
# )

# -----------------------------
# Yearly satellite heat-stress summary
# 1998, 2010, 2019, 2024
# Uses max site/pixel values per day for NOAA-style peak exposure
# -----------------------------

heat_stress_summary_by_year <- satellite_dhw_all %>%
  filter(year %in% c(1998, 2010, 2019, 2024)) %>%
  group_by(year, date_only) %>%
  summarise(
    max_sst = max(sst, na.rm = TRUE),
    max_mmm = max(mmm, na.rm = TRUE),
    max_hotspot = max(hotspot, na.rm = TRUE),
    mean_hotspot = mean(hotspot, na.rm = TRUE),
    max_dhw = max(dhw, na.rm = TRUE),
    
    above_mmm = any(sst > mmm, na.rm = TRUE),
    above_mmm_plus_1 = any(sst >= mmm + 1, na.rm = TRUE),
    hotspot_day = any(hotspot > 0, na.rm = TRUE),
    dhw_over_4 = any(dhw >= 4, na.rm = TRUE),
    
    .groups = "drop"
  ) %>%
  group_by(year) %>%
  summarise(
    days_above_mmm = sum(above_mmm, na.rm = TRUE),
    weeks_above_mmm = days_above_mmm / 7,
    
    days_above_mmm_plus_1 = sum(above_mmm_plus_1, na.rm = TRUE),
    weeks_above_mmm_plus_1 = days_above_mmm_plus_1 / 7,
    
    hotspot_days = sum(hotspot_day, na.rm = TRUE),
    hotspot_weeks = hotspot_days / 7,
    
    days_dhw_over_4 = sum(dhw_over_4, na.rm = TRUE),
    weeks_dhw_over_4 = days_dhw_over_4 / 7,
    
    peak_dhw = max(max_dhw, na.rm = TRUE),
    date_peak_dhw = date_only[which.max(max_dhw)],
    
    peak_mean_hotspot = max(mean_hotspot, na.rm = TRUE),
    mean_hotspot_across_event = mean(mean_hotspot, na.rm = TRUE),
    
    .groups = "drop"
  )

heat_stress_summary_by_year
