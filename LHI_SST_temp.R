##===============================
## 2024 MARINE HEATWAVE ANALYSIS – Lord Howe Island
## Max SST / In Situ & Satellite Data
##===============================

# -----------------------------
# 0️⃣ Clean environment
# -----------------------------
rm(list = ls())

# -----------------------------
# 1️⃣ Load required packages
# -----------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lubridate, janitor, zoo, data.table, 
               ggplot2, patchwork, viridis, scales, cowplot, paletteer)

# -----------------------------
# 2️⃣ Define paths & global parameters
# -----------------------------

insitu_start <- ymd("2023-12-01")
insitu_end   <- ymd("2024-05-31")
sat_years    <- c(1998, 2010, 2019, 2024)
dhw_window_days <- 84

# -----------------------------
# 3️⃣ Load site MMMs
# -----------------------------
mmms <- read_csv("LHI_mmms.csv") %>%
  clean_names() %>%
  mutate(site = str_to_lower(str_replace_all(site, " ", "_")))

# -----------------------------
# 4️⃣ Load raw in situ logger data
# -----------------------------
insitu_raw <- read_csv("LHIMP_DATA_in_situ_loggers.csv") %>%
  clean_names() %>%
  filter(date_time != "") %>%
  mutate(date_time = dmy_hms(paste0(date_time, ":00")),
         date_only = as_date(date_time)) %>%
  mutate(across(where(is.numeric), ~na_if(., 0))) # replace 0 with NA

# -----------------------------
# 5️⃣ Convert in situ data to long format
# -----------------------------
insitu_long <- insitu_raw %>%
  pivot_longer(observatory_rock:coral_gardens,
               names_to = "site",
               values_to = "temperature") %>%
  filter(!is.na(temperature)) %>%
  mutate(site = str_to_lower(str_replace_all(site, " ", "_"))) %>%
  left_join(mmms, by = "site") %>%
  filter(!is.na(mmm)) %>%           # only valid sites
  filter(date_only >= insitu_start & date_only <= insitu_end) %>% 
  filter(site != "comets_hole")     # remove faulty logger

# -----------------------------
# 6️⃣ Daily max & min temperatures
# -----------------------------
insitu_daily <- insitu_long %>%
  group_by(site, date_only, mmm) %>%
  summarise(
    max_temp = max(temperature, na.rm = TRUE),
    min_temp = min(temperature, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    anomaly = max_temp - mmm,
    hotspot = if_else(anomaly >= 1, anomaly, 0)
  )

# -----------------------------
# 7️⃣ Extract 2024 data
# -----------------------------
insitu_2024 <- insitu_daily %>% filter(year(date_only) == 2024)

# Max and min per site (separate tables)
max_temp_sites <- c("coral_gardens", "sylphs_hole", "north_bay",
                    "horseshoe_reef", "potholes_reef")

insitu_max_temp_2024 <- insitu_2024 %>%
  filter(site %in% max_temp_sites) %>%
  group_by(site) %>%
  summarise(
    max_temp = max(max_temp, na.rm = TRUE),
    date_of_max = date_only[which.max(max_temp)],
    .groups = "drop"
  )

insitu_min_temp_2024 <- insitu_2024 %>%
  group_by(site) %>%
  summarise(
    min_temp = min(max_temp, na.rm = TRUE),
    date_of_min = date_only[which.min(max_temp)],
    .groups = "drop"
  )

# -----------------------------
# 8️⃣ Compute Degree Heating Weeks (DHW)
# -----------------------------
insitu_dhw <- insitu_2024 %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(dhw = zoo::rollapply(hotspot,
                              width = dhw_window_days,
                              FUN = function(x) sum(x)/7,
                              fill = NA,
                              align = "right")) %>%
  ungroup()

# -----------------------------
# 9️⃣ Peak DHW and consecutive runs
# -----------------------------
peak_insitu_dhw <- insitu_dhw %>%
  group_by(site) %>%
  filter(!is.na(dhw)) %>%
  summarise(
    peak_dhw = max(dhw, na.rm = TRUE),
    longest_run_days = {
      r <- rle(dhw %>% dplyr::near(max(dhw, na.rm = TRUE)))
      if(any(r$values)) max(r$lengths[r$values]) else 0
    },
    start_date = {
      r <- rle(dhw %>% dplyr::near(max(dhw, na.rm = TRUE)))
      idx <- which(r$values & r$lengths == max(r$lengths[r$values]))
      if(length(idx) == 0) NA_Date_ else {
        idx <- idx[1]
        start_idx <- if(idx > 1) sum(r$lengths[1:(idx-1)]) + 1 else 1
        date_only[start_idx]
      }
    },
    end_date = {
      r <- rle(dhw %>% dplyr::near(max(dhw, na.rm = TRUE)))
      idx <- which(r$values & r$lengths == max(r$lengths[r$values]))
      if(length(idx) == 0) NA_Date_ else {
        idx <- idx[1]
        start_idx <- if(idx > 1) sum(r$lengths[1:(idx-1)]) + 1 else 1
        end_idx <- start_idx + r$lengths[idx] - 1
        date_only[end_idx]
      }
    },
    .groups = "drop"
  )

# -----------------------------
# 1️⃣0️⃣ Round DHW for plotting consecutive runs
# -----------------------------
insitu_dhw_rounded <- insitu_dhw %>%
  filter(!is.na(dhw)) %>%
  mutate(dhw_rounded = round(dhw)) %>%
  group_by(site) %>%
  mutate(run_id = data.table::rleid(dhw_rounded)) %>%
  group_by(site, dhw_rounded, run_id) %>%
  summarise(
    start_date = min(date_only),
    end_date   = max(date_only),
    run_length = n(),
    .groups = "drop"
  )

# -----------------------------
# 1️⃣1️⃣ Satellite Data Preprocessing
# -----------------------------
satellite <- read_csv("LHI_sst_locations_timeseries.csv") %>%
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
    hotspot = if_else(anomaly >= 1, anomaly, 0),
    above_mmm = anomaly > 0
  ) %>%
  filter(year == 2024, date_only >= ymd("2024-01-01"), date_only <= ymd("2024-05-31"))

# Max & min SST separate
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
# 1️⃣2️⃣ Satellite DHW
# -----------------------------
satellite_dhw <- satellite %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(dhw = zoo::rollapply(hotspot,
                              width = dhw_window_days,
                              FUN = function(x) sum(x)/7,
                              fill = NA,
                              align = "right")) %>%
  ungroup()

# Peak DHW for satellite
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

##===============================
## PLOTTING: In Situ & Satellite SST + DHW
##===============================

library(ggplot2)
library(dplyr)
library(zoo)
library(patchwork)

# -----------------------------
# 0️⃣ NOAA-style theme
# -----------------------------
theme_noaa <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(size = base_size - 2)
    )
}

# -----------------------------
# 1️⃣ Scale DHW to match SST for plotting
# -----------------------------
scale_dhw <- function(df, dhw_col = "dhw", sst_col = "max_temp", dhw_max_plot = 30) {
  df <- df %>%
    group_by(site) %>%
    mutate(
      scale_factor = max(.data[[sst_col]], na.rm = TRUE) / dhw_max_plot,
      dhw_scaled = .data[[dhw_col]] * scale_factor,
      dhw4_scaled = 4 * scale_factor,
      dhw8_scaled = 8 * scale_factor
    ) %>%
    ungroup()
  return(df)
}

insitu_plot_df    <- scale_dhw(insitu_dhw)
satellite_plot_df <- scale_dhw(satellite_dhw)

# -----------------------------
# 2️⃣ Labels for right-hand side of plot
# -----------------------------
make_labels <- function(df, dhw_col_scaled = "dhw_scaled") {
  df %>%
    group_by(site) %>%
    summarise(
      x_label = max(date_only) - 25,
      mmm_end = max(mmm, na.rm = TRUE),
      mmm1_end = max(mmm + 1, na.rm = TRUE),
      dhw4_end = max(.data[[paste0("dhw4_scaled")]], na.rm = TRUE),
      dhw8_end = max(.data[[paste0("dhw8_scaled")]], na.rm = TRUE),
      .groups = "drop"
    )
}

insitu_labels    <- make_labels(insitu_plot_df)
satellite_labels <- make_labels(satellite_plot_df)

# -----------------------------
# 3️⃣ In Situ SST + DHW plot
# -----------------------------
fig_insitu <- ggplot(insitu_plot_df, aes(x = date_only)) +
  
  # SST lines
  geom_line(aes(y = max_temp), color = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm + 1), color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = mmm), color = "pink", linetype = "dashed", linewidth = 0.8) +
  
  # DHW scaled
  geom_line(aes(y = dhw_scaled), color = "blue", linewidth = 0.8) +
  
  # HotSpot shading
  geom_rect(
    data = insitu_plot_df %>% filter(hotspot >= 1),
    aes(xmin = date_only - 0.5, xmax = date_only + 0.5, ymin = -Inf, ymax = Inf),
    fill = "orange", alpha = 0.2, inherit.aes = FALSE
  ) +
  
  # DHW thresholds
  geom_hline(aes(yintercept = dhw4_scaled), color = "lightblue", linetype = "dashed", linewidth = 0.8) +
  geom_hline(aes(yintercept = dhw8_scaled), color = "blue", linetype = "dashed", linewidth = 0.8) +
  
  # Right-hand side labels
  geom_text(data = insitu_labels, aes(x = x_label, y = mmm_end - 0.6, label = "MMM"),
            color = "pink", hjust = 0, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels, aes(x = x_label, y = mmm1_end + 1.0, label = "MMM + 1°C"),
            color = "red", hjust = 0, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels, aes(x = x_label, y = dhw4_end + 0.8, label = "DHW 4"),
            color = "lightblue", hjust = 0, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels, aes(x = x_label, y = dhw8_end + 0.8, label = "DHW 8"),
            color = "blue", hjust = 0, size = 3, fontface = "bold") +
  
  # Facets
  facet_wrap(~ site, scales = "free_y") +
  
  # Axes
  scale_y_continuous(
    name = "SST (°C)",
    sec.axis = sec_axis(~ . / max(insitu_plot_df$scale_factor, na.rm = TRUE),
                        name = "DHW (°C-weeks)")
  ) +
  
  labs(
    x = "Date",
    title = "In Situ Daily Max SST and DHW",
    subtitle = "Orange shading = HotSpot days (≥ MMM + 1°C)"
  ) +
  
  theme_noaa()

# -----------------------------
# 4️⃣ Satellite SST + DHW plot
# -----------------------------
fig_satellite <- ggplot(satellite_plot_df, aes(x = date_only)) +
  
  geom_line(aes(y = max_temp), color = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm + 1), color = "red", linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = mmm), color = "pink", linetype = "dashed", linewidth = 0.8) +
  geom_line(aes(y = dhw_scaled), color = "blue", linewidth = 0.8) +
  
  geom_rect(
    data = satellite_plot_df %>% filter(hotspot >= 1),
    aes(xmin = date_only - 0.5, xmax = date_only + 0.5, ymin = -Inf, ymax = Inf),
    fill = "orange", alpha = 0.2, inherit.aes = FALSE
  ) +
  
  geom_hline(aes(yintercept = dhw4_scaled), color = "lightblue", linetype = "dashed", linewidth = 0.8) +
  geom_hline(aes(yintercept = dhw8_scaled), color = "blue", linetype = "dashed", linewidth = 0.8) +
  
  geom_text(data = satellite_labels, aes(x = x_label, y = mmm_end - 1.3, label = "MMM"),
            color = "pink", hjust = 0, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels, aes(x = x_label, y = mmm1_end + 0.7, label = "MMM + 1°C"),
            color = "red", hjust = 0, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels, aes(x = x_label, y = dhw4_end + 1.2, label = "DHW 4"),
            color = "lightblue", hjust = 0, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels, aes(x = x_label, y = dhw8_end + 1.4, label = "DHW 8"),
            color = "blue", hjust = 0, size = 3, fontface = "bold") +
  
  facet_wrap(~ site, scales = "free_y") +
  scale_y_continuous(
    name = "SST (°C)",
    sec.axis = sec_axis(~ . / max(satellite_plot_df$scale_factor, na.rm = TRUE),
                        name = "DHW (°C-weeks)", breaks = seq(0, 30, 5))
  ) +
  labs(
    x = "Date",
    title = "Satellite SST and DHW",
    subtitle = "Orange shading = HotSpot days (≥ MMM + 1°C)"
  ) +
  theme_noaa()

# -----------------------------
# 5️⃣ Combine In Situ & Satellite vertically
# -----------------------------
fig_combined <- fig_insitu / fig_satellite + plot_layout(ncol = 1, heights = c(1, 1))
fig_combined

# Optional save
# ggsave("fig_combined_sst_dhw.jpg", fig_combined, width = 12, height = 10, dpi = 300, bg = "white")

