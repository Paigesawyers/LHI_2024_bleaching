##===============================
## 2024 MARINE HEATWAVE ANALYSIS – Lord Howe Island
## SST / In Situ & Satellite Data
##===============================

# -----------------------------
# 0️⃣ Clean environment
# -----------------------------
rm(list = ls())

setwd("/users/paigesawyers/Desktop/paige_all_files")
# -----------------------------
# 1. Load required packages
# -----------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lubridate, janitor, zoo, data.table, 
               ggplot2, patchwork, viridis, scales, cowplot, paletteer, dplyr)


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
# 4. Load raw in situ logger data
# -----------------------------
insitu_raw <- read_csv("LHIMP_DATA_in_situ_loggers.csv") %>%
  clean_names() %>%
  filter(date_time != "") %>%
  mutate(date_time = dmy_hms(paste0(date_time, ":00")),
         date_only = as_date(date_time)) %>%
  mutate(across(where(is.numeric), ~na_if(., 0))) # replace 0 with NA

# -----------------------------
# 5. Convert in situ data to long format
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
# 6. Daily max, min & mean temperatures
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
    anomaly_max  = max_temp  - mmm,
    anomaly_mean = mean_temp - mmm,
    
    hotspot_max  = if_else(anomaly_max  >= 1, anomaly_max,  0),
    hotspot_mean = if_else(anomaly_mean >= 1, anomaly_mean, 0)
  )
# -----------------------------
# 7. Extract 2024 data
# -----------------------------
insitu_2024 <- insitu_daily %>% filter(year(date_only) == 2024)

# Max and min per site (separate tables)
max_temp_sites <- c("coral_gardens", "sylphs_hole", "north_bay",
                    "horseshoe_reef", "potholes_reef")

insitu_max_temp_2024 <- insitu_2024 %>%
  filter(site %in% max_temp_sites) %>%
  group_by(site) %>%
  slice_max(max_temp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(site,
         max_temp,
         date_of_max = date_only)

min_temp_sites <- c("coral_gardens", "sylphs_hole", "north_bay",
                    "horseshoe_reef", "potholes_reef")

insitu_min_temp_2024 <- insitu_2024 %>%
  filter(site %in% min_temp_sites) %>%
  group_by(site) %>%
  slice_min(min_temp, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(site,
         min_temp,
         date_of_min = date_only)


# -----------------------------
# 8. Compute Degree Heating Weeks (DHW)
# -----------------------------
insitu_dhw <- insitu_daily %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(
    dhw = zoo::rollapply(
      hotspot_max,
      width = dhw_window_days,
      FUN = function(x) sum(x) / 7,
      align = "right",
      partial = TRUE
    ),
    dhw = replace_na(dhw, 0)
  ) %>%
  ungroup()

# -----------------------------
# 9. Peak DHW and consecutive runs
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
# 9.a Mean DHW and consecutive runs (calculated for DHWs used in manuscript)

insitu_dhw_mean <- insitu_daily %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(
    dhw_mean = zoo::rollapply(
      hotspot_mean,
      width = dhw_window_days,
      FUN = function(x) sum(x) / 7,
      fill = NA,
      align = "right",
      partial = TRUE
    )
  ) %>%
  ungroup()


peak_insitu_dhw_mean <- insitu_dhw_mean %>%
  group_by(site) %>%
  filter(!is.na(dhw_mean)) %>%
  summarise(
    peak_dhw = max(dhw_mean, na.rm = TRUE),
    .groups = "drop"
  )

## rounding DHWS

insitu_dhw_mean_rounded <- insitu_dhw_mean %>%
  filter(!is.na(dhw_mean)) %>%
  mutate(dhw_mean_rounded = round(dhw_mean)) %>%
  group_by(site) %>%
  mutate(run_id = data.table::rleid(dhw_mean_rounded)) %>%
  group_by(site, dhw_mean_rounded, run_id) %>%
  summarise(
    start_date = min(date_only),
    end_date   = max(date_only),
    run_length = n(),
    .groups = "drop"
  )

# -----------------------------
# -----------------------------
# 10. Round MAX DHW for plotting consecutive runs
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
# 11. Satellite Data Preprocessing
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
# 12. Satellite DHW
# -----------------------------
satellite_dhw <- satellite %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(
    dhw = zoo::rollapply(
      hotspot,
      width = dhw_window_days,
      FUN = function(x) sum(x) / 7,
      fill = NA,
      align = "right",
      partial = TRUE
    )
  ) %>%
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
## PLOTTING: In Situ & Satellite SST + DHW SUP FIGURES
##===============================
# -----------------------------
# NOAA-style theme
# -----------------------------
theme_noaa <- function(base_size = 12) {
  theme_classic(base_size = base_size) +
    theme(
      strip.background = element_blank(),
      strip.text = element_text(face = "bold"),
      axis.title = element_text(face = "bold"),
      plot.title = element_text(face = "bold"),
      plot.subtitle = element_text(face = "bold", size = base_size - 2)
    )
}

# -----------------------------
# Scale DHW to fixed SST range (30°C) for plotting
# -----------------------------
scale_dhw <- function(df, dhw_col = "dhw", sst_max = 30) {
  df %>%
    group_by(site) %>%
    mutate(
      scale_factor = sst_max / 20,  # scale DHW to match max SST 30°C
      dhw_scaled   = .data[[dhw_col]] * scale_factor,
      dhw4_scaled  = 4 * scale_factor,
      dhw8_scaled  = 8 * scale_factor
    ) %>%
    ungroup()
}

# -----------------------------
# Prepare plotting data
# -----------------------------
insitu_plot_df <- scale_dhw(insitu_dhw, dhw_col = "dhw", sst_max = 30)
satellite_plot_df <- satellite_dhw %>%
  rename(mean_temp = sst) %>%
  scale_dhw(dhw_col = "dhw", sst_max = 30) %>%
  filter(date_only <= as.Date("2024-04-30"))

# Fix site names
insitu_plot_df <- insitu_plot_df %>%
  mutate(site = if_else(site == "coral_gardens", "stephens_hole", site))
satellite_plot_df <- satellite_plot_df %>%
  mutate(site = if_else(site == "coral_gardens", "stephens_hole", site))

# -----------------------------
# Scale in situ DHW to fixed SST range (30°C) for plotting only
# -----------------------------

max_insitu_dhw <- max(insitu_dhw$dhw, na.rm = TRUE)

insitu_plot_df <- insitu_dhw %>%
  mutate(
    dhw_scaled = dhw / max_insitu_dhw * 30,   # scale DHW to fill 0–30°C
    dhw4_scaled = 4 / max_insitu_dhw * 30,
    dhw8_scaled = 8 / max_insitu_dhw * 30
  )


# -----------------------------
# Labels for plotting thresholds
# -----------------------------
make_labels_fixed <- function(df) {
  df %>%
    group_by(site) %>%
    summarise(
      x_label = as.Date("2024-04-30") - 10,
      mmm_y   = max(mmm, na.rm = TRUE),
      mmm1_y  = max(mmm + 1, na.rm = TRUE),
      dhw4_y  = 4,  # fixed position on right-hand axis
      dhw8_y  = 8,
      .groups = "drop"
    )
}

insitu_labels_fixed <- make_labels_fixed(insitu_plot_df)
satellite_labels_fixed <- make_labels_fixed(satellite_plot_df)

# -----------------------------
# Adjust label positions for better separation
# -----------------------------
insitu_labels_fixed <- insitu_labels_fixed %>%
  mutate(
    mmm_y  = mmm_y  - 0.5,
    mmm1_y = mmm1_y + 0.5,
    dhw4_y = dhw4_y - 0.3,
    dhw8_y = dhw8_y + 0.3
  )

satellite_labels_fixed <- satellite_labels_fixed %>%
  mutate(
    mmm_y  = mmm_y  - 0.5,
    mmm1_y = mmm1_y + 0.5,
    dhw4_y = dhw4_y - 0.3,
    dhw8_y = dhw8_y + 0.3
  )

# -----------------------------
# Scale in situ DHW to fit y-axis 0-30°C
# -----------------------------
max_insitu_dhw <- max(insitu_dhw$dhw, na.rm = TRUE)

insitu_plot_df <- insitu_dhw %>%
  mutate(
    dhw_scaled  = dhw / max_insitu_dhw * 30,   # scale DHW to fill 0–30°C
    dhw4_scaled = 4 / max_insitu_dhw * 30,
    dhw8_scaled = 8 / max_insitu_dhw * 30
  )

# -----------------------------
# In Situ plot
# -----------------------------
fig_insitu <- ggplot(insitu_plot_df, aes(x = date_only)) +
  
  # SST lines
  geom_line(aes(y = max_temp), colour = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm), colour = "pink", linetype = "dashed") +
  geom_line(aes(y = mmm + 1), colour = "red", linetype = "dashed") +
  
  # HotSpot shading
  geom_rect(
    data = insitu_plot_df %>% filter(hotspot_max >= 1),
    aes(xmin = date_only - 0.5,
        xmax = date_only + 0.5,
        ymin = -Inf, ymax = Inf),
    fill = "orange", alpha = 0.2, inherit.aes = FALSE
  ) +
  
  # DHW line
  geom_line(aes(y = dhw_scaled), colour = "blue", linewidth = 1) +
  
  # DHW thresholds
  geom_hline(aes(yintercept = dhw4_scaled), colour = "lightblue", linetype = "dashed") +
  geom_hline(aes(yintercept = dhw8_scaled), colour = "blue", linetype = "dashed") +
  
  # Labels
  geom_text(data = insitu_labels_fixed,
            aes(x = x_label, y = mmm_y, label = "MMM"),
            colour = "pink", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels_fixed,
            aes(x = x_label, y = mmm1_y, label = "MMM + 1°C"),
            colour = "red", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels_fixed,
            aes(x = x_label, y = dhw4_y, label = "eDHW 4"),
            colour = "lightblue", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = insitu_labels_fixed,
            aes(x = x_label, y = dhw8_y, label = "eDHW 8"),
            colour = "blue", hjust = 1.1, size = 3, fontface = "bold") +
  
  facet_wrap(~ site, scales = "free_y") +
  
  scale_x_date(limits = c(as.Date("2024-01-01"), max(insitu_plot_df$date_only)),
               expand = expansion(mult = c(0.01, 0.02))) +
  
  scale_y_continuous(
    name = "IST (°C)",
    sec.axis = sec_axis(
      trans = ~ . * (max_insitu_dhw / 30),  # convert scaled back to DHW units
      name = "eDHW (°C-weeks)"
    )
  ) +
  
  coord_cartesian(ylim = c(NA, 30)) +
  
  # Updated title: only "In Situ" italic, rest normal
  labs(
    x = "Date",
    title = expression(italic("In Situ") ~ "daily maximum IST and experimental Degree Heating Weeks")
  ) +
  theme_noaa()

# -----------------------------
# Satellite plot
# -----------------------------
fig_satellite <- ggplot(satellite_plot_df, aes(x = date_only)) +
  
  # SST lines
  geom_line(aes(y = mean_temp), colour = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm), colour = "pink", linetype = "dashed") +
  geom_line(aes(y = mmm + 1), colour = "red", linetype = "dashed") +
  
  # HotSpot shading
  geom_rect(
    data = satellite_plot_df %>% filter(hotspot >= 1),
    aes(xmin = date_only - 0.5,
        xmax = date_only + 0.5,
        ymin = -Inf, ymax = Inf),
    fill = "orange", alpha = 0.2, inherit.aes = FALSE
  ) +
  
  # DHW line
  geom_line(aes(y = dhw_scaled), colour = "blue", linewidth = 1, na.rm = TRUE) +
  
  # DHW thresholds
  geom_hline(aes(yintercept = dhw4_scaled), colour = "lightblue", linetype = "dashed") +
  geom_hline(aes(yintercept = dhw8_scaled), colour = "blue", linetype = "dashed") +
  
  # Labels
  geom_text(data = satellite_labels_fixed,
            aes(x = x_label, y = mmm_y, label = "MMM"),
            colour = "pink", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels_fixed,
            aes(x = x_label, y = mmm1_y, label = "MMM + 1°C"),
            colour = "red", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels_fixed,
            aes(x = x_label, y = dhw4_y, label = "eDHW 4"),
            colour = "lightblue", hjust = 1.1, size = 3, fontface = "bold") +
  geom_text(data = satellite_labels_fixed,
            aes(x = x_label, y = dhw8_y, label = "eDHW 8"),
            colour = "blue", hjust = 1.1, size = 3, fontface = "bold") +
  
  facet_wrap(~ site, scales = "free_y") +
  
  scale_x_date(limits = c(as.Date("2024-01-01"), as.Date("2024-04-30")),
               breaks = as.Date(c("2024-01-01","2024-02-01","2024-03-01","2024-04-01")),
               labels = c("Jan","Feb","Mar","Apr"),
               expand = expansion(add = 0)) +
  
  scale_y_continuous(
    name = "SST (°C)",
    sec.axis = sec_axis(
      transform = ~ . * (20 / 30),
      name = "DHWs (°C-weeks)"
    )
  ) +
  
  coord_cartesian(ylim = c(NA, 30)) +
  
  # Updated title: normal, no bold
  labs(
    x = "Date",
    title = "Satellite SST and DHWs"
  ) +
  theme_noaa() +
  theme(plot.title = element_text(face = "plain"))  # overrides bold from theme

# -----------------------------
# Combine plots
# -----------------------------
fig_combined <- fig_insitu / fig_satellite + plot_layout(heights = c(1,1))

# Display
fig_combined

# -----------------------------
# Optional save
# -----------------------------
ggsave(
  filename = "fig_combined.jpg",
  plot     = fig_combined,
  width    = 12,
  height   = 18,
  dpi      = 300,
  bg       = "white"
)

# -----------------------------
# Satellite plot sites combined
# -----------------------------

# -----------------------------
# 1. Combine sites: mean SST, max DHW
# -----------------------------
combined_df <- satellite_dhw %>%
  group_by(date_only) %>%
  summarise(
    mean_temp   = mean(sst, na.rm = TRUE),  # mean SST
    mmm         = mean(mmm, na.rm = TRUE),  # mean MMM
    dhw_scaled  = max(dhw, na.rm = TRUE),   # max DHW
    hotspot     = max(hotspot, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    dhw4_scaled = 4,
    dhw8_scaled = 8
  )

# -----------------------------
# 2. Extend dates to May 31 and carry last observation forward
# -----------------------------
combined_df_full <- combined_df %>%
  complete(date_only = seq(as.Date("2024-01-01"), as.Date("2024-05-31"), by = "day")) %>%
  fill(mean_temp, mmm, dhw_scaled, dhw4_scaled, dhw8_scaled, hotspot, .direction = "down")

# -----------------------------
# 3. Plot combined SST & DHW
# -----------------------------
fig_combined_sites <- ggplot(combined_df_full, aes(x = date_only)) +
  
  # HotSpot shading first (lighter)
  geom_rect(
    data = combined_df_full %>% filter(hotspot >= 1),
    aes(xmin = date_only - 0.5,
        xmax = date_only + 0.5,
        ymin = -Inf, ymax = Inf),
    fill = "orange", alpha = 0.1, inherit.aes = FALSE
  ) +
  
  # SST line (grey)
  geom_line(aes(y = mean_temp), colour = "grey50", linewidth = 0.8) +
  
  # MMM lines (deeper pink)
  geom_line(aes(y = mmm), colour = "#D61C59", linetype = "dashed") +
  geom_line(aes(y = mmm + 1), colour = "red", linetype = "dashed") +
  
  # DHW line (black)
  geom_line(aes(y = dhw_scaled), colour = "black", linewidth = 1) +
  
  # DHW thresholds
  geom_hline(aes(yintercept = dhw4_scaled), colour = "#5DADE2", linetype = "dashed") +  # lighter blue
  geom_hline(aes(yintercept = dhw8_scaled), colour = "#2166AC", linetype = "dashed") +  # dark blue
  
  # Labels using annotate() with adjusted positions
  annotate("text", x = as.Date("2024-05-25"), y = max(combined_df$mmm) - 0.5,
           label = "MMM", colour = "#D61C59", hjust = 1.1, size = 3, fontface = "bold") +
  annotate("text", x = as.Date("2024-05-25"), y = max(combined_df$mmm)+1.5,
           label = "MMM + 1°C", colour = "red", hjust = 1.1, size = 3, fontface = "bold") +
  annotate("text", x = as.Date("2024-05-25"), y = 3.5,
           label = "DHW 4°C", colour = "#5DADE2", hjust = 1.1, size = 3, fontface = "bold") +
  annotate("text", x = as.Date("2024-05-25"), y = 7.5,
           label = "DHW 8°C", colour = "#2166AC", hjust = 1.1, size = 3, fontface = "bold") +
  
  # X-axis with Jan label nudged slightly to the right
  scale_x_date(
    limits = c(as.Date("2024-01-01"), as.Date("2024-05-31")),
    breaks = as.Date(c("2024-01-01","2024-02-01","2024-03-01","2024-04-01","2024-05-01")),
    labels = c("Jan","Feb","Mar","Apr","May"),
    expand = expansion(add = 0)
  ) +
  
  theme(axis.text.x = element_text(margin = margin(t = 5))) +  # small shift down, can tweak
  
  # Y-axis with secondary DHW axis
  scale_y_continuous(
    name = "SST (°C)",
    sec.axis = sec_axis(
      trans = ~ .,  # already in same scale
      name = "DHWs (°C-weeks)"
    )
  ) +
  
  coord_cartesian(ylim = c(NA, 30)) +
  
  # Title
  labs(
    x = "Date",
    title = "Satellite Sea Surface Temperatures and Degree Heating Weeks"
  ) +
  
  theme_noaa() +
  theme(plot.title = element_text(face = "plain"))

fig_combined_sites

# ggsave(
#   filename = "fig_combined_sites.jpg",
#   plot     = fig_combined_sites,
#   width    = 12,
#   height   = 8,
#   dpi      = 300,
#   bg       = "white"
# )

### adding seasonal bar
# -----------------------------
# Create seasonal gradient
# -----------------------------
season_df <- data.frame(
  date_only = seq(as.Date("2024-01-01"), as.Date("2024-05-31"), by = "day")
) %>%
  mutate(
    doy = yday(date_only),
    # Linear gradient: Jan = 1 (summer), May = 0.5 (autumn transition)
    season_value = 1 - (doy - 1) / (151 - 1) * 0.5
  )

# -----------------------------
# Add seasonal gradient bar below plot with strong endpoints
# -----------------------------
fig_combined_sites_season <- fig_combined_sites +
  geom_tile(
    data = season_df,
    aes(x = date_only, y = -2, fill = season_value),
    height = 1.5,   # slightly thicker for visibility
    inherit.aes = FALSE
  ) +
  scale_fill_gradientn(
    colours = c("#053061", "#2166AC", "#E34A33", "#67001F"),,
    values = scales::rescale(c(0, 0.5, 0.8, 1)),  # pushes colours to stronger endpoints
    limits = c(0,1),
    guide = "none"
  ) +
  coord_cartesian(ylim = c(-2, 30))  # expand y-axis for seasonal bar

fig_combined_sites_season

# ggsave(
#   filename = "fig_combined_sites_season.jpg",
#   plot     = fig_combined_sites_season,
#   width    = 12,
#   height   = 8,
#   dpi      = 300,
#   bg       = "white"
# )

##===============================
## SST / In Situ & Satellite Data RESULTS FIGURE
##===============================
# -----------------------------
# Define plotting date ranges
# -----------------------------
plot_start_sst <- insitu_start
plot_end_sst   <- insitu_end
plot_start_dhw <- insitu_start
plot_end_dhw   <- insitu_end

# -----------------------------
# In situ DHW for plotting
# -----------------------------
insitu_dhw_plot_data <- insitu_daily %>%
  filter(site %in% names(insitu_site_map),
         date_only >= plot_start_dhw,
         date_only <= plot_end_dhw) %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(
    dhw = zoo::rollapply(hotspot_mean,
                         width = dhw_window_days,
                         FUN = function(x) sum(x)/7,
                         align = "right",
                         partial = TRUE)
  ) %>%
  replace_na(list(dhw = 0)) %>%
  ungroup()

# Satellite DHW (max across sites)
satellite_dhw_max <- satellite_dhw %>%
  filter(date_only >= plot_start_dhw & date_only <= plot_end_dhw) %>%
  group_by(date_only) %>%
  summarise(dhw = max(dhw, na.rm = TRUE), .groups = "drop") %>%
  mutate(site = "Satellite")

# Combine DHW
dhw_combined <- bind_rows(
  insitu_dhw_plot_data %>% select(site, date_only, dhw),
  satellite_dhw_max
) %>%
  mutate(
    site = recode(site, !!!insitu_site_map, "Satellite" = "Satellite"),
    site = factor(site, levels = site_levels)
  )

# -----------------------------
# SST plot with MMM dashed line
# -----------------------------
fig_sst_filtered <- ggplot(sst_combined, aes(x = date_only, y = temp, colour = site)) +
  geom_line(aes(size = site), na.rm = TRUE) +
  # add black dashed MMM line
  geom_hline(yintercept = 24.06, colour = "black", linetype = "dashed", linewidth = 0.8) +
  scale_color_manual(values = site_colors, breaks = site_levels) +
  scale_size_manual(values = site_linewidths, guide = "none") +
  scale_x_date(
    breaks = seq(as.Date("2023-12-01"), as.Date("2024-05-31"), by = "1 month"),
    labels = label_date("%b %Y"),
    limits = c(as.Date("2023-12-01"), as.Date("2024-05-31"))  # cut off after May
  ) +
  labs(x = "Date", y = "Daily mean SST (°C)", colour = "Site / Source") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# -----------------------------
# DHW plot
# -----------------------------
fig_dhw_filtered <- ggplot(dhw_combined, aes(x = date_only, y = dhw, colour = site)) +
  geom_line(aes(size = site), na.rm = TRUE) +
  scale_color_manual(values = site_colors, breaks = site_levels) +
  scale_size_manual(values = site_linewidths, guide = "none") +
  scale_x_date(
    breaks = seq(as.Date("2023-12-01"), as.Date("2024-05-31"), by = "1 month"),
    labels = label_date("%b %Y"),
    limits = c(as.Date("2023-12-01"), as.Date("2024-05-31"))  # cut off after May
  ) +
  scale_y_continuous(limits = c(0, 25)) +
  labs(x = "Date", y = "Degree Heating Weeks (°C-weeks)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = "none")


# -----------------------------
# Shared legend
# -----------------------------

shared_legend <- get_legend(
  ggplot(sst_combined, aes(x = date_only, y = temp, colour = site)) +
    geom_line(size = 1, na.rm = TRUE) +
    scale_color_manual(values = site_colors, breaks = site_levels) +
    labs(colour = NULL) +
    guides(colour = guide_legend(override.aes = list(size = 1.2))) +
    theme(legend.position = "top",
          legend.text = element_text(size = 12),
          legend.key = element_blank())
)

# -----------------------------
# Combine plots side-by-side with shared legend
# -----------------------------
fig_sst_combined <- plot_grid(
  shared_legend,
  plot_grid(fig_sst_filtered, fig_dhw_filtered, ncol = 2, align = "hv"),
  ncol = 1,
  rel_heights = c(0.1, 1)
)

# -----------------------------
# Show figure
# -----------------------------
fig_sst_combined

# -----------------------------
# Optional save
# -----------------------------
# ggsave(
#   filename = "fig_sst_combined.jpg",
#   plot     = fig_sst_combined,
#   width    = 12,
#   height   = 8,
#   dpi      = 300,
#   bg       = "white"
# )

# -----------------------------
# PAST BLEACHING EVENTS RESULTS FIGURE
# -----------------------------

# Load and process SST data ----
df_sst <- read_csv("LHI_sst_locations_timeseries.csv") %>%
  mutate(date_only = dmy(date)) %>%
  pivot_longer(cols = where(is.numeric), names_to = "site", values_to = "sst") %>%
  left_join(read_csv("LHI_mmms.csv")) %>%
  filter(!is.na(mmm)) %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(
    anomaly = sst - mmm,
    hotspot = if_else(anomaly >= 1, anomaly, 0),
    sdhw = rollapply(hotspot, width = 84, FUN = function(x) sum(x)/7, align = "right", fill = NA)
  ) %>%
  select(date_only, site, sst, sdhw) %>%
  ungroup()

# Filter for years and sites of interest ----
years_of_interest <- c(1998, 2010, 2019, 2024)
df_sst_filtered <- df_sst %>%
  filter(year(date_only) %in% years_of_interest, site != "le_meurthe") %>%
  mutate(year = year(date_only), day_of_year = yday(date_only)) %>%
  filter(month(date_only) <= 7) # January–July

month_starts <- c(1, 32, 60, 91, 121, 152, 182)
month_labels <- month.abb[1:7]
MMM <- 24.06

# SST plot: overlapping years with MMM ----
plot_sst_year <- df_sst_filtered %>%
  group_by(year, day_of_year) %>%
  summarise(sst = mean(sst, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_year, y = sst, color = factor(year))) +
  geom_line(size = 1) +
  geom_hline(yintercept = MMM, linetype = "dotted", color = "black", size = 0.8) +
  annotate("text", x = 5, y = MMM + 0.15, label = "MMM", color = "black", hjust = 0) +
  labs(title = "Satellite SST", x = "Month", y = "SST (°C)", color = "Year") +
  scale_x_continuous(breaks = month_starts, labels = month_labels) +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

# DHW plots: mean and max across sites ----
plot_dhw_mean <- df_sst_filtered %>%
  group_by(year, day_of_year) %>%
  summarise(dhw = mean(sdhw, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_year, y = dhw, color = factor(year))) +
  geom_line(size = 1) +
  labs(title = "Satellite DHW", x = "Month", y = "DHW (°C-weeks)", color = "Year") +
  scale_x_continuous(breaks = month_starts, labels = month_labels) +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

plot_dhw_max <- df_sst_filtered %>%
  group_by(year, day_of_year) %>%
  summarise(dhw = max(sdhw, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_year, y = dhw, color = factor(year))) +
  geom_line(size = 1) +
  labs(title = "Satellite DHW", x = "Month", y = "DHW (°C-weeks)", color = "Year") +
  scale_x_continuous(breaks = month_starts, labels = month_labels) +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))

# Combine SST + DHW plots side by side ----
combined_mean <- plot_sst_year + plot_dhw_mean +
  plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "top")

combined_max <- plot_sst_year + plot_dhw_max +
  plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "top")

# Display
combined_mean
combined_max

# Max DHW summary for reference ----
max_sdhw <- df_sst_filtered %>%
  group_by(year) %>%
  slice_max(sdhw, n = 1, with_ties = FALSE) %>%
  select(year, site_sdhw = site, max_sdhw = sdhw)

max_dhw_summary <- max_sdhw
max_dhw_summary

## using pastels

# Years and colors
years_of_interest <- c(1998, 2010, 2019, 2024)
year_colors <- c(
  "1998" = "#A6CEE3",  # pastel blue
  "2010" = "#B2DF8A",  # pastel green
  "2019" = "#FB9A99",  # pastel red
  "2024" = "#FDBF6F"   # pastel orange
)

MMM <- 24.06

# Month info
month_starts_sst <- c(1,32,60,91,121,152,182,213,244,274,305,335)
month_labels_sst <- month.abb

month_starts_dhw <- c(335,1,32,60,91,121,152,182)
month_labels_dhw <- c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul")

# SST plot: January → December
plot_sst_year <- df_sst %>%
  filter(year(date_only) %in% years_of_interest, site != "le_meurthe") %>%
  mutate(year = year(date_only), day_of_year = yday(date_only)) %>%
  group_by(year, day_of_year) %>%
  summarise(sst = mean(sst, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_year, y = sst, color = factor(year))) +
  geom_line(size = 1.2) +
  geom_hline(yintercept = MMM, linetype = "dotted", color = "black", size = 0.8) +
  annotate("text", x = 5, y = MMM + 0.15, label = "MMM", color = "black", hjust = 0) +
  scale_x_continuous(breaks = month_starts_sst, labels = month_labels_sst) +
  scale_color_manual(values = year_colors) +
  labs(title = "Satellite SST", x = "Month", y = "SST (°C)", color = "Year") +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))


# DHW plot: December → July with Dec first
plot_dhw_max <- df_sst %>%
  filter(year(date_only) %in% years_of_interest, site != "le_meurthe") %>%
  mutate(year = year(date_only),
         # Create a new "month_order" for plotting Dec first
         month_order = ifelse(month(date_only) == 12, 0, month(date_only)),
         day_of_year_adj = yday(date_only) + ifelse(month(date_only) == 12, -365, 0)) %>%
  filter(month(date_only) == 12 | month(date_only) <= 7) %>%
  group_by(year, day_of_year_adj) %>%
  summarise(dhw = max(sdhw, na.rm = TRUE), .groups = "drop") %>%
  ggplot(aes(x = day_of_year_adj, y = dhw, color = factor(year))) +
  geom_line(size = 1.2) +
  scale_x_continuous(
    breaks = c(-30,1,32,60,91,121,152,182), # Dec ~-30, Jan-Jul
    labels = c("Dec","Jan","Feb","Mar","Apr","May","Jun","Jul")
  ) +
  scale_color_manual(values = year_colors) +
  labs(title = "Satellite DHWs", x = "Month", y = "DHW (°C-weeks)", color = "Year") +
  theme_classic() +
  theme(legend.position = "top", plot.title = element_text(hjust = 0.5))


# Combine
fig_updated_combined_max <- plot_sst_year + plot_dhw_max +
  plot_layout(ncol = 2, guides = "collect") & theme(legend.position = "top")

fig_updated_combined_max

