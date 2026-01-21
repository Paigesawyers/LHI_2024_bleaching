##===============================
## 2024 MARINE HEATWAVE ANALYSIS – Lord Howe Island
## SST / In Situ & Satellite Data
##===============================

# -----------------------------
# 0️⃣ Clean environment
# -----------------------------
rm(list = ls())

# -----------------------------
# 1. Load required packages
# -----------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, lubridate, janitor, zoo, data.table, 
               ggplot2, patchwork, viridis, scales, cowplot, paletteer, dplyr)
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
## PLOTTING: In Situ & Satellite SST + DHW
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
      plot.subtitle = element_text(size = base_size - 2)
    )
}

# -----------------------------
# Scale DHW to SST range
# -----------------------------
scale_dhw <- function(df, dhw_col = "dhw", sst_col, dhw_max_plot = 30) {
  df %>%
    group_by(site) %>%
    mutate(
      scale_factor = max(.data[[sst_col]], na.rm = TRUE) / dhw_max_plot,
      dhw_scaled   = .data[[dhw_col]] * scale_factor,
      dhw4_scaled  = 4 * scale_factor,
      dhw8_scaled  = 8 * scale_factor
    ) %>%
    ungroup()
}

# -----------------------------
# Prepare plotting data
# -----------------------------
insitu_plot_df <- scale_dhw(
  insitu_dhw,
  dhw_col = "dhw",
  sst_col = "max_temp"
)

satellite_plot_df <- satellite_dhw %>%
  rename(mean_temp = sst) %>%
  scale_dhw(
    dhw_col = "dhw",
    sst_col = "mean_temp"
  ) %>%
  filter(date_only <= as.Date("2024-04-30"))

## changing site name 

insitu_plot_df <- insitu_plot_df %>%
  mutate(site = if_else(site == "coral_gardens", "stephens_hole", site))

satellite_plot_df <- satellite_plot_df %>%
  mutate(site = if_else(site == "coral_gardens", "stephens_hole", site))

# -----------------------------
# Right-hand labels
# -----------------------------
make_labels <- function(df) {
  df %>%
    group_by(site) %>%
    summarise(
      x_label = as.Date("2024-04-30") - 10,
      
      mmm_y  = max(mmm, na.rm = TRUE),
      mmm1_y = max(mmm + 1, na.rm = TRUE),
      
      dhw4_y = unique(dhw4_scaled),
      dhw8_y = unique(dhw8_scaled),
      
      .groups = "drop"
    )
}

insitu_labels    <- make_labels(insitu_plot_df)
satellite_labels <- make_labels(satellite_plot_df)

# -----------------------------
# In situ SST + DHW plot
# -----------------------------
fig_insitu <- ggplot(insitu_plot_df, aes(x = date_only)) +
  
  # SST
  geom_line(aes(y = max_temp), colour = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm + 1), colour = "red",  linetype = "dashed") +
  geom_line(aes(y = mmm),     colour = "pink", linetype = "dashed") +
  
  # HotSpot shading (FIXED column)
  geom_rect(
    data = insitu_plot_df %>% filter(hotspot_max >= 1),
    aes(
      xmin = date_only - 0.5,
      xmax = date_only + 0.5,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "orange",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  # DHW (drawn LAST so it’s visible from Jan)
  geom_line(aes(y = dhw_scaled), colour = "blue", linewidth = 1) +
  
  # DHW thresholds
  geom_hline(aes(yintercept = dhw4_scaled),
             colour = "lightblue", linetype = "dashed") +
  geom_hline(aes(yintercept = dhw8_scaled),
             colour = "blue", linetype = "dashed") +
  # Labels
  geom_text(
    data = insitu_labels,
    aes(x = x_label, y = mmm_y - 0.5, label = "MMM"),
    colour = "pink", hjust = 0, size = 3, fontface = "bold"
  ) +
  geom_text(
    data = insitu_labels,
    aes(x = x_label, y = mmm1_y + 0.8, label = "MMM + 1°C"),
    colour = "red", hjust = 0, size = 3, fontface = "bold"
  ) +
  geom_text(
    data = insitu_labels,
    aes(x = x_label, y = dhw4_y - 0.5, label = "DHW 4"),
    colour = "lightblue", hjust = 0, size = 3, fontface = "bold"
  ) +
  geom_text(
    data = insitu_labels,
    aes(x = x_label, y = dhw8_y + 0.8, label = "DHW 8"),
    colour = "blue", hjust = 0, size = 3, fontface = "bold"
  ) +
  
  facet_wrap(~ site, scales = "free_y") +
  
  scale_x_date(
    limits = c(as.Date("2024-01-01"),
               max(insitu_plot_df$date_only)),
    expand = expansion(mult = c(0.01, 0.02))
  ) +
  
  scale_y_continuous(
    name = "SST (°C)",
    sec.axis = sec_axis(
      ~ . / max(insitu_plot_df$scale_factor, na.rm = TRUE),
      name = "DHW (°C-weeks)"
    )
  ) +
  
  labs(
    x = "Date",
    title = "In Situ Daily Maximum SST and DHW",
    subtitle = "DHW calculated from daily maximum temperature"
  ) +
  theme_noaa()

# -----------------------------
# Satellite SST + DHW plot for individual sites
# -----------------------------
fig_satellite <- ggplot(satellite_plot_df, aes(x = date_only)) +
  
  # SST lines
  geom_line(aes(y = mean_temp), colour = "black", linewidth = 0.8) +
  geom_line(aes(y = mmm + 1), colour = "red",  linetype = "dashed", linewidth = 0.7) +
  geom_line(aes(y = mmm),     colour = "pink", linetype = "dashed", linewidth = 0.7) +
  
  # HotSpot shading
  geom_rect(
    data = satellite_plot_df %>% filter(hotspot >= 1),
    aes(
      xmin = date_only - 0.5,
      xmax = date_only + 0.5,
      ymin = -Inf,
      ymax = Inf
    ),
    fill = "orange",
    alpha = 0.2,
    inherit.aes = FALSE
  ) +
  
  # DHW
  geom_line(aes(y = dhw_scaled), colour = "blue", linewidth = 1, na.rm = TRUE) +
  
  # DHW thresholds
  geom_hline(aes(yintercept = dhw4_scaled),
             colour = "lightblue", linetype = "dashed") +
  geom_hline(aes(yintercept = dhw8_scaled),
             colour = "blue", linetype = "dashed") +
  
  # 🔤 Labels
  geom_text(data = satellite_labels,
            aes(x = x_label, y = mmm_y - 0.5, label = "MMM"),
            colour = "pink", hjust = 0, size = 3, fontface = "bold") +
  
  geom_text(data = satellite_labels,
            aes(x = x_label, y = mmm1_y + 0.8, label = "MMM + 1°C"),
            colour = "red", hjust = 0, size = 3, fontface = "bold") +
  
  geom_text(data = satellite_labels,
            aes(x = x_label, y = dhw4_y - 0.5, label = "DHW 4"),
            colour = "lightblue", hjust = 0, size = 3, fontface = "bold") +
  
  geom_text(data = satellite_labels,
            aes(x = x_label, y = dhw8_y + 0.8, label = "DHW 8"),
            colour = "blue", hjust = 0, size = 3, fontface = "bold") +
  
  facet_wrap(~ site, scales = "free_y") +
  
  scale_x_date(
    limits = c(as.Date("2024-01-01"), as.Date("2024-04-30")),
    breaks = as.Date(c(
      "2024-01-01",
      "2024-02-01",
      "2024-03-01",
      "2024-04-01"
    )),
    labels = c("Jan", "Feb", "Mar", "Apr"),
    expand = expansion(add = 0)
  ) +
  
  scale_y_continuous(
    name = "SST (°C)",
    sec.axis = sec_axis(
      ~ . / max(satellite_plot_df$scale_factor, na.rm = TRUE),
      name = "DHW (°C-weeks)"
    )
  ) +
  
  labs(
    x = "Date",
    title = "Satellite SST and DHW",
    subtitle = "Orange shading = HotSpot days (≥ MMM + 1°C)"
  ) +
  
  theme_noaa()


# -----------------------------
# Combine
# -----------------------------
fig_combined <- fig_insitu / fig_satellite +
  plot_layout(heights = c(1, 1))

fig_combined

# # Optional save
# ggsave(
#   filename = "fig_combined.jpg",
#   plot     = fig_combined,
#   width    = 12,
#   height   = 18,
#   dpi      = 300,
#   bg       = "white"
# )

# =============================
# SST & DHW COMBINED FIGURE Uses DAILY MEAN in situ SST
# =============================

# -----------------------------
# In situ SST (daily MEAN)
# -----------------------------
insitu_sst_plot <- insitu_daily %>%
  filter(
    site %in% c("coral_gardens","sylphs_hole","north_bay",
                "horseshoe_reef","potholes_reef"),
    date_only >= plot_start_sst,
    date_only <= plot_end_sst
  ) %>%
  select(site, date_only, temp = mean_temp)

# -----------------------------
# Satellite SST (daily mean across sites)
# -----------------------------
satellite_mean_plot <- satellite_full %>%
  filter(date_only >= plot_start_sst & date_only <= plot_end_sst) %>%
  select(site, date_only, temp = sst) %>%
  group_by(date_only) %>%
  summarise(temp = mean(temp, na.rm = TRUE), .groups = "drop") %>%
  mutate(site = "Satellite")

# -----------------------------
# Combine SST
# -----------------------------
sst_combined <- bind_rows(insitu_sst_plot, satellite_mean_plot) %>%
  mutate(
    site = factor(
      str_to_sentence(str_replace_all(site, "_", " ")),
      levels = site_levels
    )
  )

# -----------------------------
# SST plot
# -----------------------------
fig_sst_filtered <- ggplot(
  sst_combined,
  aes(x = date_only, y = temp, colour = site)
) +
  geom_line(aes(size = site)) +
  scale_color_manual(values = site_colors, breaks = site_levels) +
  scale_size_manual(values = site_linewidths, guide = "none") +
  geom_hline(
    yintercept = MMM,
    linetype = "dotted",
    colour = "black",
    linewidth = 0.8
  ) +
  scale_x_date(
    breaks = date_breaks("1 month"),
    labels = label_date("%b %Y"),
    limits = c(plot_start_sst, plot_end_sst)
  ) +
  labs(
    x = "Date",
    y = "Daily mean SST (°C)",
    colour = "Site / Source"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# -----------------------------
# In situ DHW
# -----------------------------
insitu_dhw_plot_data <- insitu_dhw %>%
  filter(
    site %in% c("coral_gardens","sylphs_hole","north_bay",
                "horseshoe_reef","potholes_reef"),
    date_only >= plot_start_dhw,
    date_only <= plot_end_dhw
  ) %>%
  arrange(site, date_only) %>%
  group_by(site) %>%
  mutate(
    dhw = sapply(seq_along(hotspot), function(i) {
      window_start <- max(1, i - dhw_window_days + 1)
      sum(hotspot[window_start:i], na.rm = TRUE) / 7
    })
  ) %>%
  ungroup()

# -----------------------------
# Satellite DHW (daily max across sites)
# -----------------------------
satellite_dhw_max <- satellite_dhw %>%
  filter(date_only >= plot_start_dhw & date_only <= plot_end_dhw) %>%
  group_by(date_only) %>%
  summarise(dhw = max(dhw, na.rm = TRUE), .groups = "drop") %>%
  mutate(site = "Satellite")

# -----------------------------
# Combine DHW
# -----------------------------
dhw_combined <- bind_rows(
  insitu_dhw_plot_data %>% select(site, date_only, dhw),
  satellite_dhw_max
) %>%
  mutate(
    site = factor(
      str_to_sentence(str_replace_all(site, "_", " ")),
      levels = site_levels
    )
  )

# -----------------------------
# DHW plot
# -----------------------------
fig_dhw_filtered <- ggplot(
  dhw_combined,
  aes(x = date_only, y = dhw, colour = site)
) +
  geom_line(aes(size = site)) +
  scale_color_manual(values = site_colors, breaks = site_levels) +
  scale_size_manual(values = site_linewidths, guide = "none") +
  scale_x_date(
    breaks = date_breaks("1 month"),
    labels = label_date("%b %Y"),
    limits = c(plot_start_dhw, plot_end_dhw)
  ) +
  scale_y_continuous(limits = c(0, 30)) +
  labs(
    x = "Date",
    y = "Degree Heating Weeks (°C-weeks)"
  ) +
  theme_classic() +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

# -----------------------------
# Shared legend (final site names)
# -----------------------------
legend_labels <- c(
  "Satellite",
  "North Bay",
  "Sylphs Hole",
  "Stephens Hole",
  "Horseshoe Reef",
  "Potholes"
)

shared_legend <- get_legend(
  ggplot(sst_combined, aes(x = date_only, y = temp, colour = site)) +
    geom_line(size = 1) +
    scale_color_manual(
      values = site_colors,
      breaks = site_levels,
      labels = legend_labels
    ) +
    labs(colour = NULL) +
    guides(colour = guide_legend(
      override.aes = list(size = c(1.5, rep(1, 5)))
    )) +
    theme(
      legend.position = "top",
      legend.text = element_text(size = 12),
      legend.key = element_blank()
    )
)

# -----------------------------
# Final combined figure
# -----------------------------
fig_sst_combined <- plot_grid(
  shared_legend,
  plot_grid(fig_sst_filtered, fig_dhw_filtered, ncol = 2, align = "hv"),
  ncol = 1,
  rel_heights = c(0.1, 1)
)

fig_sst_combined

# # Optional save
# ggsave(
#   filename = "fig_sst_combined.jpg",
#   plot     = fig_sst_combined,
#   width    = 12,
#   height   = 18,
#   dpi      = 300,
#   bg       = "white"
# )

