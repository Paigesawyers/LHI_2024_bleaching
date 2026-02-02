# -------------------------------------------------------------
# 1. Load data + required packages
# -------------------------------------------------------------

# Load cleaned data objects
source("LHI_coral_load_data.R")
setwd("/users/paigesawyers/Desktop/Re_Updated_LHI_data_and_labelset")

# Load packages (using pacman for convenience)
if (!require("pacman")) {
  install.packages("pacman")
}

pacman::p_load(
  janitor,
  patchwork,
  scales,
  tidyverse,
  viridis,
  zoo,
  readr,
  stringr
)

# -------------------------------------------------------------
# 2. Define plotting colours
# -------------------------------------------------------------
chosen_colours <- rev(c(
  "#FA5454",
  "#EC9B31",
  "#F3ED37",
  "#5EDA29",
  "#29DABA",
  "#4378E9",
  "#B54BFC"
))

# -------------------------------------------------------------
# 3. Load all model .rds files 
# -------------------------------------------------------------
rds_files <- list.files(
  path = "./models/",
  pattern = ".rds",
  recursive = TRUE
)

for (x in rds_files) {
  assign(
    str_remove(x, ".rds"),
    read_rds(paste0("./models/", x))
  )
}

gc()   # clean up memory


# -------------------------------------------------------------
# 4. All coral cover vs non-coral cover (Binomial model)
# -------------------------------------------------------------
(fig_all_coral_cover <-
   emmeans(model_all_coral_cover,
           ~ Month * Site,
           type = "response") %>%
   as.data.frame() %>%
   mutate(
     Site = factor(
       case_when(
         Site == "AG"    ~ "Acropora\nGardens",
         Site == "CH"    ~ "Comets\nHole",
         Site == "HR"    ~ "Horseshoe\nReef",
         Site == "NB"    ~ "Neds\nBeach",
         Site == "North" ~ "North\nBay",
         Site == "PH"    ~ "Potholes",
         Site == "SH"    ~ "Stephens\nHole",
         Site == "Sylphs"~ "Sylphs\nHole"
       ),
       levels = c(
         "Neds\nBeach",
         "Sylphs\nHole",
         "North\nBay",
         "Acropora\nGardens",
         "Stephens\nHole",
         "Comets\nHole",
         "Horseshoe\nReef",
         "Potholes"
       )
     )
   ) %>%
   ggplot(aes(
     x = Month,
     y = prob,
     ymax = asymp.LCL,
     ymin = asymp.UCL
   )) +
   geom_col(fill = "grey") +
   geom_errorbar(width = 0.2) +
   facet_grid(. ~ Site) +
   scale_y_continuous(
     expand = c(0, 0.01),
     limits = c(0, 1),
     breaks = seq(0, 1, 0.2)
   ) +
   labs(
     x = "Month",
     y = "Proportion of coral cover"
   ) +
   theme_classic() +
   theme(
     legend.position = "top",
     legend.text = element_text(size = 12),
     legend.title = element_text(size = 14),
     legend.box = "horizontal",
     legend.direction = "horizontal",
     axis.text.x = element_text(
       size = 12,
       angle = 90,
       vjust = 0.5,
       hjust = 1
     ),
     axis.title.x = element_text(
       margin = margin(t = 0.3, unit = "cm"),
       size = 14
     ),
     axis.title.y = element_text(
       margin = margin(r = 0.3, unit = "cm"),
       size = 14
     ),
     strip.background = element_blank(),
     strip.text = element_text(size = 12)
   )
)


# -------------------------------------------------------------
# 5. Live coral vs non-coral or dead coral (Binomial model)
# -------------------------------------------------------------
(fig_live_coral_cover <-
   emmeans(model_live_coral_cover,
           ~ Month * Site,
           type = "response") %>%
   as.data.frame() %>%
   mutate(
     Site = factor(
       case_when(
         Site == "AG"    ~ "Acropora\nGardens",
         Site == "CH"    ~ "Comets\nHole",
         Site == "HR"    ~ "Horseshoe\nReef",
         Site == "NB"    ~ "Neds\nBeach",
         Site == "North" ~ "North\nBay",
         Site == "PH"    ~ "Potholes",
         Site == "SH"    ~ "Stephens\nHole",
         Site == "Sylphs"~ "Sylphs\nHole"
       ),
       levels = c(
         "Neds\nBeach",
         "Sylphs\nHole",
         "North\nBay",
         "Acropora\nGardens",
         "Stephens\nHole",
         "Comets\nHole",
         "Horseshoe\nReef",
         "Potholes"
       )
     )
   ) %>%
   ggplot(aes(
     x = Month,
     y = prob,
     ymax = asymp.LCL,
     ymin = asymp.UCL
   )) +
   geom_col(fill = "grey") +
   geom_errorbar(width = 0.2) +
   facet_grid(. ~ Site) +
   scale_y_continuous(
     expand = c(0, 0.01),
     limits = c(0, 1),
     breaks = seq(0, 1, 0.2)
   ) +
   labs(
     x = "Month",
     y = "Proportion live coral cover"
   ) +
   theme_classic() +
   theme(
     legend.position = "top",
     legend.text = element_text(size = 12),
     legend.title = element_text(size = 14),
     legend.box = "horizontal",
     legend.direction = "horizontal",
     axis.text.x = element_text(
       size = 12,
       angle = 90,
       vjust = 0.5,
       hjust = 1
     ),
     axis.title.x = element_text(
       margin = margin(t = 0.3, unit = "cm"),
       size = 14
     ),
     axis.title.y = element_text(
       margin = margin(r = 0.3, unit = "cm"),
       size = 14
     ),
     strip.background = element_blank(),
     strip.text = element_text(size = 12)
   )
)


# -------------------------------------------------------------
# 6. Reordered version (matching NMDS site order)
# -------------------------------------------------------------
(fig_live_coral_cover <-
   emmeans(model_live_coral_cover,
           ~ Month * Site,
           type = "response") %>%
   as.data.frame() %>%
   mutate(
     Site = factor(
       case_when(
         Site == "AG"    ~ "Acropora\nGardens",
         Site == "CH"    ~ "Comets\nHole",
         Site == "HR"    ~ "Horseshoe\nReef",
         Site == "NB"    ~ "Neds\nBeach",
         Site == "North" ~ "North\nBay",
         Site == "PH"    ~ "Potholes",
         Site == "SH"    ~ "Stephens\nHole",
         Site == "Sylphs"~ "Sylphs\nHole"
       ),
       levels = c(
         "Potholes",
         "Stephens\nHole",
         "Neds\nBeach",
         "Sylphs\nHole",
         "Horseshoe\nReef",
         "Comets\nHole",
         "North\nBay",
         "Acropora\nGardens"
       )
     )
   ) %>%
   ggplot(aes(
     x = Month,
     y = prob,
     ymax = asymp.LCL,
     ymin = asymp.UCL
   )) +
   geom_col(fill = "grey") +
   geom_errorbar(width = 0.2) +
   facet_grid(. ~ Site) +
   scale_y_continuous(
     expand = c(0, 0.01),
     limits = c(0, 1),
     breaks = seq(0, 1, 0.2)
   ) +
   labs(
     x = "Month",
     y = "Proportion live coral cover"
   ) +
   theme_classic() +
   theme(
     legend.position = "top",
     legend.text = element_text(size = 12),
     legend.title = element_text(size = 14),
     legend.box = "horizontal",
     legend.direction = "horizontal",
     axis.text.x = element_text(
       size = 12,
       angle = 90,
       vjust = 0.5,
       hjust = 1
     ),
     axis.title.x = element_text(
       margin = margin(t = 0.3, unit = "cm"),
       size = 14
     ),
     axis.title.y = element_text(
       margin = margin(r = 0.3, unit = "cm"),
       size = 14
     ),
     strip.background = element_blank(),
     strip.text = element_text(size = 12)
   )
)

## genera distribution ----
### proportion ----
#### betabinomial ----

# Identify model objects and clean names
model_genera_objects <-
  ls()[str_detect(ls(), "model_genera_betabinom_")]

model_genera_names <- 
  str_remove_all(model_genera_objects, "model_genera_betabinom_")

# Collect betabinomial predictions for genera
for (i in seq_along(model_genera_objects)) {
  
  temp <-
    get(model_genera_objects[i]) %>% 
    emmeans::emmeans(~ Month * Site, type = "response") %>% 
    as.data.frame() %>% 
    mutate(Genus_simplified = model_genera_names[i])
  
  if (i == 1) {
    genera_betabinom_results <- temp
  } else {
    genera_betabinom_results <- bind_rows(genera_betabinom_results, temp)
  }
}

# Plot genera distribution
fig_genera_distribution <-
  genera_betabinom_results %>% 
  group_by(Month, Site) %>% 
  mutate(total = sum(prob)) %>% 
  ungroup() %>% 
  mutate(
    prob_normalised = prob / total,
    Site = factor(case_when(
      Site == "AG" ~ "Acropora\nGardens",
      Site == "CH" ~ "Comets\nHole",
      Site == "HR" ~ "Horseshoe\nReef",
      Site == "NB" ~ "Neds\nBeach",
      Site == "North" ~ "North\nBay",
      Site == "PH" ~ "Potholes",
      Site == "SH" ~ "Stephens\nHole",
      Site == "Sylphs" ~ "Sylphs\nHole"
    ),
    levels = c(
      "Neds\nBeach", "Sylphs\nHole", "North\nBay",
      "Acropora\nGardens", "Stephens\nHole",
      "Comets\nHole", "Horseshoe\nReef", "Potholes"
    )),
    Genus_simplified = factor(
      ifelse(Genus_simplified == "other_genera", "Others", Genus_simplified),
      levels = rev(levels(factor(ifelse(Genus_simplified == "other_genera", "Others", Genus_simplified))))
    )
  ) %>% 
  ggplot(aes(x = Month, y = prob_normalised, fill = Genus_simplified)) +
  geom_col() +
  facet_grid(. ~ Site) +
  scale_y_continuous(expand = c(0, 0.01), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = chosen_colours) +
  labs(x = "Month",
       y = "Proportion from live coral cover",
       fill = "Genus") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.box = "horizontal",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 12, angle = 90,
                               vjust = 0.5, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 0.3, unit = "cm"),
                                size = 14),
    axis.title.y = element_text(margin = margin(r = 0.3, unit = "cm"),
                                size = 14),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  guides(fill = guide_legend(nrow = 1, reverse = TRUE))


## health (ignoring genera) ----
### proportion ----
#### betabinomial ----

# Identify health model objects and names
model_health_all_objects <-
  ls()[str_detect(ls(), "model_health_all_betabinom_")]

model_health_all_names <- 
  str_remove_all(model_health_all_objects, "model_health_all_betabinom_")

# Collect predictions
for (i in seq_along(model_health_all_objects)) {
  
  temp <-
    get(model_health_all_objects[i]) %>% 
    emmeans::emmeans(~ Month * Site, type = "response") %>% 
    as.data.frame() %>% 
    mutate(Health_simplified = model_health_all_names[i])
  
  if (i == 1) {
    health_all_betabinom_results <- temp
  } else {
    health_all_betabinom_results <- bind_rows(health_all_betabinom_results, temp)
  }
}

# Plot health–all
fig_health_all <-
  health_all_betabinom_results %>% 
  tibble() %>% 
  group_by(Site, Month) %>% 
  mutate(total = sum(prob)) %>% 
  ungroup() %>% 
  mutate(
    estimate_normalised = prob / total,
    group = factor(
      Health_simplified,
      levels = rev(c("Healthy", "Pale", "Bleached", "Dead", "LTPM"))
    ),
    Site = factor(case_when(
      Site == "AG" ~ "Acropora\nGardens",
      Site == "CH" ~ "Comets\nHole",
      Site == "HR" ~ "Horseshoe\nReef",
      Site == "NB" ~ "Neds\nBeach",
      Site == "North" ~ "North\nBay",
      Site == "PH" ~ "Potholes",
      Site == "SH" ~ "Stephens\nHole",
      Site == "Sylphs" ~ "Sylphs\nHole"
    ),
    levels = c(
      "Neds\nBeach", "Stephens\nHole", "Potholes",
      "Horseshoe\nReef", "Sylphs\nHole", "Comets\nHole",
      "North\nBay", "Acropora\nGardens"
    ))
  ) %>% 
  ggplot(aes(
    x = Month,
    y = estimate_normalised,
    pattern = group,
    pattern_density = group,
    pattern_spacing = group,
    pattern_angle = group,
    alpha = group
  )) +
  geom_col_pattern(
    fill = "grey20", col = "white",
    pattern_fill = "white", pattern_color = "white"
  ) +
  facet_grid(. ~ Site) +
  scale_y_continuous(expand = c(0, 0.01), breaks = seq(0, 1, 0.2)) +
  scale_pattern_manual(values = rev(c("none", "none", "circle", "stripe", "stripe"))) +
  scale_pattern_density_manual(values = rev(c(1, 1, 0.3, 0.1, 0.1))) +
  scale_pattern_spacing_manual(values = rev(c(1, 1, 0.08, 0.03, 0.03))) +
  scale_pattern_angle_manual(values = rev(c(0, 0, 45, 0, 90))) +
  scale_alpha_manual(values = rev(c(1, rep(0.5, 4)))) +
  labs(x = "Month",
       y = "Proportion from all coral cover",
       pattern = "Health status",
       pattern_density = "Health status",
       pattern_spacing = "Health status",
       pattern_angle = "Health status",
       alpha = "Health status") +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.box = "horizontal",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 12, angle = 90,
                               vjust = 0.5, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 0.3, unit = "cm"),
                                size = 14),
    axis.title.y = element_text(margin = margin(r = 0.3, unit = "cm"),
                                size = 14),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  guides(
    pattern = guide_legend(reverse = TRUE),
    pattern_density = guide_legend(reverse = TRUE),
    pattern_spacing = guide_legend(reverse = TRUE),
    pattern_angle = guide_legend(reverse = TRUE),
    alpha = guide_legend(reverse = TRUE)
  )


## health (considering genera) ----
### proportion ----
#### betabinomial ----

# Identify model objects
model_health_specific_objects <-
  ls()[str_detect(ls(), "model_health_specific_betabinom_")]

model_genera_details <- 
  data.frame(x = str_remove_all(model_health_specific_objects,
                                "model_health_specific_betabinom_")) %>% 
  separate(x, into = c("Genus_simplified", "Health_simplified"))

# Collect predictions
for (i in seq_along(model_health_specific_objects)) {
  
  temp <-
    get(model_health_specific_objects[i]) %>% 
    emmeans::emmeans(~ Month * Site, type = "response") %>% 
    as.data.frame() %>% 
    mutate(
      Genus_simplified = model_genera_details[i, 1],
      Health_simplified = model_genera_details[i, 2]
    )
  
  if (i == 1) {
    health_specific_betabinom_results <- temp
  } else {
    health_specific_betabinom_results <-
      bind_rows(health_specific_betabinom_results, temp)
  }
}

# Plot health-by-genus
fig_health_specific <-
  health_specific_betabinom_results %>%
  group_by(Month, Site, Genus_simplified) %>% 
  mutate(total_prob = sum(prob)) %>% 
  ungroup() %>% 
  mutate(
    prob_normalised = prob / total_prob,
    Health_simplified = factor(
      Health_simplified,
      levels = rev(c("Healthy", "Pale", "Bleached", "Dead"))
    )
  ) %>% 
  filter(
    (Site == "AG" & Genus_simplified %in% c("Acropora")) | 
      (Site == "CH" & Genus_simplified %in% c("Isopora", "Pocillopora", "Porites")) | 
      (Site == "HR" & Genus_simplified %in% c("Isopora", "Pocillopora", "Porites")) |
      (Site == "NB" & Genus_simplified %in% c("Isopora", "Pocillopora", "Xenia")) |
      (Site == "North" & Genus_simplified %in% c("Acropora", "Isopora", "Porites")) |
      (Site == "PH" & Genus_simplified %in% c("Acropora", "Cladiella", "Isopora", "Xenia")) |
      (Site == "SH" & Genus_simplified %in% c("Isopora", "Cladiella", "Pocillopora", "Xenia")) |
      (Site == "Sylphs" & Genus_simplified %in% c("Pocillopora", "Porites"))
  ) %>% 
  mutate(
    Site = factor(case_when(
      Site == "AG" ~ "Acropora\nGardens",
      Site == "CH" ~ "Comets\nHole",
      Site == "HR" ~ "Horseshoe\nReef",
      Site == "NB" ~ "Neds\nBeach",
      Site == "North" ~ "North\nBay",
      Site == "PH" ~ "Potholes",
      Site == "SH" ~ "Stephens\nHole",
      Site == "Sylphs" ~ "Sylphs\nHole"
    ),
    levels = c(
      "Neds\nBeach", "Sylphs\nHole", "North\nBay",
      "Acropora\nGardens", "Stephens\nHole",
      "Comets\nHole", "Horseshoe\nReef", "Potholes"
    ))
  ) %>% 
  ggplot(aes(
    x = Month,
    y = prob_normalised,
    fill = Genus_simplified,
    pattern = Health_simplified,
    pattern_density = Health_simplified,
    pattern_spacing = Health_simplified,
    pattern_angle = Health_simplified,
    alpha = Health_simplified
  )) +
  geom_col_pattern(
    col = "white",
    pattern_fill = "white",
    pattern_color = "white"
  ) +
  facet_grid(Genus_simplified ~ Site) +
  scale_y_continuous(expand = c(0, 0.01), breaks = seq(0, 1, 0.2)) +
  scale_fill_manual(values = rev(chosen_colours[-4])) +
  scale_pattern_manual(values = rev(c("none", "none", "circle", "stripe"))) +
  scale_pattern_density_manual(values = rev(c(1, 1, 0.3, 0.1))) +
  scale_pattern_spacing_manual(values = rev(c(1, 1, 0.08, 0.03))) +
  scale_pattern_angle_manual(values = rev(c(0, 0, 45, 0))) +
  scale_alpha_manual(values = rev(c(1, 0.5, 0.5, 0.5))) +
  labs(
    x = "Month",
    y = "Proportion from all coral cover",
    fill = "Genus",
    pattern = "Health status",
    pattern_density = "Health status",
    pattern_spacing = "Health status",
    pattern_angle = "Health status",
    alpha = "Health status"
  ) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    legend.box = "horizontal",
    legend.direction = "horizontal",
    axis.text.x = element_text(size = 12, angle = 90,
                               vjust = 0.5, hjust = 1),
    axis.title.x = element_text(margin = margin(t = 0.3, unit = "cm"),
                                size = 14),
    axis.title.y = element_text(margin = margin(r = 0.3, unit = "cm"),
                                size = 14),
    strip.background = element_blank(),
    strip.text = element_text(size = 12)
  ) +
  guides(
    pattern = guide_legend(nrow = 2, reverse = TRUE),
    pattern_density = guide_legend(nrow = 2, reverse = TRUE),
    pattern_spacing = guide_legend(nrow = 2, reverse = TRUE),
    pattern_angle = guide_legend(nrow = 2, reverse = TRUE),
    alpha = guide_legend(nrow = 2, reverse = TRUE)
  )

fig_health_specific


# =============================================================
#  HEALTHY CORAL
# =============================================================

# Compute estimated marginal means for Healthy coral
Healthy_coral_cover <- emmeans(
  model_health_all_betabinom_Healthy,
  ~ Month * Site,
  type = "response"
)

# Convert EMM results to data frame for plotting
df_healthy_coral_cover <- as.data.frame(Healthy_coral_cover)

# Plot Healthy coral cover through time
Healthy_facetted_line_plot <- ggplot(
  df_healthy_coral_cover,
  aes(x = Month, y = prob, group = Site, colour = Site)
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  facet_wrap(~ Site) +
  labs(
    title = "Estimated Healthy Coral Cover Over Time",
    x = "Month",
    y = "Estimated Proportion of Healthy Coral",
    colour = "Site"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

Healthy_facetted_line_plot



# =============================================================
# PALE CORAL
# =============================================================

# Compute EMMs
Pale_coral_cover <- emmeans(
  model_health_all_betabinom_Pale,
  ~ Month * Site,
  type = "response"
)

# Convert to data frame
df_pale_coral_cover <- as.data.frame(Pale_coral_cover)

# Plot Pale coral cover
Pale_facetted_line_plot <- ggplot(
  df_pale_coral_cover,
  aes(x = Month, y = prob, group = Site, colour = Site)
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  facet_wrap(~ Site) +
  labs(
    title = "Estimated Pale Coral Cover Over Time",
    x = "Month",
    y = "Estimated Proportion of Pale Coral",
    colour = "Site"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

Pale_facetted_line_plot



# =============================================================
# BLEACHED CORAL
# =============================================================

# Compute EMMs
Bleached_coral_cover <- emmeans(
  model_health_all_betabinom_Bleached,
  ~ Month * Site,
  type = "response"
)

# Convert to data frame
df_bleached_coral_cover <- as.data.frame(Bleached_coral_cover)

# Plot Bleached coral cover
Bleached_facetted_line_plot <- ggplot(
  df_bleached_coral_cover,
  aes(x = Month, y = prob, group = Site, colour = Site)
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  facet_wrap(~ Site) +
  labs(
    title = "Estimated Bleached Coral Cover Over Time",
    x = "Month",
    y = "Estimated Proportion of Pale Coral",
    colour = "Site"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

Bleached_facetted_line_plot



# =============================================================
# DEAD CORAL
# =============================================================

# Compute EMMs
Dead_coral_cover <- emmeans(
  model_health_all_betabinom_Dead,
  ~ Month * Site,
  type = "response"
)

# Convert to data frame
df_dead_coral_cover <- as.data.frame(Dead_coral_cover)

# Plot Dead coral cover
Dead_facetted_line_plot <- ggplot(
  df_dead_coral_cover,
  aes(x = Month, y = prob, group = Site, colour = Site)
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  facet_wrap(~ Site) +
  labs(
    title = "Estimated Dead Coral Cover Over Time",
    x = "Month",
    y = "Estimated Proportion of Pale Coral",
    colour = "Site"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

Dead_facetted_line_plot



# =============================================================
# LTPM CORAL
# =============================================================

# Compute EMMs
LTPM_coral_cover <- emmeans(
  model_health_all_betabinom_LTPM,
  ~ Month * Site,
  type = "response"
)

# Convert to data frame
df_LTPM_coral_cover <- as.data.frame(LTPM_coral_cover)

# Plot LTPM coral cover
LTPM_facetted_line_plot <- ggplot(
  df_LTPM_coral_cover,
  aes(x = Month, y = prob, group = Site, colour = Site)
) +
  geom_line(size = 1) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL, ymax = asymp.UCL), width = 0.2) +
  facet_wrap(~ Site) +
  labs(
    title = "Estimated LTPM Coral Cover Over Time",
    x = "Month",
    y = "Estimated Proportion of LTPM Coral",
    colour = "Site"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text  = element_text(size = 10),
    axis.title = element_text(size = 12)
  )

LTPM_facetted_line_plot



# =============================================================
# COMBINE ALL HEALTH CATEGORIES
# =============================================================

# Add health category labels
df_healthy_coral_cover  <- df_healthy_coral_cover  %>% mutate(Health = "Healthy")
df_pale_coral_cover     <- df_pale_coral_cover     %>% mutate(Health = "Pale")
df_bleached_coral_cover <- df_bleached_coral_cover %>% mutate(Health = "Bleached")
df_dead_coral_cover     <- df_dead_coral_cover     %>% mutate(Health = "Dead")
df_LTPM_coral_cover     <- df_LTPM_coral_cover     %>% mutate(Health = "LTPM")

# Combine into single table
df_all_health <- bind_rows(
  df_healthy_coral_cover,
  df_pale_coral_cover,
  df_bleached_coral_cover,
  df_dead_coral_cover,
  df_LTPM_coral_cover
)



# =============================================================
# FILTER HEALTH CATEGORIES BY SITE 
# =============================================================

df_plot <- df_all_health %>%
  # Dead: keep only for NB
  filter(!(Health == "Dead" & Site != "NB")) %>%
  # LTPM: remove for AG and NB
  filter(!(Health == "LTPM" & Site %in% c("AG", "NB"))) %>%
  # Bleached: remove for AG
  filter(!(Health == "Bleached" & Site == "AG")) %>%
  mutate(
    Health = factor(Health,
                    levels = c("Healthy", "Pale", "Bleached", "Dead", "LTPM"))
  )



# =============================================================
# DEFINE COLOUR PALETTE
# =============================================================
pastel_colors <- c(
  "Healthy"  = "#A8E6CF",
  "Pale"     = "#FFD3B6",
  "Bleached" = "#FF8B94",
  "Dead"     = "#D3D3D3",
  "LTPM"     = "#CBA0FF"
)



# =============================================================
# FINAL MULTI-PANEL PLOT OF ALL HEALTH STATES
# =============================================================
combined_health_plot <-
  ggplot(df_plot,
         aes(x = Month,
             y = prob,
             group = Health,
             colour = Health)) +
  geom_line(size = 1.5) +
  geom_point(size = 2) +
  geom_errorbar(aes(ymin = asymp.LCL,
                    ymax = asymp.UCL),
                width = 0.15,
                size = 0.8) +
  facet_wrap(
    ~ fct_relevel(
      Site,
      "PH", "SH", "NB", "Sylphs",
      "HR", "CH", "North", "AG"
    ),
    labeller = as_labeller(
      c(
        "PH" = "Potholes",
        "SH" = "Stephens Hole",
        "NB" = "Neds Beach",
        "Sylphs" = "Sylphs Hole",
        "HR" = "Horseshoe Reef",
        "CH" = "Comets Hole",
        "North" = "North Bay",
        "AG" = "Acropora Gardens"
      )
    )
  ) +
  scale_colour_manual(values = pastel_colors) +
  labs(
    x = "Month",
    y = "Proportion of Coral Cover",
    colour = "Health Status"
  ) +
  theme_minimal() +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 12),
    legend.position = "top",
    legend.direction = "horizontal",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    legend.key.width = unit(1.2, "cm")
  )

combined_health_plot

# Save figure
ggsave(
  filename = "combined_health_plot.png",
  plot = combined_health_plot,
  width = 8,
  height = 8,
  dpi = 300,
  bg = "white"
)

# ------------------------------------------------------------
# NMDS on relative abundance (Bray–Curtis)
# ------------------------------------------------------------

lhi_nmds <- metaMDS(
  df_nmds %>% select(Acropora:Xenia),
  distance = "bray",
  k = 2,
  trymax = 100
)

# -----------------------------
# 2. Extract NMDS site scores
# -----------------------------
lhi_nmds_results <- as.data.frame(scores(lhi_nmds,
                                         display = "sites",
                                         tidy = TRUE)) %>%
  select(starts_with("NMDS")) %>%
  bind_cols(df_nmds %>% select(Month, Site)) %>%
  filter(!is.na(Site))  # remove any NAs

# -----------------------------
# 3. Set site factor levels for custom legend order
# -----------------------------
site_order <- c(
  "Acropora Gardens", "North Bay", "Potholes",       # top row
  "Stephens Hole", "Neds Beach", "Sylphs Hole",     # middle row
  "Horseshoe Reef", "Comets Hole"                   # bottom row
)

lhi_nmds_results <- lhi_nmds_results %>%
  mutate(Site = factor(Site, levels = site_order))

# -----------------------------
# 4. Envfit for vectors
# -----------------------------
lhi_nmds_fit <- envfit(
  lhi_nmds,
  df_nmds %>% select(Acropora:Xenia),
  permutations = 999
)

lhi_nmds_scores <- as.data.frame(scores(lhi_nmds_fit, display = "vectors")) %>%
  mutate(
    genus  = rownames(.),
    p_val  = lhi_nmds_fit$vectors$pvals,
    genus  = ifelse(genus == "other_genera", "Other genera", paste0("*", genus, "*")),
    x_label = c(-1.175, 0.105, 0.493, 0.22, 0.89, 0.355, 0.32),
    y_label = c(-0.049, 0.647, 0.538, -0.425, -0.072, -0.96, 0.863)
  )

# -----------------------------
# 5. Stress values
# -----------------------------
# Extract stress value rounded to 1 decimal
stress_val <- round(lhi_nmds$stress, 1)

# -----------------------------
# 6. Plot
# -----------------------------

# NMDS plot with stress annotation
fig_lhi_nmds <- ggplot(lhi_nmds_results, aes(x = NMDS1, y = NMDS2)) +
  geom_point(aes(colour = Site, shape = Month), size = 3, alpha = 0.7) +
  geom_mark_hull(aes(group = Site, colour = Site, fill = Site),
                 alpha = 0.3, concavity = 5,
                 expand = unit(1, "mm"), radius = unit(1, "mm")) +
  geom_segment(data = lhi_nmds_scores,
               aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.3, "cm")), colour = "black") +
  geom_richtext(data = lhi_nmds_scores,
                aes(x = x_label, y = y_label, label = genus),
                colour = "black", alpha = 0.8, size = 5,
                vjust = 0.5, label.r = unit(0.4, "lines"),
                label.padding = unit(0.2, "lines"),
                fill = "white", label.colour = "black") +
  scale_colour_manual(values = c(
    "#B54BFC", "#4378E9", "#603A40",
    "#FA5454", "#29DABA", "#191716",
    "#5EDA29", "#EC9B31"
  )) +
  scale_fill_manual(values = c(
    "#B54BFC", "#4378E9", "#603A40",
    "#FA5454", "#29DABA", "#191716",
    "#5EDA29", "#EC9B31"
  )) +
  scale_x_continuous(limits = c(-1.55, 1)) +
  theme_classic() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14, margin = margin(t = 0.3, unit = "cm"))
  ) +
  guides(
    colour = guide_legend(nrow = 3, byrow = TRUE, order = 2),
    fill = "none",
    shape = guide_legend(nrow = 3, byrow = TRUE, order = 1)
  ) +
  # Stress annotation in bottom-right corner
  annotate("text",
           x = 0.95, y = -1.45,  # adjust as needed
           label = paste0("Stress Value = ", stress_val),
           hjust = 1, vjust = 0,
           size = 6, fontface = "bold",
           colour = "black")

fig_lhi_nmds

# Save figure
ggsave(
  filename = "fig_lhi_nmds.png",
  plot = fig_lhi_nmds,
  width = 10,
  height = 8,
  dpi = 300,
  bg = "white"
)
