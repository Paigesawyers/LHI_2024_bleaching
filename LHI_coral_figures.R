# 1. Load cleaned data objects ----
source("LHI_coral_load_data.R")

# 2. Define site order, colours, and ggplot theme ----
chosen_site_order <- 
  c("Neds\nBeach",
    "Sylphs\nHole",
    "North\nBay",
    "Acropora\nGardens",
    "Stephens\nHole",
    "Comets\nHole",
    "Horseshoe\nReef",
    "Potholes")

# mind the chosen site order above (each colour for a different site)
chosen_colours <- 
  as.character(paletteer::paletteer_d("jcolors::rainbow"))[c(1:8)]

show_col(chosen_colours)

# alternative colour palette (for health plot)
pastel_colors <- 
  c("Healthy"  = "#A8E6CF",
    "Pale"     = "#FFD3B6",
    "Bleached" = "#FF8B94",
    "Dead"     = "#D3D3D3",
    "LTPM"     = "#CBA0FF")

show_col(pastel_colors)

# ggplot2 theme
theme_set(theme_minimal() +
            theme(axis.title.x = element_text(margin = margin(t = 0.3, 
                                                              unit = "cm"),
                                              size = 14),
                  axis.title.y = element_text(margin = margin(r = 0.3, 
                                                              unit = "cm"),
                                              size = 14),
                  axis.text = element_text(size = 12),
                  legend.title = element_text(size = 14),
                  legend.text = element_text(size = 12),
                  legend.position = "top",
                  legend.box = "horizontal",
                  legend.direction = "horizontal",
                  plot.tag = element_text(face = "bold"),
                  strip.background = element_blank(),
                  strip.text = element_text(size = 14)))

# 3. Load all model .rds files ----
rds_files <- 
  list.files(path = "./Models/",
             pattern = ".rds",
             recursive = TRUE)

for (x in rds_files) {
  assign(str_remove(x,
                    ".rds"),
         read_rds(paste0("./Models/", 
                         x)))
}

gc()   # clean up memory

# 4. Live coral vs non-coral or dead coral ----
(fig_live_coral_cover <-
   emmeans(model_live_coral_cover,
           ~ Month * Site,
           type = "response") %>%
   as.data.frame() %>%
   mutate(Site = factor(case_when(Site == "AG"    ~ "Acropora\nGardens",
                                  Site == "CH"    ~ "Comets\nHole",
                                  Site == "HR"    ~ "Horseshoe\nReef",
                                  Site == "NB"    ~ "Neds\nBeach",
                                  Site == "North" ~ "North\nBay",
                                  Site == "PH"    ~ "Potholes",
                                  Site == "SH"    ~ "Stephens\nHole",
                                  Site == "Sylphs"~ "Sylphs\nHole"),
                        levels = chosen_site_order)) %>%
   ggplot(aes(x = Month,
              y = prob,
              ymax = asymp.LCL,
              ymin = asymp.UCL)) +
   geom_col(fill = "darkgrey") +
   geom_errorbar(width = 0.2) +
   facet_grid(. ~ Site) +
   scale_y_continuous(expand = c(0,
                                 0.01),
                      limits = c(0,
                                 1),
                      breaks = seq(0,
                                   1,
                                   0.2)) +
   labs(x = "Period",
        y = "Proportion live coral cover") +
   theme(axis.text.x = element_text(size = 12,
                                    angle = 90,
                                    vjust = 0.5,
                                    hjust = 1)))

# ggsave(filename = "fig_live_coral_cover.jpg",
#        plot = fig_live_coral_cover,
#        width = 10,
#        height = 5,
#        dpi = 300)

# 5. Genera distribution ----
## NMDS ----
(fig_nmds_genera <-
   nmds_genera_results %>% 
   ggplot(aes(x = NMDS1,
              y = NMDS2)) +
   geom_point(aes(col = Site,
                  shape = Month),
              size = 3,
              alpha = 0.7) +
   geom_mark_hull(aes(col = Site,
                      fill = Site),
                  alpha = 0.3,
                  concavity = 5,
                  expand = unit(1, 
                                "mm"),
                  radius = unit(1, 
                                "mm")) +
   geom_segment(data = nmds_genera_scores,
                aes(x = 0, 
                    y = 0, 
                    xend = NMDS1, 
                    yend = NMDS2),
                arrow = arrow(length = unit(0.3, 
                                            "cm")), 
                color = "black") +
   geom_richtext(data = nmds_genera_scores,
                 aes(x = x_label, 
                     y = y_label, 
                     label = genus),
                 color = "black", 
                 alpha = 0.8,
                 size = 5,
                 vjust = 0.5,
                 label.r = unit(0.4, 
                                "lines"),
                 label.padding = unit(0.2,
                                      "lines"),
                 fill = "white",
                 label.colour = "black") +
   scale_colour_manual(values = chosen_colours) +
   scale_fill_manual(values = chosen_colours) +
   scale_x_continuous(limits = c(-1.55,
                                 1)) +
   scale_linetype_manual(values = c("solid",
                                    "longdash",
                                    "dotted")) +
   guides(col = guide_legend(nrow = 3, 
                             byrow = TRUE,
                             order = 2),
          fill = "none",
          shape = guide_legend(nrow = 3, 
                               byrow = TRUE,
                               order = 1),
          linetype = "none"))

ggsave("fig_nmds_genera.jpg",
       fig_nmds_genera,
       bg = "white",
       dpi = 300,
       width = 7 * 1.41 * 0.8,
       height = 7 * 1.41 * 0.8,
       units = "in")

## GLMMs ----
# Identify model objects and clean names
model_genera_objects <-
  ls()[str_detect(ls(), 
                  "model_genera_betabinom_")]

model_genera_names <- 
  str_remove_all(model_genera_objects, 
                 "model_genera_betabinom_")

# Collect betabinomial predictions for genera
for (i in seq_along(model_genera_objects)) {
  
  temp <-
    get(model_genera_objects[i]) %>% 
    emmeans::emmeans(~ Month * Site, 
                     type = "response") %>% 
    as.data.frame() %>% 
    mutate(Genus_simplified = model_genera_names[i])
  
  if (i == 1) {
    genera_betabinom_results <- 
      temp
  } else {
    genera_betabinom_results <- 
      bind_rows(genera_betabinom_results,
                temp)
  }
}

# Plot genera distribution
(fig_genera_distribution <-
    genera_betabinom_results %>% 
    group_by(Month, 
             Site) %>% 
    mutate(total = sum(prob)) %>% 
    ungroup() %>% 
    mutate(prob_normalised = prob / total,
           Site = factor(case_when(
             Site == "AG" ~ "Acropora\nGardens",
             Site == "CH" ~ "Comets\nHole",
             Site == "HR" ~ "Horseshoe\nReef",
             Site == "NB" ~ "Neds\nBeach",
             Site == "North" ~ "North\nBay",
             Site == "PH" ~ "Potholes",
             Site == "SH" ~ "Stephens\nHole",
             Site == "Sylphs" ~ "Sylphs\nHole"),
             levels = chosen_site_order),
           Genus_simplified = factor(ifelse(Genus_simplified == "other_genera", 
                                            "Others", 
                                            paste0("*",
                                                   Genus_simplified,
                                                   "*")),
                                     levels = rev(levels(factor(ifelse(Genus_simplified == "other_genera",
                                                                       "Others",
                                                                       paste0("*",
                                                                              Genus_simplified,
                                                                              "*"))))))) %>% 
    ggplot(aes(x = Month,
               y = prob_normalised, 
               fill = Genus_simplified)) +
    geom_col() +
    facet_grid(. ~ Site) +
    # scale_x_discrete(labels = function(x) paste(x, 
    #                                             "2024")) +
    scale_fill_manual(values = chosen_colours[8:1]) +
    scale_y_continuous(expand = c(0,
                                  0.01), 
                       breaks = seq(0,
                                    1,
                                    0.2)) +
    labs(x = "Period",
         y = "Proportion from live coral cover",
         fill = "Genus") +
    theme(legend.text = element_markdown(size = 12),
          axis.text.x = element_text(size = 12,
                                     angle = 90,
                                     vjust = 0.5,
                                     hjust = 1)) +
    guides(fill = guide_legend(nrow = 1, 
                               reverse = TRUE)))

ggsave(filename = "fig_genera_distribution.png",
       plot = fig_genera_distribution,
       width = 10,
       height = 5,
       dpi = 300,
       bg = "white")

# 6. Health (across genera) ----
## NMDS ----
(fig_nmds_health <-
   nmds_health_results %>% 
   ggplot(aes(x = NMDS1,
              y = NMDS2)) +
   geom_point(aes(
     col = Site,
     shape = Month,
     group = Site),
     size = 3,
     alpha = 0.7) +
   geom_mark_hull(aes(linetype = Month,
                      fill = Site,
                      col = Site),
                  alpha = 0.3,
                  concavity = 5,
                  expand = unit(1, 
                                "mm"),
                  radius = unit(1, 
                                "mm")) +
   geom_segment(data = nmds_health_scores,
                aes(x = 0, 
                    y = 0, 
                    xend = NMDS1, 
                    yend = NMDS2),
                arrow = arrow(length = unit(0.3, 
                                            "cm")), 
                color = "black") +
   geom_richtext(data = nmds_health_scores,
                 aes(x = NMDS1, 
                     y = NMDS2, 
                     label = genus),
                 color = "black", 
                 alpha = 0.8,
                 size = 5,
                 vjust = 0.5,
                 label.r = unit(0.4, 
                                "lines"),
                 label.padding = unit(0.2,
                                      "lines"),
                 fill = "white",
                 label.colour = "black") +
   scale_colour_manual(values = chosen_colours) +
   scale_fill_manual(values = chosen_colours) +
   scale_linetype_manual(values = c("solid",
                                    "longdash",
                                    "dotted")) +
   # scale_x_continuous(limits = c(-1.55,
   #                               1)) +
   guides(col = guide_legend(nrow = 3, 
                             byrow = TRUE,
                             order = 2),
          fill = "none",
          shape = guide_legend(nrow = 3, 
                               byrow = TRUE,
                               order = 1),
          linetype = "none"))

ggsave("fig_nmds_health.jpg",
       fig_nmds_health,
       bg = "white",
       dpi = 300,
       width = 7 * 1.41 * 0.8,
       height = 7 * 1.41 * 0.8,
       units = "in")

## GLMMS ----
# Identify health model objects and names
model_health_all_objects <-
  ls()[str_detect(ls(), 
                  "model_health_all_betabinom_")]

model_health_all_names <- 
  str_remove_all(model_health_all_objects, 
                 "model_health_all_betabinom_")

# Collect predictions
for (i in seq_along(model_health_all_objects)) {
  temp <-
    get(model_health_all_objects[i]) %>% 
    emmeans::emmeans(~ Month * Site, 
                     type = "response") %>% 
    as.data.frame() %>% 
    mutate(Health_simplified = model_health_all_names[i])
  
  if (i == 1) {
    health_all_betabinom_results <- temp
  } else {
    health_all_betabinom_results <- bind_rows(health_all_betabinom_results, temp)
  }
}

df_all_health <-
  map_df(mget(model_health_all_objects),
         ~ emmeans(.,
                   ~ Month * Site,
                   type = "response") %>% 
           as.data.frame()) %>% 
  mutate(Health = factor(rep(c("Bleached",
                               "Dead",
                               "Healthy",
                               "LTPM",
                               "Pale"),
                             each = 24),
                         levels = c("Healthy", 
                                    "Pale", 
                                    "Bleached", 
                                    "Dead", 
                                    "LTPM")),
         Site = factor(case_when(
           Site == "AG" ~ "Acropora\nGardens",
           Site == "CH" ~ "Comets\nHole",
           Site == "HR" ~ "Horseshoe\nReef",
           Site == "NB" ~ "Neds\nBeach",
           Site == "North" ~ "North\nBay",
           Site == "PH" ~ "Potholes",
           Site == "SH" ~ "Stephens\nHole",
           Site == "Sylphs" ~ "Sylphs\nHole"),
           levels = chosen_site_order),
         asymp.UCL = ifelse(asymp.UCL > 0.9999,
                            NaN,
                            asymp.UCL))

(combined_health_plot <-
    ggplot(df_all_health,
           aes(x = Month,
               y = prob,
               group = Health,
               colour = Health)) +
    geom_line(size = 1.2) +
    geom_point(size = 1.2) +
    geom_errorbar(aes(ymin = asymp.LCL,
                      ymax = asymp.UCL,),
                  width = 0.15,
                  size = 0.8) +
    facet_wrap(vars(Site)) +
    scale_colour_manual(values = rev(pastel_colors)) +
    labs(x = "Period",
         y = "Proportion of Coral Cover",
         linetype = "Health status",
         colour = "Site") +
    theme(axis.text.x = element_text(size = 12,
                                     angle = 90,
                                     hjust = 1,
                                     vjust = 0.5)))

ggsave(filename = "combined_health_plot.png",
       plot = combined_health_plot,
       width = 8,
       height = 8,
       dpi = 300,
       bg = "white")

# 7. Health status (within genera) ----
# Identify model objects
model_health_specific_objects <-
  ls()[str_detect(ls(), 
                  "model_health_specific_betabinom_")]

model_genera_details <- 
  data.frame(x = str_remove_all(model_health_specific_objects,
                                "model_health_specific_betabinom_")) %>% 
  separate(x,
           into = c("Genus_simplified", 
                    "Health_simplified"))

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

# plot health by genus
(fig_health_specific <-
  health_specific_betabinom_results %>%
  group_by(Month, 
           Site, 
           Genus_simplified) %>% 
  mutate(total_prob = sum(prob)) %>% 
  ungroup() %>% 
  mutate(prob_normalised = prob / total_prob,
         Health_simplified = factor(
           Health_simplified,
           levels = rev(c("Healthy",
                          "Pale",
                          "Bleached",
                          "Dead")))) %>% 
  filter((Site == "AG" & Genus_simplified %in% c("Acropora")) | 
           (Site == "CH" & Genus_simplified %in% c("Isopora", "Pocillopora", "Porites")) | 
           (Site == "HR" & Genus_simplified %in% c("Isopora", "Pocillopora", "Porites")) |
           (Site == "NB" & Genus_simplified %in% c("Isopora", "Pocillopora", "Xenia")) |
           (Site == "North" & Genus_simplified %in% c("Acropora", "Isopora", "Porites")) |
           (Site == "PH" & Genus_simplified %in% c("Acropora", "Cladiella", "Isopora", "Xenia")) |
           (Site == "SH" & Genus_simplified %in% c("Isopora", "Cladiella", "Pocillopora", "Xenia")) |
           (Site == "Sylphs" & Genus_simplified %in% c("Pocillopora", "Porites"))) %>% 
  mutate(Site = factor(case_when(Site == "AG" ~ "Acropora\nGardens",
                                 Site == "CH" ~ "Comets\nHole",
                                 Site == "HR" ~ "Horseshoe\nReef",
                                 Site == "NB" ~ "Neds\nBeach",
                                 Site == "North" ~ "North\nBay",
                                 Site == "PH" ~ "Potholes",
                                 Site == "SH" ~ "Stephens\nHole",
                                 Site == "Sylphs" ~ "Sylphs\nHole"),
                       levels = c("Neds\nBeach", "Sylphs\nHole", "North\nBay",
                                  "Acropora\nGardens", "Stephens\nHole",
                                  "Comets\nHole", "Horseshoe\nReef", "Potholes"))) %>% 
  ggplot(aes(x = Month,
             y = prob_normalised,
             fill = Site,
             pattern = Health_simplified,
             pattern_density = Health_simplified,
             pattern_spacing = Health_simplified,
             pattern_angle = Health_simplified,
             alpha = Health_simplified)) +
  geom_col_pattern(col = "white",
                   pattern_fill = "white",
                   pattern_color = "white") +
  facet_grid(Genus_simplified ~ Site) +
  # scale_x_discrete(labels = function(x) paste(x, "2024")) +
  scale_y_continuous(expand = expansion(mult = c(0.08, 0.10)),
                     breaks = c(0, 0.5, 1),
                     limits = c(0, 1)) +
scale_fill_manual(values = chosen_colours) +
  scale_pattern_manual(values = rev(c("none", "none", "circle", "stripe"))) +
  scale_pattern_density_manual(values = rev(c(1, 1, 0.3, 0.1))) +
  scale_pattern_spacing_manual(values = rev(c(1, 1, 0.08, 0.03))) +
  scale_pattern_angle_manual(values = rev(c(0, 0, 45, 0))) +
  scale_alpha_manual(values = rev(c(1, 0.5, 0.5, 0.5))) +
  labs(x = "Period",
       y = "Proportion from all coral cover",
       fill = "Genus",
       pattern = "Health status",
       pattern_density = "Health status",
       pattern_spacing = "Health status",
       pattern_angle = "Health status",
       alpha = "Health status") +
  theme(axis.text.x = element_text(size = 12,
                                   angle = 90,
                                   vjust = 0.5,
                                   hjust = 1),
        strip.text.y = element_text(size = 12,
                                    face = "italic")) +
  guides(pattern = guide_legend(reverse = TRUE),
         pattern_density = guide_legend(reverse = TRUE),
         pattern_spacing = guide_legend(reverse = TRUE),
         pattern_angle = guide_legend(reverse = TRUE),
         alpha = guide_legend(reverse = TRUE),
         fill = "none"))

ggsave(filename = "fig_health_specific.png",
       plot = fig_health_specific,
       width = 10,
       height = 10,
       dpi = 300,
       bg = "white")