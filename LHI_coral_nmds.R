# Data wrangling ----
source("LHI_coral_load_data.R")

if (!require("pacman")) {install.packages("pacman")}
pacman::p_load(bayesplot,
               brms,
               DHARMa,
               emmeans,
               ggforce,
               ggpattern,
               ggtext,
               glmmTMB,
               fastDummies,
               janitor,
               lme4,
               lmtest,
               marginaleffects,
               patchwork,
               rstan,
               tidybayes,
               tidyverse,
               vegan,
               viridis)
# Relative abundance with Bray-Curtis dissimilarity ----
lhi_nmds <-
  metaMDS(df_nmds %>% 
            select(Acropora:Xenia), 
          distance = "bray",
          k = 2, 
          trymax = 100)

lhi_nmds_results <- 
  as.data.frame(scores(lhi_nmds,
                       display = "sites",
                       tidy = T)) %>%
  select(starts_with("nmds")) %>% 
  bind_cols(df_nmds %>% 
              select(Month,
                     Site))

lhi_nmds_fit <- envfit(lhi_nmds, 
                       df_nmds %>% 
                         select(Acropora:Xenia),
                       permutations = 999)

lhi_nmds_scores <- 
  as.data.frame(scores(lhi_nmds_fit, 
                       display = "vectors")) %>% 
  mutate(genus  = rownames(.),
         p_val = lhi_nmds_fit$vectors$pvals) %>% 
  mutate(genus = ifelse(genus == "other_genera",
                        "Other genera",
                        paste0("*",
                               genus,
                               "*")),
         x_label = c(- 1.175,
                     0.105,
                     0.493,
                     0.22,
                     0.89,
                     0.355,
                     0.32),
         y_label = c(- 0.049,
                     0.647,
                     0.538,
                     - 0.425,
                     - 0.072,
                     - 0.96,
                     0.863))

fig_lhi_nmds <-
  lhi_nmds_results %>% 
  ggplot(aes(x = NMDS1,
             y = NMDS2)) +
  geom_point(aes(col = Site,
                 shape = Month,
                 group = Site),
             size = 3,
             alpha = 0.7) +
  geom_mark_hull(aes(group = Site,
                     col = Site,
                     fill = Site),
                 alpha = 0.3,
                 concavity = 5,
                 expand = unit(1, 
                               "mm"),
                 radius = unit(1, 
                               "mm")) +
  geom_segment(data = lhi_nmds_scores,
               aes(x = 0, 
                   y = 0, 
                   xend = NMDS1, 
                   yend = NMDS2),
               arrow = arrow(length = unit(0.3, 
                                           "cm")), 
               color = "black") +
  geom_richtext(data = lhi_nmds_scores,
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
  scale_colour_manual(values = c("#B54BFC",
                                 "#4378E9",
                                 "#603A40",
                                 "#FA5454",
                                 "#29DABA",
                                 "#191716",
                                 "#5EDA29",
                                 "#EC9B31")) +
  scale_fill_manual(values = c("#B54BFC",
                               "#4378E9",
                               "#603A40",
                               "#FA5454",
                               "#29DABA",
                               "#191716",
                               "#5EDA29",
                               "#EC9B31")) +
  scale_x_continuous(limits = c(-1.55,
                                1)) +
  theme_classic() +
  theme(legend.position = "top",
        legend.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.box = "horizontal",
        legend.direction = "horizontal",
        axis.text = element_text(size = 12),
        axis.title.x = element_text(margin = margin(t = 0.3, 
                                                    unit = "cm"),
                                    size = 14),
        axis.title.y = element_text(margin = margin(r = 0.3, 
                                                    unit = "cm"),
                                    size = 14)) +
  guides(col = guide_legend(nrow = 3, 
                            byrow = TRUE,
                            order = 2),
         fill = "none",
         shape = guide_legend(nrow = 3, 
                              byrow = TRUE,
                              order = 1))

# ggsave("fig_lhi_nmds.jpg",
#        fig_lhi_nmds,
#        bg = "white",
#        dpi = 300,
#        width = 7 * 1.41 * 0.8,
#        height = 7 * 1.41 * 0.8,
#        units = "in")
