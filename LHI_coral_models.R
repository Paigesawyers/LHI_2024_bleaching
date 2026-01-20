# -------------------------------------------------------------------------
# Load cleaned data and processed objects
# -------------------------------------------------------------------------
source("LHI_coral_load_data.R")


# =========================================================================
# MODELS
# =========================================================================

# -------------------------------------------------------------------------
# 1. Coral cover vs non-coral cover (binomial GLMM)
#    Response: proportion of coral vs all other benthos
#    Fixed effects: Month, Site, and their interaction
#    Random effect: Site_transect to account for repeated measures
# -------------------------------------------------------------------------

model_all_coral_cover <-
  glmmTMB(
    data = df_all_coral_cover,
    cbind(Coral, Other) ~
      Month * Site +
      (1 | Site_transect),
    family = "binomial"
  )

# Save model output
write_rds(
  model_all_coral_cover,
  "./models/model_all_coral_cover.rds"
)


# -------------------------------------------------------------------------
# 2. Live coral vs dead or non-coral cover (binomial GLMM)
#    Response: proportion of live coral vs all non-live components
#    Same fixed + random structure as above
# -------------------------------------------------------------------------

model_live_coral_cover <-
  glmmTMB(
    data = df_live_coral_cover,
    cbind(Live_coral, Other) ~
      Month * Site +
      (1 | Site_transect),
    family = "binomial"
  )

# Save model output
write_rds(
  model_live_coral_cover,
  "./models/model_live_coral_cover.rds"
)


# =========================================================================
# 3. Genera distribution models (beta-binomial GLMMs)
# =========================================================================
# These models estimate the proportion of each coral genus at each site
# Counts are modelled with a beta-binomial to handle overdispersion
# Loop runs one model per genus and saves each model as an .rds file
# -------------------------------------------------------------------------

for (selected in c(genera_to_analyse, "other_genera")) {
  
  df_genera_count %>% 
    rename(chosen = selected) %>% 
    glmmTMB(
      data = .,
      cbind(chosen, Cover_total - chosen) ~
        Month * Site +
        (1 | Site_transect),
      family = "betabinomial"
    ) %>% 
    write_rds(
      paste0(
        "./models/model_genera_betabinom_",
        selected,
        ".rds"
      )
    )
}


# =========================================================================
# 4. Health categories (all genera combined)
# =========================================================================
# For each health category (Bleached, Dead, Healthy, etc.)
# Fit a beta-binomial GLMM modelling that category vs all other states
# -------------------------------------------------------------------------

for (selected in c("Bleached", "Dead", "Healthy", "LTPM", "Pale")) {
  
  df_health_all_count %>% 
    rename(chosen = selected) %>% 
    glmmTMB(
      data = .,
      cbind(chosen, Cover_total - chosen) ~
        Month * Site +
        (1 | Site_transect),
      family = "betabinomial"
    ) %>% 
    write_rds(
      paste0(
        "./models/model_health_all_betabinom_",
        selected,
        ".rds"
      )
    )
}


# =========================================================================
# 5. Health categories by genus (beta-binomial GLMMs)
# =========================================================================
# For each genus:
#   - Filter the dataset to that genus
#   - Model each health state separately
# This captures health patterns within specific coral groups
# -------------------------------------------------------------------------

run.model_health_specific_betabinom <- function(genus) {
  
  # Subset data for the focal genus
  df_selected <-
    df_health_specific_count %>% 
    filter(Genus_simplified == genus)
  
  # Loop through health categories
  for (selected in c("Bleached", "Dead", "Healthy", "Pale")) {
    
    df_selected %>% 
      rename(chosen = selected) %>% 
      glmmTMB(
        data = .,
        cbind(chosen, Cover_total - chosen) ~
          Month * Site +
          (1 | Site_transect),
        family = "betabinomial"
      ) %>% 
      write_rds(
        paste0(
          "./models/model_health_specific_betabinom_",
          genus,
          "_",
          selected,
          ".rds"
        )
      )
  }
}

# Run the above function across all genera to analyse
walk(
  genera_to_analyse,
  run.model_health_specific_betabinom
)
