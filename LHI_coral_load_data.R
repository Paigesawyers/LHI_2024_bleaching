# 1. Clear the environment ----
# Removes objects from memory and runs garbage collection to free RAM.
rm(list = ls())
gc()

# 2. Load required libraries ----
#    'pacman' manages installation and loading of packages.
if (!require("pacman")) install.packages("pacman")
pacman::p_load(DHARMa,
               emmeans,
               ggforce,
               ggpattern,
               ggtext,
               glmmTMB,
               fastDummies,
               janitor,
               lme4,
               lmtest,
               paletteer,
               patchwork,
               scales,
               tidyverse, 
               vegan, 
               viridis,
               zoo)


# 3. Load raw data ----
#    Reads in original percent-cover annotation file.
df_raw <- read.csv("LHI_raw_data.csv")

# 4. Convert raw data to wide format ----
#    - Remove .jpg extensions
#    - Keep only confirmed annotations
#    - Extract Month, Year, Site, Transect, Quadrat
#    - Harmonize month labels and create site_transect ID

df_wide <-
  df_raw %>% 
  mutate(Image.name = str_remove_all(Image.name, 
                                     ".JPG|.jpg")) %>% 
  filter(Annotation.status == "Confirmed") %>% 
  separate(Image.name,
           into = c("Month", 
                    "Year", 
                    "Transect", 
                    "Quadrat"),
           sep = "_", extra = "merge", 
           fill = "right") %>% 
  mutate(Month = case_when(Month %in% c("Jan",
                                        "Feb") ~ "January",
                           Month %in% c("Feb(2)",
                                        "March") ~ "March",
                           TRUE ~ Month)) %>%
  mutate(Month = factor(Month, 
                        levels = c("January",
                                   "March",
                                   "May"))) %>%
  separate(Transect,
           into = c("Site", 
                    "Transect"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           fill = "right") %>% 
  unite("Site_transect", 
        c(Site, 
          Transect), 
        sep = "_", 
        remove = FALSE)

# 5. Convert wide format to long format ----
#    Creates rows for each annotation category per image.
df_long <- 
  df_wide %>%
  pivot_longer(cols = ABB:Turf, 
               names_to = "Category", 
               values_to = "Cover") %>% 
  
  ## 5a. Assign genus from annotation codes ----
mutate(Genus = case_when(Category %in% c("ABB","ABH","ABMB","ABMH","ABMOD","ABMP","ABMRD",
                                         "ABOD","ABP","ANB","ANH","ANOD","ANP","ANPPH","ANRD",
                                         "ABRD") ~ "Acropora",
                         Category %in% c("ACBB","ACDD","ACH","ACPP","ACRD") ~ "Acanthastrea",
                         Category %in% c("CYB","CYD","CYH","CYP","CYRD") ~ "Cyphastrea",
                         Category == "Astrea" ~ "Astrea",
                         Category %in% c("ISB","ISH","ISOOD","ISP","ISRD","ISORD") ~ "Isopora",
                         Category %in% c("MEB","MEDD","MEF","MEH","MEP","MOEB","MOERD","MOPRD",
                                         "MPB","MPD","MPF","MPH","MPP") ~ "Montipora",
                         Category == "OHC" ~ "Other hard coral",
                         Category %in% c("PDB","PDH","PDP","PDPPH","POF","POOD","POPHB","PORD",
                                         "PPMB","PPMF","PPMH","PPMOD","PPMP","PPMRD") ~ "Pocillopora",
                         Category %in% c("PAB","PAD","PAH","PARAP") ~ "Paragoniastra",
                         Category %in% c("POB","POH","POPPB","POPPH","PORF","PORHH","PORIBH",
                                         "POROD","PORP","PORPH","PORRD") ~ "Porites",
                         Category %in% c("SEB","SEH","SEOD","SEPL","SERD") ~ "Seriatopora",
                         Category == "OSC" ~ "Other Soft Coral",
                         Category %in% c("XccB","XccD","XccH","XccP") ~ "Xenia cf Crassa",
                         Category %in% c("CsB","CsD","CsH","CsP") ~ "Cladiella sp",
                         Category %in% c("Anemone","Ascidians","Spgs","Sand","BARK","CYANOB",
                                         "Holothuria","Rubble","CCA","M","ECH","Turf") ~ "other_cover",
                         Category %in% c("LTPM","PPMLPM") ~ "unknown",
                         TRUE ~ NA)) %>%
  
  ## 5b. Remove unassigned categories & simplify genera groups ----
filter(!is.na(Genus)) %>% 
  mutate(Genus_simplified = case_when(Genus %in% c("Acanthastrea",
                                                   "Cyphastrea",
                                                   "Astrea",
                                                   "Montipora",
                                                   "Other hard coral",
                                                   "Paragoniastra",
                                                   "Seriatopora") ~ "other_genera",
                                      Genus == "Cladiella sp" ~ "Cladiella",
                                      Genus %in% c("Xenia cf Crassa",
                                                   "Other Soft Coral") ~ "Xenia",
                                      .default = Genus),
         
         ## 5c. Assign health status (Healthy, Bleached, Dead, etc.) ----
         Health = case_when(Category %in% c("ABH","ABMH","ACH","ANH","CYH","GOH","ISH","MEH",
                                            "MPH","PAH","PDH","PPMH","SEH","CsH","XccH","POH",
                                            "Astrea","OHC","OSC","PORHH") ~ "Healthy",
                            Category %in% c("ABB","ANB","ABMB","ACBB","CYB","ISB","MEB","MPB",
                                            "PAB","PDB","POB","PPMB","SEB","XccB","CsB","MOEB",
                                            "MEF","POF","PPMF","PORF","MPF","PORIBH") ~ "Bleached",
                            Category %in% c("ABP","ANP","ABMP","ACPP","CYP","ISP","MPP","PARAP",
                                            "PDP","PORP","PPMP","SEPL","CsP","XccP","MEP","MPP",
                                            "ANPPH","MPPPH","PLPHP","PORPH") ~ "Pale",
                            Category == "LTPM" ~ "LTPM",
                            Category %in% c("CYD","AND","GOD","PAD","CsD","XccD","ACDD","MEDD","MPD")
                            ~ "Dead",
                            Category %in% c("ABMOD","ANOD","ISOOD","POOD","POROD","PPMOD","SEOD",
                                            "PDP","ABOD") ~ "Old_dead",
                            Category %in% c("ABMRD","ANRD","ISRD","MOERD","MOPRD","PORD","PORRD",
                                            "PDP","PPMRD","SERD","ABRD","ISORD","CYRD","ACRD")
                            ~ "Recent_dead",
                            .default = NA)) %>%
  
  ## 5d. Derive simplified health & coral presence classes ----
mutate(Health_simplified = ifelse(str_detect(Health, "ead"), "Dead", Health),
       Live_coral = ifelse(Health %in% c("Healthy", "Bleached", "Pale"),
                           "Live_coral", "Other"),
       Any_coral = ifelse(Genus != "other_cover", "Coral", "Other"))


# 6. Coral vs non-coral cover ----

## 6a. All coral vs other cover (wide format) ----
df_all_coral_cover <- 
  df_long %>%
  filter(Cover > 0) %>%
  group_by(Image.ID, Month, Site, Site_transect, Quadrat, Any_coral) %>%
  summarise(Cover = sum(Cover)) %>%
  ungroup() %>%
  complete(nesting(Image.ID, Month, Site, Site_transect, Quadrat),
           Any_coral, fill = list(Cover = 0)) %>%
  pivot_wider(names_from = Any_coral, values_from = Cover)

## 6b. Live coral vs dead/other cover ----
df_live_coral_cover <- 
  df_long %>%
  filter(Cover > 0) %>%
  group_by(Image.ID, Month, Site, Site_transect, Quadrat, Live_coral) %>%
  summarise(Cover = sum(Cover)) %>%
  ungroup() %>%
  complete(nesting(Image.ID, Month, Site, Site_transect, Quadrat),
           Live_coral, fill = list(Cover = 0)) %>%
  pivot_wider(names_from = Live_coral, values_from = Cover)

# 7. Genera Distribution ----
df_genera_count <-
  df_long %>%
  filter(Live_coral == "Live_coral", 
         Cover > 0) %>%
  group_by(Image.ID,
           Month,
           Site,
           Site_transect,
           Quadrat,
           Genus_simplified) %>%
  summarise(Cover = sum(Cover)) %>%
  ungroup() %>%
  group_by(Image.ID) %>% 
  mutate(Cover_total = sum(Cover)) %>%
  ungroup() %>%
  complete(nesting(Image.ID,
                   Month,
                   Site,
                   Site_transect,
                   Quadrat,
                   Cover_total),
           Genus_simplified, 
           fill = list(Cover = 0)) %>%
  pivot_wider(names_from = Genus_simplified,
              values_from = Cover)

## NMDS ----
#     Summarises genus counts by site-month for multivariate analyses.
df_nmds_genera <- 
  df_genera_count %>%
  group_by(Month, 
           Site, 
           Site_transect) %>%
  summarise(across(Cover_total:Xenia, 
                   ~ sum(.))) %>%
  ungroup() %>%
  mutate(Site = factor(case_when(Site == "AG"     ~ "Acropora Gardens",
                                 Site == "CH"     ~ "Comets Hole",
                                 Site == "HR"     ~ "Horseshoe Reef",
                                 Site == "NB"     ~ "Neds Beach",
                                 Site == "North"  ~ "North Bay",
                                 Site == "PH"     ~ "Potholes",
                                 Site == "SH"     ~ "Stephens Hole",
                                 Site == "Sylphs" ~ "Sylphs Hole"), 
                       levels = c("Neds Beach",
                                  "Sylphs Hole",
                                  "North Bay",
                                  "Acropora Gardens",
                                  "Stephens Hole",
                                  "Comets Hole",
                                  "Horseshoe Reef",
                                  "Potholes")),
         # Convert cover values to proportions within sample
         across(Acropora:Xenia, ~ . / Cover_total))

nmds_genera <-
  metaMDS(df_nmds_genera %>% 
            select(Acropora:Xenia), 
          distance = "bray",
          k = 2, 
          trymax = 100)

nmds_genera_results <- 
  as.data.frame(scores(nmds_genera,
                       display = "sites",
                       tidy = T)) %>%
  select(starts_with("nmds")) %>% 
  bind_cols(df_nmds_genera %>% 
              select(Month,
                     Site))

nmds_genera_fit <- envfit(nmds_genera, 
                          df_nmds_genera %>% 
                            select(Acropora:Xenia),
                          permutations = 999)

nmds_genera_scores <- 
  as.data.frame(scores(nmds_genera_fit, 
                       display = "vectors")) %>% 
  mutate(genus  = rownames(.),
         p_val = nmds_genera_fit$vectors$pvals) %>% 
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

stress_val <- round(nmds_genera$stress, 
                    2)

# 8. Health status (across genera) ----
df_health_all_count <- 
  df_long %>%
  filter(!is.na(Health_simplified), 
         Cover > 0) %>%
  group_by(Image.ID,
           Month,
           Site,
           Site_transect,
           Quadrat,
           Health_simplified) %>%
  summarise(Cover = sum(Cover)) %>%
  ungroup() %>%
  group_by(Image.ID) %>% 
  mutate(Cover_total = sum(Cover)) %>% 
  ungroup() %>%
  complete(nesting(Image.ID,
                   Month,
                   Site,
                   Site_transect,
                   Quadrat,
                   Cover_total),
           Health_simplified, 
           fill = list(Cover = 0)) %>%
  pivot_wider(names_from = Health_simplified, 
              values_from = Cover)

## NMDS ----
df_nmds_health <- 
  df_health_all_count %>%
  group_by(Month, 
           Site, 
           Site_transect) %>%
  summarise(across(Cover_total:Pale, 
                   ~ sum(.))) %>%
  ungroup() %>%
  mutate(Site = factor(case_when(Site == "AG"     ~ "Acropora Gardens",
                                 Site == "CH"     ~ "Comets Hole",
                                 Site == "HR"     ~ "Horseshoe Reef",
                                 Site == "NB"     ~ "Neds Beach",
                                 Site == "North"  ~ "North Bay",
                                 Site == "PH"     ~ "Potholes",
                                 Site == "SH"     ~ "Stephens Hole",
                                 Site == "Sylphs" ~ "Sylphs Hole"), 
                       levels = c("Neds Beach",
                                  "Sylphs Hole",
                                  "North Bay",
                                  "Acropora Gardens",
                                  "Stephens Hole",
                                  "Comets Hole",
                                  "Horseshoe Reef",
                                  "Potholes")),
         # Convert cover values to proportions within sample
         across(Bleached:Pale, 
                ~ . / Cover_total))

nmds_health <-
  metaMDS(df_nmds_health %>% 
            select(Bleached:Pale), 
          distance = "bray",
          k = 2, 
          trymax = 100)

nmds_health_results <- 
  as.data.frame(scores(nmds_health,
                       display = "sites",
                       tidy = T)) %>%
  select(starts_with("nmds")) %>% 
  bind_cols(df_nmds_health %>% 
              select(Month,
                     Site))

nmds_health_fit <- envfit(nmds_health, 
                          df_nmds_health %>% 
                            select(Bleached:Pale),
                          permutations = 999)

nmds_health_scores <- 
  as.data.frame(scores(nmds_health_fit, 
                       display = "vectors")) %>% 
  mutate(genus  = rownames(.),
         p_val = nmds_health_fit$vectors$pvals)

stress_val <- round(nmds_health$stress, 
                    2)


# 9. Health status (within genera) ----

## 9a. Filter for genera with sufficient data and site relevance ----
df_health_specific <- 
  df_long %>%
  filter(!is.na(Health_simplified),
         Health_simplified != "LTPM",
         Genus_simplified != "other_genera",
         Cover > 0,
         # Site-specific genus exclusion rules
         !(Site == "AG"   & Genus_simplified != "Acropora"),
         !(Site == "CH"   & Genus_simplified %in% c("Acropora",
                                                    "Cladiella",
                                                    "Xenia")),
         !(Site == "HR"   & Genus_simplified == "Acropora"),
         !(Site == "NB"   & Genus_simplified %in% c("Acropora",
                                                    "Cladiella",
                                                    "Porites")),
         !(Site == "North" & Genus_simplified %in% c("Cladiella",
                                                     "Xenia")),
         !(Site == "PH"   & Genus_simplified == "Porites"),
         !(Site == "SH"   & Genus_simplified == "Porites"),
         !(Site == "Sylphs" & Genus_simplified == "Xenia")) %>%
  group_by(Image.ID, 
           Month, 
           Site, 
           Site_transect, 
           Quadrat,
           Health_simplified, 
           Genus_simplified) %>%
  summarise(Cover = sum(Cover)) %>% 
  ungroup() %>%
  group_by(Image.ID, Genus_simplified) %>%
  mutate(Cover_total = sum(Cover)) %>% 
  ungroup() %>%
  complete(nesting(Image.ID,
                   Month,
                   Site,
                   Site_transect,
                   Quadrat,
                   Genus_simplified,
                   Cover_total),
           Health_simplified, 
           fill = list(Cover = 0))

## 9b. List of genera included in subsequent analyses ----
genera_to_analyse <- 
  df_long %>%
  filter(!Genus_simplified %in% 
           c("other_genera",
             "unknown",
             "other_cover")) %>%
  pull(Genus_simplified) %>% 
  unique()


## 9c. Count (raw cover) by genus and health ----
df_health_specific_count <- 
  df_health_specific %>%
  pivot_wider(names_from = Health_simplified, 
              values_from = Cover)

