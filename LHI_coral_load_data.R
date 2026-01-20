# ============================================================
# CLEAN DATA SCRIPT
# Purpose: Load raw point-annotation data, wrangle into 
#          wide/long formats, classify coral genus & health, 
#          and produce derived datasets for cover & community analyses.
# ============================================================


# ------------------------------------------------------------
# 1. Clear the environment
#    Removes objects from memory and runs garbage collection to free RAM.
# ------------------------------------------------------------
rm(list = ls())
gc()

# ------------------------------------------------------------
# 2. Set working directory
#    Ensures file paths for reading/writing data work correctly.
# ------------------------------------------------------------

# ------------------------------------------------------------
# 3. Load required libraries
#    'pacman' manages installation and loading of packages.
# ------------------------------------------------------------
if (!require("pacman")) install.packages("pacman")

pacman::p_load(
  bayesplot, brms, DHARMa, emmeans, ggforce, ggpattern, ggtext,
  glmmTMB, fastDummies, janitor, lme4, lmtest, marginaleffects,
  patchwork, rstan, tidybayes, tidyverse, vegan, viridis
)

# Load tidyverse explicitly (redundant but harmless)
pacman::p_load(tidyverse)


# ------------------------------------------------------------
# 4. Load raw data
#    Reads in original percent-cover annotation file.
# ------------------------------------------------------------
df_raw <- read.csv("LHI_raw_data.csv")


# ============================================================
# DATA WRANGLING
# ============================================================

# ------------------------------------------------------------
# 5. Convert raw data to wide format
#    - Remove .jpg extensions
#    - Keep only confirmed annotations
#    - Extract Month, Year, Site, Transect, Quadrat
#    - Harmonize month labels and create site_transect ID
# ------------------------------------------------------------
df_wide <-
  df_raw %>% 
  mutate(Image.name = str_remove_all(Image.name, ".JPG|.jpg")) %>% 
  filter(Annotation.status == "Confirmed") %>% 
  separate(Image.name,
           into = c("Month", "Year", "Transect", "Quadrat"),
           sep = "_", extra = "merge", fill = "right") %>% 
  mutate(Month = case_when(
    Month %in% c("Jan", "Feb") ~ "January",
    Month %in% c("Feb(2)", "March") ~ "March",
    TRUE ~ Month
  )) %>%
  mutate(Month = factor(Month, levels = c("January", "March", "May"))) %>%
  separate(Transect,
           into = c("Site", "Transect"),
           sep = "(?<=[A-Za-z])(?=[0-9])",
           fill = "right") %>% 
  unite("Site_transect", c(Site, Transect), sep = "_", remove = FALSE)


# ------------------------------------------------------------
# 6. Convert wide format to long format
#    Creates rows for each annotation category per image.
# ------------------------------------------------------------
df_long <- 
  df_wide %>%
  pivot_longer(cols = ABB:Turf, 
               names_to = "Category", 
               values_to = "Cover") %>% 
  
  # --------------------------------------------------------
# 6a. Assign genus from annotation codes
#     (Extensive lookup table based on category abbreviations)
# --------------------------------------------------------
mutate(Genus = case_when(
  Category %in% c("ABB","ABH","ABMB","ABMH","ABMOD","ABMP","ABMRD",
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
  TRUE ~ NA
)) %>%
  
  # --------------------------------------------------------
# 6b. Remove unassigned categories & simplify genera groups
# --------------------------------------------------------
filter(!is.na(Genus)) %>% 
  mutate(
    Genus_simplified = case_when(
      Genus %in% c("Acanthastrea","Cyphastrea","Astrea","Montipora",
                   "Other hard coral","Paragoniastra","Seriatopora")
      ~ "other_genera",
      Genus == "Cladiella sp" ~ "Cladiella",
      Genus %in% c("Xenia cf Crassa","Other Soft Coral") ~ "Xenia",
      .default = Genus
    ),
    
    # ------------------------------------------------------
    # 6c. Assign health status (Healthy, Bleached, Dead, etc.)
    # ------------------------------------------------------
    Health = case_when(
      Category %in% c("ABH","ABMH","ACH","ANH","CYH","GOH","ISH","MEH",
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
      .default = NA
    )
  ) %>%
  
  # --------------------------------------------------------
# 6d. Derive simplified health & coral presence classes
# --------------------------------------------------------
mutate(
  Health_simplified = ifelse(str_detect(Health, "ead"), "Dead", Health),
  Live_coral = ifelse(Health %in% c("Healthy", "Bleached", "Pale"),
                      "Live_coral", "Other"),
  Any_coral = ifelse(Genus != "other_cover", "Coral", "Other")
)


# ============================================================
# 7. Summaries: Coral vs Non-coral Cover
# ============================================================

# ------------------------------------------------------------
# 7a. All coral vs other cover (wide format)
# ------------------------------------------------------------
df_all_coral_cover <- df_long %>%
  filter(Cover > 0) %>%
  group_by(Image.ID, Month, Site, Site_transect, Quadrat, Any_coral) %>%
  summarise(Cover = sum(Cover)) %>%
  ungroup() %>%
  complete(nesting(Image.ID, Month, Site, Site_transect, Quadrat),
           Any_coral, fill = list(Cover = 0)) %>%
  pivot_wider(names_from = Any_coral, values_from = Cover)


# ------------------------------------------------------------
# 7b. Live coral vs dead/other cover
# ------------------------------------------------------------
df_live_coral_cover <- df_long %>%
  filter(Cover > 0) %>%
  group_by(Image.ID, Month, Site, Site_transect, Quadrat, Live_coral) %>%
  summarise(Cover = sum(Cover)) %>%
  ungroup() %>%
  complete(nesting(Image.ID, Month, Site, Site_transect, Quadrat),
           Live_coral, fill = list(Cover = 0)) %>%
  pivot_wider(names_from = Live_coral, values_from = Cover)


# ============================================================
# 8. Genera Distribution (Proportional and Count)
# ============================================================

# ------------------------------------------------------------
# 8a. Proportional cover by simplified genus
# ------------------------------------------------------------
df_genera_prop <-
  df_long %>%
  filter(Live_coral == "Live_coral", Cover > 0) %>%
  group_by(Image.ID, Month, Site, Site_transect, Quadrat, Genus_simplified) %>%
  summarise(Cover = sum(Cover)) %>%
  ungroup() %>%
  group_by(Image.ID) %>%
  mutate(Cover_total = sum(Cover)) %>%
  ungroup() %>%
  complete(nesting(Image.ID, Month, Site, Site_transect, Quadrat, Cover_total),
           Genus_simplified, fill = list(Cover = 0)) %>%
  mutate(Prop = Cover / Cover_total,
         Prop = ifelse(Prop == 0, 1e-5, Prop)) %>%
  group_by(Image.ID) %>% mutate(Prop_total = sum(Prop)) %>% ungroup() %>%
  mutate(Prop = Prop / Prop_total) %>%
  pivot_wider(names_from = Genus_simplified, values_from = Prop)


# ------------------------------------------------------------
# 8b. Count (raw cover sums) by simplified genus
# ------------------------------------------------------------
df_genera_count <-
  df_long %>%
  filter(Live_coral == "Live_coral", Cover > 0) %>%
  group_by(Image.ID, Month, Site, Site_transect, Quadrat, Genus_simplified) %>%
  summarise(Cover = sum(Cover)) %>%
  ungroup() %>%
  group_by(Image.ID) %>% mutate(Cover_total = sum(Cover)) %>% ungroup() %>%
  complete(nesting(Image.ID, Month, Site, Site_transect, Quadrat, Cover_total),
           Genus_simplified, fill = list(Cover = 0)) %>%
  pivot_wider(names_from = Genus_simplified, values_from = Cover)


# ============================================================
# 9. Health Summaries (Ignoring Genus)
# ============================================================

# ------------------------------------------------------------
# 9a. Proportional health cover
# ------------------------------------------------------------
df_health_all_prop <- df_long %>%
  filter(!is.na(Health_simplified), Cover > 0) %>%
  group_by(Image.ID, Month, Site, Site_transect, Quadrat, Health_simplified) %>%
  summarise(Cover = sum(Cover)) %>% ungroup() %>%
  group_by(Image.ID) %>% mutate(Cover_total = sum(Cover)) %>% ungroup() %>%
  complete(nesting(Image.ID, Month, Site, Site_transect, Quadrat, Cover_total),
           Health_simplified, fill = list(Cover = 0)) %>%
  mutate(Prop = Cover / Cover_total,
         Prop = ifelse(Prop == 0, 1e-5, Prop)) %>%
  group_by(Image.ID) %>% mutate(Prop_total = sum(Prop)) %>% ungroup() %>%
  mutate(Prop = Prop / Prop_total) %>%
  pivot_wider(names_from = Health_simplified, values_from = Prop)


# ------------------------------------------------------------
# 9b. Count (raw cover) by health category
# ------------------------------------------------------------
df_health_all_count <- df_long %>%
  filter(!is.na(Health_simplified), Cover > 0) %>%
  group_by(Image.ID, Month, Site, Site_transect, Quadrat, Health_simplified) %>%
  summarise(Cover = sum(Cover)) %>% ungroup() %>%
  group_by(Image.ID) %>% mutate(Cover_total = sum(Cover)) %>% ungroup() %>%
  complete(nesting(Image.ID, Month, Site, Site_transect, Quadrat, Cover_total),
           Health_simplified, fill = list(Cover = 0)) %>%
  pivot_wider(names_from = Health_simplified, values_from = Cover)


# ============================================================
# 10. Health Summaries by Genus (Filtered and Site-Specific)
# ============================================================

# ------------------------------------------------------------
# 10a. Filter for genera with sufficient data and site relevance
# ------------------------------------------------------------
df_health_specific <- df_long %>%
  filter(
    !is.na(Health_simplified),
    Health_simplified != "LTPM",
    Genus_simplified != "other_genera",
    Cover > 0,
    
    # Site-specific genus exclusion rules
    !(Site == "AG"   & Genus_simplified != "Acropora"),
    !(Site == "CH"   & Genus_simplified %in% c("Acropora","Cladiella","Xenia")),
    !(Site == "HR"   & Genus_simplified == "Acropora"),
    !(Site == "NB"   & Genus_simplified %in% c("Acropora","Cladiella","Porites")),
    !(Site == "North" & Genus_simplified %in% c("Cladiella","Xenia")),
    !(Site == "PH"   & Genus_simplified == "Porites"),
    !(Site == "SH"   & Genus_simplified == "Porites"),
    !(Site == "Sylphs" & Genus_simplified == "Xenia")
  ) %>%
  group_by(Image.ID, Month, Site, Site_transect, Quadrat,
           Health_simplified, Genus_simplified) %>%
  summarise(Cover = sum(Cover)) %>% ungroup() %>%
  group_by(Image.ID, Genus_simplified) %>%
  mutate(Cover_total = sum(Cover)) %>% ungroup() %>%
  complete(nesting(Image.ID, Month, Site, Site_transect, Quadrat,
                   Genus_simplified, Cover_total),
           Health_simplified, fill = list(Cover = 0))


# ------------------------------------------------------------
# 10b. List of genera included in subsequent analyses
# ------------------------------------------------------------
genera_to_analyse <- df_long %>%
  filter(!Genus_simplified %in% c("other_genera","unknown","other_cover")) %>%
  pull(Genus_simplified) %>% unique()


# ------------------------------------------------------------
# 10c. Proportional health cover by genus
# ------------------------------------------------------------
df_health_specific_prop <- df_health_specific %>%
  mutate(Prop = Cover / Cover_total,
         Prop = ifelse(Prop == 0, 1e-5, Prop)) %>%
  group_by(Image.ID, Genus_simplified) %>%
  mutate(Prop_total = sum(Prop)) %>% ungroup() %>%
  mutate(Prop = Prop / Prop_total) %>%
  pivot_wider(names_from = Health_simplified, values_from = Prop)


# ------------------------------------------------------------
# 10d. Count (raw cover) by genus and health
# ------------------------------------------------------------
df_health_specific_count <- df_health_specific %>%
  pivot_wider(names_from = Health_simplified, values_from = Cover)


# ============================================================
# 11. NMDS (Non-metric Multidimensional Scaling) Data
#     Summarises genus counts by site-month for multivariate analyses.
# ============================================================
df_nmds <- df_genera_count %>%
  group_by(Month, Site, Site_transect) %>%
  summarise(across(Cover_total:Xenia, ~ sum(.))) %>%
  ungroup() %>%
  mutate(
    Site = factor(case_when(
      Site == "AG"     ~ "Acropora Gardens",
      Site == "CH"     ~ "Comets Hole",
      Site == "HR"     ~ "Horseshoe Reef",
      Site == "NB"     ~ "Neds Beach",
      Site == "North"  ~ "North Bay",
      Site == "PH"     ~ "Potholes",
      Site == "SH"     ~ "Stephens Hole",
      Site == "Sylphs" ~ "Sylphs Hole"
    ), levels = c(
      "Neds Beach","Sylphs Hole","North Bay","Acropora Gardens",
      "Stephens Hole","Comets Hole","Horseshoe Reef","Potholes"
    )),
    
    # Convert cover values to proportions within sample
    across(Acropora:Xenia, ~ . / Cover_total)
  )
