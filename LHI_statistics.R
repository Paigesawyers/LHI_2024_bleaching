# ============================================================
# SCRIPT: Coral Cover, Genera Proportions, and Coral Health
# PURPOSE: Load fitted GLMMs, extract estimated marginal means,
#          and compute statistical contrasts for coral condition
#          and taxon-specific proportional cover.
# NOTES:
# - All response variables are counts of category-specific cover.
# - Cover_total is the total number of points per transect.
# - Proportions are modelled using betabinomial GLMMs to allow
#   overdispersion relative to a binomial model.
# - Site_transect is treated as a random effect, accommodating
#   repeated sampling of transects within sites.
# ============================================================


# ------------------------------------------------------------
# 1. LOAD STORED MODEL OBJECTS (.rds FILES)
# ------------------------------------------------------------

# Identify all RDS model files saved previously.
# These contain fitted glmmTMB model objects for later extraction.
rds_files <- list.files(
  path = "./models/",
  pattern = ".rds",
  recursive = TRUE
)

# Read each model file and assign into workspace using filename
# (minus ".rds") as the object name.
# This allows automated loading without manual naming.
for (x in rds_files) {
  assign(
    str_remove(x, ".rds"),
    read_rds(file.path("./models/", x))
  )
}


# ------------------------------------------------------------
# 2. LIVE CORAL COVER — PROPORTIONAL RESPONSE
# ------------------------------------------------------------

# Calculate Estimated Marginal Means (EMMs) for each Month × Site
# combination. EMMs represent model-adjusted predicted proportions
# of live coral cover, averaged over random effects.
live_coral_cover <- emmeans(
  model_live_coral_cover,
  ~ Month * Site,
  type = "response"       # back-transforms to proportion scale
)

# Compare Month × Site combinations using pairwise contrasts.
# This tests differences in predicted live coral cover proportions.
live_coral_cover_contrasts <- contrast(live_coral_cover, method = "pairwise")

# Convert contrast results to a table for export or reporting.
df_live_coral_cover_contrasts <- as.data.frame(live_coral_cover_contrasts)

# Calculate Estimated Marginal Means (EMMs) for each Month × Site
# combination. EMMs represent model-adjusted predicted proportions
# of live coral cover, averaged over random effects.
live_coral_cover_month <- emmeans(
  model_live_coral_cover,
  ~ Month,
  type = "response"       # back-transforms to proportion scale
)

live_coral_cover_month
# Compare Month × Site combinations using pairwise contrasts.
# This tests differences in predicted live coral cover proportions.
live_coral_cover_month_contrasts <- contrast(live_coral_cover_month, method = "pairwise")

# Convert contrast results to a table for export or reporting.
df_live_coral_cover_month_contrasts <- as.data.frame(live_coral_cover_month_contrasts)
# ------------------------------------------------------------
# 3. GENERA PROPORTION MODELS
# ------------------------------------------------------------

## --- Acropora ---
# Fit a betabinomial GLMM estimating proportion of Acropora cover.
# Response uses successes (Acropora cover) vs failures.
# Month is a fixed effect; transects nested within sites are random.
Model_genera_betabinom_Acropora_modified <- glmmTMB(
  cbind(Acropora, Cover_total - Acropora) ~ Month + (1 | Site_transect),
  data = df_genera_count,
  family = "betabinomial"
)

# Extract predicted Acropora proportion by Month.
emmeans(Model_genera_betabinom_Acropora_modified, ~ Month, type = "response")


## --- Pocillopora ---
# Get predicted Pocillopora proportion across Month × Site.
emmeans(model_genera_betabinom_Pocillopora, ~ Month * Site, type = "response")


## --- Porites ---
# Same as above for Porites proportional cover.
emmeans(model_genera_betabinom_Porites, ~ Month * Site, type = "response")


# ------------------------------------------------------------
# 4. CORAL HEALTH — HEALTHY
# ------------------------------------------------------------

# Extract predicted proportions of "healthy" coral per Month × Site.
Healthy_coral_cover <- emmeans(
  model_health_all_betabinom_Healthy,
  ~ Month * Site,
  type = "response"
)

# Convert EMMs to a table for summaries/plots.
df_healthy_coral_cover <- as.data.frame(Healthy_coral_cover)

# Pairwise comparisons across Month × Site.
Healthy_coral_cover_contrasts <- contrast(Healthy_coral_cover, method = "pairwise")
df_healthy_coral_cover_contrasts <- as.data.frame(Healthy_coral_cover_contrasts)


# ------------------------------------------------------------
# 5. CORAL HEALTH — BLEACHED
# ------------------------------------------------------------

# Predicted proportions of bleached coral across Month × Site.
Bleached_coral_cover <- emmeans(
  model_health_all_betabinom_Bleached,
  ~ Month * Site,
  type = "response"
)

# Pairwise contrasts test bleaching differences by Month × Site.
contrast(Bleached_coral_cover, method = "pairwise")


# ------------------------------------------------------------
# 6. SITE-SPECIFIC BLEACHING MODELS
# These isolated models estimate bleaching proportions only at
# specific sites where taxon/site combinations warranted separate
# modelling (e.g., limited occurrences, site-specific patterns).
# ------------------------------------------------------------

## --- Sylphs only ---
model_health_all_betabinom_Sylphs_Bleached_simplified <- df_health_specific_count %>%
  filter(Site == "Sylphs") %>%                         # isolate site
  glmmTMB(
    cbind(Bleached, Cover_total - Bleached) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Sylphs_Bleached_simplified, ~ Month, type = "response")


## --- Acropora bleaching (AG only) ---
model_health_all_betabinom_Acropora_Bleached_simplified <- df_health_specific_count %>%
  filter(Site == "AG") %>%
  glmmTMB(
    cbind(Bleached, Cover_total - Bleached) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Acropora_Bleached_simplified, ~ Month, type = "response")


## --- Acropora Pale (PH only) ---
model_health_all_betabinom_Acropora_Pale_simplified <- df_health_specific_count %>%
  filter(Site == "PH") %>%
  glmmTMB(
    cbind(Pale, Cover_total - Pale) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Acropora_Pale_simplified, ~ Month, type = "response")


# ------------------------------------------------------------
# 7. POCILLOPORA DEAD — NEDS BEACH ONLY
# ------------------------------------------------------------

# Separate model for Pocillopora mortality specific to Neds Beach.
model_health_specific_betabinom_Pocillopora_Dead_simplified <- df_health_specific_count %>%
  filter(Genus_simplified == "Pocillopora", Site == "NB") %>%
  glmmTMB(
    cbind(Dead, Cover_total - Dead) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

# Extract predicted proportions and contrasts.
Poc_dead <- emmeans(model_health_specific_betabinom_Pocillopora_Dead_simplified,
                    ~ Month, type = "response")

Poc_dead_contrasts <- contrast(Poc_dead, method = "pairwise")
df_Poc_dead_contrasts <- as.data.frame(Poc_dead_contrasts)


# ------------------------------------------------------------
# 8. CORAL MORTALITY (Low Tide Partial Mortality)
# ------------------------------------------------------------

# Predicted LTPM proportions across Month × Site.
emmeans(model_health_all_betabinom_LTPM, ~ Month * Site, type = "response")

# Compare Sites specifically within May.
LTPM_coral_cover_May <- emmeans(
  model_health_all_betabinom_LTPM,
  ~ Month,
  at = list(Month = "May"),
  type = "response"
)

LTPM_coral_cover_May_contrasts <- contrast(LTPM_coral_cover_May, method = "pairwise")
df_LTPM_coral_cover_May_contrasts <- as.data.frame(LTPM_coral_cover_May_contrasts)


# Modified model for May only, excluding sites with low counts (AG, NB).
model_health_all_betabinom_LTPM_modified <- glmmTMB(
  cbind(LTPM, Cover_total - LTPM) ~ Site + (1 | Site_transect),
  data = df_health_all_count %>% filter(!Site %in% c("AG", "NB"), Month == "May"),
  family = "betabinomial"
)

emmeans(model_health_all_betabinom_LTPM_modified, ~ Site, type = "response")


# ------------------------------------------------------------
# 9. CORAL HEALTH — DEAD
# ------------------------------------------------------------

Dead_coral_cover <- emmeans(
  model_health_all_betabinom_Dead,
  ~ Month * Site,
  type = "response"
)

Dead_coral_cover_contrasts <- contrast(Dead_coral_cover, method = "pairwise")
df_dead_coral_cover_contrasts <- as.data.frame(Dead_coral_cover_contrasts)


# ------------------------------------------------------------
# 10. POCILLOPORA BLEACHED
# ------------------------------------------------------------

Poc_bleached <- emmeans(
  model_health_specific_betabinom_Pocillopora_Bleached,
  ~ Month * Site,
  type = "response"
)

Poc_bleached_contrasts <- contrast(Poc_bleached, method = "pairwise")
df_Poc_bleached_contrasts <- as.data.frame(Poc_bleached_contrasts)

Poc_bleached <- emmeans(
  model_health_specific_betabinom_Pocillopora_Bleached,
  ~ Month,
  type = "response"
)

Poc_bleached_month_contrasts <- contrast(Poc_bleached, method = "pairwise")
df_Poc_bleached_month_contrasts <- as.data.frame(Poc_bleached_month_contrasts)

# Site-averaged bleaching per Month.
Poc_bleached_month <- emmeans(
  model_health_specific_betabinom_Pocillopora_Bleached,
  ~ Month,
  type = "response"
)

Poc_bleached_month
# ------------------------------------------------------------
# 11. PORITES BLEACHED
# ------------------------------------------------------------

Por_bleached <- emmeans(
  model_health_specific_betabinom_Porites_Bleached,
  ~ Month * Site,
  type = "response"
)
Por_bleached

Por_bleached_contrasts <- contrast(Por_bleached, method = "pairwise")
df_Por_bleached_contrasts <- as.data.frame(Por_bleached_contrasts)

Por_bleached_month <- emmeans(
  model_health_specific_betabinom_Porites_Bleached,
  ~ Month,
  type = "response"
)

Por_bleached_month_contrasts <- contrast(Por_bleached_month, method = "pairwise")
df_Por_bleached_month_contrasts <- as.data.frame(Por_bleached_month_contrasts)
# ------------------------------------------------------------
# 12. ACROPORA PALING — PH ONLY
# ------------------------------------------------------------

filtered_data <- df_health_specific_count %>%
  filter(Site == "PH") %>%
  mutate(Site = factor(Site))

model_acro_paling <- glmmTMB(
  cbind(Pale, Cover_total - Pale) ~ Month + Site + (1 | Site_transect),
  data = filtered_data,
  family = "betabinomial"
)

Acro_paling <- emmeans(model_acro_paling, ~ Month * Site, type = "response")


# ------------------------------------------------------------
# 13. CLADIELLA PALING FOR SITES PH AND SH
# ------------------------------------------------------------

model_health_all_betabinom_Cladiella_Pale_simplified <- df_health_specific_count %>%
  filter(Site == "PH", Genus_simplified == "Cladiella") %>%
  glmmTMB(
    cbind(Pale, Cover_total - Pale) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Cladiella_Pale_simplified, ~ Month, type = "response")

model_health_all_betabinom_Cladiella_Pale_simplified <- df_health_specific_count %>%
  filter(Site == "SH", Genus_simplified == "Cladiella") %>%
  glmmTMB(
    cbind(Pale, Cover_total - Pale) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Cladiella_Pale_simplified, ~ Month, type = "response")
# ------------------------------------------------------------
# 15. XENIA PALING
# ------------------------------------------------------------
model_health_all_betabinom_Xenia_Pale_simplified <- df_health_specific_count %>%
  filter(Site == "PH", Genus_simplified == "Xenia") %>%
  glmmTMB(
    cbind(Pale, Cover_total - Pale) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Xenia_Pale_simplified, ~ Month, type = "response")

model_health_all_betabinom_Xenia_Pale_simplified <- df_health_specific_count %>%
  filter(Site == "SH", Genus_simplified == "Xenia") %>%
  glmmTMB(
    cbind(Pale, Cover_total - Pale) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Xenia_Pale_simplified, ~ Month, type = "response")
