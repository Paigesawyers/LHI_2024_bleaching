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

df_live_coral_cover_contrasts$p.value <- ifelse(
  df_live_coral_cover_contrasts$p.value < 0.001,
  "<0.001",
  round(df_live_coral_cover_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(model_live_coral_cover)

model_fit <- data.frame(
  AIC = round(AIC(model_live_coral_cover), 3),
  BIC = round(BIC(model_live_coral_cover), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = live_coral_cover,
    "Pairwise_comparisons" = df_live_coral_cover_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Live_coral_cover_results_inc_site.xlsx",
  rowNames = FALSE
)

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

df_live_coral_cover_month_contrasts$p.value <- ifelse(
  df_live_coral_cover_month_contrasts$p.value < 0.001,
  "<0.001",
  round(df_live_coral_cover_month_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(model_live_coral_cover)

model_fit <- data.frame(
  AIC = round(AIC(model_live_coral_cover), 3),
  BIC = round(BIC(model_live_coral_cover), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = live_coral_cover_month,
    "Pairwise_comparisons" = df_live_coral_cover_month_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Live_coral_cover_results.xlsx",
  rowNames = FALSE
)
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
Acro_cover <- emmeans(Model_genera_betabinom_Acropora_modified, ~ Month, type = "response")


# Convert EMMs to a table for summaries/plots.
df_Acro_cover  <- as.data.frame(Acro_cover)

# Pairwise comparisons across Month × Site.
Acro_cover_contrasts <- contrast(Acro_cover, method = "pairwise")
df_Acro_cover_contrasts <- as.data.frame(Acro_cover_contrasts)


Model_genera_betabinom_Acropora<- glmmTMB(
  cbind(Acropora, Cover_total - Acropora) ~ Month * Site + (1 | Site_transect),
  data = df_genera_count,
  family = "betabinomial"
)

# Extract predicted Acropora proportion by Month.
Acro_cover <- emmeans(Model_genera_betabinom_Acropora, ~ Month * Site, type = "response")

# Convert EMMs to a table for summaries/plots.
df_Acro_cover  <- as.data.frame(Acro_cover)

# Pairwise comparisons across Month × Site.
Acro_cover_contrasts <- contrast(Acro_cover, method = "pairwise")
df_Acro_cover_contrasts <- as.data.frame(Acro_cover_contrasts)

#-------------------------------------------------------------
#SITE ONLY 
# used to determine abundance and differences between site
#----------- -------------------------------------------------

Model_genera_betabinom_Acropora<- glmmTMB(
  cbind(Acropora, Cover_total - Acropora) ~ Site + (1 | Site_transect),
  data = df_genera_count,
  family = "betabinomial"
)

# Extract predicted Acropora proportion by Month.
Acro_cover <- emmeans(Model_genera_betabinom_Acropora, ~ Site, type = "response")

# Convert EMMs to a table for summaries/plots.
df_Acro_cover  <- as.data.frame(Acro_cover)

# Pairwise comparisons across Month × Site.
Acro_cover_contrasts <- contrast(Acro_cover, method = "pairwise")
df_Acro_cover_contrasts <- as.data.frame(Acro_cover_contrasts)
df_Acro_cover_contrasts$p.value <- ifelse(
  df_Acro_cover_contrasts$p.value < 0.001,
  "<0.001",
  round(df_Acro_cover_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(Model_genera_betabinom_Acropora)

model_fit <- data.frame(
  AIC = round(AIC(Model_genera_betabinom_Acropora), 3),
  BIC = round(BIC(Model_genera_betabinom_Acropora), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = df_Acro_cover,
    "Pairwise_comparisons" = df_Acro_cover_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Acropora_cover_results.xlsx",
  rowNames = FALSE
)

## --- Pocillopora ---
# Get predicted Pocillopora proportion across Month × Site.
emmeans(model_genera_betabinom_Pocillopora, ~ Month * Site, type = "response")

Model_genera_betabinom_Poc<- glmmTMB(
  cbind(Pocillopora, Cover_total - Pocillopora) ~ Site + (1 | Site_transect),
  data = df_genera_count,
  family = "betabinomial"
)

# Extract predicted Acropora proportion by Month.
poc_cover <- emmeans(Model_genera_betabinom_Poc, ~ Site, type = "response")

# Convert EMMs to a table for summaries/plots.
df_poc_cover  <- as.data.frame(poc_cover)

# Pairwise comparisons across Month × Site.
poc_cover_contrasts <- contrast(poc_cover, method = "pairwise")
df_poc_cover_contrasts <- as.data.frame(poc_cover_contrasts)
df_poc_cover_contrasts$p.value <- ifelse(
  df_poc_cover_contrasts$p.value < 0.001,
  "<0.001",
  round(df_poc_cover_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(Model_genera_betabinom_Poc)

model_fit <- data.frame(
  AIC = round(AIC(Model_genera_betabinom_Poc), 3),
  BIC = round(BIC(Model_genera_betabinom_Poc), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = df_poc_cover,
    "Pairwise_comparisons" = df_poc_cover_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Pocillopora_cover_results.xlsx",
  rowNames = FALSE
)


## --- Porites ---
# Same as above for Porites proportional cover.
emmeans(model_genera_betabinom_Porites, ~ Month * Site, type = "response")

Model_genera_betabinom_Por<- glmmTMB(
  cbind(Porites, Cover_total - Porites) ~ Site + (1 | Site_transect),
  data = df_genera_count,
  family = "betabinomial"
)

# Extract predicted Acropora proportion by Month.
por_cover <- emmeans(Model_genera_betabinom_Por, ~ Site, type = "response")

# Convert EMMs to a table for summaries/plots.
df_por_cover  <- as.data.frame(por_cover)

# Pairwise comparisons across Month × Site.
por_cover_contrasts <- contrast(por_cover, method = "pairwise")
df_por_cover_contrasts <- as.data.frame(por_cover_contrasts)
df_por_cover_contrasts$p.value <- ifelse(
  df_por_cover_contrasts$p.value < 0.001,
  "<0.001",
  round(df_por_cover_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(Model_genera_betabinom_Por)

model_fit <- data.frame(
  AIC = round(AIC(Model_genera_betabinom_Por), 3),
  BIC = round(BIC(Model_genera_betabinom_Por), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = df_por_cover,
    "Pairwise_comparisons" = df_por_cover_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Porites_cover_results.xlsx",
  rowNames = FALSE
)


## --- Isopora ---

Model_genera_betabinom_iso<- glmmTMB(
  cbind(Isopora, Cover_total - Isopora) ~ Site + (1 | Site_transect),
  data = df_genera_count,
  family = "betabinomial"
)

# Extract predicted Acropora proportion by Month.
iso_cover <- emmeans(Model_genera_betabinom_iso, ~ Site, type = "response")

# Convert EMMs to a table for summaries/plots.
df_iso_cover  <- as.data.frame(iso_cover)

# Pairwise comparisons across Month × Site.
iso_cover_contrasts <- contrast(iso_cover, method = "pairwise")
df_iso_cover_contrasts <- as.data.frame(iso_cover_contrasts)
df_iso_cover_contrasts$p.value <- ifelse(
  df_iso_cover_contrasts$p.value < 0.001,
  "<0.001",
  round(df_iso_cover_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(Model_genera_betabinom_iso)

model_fit <- data.frame(
  AIC = round(AIC(Model_genera_betabinom_iso), 3),
  BIC = round(BIC(Model_genera_betabinom_iso), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = df_iso_cover,
    "Pairwise_comparisons" = df_iso_cover_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Isopora_cover_results.xlsx",
  rowNames = FALSE
)

## --- Cladiella---

model_genera_betabinom_Cladiella<- glmmTMB(
  cbind(Cladiella, Cover_total - Cladiella) ~ Site + (1 | Site_transect),
  data = df_genera_count,
  family = "betabinomial"
)

# Extract predicted cladiella proportion by Site.
cladiella_cover <- emmeans(model_genera_betabinom_Cladiella, ~ Site, type = "response")

# Convert EMMs to a table for summaries/plots.
df_cladiella_cover  <- as.data.frame(cladiella_cover)

# Pairwise comparisons across Month × Site.
cladiella_cover_contrasts <- contrast(cladiella_cover, method = "pairwise")
df_cladiella_cover_contrasts <- as.data.frame(cladiella_cover_contrasts)
df_cladiella_cover_contrasts$p.value <- ifelse(
  df_cladiella_cover_contrasts$p.value < 0.001,
  "<0.001",
  round(df_cladiella_cover_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(model_genera_betabinom_Cladiella)

model_fit <- data.frame(
  AIC = round(AIC(model_genera_betabinom_Cladiella), 3),
  BIC = round(BIC(model_genera_betabinom_Cladiella), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = df_cladiella_cover,
    "Pairwise_comparisons" = df_cladiella_cover_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Cladiella_cover_results.xlsx",
  rowNames = FALSE
)

## --- XENIA---

model_genera_betabinom_Xenia<- glmmTMB(
  cbind(Xenia, Cover_total - Xenia) ~ Site + (1 | Site_transect),
  data = df_genera_count,
  family = "betabinomial"
)

# Extract predicted cladiella proportion by Site.
Xenia_cover <- emmeans(model_genera_betabinom_Xenia, ~ Site, type = "response")

# Convert EMMs to a table for summaries/plots.
df_Xenia_cover  <- as.data.frame(Xenia_cover)

# Pairwise comparisons across Month × Site.
Xenia_cover_contrasts <- contrast(Xenia_cover, method = "pairwise")
df_Xenia_cover_contrasts <- as.data.frame(Xenia_cover_contrasts)
df_Xenia_cover_contrasts$p.value <- ifelse(
  df_Xenia_cover_contrasts$p.value < 0.001,
  "<0.001",
  round(df_Xenia_cover_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(model_genera_betabinom_Xenia)

model_fit <- data.frame(
  AIC = round(AIC(model_genera_betabinom_Xenia), 3),
  BIC = round(BIC(model_genera_betabinom_Xenia), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = df_Xenia_cover,
    "Pairwise_comparisons" = df_Xenia_cover_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Xenia__cover_results.xlsx",
  rowNames = FALSE
)
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

df_healthy_coral_cover_contrasts$p.value <- ifelse(
  df_healthy_coral_cover_contrasts$p.value < 0.001,
  "<0.001",
  round(df_healthy_coral_cover_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(model_health_all_betabinom_Healthy)

model_fit <- data.frame(
  AIC = round(AIC(model_health_all_betabinom_Healthy), 3),
  BIC = round(BIC(model_health_all_betabinom_Healthy), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = df_healthy_coral_cover,
    "Pairwise_comparisons" = df_healthy_coral_cover_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Healthy_cover_results.xlsx",
  rowNames = FALSE
)

## getting citations and versions for the following packages
citation("emmeans")
packageVersion("emmeans")
citation("vegan")
packageVersion("vegan")
               
# ------------------------------------------------------------
# 5. CORAL HEALTH — BLEACHED
# testing differences across sites
# ------------------------------------------------------------

model_health_all_betabinom_Bleached <- glmmTMB(
  cbind(Bleached, Cover_total - Bleached) ~
    Month * Site + (1 | Site_transect),
  data = df_health_specific_count,
  family = "betabinomial"
)

Bleached_coral_cover <- emmeans(
  model_health_all_betabinom_Bleached,
  ~ Month * Site,
  type = "response"
)

contrast(Bleached_coral_cover, method = "pairwise")

# Pairwise comparisons across Month × Site.
Bleached_coral_cover_contrasts <- contrast(Bleached_coral_cover, method = "pairwise")
df_Bleached_coral_cover_contrasts <- as.data.frame(Bleached_coral_cover_contrasts)

df_Bleached_coral_cover_contrasts$p.value <- ifelse(
  df_Bleached_coral_cover_contrasts$p.value < 0.001,
  "<0.001",
  round(df_Bleached_coral_cover_contrasts$p.value, 3)
)

# Goodness of fit
r2_vals <- performance::r2(model_health_all_betabinom_Bleached)

model_fit <- data.frame(
  AIC = round(AIC(model_health_all_betabinom_Bleached), 3),
  BIC = round(BIC(model_health_all_betabinom_Bleached), 3),
  Marginal_R2 = round(r2_vals$R2_marginal, 3),
  Conditional_R2 = round(r2_vals$R2_conditional, 3)
)

# Export all results
write.xlsx(
  list(
    "Model_fit" = model_fit,
    "Estimated_means" = Bleached_coral_cover,
    "Pairwise_comparisons" = df_Bleached_coral_cover_contrasts
  ),
  file = "/Users/paigesawyers/Desktop/Bleached_cover_results.xlsx",
  rowNames = FALSE
)
# ------------------------------------------------------------
# 6. SITE-SPECIFIC BLEACHING/MORTALITY MODELS
# These isolated models estimate bleaching proportions only at
# specific sites where taxon/site combinations warranted separate
# modelling (e.g., limited occurrences, site-specific patterns).
# use this to get bleaching prevalence
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

## bleaching mortality
model_health_all_betabinom_Sylphs_Dead_simplified <- df_health_specific_count %>%
  filter(Site == "Sylphs") %>%                         # isolate site
  glmmTMB(
    cbind(Dead, Cover_total - Dead) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Sylphs_Dead_simplified, ~ Month, type = "response")

## --- Neds only ---
model_health_all_betabinom_Neds_Bleached_simplified <- df_health_specific_count %>%
  filter(Site == "NB") %>%                         # isolate site
  glmmTMB(
    cbind(Bleached, Cover_total - Bleached) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Neds_Bleached_simplified, ~ Month, type = "response")

## bleaching mortality
model_health_all_betabinom_Neds_Dead_simplified <- df_health_specific_count %>%
  filter(Site == "NB") %>%                         # isolate site
  glmmTMB(
    cbind(Dead, Cover_total - Dead) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Neds_Dead_simplified, ~ Month, type = "response")

## --- North Bay only ---
model_health_all_betabinom_North_Bleached_simplified <- df_health_specific_count %>%
  filter(Site == "North") %>%                         # isolate site
  glmmTMB(
    cbind(Bleached, Cover_total - Bleached) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_North_Bleached_simplified, ~ Month, type = "response")

## bleaching mortality
model_health_all_betabinom_North_Dead_simplified <- df_health_specific_count %>%
  filter(Site == "North") %>%                         # isolate site
  glmmTMB(
    cbind(Dead, Cover_total - Dead) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_North_Dead_simplified, ~ Month, type = "response")

## --- Stephens Hole only ---
model_health_all_betabinom_SH_Bleached_simplified <- df_health_specific_count %>%
  filter(Site == "SH") %>%                         # isolate site
  glmmTMB(
    cbind(Bleached, Cover_total - Bleached) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_SH_Bleached_simplified, ~ Month, type = "response")

## bleaching mortality
model_health_all_betabinom_SH_Dead_simplified <- df_health_specific_count %>%
  filter(Site == "SH") %>%                         # isolate site
  glmmTMB(
    cbind(Dead, Cover_total - Dead) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_SH_Dead_simplified, ~ Month, type = "response")

## --- Comets Hole only ---
model_health_all_betabinom_CH_Bleached_simplified <- df_health_specific_count %>%
  filter(Site == "CH") %>%                         # isolate site
  glmmTMB(
    cbind(Bleached, Cover_total - Bleached) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_CH_Bleached_simplified, ~ Month, type = "response")

## bleaching mortality
model_health_all_betabinom_CH_Dead_simplified <- df_health_specific_count %>%
  filter(Site == "CH") %>%                         # isolate site
  glmmTMB(
    cbind(Dead, Cover_total - Dead) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_CH_Dead_simplified, ~ Month, type = "response")

## --- Horesehoe only ---
model_health_all_betabinom_HR_Bleached_simplified <- df_health_specific_count %>%
  filter(Site == "HR") %>%                         # isolate site
  glmmTMB(
    cbind(Bleached, Cover_total - Bleached) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_HR_Bleached_simplified, ~ Month, type = "response")

## bleaching mortality
model_health_all_betabinom_HR_Dead_simplified <- df_health_specific_count %>%
  filter(Site == "HR") %>%                         # isolate site
  glmmTMB(
    cbind(Dead, Cover_total - Dead) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_HR_Dead_simplified, ~ Month, type = "response")

## --- Potholes only ---
model_health_all_betabinom_PH_Bleached_simplified <- df_health_specific_count %>%
  filter(Site == "PH") %>%                         # isolate site
  glmmTMB(
    cbind(Bleached, Cover_total - Bleached) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_PH_Bleached_simplified, ~ Month, type = "response")

## bleaching mortality
model_health_all_betabinom_PH_Dead_simplified <- df_health_specific_count %>%
  filter(Site == "PH") %>%                         # isolate site
  glmmTMB(
    cbind(Dead, Cover_total - Dead) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_PH_Dead_simplified, ~ Month, type = "response")

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

emmeans(model_health_specific_betabinom_Pocillopora_Dead_simplified, ~ Month, type = "response")

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

df_dead_coral_cover_contrasts 

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
# 13. CLADIELLA PALING FOR SITES PH AND SH
# ------------------------------------------------------------

model_health_all_betabinom_Cladiella_Pale_PH_simplified <- df_health_specific_count %>%
  filter(Site == "PH", Genus_simplified == "Cladiella") %>%
  glmmTMB(
    cbind(Pale, Cover_total - Pale) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Cladiella_Pale_PH_simplified, ~ Month, type = "response")

model_health_all_betabinom_Cladiella_Pale_SH_simplified <- df_health_specific_count %>%
  filter(Site == "SH", Genus_simplified == "Cladiella") %>%
  glmmTMB(
    cbind(Pale, Cover_total - Pale) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Cladiella_Pale_SH_simplified, ~ Month, type = "response")
# ------------------------------------------------------------
# 15. XENIA PALING
# ------------------------------------------------------------
model_health_all_betabinom_Xenia_Pale_PH_simplified <- df_health_specific_count %>%
  filter(Site == "PH", Genus_simplified == "Xenia") %>%
  glmmTMB(
    cbind(Pale, Cover_total - Pale) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Xenia_Pale_PH_simplified, ~ Month, type = "response")

model_health_all_betabinom_Xenia_Pale_SH_simplified <- df_health_specific_count %>%
  filter(Site == "SH", Genus_simplified == "Xenia") %>%
  glmmTMB(
    cbind(Pale, Cover_total - Pale) ~ Month + (1 | Site_transect),
    data = .,
    family = "betabinomial"
  )

emmeans(model_health_all_betabinom_Xenia_Pale_SH_simplified, ~ Month, type = "response")

#-------------------------------------
#EXPORTING GLMM RESULT EXCEL DATASHEETS FOR THE MODELS ABOVE
#---------------------------------------
# ------------------------------------------------------------
# EXPORT SITE-SPECIFIC / HEALTH GLMM RESULTS
# Model fit + estimated means + pairwise comparisons
# ------------------------------------------------------------

library(emmeans)
library(openxlsx)
library(performance)
library(dplyr)

# Simple p-value cleaner
clean_pvalues <- function(df) {
  if ("p.value" %in% names(df)) {
    df$p.value <- ifelse(
      df$p.value < 0.001,
      "<0.001",
      round(df$p.value, 3)
    )
  }
  df
}

# Simple goodness-of-fit table
get_model_fit <- function(model) {
  r2_vals <- performance::r2(model)
  
  data.frame(
    AIC = round(AIC(model), 3),
    BIC = round(BIC(model), 3),
    Marginal_R2 = round(r2_vals$R2_marginal, 3),
    Conditional_R2 = round(r2_vals$R2_conditional, 3)
  )
}

# Simple export function
export_emmeans_results <- function(model, emm_formula, file_name) {
  
  emm <- emmeans(
    model,
    emm_formula,
    type = "response"
  )
  
  contrasts <- contrast(
    emm,
    method = "pairwise"
  )
  
  df_emm <- as.data.frame(emm)
  df_contrasts <- clean_pvalues(as.data.frame(contrasts))
  model_fit <- get_model_fit(model)
  
  write.xlsx(
    list(
      "Model_fit" = model_fit,
      "Estimated_means" = df_emm,
      "Pairwise_comparisons" = df_contrasts
    ),
    file = file_name,
    rowNames = FALSE
  )
}

# ------------------------------------------------------------
# SITE-SPECIFIC BLEACHING / DEAD MODELS
# ------------------------------------------------------------

export_emmeans_results(
  model_health_all_betabinom_Sylphs_Bleached_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Sylphs_Bleached_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_Sylphs_Dead_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Sylphs_Dead_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_Neds_Bleached_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Neds_Bleached_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_Neds_Dead_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Neds_Dead_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_North_Bleached_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/North_Bleached_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_North_Dead_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/North_Dead_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_SH_Bleached_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Stephens_Hole_Bleached_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_SH_Dead_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Stephens_Hole_Dead_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_CH_Bleached_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Comets_Hole_Bleached_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_CH_Dead_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Comets_Hole_Dead_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_HR_Bleached_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Horseshoe_Reef_Bleached_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_HR_Dead_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Horseshoe_Reef_Dead_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_PH_Bleached_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Potholes_Bleached_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_PH_Dead_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Potholes_Dead_results.xlsx"
)

# ------------------------------------------------------------
# ACROPORA
# ------------------------------------------------------------

export_emmeans_results(
  model_health_all_betabinom_Acropora_Bleached_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Acropora_Bleached_AG_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_Acropora_Pale_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Acropora_Pale_PH_results.xlsx"
)

# ------------------------------------------------------------
# POCILLOPORA / PORITES
# ------------------------------------------------------------

export_emmeans_results(
  model_health_specific_betabinom_Pocillopora_Dead_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Pocillopora_Dead_Neds_results.xlsx"
)

export_emmeans_results(
  model_health_specific_betabinom_Pocillopora_Bleached,
  ~ Month * Site,
  "/Users/paigesawyers/Desktop/Pocillopora_Bleached_results.xlsx"
)

export_emmeans_results(
  model_health_specific_betabinom_Porites_Bleached,
  ~ Month * Site,
  "/Users/paigesawyers/Desktop/Porites_Bleached_results.xlsx"
)

# ------------------------------------------------------------
# LTPM / DEAD / CLADIELLA / XENIA
# ------------------------------------------------------------

export_emmeans_results(
  model_health_all_betabinom_LTPM,
  ~ Month * Site,
  "/Users/paigesawyers/Desktop/LTPM_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_LTPM_modified,
  ~ Site,
  "/Users/paigesawyers/Desktop/LTPM_modified_May_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_Dead,
  ~ Month * Site,
  "/Users/paigesawyers/Desktop/Dead_coral_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_Cladiella_Pale_PH_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Cladiella_Pale_PH_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_Cladiella_Pale_SH_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Cladiella_Pale_SH_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_Xenia_Pale_PH_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Xenia_Pale_PH_results.xlsx"
)

export_emmeans_results(
  model_health_all_betabinom_Xenia_Pale_SH_simplified,
  ~ Month,
  "/Users/paigesawyers/Desktop/Xenia_Pale_SH_results.xlsx"
)
