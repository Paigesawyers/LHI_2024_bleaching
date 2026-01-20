# LHI_2024_bleaching
This repository contains all scripts, figures, raw data, licence, and documentation for the study Bleaching at the World Heritage–listed Lord Howe Island Marine Park, the southernmost coral reef ecosystem, during the 2024 global bleaching event.

# Authors
Paige Sawyers

# Contact
For questions or collaboration requests, contact Paige Sawyers at [p.sawyers@unsw.edu.au].


## Overview
This project analyses the 2024 marine heatwave at Lord Howe Island (LHI) 
---

## Overview
This project analyses the 2024 marine heatwave at Lord Howe Island (LHI).

## Repository Structure

/LHI_raw_data                  - Original raw CSV files from CoralNet  
/LHI_Coralnet_labelset        - Annotation codes
/Models                       - Saved GLMM and beta-binomial model objects (.rds)  
/Figures                      - Generated plots and visualisations  

Scripts:  
- LHI_clean_load_data_cleaned.R  
- LHI_coral_models_cleaned.R  
- LHI_statistics_cleaned.R  
- LHI_coral_figures_cleaned.R  
- LHI_SST_temp_cleaned.R  

README.md                  - This file

# Instructions
 To run this repository, first download and unzip the entire repository.

 To reproduce the study:
 1. Install all required R packages listed below.
   - You can install missing packages using install.packages("package_name") or pacman::p_load().
   
 2. Set your working directory to the main project folder:
      setwd("path/to/repository")

 3. Run the scripts in order to reproduce the full analysis pipeline:
    a) Data Cleaning & Wrangling:
       - Script: LHI_clean_load_data_cleaned.R
       - Produces cleaned datasets for coral cover, genera proportions, and health.

    b) Modelling:
       - Script: LHI_coral_models_cleaned.R
       - Fits GLMMs and beta-binomial models for coral cover and health.
       - Saves model objects in ./models/ for downstream use.

    c) Statistical Analysis:
       - Script: LHI_statistics_cleaned.R
       - Loads fitted models, extracts EMMs, computes contrasts, and generates summary tables.

    d) Figures:
       - Script: LHI_coral_figures_cleaned.R
       - Generates visualisations of coral cover, genera distributions, health categories, and NMDS ordinations.

    e) SST & DHW Analysis:
       - Script: LHI_SST_temp_cleaned.R
       - Processes in situ and satellite SST data, computes anomalies, DHW, and plots NOAA-style facetted figures.

 4. Optional:
    - Export plots or summary tables by uncommenting the ggsave or write_csv lines in each script.

# Notes:
 - Scripts are designed to be run sequentially; intermediate outputs are used by downstream scripts.
 - Ensure raw data files are in the correct directories as expected by each script.
 - All figures and model outputs are reproducible from the cleaned datasets and stored models.

# Packages Used
The following R packages are required across all LHI analysis scripts:
1. tidyverse  : Data wrangling, manipulation, and plotting (dplyr, ggplot2, tidyr, stringr, readr, etc.)
2. janitor    : Clean and standardise column names
3. lubridate  : Date/time parsing and handling
4. zoo        : Rolling calculations (Degree Heating Weeks)
5. data.table : Fast data manipulation and run-length encoding
6. patchwork  : Combine multiple ggplot figures
7. viridis    : Colour palettes for plots
8. scales     : Axis scaling and formatting in ggplot
9. cowplot    : Plot annotations and combined layouts
10. paletteer : Access extended colour palettes
11. glmmTMB   : Fit GLMMs including beta-binomial models
12. emmeans   : Estimated Marginal Means and pairwise contrasts

# Notes:
- tidyverse already loads core packages like ggplot2, dplyr, tidyr, stringr, and readr.
- Some packages are only used in specific scripts:
- zoo & data.table: SST/DHW rolling calculations
- patchwork, viridis, scales, cowplot, paletteer: Figure plotting
- glmmTMB & emmeans: Models and statistical analyses
- Ensure all packages are installed before running the full pipeline.

## Data
This repository uses two primary CSV files from CoralNet annotations:

### `LHI_raw_data.csv`  
Raw point-annotation data exported from CoralNet. Each row represents a single point annotation per image.

- `Site`                 - Reef site name  
- `Transect`             - Transect number within the site  
- `Quadrat`              - Quadrat number within the transect  
- `Image_ID`             - Identifier for each photo  
- `Month`                - Month of survey  
- `Year`                 - Year of survey  
- `Annotation_code`      - CoralNet annotation code for the point  

### `LHI_Coralnet_labelset.csv`  
Mapping of CoralNet annotation codes to coral genera and health categories.  
- `Code` - CoralNet annotation code  
- `Hard_coral`       - Full label of the annotation 
- `Soft_coral`       - Full label of the annotation 
- `Other`            - Full label of the annotation 
 
# =============================================================================
# Data Cleaning and Wrangling Script
# =============================================================================
 Script: LHI_clean_load_data_cleaned.R
  Purpose:
  Process raw coral point-annotation data, standardise formats, classify coral
  genus and health, and produce derived datasets for coral cover, community,
  and health analyses.

# Key Steps:
 1. Environment Setup
   - Clears workspace and runs garbage collection.
   - Sets working directory.
   - Loads required R packages using pacman.

 2. Raw Data Loading
   - Reads original point-annotation file: LHI_percent_covers_V10.csv
   
 3. Data Wrangling
    - Wide Format: removes file extensions, filters confirmed annotations,
      extracts Month, Year, Site, Transect, Quadrat, harmonises month labels,
      creates Site_transect identifier.
    - Long Format: converts to row-per-category format.
    - Genus Assignment: maps annotation codes to coral genera.
    - Health Classification: assigns each annotation to Healthy, Bleached, Pale,
      Dead, Old_dead, Recent_dead, LTPM, or Unknown.
    - Simplified Categories:
       - Genus_simplified: groups minor taxa under other_genera, merges soft corals
       - Health_simplified: collapses multiple dead classes into Dead
       - Live_coral & Any_coral flags for presence/absence

 4. Derived Summaries
    - Coral vs Non-Coral Cover: wide-format datasets per image
    - Genera Summaries:
       - Proportional cover: df_genera_prop
       - Raw counts: df_genera_count
    - Health Summaries:
       - Overall proportional: df_health_all_prop
       - Count-based: df_health_all_count
       - Genus-specific, site-filtered: df_health_specific_prop & df_health_specific_count

 5. NMDS Data Preparation
    - Aggregates genus counts by Site × Month × Transect for multivariate analyses
    - Converts counts to relative proportions per sample

# Outputs:
   df_wide                    - Cleaned wide-format annotation table
   df_long                    - Long-format table with genus and health
   df_all_coral_cover         - Coral vs non-coral cover per image
   df_live_coral_cover        - Live coral vs dead/other per image
   df_genera_prop             - Proportional cover by simplified genus
   df_genera_count            - Raw counts by simplified genus
   df_health_all_prop         - Proportional health cover (all genera)
   df_health_all_count        - Count-based health cover (all genera)
   df_health_specific_prop    - Proportional health cover by genus (site-filtered)
   df_health_specific_count   - Count-based health cover by genus (site-filtered)
   df_nmds                    - NMDS-ready site × month community proportions

# Notes:
 - Categories with no assigned genus or insufficient data are removed.
 - Health and genus simplifications are designed for downstream GLMMs and
   community analyses.
 - Outputs are used in subsequent scripts for coral cover modelling, genera-
   specific analyses, and health contrasts.

# =============================================================================
# Models
# =============================================================================
 Script: LHI_coral_models_cleaned.R

 Purpose:
   Fit statistical models to processed coral cover and health data, including:
   - Coral vs non-coral cover
   - Live coral vs dead/other cover
   - Genus-specific proportional cover
   - Health categories (Bleached, Dead, Healthy, Pale, LTPM) overall
   - Health categories by genus
   Models are saved as .rds files for downstream extraction of estimated
   marginal means and contrasts.

# Workflow:
 1. Load cleaned data
    - Source "LHI_coral_load_data.R" to import processed datasets
    - Datasets include df_all_coral_cover, df_live_coral_cover, df_genera_count,
      df_health_all_count, df_health_specific_count, and genera_to_analyse.

 2. Coral Cover Models
    a) All coral vs non-coral (binomial GLMM)
       - Response: proportion of coral points
       - Fixed effects: Month, Site, Month*Site interaction
       - Random effect: Site_transect
    b) Live coral vs non-live (binomial GLMM)
       - Same structure as above

 3. Genera Distribution Models (beta-binomial GLMMs)
    - Estimate proportional cover of each genus across Month × Site
    - Handle overdispersion using beta-binomial family
    - One model per genus; outputs saved individually

 4. Health Category Models (all genera combined)
    - Fit separate beta-binomial GLMM for each health state:
      Bleached, Dead, Healthy, Pale, LTPM
    - Models predict proportion of that category vs all others

 5. Health by Genus (beta-binomial GLMMs)
    - Filter dataset to a specific genus
    - Fit separate model per health category
    - Captures genus-specific patterns across Month × Site
    - Outputs saved as individual .rds files per genus and category

# Outputs:
   Models are saved in "./models/" folder with descriptive filenames:
   - model_all_coral_cover.rds
   - model_live_coral_cover.rds
   - model_genera_betabinom_<Genus>.rds
   - model_health_all_betabinom_<Health>.rds
   - model_health_specific_betabinom_<Genus>_<Health>.rds

# Notes:
 - All models use Site_transect as a random effect to account for repeated sampling.
 - Beta-binomial models are used for proportional cover to accommodate overdispersion.
 - Outputs are later used for estimated marginal means (emmeans) extraction
   and pairwise contrasts in downstream analysis scripts.

# =============================================================================
# Figures
# =============================================================================
 Script: LHI_coral_figures_cleaned.R

 Purpose:
   Generate visualisations from previously fitted coral cover and health models,
   and perform NMDS on relative abundances of coral genera. Includes:
   - Estimated marginal means (EMMs) for all coral cover, live coral, and health
     categories.
   - Proportional cover visualisations by genus and health status.
   - Multi-panel plots across months and sites.
   - NMDS ordination and environmental fitting (Bray–Curtis dissimilarity).

# Workflow:
 1. Load data and packages:
    - Cleaned data objects from "LHI_coral_load_data.R".
    - Required packages: janitor, patchwork, scales, tidyverse, viridis, zoo, readr, stringr.

 2. Define colour palettes:
    - Coral genus colours: chosen_colours
    - Health categories: pastel_colors

 3. Load all model outputs dynamically from "./models/".

 4. Coral cover visualisations:
    a) All coral vs non-coral cover (Binomial GLMM)
    b) Live coral vs non-coral or dead coral (Binomial GLMM)
    c) Reordered plots matching NMDS site order

 5. Genus distribution:
    - Extract EMMs from beta-binomial GLMMs for each genus.
    - Normalise probabilities to sum to 1 per site/month.
    - Plot proportional cover per genus with chosen colour palette.

 6. Health category visualisations:
    a) Health categories across all genera:
       - Bleached, Dead, Healthy, LTPM, Pale.
       - Normalised proportional cover.
       - Patterned bar plots and line plots.
    b) Health by genus:
       - Extract predictions from genus-specific health models.
       - Normalised by genus cover per site/month.
       - Filtered to relevant genera per site.

 7. Individual health category line plots:
    - Healthy, Pale, Bleached, Dead, LTPM.
    - Facetted by site.

 8. Combined health states plot:
    - Merge individual health category data frames.
    - Filter by site-specific exclusions (Dead/LTPM/Bleached).
    - Plot multi-panel line plot facetted by site.

 9. NMDS analysis on relative abundance of genera:
    - Bray–Curtis dissimilarity.
    - metaMDS with k=2 dimensions.
    - Environmental fitting with envfit.
    - Plot site scores, genus vectors, and convex hulls per site.

# Outputs:
   - ggplot objects for each coral cover and health category:
     fig_all_coral_cover, fig_live_coral_cover, fig_genera_distribution,
     fig_health_all, fig_health_specific
   - Line plots for individual health states:
     Healthy_facetted_line_plot, Pale_facetted_line_plot,
     Bleached_facetted_line_plot, Dead_facetted_line_plot,
     LTPM_facetted_line_plot
   - Combined multi-panel health plot: combined_health_plot
   - NMDS ordination plot: fig_lhi_nmds

# Notes:
   - All plots use EMMs from previously fitted models to show estimated proportions.
   - Colours and patterns are carefully chosen for clarity across sites and categories.
   - NMDS vectors are annotated with genus names and adjusted for readability.
   - Filtering applied to health categories per site to reflect biologically
     meaningful occurrences.

# =============================================================================
# SST & DHW Analysis
# =============================================================================
 Script: LHI_SST_temp_cleaned.R

 Purpose:
   Process in situ and satellite sea surface temperature (SST) data for Lord Howe
   Island to characterise the 2024 marine heatwave. Includes:
   - Cleaning and formatting of raw in situ logger data.
   - Calculation of daily maximum and minimum SST per site.
   - Computation of temperature anomalies relative to Maximum Monthly Mean (MMM).
   - Identification of hotspot days (SST ≥ MMM + 1°C).
   - Calculation of Degree Heating Weeks (DHW) from both in situ and satellite data.
   - Identification of peak DHW events and longest consecutive heat stress periods.
   - NOAA-style facetted visualisations for both in situ and satellite SST/DHW data.

# Workflow:
 1. Clean environment and load required packages:
    - Packages: tidyverse, lubridate, janitor, zoo, data.table, ggplot2,
      patchwork, viridis, scales, cowplot, paletteer

 2. Define paths and global parameters:
    - Working directory
    - Date ranges for in situ data
    - Satellite years for comparison
    - DHW rolling window (default 84 days / 12 weeks)

 3. Load site MMMs:
    - Clean site names and join with SST data

 4. Load raw in situ logger data:
    - Convert timestamps
    - Replace 0 values with NA

 5. Convert in situ data to long format:
    - Filter invalid sites and faulty loggers
    - Merge with MMMs

 6. Calculate daily max/min temperatures, anomalies, and hotspots

 7. Filter dataset to 2024

 8. Extract site-specific max and min temperatures and corresponding dates

 9. Compute Degree Heating Weeks (DHW):
    - Rolling 84-day window for both in situ and satellite datasets

 10. Identify peak DHW and consecutive runs:
     - Record start/end dates and duration of peak DHW events

 11. Round DHW for plotting consecutive runs

 12. Satellite data preprocessing:
     - Load satellite SST data
     - Calculate anomalies, hotspots, and filter to 2024

 13. Satellite DHW computation and peak DHW extraction

 14. NOAA-style plotting:
     - Scale DHW to match SST for plotting
     - Facetted line plots of SST, MMM, MMM+1°C, hotspots, and scaled DHW
     - Hotspot shading for SST ≥ MMM + 1°C
     - Combine in situ and satellite plots vertically using patchwork

# Outputs:
   - Data frames:
       insitu_daily, insitu_dhw, insitu_max_temp_2024, insitu_min_temp_2024
       satellite, satellite_dhw, satellite_max_temp_2024, satellite_min_temp_2024
   - Peak DHW summaries:
       peak_insitu_dhw, satellite_peak_dhw_duration
   - Plots:
       fig_insitu, fig_satellite, fig_combined (NOAA-style facetted plots)

# Notes:
   - Hotspot days are defined as SST ≥ MMM + 1°C.
   - DHW is scaled for plotting alongside SST (max 30°C-weeks).
   - Some sites may be excluded due to missing data or faulty loggers (e.g., comets_hole).
   - Rolling DHW window defaults to 84 days (12 weeks) to match NOAA methodology.
   - Facetted plots allow visual comparison across multiple reef sites.

# =============================================================================
 Statistical Analysis
# =============================================================================
 Script: LHI_statistics_cleaned.R

 Purpose:
   Load previously fitted GLMMs and extract Estimated Marginal Means (EMMs) to:
   - Quantify coral cover, live coral, and health categories across sites and months.
   - Compute proportional cover for key coral genera (Acropora, Pocillopora, Porites).
   - Analyse health states including Healthy, Pale, Bleached, Dead, and Low Tide Partial Mortality (LTPM).
   - Perform pairwise contrasts to assess statistical differences across Months and Sites.

# Workflow:
 1. Load stored model objects (.rds files) from "./models/":
    - Fitted GLMMs for coral cover, genera proportions, and health categories.
    - Automated loading assigns each model to a workspace object named after the RDS file.

 2. Live coral cover:
    - Extract EMMs for Month × Site combinations.
    - Perform pairwise contrasts for predicted proportions of live coral.

 3. Genera proportion models:
    - Acropora, Pocillopora, and Porites modelled using betabinomial GLMMs.
    - Month (and Site, if relevant) as fixed effects; Site_transect as random effect.
    - Extract predicted proportional cover via emmeans.

 4. Coral health categories:
    - Healthy: EMMs by Month × Site; pairwise contrasts.
    - Bleached: EMMs by Month × Site; site-specific models for low-occurrence taxa.
    - Low Tide Partial Mortality (LTPM): Month × Site predictions; modified May-only model for certain sites.
    - Dead: Month × Site predictions and contrasts.
    - Pale (Acropora-specific): PH-only and North-only single-site models.

 5. Site-specific and genus-specific models:
    - Isolated models for rare combinations or site-specific patterns:
       - Sylphs Bleached
       - Acropora Bleached (AG only)
       - Acropora Pale (PH only)
       - Pocillopora Dead (Neds Beach only)
       - Acropora Paling (PH and North only)

 6. Pairwise contrasts:
    - Conducted on all relevant EMMs to test differences across Month × Site combinations.
    - Exported as data frames for further plotting or reporting.

# Outputs:
   - EMM objects:
       live_coral_cover, Healthy_coral_cover, Bleached_coral_cover, Dead_coral_cover
       Acro_paling, Poc_bleached, Por_bleached, LTPM predictions
   - Data frames of contrasts:
       df_live_coral_cover_contrasts, df_healthy_coral_cover_contrasts,
       df_Poc_dead_contrasts, df_Poc_bleached_contrasts, df_Por_bleached_contrasts,
       df_LTPM_coral_cover_May_contrasts, df_dead_coral_cover_contrasts
   - Site- and genus-specific model outputs for rare combinations

# Notes:
   - Response variables are counts; proportions modelled using betabinomial GLMMs to allow overdispersion.
   - Site_transect is included as a random effect to account for repeated transect sampling.
   - Predictions (EMMs) are back-transformed to proportion scale.
   - Filtering applied for single-site or rare taxa models to ensure biologically meaningful estimates.
   - All models assume Month as fixed effect; some models include Site where relevant.


