  #'  --------------------------------
  #'  3. Relative density indices (RDI)
  #'  
  #'  R code from Bassing, S.B., D.E. Ausband, M.A. Mumma, J. Baumgardt, S. Thompson, 
  #'  M.A. Hurley, and M. Falcy. 2025. "Disentangling the web of species interactions 
  #'  in a multi-predator, multi-prey community"
  #'  --------------------------------
  #'  Script to generate relative density indices (RDI) for each species, season, 
  #'  and cluster polygon using camera-specific relative abundance indices. Then
  #'  formats data for 1-year time lag and explores/visualizes data. This script 
  #'  can be run independently but it is also sourced by 4_Structural_Equation_Models.R 
  #'  to create input data for competing SEMs.
  #'  
  #'  SPATIALLY ANONYMIZED OUTPUTS FROM 2_CLUSTER_CAMERAS.R ARE PROVIDED TO ENSURE 
  #'  CODE AND ANALYSIS ARE TRANSPARENT AND REPRODUCIBLE. CAMERA LOCATIONS USED 
  #'  TO GENERATE THESE OUTPUTS WITH 2_CLUSTER_CAMERAS.R SCRIPT ARE SENSITIVE BUT 
  #'  AVAILABLE UPON REQUEST FROM IDAHO DEPT OF FISH & GAME FOR QUALIFIED RESEARCHERS.
  #'  
  #'  Requires: 
  #'    1. Relative abundance data (RN_abundance.RData) estimated by Royle-Nichols 
  #'    abundance model in 1_Relative_abundance_Royle-Nichols_model.R
  #'    2. Cameras assigned to clusters generated in 2_Cluster_cameras.R housed in
  #'    "./Data/" folder. Data include camera ID, season, GMU, and camera setup type, 
  #'    plus GMU-specific cluster and area of cluster each camera was grouped in.
  #'        1. stations_smr2020.RData
  #'        2. stations_smr2021.RData
  #'        3. stations_smr2022.RData 
  #'    3. Shapefiles of cluster polygons generated in 2_Cluster_cameras.R housed 
  #'    in "./Outputs/Spatial_outputs/" folder
  #'        1. cluster_polygons_gmu1
  #'        2. cluster_polygons_gmu6
  #'        3. cluster_polygons.gmu10a 
  #'    4. Covariate data housed in "./Data/" folder
  #'       Covariate data includes:
  #'        1. GMU: Game Management Unit where each camera cluster occurred
  #'        2. ClusterID: Unique number of clusters within a GMU 
  #'        (ClusterID = 1 ... n for n clusters in each GMU)
  #'        3. Year: Year of data collection
  #'        4. DisturbedForest_last20Yrs: percent of cluster comprising forested 
  #'        habitat that was disturbed within the 20 years prior to year of data 
  #'        collection. Covariate generated using canopy loss data from the Global 
  #'        Forest Cover Change dataset extracted from Google Earth Engine. 
  #'        5. Annual harvest: total number of wolves harvested (trapped, hunted, 
  #'        or removed for management) reported within 1 year prior to start of 
  #'        annual data collection. Data provided by Idaho Department of Fish & Game.
  #'  --------------------------------
  
  #'  Load libraries
  library(piecewiseSEM)
  library(labelled)
  library(lme4)
  library(MASS)
  library(tidyverse)
  library(sf)
  library(terra)
  library(mapview)
  library(ggplot2)
  
  #'  Load RN model local abundance estimates
  load("./Outputs/RN_abundance.RData")
  
  #'  Load camera station data, including assigned cluster (clusters based on output
  #'  from 2_Cluster_cameras.R)
  load("./Data/stations_smr2020.RData")
  load("./Data/stations_smr2021.RData")
  load("./Data/stations_smr2022.RData")
  
  #'  Merge camera cluster data together and reduce to one observation per camera
  clusters_all <- bind_rows(stations_smr2020, stations_smr2021, stations_smr2022) %>%
    group_by(CamID) %>%
    slice(1L) %>%
    ungroup() %>%
    dplyr::select(-Season) %>%
    #'  Change GMU and Setup to characters that match thos in RN_abundance
    mutate(GMU = ifelse(GMU == 1, "GMU10A", GMU),
           GMU = ifelse(GMU == 2, "GMU6", GMU),
           GMU = ifelse(GMU == 3, "GMU1", GMU),
           Setup = ifelse(Setup == 1, "U", "P"))

  #'  Read in cluster polygon shapefiles (cluster polygons based on output from 
  #'  2_Cluster_cameras.R)
  gmu1_poly <- st_read("./Outputs/Spatial_outputs/cluster_polygons_gmu1.shp") %>% mutate(GMU = "GMU1")
  gmu6_poly <- st_read("./Outputs/Spatial_outputs/cluster_polygons_gmu6.shp") %>% mutate(GMU = "GMU6")
  gmu10a_poly <- st_read("./Outputs/Spatial_outputs/cluster_polygons_gmu10a.shp") %>% mutate(GMU = "GMU10A")

  #'  Merge cluster polygons across GMUs
  cluster_poly <- bind_rows(gmu1_poly, gmu6_poly, gmu10a_poly) %>%
    rename(ClusterID = Clusters)
  mapview::mapview(cluster_poly, zcol = "ClusterID")
  
  #'  Load GMU spatial data
  gmu_wgs84 <- st_read("./Data/Spatial_data/Study_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Load covariate data 
  load("./Data/covariates.RData")
  
  #'  ---------------------------------------------------
  ####  Calculate index of relative density per species  ####
  #'  ---------------------------------------------------
  #'  Join local abundance estimates with spatial cluster data
  cluster_RAI <- function(rai, clusters) {
    clustered_rai <- rai %>%
      left_join(clusters, by = c("CamID", "GMU", "Setup"))
    return(clustered_rai)
  }
  RN_abundance_cluster <- lapply(RN_abundance, cluster_RAI, clusters = clusters_all)
  
  #'  Calculate index of relative density per species per cluster per year
  density_per_cluster <- function(rai, yr) {
    relative_density <- rai %>%
      group_by(GMU, Clusters, Species) %>%
      reframe(SppN = sum(RN.n),
              SppDensity.km2 = SppN/Area_km2,
              SppDensity.100km2 = SppDensity.km2*100,
              SppN.r = sum(round(RN.n, 0)),
              SppDensity.km2.r = SppN.r/Area_km2,
              SppDensity.100km2.r = SppDensity.km2.r*100) %>% 
      unique() %>%
      mutate(Year = yr) %>%
      #'  Join spatial polygon data to relative density estimates
      left_join(cluster_poly, by = c("GMU", "Clusters" = "ClusterID")) %>%
      st_as_sf()
    return(relative_density)
  }
  yr <- list("2020", "2021", "2022")
  cluster_density <- mapply(density_per_cluster, rai = RN_abundance_cluster, yr = yr, SIMPLIFY = FALSE)
  
  save(cluster_density, file = "./Outputs/RelativeDensityIndex_per_SppCluster.RData")
  
  #'  ---------------------------------
  ####  Explore density relationships  ####
  #'  ---------------------------------
  #'  Visualize with ggplot
  dat <- bind_rows(cluster_density[[1]], cluster_density[[2]], cluster_density[[3]])
  (wolf_density <- dat %>% filter(Species == "wolf") %>%
      ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "blue") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (wolves/100km^2)"))
  (bear_density <- dat %>% filter(Species == "bear_black") %>%
      ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "red") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (black bear/100km^2)"))
  (lion_density <- dat %>% filter(Species == "mountain_lion") %>%
      ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "green") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (mountain lion/100km^2)"))
  (coy_density <- dat %>% filter(Species == "coyote") %>%
      ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "orange") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (coyote/100km^2)"))
  (elk_density <- dat %>% filter(Species == "elk") %>%
      ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "purple") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (elk/100km^2)"))
  (moose_density <- dat %>% filter(Species == "moose") %>%
      ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "gold") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (moose/100km^2)"))
  (wtd_density <- dat %>% filter(Species == "whitetailed_deer") %>%
      ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) +
      geom_sf(aes(fill = SppDensity.100km2.r, group = Year), alpha = 0.2) +
      scale_fill_gradient(low = "white", high = "black") +
      labs(x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) + 
      facet_wrap(~Year) +
      ggtitle("Relative density index (white-tailed deer/100km^2)"))
  
  #' #'  Save maps of relative density indices
  #' ggsave("./Outputs/cluster_wolf_density.tiff", wolf_density,
  #'        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  #' ggsave("./Outputs/cluster_bear_density.tiff", bear_density,
  #'        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  #' ggsave("./Outputs/cluster_lion_density.tiff", lion_density,
  #'        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  #' ggsave("./Outputs/cluster_coy_density.tiff", coy_density,
  #'        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  #' ggsave("./Outputs/cluster_elk_density.tiff", elk_density,
  #'        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  #' ggsave("./Outputs/cluster_moose_density.tiff", moose_density,
  #'        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  #' ggsave("./Outputs/cluster_wtd_density.tiff", wtd_density,
  #'        units = "in", width = 6, height = 4, dpi = 300, device = "tiff", compression = "lzw")
  
  #'  ---------------------------------------
  ####  Format for SEM with 1-year time lag  ####
  #'  ---------------------------------------
  #'  Long to wide data structure
  wide_data <- function(dat) {
    pivot_data_wide <- as.data.frame(dat) %>%
      filter(!is.na(Clusters)) %>%
      #'  Create categorical year variable
      mutate(Year = as.numeric(Year),
             year = ifelse(Year == 2020, "yr1", Year),
             year = ifelse(Year == 2021, "yr2", year),
             year = ifelse(Year == 2022, "yr3", year)) %>% 
      dplyr::select(c(GMU, Clusters, Year, year, Species, SppDensity.100km2.r)) %>%
      #'  Create column per species with their site-specific local abundance 
      pivot_wider(names_from = "Species",
                  values_from = c("SppDensity.100km2.r")) %>%
      left_join(covariates, by = c("GMU", "Year", "Clusters" = "ClusterID"))
    #'  Review to make sure it time lags look right given these observations are 
    #'  hypothesized to affect the next set of observations
    print(pivot_data_wide)
    return(pivot_data_wide)
  }
  density_wide <- lapply(cluster_density, wide_data)
  
  #'  Group data across years so t-1 data affects t data, regardless of year
  #'  Stack year 1 & year 2 data
  dat_t_minus_1 <- bind_rows(density_wide[[1]], density_wide[[2]])
  #'  Stack year 2 & year 3 data
  dat_t <- bind_rows(density_wide[[2]], density_wide[[3]])
  #'  List stacked data
  dat_stack_list <- list(dat_t_minus_1, dat_t)
  
  #'  Group all years together (stacking Year1, Year2, and Year3 together)
  density_wide_allyrs <- rbind(density_wide[[1]], density_wide[[2]], density_wide[[3]]) 
  
  #'  Visualize spread of relative density estimates across clusters per year
  plot_histogram <- function(dat, yr) {
    hist(dat$bear_black, main = paste("Relative bear density across clusters, \nsummer", yr))
    hist(dat$coyote, main = paste("Relative coyote density across clusters, \nsummer", yr))
    hist(dat$mountain_lion, main = paste("Relative lion density across clusters, \nsummer", yr))
    hist(dat$wolf, main = paste("Relative wolf density across clusters, \nsummer", yr))
    hist(dat$elk, main = paste("Relative elk density across clusters, \nsummer", yr))
    hist(dat$moose, main = paste("Relative moose density across clusters, \nsummer", yr))
    hist(dat$whitetailed_deer, main = paste("Relative wtd density across clusters, \nsummer", yr))
    hist(dat$DisturbedForest_last20Yrs, main = paste("Percent disturbed forest, \nsummer", yr))
    hist(dat$annual_harvest, main = paste("Wolves harvested over past 12 months, \nstarting June", yr))
  }
  plot_histogram(density_wide[[1]], yr = "2020")
  plot_histogram(density_wide[[2]], yr = "2021")
  plot_histogram(density_wide[[3]], yr = "2022")
  
  #'  ----------------------------------------------
  #####  Time period format: t-1 vs t across years  #####
  #'  ----------------------------------------------
  #'  Grouping data across years so t-1 data affects t data, regardless of year,
  #'  increases sample size and eliminates annual variation in estimated relationships
  #'  Stack year 1 & year 2 data
  dat_t_minus_1 <- dat_stack_list[[1]] %>% 
    #'  Indicate whether observation was from first or second year in grouped data
    mutate(GroupYear = ifelse(year == "yr1", "first", "second")) %>%
    #'  Add time period identifier to each column name
    rename_with(.cols = bear_black:annual_harvest, function(x){paste0(x, ".Tminus1")}) %>%
    dplyr::select(-year)
  #'  Stack year 2 & year 3 data
  dat_t <- dat_stack_list[[2]] %>% 
    #'  Indicate whether observation was from first or second year in grouped data
    mutate(GroupYear = ifelse(year == "yr2", "first", "second")) %>%
    #'  Add time period identifier to each column name
    rename_with(.cols = bear_black:annual_harvest, function(x){paste0(x, ".T")}) %>%
    dplyr::select(-year)
  
  #'  Join t-1 and t data based on camera location
  density_wide_1YrLag_20s_22s <- full_join(dat_t_minus_1, dat_t, by = c("GMU", "Clusters", "GroupYear")) %>%
    relocate(GroupYear, .after = Clusters) %>%
    #'  Drop sites with NAs (missing 1+ years of data)
    na.omit(.) %>%
    mutate(Clusters = as.character(Clusters),
           Year.x = as.character(Year.x),
           Year.y = as.character(Year.y)) 
  
  #'  Z-transform local abundance estimates (not actually needed - piecewiseSEM takes care of this)
  localN_z_1YrLag <- density_wide_1YrLag_20s_22s %>%
    mutate(across(where(is.numeric), ~(.x - mean(.x, na.rm = TRUE))/sd(.x, na.rm = TRUE)))
    
  #'  Create correlation matrix for all continuous covariates at once
  cov_correlation <- function(dat) {
    covs <- dat %>%
      dplyr::select(contains(c(".Tminus1", ".T")))
    cor_matrix <- cor(covs, use = "complete.obs")
    return(cor_matrix)
  }
  cov_correlation(localN_z_1YrLag) 
  
  #'  Visualize data
  plot_histograms <- function(dat) {
    ndat <- dplyr::select(dat, contains(c(".Tminus1", ".T"))) %>% as.matrix(.)
    for(i in 1:ncol(ndat)) {
      hist(ndat[,i])
    }
  }
  plot_histograms(localN_z_1YrLag)
  
  #'  Save outputs
  save(density_wide_1YrLag_20s_22s, file = "./Outputs/data_for_SEM_1YrLag_n.RData")
  
  #'  Reminder to ignore warnings if sourcing this script with 4_Structural_Equation_Models.R
  print("IGNORE WARNINGS!")
  
  #'  Next step:
  #'  Run structural equation models (script 4_Structural_Equation_Models.R)
  
  
