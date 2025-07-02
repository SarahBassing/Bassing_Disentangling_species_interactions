  #'  -------------------------
  #'  6. Figures for publication
  #'  
  #'  R code from Bassing, S.B., D.E. Ausband, M.A. Mumma, J. Baumgardt, S. Thompson, 
  #'  M.A. Hurley, and M. Falcy. 2025. "Disentangling the web of species interactions 
  #'  in a multi-predator, multi-prey community"
  #'  -------------------------
  #'  Script to generate maps used in publication.
  #'  
  #'  THIS SCRIPT IS PROVIDED FOR TRANSPARENCY BUT ONLY THE STUDY AREA MAP CAN BE
  #'  REPRODUCED WITH PUBLICLY AVAILABLE DATA. CAMERA LOCATIONS ARE SENSITIVE BUT 
  #'  AVAILABLE UPON REQUEST FROM IDAHO DEPT OF FISH & GAME FOR QUALIFIED RESEARCHERS.
  #'  
  #'  Requires:
  #'    1. Spatial data from camera clusters generated with 2_Cluster_cameras.R
  #'    2. Spatial data for Royle-Nichols relative abundance indices generated with
  #'    1_Relative_abundance_Royle-Nichols_model.R
  #'    3. Relative density indices per cluster generated with 3_Relative_density_indices.R
  #'    4. Spatial data for camera locations and water bodies saved in Data folder
  #'  -------------------------
  
  #'  Load libraries
  library(sf)
  library(terra)
  library(ggplot2)
  library(cowplot)
  library(ggh4x)
  library(patchwork)
  library(tidyverse)
  
  #'  Load spatial data that is publicly available or generated with previous scripts
  #'  Cluster polygons
  gmu1_poly <- st_read("./Outputs/Spatial_outputs/cluster_polygons_gmu1.shp") %>% mutate(GMU = "GMU1")
  gmu6_poly <- st_read("./Outputs/Spatial_outputs/cluster_polygons_gmu6.shp") %>% mutate(GMU = "GMU6")
  gmu10a_poly <- st_read("./Outputs/Spatial_outputs/cluster_polygons_gmu10a.shp") %>% mutate(GMU = "GMU10A")
  
  #'  Relative Density Index per cluster (loads as cluster_density)
  load("./Outputs/RelativeDensityIndex_per_SppCluster.RData")
  head(cluster_density[[1]]); head(cluster_density[[2]]); head(cluster_density[[3]])
  
  #'  Idaho GMUs and water
  usa_wgs84 <- st_read("./Data/Spatial_data/US_states.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  Idaho <- st_read("./Data/Spatial_data/IdahoState.shp")
  all_gmu_wgs84 <- st_read("./Data/Spatial_data/Game_Management_Units.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  gmu_wgs84 <- st_read("./Data/Spatial_data/Study_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  bigwater_wgs84 <- st_read("./Data/Spatial_data/Idaho_waterbodies_1km2.shp") %>%
    st_transform(crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
    filter(gnis_name == "Priest Lake" | gnis_name == "Upper Priest Lake" | 
             gnis_name == "Lake Pend Oreille" | gnis_name == "Chatcolet Lake" | 
             gnis_name == "Dworshak Reservoir") %>%
    dplyr::select(gnis_name)
  priestrivers <- st_read("./Data/Spatial_data/PriestRiver_flowline.shp")
  pendoreille <- st_read("./Data/Spatial_data/Pendoreille_flowline.shp")
  kootenairiver <- st_read("./Data/Spatial_data/KootenaiRiver_flowline.shp")
  clarkfork <- st_read("./Data/Spatial_data/ClarkForkRiver_flowline.shp") 
  clearwater <- st_read("./Data/Spatial_data/NorthForkClearwater_flowline.shp")
  rivers <- bind_rows(priestrivers, pendoreille, kootenairiver, clarkfork, clearwater) 
  rivers_clip <- st_intersection(rivers, gmu_wgs84)
  
  #'  ------------------
  ####  Study area map  ####
  #'  ------------------
  #####  USA with Idaho highlighted  #####
  usa_map <- ggplot() +
    geom_sf(data = usa_wgs84, fill = "white", color = "gray50", linewidth = 0.2) +
    geom_sf(data = Idaho, fill = "gray20", color = "gray20") +
    #'  Constrain plot to the lower 48 
    coord_sf(xlim = c(-124.2, -68.2), ylim = c(25, 49)) +
    theme_minimal() +  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          axis.text.x=element_blank(), axis.ticks.x=element_blank(), axis.title.x=element_blank(),
          axis.text.y=element_blank(), axis.ticks.y=element_blank(), axis.title.y=element_blank(),
          #'  No margins around figure
          plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt"))
  
  #####  Map #2  #####
  #'  Study area colors
  gmu_wgs84 <- mutate(gmu_wgs84, NAME = as.factor(NAME),
                      NAME = factor(NAME, levels = c("1", "6", "10A")))
  gmu_color <- c("darkcyan", "lightcoral", "darkgoldenrod3")
  
  #'  Plot study areas within Idaho
  ID_study_areas <- ggplot() +
    geom_sf(data = all_gmu_wgs84, fill = "gray95", color="gray50", linewidth = 0.2) + 
    geom_sf(data = gmu_wgs84, aes(fill = NAME, color = NAME), linewidth = 0.2, alpha = 0.7) + 
    scale_fill_manual(values = gmu_color) +
    scale_color_manual(values = gmu_color) +
    labs(fill = "GMU", color = "GMU") +
    geom_rect(aes(xmin = -117.3, xmax = -114.15, ymin = 45.75, ymax = 49.25), color = "black", fill = NA, linewidth = 0.2) + 
    theme_minimal() +
    theme(axis.text = element_blank(),
          axis.ticks = element_blank(),
          panel.grid.minor = element_blank(),
          panel.grid.major = element_blank(),
          legend.title = element_text(size = 12), 
          legend.text = element_text(size = 12), 
          legend.key.size = unit(0.5, "lines")) +
    theme(legend.justification = c(1, 0)) 
  
  #'  Draw map of study area with inset of USA & Idaho
  StudyArea_Map <- ggdraw(ID_study_areas + 
                            theme(legend.box.margin = margin(0,-8,0,-10))) + 
    draw_plot(
      {
        usa_map
      },
      #'  Distance along a (0,1) x-axis to draw the left edge of the plot
      x = 0.40, 
      #'  Distance along a (0,1) y-axis to draw the bottom edge of the plot
      y = 0.50, 
      #'  Width & height of the plot expressed as proportion of the entire ggdraw object
      #'  THIS DEPENDS ON HOW BIG YOUR PLOT WINDOW IS TOO!!!!
      width = 0.55, 
      height = 0.55) +
    labs(subtitle = "a") + 
    theme(plot.subtitle = element_text(hjust = 0))
  
  plot(StudyArea_Map) 
  #'  NOTE: placement of inset and margins will vary with plot window. Values 
  #'  for xy,wh are based on exported figures for publication
  
  
  #'  --------------------------------------------
  ####  Map Relative Abundance & Density Indices  ####
  #'  --------------------------------------------
  #'  THE FOLLOWING CODE IS PROVIDED FOR TRANSPARENCY BUT CANNOT BE RUN WITH 
  #'  PUBLICLY AVAILABLE DATA. CAMERA LOCATIONS ARE SENSITIVE BUT AVAILABLE
  #'  UPON REQUEST FROM IDAHO DEPARTMENT OF FISH & GAME FOR QUALIFIED RESEARCHERS.
  #'  Camera clusters
  
  #'  Load camera cluster data and spatial data for RAI
  clusters_gmu1 <- st_read("./Outputs/Spatial_outputs/cam_clusters_gmu1.shp") 
  clusters_gmu6 <- st_read("./Outputs/Spatial_outputs/cam_clusters_gmu6.shp") 
  clusters_gmu10a <- st_read("./Outputs/Spatial_outputs/cam_clusters_gmu10a.shp") 
  cameras_per_cluster <- bind_rows(clusters_gmu1, clusters_gmu6, clusters_gmu10a) %>%
    group_by(CamID) %>%
    slice(1L) %>%
    ungroup()
  
  #'  Royle-Nichols estimates spatial data
  load("./Outputs/Spatial_outputs/spatial_RN_list.RData")
  
    #'  Camera locations
  cams_20s_wgs84 <- st_read("./Data/Spatial_data/cams_2020_wgs84.shp")
  cams_21s_wgs84 <- st_read("./Data/Spatial_data/cams_2021_wgs84.shp")
  cams_22s_wgs84 <- st_read("./Data/Spatial_data/cams_2022_wgs84.shp")
  cam_list <- list(cams_20s_wgs84, cams_21s_wgs84, cams_22s_wgs84)
  
  
  #'  -----------------------
  ####  Format spatial data  ####
  #'  -----------------------
  #'  Define study years
  year_list <- list("2020", "2021", "2022")
  
  #'  Add year to each RN abundance dataframe and unlist per species
  add_yr <- function(dat, yr) {
    dat$Year <- yr
    return(dat)
  }
  rn_bear_all <- mapply(add_yr, dat = spatial_RN_list[[1]], yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_coy_all <- mapply(add_yr, dat = spatial_RN_list[[3]], yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_lion_all <- mapply(add_yr, dat = spatial_RN_list[[4]], yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_wolf_all <- mapply(add_yr, dat = spatial_RN_list[[5]], yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_elk_all <- mapply(add_yr, dat = spatial_RN_list[[6]], yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_moose_all <- mapply(add_yr, dat = spatial_RN_list[[8]], yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_wtd_all <- mapply(add_yr, dat = spatial_RN_list[[9]], yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  #'  Merge cluster polygons across GMUs
  cluster_poly <- bind_rows(gmu1_poly, gmu6_poly, gmu10a_poly) %>%
    rename(ClusterID = Clusters)
  
  #'  Merge spatial camera cluster data together and reduce to one observation per camera
  clusters_all <- bind_rows(clusters_gmu1, clusters_gmu6, clusters_gmu10a) %>%
    dplyr::select(c(CamID, CellID, GMU, Setup, Clustrs, are_km2, geometry)) %>%
    rename(ClusterID = Clustrs) %>%
    rename(area_km2 = are_km2) %>%
    group_by(CamID) %>%
    slice(1L) %>%
    ungroup()
    
  #'  Merge RDI into single document
  rdi <- bind_rows(cluster_density) %>%
    filter(!is.na(ClusterID))
    
  
  #'  ------------------------------------------
  ####  Plot Royle-Nichols abundance estimates  ####
  #'  ------------------------------------------
  #'  Make one giant faceted plot where rows represent GMU and columns represent years 
  #'  to ensure that the dot sizes are all consistent for at least a single species
  map_RoyleNichols <- function(sf_rn, spp) {
    #'  Define size of circles
    size_breaks <- c(0, 1, 2, 3, 5, 7, 9, 12)
    
    sf_rn <- mutate(sf_rn, GMU = factor(GMU, levels = c("GMU1", "GMU6", "GMU10A")))
    pal <- c("darkcyan", "lightcoral", "darkgoldenrod3")
    
    #'  Create figure
    spp_rn <- ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) + 
      geom_sf(data = rivers_clip, color = "lightskyblue3") + 
      geom_sf(data = bigwater_wgs84, fill = "lightskyblue2") +
      geom_sf(data = sf_rn, aes(size = RN.n.rounded, colour = GMU, fill = GMU), shape = 21, alpha = 3/10) +
      scale_size_continuous(breaks = size_breaks, range = c(0,12)) +
      scale_color_manual(values = pal) +
      scale_fill_manual(values = pal) +
      labs(size = "Predicted \nRAI", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 18)) +
      labs(title = paste("Predicted relative abundance indicies for", spp)) + 
      # facet_grid(rows = vars(Species), cols = vars(Year)) +
      # force_panelsizes(rows = 2, cols = 1, TRUE)
      facet_wrap(~Year) 
    
    #'  Plot each map
    print(spp_rn)
    
    return(spp_rn)
  }
  rn_maps_bear <- map_RoyleNichols(rn_bear_all, spp = "black bears")
  rn_maps_coy <- map_RoyleNichols(rn_coy_all, spp = "coyotes")
  rn_maps_lion <- map_RoyleNichols(rn_lion_all, spp = "mountain lions")
  rn_maps_wolf <- map_RoyleNichols(rn_wolf_all, spp = "wolves") 
  rn_maps_elk <- map_RoyleNichols(rn_elk_all, spp = "elk")
  rn_maps_moose <- map_RoyleNichols(rn_moose_all, spp = "moose")
  rn_maps_wtd <- map_RoyleNichols(rn_wtd_all, spp = "white-tailed deer")
  
  
  (rn_maps_wolf_refined <- rn_maps_wolf + 
      scale_size_continuous(breaks = c(0,1, 2, 5, 7), range = c(0,7)) +
      guides(colour = "none", fill = "none") + ggtitle(NULL) + 
      labs(subtitle = "b") +
      theme(legend.title = element_text(size = 12), 
            legend.text = element_text(size = 12), 
            legend.key.size = unit(0.5, "lines"),
            legend.justification = c(1, 0),
            plot.subtitle = element_text(vjust = 0),
            #'  No margins around figure
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")))

  ggsave("./Outputs/RN_map_blackbear.tiff", rn_maps_bear,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RN_map_coyote.tiff", rn_maps_coy,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RN_map_lion.tiff", rn_maps_lion,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RN_map_wolf.tiff", rn_maps_wolf,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RN_map_elk.tiff", rn_maps_elk,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RN_map_moose.tiff", rn_maps_moose,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RN_map_wtd.tiff", rn_maps_wtd,
         units = "in", width = 13, height = 11, dpi = 400, device = "tiff", compression = "lzw")
  
  
  #'  -------------------------------------------
  ####  Plot Relative Density Index per Cluster  ####
  #'  -------------------------------------------
  map_rdi <- function(sf_rdi, q) {
    
    #'  Define 5%, 25%, 50%, 75%, and 95% quantiles of RDI
    q0 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0), 0))
    q5 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.05), 0))
    q25  <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.25), 0))
    q50 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.5), 0))
    q75 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.75), 0))
    q95 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.95), 0))
    q100 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 1), 0))
    
    #'  Define breaks in continuous fill scale (if using this version)
    size_breaks <- c(q0, q5, q25, q50, q75, q95, q100)
    print(size_breaks)
    
    #'  Define spatial extent of study areas
    bbox <- st_bbox(gmu_wgs84)
    
    #'  Create figure
    spp_rdi <- ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) +
      geom_sf(data = sf_rdi, aes(fill = SppDensity.100km2.r, group = Year)) +  
      scale_fill_gradient2(low = "#f5f5f5", mid = "#54a1a1", high = "#1f6f6f", midpoint = q) +
      geom_sf(data = rivers_clip, color = "lightskyblue3") +
      geom_sf(data = bigwater_wgs84, fill = "lightskyblue2") +
      coord_sf(xlim = c(bbox[1], bbox[3]), ylim = c(bbox[2], bbox[4])) +
      labs(fill = "Relative \nDensity \nIndex (RDI)", x = "Longitude", y = "Latitude") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            text = element_text(size = 18)) +
      facet_wrap(~Year)
    
    return(spp_rdi)
  }
  (rdi_map_bear_2020 <- map_rdi(rdi[rdi$Species == "bear_black" & rdi$Year == "2020",], q = 9) + theme(axis.title.x = element_blank())) 
  (rdi_map_bear_2021 <- map_rdi(rdi[rdi$Species == "bear_black" & rdi$Year == "2021",], q = 8) + 
      guides(fill="none") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_bear_2022 <- map_rdi(rdi[rdi$Species == "bear_black" & rdi$Year == "2022",], q = 9) + 
      guides(fill="none") + theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_bear <- rdi_map_bear_2020 + rdi_map_bear_2021 + rdi_map_bear_2022 + plot_layout(ncol = 3, guides = 'collect') + 
    plot_annotation(title = "Relative density indices for black bears", theme = theme(plot.title = element_text(size = 18))))
  
  (rdi_map_coy_2020 <- map_rdi(rdi[rdi$Species == "coyote" & rdi$Year == "2020",], q = 11) + guides(fill="none") + theme(axis.title.x = element_blank()))
  (rdi_map_coy_2021 <- map_rdi(rdi[rdi$Species == "coyote" & rdi$Year == "2021",], q = 11) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_coy_2022 <- map_rdi(rdi[rdi$Species == "coyote" & rdi$Year == "2022",], q = 11) + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_coy <- rdi_map_coy_2020 + rdi_map_coy_2021 + rdi_map_coy_2022 + plot_layout(ncol = 3, guides = 'collect') + 
    plot_annotation(title = "Relative density indices for coyotes", theme = theme(plot.title = element_text(size = 18))))
  
  (rdi_map_lion_2020 <- map_rdi(rdi[rdi$Species == "mountain_lion" & rdi$Year == "2020",], q = 7) + theme(axis.title.x = element_blank()))
  (rdi_map_lion_2021 <- map_rdi(rdi[rdi$Species == "mountain_lion" & rdi$Year == "2021",], q = 7) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_lion_2022 <- map_rdi(rdi[rdi$Species == "mountain_lion" & rdi$Year == "2022",], q = 7) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_lion <- rdi_map_lion_2020 + rdi_map_lion_2021 + rdi_map_lion_2022 + plot_layout(ncol = 3, guides = 'collect') + 
    plot_annotation(title = "Relative density indices for mountain lions", theme = theme(plot.title = element_text(size = 18))))
  
  (rdi_map_elk_2020 <- map_rdi(rdi[rdi$Species == "elk" & rdi$Year == "2020",], q = 16) + theme(axis.title.x = element_blank()))
  (rdi_map_elk_2021 <- map_rdi(rdi[rdi$Species == "elk" & rdi$Year == "2021",], q = 16) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_elk_2022 <- map_rdi(rdi[rdi$Species == "elk" & rdi$Year == "2022",], q = 16) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_elk <- rdi_map_elk_2020 + rdi_map_elk_2021 + rdi_map_elk_2022 + plot_layout(ncol = 3, guides = 'collect') + 
      plot_annotation(title = "Relative density indices for elk", theme = theme(plot.title = element_text(size = 18))))
  
  (rdi_map_moose_2020 <- map_rdi(rdi[rdi$Species == "moose" & rdi$Year == "2020",], q = 7) + guides(fill="none") + 
      theme(axis.title.x = element_blank()))
  (rdi_map_moose_2021 <- map_rdi(rdi[rdi$Species == "moose" & rdi$Year == "2021",], q = 7) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_moose_2022 <- map_rdi(rdi[rdi$Species == "moose" & rdi$Year == "2022",], q = 7) + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_moose <- rdi_map_moose_2020 + rdi_map_moose_2021 + rdi_map_moose_2022 + plot_layout(ncol = 3, guides = 'collect') + 
      plot_annotation(title = "Relative density indices for moose", theme = theme(plot.title = element_text(size = 18))))
  
  (rdi_map_wtd_2020 <- map_rdi(rdi[rdi$Species == "whitetailed_deer" & rdi$Year == "2020",], q = 20) + guides(fill="none") +
      theme(axis.title.x = element_blank()))
  (rdi_map_wtd_2021 <- map_rdi(rdi[rdi$Species == "whitetailed_deer" & rdi$Year == "2021",], q = 20) + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank()))
  (rdi_map_wtd_2022 <- map_rdi(rdi[rdi$Species == "whitetailed_deer" & rdi$Year == "2022",], q = 20) + guides(fill="none") + 
      theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.title.y = element_blank(), axis.title.x = element_blank()))
  (rdi_map_wtd <- rdi_map_wtd_2020 + rdi_map_wtd_2021 + rdi_map_wtd_2022 + plot_layout(ncol = 3, guides = 'collect') + 
      plot_annotation(title = "Relative density indices for white-tailed deer", theme = theme(plot.title = element_text(size = 18))))
  
  #'  Plot wolf RDI across years to help with color gradient
  #'  Scale of RDI is so different in 2021 that it makes things hard to compare if
  #'  each year is plotted separately
  wolf_rdi <- function(sf_rdi, q) {
    
    #'  Define 5%, 25%, 50%, 75%, and 95% quantiles of RDI
    q0 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0), 0))
    q5 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.05), 0))
    q25  <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.25), 0))
    q50 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.5), 0))
    q75 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.75), 0))
    q95 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 0.95), 0))
    q100 <- as.numeric(round(quantile(sf_rdi$SppDensity.100km2.r, 1), 0))
    
    #'  Define breaks in continuous fill scale (if using this version)
    size_breaks <- c(q0, q5, q25, q50, q75, q95, q100)
    print(size_breaks)
    
    #'  Define spatial extent of study areas
    bbox <- st_bbox(gmu_wgs84)
    
    #'  Create figure
    spp_rdi <- ggplot() +
      geom_sf(data = gmu_wgs84, fill = NA) +
      geom_sf(data = sf_rdi, aes(fill = SppDensity.100km2.r, group = Year)) + 
      scale_fill_gradient2(low = "#f5f5f5", mid = "#54a1a1", high = "#1f6f6f", midpoint = q) +
      geom_sf(data = rivers_clip, color = "lightskyblue3") +
      geom_sf(data = bigwater_wgs84, fill = "lightskyblue2") +
      labs(fill = "Relative \nDensity \nIndex (RDI)", x = "Longitude", y = "Latitude",
           subtitle = "c") +
      theme_classic() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
            legend.title = element_text(size = 12), 
            legend.text = element_text(size = 12), 
            legend.key.size = unit(0.5, "lines"),
            legend.justification = c(1, 0),
            plot.subtitle = element_text(vjust = 0),
            #'  No margins around figure
            plot.margin = margin(t = 0, r = 0, b = 0, l = 0, unit = "pt")) +
      facet_wrap(~Year)
    
    return(spp_rdi)
  }
  (rdi_map_wolf <- wolf_rdi(rdi[rdi$Species == "wolf",], q = 4) +
      theme(axis.title.y = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank())) 
    
  #'  Save
  ggsave("./Outputs/RDI_bear.tiff", rdi_map_bear,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RDI_coyote.tiff", rdi_map_coy,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RDI_lion.tiff", rdi_map_lion,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RDI_wolf_with_cams.tiff", rdi_map_wolf,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RDI_elk.tiff", rdi_map_elk,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RDI_moose.tiff", rdi_map_moose,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  ggsave("./Outputs/RDI_wtd.tiff", rdi_map_wtd,
         units = "in", width = 12, height = 10, dpi = 400, device = "tiff", compression = "lzw")
  
  
  #'  -------------------------------------------------------------
  #### Figure mapping study area, Wolf RN estimates, and Wolf RDI  ####
  #'  -------------------------------------------------------------
  #'  Define lay out of three plots in relation to each other
  layout <- c(
    area(t = 1, b = 3, l = 1, r = 3),
    area(t = 3, b = 12, l = 3.5, r = 7),
    area(t = 3, b = 12, l = 8.5, r = 12)
  )
  
  #'  Create paneled figure
  Full_figure <- StudyArea_Map + rn_maps_wolf_refined + rdi_map_wolf + 
    plot_layout(ncol = 3, guides = 'collect') + 
    plot_layout(widths = c(0.85, 1.75, 1.75)) +
    plot_annotation(title = "Predicted wolf relative abundance and density in each study area") & 
    theme(legend.justification = c(1, 0),
          plot.title = element_text(size = 16, vjust = -5),
          text = element_text(size = 14),
          plot.tag = element_text(size = 14),
          plot.margin = unit(c(t = 1, r = 0, b = 0, l = 0), "cm")) & 
    labs(title = NULL)
  Full_figure
  #'  Note: proportions will look terrible in plot viewer - this is formatted for
  #'  exported figures
  
  
  #'  Save
  ggsave("./Outputs/StudyArea_WolfResults.tiff", Full_figure, 
         units = "cm", width = 36, height = 17, dpi = 600, device = 'tiff', compression = 'lzw')
  
  
  #'  Fin!
  