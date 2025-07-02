  #'  ----------------------------------
  #'  5. Result tables for RN model and relative density indices
  #'  
  #'  R code from Bassing, S.B., D.E. Ausband, M.A. Mumma, J. Baumgardt, S. Thompson, 
  #'  M.A. Hurley, and M. Falcy. 2025. "Disentangling the web of species interactions 
  #'  in a multi-predator, multi-prey community"
  #'  ----------------------------------
  #'  Script to summarize results from species and year specific Royle-Nichols
  #'  abundance models and derived relative density indices.
  #'  
  #'  Requires:
  #'    1. Saved RN model outputs generated with 1_Relative_abundance_Royle-Nichols_model.R
  #'    2. Saved relative density indices (RDI) ouptuts generated with 3_Relative_density_indices.R
  #'  ----------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(stringr)
  library(tidyverse)
    
  #'  Load model outputs (make sure only one run per model is in this folder)
  filenames <- list.files("./Outputs/JAGS_out", pattern="*.RData", full.names=TRUE)
  lapply(filenames, load, environment())
  
  #'  List models based on species and year
  yr_list <- list("2020", "2021", "2022")
  bear_list <- list(RN_bear_2020, RN_bear_2021, RN_bear_2022)
  coy_list <- list(RN_coy_2020, RN_coy_2021, RN_coy_2022)
  lion_list <- list(RN_lion_2020, RN_lion_2021, RN_lion_2022)
  wolf_list <- list(RN_wolf_2020, RN_wolf_2021, RN_wolf_2022)
  elk_list <- list(RN_elk_2020, RN_elk_2021, RN_elk_2022)
  moose_list <- list(RN_moose_2020, RN_moose_2021, RN_moose_2022)
  wtd_list <- list(RN_wtd_2020, RN_wtd_2021, RN_wtd_2022)
  
  #'  Define number of significant digits to round to
  rounddig <- 3
  
  #'  ----------------------------------------------
  ####  Average parameter estimates in wide format  ####
  #'  ----------------------------------------------
  #'  Model outputs (mean lambda, p, and r) per species
  params <- function(mod, spp, yr) {
    #'  Grab parameter estimates and SD per site
    lambda <- round(mod$mean$mu.lambda, 2)
    lambda.ll <- round(mod$q2.5$mu.lambda, rounddig)
    lambda.ul <- round(mod$q97.5$mu.lambda, rounddig)
    r <- round(mod$mean$mu.r, rounddig)
    r.ll <- round(mod$q2.5$mu.r, rounddig)
    r.ul <- round(mod$q97.5$mu.r, rounddig)
    p <- round(mod$mean$mean.p, rounddig)
    p.ll <- round(mod$q2.5$mean.p, rounddig)
    p.ul <- round(mod$q97.5$mean.p, rounddig)
    psi <- round(mod$mean$mean.psi, rounddig)
    psi.ll <- round(mod$q2.5$mean.psi, rounddig)
    psi.ul <- round(mod$q97.5$mean.psi, rounddig)
    #'  Merge and format into single data frame
    out <- cbind(lambda, lambda.ll, lambda.ul, r, r.ll, r.ul, p, p.ll, p.ul, psi, psi.ll, psi.ul)
    param_est <- as.data.frame(out) %>%
      mutate(Species = spp,
             Year = yr,
             # lambda = as.numeric(lambda),
             lambda.cri = paste0(lambda, " (", lambda.ll, " - ", lambda.ul, ")"),
             # r = as.numeric(r),
             r.cri = paste0(r, " (", r.ll, " - ", r.ul, ")"),
             # p = as.numeric(p),
             p.cri = paste0(p, " (", p.ll, " - ", p.ul, ")"),
             # psi = as.numeric(psi),
             psi.cri = paste0(psi, " (", psi.ll, " - ", psi.ul, ")")) %>%
      dplyr::select(c(Species, Year, lambda.cri, r.cri, p.cri, psi.cri))
    return(param_est)
  }
  study_yr <- list("2020", "2021", "2022")
  bear_mean_params <- mapply(params, bear_list, spp = "Black bear", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  coy_mean_params <- mapply(params, coy_list, spp = "Coyote", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  lion_mean_params <- mapply(params, lion_list, spp = "Mountain lion", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  wolf_mean_params <- mapply(params, wolf_list, spp = "Wolf", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  elk_mean_params <- mapply(params, elk_list, spp = "Elk", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  moose_mean_params <- mapply(params, moose_list, spp = "Moose", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  wtd_mean_params <- mapply(params, wtd_list, spp = "White-tailed deer", yr = study_yr, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  mean_params <- bind_rows(bear_mean_params, coy_mean_params, lion_mean_params, wolf_mean_params, elk_mean_params, moose_mean_params, wtd_mean_params)
  write_csv(mean_params, "./Outputs/RN_mean_lambda_table.csv")
  
  #'  --------------------------------------------
  ####  All coefficient estimates in long format  ####
  #'  --------------------------------------------
  mod_out_summary_table <- function(mod, spp, yr) {
    #'  Retain & reformat model coefficients and mean derived parameters (r, lambda, psi)
    out <- as.data.frame(mod$summary) %>%
      rownames_to_column(., "Parameter") %>%
      transmute(
        Species = spp,
        Year = yr,
        Parameter = Parameter,
        Mean = round(mean, rounddig),
        Lower_CRI = round(`2.5%`, rounddig),
        Upper_CRI = round(`97.5%`, rounddig),
        Mean = as.numeric(Mean),
        Lower_CRI = as.numeric(Lower_CRI),
        Upper_CRI = as.numeric(Upper_CRI),
        `95% CRI` = paste0(" ", Lower_CRI, " - ", Upper_CRI),
        overlap0 = ifelse(overlap0 == 0, FALSE, TRUE)) %>%
      dplyr::select(-c(Lower_CRI, Upper_CRI)) %>%
      #'  Give parameters more meaningful names
      mutate(
        Parameter = ifelse(Parameter == "beta0", "Intercept [GMU10A]", Parameter),
        Parameter = ifelse(Parameter == "beta1[2]", "GMU6", Parameter),              #beta1[1] set to 0 in JAGS
        Parameter = ifelse(Parameter == "beta1[3]", "GMU1", Parameter),
        Parameter = ifelse(Parameter == "alpha0", "Intercept [Random]", Parameter),
        Parameter = ifelse(Parameter == "alpha1[2]", "Setup [Road]", Parameter),     #alpha1[1] set to 0 in JAGS
        Parameter = ifelse(Parameter == "mu.lambda", "mean lambda", Parameter),
        Parameter = ifelse(Parameter == "mu.r", "mean Pr(ind. detection)", Parameter),
        Parameter = ifelse(Parameter == "mean.p", "mean Pr(detection)", Parameter),
        Parameter = ifelse(Parameter == "mean.psi", "mean Pr(occupancy)", Parameter),
        overlap0 = ifelse(str_detect(Parameter, "Intercept "), NA, overlap0))  %>%
      #'  Save only coefficients and mean derived parameters
      filter(Parameter == "Intercept [GMU10A]" | Parameter == "GMU6" | Parameter == "GMU1" | 
               Parameter == "Intercept [Random]" | Parameter == "Setup [Road]" | 
               Parameter == "mean lambda" | Parameter == "mean Pr(ind. detection)" | 
               Parameter == "mean Pr(detection)" | Parameter == "mean Pr(occupancy)")
    
    print(out)
    return(out)
  }
  bear_tbl <- mapply(mod_out_summary_table, mod = bear_list, spp = "Black bear", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  coy_tbl <- mapply(mod_out_summary_table, mod = coy_list, spp = "Coyote", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  lion_tbl <- mapply(mod_out_summary_table, mod = lion_list, spp = "Mountain lion", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  wolf_tbl <- mapply(mod_out_summary_table, mod = wolf_list, spp = "Wolf", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  elk_tbl <- mapply(mod_out_summary_table, mod = elk_list, spp = "Elk", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  moose_tbl <- mapply(mod_out_summary_table, mod = moose_list, spp = "Moose", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  wtd_tbl <- mapply(mod_out_summary_table, mod = wtd_list, spp = "White-tailed deer", yr = yr_list, SIMPLIFY = FALSE) %>% do.call("rbind", .)
  
  RN_result_tbl <- rbind(bear_tbl, coy_tbl, lion_tbl, wolf_tbl, elk_tbl, moose_tbl, wtd_tbl) 
  RN_result_tbl_list <- list(bear_tbl, coy_tbl, lion_tbl, wolf_tbl, elk_tbl, moose_tbl, wtd_tbl) 
  
  #'  Save
  write_csv(RN_result_tbl, file = "./Outputs/RN_result_tbl.csv")
  save(RN_result_tbl_list, file = "./Outputs/RN_result_tbl_list.R")
  
  #'  -------------------
  ####  RDI per cluster  ####
  #'  -------------------
  #'  Relative Density Index per cluster (loads as cluster_density)
  load("./Outputs/RelativeDensityIndex_per_SppCluster.RData")
  head(cluster_density[[1]]); head(cluster_density[[2]]); head(cluster_density[[3]])
  rdi <- bind_rows(cluster_density[[1]], cluster_density[[2]], cluster_density[[3]]) %>%
    as.data.frame(.) %>%
    dplyr::select(c(GMU, Year, Clusters, area_km2, Species, SppDensity.100km2.r)) %>%
    filter(!is.na(Clusters))
  
  #'  Quick visualization of spread of RDI per species across years & GMUs
  ggplot(rdi, aes(Species, SppDensity.100km2.r, fill=factor(Species))) +
    geom_boxplot()
  
  #'  RDI per cluster, year, and species in wide-format for publication
  rdi_all <- rdi %>%
    mutate(SppDensity.100km2.r = round(SppDensity.100km2.r, 2),
           area_km2 = round(area_km2, 2),
           Species = ifelse(Species == "bear_black", "Black bear RDI", Species),
           Species = ifelse(Species == "coyote", "Coyote RDI", Species),
           Species = ifelse(Species == "elk", "Elk RDI", Species),
           Species = ifelse(Species == "moose", "Moose RDI", Species),
           Species = ifelse(Species == "mountain_lion", "Mountain lion RDI", Species),
           Species = ifelse(Species == "whitetailed_deer", "White-tailed deer RDI", Species),
           Species = ifelse(Species == "wolf", "Wolf RDI", Species),
           GMU = factor(GMU, levels = c("GMU1", "GMU6", "GMU10A"))) %>%
    pivot_wider(names_from = Species, values_from = SppDensity.100km2.r) %>%
    rename("Area (km2)" = "area_km2") %>%
    rename("Cluster" = "Clusters") %>%
    arrange(GMU, Year)
  
  write_csv(rdi_all, "./Outputs/Cluster_Relative_Density_all_Spp.csv")
  
  #'  Average RDI per species per GMU, across years
  mean_rdi_per_GMU <- rdi %>%
    dplyr::select(c(Species, GMU, SppDensity.100km2.r)) %>%
    group_by(Species, GMU) %>%
    summarise(mean_RDI = round(mean(SppDensity.100km2.r), 2),
              se_RDI = round(sd(SppDensity.100km2.r)/sqrt(length(SppDensity.100km2.r)), 2)) %>%
    ungroup() %>%
    # arrange(Species, -mean_RDI, GMU) %>%
    mutate(Species = ifelse(Species == "bear_black", "Black bear", Species),
           Species = ifelse(Species == "coyote", "Coyote", Species),
           Species = ifelse(Species == "elk", "Elk", Species),
           Species = ifelse(Species == "moose", "Moose", Species),
           Species = ifelse(Species == "mountain_lion", "Mountain lion", Species),
           Species = ifelse(Species == "whitetailed_deer", "White-tailed deer", Species),
           Species = ifelse(Species == "wolf", "Wolf", Species),
           mean_RDI = paste0(mean_RDI, " (", se_RDI, ")"),
           GMU = factor(GMU, levels = c("GMU1", "GMU6", "GMU10A"))) %>%
    dplyr::select(-se_RDI) %>%
    rename("Mean RDI (SE)" = "mean_RDI") 

  write_csv(mean_rdi_per_GMU, "./Outputs/GMU_Average_Relative_Density_all_Spp.csv")
  
  #'  Average RDI per species, across years and GMUs
  (mean_rdi <- rdi %>%
    dplyr::select(c(Species, SppDensity.100km2.r)) %>%
    group_by(Species) %>%
    summarise(mean_RDI = round(mean(SppDensity.100km2.r), 2),
              se_RDI = round(sd(SppDensity.100km2.r)/sqrt(length(SppDensity.100km2.r)), 2)) %>%
    ungroup() %>%
    arrange(-mean_RDI))
  
  #'  Next step:
  #'  Generate figures for publication (script 6_Result_Figures.R)
  
  