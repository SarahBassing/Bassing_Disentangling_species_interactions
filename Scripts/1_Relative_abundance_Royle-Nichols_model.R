  #'  --------------------------------------
  #'  1. Relative abundance index estimation
  #'  
  #'  R code from Bassing, S.B., D.E. Ausband, M.A. Mumma, J. Baumgardt, S. Thompson, 
  #'  M.A. Hurley, and M. Falcy. 2025. "Disentangling the web of species interactions 
  #'  in a multi-predator, multi-prey community"
  #'  --------------------------------------
  #'  Script loads detection histories for each species, site, and summer, then 
  #'  bundles data and runs Royle-Nichols abundance models in JAGS to estimate 
  #'  relative abundance index (RAI) at each camera site. Script sources
  #'  1a_RNmodel_JAGS_code_2020.R and 1b_RNmodel_JAGS_code_2021_2022.R to fit
  #'  models with JAGS. Using 2 different scripts because surveyed different 
  #'  number of study areas in 2020 vs 2021/2022, leading to different number of
  #'  parameters included in the models.
  #'  
  #'  Anonymized data formatted as detection histories are described below. 
  #'  Qualified researchers can contact Matt Boone (idfgdatarequests@idfg.idaho.gov), 
  #'  Data Management Lead, Idaho Department of Fish and Game, Idaho Fish and Wildlife 
  #'  Information System, 600 S Walnut, Boise, ID 83712 for full data requests. 
  #'  
  #'  3 study areas (GMU 1, 6, and 10A) were surveyed June 1 - Sept 15, 2020, 2021, and 2022.
  #'  In summer 2020, only GMU 6 & 10A were surveyed, with a total of 476 camera sites in operation. 
  #'  In summer 2021, all study areas were surveyed, with a total of 694 camera sites in operation. 
  #'  In summer 2022, all study areas were surveyed, with a total of 700 camera sites in operation. 
  #'  Script loads 2 data sets comprising 6 files:
  #'    1. DH_smr2020.RData, DH_smr2021.RData, DH_smr2022.RData: 
  #'        Each contains a list of 7 data frames. Each data frame contains 1 detection 
  #'        history per species: [[1]] = black bear; [[2]] = coyote; [[3]] = mountain lion;
  #'        [[4]] = wolf; [[5]] = elk; [[6]] = moose; [[7]] = white-tailed deer
  #'    2. stations_smr2020.RData, stations_smr2021.RData, stations_smr2022.RData: 
  #'        Each is a data frame where rows represent a unique camera station and 
  #'        columns indicate each station's unique ID, season of deployment, study  
  #'        area (GMU) of deployment, and the camera setup type. Also contains 
  #'        the cluster each camera is assigned to in a later stage and the area 
  #'        of the cluster. These data are provided for a later stage of the analysis
  #'        (3_Relative_density_indices.R) and can be ignored here.
  #'            GMU has 3 factor levels: 1 = GMU10A, 2 = GMU6, 3 = GMU1 
  #'            Setup has 2 factor levels: 1 = random placement, 2 = trail/road placement
  #'  --------------------------------------
  
  #'  Clean workspace
  rm(list = ls())
  
  #'  Load libraries
  library(camtrapR)
  library(chron)
  library(lubridate)
  library(unmarked)
  library(jagsUI)
  library(mcmcplots)
  library(tidyverse)
  
  #'  -------------
  ####  Load data  ####
  #'  -------------
  #'  Detection histories - each file contains a list of 7 detection histories, 1 per species for a given year.
  load("./Data/DH_smr2020.RData") # summer 2020 detection histories
  load("./Data/DH_smr2021.RData") # summer 2021 detection histories
  load("./Data/DH_smr2022.RData") # summer 2022 detection histories
  #'  Camera station data
  load("./Data/stations_smr2020.RData") # summer 2020 sites
  load("./Data/stations_smr2021.RData") # summer 2021 sites
  load("./Data/stations_smr2022.RData") # summer 2022 sites
  
  #'  ------------------------
  ####  Setup data for JAGS  ####
  #'  ------------------------
  #'  Bundle detection histories and covariates for each species and year
  bundle_dat <- function(dh, nsite, nsurvey, cov, effort) {
    #'  Convert detection history to matrix
    dh <- as.matrix(dh)
    dimnames(dh) <- NULL
    #'  Count number of sites per GMU
    ncams_perGMU <- cov %>%
      group_by(GMU) %>%
      summarise(nsites = n()) %>%
      ungroup()
    #'  Split up covariates by GMU
    covs_GMU10A <- filter(cov, GMU == 1)
    covs_GMU6 <- filter(cov, GMU == 2)
    covs_GMU1 <- filter(cov, GMU == 3)
    #'  Bundle data for JAGS
    bundled <- list(y = dh, 
                    nsites = dim(dh)[1], 
                    nsurveys = dim(dh)[2], 
                    ngmu = max(as.numeric(cov$GMU)),
                    nsets = max(as.numeric(cov$Setup)),
                    ncams1 = as.numeric(ncams_perGMU[1,2]), # GMU10A
                    ncams2 = as.numeric(ncams_perGMU[2,2]), # GMU6
                    ncams3 = as.numeric(ifelse(is.na(ncams_perGMU[3,2]), 0, ncams_perGMU[3,2])), #GMU1
                    gmu = as.numeric(cov$GMU), 
                    setup = as.numeric(cov$Setup),
                    #'  Area of each (km2)
                    area1 = as.numeric(8527.31),
                    area2 = as.numeric(5905.44),
                    area3 = as.numeric(14648.92))
    str(bundled)
    return(bundled)
  }
  data_bundle_2020 <- lapply(DH_smr2020, bundle_dat, cov = stations_smr2020)
  data_bundle_2021 <- lapply(DH_smr2021, bundle_dat, cov = stations_smr2021)
  data_bundle_2022 <- lapply(DH_smr2022, bundle_dat, cov = stations_smr2022)
  
  #'  Initial values
  #'  Using naive occupancy as a starting point for local abundance
  initial_n <- function(dh) {
    #'  Max value per row
    ninit <- apply(dh, 1, max, na.rm = TRUE)
    ninit <- as.vector(ninit)
    return(ninit)
  }
  #'  Apply function per species for each year
  ninit_2020 <- lapply(DH_smr2020, initial_n)
  ninit_2021 <- lapply(DH_smr2021, initial_n)
  ninit_2022 <- lapply(DH_smr2022, initial_n)
  
  #'  Parameters monitored
  params <- c("beta0", "beta1", "alpha0", "alpha1", "rSetup", "mu.r", "mean.p", 
              "mu.lambda", "totalN", "occSites", "mean.psi", 
              "totalN.gmu10a", "densitykm2.gmu10a", "density100km2.gmu10a", 
              "totalN.gmu6", "densitykm2.gmu6", "density100km2.gmu6", 
              "totalN.gmu1", "densitykm2.gmu1", "density100km2.gmu1", "N")
  #'  mu.lambda = lambda averaged across all GMUs
  #'  mu.r = per-individual detection probability averaged across all sites 
  
  #'  MCMC settings
  nc <- 3
  ni <- 50000
  nb <- 10000
  nt <- 10
  na <- 5000
  
  #'  --------------------------
  ####  Run RN model with JAGS  ####
  #'  --------------------------
  #####  2020  Analyses  #####
  source("./Scripts/1a_RNmodel_JAGS_code_2020.R")
  
  ######  Black bear  ######
  start.time = Sys.time()
  inits_bear2020 <- function(){list(N = ninit_2020[[1]])}
  RN_bear_2020 <- jags(data_bundle_2020[[1]], inits = inits_bear2020, params,
                      "./Outputs/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_2020$summary)
  which(RN_bear_2020$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_2020$samples)
  save(RN_bear_2020, file = "./Outputs/JAGS_out/RN_bear_2020.RData")

  ######  Coyote  ######
  start.time = Sys.time()
  inits_coy2020 <- function(){list(N = ninit_2020[[2]])}
  RN_coy_2020 <- jags(data_bundle_2020[[2]], inits = inits_coy2020, params,
                     "./Outputs/JAGS_RNmod_2020.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_2020$summary)
  which(RN_coy_2020$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_2020$samples)
  save(RN_coy_2020, file = "./Outputs/JAGS_out/RN_coy_2020.RData")
  
  ######  Mountain lion  ######
  start.time = Sys.time()
  inits_lion2020 <- function(){list(N = ninit_2020[[3]])}
  RN_lion_2020 <- jags(data_bundle_2020[[3]], inits = inits_lion2020, params,
                      "./Outputs/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_2020$summary)
  which(RN_lion_2020$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_2020$samples)
  save(RN_lion_2020, file = "./Outputs/JAGS_out/RN_lion_2020.RData") 
  
  ######  Wolf  ######
  ni_wolf <-  100000 
  start.time = Sys.time()
  inits_wolf2020 <- function(){list(N = ninit_2020[[4]])}
  RN_wolf_2020 <- jags(data_bundle_2020[[4]], inits = inits_wolf2020, params,
                      "./Outputs/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_2020$summary)
  which(RN_wolf_2020$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_2020$samples)
  save(RN_wolf_2020, file = "./Outputs/JAGS_out/RN_wolf_2020.RData") 
  
  ######  Elk  ######
  start.time = Sys.time()
  inits_elk2020 <- function(){list(N = ninit_2020[[5]])}
  RN_elk_2020 <- jags(data_bundle_2020[[5]], inits = inits_elk2020, params,
                      "./Outputs/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_2020$summary)
  which(RN_elk_2020$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_2020$samples)
  save(RN_elk_2020, file = "./Outputs/JAGS_out/RN_elk_2020.RData")
  
  ######  Moose  ######
  start.time = Sys.time()
  inits_moose2020 <- function(){list(N = ninit_2020[[6]])}
  RN_moose_2020 <- jags(data_bundle_2020[[6]], inits = inits_moose2020, params,
                      "./Outputs/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_2020$summary)
  which(RN_moose_2020$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_2020$samples)
  save(RN_moose_2020, file = "./Outputs/JAGS_out/RN_moose_2020.RData")
  
  ######  White-tailed Deer  ######
  ni_wtd <- 100000
  start.time = Sys.time()
  inits_wtd2020 <- function(){list(N = ninit_2020[[7]])}
  RN_wtd_2020 <- jags(data_bundle_2020[[7]], inits = inits_wtd2020, params,
                      "./Outputs/JAGS_RNmod_2020.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_2020$summary)
  which(RN_wtd_2020$summary[,"Rhat"] > 1.1) 
  mcmcplot(RN_wtd_2020$samples)
  save(RN_wtd_2020, file = "./Outputs/JAGS_out/RN_wtd_2020.RData")
  
  
  #####  2021 & 2022  Analyses  #####
  source("./Scripts/1b_RNmodel_JAGS_code_2021_2022.R")
  
  ######  Black bear  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_bear2021 <- function(){list(N = ninit_2021[[1]])}
  RN_bear_2021 <- jags(data_bundle_2021[[1]], inits = inits_bear2021, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_2021$summary)
  which(RN_bear_2021$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_2021$samples)
  save(RN_bear_2021, file = "./Outputs/JAGS_out/RN_bear_2021.RData")
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_bear2022 <- function(){list(N = ninit_2022[[1]])}
  RN_bear_2022 <- jags(data_bundle_2022[[1]], inits = inits_bear2022, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_bear_2022$summary)
  which(RN_bear_2022$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_bear_2022$samples)
  save(RN_bear_2022, file = "./Outputs/JAGS_out/RN_bear_2022.RData")
  
  ######  Coyote  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_coy2021 <- function(){list(N = ninit_2021[[2]])}
  RN_coy_2021 <- jags(data_bundle_2021[[2]], inits = inits_coy2021, params,
                     "./Outputs/JAGS_RNmod_2021_2022.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_2021$summary)
  which(RN_coy_2021$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_2021$samples)
  save(RN_coy_2021, file = "./Outputs/JAGS_out/RN_coy_2021.RData")
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_coy2022 <- function(){list(N = ninit_2022[[2]])}
  RN_coy_2022 <- jags(data_bundle_2022[[2]], inits = inits_coy2022, params,
                     "./Outputs/JAGS_RNmod_2021_2022.txt",
                     n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                     n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_coy_2022$summary)
  which(RN_coy_2022$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_coy_2022$samples)
  save(RN_coy_2022, file = "./Outputs/JAGS_out/RN_coy_2022.RData")
  
  ######  Mountain lion  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_lion2021 <- function(){list(N = ninit_2021[[3]])}
  RN_lion_2021 <- jags(data_bundle_2021[[3]], inits = inits_lion2021, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_2021$summary)
  which(RN_lion_2021$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_2021$samples)
  save(RN_lion_2021, file = "./Outputs/JAGS_out/RN_lion_2021.RData")
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_lion2022 <- function(){list(N = ninit_2022[[3]])}
  RN_lion_2022 <- jags(data_bundle_2022[[3]], inits = inits_lion2022, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_lion_2022$summary)
  which(RN_lion_2022$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_lion_2022$samples)
  save(RN_lion_2022, file = "./Outputs/JAGS_out/RN_lion_2022.RData")
  
  ######  Wolf  ######
  ni_wolf <- 100000 
  #'  Summer 2021
  start.time = Sys.time()
  inits_wolf2021 <- function(){list(N = ninit_2021[[4]])}
  RN_wolf_2021 <- jags(data_bundle_2021[[4]], inits = inits_wolf2021, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_2021$summary)
  which(RN_wolf_2021$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_2021$samples)
  save(RN_wolf_2021, file = "./Outputs/JAGS_out/RN_wolf_2021.RData")
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_wolf2022 <- function(){list(N = ninit_2022[[4]])}
  RN_wolf_2022 <- jags(data_bundle_2022[[4]], inits = inits_wolf2022, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wolf, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wolf_2022$summary)
  which(RN_wolf_2022$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_wolf_2022$samples)
  save(RN_wolf_2022, file = "./Outputs/JAGS_out/RN_wolf_2022.RData")
  
  ######  Elk  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_elk2021 <- function(){list(N = ninit_2021[[5]])}
  RN_elk_2021 <- jags(data_bundle_2021[[5]], inits = inits_elk2021, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_2021$summary)
  which(RN_elk_2021$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_2021$samples)
  save(RN_elk_2021, file = "./Outputs/JAGS_out/RN_elk_2021.RData")
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_elk2022 <- function(){list(N = ninit_2022[[5]])}
  RN_elk_2022 <- jags(data_bundle_2022[[5]], inits = inits_elk2022, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_elk_2022$summary)
  which(RN_elk_2022$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_elk_2022$samples)
  save(RN_elk_2022, file = "./Outputs/JAGS_out/RN_elk_2022.RData")
  
  ######  Moose  ######
  #'  Summer 2021
  start.time = Sys.time()
  inits_moose2021 <- function(){list(N = ninit_2021[[6]])}
  RN_moose_2021 <- jags(data_bundle_2021[[6]], inits = inits_moose2021, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_2021$summary)
  which(RN_moose_2021$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_2021$samples)
  save(RN_moose_2021, file = "./Outputs/JAGS_out/RN_moose_2021.RData")
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_moose2022 <- function(){list(N = ninit_2022[[6]])}
  RN_moose_2022 <- jags(data_bundle_2022[[6]], inits = inits_moose2022, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_moose_2022$summary)
  which(RN_moose_2022$summary[,"Rhat"] > 1.1)
  mcmcplot(RN_moose_2022$samples)
  save(RN_moose_2022, file = "./Outputs/JAGS_out/RN_moose_2022.RData")
  
  ######  White-tailed Deer  ######
  ni_wtd <- 100000
  #'  Summer 2021
  start.time = Sys.time()
  inits_wtd2021 <- function(){list(N = ninit_2021[[7]])}
  RN_wtd_2021 <- jags(data_bundle_2021[[7]], inits = inits_wtd2021, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_2021$summary)
  which(RN_wtd_2021$summary[,"Rhat"] > 1.1) 
  mcmcplot(RN_wtd_2021$samples)
  save(RN_wtd_2021, file = "./Outputs/JAGS_out/RN_wtd_2021.RData")
  
  #'  Summer 2022
  start.time = Sys.time()
  inits_wtd2022 <- function(){list(N = ninit_2022[[7]])}
  RN_wtd_2022 <- jags(data_bundle_2022[[7]], inits = inits_wtd2022, params,
                      "./Outputs/JAGS_RNmod_2021_2022.txt",
                      n.adapt = na, n.chains = nc, n.thin = nt, n.iter = ni_wtd, 
                      n.burnin = nb, parallel = TRUE)
  end.time <- Sys.time(); (run.time <- end.time - start.time)
  print(RN_wtd_2022$summary)
  which(RN_wtd_2022$summary[,"Rhat"] > 1.1) 
  mcmcplot(RN_wtd_2022$samples)
  save(RN_wtd_2022, file = "./Outputs/JAGS_out/RN_wtd_2022.RData")
  
  
  #'  --------------------------------------
  ####  Extract relative abundance indices  ####
  #'  --------------------------------------
  #'  Load model outputs
  filenames <- list.files("./Outputs/JAGS_out", pattern="*.RData", full.names=TRUE)
  lapply(filenames, load, environment())
  
  #'  List annual model outputs together per species
  rn_bear_list <- list(RN_bear_2020, RN_bear_2021, RN_bear_2022)
  rn_coy_list <- list(RN_coy_2020, RN_coy_2021, RN_coy_2022)
  rn_lion_list <- list(RN_lion_2020, RN_lion_2021, RN_lion_2022)
  rn_wolf_list <- list(RN_wolf_2020, RN_wolf_2021, RN_wolf_2022)
  rn_elk_list <- list(RN_elk_2020, RN_elk_2021, RN_elk_2022)
  rn_moose_list <- list(RN_moose_2020, RN_moose_2021, RN_moose_2022)
  rn_wtd_list <- list(RN_wtd_2020, RN_wtd_2021, RN_wtd_2022)
  
  #'  List one detection history per year (doesn't matter which species it relates to)
  dh_list <- list(DH_smr2020[[1]], DH_smr2021[[1]], DH_smr2022[[1]])
  
  #'  Save estimated N per site
  estimated_N <- function(mod, dh, spp) {
    #'  Grab estimated N and SD per site
    RN.n <- mod$mean$N
    RN.sd <- mod$sd$N
    #'  Grab camera location
    locs <- rownames(dh)
    #'  Merge and format into single data frame with corresponding N & SD per site
    out <- cbind(locs, RN.n, RN.sd)
    RN_est <- as.data.frame(out) %>%
      mutate(CamID = locs, 
             Species = spp,
             RN.n = as.numeric(RN.n),
             RN.sd = as.numeric(RN.sd)) %>%
      separate(locs, c("GMU", "Setup", "site"), sep = "_") %>%
      mutate(CellID = paste0(GMU, "_", site)) %>%
      dplyr::select(-c(site, CellID)) %>%
      relocate(CamID, .before = "GMU") %>%
      relocate(Species, .after = "Setup") 
    return(RN_est)
  }
  rn_bear_out <- mapply(estimated_N, rn_bear_list, dh = dh_list, spp = "bear_black", SIMPLIFY = FALSE)
  rn_coy_out <- mapply(estimated_N, rn_coy_list, dh = dh_list, spp = "coyote", SIMPLIFY = FALSE)
  rn_lion_out <- mapply(estimated_N, rn_lion_list, dh = dh_list, spp = "mountain_lion", SIMPLIFY = FALSE)
  rn_wolf_out <- mapply(estimated_N, rn_wolf_list, dh = dh_list, spp = "wolf", SIMPLIFY = FALSE)
  rn_elk_out <- mapply(estimated_N, rn_elk_list, dh = dh_list, spp = "elk", SIMPLIFY = FALSE)
  rn_moose_out <- mapply(estimated_N, rn_moose_list, dh = dh_list, spp = "moose", SIMPLIFY = FALSE)
  rn_wtd_out <- mapply(estimated_N, rn_wtd_list, dh = dh_list, spp = "whitetailed_deer", SIMPLIFY = FALSE)
  
  #'  Merge site-specific N for all species per year
  rn_2020 <- rbind(rn_bear_out[[1]], rn_coy_out[[1]], rn_lion_out[[1]], rn_wolf_out[[1]],
                   rn_elk_out[[1]], rn_moose_out[[1]], rn_wtd_out[[1]]) %>%
    arrange(CamID, Species) %>%
    mutate(season = "Smr20") %>%
    relocate(season, .after = Setup)
  rn_2021 <- rbind(rn_bear_out[[2]], rn_coy_out[[2]], rn_lion_out[[2]], rn_wolf_out[[2]],
                   rn_elk_out[[2]], rn_moose_out[[2]], rn_wtd_out[[2]]) %>%
    arrange(CamID, Species) %>%
    mutate(season = "Smr21") %>%
    relocate(season, .after = Setup)
  rn_2022 <- rbind(rn_bear_out[[3]], rn_coy_out[[3]], rn_lion_out[[3]], rn_wolf_out[[3]],
                   rn_elk_out[[3]], rn_moose_out[[3]], rn_wtd_out[[3]]) %>%
    arrange(CamID, Species) %>%
    mutate(season = "Smr22") %>%
    relocate(season, .after = Setup)
 
  RN_abundance <- list(rn_2020, rn_2021, rn_2022)
  save(RN_abundance, file = "./Outputs/RN_abundance.RData")
  
  #'  -----------------------------------------
  #####  Visualize relative abundance indices  #####
  #'  -----------------------------------------
  #'  THE FOLLOWING CODE IS PROVIDED FOR TRANSPARENCY BUT CANNOT BE RUN WITH 
  #'  PUBLICLY AVAILABLE DATA. CAMERA LOCATIONS ARE SENSITIVE BUT AVAILABLE
  #'  UPON REQUEST FROM IDAHO DEPARTMENT OF FISH & GAME FOR QUALIFIED RESEARCHERS.
  
  #'  Map relative abundance data per species, study area, and year
  library(sf)
  library(ggplot2)
  library(patchwork)
  
  #'  Load spatial data
  #'  GMU shapefiles provided by Idaho Department of Fish and Game
  gmu_wgs84 <- st_read("./Data/Spatial_data/Study_GMUs.shp") %>%
    st_transform("+proj=longlat +datum=WGS84 +no_defs")
  
  #'  Camera site coordinates have been adjusted slightly to maintain location anonymity  
  cams_2020_wgs84 <- st_read("./Data/Spatial_data/cams_2020_wgs84.shp")
  cams_2021_wgs84 <- st_read("./Data/Spatial_data/cams_2021_wgs84.shp")
  cams_2022_wgs84 <- st_read("./Data/Spatial_data/cams_2022_wgs84.shp")
  
  #'  Original water data: National Hydrology Database Idaho State
  bigwater <- st_read("./Data/Spatial_data/Idaho_waterbodies_1km2.shp")
  bigwater_wgs84 <- st_transform(bigwater, crs = "+proj=longlat +datum=WGS84 +no_defs") %>%
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
  
  #'  List camera spatial data
  cam_list <- list(cams_2020_wgs84, cams_2021_wgs84, cams_2022_wgs84)
  
  #'  Append RN local abundance estimates to spatial data
  spatial_rn <- function(rn, spp, cams) {
    #'  Filter data to single species
    single_spp_rn <- rn %>%
      filter(Species == spp) %>%
      mutate(RN.n.rounded = round(RN.n, 0))
    
    #'  Join spatial data with rn data
    rn_shp <- full_join(cams, single_spp_rn, by = "CamID") %>%
      filter(!is.na(Species)) 
    
    return(rn_shp)
  }
  spatial_rn_bear <- mapply(rn = rn_bear_out, spatial_rn, spp = "bear_black", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_coy <- mapply(rn = rn_coy_out, spatial_rn, spp = "coyote", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_lion <- mapply(rn = rn_lion_out, spatial_rn, spp = "mountain_lion", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_wolf <- mapply(rn = rn_wolf_out, spatial_rn, spp = "wolf", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_elk <- mapply(rn = rn_elk_out, spatial_rn, spp = "elk", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_moose <- mapply(rn = rn_moose_out, spatial_rn, spp = "moose", cams = cam_list, SIMPLIFY = FALSE)
  spatial_rn_wtd <- mapply(rn = rn_wtd_out, spatial_rn, spp = "whitetailed_deer", cams = cam_list, SIMPLIFY = FALSE)
  
  #'  List spatial RN local abundance data and save
  spatial_RN_list <- list(spatial_rn_bear, spatial_rn_coy, spatial_rn_lion, spatial_rn_wolf,
                          spatial_rn_elk, spatial_rn_moose, spatial_rn_wtd)
  save(spatial_RN_list, file = "./Outputs/Spatial_outputs/spatial_RN_list.RData")
  
  #'  Merge wolf results into single large spatial data frame and save as a shapefile
  spatial_rn_wolf_locs <- bind_rows(spatial_rn_wolf[[1]], spatial_rn_wolf[[2]], spatial_rn_wolf[[3]])
  st_write(spatial_rn_wolf_locs, "./Outputs/Spatial_outputs/spatial_rn_wolf_locs.shp")
  
  year_list <- list("2020", "2021", "2022")
  
  #'  Add year to each dataframe and unlist into one single large dataframe per species
  add_yr <- function(dat, yr) {
    dat$Year <- yr
    return(dat)
  }
  rn_bear_all <- mapply(add_yr, dat = spatial_rn_bear, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_coy_all <- mapply(add_yr, dat = spatial_rn_coy, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_lion_all <- mapply(add_yr, dat = spatial_rn_lion, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_wolf_all <- mapply(add_yr, dat = spatial_rn_wolf, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_elk_all <- mapply(add_yr, dat = spatial_rn_elk, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_moose_all <- mapply(add_yr, dat = spatial_rn_moose, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  rn_wtd_all <- mapply(add_yr, dat = spatial_rn_wtd, yr = year_list, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  #'  Make one giant faceted plot where rows represent GMU and columns represent years 
  #'  to ensure that the dot sizes are all consistent for at least a single species
  library(ggh4x)
  map_rn <- function(sf_rn, spp) {
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
      facet_wrap(~Year) 
      
    #'  Plot each map
    plot(spp_rn)
    
    return(spp_rn)
  }
  rn_maps_bear <- map_rn(rn_bear_all, spp = "black bears")
  rn_maps_coy <- map_rn(rn_coy_all, spp = "coyotes")
  rn_maps_lion <- map_rn(rn_lion_all, spp = "mountain lions")
  rn_maps_wolf <- map_rn(rn_wolf_all, spp = "wolves")
  rn_maps_elk <- map_rn(rn_elk_all, spp = "elk")
  rn_maps_moose <- map_rn(rn_moose_all, spp = "moose")
  rn_maps_wtd <- map_rn(rn_wtd_all, spp = "white-tailed deer")
  
  #'  Save figures
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
  
  
  #'  Next step:
  #'  Cluster cameras based on wolf relative abundance index (script 2_Cluster_cameras.R)
  