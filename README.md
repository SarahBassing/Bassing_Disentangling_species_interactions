# Bassing_Disentangling_species_interactions
R code and data associated with Bassing, Ausband, Mumma, Baumgardt, Thompson, Hurley, &amp; Falcy. 2025. Disentangling the web of species interactions in a multi-predator, multi-prey community

Anonymized data formatted as detection histories are described below. Some portions of this repository are provided for transparency but cannot be run with the publicly available data provided in this repository.
Qualified researchers can contact the Data Management Lead, Idaho Department of Fish and Game, Idaho Fish and Wildlife Information System, 600 S Walnut, Boise, ID 83712 for full data requests: idfgdatarequests@idfg.idaho.gov 

Repository contains 8 scripts, anonymized data, and shapefiles needed for some spatial analyses.
"Scripts" folder:
  1_Relative_abundance_Royle-Nichols_models.R: script to estimate relative abundance indices (RAI)
  1a_RNmodel_JAGS_code_2020.R: JAGS code for Royle-Nichols abundance model based on 2020 data
  1b_RNmodel_JAGS_code_2021_2022.R: JAGS code for Royle-Nichols abundance model based on 2021 & 2022 data (differ from 2020 code because requires slightly different model parameterization with addition of GMU1 study area)
  2_Cluster_cameras.R: script to group cameras based on wolf RAI and geographic distance (CANNOT BE RUN WITH PUBLICLY AVAILABLE DATA PROVIDED IN THIS REPOSITORY)
  3_Relative_density_indices.R: script to generate relative density indices using RAIs
  4_Structural_Equation_Models.R: script to test hypotheses using structural equation models (SEM)
  5_Result_Tables_and_Summaries.R: script to generate result tables for publication
  6_Result_Figures.R: script to generate figures for publication (MOST CANNOT BE RUN WITH PUBLICALLY AVAILABLE DATA PROVIDED IN THIS REPOSITORY)
  
"Data.zip" folder: Once unzipped, will contain
  1. "./Data/Spatial_data/"
    1. ClarkForkRiver_flowline: spatial line data for Clark Fork River flowing through study area
    2. Game_Management_Units: spatial polygon data for Idaho Dept. Fih & Game GMUs
    3. Idaho_waterbodies_1km2: spatial polygon data for waterbodies in study areas >=1 sq-km
    4. IdahoState: spatial polygon data for Idaho State administrative boundary
    5. KootenaiRiver_flowline: spatial line data for Kootenai River flowing through study area
    6. NorthForkClearwater_flowline: spatial line data for N. Fork Clearwater River flowing 
        through study area
    7. Pendoreille_flowline: spatial line data for Pendoreille River flowing through study area
    8. Study_GMUs: spatial polygon data for GMUs 1, 6, and 10A
    9. US_states: spatial polygon data for USA administrative boundary
  2. covariates.RData
      Covaraite data included in Structural Equation Models
  3. DH_smr2020.RData, DH_smr2021.RData, DH_smr2022.Rdata
      Species and season specific detection histories used for Royle-Nichols abundance models 
  4. stations_smr2020.RData, stations_srm2021.RData, station_smr2022.RData
      Camera site information used for Royle-Nichols abundance models and to calculate
      Relative Density Indices for each cluster based on cameras assigned to each cluster
  
"Outputs.zip" folder: Once unzipped, will contain
  1. "./Outputs/JAGS_out/"
      Currently empty but is automatically populated with model outputs when Royle-Nichols 
      abundance models are run using 1_Relative_abundance_Royle-Nichols_models.R
  2. "./Outputs/Spatial_outputs/"
      Currently contains cluster_polygons_gmu1, cluster_polygons_gmu6, and cluster_polygons_gmu10a, 
      spatial polygon data for GMU-specific cluster polygons created in 2_Cluster_cameras.R. 
      These cannot be generated with the publicly available dataso are provided here. 
      1_Relative_abundance_Royle-Nichols_models.R and 2_Cluster_cameras.R save additional spatial 
      outputs (if scripts were run in their entirety)
  3. "./Outputs" 
      All scripts will automatically save various output files to this folder

