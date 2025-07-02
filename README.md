# Bassing_Disentangling_species_interactions
R code and data associated with Bassing, Ausband, Mumma, Baumgardt, Thompson, Hurley, &amp; Falcy. 2025. Disentangling the web of species interactions in a multi-predator, multi-prey community

Anonymized data formatted as detection histories are described below. Some portions of this repository are provided for transparency but cannot be run with the publicly available data provided in this repository.
Qualified researchers can contact the Data Management Lead, Idaho Department of Fish and Game, Idaho Fish and Wildlife Information System, 600 S Walnut, Boise, ID 83712 for full data requests: idfgdatarequests@idfg.idaho.gov 

Repository contains 8 scripts, anonymized data, and shapefiles needed for some spatial analyses.

_"Scripts" folder:_
1. 1_Relative_abundance_Royle-Nichols_models.R: script to estimate relative abundance indices (RAI)
2. 1a_RNmodel_JAGS_code_2020.R: JAGS code for Royle-Nichols abundance model based on 2020 data
3. 1b_RNmodel_JAGS_code_2021_2022.R: JAGS code for Royle-Nichols abundance model based on 2021 & 2022 data (differ from 2020 code because requires slightly different model parameterization with addition of GMU1 study area)
4. 2_Cluster_cameras.R: script to group cameras based on wolf RAI and geographic distance (CANNOT BE RUN WITH PUBLICLY AVAILABLE DATA PROVIDED IN THIS REPOSITORY)
5. 3_Relative_density_indices.R: script to generate relative density indices using RAIs
6. 4_Structural_Equation_Models.R: script to test hypotheses using structural equation models (SEM)
7. 5_Result_Tables_and_Summaries.R: script to generate result tables for publication
8. 6_Result_Figures.R: script to generate figures for publication (MOST CANNOT BE RUN WITH PUBLICALLY AVAILABLE DATA PROVIDED IN THIS REPOSITORY)
  
_"Data.zip" folder: 
Once unzipped, will contain_
  1. "./Data/Spatial_data/"
     
     **ClarkForkRiver_flowline**: spatial line data for Clark Fork River flowing through study area;
     **Game_Management_Units**: spatial polygon data for Idaho Dept. Fih & Game GMUs;
     **Idaho_waterbodies_1km2**: spatial polygon data for waterbodies in study areas >=1 sq-km;
     **IdahoState**: spatial polygon data for Idaho State administrative boundary;
     **KootenaiRiver_flowline**: spatial line data for Kootenai River flowing through study area;
     **NorthForkClearwater_flowline**: spatial line data for N. Fork Clearwater River flowing through study area;
     **Pendoreille_flowline**: spatial line data for Pendoreille River flowing through study area;
     **Study_GMUs**: spatial polygon data for GMUs 1, 6, and 10A;
     **US_states**: spatial polygon data for USA administrative boundary
  2. "./Data/"
     
      **covariates.RData**: Covaraite data included in Structural Equation Models,
      **DH_smr2020.RData, DH_smr2021.RData, DH_smr2022.Rdata**: Species and season specific detection histories used for Royle-Nichols abundance models,
      **stations_smr2020.RData, stations_srm2021.RData, station_smr2022.RData**: Camera site information used for Royle-Nichols abundance models and to calculate Relative Density Indices for each cluster based on cameras assigned to each cluster
  
_"Outputs.zip" folder: 
Once unzipped, will contain_
  1. "./Outputs/JAGS_out/"
     
      Currently empty but is automatically populated with model outputs when Royle-Nichols 
      abundance models are run using 1_Relative_abundance_Royle-Nichols_models.R
  2. "./Outputs/Spatial_outputs/"
     
      Currently contains **cluster_polygons_gmu1, cluster_polygons_gmu6, and cluster_polygons_gmu10a** shapefiles- 
      spatial polygon data for GMU-specific cluster polygons created in 2_Cluster_cameras.R. 
      These cannot be generated with the publicly available dataso are provided here. 
      1_Relative_abundance_Royle-Nichols_models.R and 2_Cluster_cameras.R save additional spatial 
      outputs (if scripts were run in their entirety)
  3. "./Outputs/"
     
      All scripts will automatically save various output files to this folder

