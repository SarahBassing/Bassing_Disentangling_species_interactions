  #'  --------------------------------
  #'  4. Structural Equation Models
  #'  
  #'  R code from Bassing, S.B., D.E. Ausband, M.A. Mumma, J. Baumgardt, S. Thompson, 
  #'  M.A. Hurley, and M. Falcy. 2025. "Disentangling the web of species interactions 
  #'  in a multi-predator, multi-prey community"
  #'  --------------------------------
  #'  Script to run structural equation models. Script sources3_Relative_density_indices.R 
  #'  script and runs competing structural equation models (SEM) to test hypotheses 
  #'  about how predator-prey and predator-predator interactions influence the
  #'  relative densities of large mammals in a multi-predator, multi-prey community.
  #'  Initial model representing a hypothesis is first presented, then a second
  #'  version of the model is presented that was adjusted based on tests of 
  #'  d-separation and add/removing relationships between nodes.
  #'  
  #'  Required: 
  #'    1. Source 3_Relative_density_indices.R which includes relative density 
  #'    indices and covariates formatted for SEM (density_wide_1YrLag_20s_22s object)
  #'  --------------------------------
  
  #'  Clean workspace
  rm(list = ls())

  library(piecewiseSEM)
  library(semEff)
  library(labelled)
  library(DiagrammeR)
  library(lme4)
  library(tidyverse)
  
  #'  Run script that formats density data for SEMs
  source("./Scripts/3_Relative_density_indices.R")
  
  #'  Take a quick look
  head(density_wide_1YrLag_20s_22s)
  
  #'  Set options so all no rows are omitted in model output
  options(max.print = 9999)
  
  
  #'  -------------------------------------------------
  ####  SEM with 1-year time lag, no annual variation  ####
  #'  -------------------------------------------------
  #'  ---------------------------------------
  #####  Top down, interference competition  #####
  #'  ---------------------------------------
  top_down_inter <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + annual_harvest.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_inter)
  AIC_psem(top_down_inter, AIC.type = "loglik")
  dSep(top_down_inter)

  #'  Reduced model based on tests of d-separation
  #'  For all linear models at once-
  #'    1. Remove non-significant relationships that were initially hypothesized, rerun
  #'    2. Add newly significant relationships identified by d-sep test that were initially hypothesized, rerun one at a time
  #'    3. Add newly significant relationships identified by d-sep test that were not initially hypothesized but biologically plausible, rerun one at a time
  #'    4. Add newly significant relationships identified by d-sep test that were not initially expected to have a causal relationship, rerun
  top_down_inter.dsep <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.T, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1 + moose.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(top_down_inter.dsep)
  
  
  #'  ---------------------------------------
  #####  Top-down, exploitation competition  #####
  #'  ---------------------------------------
  top_down_exploit <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + wolf.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + annual_harvest.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  ) 
  summary(top_down_exploit)
  
  #'  Reduced model based on tests of d-separation
  top_down_exploit.dsep <- psem(
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + mountain_lion.Tminus1 + bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(wolf.T ~ wolf.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + mountain_lion.Tminus1 + elk.T, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T, data = density_wide_1YrLag_20s_22s),
    lm(coyote.T ~ coyote.Tminus1 + moose.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  ) 
  summary(top_down_exploit.dsep)
  
  
  #'  ----------------------------------------
  #####  Bottom up, interference competition  #####
  #'  ----------------------------------------
  bottom_up_inter <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + bear_black.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1 + wolf.Tminus1 + mountain_lion.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s), 
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_inter)
  AIC_psem(bottom_up_inter, AIC.type = "loglik")
  
  #'  Reduced model based on tests of d-separation
  bottom_up_inter.dsep <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + bear_black.T, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + bear_black.Tminus1 + wolf.T, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s), 
    coyote.T %~~% mountain_lion.T,
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_inter.dsep)
  
  
  #'  --------------
  #####  Bottom up, exploitation competition  #####
  #'  --------------
  bottom_up_exploit <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1 + elk.Tminus1 + whitetailed_deer.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1 + whitetailed_deer.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1 + DisturbedForest_last20Yrs.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_exploit)
  AIC_psem(bottom_up_exploit, AIC.type = "loglik")
  
  #'  Reduced model based on tests of d-separation
  bottom_up_exploit.dsep <- psem(
    lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s),
    lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s), 
    lm(coyote.T ~ coyote.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s),
    lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + bear_black.T, data = density_wide_1YrLag_20s_22s),
    lm(elk.T ~ elk.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + bear_black.Tminus1 + wolf.T, data = density_wide_1YrLag_20s_22s),
    lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s),
    data = density_wide_1YrLag_20s_22s
  )
  summary(bottom_up_exploit.dsep)
  
  #'  -----------------------------
  ####  Model selection using AIC  ####
  #'  -----------------------------
  modSelect <- AIC(top_down_inter.dsep, top_down_exploit.dsep, bottom_up_inter.dsep, bottom_up_exploit.dsep)  
  mod_names <- c("top_down_inter.dsep", "top_down_exploit.dsep", "bottom_up_inter.dsep", "bottom_up_exploit.dsep") 
  modSelect <- bind_cols(mod_names, modSelect)
  names(modSelect) <- c("Model", "AIC", "K", "n")
  (sem_aic <- arrange(modSelect, AIC, decreasing = TRUE))
  
  #'  Review Fisher's C and Chi-square test statistics for each model
  sem_gof <- function(mod, gof, modname) {
    model <- modname
    chisq <- LLchisq(mod)
    chisq_stat <- round(chisq[[1]], 2)
    chisq_df <- chisq[[2]]
    chisq_pval <- round(chisq[[3]], 2)
    fishers <- fisherC(mod)
    fishers_stat <- round(fishers[[1]], 2)
    fishers_df <- fishers[[2]]
    fishers_pval <- round(fishers[[3]], 2)
    GoF <- data.frame(model, chisq_stat, chisq_df, chisq_pval, fishers_stat, fishers_df, fishers_pval)
    names(GoF) <- c("Model", "Chi2", "Chi2 df", "Chi2 p-value", "Fisher's C", "Fisher's C df", "Fisher's C p-value")
    return(GoF)
  }
  mod_list <- list(top_down_inter.dsep, top_down_exploit.dsep, bottom_up_inter.dsep, bottom_up_exploit.dsep) 
  mod_names <- c("top_down_inter.dsep", "top_down_exploit.dsep", "bottom_up_inter.dsep", "bottom_up_exploit.dsep")
  sem_GoF <- mapply(sem_gof, mod_list, modname = mod_names, SIMPLIFY = FALSE) %>% bind_rows()
  
  modSel_table <- full_join(sem_aic, sem_GoF, by = "Model") %>%
    mutate(AIC = round(AIC, 2),
           deltaAIC = AIC - first(x = AIC),
           Model = ifelse(Model == "top_down_inter.dsep", "Top-down, interference competition", Model),
           Model = ifelse(Model == "top_down_exploit.dsep", "Top-down, exploitation competition", Model),
           Model = ifelse(Model == "bottom_up_inter.dsep", "Bottom-up, interference competition", Model),
           Model = ifelse(Model == "bottom_up_exploit.dsep", "Bottom-up, exploitation competition", Model)) %>%
    relocate(deltaAIC, .after = AIC) %>%
    dplyr::select(-n)
  write_csv(modSel_table, "./Outputs/ModSel_GoF_table.csv")
  
  #'  ----------------------
  ####  Dig into top model  ####
  #'  ----------------------
  #'  Dig into top model based on AIC, Fisher's C
  #'  Check for multicollinearity
  RVIF(bottom_up_exploit.dsep[[1]]) 
  RVIF(bottom_up_exploit.dsep[[2]]) 
  RVIF(bottom_up_exploit.dsep[[3]]) 
  RVIF(bottom_up_exploit.dsep[[4]]) 
  RVIF(bottom_up_exploit.dsep[[5]]) 
  RVIF(bottom_up_exploit.dsep[[6]]) 
  RVIF(bottom_up_exploit.dsep[[7]]) 
  
  RVIF(bottom_up_inter.dsep[[1]]) 
  RVIF(bottom_up_inter.dsep[[2]]) 
  RVIF(bottom_up_inter.dsep[[3]]) 
  RVIF(bottom_up_inter.dsep[[4]]) 
  RVIF(bottom_up_inter.dsep[[5]]) 
  RVIF(bottom_up_inter.dsep[[6]]) 
  RVIF(bottom_up_inter.dsep[[7]]) 
  
  #'  Run some basic model diagnostics on individuals models keeping in mind
  #'  these linear regressions ignore important indirect effects the SEM identified
  wtd_mod <- lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + bear_black.T, data = density_wide_1YrLag_20s_22s)
  plot(wtd_mod)
  elk_mod <- lm(elk.T ~ elk.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + bear_black.Tminus1 + wolf.T, data = density_wide_1YrLag_20s_22s)
  plot(elk_mod)
  moose_mod <- lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(moose_mod)
  lion_mod <- lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(lion_mod)
  bear_mod <- lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(bear_mod)
  wolf_mod <- lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  plot(wolf_mod)
  coy_mod <- lm(coyote.T ~ coyote.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s)
  plot(coy_mod)
  
  #'  ------------------
  #####  Visualize SEM  #####
  #'  ------------------
  #'  Couple ways to visualize SEMs with network diagrams
  piecewiseSEM:::plot.psem(bottom_up_exploit.dsep, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"),
                           layout = "circle",
                           output = "visNetwork") 
  piecewiseSEM:::plot.psem(bottom_up_exploit.dsep, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  
  piecewiseSEM:::plot.psem(bottom_up_inter.dsep, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(bottom_up_inter.dsep, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  
  piecewiseSEM:::plot.psem(top_down_inter.dsep, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"),
                           layout = "circle",
                           output = "visNetwork") 
  piecewiseSEM:::plot.psem(top_down_inter.dsep, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))
  
  piecewiseSEM:::plot.psem(top_down_exploit.dsep, 
                           node_attrs = data.frame(shape = "rectangle", color = "orange", fontcolor = "black"), 
                           layout = "circle",
                           output = "visNetwork")
  piecewiseSEM:::plot.psem(top_down_exploit.dsep, 
                           node_attrs = data.frame(shape = "rectangle", fontcolor = "black", fillcolor = "orange"),
                           edge_attrs = data.frame(style = "solid", color = "black"))

  
  #'  -------------------------------
  #####  DIRECT vs INDIRECT EFFECTS  #####
  #'  -------------------------------
  #'  Calculate direct, indirect, total, and mediator effects (SE & 95% CI) for 
  #'  all endogenous (response) variables from top SEM using semEff package 
  #'  https://murphymv.github.io/semEff/articles/semEff.html
  #'  
  #'  First: bootstrap standardized model coefficients (necessary for calculating SEs)
  #'  THIS TAKES AWHILE! Especially with R = 10000
  #'  
  #'  Second: calculate standardized effects for all casual pathways
  #'  Note, these are standardized unique effects (i.e., adjusted for multicollinearity;
  #'  i.e., semipartial correlations), allowing us to fully partition effects in the system
  #'  These tend to be smaller than the unadjusted standardized coefficients
  #'  If there are large differences btwn the two then consideration should be given
  #'  to the impact and relevance of multicollinearity in the system (check with RVIF())
  bottom_up_exploit.dsep_bootEff <- bootEff(bottom_up_exploit.dsep, R = 10000, seed = 13, 
                                            type = "nonparametric", parallel = "multicore", ncpus = 5) 
  (bottom_up_exploit.dsep_semEff <- semEff(bottom_up_exploit.dsep_bootEff))
  summary(bottom_up_exploit.dsep_semEff)
  
  #'  Save since this takes forever
  save(bottom_up_exploit.dsep_bootEff, file = "./Outputs/bottom_up_exploit.dsep_bootEff.RData")
  save(bottom_up_exploit.dsep_semEff, file = "./Outputs/bottom_up_exploit.dsep_semEff.RData")
  
  #'  ------------------
  #####  Result tables  #####
  #'  ------------------
  ######  Direct & indirect effects table  ######
  #'  -------------------------------------
  #'  Extract standardized direct, indirect, etc. effects from SEMs after bootstrapping
  #'  These values are what go into the final network diagram figure showing
  #'  direct and indirect relationships and 95% CIs
  library(stringi)
  #'  Create results tables
  tbl_std_est <- function(mod, mod_name) {
    #'  Extract standardized results
    predictor_wtd <- mod$Summary$whitetailed.deer.T
    predictor_elk <- mod$Summary$elk.T
    predictor_moose <- mod$Summary$moose.T
    predictor_wolf <- mod$Summary$wolf.T
    predictor_lion <- mod$Summary$mountain.lion.T
    predictor_bear <- mod$Summary$bear.black.T
    predictor_coy <- mod$Summary$coyote.T
    
    #'  Bind into a single data frame and add column for endogenous variable
    tbl_wtd <- as.data.frame(predictor_wtd) %>% bind_cols("White-tailed deer t")
    tbl_elk <- as.data.frame(predictor_elk) %>% bind_cols("Elk t")
    tbl_moose <- as.data.frame(predictor_moose) %>% bind_cols("Moose t")
    tbl_wolf <- as.data.frame(predictor_wolf) %>% bind_cols("Wolf t")
    tbl_lion <- as.data.frame(predictor_lion) %>% bind_cols("Mountain lion t")
    tbl_bear <- as.data.frame(predictor_bear) %>% bind_cols("Black bear t")
    tbl_coy <- as.data.frame(predictor_coy) %>% bind_cols("Coyote t")
    
    #'  Rename columns (this is weird thanks to odd formatting from semEff output)
    col_names <- c("Effect type", "Exogenous_variable", "space", "Standardized_effect", 
                   "space", "Bias", "space", "Std.Error", "space", "lower_CI", "upper_CI", 
                   "space", "Signif", "Endogenous_variable")
    names(tbl_wtd) <- col_names
    names(tbl_elk) <- col_names
    names(tbl_moose) <- col_names
    names(tbl_wolf) <- col_names
    names(tbl_lion) <- col_names
    names(tbl_bear) <- col_names
    names(tbl_coy) <- col_names
    
    #'  Bind results together and clean up final table
    full_tbl <- bind_rows(tbl_wtd, tbl_elk, tbl_moose, tbl_wolf, tbl_lion, tbl_bear, tbl_coy) %>%
      relocate(Endogenous_variable, .before = Exogenous_variable) %>%
      dplyr::select(c("Endogenous_variable", "Exogenous_variable", "Standardized_effect", 
                      "Std.Error", "lower_CI", "upper_CI", "Signif")) %>%
      #'  Remove extra spaces before or end of words
      mutate(across(everything(), ~stri_trim(.))) %>%
      #'  Replace empty cells with NA
      mutate(across(everything(), ~na_if(., ""))) %>%
      #'  Remove any rows with `NA`
      # filter(!if_any(everything(), is.na)) %>%
      filter(!is.na(Exogenous_variable)) %>%
      filter(Exogenous_variable != "n/a") %>%
      #'  Remove duplicate information
      distinct(.) %>%
      #'  Create single column for confidence intervals
      mutate("95% CI" = paste(" ", lower_CI, "-", upper_CI)) %>%
      dplyr::select(-c(lower_CI, upper_CI)) %>%
      arrange(Endogenous_variable, Exogenous_variable) %>%
      #'  Add column reporting model name
      mutate(Model = mod_name) %>%
      relocate(Model, .before = Endogenous_variable) %>%
      relocate(Signif, .after = "95% CI") %>%
      #'  Remove "." and adjust time period indicator
      mutate(Exogenous_variable = ifelse(Exogenous_variable == "bear.black.Tminus1", "bear black.Tminus1", Exogenous_variable),
             Exogenous_variable = ifelse(Exogenous_variable == "mountain.lion.Tminus1", "mountain lion.Tminus1", Exogenous_variable),
             Exogenous_variable = ifelse(Exogenous_variable == "whitetailed.deer.Tminus1", "white-tailed deer.Tminus1", Exogenous_variable),
             Exogenous_variable = gsub(".Tminus1", " t-1", Exogenous_variable),
             Exogenous_variable = str_to_sentence(Exogenous_variable))
    
    return(full_tbl)
  }
  result_tbl_bottom_up_exploit.dsep <- tbl_std_est(bottom_up_exploit.dsep_semEff, mod_name = "Bottom-up, Exploit")
  write_csv(result_tbl_bottom_up_exploit.dsep, "./Outputs/result_table_top_model_std_effects.csv")
  
  #'  --------------------------------------------
  ######  Unstandardized regression coefficients  ######
  #'  --------------------------------------------
  #'  Run linear regression models from best supported SEM
  wtd_mod <- lm(whitetailed_deer.T ~ whitetailed_deer.Tminus1 + bear_black.T, data = density_wide_1YrLag_20s_22s)
  elk_mod <- lm(elk.T ~ elk.Tminus1 + whitetailed_deer.Tminus1 + whitetailed_deer.T + bear_black.Tminus1 + wolf.T, data = density_wide_1YrLag_20s_22s)
  moose_mod <- lm(moose.T ~ moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  wolf_mod <- lm(wolf.T ~ wolf.Tminus1 + moose.Tminus1, data = density_wide_1YrLag_20s_22s)
  lion_mod <- lm(mountain_lion.T ~ mountain_lion.Tminus1 + elk.Tminus1, data = density_wide_1YrLag_20s_22s)
  bear_mod <- lm(bear_black.T ~ bear_black.Tminus1, data = density_wide_1YrLag_20s_22s)
  coy_mod <- lm(coyote.T ~ coyote.Tminus1 + mountain_lion.T, data = density_wide_1YrLag_20s_22s)
  
  #'  List models
  regression_list <- list(wtd_mod, elk_mod, moose_mod, wolf_mod, lion_mod, bear_mod, coy_mod)
  response_var <- list("White-tailed deer t", "Elk t", "Moose t", "Wolf t", "Mountain lion t", "Black bear t", "Coyote t")
  
  #'  Create results table of unstandardized regression coefficients and R2 values
  SEM_unstandardized_tbl <- function(mod, response_var) {
    #'  Grab coefficient estimates, SE, and p-value from linear model
    Params <- rownames(as.data.frame(mod$coefficients))
    mod_summary <- summary(mod)
    coefEst <- round(mod_summary$coefficients[,1], 2)
    stdErr <- round(mod_summary$coefficients[,2], 2)
    pVal <- round(mod_summary$coefficients[,4], 5)
    r2 <- round(mod_summary$r.squared, 2)
    
    #'  Merge into a single data frame
    lm_out <- bind_cols(response_var, Params, coefEst, stdErr, pVal, r2) 
    names(lm_out) <- c("Response", "Parameter", "Estimate", "SE", "p-value", "R2")
    #'  Clean up parameter names - remove "." and adjust time period indicator
    lm_out <- lm_out %>%
      mutate(Parameter = ifelse(Parameter == "bear_black.Tminus1", "bear black.Tminus1", Parameter),
             Parameter = ifelse(Parameter == "bear_black.T", "bear black.T", Parameter),
             Parameter = ifelse(Parameter == "mountain_lion.Tminus1", "mountain lion.Tminus1", Parameter),
             Parameter = ifelse(Parameter == "mountain_lion.T", "mountain lion.T", Parameter),
             Parameter = ifelse(Parameter == "whitetailed_deer.Tminus1", "white-tailed deer.Tminus1", Parameter),
             Parameter = ifelse(Parameter == "whitetailed_deer.T", "white-tailed deer.T", Parameter),
             Parameter = gsub(".Tminus1", " t-1", Parameter),
             Parameter = gsub(".T", " t", Parameter),
             Parameter = str_to_sentence(Parameter))
   
    return(lm_out) 
  }
  sem_regression_tbl <- mapply(SEM_unstandardized_tbl, regression_list, response_var = response_var, SIMPLIFY = FALSE) %>% bind_rows(.)
  
  #'  Save for publication
  write_csv(sem_regression_tbl, "./Outputs/results_table_regres_coefs_tbl.csv")
  
  
  #'  Next step:
  #'  Generate summary and result tables for publication (script 5_Result_Tables_and_Summaries.R)
  