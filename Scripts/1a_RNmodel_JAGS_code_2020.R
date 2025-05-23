  #'  -------------------------------
  #'  1a. Royle-Nichols abundance model (2020 data)
  #'  
  #'  R code from Bassing, S.B., D.E. Ausband, M.A. Mumma, J. Baumgardt, S. Thompson, 
  #'  M.A. Hurley, and M. Falcy. 2025. "Disentangling the web of species interactions 
  #'  in a multi-predator, multi-prey community"
  #'  
  #'  -------------------------------
  #'  RN model to estimate relative abundance from binary detection/non-detection 
  #'  data, assuming heterogeneous abundance affects detection probability.
  #'  
  #'  Relevant parameters and data:
  #'  lambda: The number of animals available for detection at site i, N[i], is a 
  #'  Poisson-distributed random variable with mean lambda[i]. N is latent.
  #'  
  #'  y: The observed detection/non-detection data for each site i and survey 
  #'  occasion j, y[i,j], is a Bernoulli-distributed random variable with  
  #'  probability p[i,j]. y[i,j] is conditional on N[i].
  #'  
  #'  p & r: The probability of detecting occupancy at site i during occasion j, 
  #'  p[i,j], is a function of the per-individual detection probability, r[i,j], 
  #'  a binomial sampling probability that a particular individual is detected at 
  #'  site i during occasion j, and the number of individuals at site i, N[i].
  #'  -------------------------------

  cat(file = './Outputs/JAGS_RNmod_2020.txt', "
      model{
      
      #'  Define priors
      #'  -------------
      #'  Abundance priors
      beta0 ~ dunif(-10, 10)      # Abundance intercept
      mean.lambda <- exp(beta0)   # Mean lambda for GMU10A
      
      #'  Categorical effect for GMU needs multiple beta1 coefficients
      beta1[1] <- 0
      for(gmu in 2:ngmu) {
        beta1[gmu] ~ dnorm(0, 0.001)
      }
      
      #'  Detection priors
      mean.r ~ dunif(0, 1)        # Detection intercept (on probability scale)
      alpha0 <- logit(mean.r)     # Detection intercept (on logit scale)
      
      #'  Categorical effect for camera setup needs multiple alpha1 coefficients
      alpha1[1] <- 0
      for(cam in 2:nsets) {
        alpha1[cam] ~ dnorm(0,0.001)
      }
      
      
      #'  Define likelihood
      #'  -----------------
      #'  Latent state (abundance)
      for(i in 1:nsites){
        N[i] ~ dpois(lambda[i])
        lambda[i] <- exp(beta0 + beta1[gmu[i]])
        
        #'  Detection state
        for(j in 1:nsurveys){
          y[i,j] ~ dbern(p[i,j])
          p[i,j] <- 1 - pow((1 - r[i,j]), N[i])
          logit(r[i,j]) <- alpha0 + alpha1[setup[i]]  
        }
      }
      
      #'  Derived parameters
      #'  ------------------
      #'  Mean lambda per GMU
      for(gmu in 1:ngmu) {
        lambdaGMU[gmu] <- exp(beta0 + beta1[gmu])
      }
      
      #'  Mean lambda averaged across GMUs
      mu.lambda <- mean(lambdaGMU[])
      
      #'  Predicted site-level abundance per GMU 
      for(i in 1:ncams1) {
        Ngmu1[i] <- exp(beta0 + beta1[gmu[1]])
      }
      for(i in 1:ncams2) {
        Ngmu2[i] <- exp(beta0 + beta1[gmu[2]])
      }
      
      #'  Total abundance across camera sites
      totalN <- sum(N[])
      
      #'  Mean density per GMU
      totalN.gmu10a <- sum(Ngmu1[])
      densitykm2.gmu10a <- totalN.gmu10a/area1
      density100km2.gmu10a <- densitykm2.gmu10a * 100

      totalN.gmu6 <- sum(Ngmu2[])
      densitykm2.gmu6 <- totalN.gmu6/area2
      density100km2.gmu6 <- densitykm2.gmu6 * 100

      #'  Total sites occupied (N > 0)
      for(i in 1:nsites) {
        occupied[i] <- ifelse(N[i] > 0, 1, 0)
      }
      occSites <- sum(occupied[])

      #'  Mean occupancy probability
      mean.psi <- 1 - exp(-mu.lambda)
 
      #'  Mean per-individual detection probability (r) per camera setup
      for(cam in 1:nsets) {
        rSetup[cam] <- 1/(1 + exp(-(alpha0 + alpha1[cam])))
      }
      #'  per-individual detection probability (r) averaged across all camera setups 
      mu.r <- mean(rSetup[])

      #'  Mean detection probability (p)
      for(i in 1:nsites) {
        p.occasion[i] <- mean(p[i,])
      }
      mean.p <- mean(p.occasion[])
      
      }
      ")
