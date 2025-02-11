###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####         HYBRID  ALGORITHM to generate choice and RTs samples
###############################################################################
##### This hybrid approach combines:
##### 1. Direct sampling from von Mises distribution for angular choices
##### 2. Metropolis-Hastings algorithm for response times
###############################################################################
########################################################   by Adriana F. Chávez 
library(circular)   # For von Mises distribution (rvonmises)

rCDDM_Hybrid <- function(n, par, plot=FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Parameter extraction and validation
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extract model parameters
  mu1 <- par$mu1;   mu2 <- par$mu2                 # Drift vector in Cartesian coordinates
  drift <- par$drift;   theta <- par$theta         # Drift direction in polar coordinates
  tzero <- par$tzero;   boundary <- par$boundary   # Non-decision time and boundary
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Defensive coding: ensure we have polar coordinates
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  noPolar <- is.null(theta) & is.null(drift)
  noRect <- is.null(mu1) & is.null(mu2)
  if(noPolar){
    if(noRect){ 
      stop("Provide Cartesian or Polar coordinates", call. = FALSE)
    }else{
      Mu <- rectToPolar(mu1, mu2)
      drift <- Mu$dLength;    par$drift <- drift
      theta <- Mu$dAngle;     par$theta <- theta
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Initialize sampling parameters
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get density characteristics from CDDM
  test.Density <- keyDensityPoints(par)
  min.RT <- test.Density$min.RT          # Minimum possible RT
  max.RT <- test.Density$max.RT          # Maximum considered RT
  max.Density <- test.Density$max.Density # Peak density
  mean.rt <- ezcddm_MRT(drift, boundary, tzero)  # Theoretical mean RT  
  # Initialize proposal parameters for RT
  mu_rt <- mean.rt    # ex-Gaussian mean (centered at theoretical mean)
  sigma_rt <- 0.5     # ex-Gaussian SD (will be tuned)
  tau_rt <- 0.5       # ex-Gaussian exponential component (will be tuned)
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Direct sampling of angular choices
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set von Mises concentration based on drift and boundary
  # Higher drift*boundary = more concentrated around the mean
  kappa <- drift * boundary  
  # Generate all angular choices directly from von Mises
  choices <- as.numeric(rvonmises(n, mu=circular(theta), kappa=kappa)) %% (2*pi)  
  # Initialize samples dataframe with choices (RTs will be filled later)
  samples <- data.frame(Choice = choices, RT = numeric(n))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Warmup phase: Tune RT proposal distribution
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Initialize M-H parameters
  M <- 100                # Number of warmup iterations
  ARate_des <- 0.4       # Desired acceptance rate (typical for M-H)
  ARate_obs <- ARate_des # Observed acceptance rate
  current_rt <- mean.rt   # Initial RT state
  
  # Begin warmup phase
  warm_up <- TRUE
  while(warm_up){
    # Generate M warmup samples
    u <- runif(M,0,1)  # Random numbers for acceptance decisions
    change <- rep(NA,M)  # Track which proposals are accepted
    
    for(i in 1:M){
      # Generate candidate RT from proposal distribution
      cand_rt <- rexGAUS(1, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
      cand_rt <- max(cand_rt, min.RT)  # Ensure RT exceeds minimum
      
      # Calculate acceptance ratio using first choice as representative
      # (any choice would work since we're just tuning the RT proposal)
      target_num <- dCDDM(c(choices[1], cand_rt), drift, theta, tzero, boundary)
      target_den <- dCDDM(c(choices[1], current_rt), drift, theta, tzero, boundary)
      
      prop_num <- dexGAUS(current_rt, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
      prop_den <- dexGAUS(cand_rt, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
      
      ratio <- (target_num * prop_num) / (target_den * prop_den)
      
      # Accept/reject
      change[i] <- ratio > u[i]
      if(change[i]){ current_rt <- cand_rt }
    }
    
    # Update acceptance rate and tune proposal parameters
    ARate_obs <- mean(change)
    if(ARate_obs < ARate_des){
      sigma_rt <- sigma_rt * 1.1  # Broaden ex-Gaussian
      tau_rt <- tau_rt * 1.1
    } else {
      sigma_rt <- sigma_rt * 0.9  # Narrow ex-Gaussian
      tau_rt <- tau_rt * 0.9
    }
    
    # Check if acceptance rate is within desired range
    warm_up <- ARate_obs < (ARate_des * 0.8) || ARate_obs > (ARate_des * 1.2)
    if(ARate_obs == 0){ ARate_obs <- ARate_des }  # Prevent division by zero
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Main sampling phase: Generate RTs for each choice
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  u <- runif(n,0,1)  # Random numbers for acceptance decisions
  
  for(i in 1:n){
    # Generate candidate RT
    cand_rt <- rexGAUS(1, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
    cand_rt <- max(cand_rt, min.RT)
    
    # Calculate acceptance ratio for current choice
    target_num <- dCDDM(c(samples$Choice[i], cand_rt), drift, theta, tzero, boundary)
    target_den <- dCDDM(c(samples$Choice[i], current_rt), drift, theta, tzero, boundary)
    
    prop_num <- dexGAUS(current_rt, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
    prop_den <- dexGAUS(cand_rt, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
    
    ratio <- (target_num * prop_num) / (target_den * prop_den)
    
    # Store sample
    if(ratio > u[i]) {
      current_rt <- cand_rt
      samples$RT[i] <- cand_rt
    } else {
      samples$RT[i] <- current_rt
    }
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Post-processing: Rotate choices back to original angle
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  samples$Choice <- (samples$Choice - rotation_angle) %% (2*pi)
  
  return(samples)
} 