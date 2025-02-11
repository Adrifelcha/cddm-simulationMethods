###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                      REJECTION  ALGORITHM
###############################################################################
# Function to generate samples from CDDM using independent Metropolis algorithms
# for choices and RTs. This version:
# 1. Uses von Mises proposals for choices
# 2. Uses ex-Gaussian proposals for RTs
# 3. Samples each dimension independently
###############################################################################
########################################################   by Adriana F. Chávez 
library(circular)   # For von Mises distribution (rvonmises)

rCDDM_MetropolisHastings <- function(n, par, plot=FALSE){
  # Extract model parameters
  mu1 <- par$mu1;   mu2 <- par$mu2
  drift <- par$drift;   theta <- par$theta
  tzero <- par$tzero;   boundary <- par$boundary
  
  # Defensive coding
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
  
  # Store original theta and compute rotation angle to center at pi
  theta_original <- theta
  rotation_angle <- pi - theta
  theta <- pi  # Center the distribution at pi
  
  # Get density characteristics
  test.Density <- keyDensityPoints(par)
  min.RT <- test.Density$min.RT
  max.RT <- test.Density$max.RT
  max.Density <- test.Density$max.Density
  predRT <- ezcddm_MRT(drift, boundary, tzero)
  
  # Initialize proposal parameters
  kappa <- 0.5        # von Mises concentration
  mu_rt <- predRT     # ex-Gaussian mean
  sigma_rt <- 0.5     # ex-Gaussian SD
  tau_rt <- 0.5       # ex-Gaussian exponential component
  
  # Initialize Metropolis parameters
  M <- 300               # Warmup iterations
  ARate_des <- 0.4       # Target acceptance rate
  ARate_obs_choice <- ARate_des # Observed rates for each dimension
  ARate_obs_rt <- ARate_des
  current_choice <- pi    # Current states
  current_rt <- predRT
  
  # Begin warmup phase
  warm_up <- TRUE
  while(warm_up){
    # Generate M warmup samples
    u <- runif(M,0,1)
    change_choice <- rep(NA,M)
    change_rt <- rep(NA,M)
    
    for(i in 1:M){
      # Generate candidates independently
      cand_choice <- as.numeric(rvonmises(1, mu=circular(current_choice), kappa=kappa)) %% (2*pi)
      cand_rt <- rexGAUS(1, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
      cand_rt <- max(cand_rt, min.RT)
      
      # Calculate acceptance ratios separately for each dimension
      # For choice (keeping current RT)
      target_num_choice <- dCDDM(c(cand_choice, current_rt), drift, theta, tzero, boundary)
      target_den_choice <- dCDDM(c(current_choice, current_rt), drift, theta, tzero, boundary)
      prop_num_choice <- dvonmises(circular(current_choice), mu=circular(cand_choice), kappa=kappa)
      prop_den_choice <- dvonmises(circular(cand_choice), mu=circular(current_choice), kappa=kappa)
      ratio_choice <- (target_num_choice * prop_num_choice) / (target_den_choice * prop_den_choice)
      
      # For RT (keeping current choice)
      target_num_rt <- dCDDM(c(current_choice, cand_rt), drift, theta, tzero, boundary)
      target_den_rt <- dCDDM(c(current_choice, current_rt), drift, theta, tzero, boundary)
      prop_num_rt <- dexGAUS(current_rt, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
      prop_den_rt <- dexGAUS(cand_rt, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
      ratio_rt <- (target_num_rt * prop_num_rt) / (target_den_rt * prop_den_rt)
      
      # Accept/reject independently
      change_choice[i] <- ratio_choice > u[i]
      if(change_choice[i]){ current_choice <- cand_choice }
      
      change_rt[i] <- ratio_rt > u[i]
      if(change_rt[i]){ current_rt <- cand_rt }
    }
    
    # Update acceptance rates and tune parameters independently
    ARate_obs_choice <- mean(change_choice)
    ARate_obs_rt <- mean(change_rt)
    
    # Tune choice proposal
    if(ARate_obs_choice < ARate_des){
      kappa <- kappa * 0.9  # Broaden von Mises
    } else {
      kappa <- kappa * 1.1  # Narrow von Mises
    }
    
    # Tune RT proposal
    if(ARate_obs_rt < ARate_des){
      sigma_rt <- sigma_rt * 1.1  # Broaden ex-Gaussian
      tau_rt <- tau_rt * 1.1
    } else {
      sigma_rt <- sigma_rt * 0.9  # Narrow ex-Gaussian
      tau_rt <- tau_rt * 0.9
    }
    
    # Check convergence for both dimensions
    warm_up <- ARate_obs_choice < (ARate_des * 0.8) || ARate_obs_choice > (ARate_des * 1.2) ||
               ARate_obs_rt < (ARate_des * 0.8) || ARate_obs_rt > (ARate_des * 1.2)
    
    # Prevent division by zero
    if(ARate_obs_choice == 0){ ARate_obs_choice <- ARate_des }
    if(ARate_obs_rt == 0){ ARate_obs_rt <- ARate_des }
  }
  
  # Main sampling phase
  samples <- data.frame(Choice = numeric(n), RT = numeric(n))
  u <- runif(n,0,1)
  
  for(i in 1:n){
    # Generate candidates independently
    cand_choice <- as.numeric(rvonmises(1, mu=circular(current_choice), kappa=kappa)) %% (2*pi)
    cand_rt <- rexGAUS(1, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
    cand_rt <- max(cand_rt, min.RT)
    
    # Calculate acceptance ratios separately
    # For choice
    target_num_choice <- dCDDM(c(cand_choice, current_rt), drift, theta, tzero, boundary)
    target_den_choice <- dCDDM(c(current_choice, current_rt), drift, theta, tzero, boundary)
    prop_num_choice <- dvonmises(circular(current_choice), mu=circular(cand_choice), kappa=kappa)
    prop_den_choice <- dvonmises(circular(cand_choice), mu=circular(current_choice), kappa=kappa)
    ratio_choice <- (target_num_choice * prop_num_choice) / (target_den_choice * prop_den_choice)
    
    # For RT
    target_num_rt <- dCDDM(c(current_choice, cand_rt), drift, theta, tzero, boundary)
    target_den_rt <- dCDDM(c(current_choice, current_rt), drift, theta, tzero, boundary)
    prop_num_rt <- dexGAUS(current_rt, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
    prop_den_rt <- dexGAUS(cand_rt, mu=mu_rt, sigma=sigma_rt, tau=tau_rt)
    ratio_rt <- (target_num_rt * prop_num_rt) / (target_den_rt * prop_den_rt)
    
    # Accept/reject and store samples
    if(ratio_choice > u[i]){
      current_choice <- cand_choice
      samples$Choice[i] <- cand_choice
    } else {
      samples$Choice[i] <- current_choice
    }
    
    if(ratio_rt > u[i]){
      current_rt <- cand_rt
      samples$RT[i] <- cand_rt
    } else {
      samples$RT[i] <- current_rt
    }
  }
  
  # Rotate choices back to original theta
  samples$Choice <- (samples$Choice - rotation_angle) %% (2*pi)
  
  return(samples)
}

x <- rCDDM_MetropolisHastings(10, par=list(mu1=1, mu2=1, tzero=0, boundary=1))