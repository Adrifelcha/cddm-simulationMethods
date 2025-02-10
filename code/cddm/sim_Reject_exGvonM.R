###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                      REJECTION  ALGORITHM
###############################################################################
#####      This function simulates data from the CDDM using a rejection
#####      sampling algorithm with an ex-Gaussian proposal distribution for the
#####      response time and a von Mises distribution for the choice.
###############################################################################
########################################################   by Adriana F. Chavez 
library(circular)  # For von Mises distribution (rvm)

rCDDM_Reject_exGvonM <- function(n, par, plot=FALSE){
  # Extract parameters and validate
  mu1 <- par$mu1;       mu2 <- par$mu2        # Cartesian drift components
  drift <- par$drift;   theta <- par$theta     # Polar drift components
  tzero <- par$tzero;   boundary <- par$boundary # Non-decision time and boundary
  
  # Parameter validation: Check if we have either Cartesian or Polar coordinates
  noPolar <- is.null(theta) & is.null(drift)
  noRect <- is.null(mu1) & is.null(mu2)
  if(noPolar){
      if(noRect){  stop("Provide Cartesian or Polar coordinates", call. = FALSE)
      }else{
        # Convert Cartesian to polar coordinates if needed
        Mu <- rectToPolar(mu1, mu2)
        drift <- Mu$dLength;    par$drift <- drift  # Magnitude of drift
        theta <- Mu$dAngle;     par$theta <- theta  # Direction of drift
      }
  }
  
  # Get density characteristics for parameter scaling
  test.Density <- keyDensityPoints(par)
  min.RT <- test.Density$min.RT           # Minimum response time
  max.RT <- test.Density$max.RT           # Maximum response time
  pred.RT <- test.Density$pred.RT         # Most likely response time
  max.Density <- test.Density$max.Density*1.1 # Maximum density (with 10% buffer)
  
  # Initialize sampling variables
  n.keep <- 0  # Counter for accepted samples
  samples <- data.frame(Choice = numeric(), RT = numeric())  # Store accepted samples
  
  # Main rejection sampling loop
  while(n.keep < n){
    # Generate batch of candidates from proposal distributions
    kappa <- drift * boundary  # Concentration parameter for von Mises
    cand <- matrix(NA, nrow=n, ncol=2)
    cand[,1] <- rvonmises(n, theta, kappa)  # Sample choices from von Mises
    
    # Sample RTs from ex-Gaussian
    mu <- pred.RT
    sigma <- boundary/drift
    tau <- boundary/(2*drift)
    cand[,2] <- rexGAUS(n, mu, sigma, tau)
    
    # Filter out RTs below minimum
    valid_rt <- cand[,2] >= min.RT
    if(sum(valid_rt) == 0) next
    
    # Calculate densities for valid candidates
    valid_cand <- cand[valid_rt,]
    
    # Calculate proposal density
    prop_dens_choice <- dvonmises(valid_cand[,1], theta, kappa)
    prop_dens_rt <- dexGAUS(valid_cand[,2], mu, sigma, tau)
    prop_density <- prop_dens_choice * prop_dens_rt
    
    # Calculate target CDDM density
    target_density <- dCDDM(valid_cand, drift, theta, tzero, boundary)
    
    # Generate rejection thresholds
    rej.crit <- runif(nrow(valid_cand), 0, max.Density * prop_density)
    
    # Determine which candidates to keep
    keep <- (target_density >= rej.crit)
    
    # Visualize accepted and rejected samples if plotting
    if(plot){
      a$points3d(valid_cand[!keep,1], valid_cand[!keep,2], rej.crit[!keep],
                 col = rgb(204/255,0,0,0.2), pch = 16, cex = 0.3)  # Rejected (red)
      a$points3d(valid_cand[keep,1], valid_cand[keep,2], rej.crit[keep],
                 col = rgb(0,204/255,0,0.8), pch = 16, cex = 0.3)  # Accepted (green)
    }
    
    # Store accepted samples
    if(sum(keep) > 0) {
      new_samples <- data.frame(
        Choice = valid_cand[keep,1],
        RT = valid_cand[keep,2]
      )
      samples <- rbind(samples, new_samples)
    }
    
    n.keep <- nrow(samples)  # Update counter
  }
  
  # Return exactly n samples
  samples <- samples[1:n,]
  
  return(samples)
}