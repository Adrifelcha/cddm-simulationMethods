###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                      REJECTION  ALGORITHM
###############################################################################
#####      This function simulates data from the CDDM using a rejection
#####      sampling algorithm with a bivariate normal proposal distribution 
#####      for the choice and response time.
###############################################################################
########################################################   by Adriana F. Chávez
library(mvtnorm)  # For multivariate normal distribution


sample.Reject.2DNormal.cddm <- function(n, par, plot=FALSE){
  # Extract parameters and validate
  mu1 <- par$mu1;       mu2 <- par$mu2        # Cartesian drift components
  drift <- par$drift;   theta <- par$theta     # Polar drift components
  tzero <- par$tzero;   boundary <- par$boundary # Non-decision time and boundary
  
  # Parameter validation: Check if we have either Cartesian or Polar coordinates
  noPolar <- is.null(theta) & is.null(drift)
  noRect <- is.null(mu1) & is.null(mu2)
  if(noPolar){
    if(noRect){
      stop("Provide Cartesian or Polar coordinates", call. = FALSE)
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
  
  # Set up bivariate normal proposal parameters
  proposal_mean <- c(theta, pred.RT)  # Center at predicted angle and RT
  proposal_sigma <- matrix(c(
    (2*pi/6)^2,     0,        # Variance for angle (covers ~1/6 of circle)
    0,        (max.RT/4)^2    # Variance for RT (covers ~1/4 of range)
  ), nrow=2)
  
  # Optional visualization setup
  if(plot){
    par(mar = c(0, 0, 0, 0))
    nLines <- 40
    x.C <- seq(0, 2*pi, length.out=nLines)
    y.RT <- seq(min.RT, max.RT, length.out=nLines)
    
    # Compute CDDM density for visualization
    z.Dens <- matrix(NA, nrow=nLines, ncol=nLines)
    for(c in 1:nLines){ for(t in 1:nLines){
      z.Dens[c,t] <- dCDDM(c(x.C[c], y.RT[t]), drift, theta, tzero, boundary)
    }}
    
    # Create 3D scatter plot
    a <- scatterplot3d(x.C, y.RT, diag(z.Dens), 
                      type="l", color="blue",
                      xlab="Choices", ylab="RT", zlab="Density",
                      zlim=c(0, max.Density))
  }
  
  # Initialize sampling variables
  n.keep <- 0
  samples <- data.frame(Choice = numeric(), RT = numeric())
  
  # Main rejection sampling loop
  while(n.keep < n){
    # Generate batch of candidates from bivariate normal
    cand <- rmvnorm(n, proposal_mean, proposal_sigma)
    
    # Wrap angles to [0, 2π] and filter valid RTs
    cand[,1] <- cand[,1] %% (2*pi)
    valid_rt <- cand[,2] >= min.RT
    if(sum(valid_rt) == 0) next
    
    valid_cand <- cand[valid_rt,]
    
    # Calculate proposal density (bivariate normal)
    prop_density <- dmvnorm(valid_cand, proposal_mean, proposal_sigma)
    
    # Calculate target CDDM density
    target_density <- dCDDM(valid_cand, drift, theta, tzero, boundary)
    
    # Generate rejection thresholds
    rej.crit <- runif(nrow(valid_cand), 0, max.Density * max(prop_density))
    
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
    
    n.keep <- nrow(samples)
  }
  
  # Return exactly n samples
  samples <- samples[1:n,]
  
  return(samples)
}
