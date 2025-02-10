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

rCDDM_Reject_2DNormal <- function(n, par, plot=FALSE, storePlotToFile=FALSE){
  # Extract parameters and validate
  mu1 <- par$mu1;       mu2 <- par$mu2           # Cartesian drift components
  drift <- par$drift;   theta <- par$theta       # Polar drift components
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

  # Calculate rotation needed to center theta at pi
  rotation <- pi - theta
  rotated_theta <- pi  # This will be our sampling center
  
  # Get density characteristics for parameter scaling
  test_density <- keyDensityPoints(par)
  min_rt <- test_density$min.RT
  max_rt <- test_density$max.RT
  max_density <- (test_density$max.Density) * (1.1)

  # Get theoretical moments for RT distribution
  mean_rt <- ezcddm_MRT(drift, boundary, tzero)
  rt_var <- ezcddm_VRT(drift, boundary)

  # Set up bivariate normal proposal parameters
  proposal_mean <- c(rotated_theta, mean_rt)  # Center at pi and theoretical mean RT
  proposal_sigma <- matrix(c(
    (2*pi/6)^2,     0,        # Variance for angle (covers ~1/6 of circle)
    0,              rt_var     # Use theoretical variance for RT
  ), nrow=2)
  
  # Optional visualization setup
  if(plot|storePlotToFile){
    if(storePlotToFile){
      pdf(file = here("results", "show_Reject_2DNormal.pdf"), width = 10, height = 10)
    }
    # Set up grid for plotting
    nSupp <- 100; nLines <- 30
    support.C <- seq(0, 2*pi, length.out=nSupp)
    support.RT1 <- seq(min_rt, max_rt, length.out=nSupp)
    support.RT2 <- rev(support.RT1)
    support.theta <- rep(theta, nSupp)
    
    # Calculate density values for different slices
    z.diag1 <- dCDDM(cbind(support.C, support.RT1), drift, theta, tzero, boundary)
    z.diag2 <- dCDDM(cbind(support.C, support.RT2), drift, theta, tzero, boundary)
    z.RT_at_theta <- dCDDM(cbind(support.theta, support.RT1), drift, theta, tzero, boundary)
    
    # Create 3D scatter plot
    a <- scatterplot3d(support.C, support.RT1, z.diag1, 
                      xlim=c(0, 2*pi), ylim=c(min_rt, max_rt), 
                      zlim=c(0, max_density),
                      xlab="Choices", ylab="RT", zlab="Density",
                      color="blue", type="l", bg="white")
    
    # Add additional density curves
    a$points3d(support.C, support.RT2, z.diag2, col="blue", type="l")
    a$points3d(support.theta, support.RT1, z.RT_at_theta, col="blue", type="l")
    
    # Add grid lines
    L <- round(nSupp/nLines, 0)
    for(i in 1:nLines){
      choose.RT <- rep(support.RT1[i*L], nSupp)
      choose.C <- rep(support.C[i*L], nSupp)
      z.overRT <- dCDDM(cbind(choose.C, support.RT1), drift, theta, tzero, boundary)
      z.overC <- dCDDM(cbind(support.C, choose.RT), drift, theta, tzero, boundary)
      a$points3d(support.C, choose.RT, z.overC, col="blue", type="l")
      a$points3d(choose.C, support.RT1, z.overRT, col="blue", type="l")
    }
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
    valid_rt <- cand[,2] >= min_rt
    if(sum(valid_rt) == 0) next
    
    valid_cand <- cand[valid_rt,]
    
    # Rotate angles back to original theta center
    valid_cand[,1] <- (valid_cand[,1] - rotation) %% (2*pi)    
   
    # Calculate target CDDM density
    target_density <- dCDDM(valid_cand, drift, theta, tzero, boundary)
    
    # Generate rejection thresholds (using only max_density)
    rej.crit <- runif(nrow(valid_cand), 0, max_density)
    
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
  
  if(storePlotToFile){    
    dev.off()    
  }

  # Return exactly n samples
  samples <- samples[1:n,]
  
  return(samples)
}
