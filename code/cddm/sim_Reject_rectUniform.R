###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                      REJECTION  ALGORITHM
###############################################################################
#####      This function simulates data from the CDDM using a rejection
#####      sampling algorithm with a uniform proposal distribution for the
#####      choice and response time.
###############################################################################
########################################################   by Adriana F. Chavez 


rCDDM_Reject_rectUniform <- function(n, par, plot=FALSE){
  # Extract parameters from input list
  mu1 <- par$mu1;       mu2 <- par$mu2        # Cartesian drift components
  drift <- par$drift;   theta <- par$theta     # Polar drift components
  tzero <- par$tzero;   boundary <- par$boundary # Non-decision time and boundary

  # Parameter Validation: Check if we have either Cartesian or Polar coordinates
  noPolar <- is.null(theta) & is.null(drift)   # Check if polar coordinates missing
  noRect <- is.null(mu1) & is.null(mu2)        # Check if Cartesian coordinates missing
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

  # Get key density points for setting up the sampling space
  test.Density <- keyDensityPoints(par)
  min.RT <- test.Density$min.RT           # Minimum response time
  max.RT <- test.Density$max.RT           # Maximum response time
  max.Density <- test.Density$max.Density*1.1  # Maximum density (with 10% buffer)
  base.C <- c(0, 2*pi)                    # Choice range (0 to 2π)
  base.RT <- c(min.RT, max.RT)            # RT range
  
  # Optional visualization of the density function
  if(plot){
    # Set up grid for plotting
    nSupp <- 100; nLines <- 30
    support.C <- seq(base.C[1],base.C[2],length.out=nSupp)
    support.RT1 <- seq(base.RT[1],base.RT[2],length.out=nSupp)
    support.RT2 <- rev(support.RT1)
    support.theta <- rep(theta,nSupp)
    
    # Calculate density values for different slices
    z.diag1 <- dCDDM(cbind(support.C,support.RT1),drift, theta, tzero, boundary)
    z.diag2 <- dCDDM(cbind(support.C,support.RT2),drift, theta, tzero, boundary)
    z.RT_at_theta <- dCDDM(cbind(support.theta,support.RT1),drift, theta, tzero, boundary)
    
    # Create 3D scatter plot
    a <- scatterplot3d(support.C, support.RT1, z.diag1, 
                       xlim=base.C, ylim=base.RT, zlim=c(0,max.Density),
                       xlab="Choices", ylab="RT", zlab="Density",
                       color="blue", type="l", bg="green")
    
    # Add additional density curves
    a$points3d(support.C, support.RT2, z.diag2, col = "blue", type="l")
    a$points3d(support.theta, support.RT1, z.RT_at_theta, col = "blue", type="l")
    
    # Add grid lines
    L <- round(nSupp/nLines,0)
    for(i in 1:nLines){
      choose.RT <- rep(support.RT1[i*L],nSupp)
      choose.C  <- rep(support.C[i*L],nSupp)
      z.overRT <- dCDDM(cbind(choose.C,support.RT1),drift, theta, tzero, boundary)
      z.overC <-  dCDDM(cbind(support.C,choose.RT),drift, theta, tzero, boundary)
      a$points3d(support.C, choose.RT, z.overC, col = "blue", type="l")
      a$points3d(choose.C, support.RT1, z.overRT, col = "blue", type="l")
    }
  }
  
  # Initialize sampling variables
  n.keep <- 0  # Counter for accepted samples
  samples <- data.frame(Choice = numeric(), RT = numeric())  # Store accepted samples
  
  # Main rejection sampling loop
  while(n.keep < n){
        # Generate candidate samples
        cand <- matrix(NA, nrow=n, ncol=2)
        cand[,1] <- runif(n,base.C[1],base.C[2])  # Sample choices uniformly
        cand[,2] <- runif(n,base.RT[1],base.RT[2])  # Sample RTs uniformly
        
        # Evaluate density at candidate points
        eval <- dCDDM(cand,drift, theta, tzero, boundary)
        rej.crit <- runif(n,0,max.Density)  # Generate rejection threshold
        keep <- (eval >= rej.crit)  # Accept samples above threshold
        
        # Visualize accepted and rejected samples if plotting
        if(plot){
          a$points3d(cand[!keep,1], cand[!keep,2], rej.crit[!keep],
                     col = rgb(204/255,0,0,0.2), pch = 16, cex = 0.3)  # Rejected (red)
          a$points3d(cand[keep,1], cand[keep,2], rej.crit[keep],
                     col = rgb(0,204/255,0,0.8), pch = 16, cex = 0.3)  # Accepted (green)
        }
        
        # Store accepted samples
        if(sum(keep) > 0) {
            new_samples <- data.frame(
                Choice = cand[keep,1],
                RT = cand[keep,2]
            )
            samples <- rbind(samples, new_samples)
        }
        
        n.keep <- n.keep + sum(keep)  # Update counter
  }
  
  # Return exactly n samples (may have generated more in last iteration)
  samples <- samples[1:n,]
  
  return(samples)
}