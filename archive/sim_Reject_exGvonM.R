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

rCDDM_Reject_exGvonM <- function(n, par, plot=FALSE, storePlotToFile=FALSE){
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
  
  # Optional visualization setup
  if(plot|storePlotToFile){
      if(storePlotToFile){
          pdf(file = here("results", "show_Reject_exGvonM.pdf"), width = 10, height = 10)
      }
      # Set up grid for plotting
      nSupp <- 100; nLines <- 30
      support.C <- seq(0, 2*pi, length.out=nSupp)
      support.RT1 <- seq(min.RT, max.RT, length.out=nSupp)
      support.RT2 <- rev(support.RT1)
      support.theta <- rep(theta, nSupp)
      
      # Calculate density values for different slices
      z.diag1 <- dCDDM(cbind(support.C, support.RT1), drift, theta, tzero, boundary)
      z.diag2 <- dCDDM(cbind(support.C, support.RT2), drift, theta, tzero, boundary)
      z.RT_at_theta <- dCDDM(cbind(support.theta, support.RT1), drift, theta, tzero, boundary)
      
      # Create 3D scatter plot
      a <- scatterplot3d(support.C, support.RT1, z.diag1, 
                        xlim=c(0, 2*pi), ylim=c(min.RT, max.RT), 
                        zlim=c(0, max.Density),
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
  
    # Calculate target CDDM density
    target_density <- dCDDM(valid_cand, drift, theta, tzero, boundary)
    
    # Generate rejection thresholds
    rej.crit <- runif(nrow(valid_cand), 0, max.Density)
    
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
  
  if(storePlotToFile){    
      dev.off()    
  }
  
  # Return exactly n samples
  samples <- samples[1:n,]
  
  return(samples)
}

#rCDDM_Reject_exGvonM(n=20000, par=list(mu1=1, mu2=1, tzero=0, boundary=1), plot=TRUE, storePlotToFile=TRUE)