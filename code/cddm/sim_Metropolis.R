###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                     Metropolis Sampling Algorithm
###############################################################################
########################################################   by Adriana F. Chávez
library(mvtnorm)


# Helper function to generate valid candidate
generate_valid_candidate <- function(current, Sigma, logRT, min.RT, max.RT) {
  repeat {
          cand <- rmvnorm(1, mean=current, sigma=Sigma)
          cand[1] <- cand[1] %% (2*pi)
          if(logRT) {
                      cand[2] <- exp(cand[2])
          }
          if(cand[2] > min.RT && cand[2] < max.RT) {
              return(cand)
          }
  }
}


# Function to generate samples from CDDM using Metropolis algorithm
rCDDM_Metropolis <- function(n, par, plot = FALSE, logRT = FALSE){
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extract model parameters from the parameter list
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  mu1 <- par$mu1;       mu2 <- par$mu2           # Drift vector in Cartesian coordinates
  drift <- par$drift;   theta <- par$theta       # Drift direction in polar coordinates
  tzero <- par$tzero;   boundary <- par$boundary # Decision boundary and non-decision time

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #               Defensive Coding
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  noPolar <- is.null(theta) & is.null(drift)
  noRect <- is.null(mu1) & is.null(mu2)
  if(noPolar){
    if(noRect){ stop("Provide Cartesian or Polar coordinates", call. = FALSE)
    }else{
      # Convert Cartesian to polar coordinates if needed
      Mu <- rectToPolar(mu1, mu2)
      drift <- Mu$dLength;    par$drift <- drift  # Magnitude of drift
      theta <- Mu$dAngle;     par$theta <- theta  # Direction of drift
    }
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Store original theta and compute rotation angle to center at pi
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  theta_original <- theta
  rotation_angle <- pi - theta
  theta <- pi  # Center the distribution at pi  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get density characteristics for tuning the proposal distribution
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  test.Density <- keyDensityPoints(par)  # Calculate key density points
  min.RT <- test.Density$min.RT          # Minimum possible reaction time
  max.RT <- test.Density$max.RT          # Maximum considered reaction time
  max.Density <- test.Density$max.Density # Peak density value
  if(logRT){
              # Log-transform RTs and density
              log_min.RT <- log(min.RT)
              log_max.RT <- log(max.RT)
              log_max.Density <- log(max.Density)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define predicted choice and RT to center the proposal distribution
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  predRT <- ezcddm_MRT(drift, boundary, tzero)  # Expected RT
  if(logRT){
      # Log-transform RTs
      predRT <- log(predRT)
  }
  predChoice <- theta                           # Expected choice
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Plot the proposal distribution
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(plot){
            par(mar = c(0, 0, 0, 0))  # Set margins to zero
            # Create grid for density visualization
            nLines <- 40  # Number of grid lines
            x.C <- seq(0,2*pi,length.out=nLines)            # Grid points for choices
            y.RT <- seq(min.RT,max.RT,length.out=nLines)    # Grid points for RT
            x.theta <- rep(theta,nLines)                    # Vector of mean directions

            # Compute density values for visualization
            z.Dens <- matrix(NA, nrow=nLines, ncol=nLines)  # Initialize density matrix
            for(c in 1:nLines){ 
                for(t in 1:nLines){  # Loop through grid
                z.Dens[c,t] <- dCDDM(c(x.C[c],y.RT[t]),drift,theta,tzero,boundary)  # Compute density
                }
            }
            # Compute density along mean direction
            theta.Dens <- dCDDM(cbind(x.theta,y.RT),drift,theta,tzero,boundary)
            
            # Initialize 3D scatter plot - Blank plot using diag(z.Dens) as reference
            a <- scatterplot3d(x.C, y.RT, diag(z.Dens), zlab="Density", color = "white", type="l",
                              xlab="Choices", ylab="RT", zlim = c(0, max(z.Dens, theta.Dens)))
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Initialize Metropolis-Hastings algorithm parameters
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  M <- 300               # Number of iterations for warmup phase
  ARate_des <- 0.4       # Target acceptance rate (40% is typical for M-H)
  ARate_obs <- ARate_des # Default observed acceptance rate (Scale ARate_obs/ARate_des is 1)
  Mu <- c(predChoice, predRT)    # Starting point: predicted choice and RT
  
  # Smarter initial proposal variances based on parameter values
  if(logRT) {
    var_RT <- (log(boundary))^2/4  # Scale with log boundary
  } else {
    var_RT <- (boundary/2)^2       # Scale with boundary
  }
  var_choice <- min((2*pi)^2/4, (drift^2 + 1))  # Scale with drift, capped at π²
  
  # Initialize covariance matrix
  Sigma <- diag(c(var_choice, var_RT))
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Begin warmup phase to tune proposal distribution
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  warm_up <- TRUE  # Flag to control warmup loop
  if(plot){
            nLines2 <- 15    # Number of grid lines for proposal distribution visualization
            x.C2 <- seq(0,2*pi,length.out=nLines2)            # Grid points for choices
            y.RT2 <- seq(min.RT,max.RT,length.out=nLines2)    # Grid points for RT    
  }
  
  # Warmup continues until the desired acceptance rate is achieved
  while(warm_up){
    #>>>>>>>>>>>>>> Adjust proposal covariance based on acceptance rate
    scalingFactor <- ((ARate_obs/ARate_des)^(0.5))
    newSigma <- Sigma*scalingFactor  # Scale covariance matrix
    Sigma <- newSigma  # Update proposal covariance
    max.DensityMVN <- dmvnorm(Mu,Mu,Sigma)  # Maximum density of proposal distribution
    

    #>>>>>>>>>>>>>> Plot the proposal distribution
    if(plot){  
            # Generate random RGB values for plotting
            red <- runif(1,0,1);   green <- runif(1,0,1);   blue <- runif(1,0,1)   
            baseColor <- rgb(red,green,blue,0.05)   # Create transparent color for density
            meanColor <- rgb(red,green,blue,0.4)    # Create more opaque color for points

            # Compute proposal distribution density
            z.Dens2 <- matrix(NA, nrow=nLines2, ncol=nLines2)
            for(c in 1:nLines2){ 
                for(t in 1:nLines2){
                    z.Dens2[c,t] <- dmvnorm(c(x.C2[c],y.RT2[t]),mean=Mu, sigma = Sigma)
                }
            }
            z.Dens2 <- (z.Dens2*max.Density)/max.DensityMVN  # Scale density values
            
            # Plot proposal distribution
            for(i in 1:nLines2){  
              a$points3d(rep(x.C2[i],nLines2), y.RT2, z.Dens2[i,], type="l", col = baseColor)
              a$points3d(x.C2, rep(y.RT2[i],nLines2), z.Dens2[,i],  col = baseColor, type="l")
            }
            a$points3d(Mu[1],Mu[2],max.Density, col = meanColor, cex=0.5, pch=16)
    }
    
    #>>>>>>>>>>>>>> Generate samples from proposal distribution
    u <- runif(M,0,1)      # Generate M uniform criteria
    change <- rep(NA,M)    # Variable tracking which candidates are accepted
    for(i in 1:M){  
          cand <- rmvnorm(1, Mu, Sigma)  # Generate candidates
          cand[1] <- cand[1] %% (2*pi)   # Wrap choice angle to [0, 2π]
          
          # Handle RT bounds and transformations
          if(logRT){
              cand[2] <- exp(cand[2])    # Log-transform RTs if needed
          }
          
          # Only proceed with valid RTs
          if(cand[2] > min.RT && cand[2] < max.RT) {
              # Compute acceptance ratio
              ratio.num <- max(dCDDM(cand, drift, theta, tzero, boundary), 0, na.rm=TRUE)
              ratio.den <- dCDDM(Mu, drift, theta, tzero, boundary)
              ratio <- ratio.num/ratio.den  
              
              # Record acceptance/rejection
              change[i] <- ratio > u[i]
              
              # Update if accepted
              if(change[i]){  
                  if(logRT){
                      cand[2] <- log(cand[2])  # Transform back for proposal distribution
                  }
                  Mu <- as.vector(cand)
                  
                  # Plot updates if needed
                  if(plot){ 
                      z.Dens2 <- matrix(NA, nrow=nLines2, ncol=nLines2)
                      for(c in 1:nLines2){ 
                          for(t in 1:nLines2){
                              z.Dens2[c,t] <- dmvnorm(c(x.C2[c],y.RT2[t]),mean=Mu, sigma = Sigma)
                          }
                      }
                      z.Dens2 <- (z.Dens2*max.Density)/max.DensityMVN
                      for(i in 1:nLines2){  
                        a$points3d(rep(x.C2[i],nLines2), y.RT2, z.Dens2[i,], type="l", col = baseColor)
                        a$points3d(x.C2, rep(y.RT2[i],nLines2), z.Dens2[,i],  col = baseColor, type="l")
                      }
                      a$points3d(Mu[1],Mu[2],max.Density, col = meanColor, cex=0.5, pch=16)
                  }
              }
          } else {
              change[i] <- FALSE  # Record rejection for invalid RT
          }
    }
    
    # Evaluate the acceptance rate
    ARate_obs <- mean(change)       # Calculate observed acceptance rate
    warm_up <- ARate_obs<ARate_des  # Determine if warmup is complete
    if(ARate_obs==0){ ARate_obs <- ARate_des}  # Prevent division by zero
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Add the prescribed CDDM density plot to the 3D scatter plot
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(plot){  
            # Color for CDDM density visualization
            baseColorCDDM <- rgb(0,0,0,0.5)  
            for(i in 1:nLines){  # Plot density contours
                a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l", col = baseColorCDDM)
                a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i],  col = baseColorCDDM, type="l")
            }
            a$points3d(x.theta, y.RT, theta.Dens, col = baseColorCDDM, type="l")
  }  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Main sampling phase
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # Initialize sampling
  samples <- matrix(NA, nrow=n, ncol=2)
  u <- runif(n,0,1)
  current <- Mu
  
  # Main sampling loop
  for(i in 1:n){
      repeat{
              # Generate valid candidate
              cand <- generate_valid_candidate(current, Sigma, logRT, min.RT, max.RT)
              
              # Compute acceptance ratio
              ratio.num <- max(dCDDM(cand, drift, theta, tzero, boundary), 0, na.rm=TRUE)
              ratio.den <- dCDDM(current, drift, theta, tzero, boundary)
              ratio <- ratio.num/ratio.den
              
              # Accept if ratio exceeds rejection criterion
              if(ratio > u[i]) {
                current <- cand
                break
              }
      }

      # Store current state
      samples[i,] <- current
      if(logRT){
        samples[i,2] <- exp(current[2])
      }
  }
  
  # Rotate samples back to original theta
  samples[,1] <- (samples[,1] - rotation_angle) %% (2*pi)
  
  if(plot){   # Plot final samples
              # Rotate samples temporarily back for plotting
              plot_samples <- samples
              plot_samples[,1] <- (plot_samples[,1] + rotation_angle) %% (2*pi)
              a$points3d(plot_samples[,1], plot_samples[,2],
                        dCDDM(as.matrix(plot_samples), drift, pi, tzero, boundary),
                        col = rgb(b,e,d,0.5), pch = 16, cex = 0.5)
  }
  
  # Convert matrix to data frame and add column names
  samples <- as.data.frame(samples)
  colnames(samples) <- c("Choice", "RT")
  
  return(samples)  # Return data frame of samples
}


#x  <- rCDDM_Metropolis(1000, par = list(mu1 = 0.5,  mu2 = 0.5, boundary = 5, tzero = 0.1), logRT = TRUE)