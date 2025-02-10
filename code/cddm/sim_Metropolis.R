###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                     Metropolis Sampling Algorithm
###############################################################################
########################################################   by Adriana F. Ch?vez 
#if(!exists("superCalled")){superCalled <- FALSE}
#if(!superCalled){ 
#  source("./dCDDM.R") 
#  source("./sim_auxiliarFunctions.R")
#}
library(mvtnorm)

# Function to generate samples from CDDM using Metropolis-Hastings algorithm
sample.Metropolis.cddm <- function(n, par, plot=FALSE){  # n=number of samples, par=parameters, plot=visualization flag
  # Extract model parameters from the parameter list
  mu1 <- par$mu1
  mu2 <- par$mu2
  drift <- par$drift      # Drift rate: strength/speed of information accumulation
  theta <- par$theta      # Mean direction: predicted/target angle in radians
  tzero <- par$tzero      # Non-decision time: time for non-decision components
  boundary <- par$boundary # Decision boundary: threshold for decision
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  #               Defensive Coding                                         #
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
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  
  # Initialize predicted choice angle from theta
  predChoice <- theta     # Predicted choice angle (initialized to mean direction)      
  
  # Get density characteristics for tuning the proposal distribution
  test.Density <- keyDensityPoints(par)  # Calculate key points of the density function
  min.RT <- test.Density$min.RT          # Minimum possible reaction time
  max.RT <- test.Density$max.RT          # Maximum considered reaction time
  max.Density <- test.Density$max.Density # Peak density value
  predRT <- test.Density$pred.RT         # Most likely reaction time
  
  if(plot){          # If plotting is enabled
    par(mar = c(0, 0, 0, 0))  # Set margins to zero
    # Create grid for density visualization
    nLines <- 40  # Number of grid lines
    x.C <- seq(0,2*pi,length.out=nLines)            # Grid points for choices
    y.RT <- seq(min.RT,max.RT,length.out=nLines)    # Grid points for RT
    x.theta <- rep(theta,nLines)                    # Vector of mean directions
    
    # Compute density values for visualization
    z.Dens <- matrix(NA, nrow=nLines, ncol=nLines)  # Initialize density matrix
    for(c in 1:nLines){ for(t in 1:nLines){  # Loop through grid
      z.Dens[c,t] <- dCDDM(c(x.C[c],y.RT[t]),drift,theta,tzero,boundary)  # Compute density
    }}
    # Compute density along mean direction
    theta.Dens <- dCDDM(cbind(x.theta,y.RT),drift,theta,tzero,boundary)
    
    # Create 3D scatter plot
    a <- scatterplot3d(x.C, y.RT, diag(z.Dens), zlab="Density", color = "white", type="l",
                       xlab="Choices", ylab="RT", zlim = c(0, max(z.Dens, theta.Dens)))
  }
  
  # Initialize Metropolis-Hastings algorithm parameters
  M <- 300                # Number of iterations for warmup phase
  ARate_des <- 0.4       # Target acceptance rate (40% is typical for M-H)
  ARate_obs <- ARate_des # Initialize observed acceptance rate
  Mu <- c(predChoice, predRT)  # Starting point: predicted choice and RT
  Sigma <- diag(c((2*pi)^2,50))  # Initial proposal covariance (choice variance, RT variance)
  
  # Begin warmup phase to tune proposal distribution
  warm_up <- TRUE  # Flag to control warmup loop
  nLines2 <- 15    # Number of grid lines for proposal distribution visualization
  x.C2 <- seq(0,2*pi,length.out=nLines2)            # Grid points for choices
  y.RT2 <- seq(min.RT,max.RT,length.out=nLines2)    # Grid points for RT
  
  while(warm_up){  # Continue until desired acceptance rate is achieved
    # Generate random RGB values for plotting
    b <- runif(1,0,1)   # Random blue component
    e <- runif(1,0,1)   # Random green component
    d <- runif(1,0,1)   # Random red component
    baseColor <- rgb(b,e,d,0.05)   # Create transparent color for density
    meanColor <- rgb(b,e,d,0.4)    # Create more opaque color for points
    
    # Adjust proposal covariance based on acceptance rate
    newSigma <- Sigma*((ARate_obs/ARate_des)^(0.5))  # Scale covariance matrix
    Sigma <- newSigma  # Update proposal covariance
    max.DensityMVN <- dmvnorm(Mu,Mu,Sigma)  # Maximum density of proposal distribution
    
    if(plot){  # Visualize proposal distribution
      # Compute proposal distribution density
      z.Dens2 <- matrix(NA, nrow=nLines2, ncol=nLines2)
      for(c in 1:nLines2){ for(t in 1:nLines2){
        z.Dens2[c,t] <- dmvnorm(c(x.C2[c],y.RT2[t]),mean=Mu, sigma = Sigma)
      }}
      z.Dens2 <- (z.Dens2*max.Density)/max.DensityMVN  # Scale density values
      
      # Plot proposal distribution
      for(i in 1:nLines2){  
        a$points3d(rep(x.C2[i],nLines2), y.RT2, z.Dens2[i,], type="l", col = baseColor)
        a$points3d(x.C2, rep(y.RT2[i],nLines2), z.Dens2[,i],  col = baseColor, type="l")
      }
      a$points3d(Mu[1],Mu[2],max.Density, col = meanColor, cex=0.5, pch=16)
    }
    
    # Generate M warmup samples
    u <- runif(M,0,1)      # Random uniforms for acceptance decisions
    change <- rep(NA,M)    # Track which proposals are accepted
    for(i in 1:M){  # Loop through warmup iterations
      # Generate candidate from proposal distribution
      cand <- rmvnorm(1,Mu,Sigma)  # Draw from bivariate normal
      cand[1] <- cand[1] %% (2*pi)  # Wrap choice angle to [0, 2π]
      
      # Calculate acceptance ratio
      ratio.num <- max(dCDDM(cand,drift, theta, tzero, boundary),0,na.rm = TRUE)  # Candidate density
      ratio.den <- dCDDM(Mu,drift, theta, tzero, boundary)  # Current state density
      ratio <- ratio.num/ratio.den  # Compute acceptance ratio
      
      # Accept/reject step
      change[i] <- ratio>u[i]  # Accept if ratio exceeds random uniform
      if(change[i]){  # If accepted
        Mu <- as.vector(cand)  # Update current state
        if(plot){  # Update visualization if plotting enabled
          z.Dens2 <- matrix(NA, nrow=nLines2, ncol=nLines2)
          for(c in 1:nLines2){ for(t in 1:nLines2){
            z.Dens2[c,t] <- dmvnorm(c(x.C2[c],y.RT2[t]),mean=Mu, sigma = Sigma)
          }}
          z.Dens2 <- (z.Dens2*max.Density)/max.DensityMVN
          for(i in 1:nLines2){  
            a$points3d(rep(x.C2[i],nLines2), y.RT2, z.Dens2[i,], type="l", col = baseColor)
            a$points3d(x.C2, rep(y.RT2[i],nLines2), z.Dens2[,i],  col = baseColor, type="l")
          }
          a$points3d(Mu[1],Mu[2],max.Density, col = meanColor, cex=0.5, pch=16)
        }
      }
    }
    
    # Update acceptance rate and check convergence
    ARate_obs <- mean(change)  # Calculate observed acceptance rate
    warm_up <- ARate_obs<ARate_des  # Continue if acceptance rate too low
    if(ARate_obs==0){ ARate_obs <- ARate_des}  # Prevent division by zero
  }
  
  if(plot){  # Plot final CDDM density
    baseColorCDDM <- rgb(0,0,0,0.5)  # Color for CDDM density visualization
    for(i in 1:nLines){  # Plot density contours
      a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l", col = baseColorCDDM)
      a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i],  col = baseColorCDDM, type="l")
    }
    a$points3d(x.theta, y.RT, theta.Dens, col = baseColorCDDM, type="l")
  }
  
  # Main sampling phase
  samples <- matrix(NA, nrow=n, ncol=2)  # Initialize matrix for samples
  u <- runif(n,0,1)  # Generate random uniforms for acceptance decisions
  for(i in 1:n){  # Loop through desired number of samples
    # Generate valid candidate (ensure RT > minimum)
    cand <- c(0,-1)  # Initialize invalid candidate
    while(cand[2]<=min.RT){  # Keep trying until valid RT
      cand <- rmvnorm(1,Mu,Sigma)  # Generate candidate
    }
    samples[i,1] <- cand[1] %% (2*pi)  # Store wrapped choice angle
    samples[i,2] <- cand[2]            # Store reaction time
    
    # Perform Metropolis acceptance step
    ratio.num <- max(dCDDM(cand,drift, theta, tzero, boundary),0,na.rm = TRUE)  # Candidate density
    ratio.den <- dCDDM(Mu,drift, theta, tzero, boundary)  # Current state density
    ratio <- ratio.num/ratio.den  # Compute acceptance ratio
    if(ratio>u[i]){  # Accept if ratio exceeds random uniform
      Mu <- as.vector(cand)  # Update current state
    }
  }
  
  if(plot){  # Plot final samples
    a$points3d(samples[,1],samples[,2],dCDDM(samples,drift,theta,tzero,boundary),
               col = rgb(b,e,d,0.5), pch = 16, cex = 0.5)
  }
  
  # Convert matrix to data frame and add column names
  samples <- as.data.frame(samples)
  colnames(samples) <- c("Choice", "RT")
  
  return(samples)  # Return data frame of samples
}