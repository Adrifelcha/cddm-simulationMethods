###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                      REJECTION  ALGORITHM
###############################################################################
#####      Our implementation considers three types of proposal distributions:
#####      1. Bivariate normal  (type="2DNormal")
#####      2. Ex-Gaussian and von Mises (type="exGvonM")
#####      3. Uniform (type="Uniform")
###############################################################################
########################################################   by Adriana F. Chávez 
library(mvtnorm)  # For multivariate normal distribution
library(circular)  # For von Mises distribution (rvm)

rCDDM_Reject <- function(n, par, type="2DNormal", plot=FALSE, createPDF=FALSE){  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Input validation
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Check if the specified proposal distribution type is valid
    if(!(type %in% c("2DNormal", "exGvonM", "Uniform"))){
        stop("Invalid type. Type must be one of: '2DNormal', 'exGvonM', or 'Uniform'.", call. = FALSE)
    }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Prepare parameter values
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~    
  # Extract parameter values in parameter list
  mu1 <- par$mu1;       mu2 <- par$mu2        # Cartesian drift components
  drift <- par$drift;   theta <- par$theta     # Polar drift components
  tzero <- par$tzero;   boundary <- par$boundary # Non-decision time and boundary  
  # Parameter validation: Check that drift vector is correctly specified
  noPolar <- is.null(theta) & is.null(drift)
  noRect <- is.null(mu1) & is.null(mu2)
  if(noPolar){
      if(noRect){  stop("Provide Cartesian or Polar coordinates", call. = FALSE)
      }else{
        # This algorithm requires polar coordinates
        Mu <- rectToPolar(mu1, mu2)
        drift <- Mu$dLength;    par$drift <- drift  # Magnitude of drift
        theta <- Mu$dAngle;     par$theta <- theta  # Direction of drift
      }
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Get key Choice and RT values based on the CDDM density
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  test_Density <- keyDensityPoints(par)
  min_RT <- test_Density$min.RT           # Minimum response time
  max_RT <- test_Density$max.RT           # Maximum response time
  pred_RT <- test_Density$pred.RT         # Most likely response time
  max_Density <- test_Density$max.Density*1.1 # Maximum density (with 10% buffer)  
# Get theoretical moments for RT distribution
  mean_rt <- ezcddm_MRT(drift, boundary, tzero)   # Expected RT
  rt_var <- ezcddm_VRT(drift, boundary)           # RT variance
  choice_var <- ezcddm_VCA(drift, boundary, tzero)   # Choice variance

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Set up parameters relevant for proposal distributions considered
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(type=="2DNormal"){
  rotation <- pi - theta
  rotated_theta <- pi
  proposal_mean <- c(rotated_theta, mean_rt)  # Center at pi and theoretical mean RT
  # Make sure proposal distribution is wider (use theoretical variances as reference)
  proposal_sigma <- matrix(c(max(c((2*pi/6)^2, choice_var),0),  
                             0,rt_var*2), nrow=2)
  }else if(type=="exGvonM"){
    kappa <- drift * boundary  # Concentration parameter for von Mises
    mu <- mean_rt
    sigma <- boundary/drift
    tau <- boundary/(2*drift)            
  }else if(type=="Uniform"){
    base.C <- c(0, 2*pi)
    base.RT <- c(min_RT, max_RT)
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Optional visualization setup
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  if(plot|createPDF){
      if(createPDF){
          if(type=="2DNormal"){
              pdf(file = here("results", "show_Reject_2DNormal.pdf"), width = 10, height = 10)
          }else if(type=="exGvonM"){
              pdf(file = here("results", "show_Reject_exGvonM.pdf"), width = 10, height = 10)
          }else if(type=="Uniform"){
              pdf(file = here("results", "show_Reject_Uniform.pdf"), width = 10, height = 10)
          }else{
              stop("Invalid type", call. = FALSE)
          }
      }
      # Set up grid for plotting
      nSupp <- 100; nLines <- 30
      support.C <- seq(0, 2*pi, length.out=nSupp)
      support.RT1 <- seq(min_RT, max_RT, length.out=nSupp)
      support.RT2 <- rev(support.RT1)
      support.theta <- rep(theta, nSupp)
      
      # Calculate density values for different slices
      z.diag1 <- dCDDM(cbind(support.C, support.RT1), drift, theta, tzero, boundary)
      z.diag2 <- dCDDM(cbind(support.C, support.RT2), drift, theta, tzero, boundary)
      z.RT_at_theta <- dCDDM(cbind(support.theta, support.RT1), drift, theta, tzero, boundary)
      
      # Create 3D scatter plot
      a <- scatterplot3d(support.C, support.RT1, z.diag1, 
                        xlim=c(0, 2*pi), ylim=c(min_RT, max_RT), 
                        zlim=c(0, max_Density),
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

    if(type=="exGvonM"){
        # Generate batch of candidates from proposal distributions        
        cand <- matrix(NA, nrow=n, ncol=2)
        cand[,1] <- rvonmises(n, theta, kappa)  # Sample choices from von Mises
        
        # Sample RTs from ex-Gaussian
        cand[,2] <- rexGAUS(n, mu, sigma, tau)
        
        # Filter out RTs below minimum
        valid_rt <- cand[,2] >= min_RT
        if(sum(valid_rt) == 0) next
        
        # Calculate densities for valid candidates
        valid_cand <- cand[valid_rt,]
    }else if(type=="2DNormal"){
        # Generate batch of candidates from bivariate normal
        cand <- rmvnorm(n, proposal_mean, proposal_sigma)
        
        # Wrap angles to [0, 2π] and filter valid RTs
        cand[,1] <- cand[,1] %% (2*pi)
        valid_rt <- cand[,2] >= min_RT
        if(sum(valid_rt) == 0) next
        
        valid_cand <- cand[valid_rt,]
        
        # Rotate angles back to original theta center
        valid_cand[,1] <- (valid_cand[,1] - rotation) %% (2*pi)        
    }else if(type=="Uniform"){
        # Generate candidate samples
        cand <- matrix(NA, nrow=n, ncol=2)
        cand[,1] <- runif(n,base.C[1],base.C[2])  # Sample choices uniformly
        cand[,2] <- runif(n,base.RT[1],base.RT[2])  # Sample RTs uniformly
        valid_cand <- cand
    }
  
    # Calculate target CDDM density
    target_density <- dCDDM(valid_cand, drift, theta, tzero, boundary)
    
    # Generate rejection thresholds
    rej.crit <- runif(nrow(valid_cand), 0, max_Density)
    
    # Determine which candidates to keep
    keep <- (target_density >= rej.crit)
    
    # Visualize accepted and rejected samples if plotting
    if(plot|createPDF){
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
  
  if(createPDF){    
      dev.off()    
  }
  
  # Return exactly n samples
  samples <- samples[1:n,]
  
  return(samples)
}
