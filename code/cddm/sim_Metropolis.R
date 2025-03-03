###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                     Metropolis Sampling Algorithm
###############################################################################
########################################################   by Adriana F. Chávez
library(mvtnorm)

###############################################################################
# Helper function for the sampling phase (ensure that the candidate is valid)
###############################################################################
generate_validMVN_candidate <- function(current, Sigma, logRT, tzero, max.RT) {
  repeat {
          # Get bivariate normal candidate
          cand <- rmvnorm(1, mean=current, sigma=Sigma)
          # Wrap choice angle to [0, 2π]
          cand[1] <- cand[1] %% (2*pi)
          # Log-transform RTs if needed
          if(logRT) {
               cand[2] <- exp(cand[2])               
          }
          # Ensure candidate is within RT bounds
          if(cand[2] > tzero && cand[2] < max.RT) {
              return(cand)
          }          
  }
}

###############################################################################
# MAIN FUNCTION: Generate samples from CDDM using Metropolis algorithm
###############################################################################
rCDDM_Metropolis <- function(n, par, plot = FALSE, logRT = FALSE, plot_warmup = FALSE, 
                             n_chains = 1, debug = FALSE){
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Prepare model parameters
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Extract model parameters from the parameter list
  mu1 <- par$mu1;       mu2 <- par$mu2           # Drift vector in Cartesian coordinates
  drift <- par$drift;   theta <- par$theta       # Drift direction in polar coordinates
  tzero <- par$tzero;   boundary <- par$boundary # Decision boundary and non-decision time

  # Ensure we have drift and theta  
  noPolar <- is.null(theta) & is.null(drift)
  noRect <- is.null(mu1) & is.null(mu2)
  if(noPolar){
    # Defensive coding: Ensure we have either polar or Cartesian coordinates
    if(noRect){ stop("Provide Cartesian or Polar coordinates", call. = FALSE)
    }else{
      # Convert Cartesian to polar coordinates if needed
      Mu <- rectToPolar(mu1, mu2)
      drift <- Mu$dLength;    par$drift <- drift  # Magnitude of drift
      theta <- Mu$dAngle;     par$theta <- theta  # Direction of drift
    }
  }
  
  # Get density characteristics for tuning the proposal distribution
  test.Density <- keyDensityPoints(par)  # Calculate key density points
  min.RT <- test.Density$min.RT          # Minimum possible reaction time
  max.RT <- test.Density$max.RT          # Maximum considered reaction time
  max.Density <- test.Density$max.Density # Peak density value  

  # Center theta at pi  
  theta_original <- theta
  rotation_angle <- pi - theta
  theta <- pi  # Center the distribution at pi  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Prepare proposal distribution
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Define predicted choice and RT to center the proposal distribution  
  predRT <- ezcddm_MRT(drift, boundary, tzero)  # Expected RT  
  predChoice <- theta                           # Expected choice

  # Determine reasonable variance of choice and RT
  predVarChoice <- ezcddm_VCA(drift, boundary)
  predVarRT <- ezcddm_VRT(drift, boundary)

  # If logRT is TRUE, transform RT to log-space
  if(logRT){   
                logRT_stats <- logRT_stats(predRT, predVarRT)
                predVarRT <- logRT_stats$sigma_logRT 
                predRT <- logRT_stats$mu_logRT
  }

  # Increase variance of choice and RT
  increase_var <- 4
  var_choice <- predVarChoice * increase_var
  var_RT <- predVarRT * increase_var

  # Initial proposal distribution parameters
  Sigma <- diag(c(var_choice,var_RT)) # Initial covariance matrix
  Mu <- c(predChoice, predRT)         # Initial mean vector

  # Start plotting space
  if(plot){
            # Determine best layout for plotting
            if(plot_warmup) {     # IF we are also plotting the warmup phase
                    # Create a layout matrix for 3x2 arrangement
                    # Value 1 is for main plot spanning both columns in top row
                    # Values 2-5 are for warmup plots in bottom two rows
                    layout_mat <- matrix(c(1,1,    # top row: main plot
                                          2,3,      # middle row: first two warmup plots
                                          4,5),     # bottom row: last two warmup plots
                                        nrow=3, ncol=2, byrow=TRUE)
                    layout(layout_mat, heights=c(1,1,1))  # Make the top row slightly taller
            } else {                # If we are not plotting the warmup phase
              par(mfrow=c(1,1), mar = c(0, 0, 0, 0))  # Set margins to zero
            }

            # Create grid for density visualization
            nLines <- 40  # Number of grid lines
            x.C <- seq(0,2*pi,length.out=nLines)            # Grid points for choices
            y.RT <- seq(min.RT,max.RT,length.out=nLines)    # Grid points for RT
            x.theta <- rep(theta,nLines)                    # Vector of mean directions

            # Compute CDDM density values to be used for plotting
            z.Dens <- matrix(NA, nrow=nLines, ncol=nLines)  # Initialize density matrix
            for(c in 1:nLines){     # Loop through grid points for choices
                for(t in 1:nLines){  # Loop through grid points for RT
                    z.Dens[c,t] <- dCDDM(c(x.C[c],y.RT[t]),drift,theta,tzero,boundary)
                }
            }            
            theta.Dens <- dCDDM(cbind(x.theta,y.RT),drift,theta,tzero,boundary)
            
            # Initialize empty 3D scatter plot 
            yrange <- c(min.RT,max.RT)  # RT range
            if(logRT){ yrange <- c(log(min.RT),log(max.RT))  # Log-transformed RT range
                       y.RT <- log(y.RT) }                   # Log-transformed RT values
            a <- scatterplot3d(x.C, y.RT, rep(0,nLines), zlab="Density", color = "white", type="l",
                               xlab="Choices", ylab="RT", zlim = c(0, max(z.Dens, theta.Dens)),
                               xlim = c(0,2*pi), ylim = yrange)

            # Initialize variables for proposal distribution visualization (2DNormal)
            nLines2 <- 15    # Number of grid lines for proposal distribution visualization
            x.C2 <- seq(0,2*pi,length.out=nLines2)            # Grid points for choices
            y.RT2 <- seq(min.RT,max.RT,length.out=nLines2)    # Grid points for RT    
            if(logRT){ y.RT2 <- log(y.RT2) }                  # Log-transformed RT values
            
            # Empty variables to store the colors used for each proposal distribution
            red <- NULL; green <- NULL; blue <- NULL
  }
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Initialize Metropolis algorithm: Warmup phase
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  M <- 400            # Number of iterations for warmup phase
  ARate_des <- 0.4    # Target acceptance rate (40% is typical for M-H)
    
   # Initialize warmup tracking
  if(plot_warmup | n_chains > 1) {
        # Initialize as list of n_chains empty data frames
        warmup_history <- vector("list", n_chains)
        for(i in 1:n_chains) {
              warmup_history[[i]] <- data.frame(iteration = integer(),
                                                sigma_choice = numeric(),
                                                sigma_rt = numeric(),
                                                acceptance_rate = numeric())
        }
  }
  
  # Initialize chain states with different starting points for each chain
  chain_states <- list()
  jitter <- runif(n_chains, 0.5, 1.5)
  for(i in 1:n_chains) {
        # For chain 1, use original starting point
        if (i == 1) {
            chain_states[[i]] <- list(mu = Mu, sigma = Sigma, ARate_obs = ARate_des)
        } else {
            # For other chains, perturb the initial proposal covariance
            jitter.sigma <- Sigma * matrix(c(jitter[i], 0,  # Scale choice variance
                                           0, jitter[i]),  # Scale RT variance
                                           nrow = 2, ncol = 2)
            chain_states[[i]] <- list(mu = Mu, sigma = jitter.sigma, ARate_obs = ARate_des)
        }
  }

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # B E G I N   W A R M U P   P H A S E   T O   T U N E   P R O P O S A L   D I S T R I B U T I O N
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  for(chain in 1:n_chains) {  # For each chain....
      # Load initial state
      current_state <- chain_states[[chain]]
      # Load initial variables for warmup phase
      warm_up <- TRUE   
      iteration <- 0
      
      # Start chain-specific warmup loop
      while(warm_up) { # Continue until this chain is warmed up
            # Update iteration counter
            iteration <- iteration + 1
            
            #>>>>>>>>>>>>>> Check if chain has achieved desired acceptance rate
            if(current_state$ARate_obs > ARate_des) {
                warm_up <- FALSE
                next
            }

            #>>>>>>>>>>>>>> Adjust proposal covariance based on acceptance rate
            scalingFactor <- ((current_state$ARate_obs/ARate_des)^(0.5))
            newSigma <- current_state$sigma * scalingFactor  # Scale covariance matrix
            current_state$sigma <- newSigma  # Update proposal covariance
            
            #>>>>>>>>>>>>>> Plot the proposal distribution
            if(plot && chain == n_chains){  # Only plot for the last chain
                    # Generate random RGB values for plotting
                    red <- runif(1,0,1); green <- runif(1,0,1); blue <- runif(1,0,1)   
                    baseColor <- rgb(red, green, blue, 0.05)
                    meanColor <- rgb(red, green, blue, 0.4)
                    
                    # Compute proposal distribution density
                    z.Dens2 <- matrix(NA, nrow=nLines2, ncol=nLines2)
                    for(c in 1:nLines2){ 
                        for(t in 1:nLines2){
                            z.Dens2[c,t] <- dmvnorm(c(x.C2[c],y.RT2[t]),mean=current_state$mu, sigma = current_state$sigma)
                        }
                    }
                    z.Dens2 <- (z.Dens2*max.Density)/dmvnorm(current_state$mu, current_state$mu, current_state$sigma)  # Scale density values
                    
                    # Plot proposal distribution
                    for(i in 1:nLines2){  
                        a$points3d(rep(x.C2[i],nLines2), y.RT2, z.Dens2[i,], type="l", col = baseColor)
                        a$points3d(x.C2, rep(y.RT2[i],nLines2), z.Dens2[,i],  col = baseColor, type="l")
                    }
                    a$points3d(current_state$mu[1],current_state$mu[2],max.Density, col = meanColor, cex=0.5, pch=16)
            }
            
            #>>>>>>>>>>>>>> Generate samples from proposal distribution
            u <- runif(M,0,1)      # Generate M uniform criteria
            change <- rep(NA,M)    # Variable tracking which candidates are accepted
            for(i in 1:M){  
                  cand <- rmvnorm(1, current_state$mu, current_state$sigma)  # Generate candidates
                  cand[1] <- cand[1] %% (2*pi)          # Wrap choice angle to [0, 2π]
                  if(logRT){ cand[2] <- exp(cand[2]) }  # Log-transform RTs if needed
                  
                  # Evaluate valid RTs
                  if(cand[2] > min.RT && cand[2] < max.RT) {
                      # Compute acceptance ratio
                      ratio.num <- max(dCDDM(cand, drift, theta, tzero, boundary), 0, na.rm=TRUE)
                      ratio.den <- dCDDM(current_state$mu, drift, theta, tzero, boundary)
                      ratio <- ratio.num/ratio.den  
                      
                      # Record acceptance/rejection
                      change[i] <- ratio > u[i]
                      
                      # Update if accepted
                      if(change[i]){  
                          if(logRT){    cand[2] <- log(cand[2])  }   # Log-transform RTs if needed
                          current_state$mu <- as.vector(cand)
                          
                          # Plot updated proposal distribution (only for the last chain)
                          if(plot && chain == n_chains){ 
                                  z.Dens2 <- matrix(NA, nrow=nLines2, ncol=nLines2)
                                  for(c in 1:nLines2){ 
                                      for(t in 1:nLines2){
                                          z.Dens2[c,t] <- dmvnorm(c(x.C2[c],y.RT2[t]),mean=current_state$mu, sigma = current_state$sigma)
                                      }
                                  }
                                  z.Dens2 <- (z.Dens2*max.Density)/dmvnorm(current_state$mu, current_state$mu, current_state$sigma)
                                  for(i in 1:nLines2){  
                                      a$points3d(rep(x.C2[i],nLines2), y.RT2, z.Dens2[i,], type="l", col = baseColor)
                                      a$points3d(x.C2, rep(y.RT2[i],nLines2), z.Dens2[,i],  col = baseColor, type="l")
                                  }
                                  a$points3d(current_state$mu[1],current_state$mu[2],max.Density, col = meanColor, cex=0.5, pch=16)
                          }
                      }
                  
                  } else {    # If RT is invalid, record rejection
                      change[i] <- FALSE  
                  }
            }
            
            # Evaluate the acceptance rate
            current_state$ARate_obs <- mean(change)
            
            # Track warmup progress
            if(plot_warmup | n_chains > 1) {
                warmup_history[[chain]] <- rbind(warmup_history[[chain]],
                                                 data.frame(iteration = iteration,
                                                            sigma_choice = sqrt(current_state$sigma[1,1]),
                                                            sigma_rt = sqrt(current_state$sigma[2,2]),
                                                            acceptance_rate = current_state$ARate_obs
                                                            )
                                                )
            }
      }# Close warmup loop
      
      # Store final state for this chain
      chain_states[[chain]] <- current_state
  }
  
  
  # Plot the theoretical CDDM density on the 3D scatter plot for reference  
  if(plot){ 
            baseColorCDDM <- rgb(0,0,0,0.5)  
            for(i in 1:nLines){  # Plot density contours
                a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l", col = baseColorCDDM)
                a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i],  col = baseColorCDDM, type="l")
            }
            a$points3d(x.theta, y.RT, theta.Dens, col = baseColorCDDM, type="l")
  }  

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Integrate WARMUP results across Chains
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # Compute average proposal distribution across chains
  final_sigma <- matrix(0, 2, 2)  # Start at 0
  final_mu <- c(0, 0)
  
  for(i in 1:n_chains) { # Add up the proposal parameters per chain
      final_sigma <- final_sigma + chain_states[[i]]$sigma
      final_mu <- final_mu + chain_states[[i]]$mu
  }
  final_sigma <- final_sigma / n_chains
  final_mu <- final_mu / n_chains

  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  # B E G I N   M A I N   S A M P L I N G   P H A S E
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
  samples <- matrix(NA, nrow=n, ncol=2)  # Initialize matrix to store samples
  u <- runif(n,0,1)                      # Generate n uniform criteria
  current <- final_mu                    # Start at the mean of the proposal distribution
  Sigma <- final_sigma                   # Use the integrated proposal distribution
    
  for(i in 1:n){      
        # Custom debugging output
        if(debug){  #rep_count <- 1
                    cat(sprintf("\rProgress: %d/%d (%.1f%%)",i, n, 100 * i/n))    
                 }        
        a <- 0
        repeat{ 
                #a <- a+1
                #cat(a)
                # Generate valid candidate
                cand <- generate_validMVN_candidate(current, Sigma, logRT, tzero, max.RT)          
                if(logRT){  current[2] <- exp(current[2])  }  # Log-transform RTs if needed

                # Compute acceptance ratio
                ratio.num <- max(dCDDM(cand, drift, theta, tzero, boundary), 0, na.rm=TRUE)
                ratio.den <- dCDDM(current, drift, theta, tzero, boundary)
                ratio <- ratio.num/ratio.den
                
                # Print debugging information
                if(debug){  cat(sprintf("  ratio: %.4f, criterion: %.4f\n", ratio, u[i]))  }            
                
                # Go back to proposal distribution logRT scale if needed
                if(logRT){  cand[2] <- log(cand[2])  }  # Log-transform RTs if needed

                # Take step if ratio exceeds rejection criterion          
                if(ratio > u[i]) {
                    current <- cand
                    break
                }
                #if(debug){ rep_count <- rep_count + 1  }
              }
      #cat("i=",i, "ratio=", ratio, "u=", u[i], "\n", "cand=", cand, "\n", "current=", current, "\n")
      # Store current state
      samples[i,] <- current      
  }
  
  if(logRT){  samples[,2] <- exp(samples[,2])  }
  # Rotate samples back to original theta
  samples[,1] <- (samples[,1] - rotation_angle) %% (2*pi)  
  
  if(plot){   # Plot final samples
              plot_samples <- samples
              plot_samples[,1] <- (plot_samples[,1] + rotation_angle) %% (2*pi)
              a$points3d(plot_samples[,1], plot_samples[,2],
                        dCDDM(as.matrix(plot_samples), drift, pi, tzero, boundary),
                        col = "white", 
                        pch = 16, cex = 0.5)
  }

  if(!plot && plot_warmup) {
              layout_mat <- matrix(c(1,2,    # top row: main plot
                                     3,4),   # bottom row: last two warmup plots
                                  nrow=2, ncol=2, byrow=TRUE)
              layout(layout_mat, heights=c(1,1), widths = c(1,1))  # Make the top row slightly taller
  }

  # Plot warmup diagnostics if requested
  if(plot_warmup && length(warmup_history) > 0) {
    # Now create the warmup diagnostic plots
    par(mar=c(4,4,2,1))  # Reset margins for the diagnostic plots
    
    # Find convergence iterations for each chain
    convergence_iters <- sapply(warmup_history, function(x) {
      which(x$acceptance_rate >= ARate_des)[1]
    })
    
    # Plot sigma evolution for choices - one line per chain
    plot(NULL, 
         xlim=range(sapply(warmup_history, function(x) x$iteration)),
         ylim=range(sapply(warmup_history, function(x) x$sigma_choice)),
         xlab="Iteration", ylab="Sigma (Choice)",
         main="Evolution of Choice SD")
    for(ch in 1:n_chains) {
      chain_data <- warmup_history[[ch]]
      lines(chain_data$iteration, chain_data$sigma_choice, 
            col=rainbow(n_chains)[ch], lwd=2)
      # Vertical line at convergence
      abline(v=convergence_iters[ch], col=rainbow(n_chains)[ch], lty=3)
      # Horizontal line at final value
      abline(h=chain_data$sigma_choice[nrow(chain_data)], 
            col=rainbow(n_chains)[ch], lty=2)
    }
    # Add thick black line for integrated final value
    abline(h=sqrt(final_sigma[1,1]), col="black", lwd=3)
    
    # Plot sigma evolution for RT - one line per chain
    plot(NULL,
         xlim=range(sapply(warmup_history, function(x) x$iteration)),
         ylim=range(sapply(warmup_history, function(x) x$sigma_rt)),
         xlab="Iteration", ylab="Sigma (RT)",
         main="Evolution of RT SD")
    for(ch in 1:n_chains) {
      chain_data <- warmup_history[[ch]]
      lines(chain_data$iteration, chain_data$sigma_rt, 
            col=rainbow(n_chains)[ch], lwd=2)
      # Vertical line at convergence
      abline(v=convergence_iters[ch], col=rainbow(n_chains)[ch], lty=3)
      # Horizontal line at final value
      abline(h=chain_data$sigma_rt[nrow(chain_data)], 
            col=rainbow(n_chains)[ch], lty=2)
    }
    # Add thick black line for integrated final value
    abline(h=sqrt(final_sigma[2,2]), col="black", lwd=3)
    
    # Plot acceptance rate - one line per chain
    plot(NULL,
         xlim=range(sapply(warmup_history, function(x) x$iteration)),
         ylim=c(0, max(1, max(unlist(sapply(warmup_history, function(x) x$acceptance_rate))))),
         xlab="Iteration", ylab="Acceptance Rate",
         main="Acceptance Rate Evolution")
    abline(h=ARate_des, lty=2, col="gray")
    for(ch in 1:n_chains) {
      chain_data <- warmup_history[[ch]]
      lines(chain_data$iteration, chain_data$acceptance_rate, 
            col=rainbow(n_chains)[ch], lwd=2)
      abline(v=convergence_iters[ch], col=rainbow(n_chains)[ch], lty=3)
    }
    
    # Plot trajectory in 2D - one path per chain
    plot(NULL,
         xlim=range(sapply(warmup_history, function(x) x$sigma_choice)),
         ylim=range(sapply(warmup_history, function(x) x$sigma_rt)),
         xlab="Sigma (Choice)", ylab="Sigma (RT)",
         main="Proposal SD Trajectory")
    for(ch in 1:n_chains) {
      chain_data <- warmup_history[[ch]]
      points(chain_data$sigma_choice, chain_data$sigma_rt,
             col=rainbow(n_chains, alpha=0.5)[ch], pch=16)
      if(nrow(chain_data) > 1) {
        arrows(chain_data$sigma_choice[-nrow(chain_data)],
               chain_data$sigma_rt[-nrow(chain_data)],
               chain_data$sigma_choice[-1],
               chain_data$sigma_rt[-1],
               length=0.1, col=rainbow(n_chains, alpha=0.3)[ch])
      }
      # Mark start and end points
      points(chain_data$sigma_choice[1], chain_data$sigma_rt[1],
             col="green", pch=16, cex=1.5)
      points(chain_data$sigma_choice[nrow(chain_data)],
             chain_data$sigma_rt[nrow(chain_data)],
             col="red", pch=16, cex=1.5)
    }
    # Add legend
    legend("topright", legend=paste("Chain", 1:n_chains),
           col=rainbow(n_chains), lwd=2, bty="n")
  }

 
  # Convert matrix to data frame and add column names
  samples <- as.data.frame(samples)
  colnames(samples) <- c("Choice", "RT")
  
  return(samples)  # Return data frame of samples
}


TryTEST <- FALSE
if(TryTEST){
start_time <- Sys.time()
pdf(paste0(here("results", "cddm_metropolis_warmup.pdf")), width=8, height=12)
x  <- rCDDM_Metropolis(n = 1000, par = list(mu1 = 0.5,  mu2 = 0.5, boundary = 5, tzero = 0.1), plot=TRUE, logRT = FALSE, plot_warmup = TRUE, n_chains = 5)
#y  <- rCDDM_Metropolis(n = 1000, par = list(mu1 = 0.5,  mu2 = 0.5, boundary = 5, tzero = 0.1), logRT = TRUE)
dev.off()
end_time <- Sys.time()
end_time - start_time

start_time <- Sys.time()
x  <- rCDDM_Metropolis(n = 1000, par = list(mu1 = 0.5,  mu2 = 0.5, boundary = 5, tzero = 0.1), plot=FALSE, logRT = FALSE, plot_warmup = FALSE, n_chains = 2)
end_time <- Sys.time()
end_time - start_time
}