###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                SEQUENTIAL REJECTION  ALGORITHM
###############################################################################
#####      Implementation using Ex-Gaussian and von Mises distributions
###############################################################################
########################################################   by Adriana F. Chávez 
library(mvtnorm)  # For multivariate normal distribution
library(circular)  # For von Mises distribution (rvm)

rCDDM_Reject_seq <- function(n, par, plot=FALSE, createPDF=FALSE){  
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Input validation
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    
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
    # Set up parameters for proposal distributions
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    kappa <- drift * boundary  # Concentration parameter for von Mises
    mu <- mean_rt
    sigma <- boundary/drift
    tau <- boundary/(2*drift)            

    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Optional visualization setup
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    if(plot|createPDF){
        if(createPDF){
            pdf(file = here("results", "show_Reject_Sequential.pdf"), width = 10, height = 10)
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
        # Step 1: Generate and filter RT candidates using marginal density
        rt_candidates <- rexGAUS(n, mu, sigma, tau)
        valid_rt <- rt_candidates >= min_RT
        if(sum(valid_rt) == 0) next
        
        rt_valid <- rt_candidates[valid_rt]
        
        # Calculate marginal RT density for valid candidates
        rt_target_density <- dCDDM_RT(rt_valid, drift, boundary, tzero)
        rt_max_density <- max(rt_target_density) * 1.1  # Add 10% buffer
        
        # Rejection sampling for RTs
        rt_rej_crit <- runif(length(rt_valid), 0, rt_max_density)
        rt_keep <- rt_target_density >= rt_rej_crit
        
        if(sum(rt_keep) == 0) next
        accepted_rts <- rt_valid[rt_keep]
        
        # Step 2: Generate matching choice samples for accepted RTs
        n_accepted_rt <- length(accepted_rts)
        choice_candidates <- rvonmises(n_accepted_rt, theta, kappa)
        
        # Calculate conditional choice density given accepted RTs
        valid_cand <- cbind(choice_candidates, accepted_rts)
        target_density <- dCDDM(valid_cand, drift, theta, tzero, boundary)
        
        # Final rejection step for choice values
        rej_crit <- runif(n_accepted_rt, 0, max_Density)
        keep <- target_density >= rej_crit
        
        # Visualize samples if plotting enabled
        if(plot|createPDF){
            a$points3d(valid_cand[!keep,1], valid_cand[!keep,2], rej_crit[!keep],
                      col = rgb(204/255,0,0,0.2), pch = 16, cex = 0.3)  # Rejected (red)
            a$points3d(valid_cand[keep,1], valid_cand[keep,2], rej_crit[keep],
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

#x <- rCDDM_Reject_seq(1000, par = list(mu1 = 2, mu2 = 1, drift = 1, theta = 0, tzero = 0, boundary = 1), plot = TRUE, createPDF = TRUE)

