###############################################################################
###############################################################################
#####   Approximate the theoretical CDF function using the total volume
###############################################################################
########################################################   by Adriana F. Chavez 
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){     source("code/cddm/dCDDM.R")       }
library("scatterplot3d")
library("plot3D")

# Auxiliary functions needed are available in pCDDM_aux.R
# integrate_grid()
# pCDDM_grid()

###########################################################################
# Parameters and variable definitions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data ~~~~  
# data <- cbind(rad, time)
# rad:  Upper bound for choice (in radians)
# time: Upper bound for response time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Parameters ~~~~
# drift: Drift rate magnitude
# theta: Drift angle (in radians)
# tzero: Non-decision time
# boundary: Decision boundary
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Method ~~~~
# method: Method for numerical integration ("grid" or "monte_carlo")
# n_points Number of points for numerical integration
# type: Type of CDF to compute ("joint", "RT", or "rad")
###########################################################################
pCDDM <- function(data, drift, theta, tzero, boundary, method="monte_carlo", 
                  n_points=NA, show=FALSE, type="joint", max_time=10) {    
    # Convert single vector to matrix if needed
    if(is.vector(data)){        data <- matrix(data, nrow=1)        }
    
    # Handle different types of CDF calculations
    if(type == "joint") {     
        rad <- data[,1]; time <- data[,2]        
    } else if(type == "RT") { 
        rad <- rep(2*pi, nrow(data)); time <- data[,1]        
    } else if(type == "rad") { 
        rad <- data[,1]; time <- rep(max_time, nrow(data))
    } else {
        stop("Unknown type. Use 'joint', 'RT', or 'rad'")
    }
    # Initialize output vector    
    probs <- rep(0, length(rad))  
    # Input validation
    valid_idx <- which(time >= tzero & rad > 0)    
    if(length(valid_idx) == 0){ return(probs) }
    
    if(method == "monte_carlo") {
        if(is.na(n_points)) { 
            if(nrow(data) > 1) { n_points <- 50000 }
            else { n_points <- 5000 }
        }
        
        # Find the maximum boundaries from all valid observations
        max_rad <- max(rad[valid_idx])
        max_time <- max(time[valid_idx])
        
        # Generate one large set of random points for the maximum area
        points <- matrix(0, nrow = n_points, ncol = 2)
        points[,1] <- runif(n_points, 0, max_rad)
        points[,2] <- runif(n_points, tzero, max_time)
        
        # Compute densities once for all points
        densities <- dCDDM(points, drift, theta, tzero, boundary)
        
        # For each observation, compute probability using points within its region
        probs[valid_idx] <- sapply(valid_idx, function(i) {
            # Find points within the current observation's region
            in_region <- points[,1] <= rad[i] & points[,2] <= time[i] & points[,2] >= tzero
            
            # Compute area ratio to adjust for different regions
            total_area <- max_rad * (max_time - tzero)
            region_area <- rad[i] * (time[i] - tzero)
            
            # Compute mean density adjusted by area ratio
            if(sum(in_region) > 0) {
                mean(densities[in_region]) * region_area
            } else {
                0
            }
        })
        
        # Visualization if requested
        if(show && length(valid_idx) > 0) {
            i <- valid_idx[1]  # Visualize first valid point
            plot_monte_carlo_cdf(points, densities, rad[i], time[i], tzero, probs[i])
        }
    } else if(method == "grid") {
        if(is.na(n_points)) {       n_points <- 1000        }
        # Similar optimization needed for grid method
        probs[valid_idx] <- sapply(valid_idx, function(i) {
            integrate_grid(rad[i], time[i], tzero, drift, theta, boundary, n_points)
        })
    }
    
    # Ensure output is between 0 and 1
    probs <- pmax(0, pmin(1, probs))
    
    return(probs)
}




# Test Monte Carlo convergence and timing
if(!exists("test")) { test <- TRUE }
if(test) {
    set.seed(789)
    # Test parameters
    drift <- 1
    theta <- pi/4
    tzero <- 0.1
    boundary <- 4
    data <- c(pi, 2)
    
    # Different n_points to test
    n_points_seq <- seq(500, 10000, 500)  # More evenly distributed sequence
    n_repetitions <- 100  # Number of times to repeat each test
    
    # Store results
    results <- data.frame()
    
    # Run tests
    for(n in n_points_seq) {
        cat("\nTesting n_points =", n, "\n")
        
        # Store results for this n_points
        values <- numeric(n_repetitions)
        times <- numeric(n_repetitions)
        
        for(i in 1:n_repetitions) {
            start_time <- Sys.time()
            values[i] <- pCDDM(data, drift, theta, tzero, boundary, 
                              method="monte_carlo", n_points=n)
            end_time <- Sys.time()
            times[i] <- as.numeric(end_time - start_time, units="secs")
        }
        
        # Compute statistics
        results <- rbind(results, data.frame(
            n_points = n,
            mean_value = mean(values),
            sd_value = sd(values),
            mean_time = mean(times),
            sd_time = sd(times)
        ))
    }
 
    # Print results
    print("Monte Carlo Results:")
    print(results)
    
    # Create PDF file
    pdf("results/pCDDM_testing.pdf", width=10, height=5)
    
    # Set up side-by-side plots
    par(mfrow=c(1,2))
    
    # Plot 1: Convergence of values
    plot(results$n_points, results$mean_value, type="l",
         xlab="Number of Points", ylab="CDF Value",
         main="Monte Carlo Convergence", 
         ylim=range(c(results$mean_value + c(-results$sd_value, results$sd_value))))
    # Add error bands
    polygon(c(results$n_points, rev(results$n_points)),
           c(results$mean_value + results$sd_value,
             rev(results$mean_value - results$sd_value)),
           col=rgb(0,0,0,0.2), border=NA)
    lines(results$n_points, results$mean_value, lwd=2)
    
    # Plot 2: Execution time
    plot(results$n_points, results$mean_time, type="l",
         xlab="Number of Points", ylab="Time (seconds)",
         main="Execution Time")
    # Add error bands
    polygon(c(results$n_points, rev(results$n_points)),
           c(results$mean_time + results$sd_time,
             rev(results$mean_time - results$sd_time)),
           col=rgb(0,0,0,0.2), border=NA)
    lines(results$n_points, results$mean_time, lwd=2)
    
    # Reset plot parameters and close PDF
    par(mfrow=c(1,1))
    dev.off()
}

