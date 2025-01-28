###############################################################################
###############################################################################
#####   Approximate the theoretical CDF function using the total volume
###############################################################################
########################################################   by Adriana F. Chavez 
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){     source("code/cddm/dCDDM.R")       }
library("scatterplot3d")
library("plot3D")


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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pCDDM <- function(data, drift, theta, tzero, boundary, 
                  method="monte_carlo", n_points=NA, show=FALSE,
                  type="joint", max_time=10) {
    
    # Convert single vector to matrix if needed
    if(is.vector(data)) {
        data <- matrix(data, nrow=1)
    }
    
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
    
    # Input validation
    probs <- rep(0, length(rad))  # Initialize output vector
    valid_idx <- which(time >= tzero & rad > 0)
    
    if(length(valid_idx) == 0) return(probs)
    
    if(method == "monte_carlo") {
        if(is.na(n_points)) {       n_points <- 4000        }
        n_valid <- length(valid_idx)
        
        # Generate all points at once
        all_points <- matrix(0, nrow = n_points * n_valid, ncol = 2)
        
        # Create indices for each observation
        obs_indices <- rep(1:n_valid, each = n_points)
        
        # Generate random points for all observations at once
        all_points[,1] <- runif(n_points * n_valid) * rep(rad[valid_idx], each = n_points)
        all_points[,2] <- runif(n_points * n_valid) * rep(time[valid_idx] - tzero, each = n_points) + tzero
        
        # Compute all densities at once
        all_densities <- dCDDM(all_points, drift, theta, tzero, boundary)
        
        # Reshape and compute means
        densities_matrix <- matrix(all_densities, nrow = n_points, ncol = n_valid)
        areas <- rad[valid_idx] * (time[valid_idx] - tzero)
        probs[valid_idx] <- areas * colMeans(densities_matrix)
        
        # Visualization if requested
        if(show && length(valid_idx) > 0) {
            par(mfrow=c(1,1), mar = c(5, 3, 6, 3), oma = c(0, 0, 0.3, 0))
            # Load required package
            if (!require("scatterplot3d")) {
                install.packages("scatterplot3d")
                library(scatterplot3d)
            }
                                    
            i <- valid_idx[1]  # Visualize first valid point
            points_i <- all_points[obs_indices == 1, ]
            densities_i <- all_densities[obs_indices == 1]
            
            # Create base 3D scatterplot
            s3d <- scatterplot3d(points_i[,1], points_i[,2], densities_i,
                               xlab = "Choice (radians)",
                               ylab = "Response Time",
                               zlab = "Density",
                               main = paste("Monte Carlo Approximation to CDF\n", 
                                          "Choice =", round(rad[i], 3), 
                                          ", RT =", round(time[i], 3),
                                          "\nEstimated CDF =", round(probs[i], 4)),
                               color = rgb(0, 0, 1, 0.3),
                               pch = 19,
                               angle = 45,
                               cex.main = 1.5,    
                               cex.lab = 1.2)
            
            # Add points on the bottom plane
            s3d$points3d(points_i[,1], points_i[,2], 
                        rep(0, nrow(points_i)),
                        col = rgb(0.7, 0.7, 0.7, 0.5),
                        pch = 19)
            
            # Add vertical lines connecting base points to density heights
            for(j in 1:nrow(points_i)) {
                s3d$points3d(rep(points_i[j,1], 2),
                            rep(points_i[j,2], 2),
                            c(0, densities_i[j]),
                            type = "l",
                            col = rgb(0, 0, 1, 0.1))
            }
            
            # Add integration region boundaries
            s3d$points3d(c(0, rad[i], rad[i], 0, 0), 
                        c(tzero, tzero, time[i], time[i], tzero), 
                        rep(0, 5), 
                        type = "l", 
                        col = "red",
                        lwd = 2)
            
            # Add target point (choice, rt) with projections
            s3d$points3d(rad[i], time[i], 0,
                        col = "darkred",
                        pch = 19,
                        cex = 2)
            
            # Projection lines to axes
            s3d$points3d(c(rad[i], rad[i]), c(tzero, time[i]), c(0, 0),
                        type = "l", col = "darkred", lty = 2)
            s3d$points3d(c(0, rad[i]), c(time[i], time[i]), c(0, 0),
                        type = "l", col = "darkred", lty = 2)
            
            # Add legend
            legend("topright", 
                   c("Density Points", "Random Samples", 
                     "Integration Region", "Data Point"), 
                   col = c(rgb(0, 0, 1, 0.3), 
                          rgb(0.7, 0.7, 0.7, 0.5),
                          "red",
                          "darkred"), 
                   pch = c(19, 19, NA, 19),
                   lty = c(NA, NA, 1, NA),
                   cex = 1.4,      # Increased legend text size
                   pt.cex = 1.5)   # Increased legend symbol size
        }
    } else if(method == "grid") {
        if(is.na(n_points)) {       n_points <- 1000        }
        # Similar optimization needed for grid method
        probs[valid_idx] <- sapply(valid_idx, function(i) {
            rad_grid <- seq(0, rad[i], length.out=n_points)
            time_grid <- seq(tzero, time[i], length.out=n_points)
            
            da <- rad_grid[2] - rad_grid[1]
            dt <- time_grid[2] - time_grid[1]
            
            grid_points <- as.matrix(expand.grid(rad=rad_grid, time=time_grid))
            
            densities <- sapply(1:nrow(grid_points), function(j) {
                dCDDM(grid_points[j,], drift, theta, tzero, boundary)
            })
            sum(densities) * da * dt
        })
    }
    
    # Ensure output is between 0 and 1
    probs <- pmax(0, pmin(1, probs))
    
    return(probs)
}

# Example usage:
if(!exists("test")) { test <- TRUE }
if(test) {
    # Test parameters
    drift <- 1;  theta <- pi/4; tzero <- 0.1; boundary <- 4

    # Test full CDF
    start_time <- Sys.time()
    p1 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="grid")
    end_time <- Sys.time()
    print(end_time - start_time)
    start_time2 <- Sys.time()
    p2 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="monte_carlo", n_points=2000)
    end_time2 <- Sys.time()
    print(end_time2 - start_time2)
    start_time3 <- Sys.time()
    p3 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="monte_carlo", n_points=2000, show=TRUE)
    end_time3 <- Sys.time()
    print(end_time3 - start_time3)
    
    # Test marginal CDFs
    p_rt <- pCDDM(c(2, 1), drift, theta, tzero, boundary, type="RT")
    p_rad <- pCDDM(c(pi, 0.1), drift, theta, tzero, boundary, type="rad")
    
    # Print results
    cat("Grid method:", p1, "\n")
    cat("Monte Carlo method:", p2, "\n")
    cat("RT marginal CDF:", p_rt, "\n")
    cat("Rad marginal CDF:", p_rad, "\n")
}

# Test Monte Carlo convergence and timing
if(!exists("test")) { test <- TRUE }
if(test) {
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
    pdf("tests/pCDDM_testing3.pdf", width=10, height=5)
    
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

