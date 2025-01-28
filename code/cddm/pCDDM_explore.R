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
    
    # Handle different types of CDF calculations
    if(type == "joint") {     rad <- data[1]; time <- data[2]        
    } else if(type == "RT") { rad <- 2*pi; time <- data[1]        
    } else if(type == "rad") { rad <- data[1]; time <- max_time } else {
        stop("Unknown type. Use 'joint', 'RT', or 'rad'")
    }
    
    # Input validation
    if(time < tzero) return(0)  # No probability mass before tzero
    if(rad <= 0) return(0)    # No probability mass at or below 0 rad
    
    if(method == "grid") {
        if(is.na(n_points)) {
            n_points <- 1000
        }
        # Grid-based numerical integration
        rad_grid <- seq(0, rad, length.out=n_points)
        time_grid <- seq(tzero, time, length.out=n_points)
        
        # Create meshgrid
        da <- rad_grid[2] - rad_grid[1]
        dt <- time_grid[2] - time_grid[1]
        
        # Create all combinations as a matrix
        grid_points <- as.matrix(expand.grid(rad=rad_grid, time=time_grid))
        
        # Vectorized computation
        densities <- sapply(1:nrow(grid_points), function(i) {
            dCDDM(grid_points[i,], drift, theta, tzero, boundary)
        })
        total_prob <- sum(densities) * da * dt
        
        return(total_prob)
        
    } else if(method == "monte_carlo") {        
        if(is.na(n_points)) {
            n_points <- 2000
        }
        # Monte Carlo integration
        points <- matrix(
            c(runif(n_points, 0, rad),     # rads
              runif(n_points, tzero, time)),  # times
            ncol=2
        )
        
        # Evaluate densities
        densities <- sapply(1:nrow(points), function(i) {
            dCDDM(points[i,], drift, theta, tzero, boundary)
        })
        
        # Monte Carlo integration formula
        area <- rad * (time - tzero)
        total_prob <- area * mean(densities)      
        
        return(total_prob)
    } else {
        stop("Unknown method. Use 'grid' or 'monte_carlo'")
    }
    
    # Ensure output is between 0 and 1
    total_prob <- pmax(0, pmin(1, total_prob))
    
    return(total_prob)
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
    
    # Test marginal CDFs
    p_rt <- pCDDM_RT(2, drift, theta, tzero, boundary)
    p_rad <- pCDDM_rad(pi, drift, theta, tzero, boundary)
    
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
    pdf("tests/pCDDM_testing.pdf", width=10, height=5)
    
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

