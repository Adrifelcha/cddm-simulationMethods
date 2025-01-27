###############################################################################
###############################################################################
#####   Approximate the theoretical CDF function using the total volume
###############################################################################
########################################################   by Adriana F. Chavez 
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){     source("./dCDDM.R")       }
library("scatterplot3d")
library("plot3D")


# Parameters and variable definitions:
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Data ~~~~
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
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
pCDDM <- function(rad, time, drift, theta, tzero, boundary, 
                  method="grid", n_points=1000, show=FALSE) {
    
    # Input validation
    if(time < tzero) return(0)  # No probability mass before tzero
    if(rad <= 0) return(0)    # No probability mass at or below 0 rad
    
    if(method == "grid") {
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

# Calculate marginal CDF for response times only
pCDDM_RT <- function(time, drift, theta, tzero, boundary, 
                     method="monte_carlo", n_points=1000) {
    return(pCDDM(2*pi, time, drift, theta, tzero, boundary, 
                 method, n_points))
}

# Calculate marginal CDF for rads only
pCDDM_rad <- function(rad, drift, theta, tzero, boundary, 
                        method="monte_carlo", n_points=1000, max_time=10) {
    # Use a large but finite time instead of Inf
    return(pCDDM(rad, max_time, drift, theta, tzero, boundary, 
                 method, n_points))
}

# Example usage:
if(!exists("test")) { test <- TRUE }
if(test) {
    # Test parameters
    drift <- 1;  theta <- pi/4; tzero <- 0.1; boundary <- 4

    # Test full CDF
    start_time <- Sys.time()
    p1 <- pCDDM(pi, 2, drift, theta, tzero, boundary)
    end_time <- Sys.time()
    print(end_time - start_time)
    start_time2 <- Sys.time()
    p2 <- pCDDM(pi, 2, drift, theta, tzero, boundary, method="monte_carlo")
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

# Example usage with visualization:
if(test) {
    # Test parameters
    drift <- 1
    theta <- pi/4
    tzero <- 0.1
    boundary <- 5
    
    # Test with visualization
    p2 <- pCDDM(pi, 2, drift, theta, tzero, boundary, 
                method="monte_carlo", visualize=TRUE)
}