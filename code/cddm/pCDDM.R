###############################################################################
###############################################################################
#####   Approximate the theoretical CDF function using the total volume
###############################################################################
########################################################   by Adriana F. Chavez 
##############################################################################
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