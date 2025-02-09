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
    
    # Create data matrix
    data_matrix <- cbind(rad, time)
    
    # Call appropriate method
    if(method == "monte_carlo") {
        approximate_cdf <- pCDDM_monte_carlo(data_matrix, drift, theta, tzero, boundary, probs, valid_idx, n_points, show)
    } else if(method == "grid") {
        approximate_cdf <- pCDDM_grid(data_matrix, drift, theta, tzero, boundary, probs, valid_idx, n_points, show)
    }
    
    return(approximate_cdf)
}

