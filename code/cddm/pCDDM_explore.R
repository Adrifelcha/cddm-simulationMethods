#' Calculate the cumulative distribution function for CDDM
#' 
#' @param angle Upper bound for angle (in radians)
#' @param time Upper bound for response time
#' @param drift Drift rate magnitude
#' @param theta Drift angle (in radians)
#' @param tzero Non-decision time
#' @param boundary Decision boundary
#' @param method Method for numerical integration ("grid" or "monte_carlo")
#' @param n_points Number of points for numerical integration
#' @return Probability P(C ≤ angle, T ≤ time) between 0 and 1
pCDDM <- function(angle, time, drift, theta, tzero, boundary, 
                  method="grid", n_points=1000) {
    
    # Input validation
    if(time < tzero) return(0)  # No probability mass before tzero
    if(angle <= 0) return(0)    # No probability mass at or below 0 angle
    
    if(method == "grid") {
        # Grid-based numerical integration
        angle_grid <- seq(0, angle, length.out=n_points)
        time_grid <- seq(tzero, time, length.out=n_points)
        
        # Create meshgrid
        da <- angle_grid[2] - angle_grid[1]
        dt <- time_grid[2] - time_grid[1]
        
        # Create all combinations
        grid_points <- expand.grid(angle=angle_grid, time=time_grid)
        
        # Vectorized computation
        densities <- dCDDM(grid_points, drift, theta, tzero, boundary)
        total_prob <- sum(densities) * da * dt
        
    } else if(method == "monte_carlo") {
        # Monte Carlo integration
        points <- cbind(
            runif(n_points, 0, angle),     # angles
            runif(n_points, tzero, time)   # times
        )
        
        # Vectorized density computation
        densities <- dCDDM(points, drift, theta, tzero, boundary)
        
        # Monte Carlo integration formula
        area <- angle * (time - tzero)
        total_prob <- area * mean(densities)
        
    } else {
        stop("Unknown method. Use 'grid' or 'monte_carlo'")
    }
    
    # Ensure output is between 0 and 1
    total_prob <- pmax(0, pmin(1, total_prob))
    
    return(total_prob)
}

#' Calculate marginal CDF for response times only
#' 
#' @param time Upper bound for response time
#' @param drift Drift rate magnitude
#' @param theta Drift angle (in radians)
#' @param tzero Non-decision time
#' @param boundary Decision boundary
#' @param method Integration method
#' @param n_points Number of integration points
#' @return Probability P(T ≤ time)
pCDDM_RT <- function(time, drift, theta, tzero, boundary, 
                     method="grid", n_points=1000) {
    return(pCDDM(2*pi, time, drift, theta, tzero, boundary, 
                 method, n_points))
}

#' Calculate marginal CDF for angles only
#' 
#' @param angle Upper bound for angle (in radians)
#' @param drift Drift rate magnitude
#' @param theta Drift angle (in radians)
#' @param tzero Non-decision time
#' @param boundary Decision boundary
#' @param method Integration method
#' @param n_points Number of integration points
#' @return Probability P(C ≤ angle)
pCDDM_angle <- function(angle, drift, theta, tzero, boundary, 
                        method="grid", n_points=1000) {
    return(pCDDM(angle, Inf, drift, theta, tzero, boundary, 
                 method, n_points))
}

# Example usage:
if(!exists("test")) { test <- TRUE }
if(test) {
    # Test parameters
    drift <- 1
    theta <- pi/4
    tzero <- 0.1
    boundary <- 5
    
    # Test full CDF
    p1b <- pCDDM(pi, 2, drift, theta, tzero, boundary)
    p2 <- pCDDM(pi, 2, drift, theta, tzero, boundary, method="monte_carlo")
    
    # Test marginal CDFs
    p_rt <- pCDDM_RT(2, drift, theta, tzero, boundary)
    p_angle <- pCDDM_angle(pi, drift, theta, tzero, boundary)
    
    # Print results
    cat("Grid method:", p1, "\n")
    cat("Monte Carlo method:", p2, "\n")
    cat("RT marginal CDF:", p_rt, "\n")
    cat("Angle marginal CDF:", p_angle, "\n")
}