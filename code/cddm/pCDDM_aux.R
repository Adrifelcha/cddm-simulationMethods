
###################################################################################
####   Different approximations for the CDF of the CDDM
###################################################################################

integrate_grid <- function(rad, time, tzero, drift, theta, boundary, n_points=1000) {
    rad_grid <- seq(0, rad, length.out=n_points)
    time_grid <- seq(tzero, time, length.out=n_points)
    
    da <- rad_grid[2] - rad_grid[1]
    dt <- time_grid[2] - time_grid[1]
    
    grid_points <- as.matrix(expand.grid(rad=rad_grid, time=time_grid))
    
    densities <- sapply(1:nrow(grid_points), function(j) {
        dCDDM(grid_points[j,], drift, theta, tzero, boundary)
    })
    sum(densities) * da * dt
}

pCDDM_grid <- function(data, drift, theta, tzero, boundary, n_points=1000) {
    # Initialize output vector    
    probs <- rep(0, length(data[,1]))
    
    # Input validation
    valid_idx <- which(data[,2] >= tzero & data[,1] > 0)    
    if(length(valid_idx) == 0){ return(probs) }
    
    # Compute probabilities for valid indices
    probs[valid_idx] <- sapply(valid_idx, function(i) {
        integrate_grid(data[i,1], data[i,2], tzero, n_points, drift, theta, boundary)
    })
    
    # Ensure output is between 0 and 1
    pmax(0, pmin(1, probs))
}