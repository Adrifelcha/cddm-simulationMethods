rCDDM_inverse <- function(n, drift, theta, tzero, boundary, 
                         n_angle = 50, n_rt = 50, max_rt = 2) {
    
    # 1. Generate grid of angles and RTs
    angles <- seq(0, 2*pi, length.out = n_angle)
    rts <- seq(tzero, max_rt, length.out = n_rt)
    
    # Create all combinations
    grid_points <- expand.grid(angle = angles, rt = rts)
    
    # 2. Calculate CDF values for all grid points at once
    cdf_values <- pCDDM(as.matrix(grid_points), drift, theta, tzero, boundary,
                        method = "monte_carlo", n_points = 1000)
    
    # Store grid and CDF values in a matrix for easier searching
    lookup_table <- cbind(grid_points, cdf = cdf_values)
    
    # 3. Generate uniform random numbers
    u <- runif(n)
    
    # 4. Find closest CDF values and interpolate
    samples <- matrix(0, nrow = n, ncol = 2)
    for(i in 1:n) {
        # Find closest CDF value
        idx <- which.min(abs(lookup_table$cdf - u[i]))
        
        # Simple interpolation if not at edges
        if(idx > 1 && idx < nrow(lookup_table)) {
            # Linear interpolation between closest points
            w1 <- abs(lookup_table$cdf[idx] - u[i])
            w2 <- abs(lookup_table$cdf[idx-1] - u[i])
            total_weight <- w1 + w2
            
            samples[i,1] <- (lookup_table$angle[idx] * w2 + 
                            lookup_table$angle[idx-1] * w1) / total_weight
            samples[i,2] <- (lookup_table$rt[idx] * w2 + 
                            lookup_table$rt[idx-1] * w1) / total_weight
        } else {
            # If at edges, use exact values
            samples[i,1] <- lookup_table$angle[idx]
            samples[i,2] <- lookup_table$rt[idx]
        }
    }
    
    return(samples)
}

# Test the function
test_simulation <- function(n = 1000) {
    start_time <- Sys.time()
    
    samples <- rCDDM_inverse(n, 
                            drift = 1, 
                            theta = pi/4, 
                            tzero = 0.1, 
                            boundary = 4)
    
    end_time <- Sys.time()
    
    # Print timing
    difftime(end_time, start_time, units = "secs")
    
    # Plot results
    par(mfrow = c(2,2))
    
    # Histogram of angles
    hist(samples[,1], main = "Angles", xlab = "Angle (radians)",
         breaks = 30, probability = TRUE)
    
    # Histogram of RTs
    hist(samples[,2], main = "Response Times", xlab = "RT (seconds)",
         breaks = 30, probability = TRUE)
    
    # Scatterplot
    plot(samples[,1], samples[,2], 
         main = "Joint Distribution",
         xlab = "Angle (radians)", 
         ylab = "RT (seconds)",
         pch = 19, col = rgb(0,0,1,0.3))
    
    # Return samples invisibly
    invisible(samples)
}