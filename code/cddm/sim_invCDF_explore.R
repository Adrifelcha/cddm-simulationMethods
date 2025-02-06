rCDDM_inverse <- function(n, par,
                          n_angle = 100, n_rt = 100, max_rt = 10,
                          grid_method = "regular") {    
    # Extract parameters from input list
    mu1 <- par$mu1;       mu2 <- par$mu2        # Cartesian drift components
    drift <- par$drift;   theta <- par$theta     # Polar drift components
    tzero <- par$tzero;   boundary <- par$boundary # Non-decision time and boundary

    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
    # Parameter Validation: Check if we have either Cartesian or Polar coordinates
    noPolar <- is.null(theta) & is.null(drift)   # Check if polar coordinates missing
    noRect <- is.null(mu1) & is.null(mu2)        # Check if Cartesian coordinates missing
    if(noPolar){
        if(noRect){
        stop("Provide Cartesian or Polar coordinates", call. = FALSE)
        }else{
        # Convert Cartesian to polar coordinates if needed
        Mu <- rectToPolar(mu1, mu2)
        drift <- Mu$dLength;    par$drift <- drift  # Magnitude of drift
        theta <- Mu$dAngle;     par$theta <- theta  # Direction of drift
        }
    }

    # 1. Generate grid of angles and RTs based on specified method
    if(grid_method == "regular") {
        angles <- seq(0, 2*pi, length.out = n_angle)
        rts <- seq(tzero, max_rt, length.out = n_rt)
    } else if(grid_method == "adaptive") {
        # Create denser grid points where CDF changes more rapidly
        # For angles: denser around theta
        angle_seq <- seq(0, 2*pi, length.out = 1000)
        angle_density <- dnorm(angle_seq, 
                             mean = theta %% (2*pi), 
                             sd = pi/4)
        angles <- sort(sample(angle_seq, 
                            size = n_angle, 
                            prob = angle_density))
        
        # For RTs: denser near tzero and expected peak
        expected_peak <- tzero + (boundary^2)/(2*drift)
        rt_seq <- seq(tzero, max_rt, length.out = 1000)
        rt_density <- dnorm(rt_seq, 
                          mean = expected_peak, 
                          sd = (max_rt - tzero)/6)
        rts <- sort(sample(rt_seq, 
                          size = n_rt, 
                          prob = rt_density))
    }
    
    # Create all combinations
    grid_points <- expand.grid(angle = angles, rt = rts)
    
    # 2. Calculate CDF values for all grid points
    cdf_values <- pCDDM(as.matrix(grid_points), drift, theta, tzero, boundary,
                        method = "monte_carlo", n_points = 50000)
    
    # Remove any NA values from the lookup table
    valid_idx <- !is.na(cdf_values)
    grid_points <- grid_points[valid_idx, ]
    cdf_values <- cdf_values[valid_idx]
    
    # Store grid and CDF values in a matrix
    lookup_table <- cbind(grid_points, cdf = cdf_values)
    
    # Sort by CDF values for faster interpolation
    lookup_table <- lookup_table[order(lookup_table$cdf), ]
    
    # Remove duplicate CDF values to ensure monotonicity
    lookup_table <- lookup_table[!duplicated(lookup_table$cdf), ]
    
    # Ensure we have enough valid points
    if(nrow(lookup_table) < 2) {
        stop("Not enough valid points in lookup table after removing NAs and duplicates")
    }
    
    # 3. Generate uniform random numbers
    u <- runif(n)
    
    # 4. Find closest CDF values and interpolate
    samples <- matrix(0, nrow = n, ncol = 2)
    
    # Use binary search for each uniform value
    for(i in 1:n) {
        # Binary search for closest CDF value
        left <- 1
        right <- nrow(lookup_table)
        
        while(left <= right) {
            mid <- floor((left + right)/2)
            if(lookup_table$cdf[mid] < u[i]) {
                left <- mid + 1
            } else {
                right <- mid - 1
            }
        }
        idx <- left
        
        # Ensure idx is within bounds
        idx <- min(max(idx, 2), nrow(lookup_table))
        
        # Interpolate between closest points
        lower_pt <- lookup_table[idx-1, ]
        upper_pt <- lookup_table[idx, ]
        
        # Calculate weights based on distance to target u
        cdf_diff <- upper_pt$cdf - lower_pt$cdf
        if(cdf_diff > 0) {  # Ensure we don't divide by zero
            w_upper <- (u[i] - lower_pt$cdf) / cdf_diff
            w_lower <- 1 - w_upper
            
            # Handle angle wrap-around
            if(abs(lower_pt$angle - upper_pt$angle) > pi) {
                if(lower_pt$angle > upper_pt$angle) {
                    upper_pt$angle <- upper_pt$angle + 2*pi
                } else {
                    lower_pt$angle <- lower_pt$angle + 2*pi
                }
            }
            
            # Linear interpolation for both angle and RT
            samples[i, 1] <- (w_lower * lower_pt$angle + w_upper * upper_pt$angle) %% (2*pi)
            samples[i, 2] <- w_lower * lower_pt$rt + w_upper * upper_pt$rt
        } else {
            # If we can't interpolate, use the closer point
            if(abs(u[i] - lower_pt$cdf) < abs(u[i] - upper_pt$cdf)) {
                samples[i, 1] <- lower_pt$angle
                samples[i, 2] <- lower_pt$rt
            } else {
                samples[i, 1] <- upper_pt$angle
                samples[i, 2] <- upper_pt$rt
            }
        }
    }
    
    # Convert to data frame with proper column names
    samples <- data.frame(
        Choice = samples[, 1],
        RT = samples[, 2]
    )
    
    return(samples)
}

