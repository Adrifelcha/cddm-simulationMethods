# Implements inverse CDF sampling for the Circular Drift Diffusion Model (CDDM)
# This method creates a lookup table of CDF values and uses interpolation to generate samples
rCDDM_inverse <- function(n, par,
                          n_angle = 100, n_rt = 100, max_rt = NULL,  # Make max_rt optional
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

    # Use keyDensityPoints to determine max_rt if not provided
    if (is.null(max_rt)) {
        key_points <- keyDensityPoints(par)
        max_rt <- key_points$max.RT
        # Optional: Use key points for RT grid in adaptive method
        pred_rt <- key_points$pred.RT
        min_rt <- key_points$min.RT
    } else {
        min_rt <- tzero
        pred_rt <- tzero + (boundary^2)/(2*drift)  # theoretical mode
    }

    # 1. Create sampling grid for angles and response times (RTs)
    if(grid_method == "regular") {
        # Use uniform spacing for both angles and RTs
        angles <- seq(0, 2*pi, length.out = n_angle)
        rts <- seq(min_rt, max_rt, length.out = n_rt)
    } else if(grid_method == "adaptive") {
        # Create non-uniform grid with more points in high-density regions
        # For angles: concentrate points around the drift direction (theta)
        angle_seq <- seq(0, 2*pi, length.out = 1000)
        angle_density <- dnorm(angle_seq, 
                             mean = theta %% (2*pi), 
                             sd = pi/4)
        angles <- sort(sample(angle_seq, 
                            size = n_angle, 
                            prob = angle_density))
        
        # For RTs: concentrate points around predicted RT
        rt_seq <- seq(min_rt, max_rt, length.out = 1000)
        rt_density <- dnorm(rt_seq, 
                          mean = pred_rt, 
                          sd = (max_rt - min_rt)/6)
        rts <- sort(sample(rt_seq, 
                          size = n_rt, 
                          prob = rt_density))
    }
    
    # Create grid of all angle-RT combinations
    grid_points <- expand.grid(angle = angles, rt = rts)
    
    # 2. Build CDF lookup table
    # Calculate CDF values for all grid points using Monte Carlo integration
    cdf_values <- pCDDM(as.matrix(grid_points), drift, theta, tzero, boundary,
                        method = "monte_carlo", n_points = 50000)
    
    # Clean up lookup table by removing invalid points and ensuring monotonicity
    valid_idx <- !is.na(cdf_values)
    grid_points <- grid_points[valid_idx, ]
    cdf_values <- cdf_values[valid_idx]
    lookup_table <- cbind(grid_points, cdf = cdf_values)
    lookup_table <- lookup_table[order(lookup_table$cdf), ]  # Sort by CDF values
    lookup_table <- lookup_table[!duplicated(lookup_table$cdf), ]  # Ensure monotonicity
    
    if(nrow(lookup_table) < 2) {
        stop("Not enough valid points in lookup table after removing NAs and duplicates")
    }
    
    # 3. Generate samples using inverse CDF method
    u <- runif(n)  # Generate uniform random numbers
    samples <- matrix(0, nrow = n, ncol = 2)
    
    # 4. Interpolate samples from lookup table
    for(i in 1:n) {
        # Use binary search to find closest CDF values efficiently
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
        idx <- min(max(left, 2), nrow(lookup_table))
        
        # Get bracketing points for interpolation
        lower_pt <- lookup_table[idx-1, ]
        upper_pt <- lookup_table[idx, ]
        
        # Interpolate between points based on uniform random value
        cdf_diff <- upper_pt$cdf - lower_pt$cdf
        if(cdf_diff > 0) {
            # Calculate interpolation weights
            w_upper <- (u[i] - lower_pt$cdf) / cdf_diff
            w_lower <- 1 - w_upper
            
            # Handle circular wraparound for angles (e.g., 0 and 2π are the same)
            if(abs(lower_pt$angle - upper_pt$angle) > pi) {
                if(lower_pt$angle > upper_pt$angle) {
                    upper_pt$angle <- upper_pt$angle + 2*pi
                } else {
                    lower_pt$angle <- lower_pt$angle + 2*pi
                }
            }
            
            # Linear interpolation for both dimensions
            samples[i, 1] <- (w_lower * lower_pt$angle + w_upper * upper_pt$angle) %% (2*pi)
            samples[i, 2] <- w_lower * lower_pt$rt + w_upper * upper_pt$rt
        } else {
            # If points have same CDF value, use the closer point
            if(abs(u[i] - lower_pt$cdf) < abs(u[i] - upper_pt$cdf)) {
                samples[i, ] <- c(lower_pt$angle, lower_pt$rt)
            } else {
                samples[i, ] <- c(upper_pt$angle, upper_pt$rt)
            }
        }
    }
    
    # Return samples as a data frame
    return(data.frame(Choice = samples[, 1], RT = samples[, 2]))
}

