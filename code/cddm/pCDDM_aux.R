###################################################################################
####   Auxiliary functions for approximating the CDF of the CDDM
###################################################################################
###################################################################################
################################################### By Adriana F. Chavez De la Peña
###################################################################################


###################################################################################
####   Numerical integration using an adaptive grid method
###################################################################################
# Grid-based method for approximating the Cumulative Distribution Function (CDF)
# of the Circular Drift Diffusion Model (CDDM). This method uses adaptive grid
# sampling focused on regions of high density.
pCDDM_grid <- function(data, drift, theta, tzero, boundary, probs, valid_idx, n_points=NA, show=FALSE, rotation=0) {
    if(is.na(n_points)) { 
        n_points <- if(nrow(data) > 1) 500 else 1000 
    }

    # Get maximum bounds for integration from all valid data points
    max.rad <- max(data[valid_idx,1])
    max.time <- max(data[valid_idx,2])

    # Normalize theta to ensure it's in the correct range for integration
    theta_norm <- theta %% (2*pi)
    if(theta_norm > max.rad) theta_norm <- theta_norm - 2*pi

    # Get theoretical moments and key points of the CDDM distribution
    par <- list(drift = drift, theta = theta,
                tzero = tzero, boundary = boundary)
    key_points <- keyDensityPoints(par)      # Get important points in the distribution
    mean_rt <- ezcddm_MRT(drift, boundary, tzero)  # Theoretical mean response time
    rt_var <- ezcddm_VRT(drift, boundary)    # Theoretical response time variance
    
    # Create adaptive angle grid (once for all data points)
    angle_seq <- seq(0, max.rad, length.out=n_points)
    angle_density <- dnorm(angle_seq, mean=theta_norm, sd=pi/6) + 
                    0.5 * dnorm(angle_seq, mean=0, sd=pi/8) +  
                    0.2  
    angle_density <- angle_density / sum(angle_density)
    rad_grid <- sort(sample(angle_seq, size=n_points/2, prob=angle_density, replace=FALSE))

    # Create adaptive time grid (once for all data points)
    time_seq <- seq(tzero, max.time, length.out=n_points)
    time_density <- dnorm(time_seq, mean=mean_rt, sd=sqrt(rt_var)) + 
                   2 * dnorm(time_seq, mean=key_points$pred.RT, sd=(key_points$max.RT - key_points$min.RT)/8) +
                   dnorm(time_seq, mean=tzero, sd=(max.time-tzero)/10) +
                   0.2  
    time_density <- time_density / sum(time_density)
    time_grid <- sort(sample(time_seq, size=n_points/2, prob=time_density, replace=FALSE))

    # Create grid points and compute densities once
    grid_points <- as.matrix(expand.grid(rad=rad_grid, time=time_grid))
    densities <- dCDDM(grid_points, drift, theta, tzero, boundary)
    
    # Compute importance sampling weights once
    weights <- outer(angle_density[findInterval(rad_grid, angle_seq)],
                    time_density[findInterval(time_grid, time_seq)])
    weights <- weights / sum(weights)
    
    # Optional visualization of the adaptive grid
    if(show && length(valid_idx) > 0) {
        plot_adaptive_grid(data, valid_idx, grid_points, densities, tzero, rotation)
    }
    
    # For each valid observation, compute its CDF using the pre-computed grid
    probs[valid_idx] <- sapply(valid_idx, function(i) {
        # Find points that fall within the region for this observation
        in_region <- grid_points[,1] <= data[i,1] & 
                    grid_points[,2] <= data[i,2]
        
        # Compute probability for this observation using points in region
        if(sum(in_region) > 0) {
            region_prob <- sum(densities[in_region] * weights[in_region]) * 
                          (data[i,1] * (data[i,2] - tzero))
            pmax(0, pmin(1, region_prob))
        } else {
            0
        }
    })
    
    return(probs)
}

# Monte Carlo method for CDF computation
pCDDM_monte_carlo <- function(data, drift, theta, tzero, boundary, probs, valid_idx, n_points=NA, show=FALSE, rotation=0) {
    # Set default number of Monte Carlo points based on data size
    # More points for multiple observations, fewer for single observation
    if(is.na(n_points)) { 
        n_points <- if(nrow(data) > 1) 70000 else 7000
    }
    
    # Determine the integration bounds from valid data points
    # This defines the rectangle over which we'll integrate
    max_rad <- max(data[valid_idx,1])   # Maximum angle in radians
    max_time <- max(data[valid_idx,2])   # Maximum response time
    
    # Generate uniform random points within the integration bounds
    # These points will be used to estimate the CDF through Monte Carlo integration
    points <- matrix(0, nrow = n_points, ncol = 2)
    points[,1] <- runif(n_points, 0, max_rad)      # Random angles
    points[,2] <- runif(n_points, tzero, max_time) # Random times
    
    # Calculate CDDM density at each random point
    densities <- dCDDM(points, drift, theta, tzero, boundary)
    
    # For each valid observation, compute the cumulative probability
    probs[valid_idx] <- sapply(valid_idx, function(i) {
        # Identify points that fall within the region of integration
        # i.e., points with angle <= observed angle AND time <= observed time
        in_region <- points[,1] <= data[i,1] & 
                    points[,2] <= data[i,2]
        
        # Calculate area for the region up to observation point
        region_area <- data[i,1] * (data[i,2] - tzero) # Area up to observation point
        
        # Compute probability using Monte Carlo integration:
        # 1. If we have points in the region, average their densities and multiply by region area
        # 2. If no points in region, return 0
        if(sum(in_region) > 0) {
            mean(densities[in_region]) * region_area
        } else {
            0
        }
    })
    
    # Optional visualization of the Monte Carlo approximation
    if(show && length(valid_idx) > 0) {
        i <- valid_idx[1]
        pdf("monte_carlo_cdf.pdf")
        plot_monte_carlo_cdf(points, densities, data[i,1], data[i,2], tzero, probs[i])
        dev.off()
    }
    
    # Ensure probabilities are bounded between 0 and 1
    pmax(0, pmin(1, probs))
}

################################################################################
################################################################################
# P L O T T I N G   F U N C T I O N S ##########################################
################################################################################
################################################################################
# Plot the adaptive grid using pre-computed points and densities
plot_adaptive_grid <- function(data, valid_idx, grid_points, densities, tzero, rotation=0) {
    # Load required package
    if (!require("scatterplot3d")) {
        install.packages("scatterplot3d")
        library(scatterplot3d)
    }
    
    # Create custom x-axis labels
    x_at <- seq(0, 2*pi, by=pi/2)  # Tick marks at 0, π/2, π, 3π/2, 2π
    x_labels <- sprintf("%.2f", (x_at - rotation) %% (2*pi))  # Convert to original angles
    
    # Set up plotting parameters
    par(mfrow=c(1,1), mar=c(5, 3, 6, 3), oma=c(0, 0, 0.3, 0))
    
    # Create 3D scatterplot with custom x-axis
    s3d <- scatterplot3d(grid_points[,1], grid_points[,2], densities,
                        xlab="Choice (radians)", ylab="Response Time",
                        zlab="Density", color=rgb(0, 0, 1, 0.3),
                        main="Adaptive Grid Integration\nCDDM Density Surface",
                        pch=19, angle=45, cex.main=1.5, cex.lab=1.2,
                        x.ticklabs=x_labels, xlim=c(0, 2*pi))
    
    # Add points on the bottom plane
    s3d$points3d(grid_points[,1], grid_points[,2], 
                 rep(0, nrow(grid_points)),
                 col=rgb(0.7, 0.7, 0.7, 0.5), pch=19)
    
    # Add vertical lines connecting points to their densities
    for(i in 1:nrow(grid_points)) {
        s3d$points3d(rep(grid_points[i,1], 2), 
                    rep(grid_points[i,2], 2),
                    c(0, densities[i]), 
                    type="l",
                    col=rgb(0, 0, 1, 0.1))
    }
    
    # Add markers for observed data points
    for(i in valid_idx) {
        s3d$points3d(rep(data[i,1], 2), 
                    rep(data[i,2], 2),
                    c(0, max(densities)), 
                    type="l", 
                    col="red", lwd=2)
        # Add point marker at the base
        s3d$points3d(data[i,1], data[i,2], 0,
                    col="red", pch=16, cex=1.5)
    }
    
    # Add legend
    legend("topright", 
           c("Density Points", "Grid Base Points", "Observed Data"), 
           col=c(rgb(0, 0, 1, 0.3), 
                 rgb(0.7, 0.7, 0.7, 0.5),
                 "red"), 
           pch=c(19, 19, 16), 
           cex=1.2, pt.cex=1.2)
}

plot_monte_carlo_cdf <- function(points, densities, rad, time, tzero, prob) {
    # Load required package
    if (!require("scatterplot3d")) {
        install.packages("scatterplot3d")
        library(scatterplot3d)
    }
    
    # Set up plotting parameters
    par(mfrow=c(1,1), mar = c(5, 3, 6, 3), oma = c(0, 0, 0.3, 0))
    
    # Filter points for the region
    in_region <- points[,1] <= rad & points[,2] <= time
    points_i <- points[in_region, ]
    densities_i <- densities[in_region]
    
    # Create base 3D scatterplot
    s3d <- scatterplot3d(points_i[,1], points_i[,2], densities_i,
                       xlab = "Choice (radians)", ylab = "Response Time",
                       zlab = "Density", color = rgb(0, 0, 1, 0.3),
                       main = paste("Monte Carlo Approximation to CDF\n", 
                                  "Choice =", round(rad, 3), 
                                  ", RT =", round(time, 3),
                                  "\nEstimated CDF =", round(prob, 4)),                       
                       pch = 19, angle = 45, cex.main = 1.5, cex.lab = 1.2)    
    # Add points on the bottom plane
    s3d$points3d(points_i[,1], points_i[,2], rep(0, nrow(points_i)),
                col = rgb(0.7, 0.7, 0.7, 0.5), pch = 19)
    
    # Add vertical lines connecting base points to density heights
    for(j in 1:nrow(points_i)) {
        s3d$points3d(rep(points_i[j,1], 2), rep(points_i[j,2], 2),
                    c(0, densities_i[j]), type = "l",
                    col = rgb(0, 0, 1, 0.1))
    }
    
    # Add integration region boundaries
    s3d$points3d(c(0, rad, rad, 0, 0), 
                c(tzero, tzero, time, time, tzero), 
                rep(0, 5), type = "l", 
                col = "red", lwd = 2)
    
    # Add target point (choice, rt) with projections
    s3d$points3d(rad, time, 0, col = "darkred", pch = 19, cex = 2)
    
    # Projection lines to axes
    s3d$points3d(c(rad, rad), c(tzero, time), c(0, 0),
                type = "l", col = "darkred", lty = 2)
    s3d$points3d(c(0, rad), c(time, time), c(0, 0),
                type = "l", col = "darkred", lty = 2)
    
    # Add legend
    legend("topright", 
           c("Density Points", "Random Samples", 
             "Integration Region", "Data Point"), 
           col = c(rgb(0, 0, 1, 0.3), 
                  rgb(0.7, 0.7, 0.7, 0.5),
                  "red", "darkred"), 
           pch = c(19, 19, NA, 19), lty = c(NA, NA, 1, NA),
           cex = 1.4, pt.cex = 1.5)
}
