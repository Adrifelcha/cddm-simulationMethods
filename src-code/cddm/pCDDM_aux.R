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
    max_rad <- max(data[valid_idx,1])
    max_time <- max(data[valid_idx,2])

    # Normalize theta to ensure it's in the correct range for integration
    theta_norm <- theta %% (2*pi)
    if(theta_norm > max_rad) theta_norm <- theta_norm - 2*pi

    # Get theoretical moments and key points of the CDDM distribution
    par <- list(drift = drift, theta = theta,
                tzero = tzero, boundary = boundary)
    key_points <- keyDensityPoints(par)      # Get important points in the distribution
    mean_rt <- ezcddm_MRT(drift, boundary, tzero)  # Theoretical mean response time
    rt_var <- ezcddm_VRT(drift, boundary)    # Theoretical response time variance
    
    # Create adaptive angle grid with more robust density-based sampling
    angle_seq <- seq(0, 2*pi, length.out=n_points)
    angle_density <- dnorm(angle_seq, mean=theta_norm, sd=pi/8) + 
                    0.3 * dnorm(angle_seq, mean=0, sd=pi/10) +
                    0.1  # Increased minimum baseline to ensure positive probabilities
    
    # Apply softer threshold to angle density
    density_threshold <- max(angle_density) * 0.01  # Reduced threshold to 1%
    angle_density[angle_density < density_threshold] <- density_threshold  # Set to minimum instead of 0
    angle_density <- angle_density / sum(angle_density)
    rad_grid <- sort(sample(angle_seq, size=n_points/2, prob=angle_density, replace=FALSE))

    # Create adaptive time grid with more robust density-based sampling
    time_seq <- seq(tzero, key_points$max.RT, length.out=n_points)
    time_density <- 2 * dnorm(time_seq, mean=mean_rt, sd=sqrt(rt_var)) + 
                   3 * dnorm(time_seq, mean=key_points$pred.RT, sd=(key_points$max.RT - key_points$min.RT)/10) +
                   dnorm(time_seq, mean=tzero, sd=(key_points$max.RT-tzero)/12) +
                   0.1  # Increased minimum baseline
    
    # Apply softer threshold to time density
    time_threshold <- max(time_density) * 0.01  # Reduced threshold to 1%
    time_density[time_density < time_threshold] <- time_threshold  # Set to minimum instead of 0
    time_density <- time_density / sum(time_density)
    time_grid <- sort(sample(time_seq, size=n_points/2, prob=time_density, replace=FALSE))

    # Create grid points only up to data maximums
    rad_grid_filtered <- rad_grid[rad_grid <= max_rad]
    time_grid_filtered <- time_grid[time_grid <= max_time]
    grid_points <- as.matrix(expand.grid(rad=rad_grid_filtered, time=time_grid_filtered))
    
    # Compute densities only for the filtered grid points
    densities <- dCDDM(grid_points, drift, theta, tzero, boundary)
    
    # Compute importance sampling weights once
    weights <- outer(angle_density[findInterval(rad_grid_filtered, angle_seq)],
                    time_density[findInterval(time_grid_filtered, time_seq)])
    weights <- weights / sum(weights)
    
    # Optional visualization of the adaptive grid
    if(show && length(valid_idx) > 0) {
        plot_adaptive_grid(data, valid_idx, grid_points, densities, tzero, rotation,
                          drift, theta, boundary)
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
plot_adaptive_grid <- function(data, valid_idx, grid_points, densities, tzero, rotation=0, 
                             drift, theta, boundary, n_background=50) {
    # Load required packages
    if (!require("scatterplot3d")) {
        install.packages("scatterplot3d")
        library(scatterplot3d)
    }
    
    # If no grid points, return early
    if(nrow(grid_points) == 0) {
        warning("No grid points to plot")
        return(NULL)
    }
    
    # Create background grid for full distribution
    par <- list(drift=drift, theta=theta, tzero=tzero, boundary=boundary)
    key_points <- keyDensityPoints(par)
    
    # Create regular grid for background surface
    rad_bg <- seq(0, 2*pi, length.out=n_background)
    time_bg <- seq(tzero, key_points$max.RT, length.out=n_background)
    
    # Create density matrix for surface
    density_matrix <- matrix(0, n_background, n_background)
    for(i in 1:n_background) {
        for(j in 1:n_background) {
            density_matrix[i,j] <- dCDDM(c(rad_bg[i], time_bg[j]), 
                                       drift, theta, tzero, boundary)
        }
    }
    
    # Set up plotting parameters
    par(mfrow=c(1,1), mar=c(0, 1.5, 1, 0.5), oma=c(0, 0, 0, 0))
    
    # Create surface plot for background distribution
    persp_output <- persp(rad_bg, time_bg, density_matrix,
                         theta=45, phi=30, expand=0.7,
                         col=rgb(0.8, 0.8, 1, 0.3),  # Light blue transparent surface
                         border=rgb(0.7, 0.7, 1, 0.2),  # Slightly darker grid lines
                         xlab="Choice (radians)", 
                         ylab="Response Time",
                         zlab="Density",
                         main="",
                         ticktype="detailed",
                         zlim=c(0, max(density_matrix)),
                         nticks=5,
                         cex.axis=0.8,
                         box=TRUE)
    
    mtext("Adaptive Grid Integration\nCDDM Density Surface", side=3, f=2, 
          line=-2.5, cex=1.5, adj=0)

    # Convert 3D coordinates to 2D for plotting points
    convert_coords <- function(x, y, z) {
        trans3d(x, y, z, pmat=persp_output)
    }
    
    # Get unique grid points for each dimension
    unique_rad <- sort(unique(grid_points[,1]))
    unique_time <- sort(unique(grid_points[,2]))
    
    # Draw lines following the density surface for radians
    for(rad in unique_rad) {
        time_points <- grid_points[grid_points[,1] == rad, 2]
        density_points <- numeric(length(time_points))
        for(i in seq_along(time_points)) {
            density_points[i] <- dCDDM(c(rad, time_points[i]), 
                                     drift, theta, tzero, boundary)
        }
        line_2d <- convert_coords(rep(rad, length(time_points)),
                                time_points,
                                density_points)
        lines(line_2d, col=rgb(0, 0, 1, 0.3), lwd=1.5)
    }
    
    # Draw lines following the density surface for time
    for(time in unique_time) {
        rad_points <- grid_points[grid_points[,2] == time, 1]
        density_points <- numeric(length(rad_points))
        for(i in seq_along(rad_points)) {
            density_points[i] <- dCDDM(c(rad_points[i], time), 
                                     drift, theta, tzero, boundary)
        }
        line_2d <- convert_coords(rad_points,
                                rep(time, length(rad_points)),
                                density_points)
        lines(line_2d, col=rgb(0, 0, 1, 0.3), lwd=1.5)
    }
    
    # Add markers for observed data points if any exist
    if(length(valid_idx) > 0) {
        for(i in valid_idx) {
            # Calculate actual density at the observed point
            obs_density <- dCDDM(c(data[i,1], data[i,2]), drift, theta, tzero, boundary)
            
            # Add point at the actual density height
            data_point_2d <- convert_coords(data[i,1], data[i,2], obs_density)
            points(data_point_2d, col="red", pch=16, cex=0.7)
        }
    }
    
    # Add legend in top-right corner of plot window
    par(xpd=TRUE)
    legend("topright", 
           c("Full Distribution", "Grid Coverage", "Observed Data"), 
           col=c(rgb(0.7, 0.7, 1, 1),  # Color of the distribution grid lines
                 rgb(0, 0, 1, 1),      # Color of the adaptive grid
                 "red"), 
           pch=c(NA, NA, 16),
           lty=c(1, 1, NA),
           lwd=c(2 , 2, NA),
           pt.cex=1,
           cex=1.2,
           bty="n")
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
