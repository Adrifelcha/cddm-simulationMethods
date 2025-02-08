###################################################################################
####   Different approximations for the CDF of the CDDM
###################################################################################

# Performs numerical integration using an adaptive grid method
integrate_grid <- function(rad, time, tzero, drift, theta, boundary, n_points=1000) {
    # Create parameter list for keyDensityPoints
    par <- list(
        drift = drift,
        theta = theta,
        tzero = tzero,
        boundary = boundary
    )
    
    # Get key density points
    key_points <- keyDensityPoints(par)
    
    # Create adaptive grid for angles (denser around theta)
    theta_norm <- theta %% (2*pi)
    if(theta_norm > rad) theta_norm <- theta_norm - 2*pi
    
    # Create base sequence and density weights for angles
    angle_seq <- seq(0, rad, length.out=1000)
    angle_density <- dnorm(angle_seq, mean=theta_norm, sd=pi/2) + 0.2  # Add baseline density
    angle_density <- angle_density / max(angle_density)
    
    # Sample angles with density-based weights
    rad_grid <- sort(sample(angle_seq, size=n_points, prob=angle_density, replace=TRUE))
    
    # Create adaptive grid for times using key points
    time_seq <- seq(key_points$min.RT, min(time, key_points$max.RT), length.out=1000)
    time_density <- dnorm(time_seq, mean=key_points$pred.RT, sd=(time-tzero)/4) + 0.2
    time_density <- time_density / max(time_density)
    
    time_grid <- sort(sample(time_seq, 
                            size=n_points, 
                            prob=time_density,
                            replace=TRUE))
    
    # Calculate total area for normalization
    total_area <- (rad) * (time - tzero)
    
    # Create grid points and calculate densities
    grid_points <- as.matrix(expand.grid(rad=rad_grid, time=time_grid))
    densities <- sapply(1:nrow(grid_points), function(j) {
        dCDDM(grid_points[j,], drift, theta, tzero, boundary)
    })
    
    # Normalize by total area and number of points
    mean(densities) * total_area
}

# Wrapper function for grid-based CDF computation
pCDDM_grid <- function(data, drift, theta, tzero, boundary, probs, valid_idx, n_points=NA, show=FALSE) {
    if(is.na(n_points)) { n_points <- 100 }
    
    # Only compute for valid indices
    probs[valid_idx] <- sapply(valid_idx, function(i) {
        integrate_grid(data[i,1], data[i,2], tzero, drift, theta, boundary, n_points)
    })
    
    pmax(0, pmin(1, probs))
}

# Monte Carlo method for CDF computation
pCDDM_monte_carlo <- function(data, drift, theta, tzero, boundary, probs, valid_idx, n_points=NA, show=FALSE) {
    if(is.na(n_points)) { 
        n_points <- if(nrow(data) > 1) 50000 else 7000
    }
    
    # Find the maximum boundaries from valid observations only
    max_rad <- max(data[valid_idx,1])
    max_time <- max(data[valid_idx,2])
    
    # Generate points and compute densities
    points <- matrix(0, nrow = n_points, ncol = 2)
    points[,1] <- runif(n_points, 0, max_rad)
    points[,2] <- runif(n_points, tzero, max_time)
    densities <- dCDDM(points, drift, theta, tzero, boundary)
    
    # Compute probabilities only for valid indices
    probs[valid_idx] <- sapply(valid_idx, function(i) {
        in_region <- points[,1] <= data[i,1] & 
                    points[,2] <= data[i,2] & 
                    points[,2] >= tzero
        
        total_area <- max_rad * (max_time - tzero)
        region_area <- data[i,1] * (data[i,2] - tzero)
        
        if(sum(in_region) > 0) {
            mean(densities[in_region]) * region_area
        } else {
            0
        }
    })
    
    if(show && length(valid_idx) > 0) {
        i <- valid_idx[1]
        plot_monte_carlo_cdf(points, densities, data[i,1], data[i,2], tzero, probs[i])
    }
    
    pmax(0, pmin(1, probs))
}

################################################################################
################################################################################
# P L O T T I N G   F U N C T I O N S ##########################################
################################################################################
################################################################################
# Helper function for visualization
plot_adaptive_grid <- function(rad, time, tzero, drift, theta, boundary, n_points) {
    # Generate grid points as in integrate_grid
    theta_norm <- theta %% (2*pi)
    if(theta_norm > rad) theta_norm <- theta_norm - 2*pi
    
    angle_seq <- seq(0, rad, length.out=1000)
    angle_density <- dnorm(angle_seq, mean=theta_norm, sd=pi/4) + 
                    dnorm(angle_seq, mean=theta_norm + 2*pi, sd=pi/4) +
                    dnorm(angle_seq, mean=theta_norm - 2*pi, sd=pi/4)
    angle_density <- angle_density / max(angle_density)
    rad_grid <- sort(sample(angle_seq, size=n_points, prob=angle_density, replace=TRUE))
    
    expected_peak <- tzero + (boundary^2)/(2*drift)
    time_seq <- seq(tzero, time, length.out=1000)
    time_density <- dnorm(time_seq, mean=expected_peak, sd=(time-tzero)/6) + 
                   dnorm(time_seq, mean=tzero, sd=(time-tzero)/8)
    time_density <- time_density / max(time_density)
    time_grid <- sort(sample(time_seq, size=n_points, prob=time_density, replace=TRUE))
    
    # Create plot
    plot(expand.grid(rad_grid, time_grid), 
         pch=16, cex=0.5, col=rgb(0,0,1,alpha=0.3),
         xlim=c(0, rad), ylim=c(tzero, time),
         xlab="Choice Angle (radians)", ylab="Response Time")
    
    # Add reference lines
    abline(v=theta_norm, col="red", lty=2)  # drift direction
    abline(h=expected_peak, col="green", lty=2)  # expected peak time
    
    # Add title and legend
    title("Adaptive Grid Integration")
    legend("topright", 
           legend=c("Grid Points", "Drift Direction", "Expected Peak Time"),
           pch=c(16, NA, NA),
           lty=c(NA, 2, 2),
           col=c(rgb(0,0,1,alpha=0.3), "red", "green"))
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
    in_region <- points[,1] <= rad & points[,2] <= time & points[,2] >= tzero
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
