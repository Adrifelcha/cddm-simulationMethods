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
                       xlab = "Choice (radians)",
                       ylab = "Response Time",
                       zlab = "Density",
                       main = paste("Monte Carlo Approximation to CDF\n", 
                                  "Choice =", round(rad, 3), 
                                  ", RT =", round(time, 3),
                                  "\nEstimated CDF =", round(prob, 4)),
                       color = rgb(0, 0, 1, 0.3),
                       pch = 19,
                       angle = 45,
                       cex.main = 1.5,    
                       cex.lab = 1.2)
    
    # Add points on the bottom plane
    s3d$points3d(points_i[,1], points_i[,2], 
                rep(0, nrow(points_i)),
                col = rgb(0.7, 0.7, 0.7, 0.5),
                pch = 19)
    
    # Add vertical lines connecting base points to density heights
    for(j in 1:nrow(points_i)) {
        s3d$points3d(rep(points_i[j,1], 2),
                    rep(points_i[j,2], 2),
                    c(0, densities_i[j]),
                    type = "l",
                    col = rgb(0, 0, 1, 0.1))
    }
    
    # Add integration region boundaries
    s3d$points3d(c(0, rad, rad, 0, 0), 
                c(tzero, tzero, time, time, tzero), 
                rep(0, 5), 
                type = "l", 
                col = "red",
                lwd = 2)
    
    # Add target point (choice, rt) with projections
    s3d$points3d(rad, time, 0,
                col = "darkred",
                pch = 19,
                cex = 2)
    
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
                  "red",
                  "darkred"), 
           pch = c(19, 19, NA, 19),
           lty = c(NA, NA, 1, NA),
           cex = 1.4,      # Increased legend text size
           pt.cex = 1.5)   # Increased legend symbol size
}