
# Select arbitrary parameter values
    drift <- 1;  theta <- pi/4; tzero <- 0.1; boundary <- 4


# Test full CDF
start_time <- Sys.time()
p1a <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="grid")
end_time <- Sys.time()
print(end_time - start_time)
start_time2 <- Sys.time()
p2 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="monte_carlo", n_points=2000)
end_time2 <- Sys.time()
print(end_time2 - start_time2)
start_time3 <- Sys.time()
p3 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="monte_carlo", n_points=2000, show=TRUE)
end_time3 <- Sys.time()
print(end_time3 - start_time3)
    
    # Test marginal CDFs
    p_rt <- pCDDM(c(2, 1), drift, theta, tzero, boundary, type="RT")
    p_rad <- pCDDM(c(pi, 0.1), drift, theta, tzero, boundary, type="rad")
    
    # Print results
    cat("Grid method:", p1, "\n")
    cat("Monte Carlo method:", p2, "\n")
    cat("RT marginal CDF:", p_rt, "\n")
    cat("Rad marginal CDF:", p_rad, "\n")


# Example usage:
if(!exists("test")) { test <- TRUE }
if(test) {
    # Test parameters
    drift <- 1;  theta <- pi/4; tzero <- 0.1; boundary <- 4

    # Test full CDF
    start_time <- Sys.time()
    p1 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="grid")
    end_time <- Sys.time()
    print(end_time - start_time)
    start_time2 <- Sys.time()
    p2 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="monte_carlo", n_points=2000)
    end_time2 <- Sys.time()
    print(end_time2 - start_time2)
    start_time3 <- Sys.time()
    p3 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="monte_carlo", n_points=2000, show=TRUE)
    end_time3 <- Sys.time()
    print(end_time3 - start_time3)
    
    # Test marginal CDFs
    p_rt <- pCDDM_RT(2, drift, theta, tzero, boundary)
    p_rad <- pCDDM_rad(pi, drift, theta, tzero, boundary)
    
    # Print results
    cat("Grid method:", p1, "\n")
    cat("Monte Carlo method:", p2, "\n")
    cat("RT marginal CDF:", p_rt, "\n")
    cat("Rad marginal CDF:", p_rad, "\n")
}

# Test Monte Carlo convergence and timing
if(!exists("test")) { test <- TRUE }
if(test) {
    # Test parameters
    drift <- 1
    theta <- pi/4
    tzero <- 0.1
    boundary <- 4
    data <- c(pi, 2)
    
    # Different n_points to test
    n_points_seq <- seq(500, 10000, 500)  # More evenly distributed sequence
    n_repetitions <- 100  # Number of times to repeat each test
    
    # Store results
    results <- data.frame()
    
    # Run tests
    for(n in n_points_seq) {
        cat("\nTesting n_points =", n, "\n")
        
        # Store results for this n_points
        values <- numeric(n_repetitions)
        times <- numeric(n_repetitions)
        
        for(i in 1:n_repetitions) {
            start_time <- Sys.time()
            values[i] <- pCDDM(data, drift, theta, tzero, boundary, 
                              method="monte_carlo", n_points=n)
            end_time <- Sys.time()
            times[i] <- as.numeric(end_time - start_time, units="secs")
        }
        
        # Compute statistics
        results <- rbind(results, data.frame(
            n_points = n,
            mean_value = mean(values),
            sd_value = sd(values),
            mean_time = mean(times),
            sd_time = sd(times)
        ))
    }
 
    # Print results
    print("Monte Carlo Results:")
    print(results)
    
    # Create PDF file
    pdf("tests/pCDDM_testing2.pdf", width=10, height=5)
    
    # Set up side-by-side plots
    par(mfrow=c(1,2))
    
    # Plot 1: Convergence of values
    plot(results$n_points, results$mean_value, type="l",
         xlab="Number of Points", ylab="CDF Value",
         main="Monte Carlo Convergence", 
         ylim=range(c(results$mean_value + c(-results$sd_value, results$sd_value))))
    # Add error bands
    polygon(c(results$n_points, rev(results$n_points)),
           c(results$mean_value + results$sd_value,
             rev(results$mean_value - results$sd_value)),
           col=rgb(0,0,0,0.2), border=NA)
    lines(results$n_points, results$mean_value, lwd=2)
    
    # Plot 2: Execution time
    plot(results$n_points, results$mean_time, type="l",
         xlab="Number of Points", ylab="Time (seconds)",
         main="Execution Time")
    # Add error bands
    polygon(c(results$n_points, rev(results$n_points)),
           c(results$mean_time + results$sd_time,
             rev(results$mean_time - results$sd_time)),
           col=rgb(0,0,0,0.2), border=NA)
    lines(results$n_points, results$mean_time, lwd=2)
    
    # Reset plot parameters and close PDF
    par(mfrow=c(1,1))
    dev.off()
}

