compare_grid_sizes <- function(drift=2, theta=pi/4, tzero=0.3, boundary=1) {
    # Test points to evaluate
    test_points <- data.frame(
        rad = c(pi/2, pi, 3*pi/2, 2*pi),  # Different angular bounds
        time = c(1.0, 1.5, 2.0, 2.5)      # Different time bounds
    )
    
    # Get "true" probabilities using very fine Monte Carlo integration
    true_probs <- pCDDM(as.matrix(test_points), drift, theta, tzero, boundary,
                        method="monte_carlo", n_points=100000)
    
    # Grid sizes to test
    grid_sizes <- c(10, 20, 50, 100, 200, 500)
    
    # Initialize results matrices
    results_regular <- matrix(NA, nrow=length(grid_sizes), ncol=3)
    results_adaptive <- matrix(NA, nrow=length(grid_sizes), ncol=3)
    colnames(results_regular) <- colnames(results_adaptive) <- 
        c("n_points", "mean_abs_error", "computation_time")
    
    # Test each grid size
    for(i in seq_along(grid_sizes)) {
        n <- grid_sizes[i]
        
        # Regular grid
        time_regular <- system.time({
            probs_regular <- pCDDM(as.matrix(test_points), drift, theta, tzero, boundary,
                                 method="grid", n_points=n)
        })
        error_regular <- mean(abs(probs_regular - true_probs))
        
        # Adaptive grid
        time_adaptive <- system.time({
            probs_adaptive <- pCDDM(as.matrix(test_points), drift, theta, tzero, boundary,
                                  method="grid", n_points=n)
        })
        error_adaptive <- mean(abs(probs_adaptive - true_probs))
        
        # Store results
        results_regular[i,] <- c(n, error_regular, time_regular["elapsed"])
        results_adaptive[i,] <- c(n, error_adaptive, time_adaptive["elapsed"])
    }
    
    # Create comparison plots
    par(mfrow=c(1,2))
    
    # Error plot
    plot(log10(results_regular[,"n_points"]), log10(results_regular[,"mean_abs_error"]),
         type="b", col="blue", pch=16,
         xlab="log10(Grid Points per Dimension)",
         ylab="log10(Mean Absolute Error)",
         main="Accuracy Comparison")
    lines(log10(results_adaptive[,"n_points"]), log10(results_adaptive[,"mean_abs_error"]),
          type="b", col="red", pch=16)
    legend("topright", c("Regular Grid", "Adaptive Grid"),
           col=c("blue", "red"), pch=16, lty=1)
    
    # Time plot
    plot(log10(results_regular[,"n_points"]), results_regular[,"computation_time"],
         type="b", col="blue", pch=16,
         xlab="log10(Grid Points per Dimension)",
         ylab="Computation Time (seconds)",
         main="Computational Efficiency")
    lines(log10(results_adaptive[,"n_points"]), results_adaptive[,"computation_time"],
          type="b", col="red", pch=16)
    legend("topleft", c("Regular Grid", "Adaptive Grid"),
           col=c("blue", "red"), pch=16, lty=1)
    
    # Return results as data frames
    list(
        regular = as.data.frame(results_regular),
        adaptive = as.data.frame(results_adaptive)
    )
}

# Run comparison
#results <- compare_grid_sizes()
#print("Regular Grid Results:")
#print(results$regular)
#print("\nAdaptive Grid Results:")
#print(results$adaptive) 