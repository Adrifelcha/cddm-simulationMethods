########################################################
# A short script to test the functionality of pCDDM()
# pCDDM() is a function that approximates the cumulative 
# distribution function (CDF) of the pCDDM model through
# numerical integration.
########################################################
cat("Running pCDDM_tests.R\n\n")

#######################################################
# Load libraries and custom function scripts
#######################################################
cat("Loading R libraries...\n")
library("scatterplot3d")
library("plot3D")

cat("Loading custom function scripts from /code/cddm...\n\n")
source("../code/cddm/dCDDM_aux.R")
source("../code/cddm/dCDDM.R")
source("../code/cddm/pCDDM_aux.R")
source("../code/cddm/pCDDM_plots.R")
source("../code/cddm/pCDDM.R")

#######################################################
# Set test parameter values and data
#######################################################
set.seed(789)
cat("Setting arbitrary parameter values:\n")
drift <- 1
theta <- pi/4
tzero <- 0.1
boundary <- 4
cat("Drift:", drift, "\n")
cat("Theta:", theta, "\n")
cat("Tzero:", tzero, "\n")
cat("Boundary:", boundary, "\n\n")
data <- c(pi, 2)

#######################################################
# Test functionality of pCDDM()
#######################################################
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Test 1: Comparing pCDDM() results for a single pair of data\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Setting test data:\n")
cat("(radians =", data[1], ", rt =", data[2], ")\n\n")
cat("Testing pCDDM() with grid method...(this may take a minute)\n")
start_time <- Sys.time()
p1 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="grid")
end_time <- Sys.time()
cat("Grid method result:", p1, "\n")
difftime(end_time, start_time, units = "secs")

cat("\nTesting pCDDM() with Monte Carlo method...\n")
start_time2 <- Sys.time()
p2 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="monte_carlo", n_points=2000)
end_time2 <- Sys.time()
cat("Monte Carlo method result:", p2, "\n")
difftime(end_time2, start_time2, units = "secs")

cat("\nTesting pCDDM() with Monte Carlo method and showing plot...\n")
start_time3 <- Sys.time()
p3 <- pCDDM(c(pi, 2), drift, theta, tzero, boundary, method="monte_carlo", n_points=2000, show=TRUE)
end_time3 <- Sys.time()
cat("Monte Carlo (with plot) result:", p3, "\n")
difftime(end_time3, start_time3, units = "secs")
    
cat("\nTesting pCDDM() with marginal CDFs...\n")
p_rt <- pCDDM(c(2, 1), drift, theta, tzero, boundary, type="RT")
cat("RT marginal CDF:", p_rt, "\n")
p_rad <- pCDDM(c(pi, 0.1), drift, theta, tzero, boundary, type="rad")
cat("Rad marginal CDF:", p_rad, "\n\n")
    
#################################################################
# Test convergence and execution time of pCDDM()
#################################################################
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Test 2: Testing convergence and execution time of pCDDM()\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
# Different n_points to test
n_points_seq <- seq(500, 10000, 500) 
# Number of times to repeat each test
n_repetitions <- 100  
    
# Store results
results <- data.frame()
# Progress counter
total_iterations <- length(n_points_seq) * n_repetitions
current_iteration <- 0
# Run tests
for(n in n_points_seq) {
    # Store results for this n_points
    values <- numeric(n_repetitions)
    times <- numeric(n_repetitions)
    
    for(i in 1:n_repetitions) {        
        # Progress indicator
        current_iteration <- current_iteration + 1            
        cat(sprintf("\rProgress: %d/%d (%.1f%%) - Running %s n_points, rep=%s",
            current_iteration, total_iterations,
            100 * current_iteration/total_iterations,
            n, i))
        # Run pCDDM()
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
                    sd_time = sd(times)))
}
 
# Print results
cat("\n\nResults:\n")
print(results)
    
testResults <- "../results/pCDDM_testing.pdf"
cat("\nStoring results in figure", testResults, "\n")
# Create PDF file
pdf(testResults, width=10, height=5)
    # Set up side-by-side plots
    par(mfrow=c(1,2))

    # Plot 1: Convergence of values
    plot(results$n_points, results$mean_value, type="l", ylab="CDF Value",
            xlab="Number of Points", main="Monte Carlo Convergence", 
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

N <- 1000
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
cat("Test 3: Comparing pCDDM() execution time for n=", N, "observations against pnorm()\n")
cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")

c <- rnorm(N,0,1)
start1 <- Sys.time()
test1b <- pnorm(c,0,1)
end1 <- Sys.time()
cat("Using pnorm() to compute CDF of n=", N, "observations\n")
difftime(end1, start1, units = "secs")

d <- cbind(runif(N,0,2*pi), runif(N,0,2))
start2 <- Sys.time()
test2d <- pCDDM(d, 1, pi/4, 0.1, 4, method="monte_carlo")
end2 <- Sys.time()
cat("Using pCDDM() to compute CDF of n=", N, "pairs of observations\n")
difftime(end2, start2, units = "secs")
