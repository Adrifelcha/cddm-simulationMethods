source("code/cddm/sim_randomWalk.R")


# Define parameter sets to test
param_sets <- list( easy = list(par = list(
                                 mu1 = 2.0,  mu2 = 2.0,  # (strong drift)
                                 boundary = 5, tzero = 0.1
                               )),
                    hard = list(par = list(
                                 mu1 = 0.5,  mu2 = 0.5,  # (weak drift)
                                 boundary = 5, tzero = 0.1
                               )))
# Different trial sizes to test
trial_sizes <- c(100, 500, 1000)
# Number of replications
n_reps <- 10

# Function to run one benchmark test
run_single_test <- function(params, n_trials) {
    # Create a new list with n and the existing parameters
    full_params <- list(
        n = n_trials,
        par = params$par
    )
    
    start_time <- Sys.time()
    result <- do.call(sample.RW.cddm, full_params)
    end_time <- Sys.time()
    
    # Get final states
    final_coords <- getFinalState(result$random.walk)
    
    # Calculate metrics
    execution_time <- as.numeric(difftime(end_time, start_time, units="secs"))
    completion <- mean(!is.na(result$bivariate.data$RT))
    
    # Calculate distance from boundary for each trial
    distances_from_boundary <- sqrt(rowSums(final_coords^2))
    circumference_precision <- mean(abs(distances_from_boundary - params$par$boundary))
    
    return(list(
        execution_time = execution_time,
        completion = completion,
        circumference_precision = circumference_precision,
        mean_rt = mean(result$bivariate.data$RT, na.rm=TRUE),
        sd_rt = sd(result$bivariate.data$RT, na.rm=TRUE)
    ))
}

# Run full benchmark
results <- data.frame()

for(param_name in names(param_sets)) {
    for(n_trials in trial_sizes) {
        for(rep in 1:n_reps) {
            cat(sprintf("\nRunning %s, trials=%d, rep=%d", param_name, n_trials, rep))
            
            bench <- run_single_test(param_sets[[param_name]], n_trials)
            
            results <- rbind(results, data.frame(
                param_set = param_name,
                n_trials = n_trials,
                replication = rep,
                execution_time = bench$execution_time,
                completion = bench$completion,
                circumference_precision = bench$circumference_precision,
                mean_rt = bench$mean_rt,
                sd_rt = bench$sd_rt
            ))
        }
    }
}

# Create ordered groups that will be used across all plots
plot_order <- paste(rep(trial_sizes, each=2), 
                   rep(c("easy", "hard"), length(trial_sizes)))
results$group <- factor(paste(results$n_trials, results$param_set),
                       levels = plot_order)

# Analyze results
summary_stats <- aggregate(
    cbind(execution_time, completion, circumference_precision) 
    ~ param_set + n_trials, 
    data = results,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
)

# Create PDF for jittered points
pdf("tests/test_RW-old.pdf", width=10, height=8)

# Set up plotting parameters
par(mfrow=c(2,1), 
    mar=c(5,5,3,2),    
    cex.main=1.8,      
    cex.lab=1.2,
    cex.axis=1.1)

x_positions <- as.numeric(results$group)

# Calculate means for each group
exec_means <- tapply(results$execution_time, results$group, mean)
circ_means <- tapply(results$circumference_precision, results$group, mean)

# Create color vector based on parameter set
point_colors <- ifelse(results$param_set == "easy", "lightblue", "lightgreen")

# Plot execution time
plot(jitter(x_positions, amount=0.2), results$execution_time,
     col=point_colors,
     pch=19,
     main="Execution Time Distribution",
     xlab="Number of Trials",
     ylab="Time (seconds)",
     xaxt="n",
     yaxt="n",  # Suppress default y-axis
     xlim=c(0.5, 6.5))
# Add custom y-axis with 0 decimal places and horizontal labels
axis(2, at=axTicks(2), labels=sprintf("%.0f", axTicks(2)), las=2)
# Add mean lines and text
segments(1:length(exec_means) - 0.25, exec_means,
         1:length(exec_means) + 0.25, exec_means,
         lwd=2, col="black")
text(1:length(exec_means) + 0.3, exec_means, 
     sprintf("%.4f", exec_means),  # Keep 4 decimal places for mean values
     adj=0, cex=0.8)
# Add custom x-axis labels
axis(1, at=seq_along(levels(results$group)), labels=levels(results$group))
# Add legend to top plot with bold title and parameter values
legend("topleft", 
       legend=c(
           expression(paste("Easy (", delta, " = 2.8, ", theta, " = ", pi/4, ", ", tau, " = 0.1)")),
           expression(paste("Hard (", delta, " = 0.7)"))
       ),
       fill=c("lightblue", "lightgreen"),
       title="Parameter Set",
       title.font=2,
       cex=1.2,
       bty="n")

# Plot circumference precision
plot(jitter(x_positions, amount=0.2), results$circumference_precision,
     col=point_colors,
     pch=19,
     main="Circumference Precision Distribution",
     xlab="Number of Trials",
     ylab="Average Distance from Boundary",
     xaxt="n",
     yaxt="n",  # Suppress y-axis
     xlim=c(0.5, 6.5))

# Calculate y-axis values
y_range <- range(results$circumference_precision)
y_mid <- mean(y_range)
y_values <- c(y_range[1], y_mid, y_range[2])

# Add custom y-axis with 3 values in scientific notation
axis(2, at=y_values, labels=sprintf("%.2e", y_values), las=2)

# Add mean lines and text
segments(1:length(circ_means) - 0.25, circ_means,
         1:length(circ_means) + 0.25, circ_means,
         lwd=2, col="black")
text(1:length(circ_means) + 0.3, circ_means, 
     sprintf("%.4f", circ_means),  # Keep 4 decimal places for mean values
     adj=0, cex=0.8)
# Add custom x-axis labels
axis(1, at=seq_along(levels(results$group)), labels=levels(results$group))

dev.off()


# Print summary
print(summary_stats)