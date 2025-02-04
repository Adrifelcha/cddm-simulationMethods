#############################################################
method_tested <- "RandomWalk"
# Possible methods:
# 1) "Metropolis"
# 2) "RandomWalk"
# 3) "inverseCDF"
# 4) "Rejection"
cat("Simulation algorithm to be tested:", method_tested, "\n\n")
#############################################################

library("here")
superCalled <- TRUE
source(here("code", "cddm", "sim_randomWalk.R"))        
# Load all R files from code/cddm folder
r_files <- list.files(path = here("code", "cddm"), 
                      pattern = "\\.R$", 
                      full.names = TRUE)
for(file in r_files) {
    source(file)
}

#############################################################
#### S E T T I N G S ########################################
#############################################################

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


#############################################################
#### R U N N I N G     T E S T S ############################
#############################################################
results <- data.frame()
# Progress counter
total_iterations <- length(names(param_sets)) * length(trial_sizes) * n_reps
current_iteration <- 0
for(param_name in names(param_sets)) {
    for(n_trials in trial_sizes) {
        for(rep in 1:n_reps) {
            current_iteration <- current_iteration + 1
            
            # Progress indicator
            cat(sprintf("\rProgress: %d/%d (%.1f%%) - Running %s, trials=%d, rep=%d",
                current_iteration, total_iterations,
                100 * current_iteration/total_iterations,
                param_name, n_trials, rep))
            
            # Run test
            bench <- single_algorithm_test(param_sets[[param_name]], n_trials, method_tested)
            # Process results
            results <- rbind(results, 
                             data.frame(param_set = param_name,
                                        n_trials = n_trials,
                                        replication = rep,
                                        execution_time = bench$execution_time,
                                        completion = bench$completion,
                                        circumference_precision = bench$circumference_precision,
                                        mean_rt = bench$mean_rt,
                                        sd_rt = bench$sd_rt))
        }
    }
}
# Fix the order in which cells will be displayed
plot_order <- paste(rep(trial_sizes, each=2), 
                    rep(c("easy", "hard"), length(trial_sizes)))
# Prepare results
results$group <- factor(paste(results$n_trials, results$param_set),
                       levels = plot_order)
# Save results!
filename <- sprintf("tests/test_%s_%s.RData", method_tested, format(Sys.Date(), "%Y%m%d"))
save(results,file = filename)
#############################################################
#### G E T     R E S U L T S ################################
#############################################################
# Analyze results
summary_stats <- aggregate(
    cbind(execution_time, completion, circumference_precision) 
    ~ param_set + n_trials, 
    data = results,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
)

#############################################################
#### P L O T T I N G     R E S U L T S ######################
#############################################################
# Open pdf
figname <- sprintf("tests/test_%s_%s.pdf", method_tested, format(Sys.Date(), "%Y%m%d"))
pdf(figname, width=10, height=8)

# Set up plotting parameters
par(mfrow=c(2,1), mar=c(5,5,3,2),    
    cex.main=1.8, cex.lab=1.2, cex.axis=1.1)
x_positions <- as.numeric(results$group)
y_range <- range(results$circumference_precision)
y_mid <- mean(y_range)
y_values <- c(y_range[1], y_mid, y_range[2])

# Calculate means for each group
exec_means <- tapply(results$execution_time, results$group, mean)
circ_means <- tapply(results$circumference_precision, results$group, mean)

# Create color vectors with transparency
transparent_blue <- rgb(173/255, 216/255, 230/255, alpha=0.3)  # lightblue with 30% opacity
transparent_green <- rgb(144/255, 238/255, 144/255, alpha=0.3)  # lightgreen with 30% opacity
point_colors <- ifelse(results$param_set == "easy", transparent_blue, transparent_green)

# Plot execution time
plot(jitter(x_positions, amount=0.2), results$execution_time, col=point_colors, 
     pch=19, main="Execution Time Distribution", xaxt="n", yaxt="n", 
     xlab="Number of Trials", ylab="Time (seconds)", xlim=c(0.5, 6.5))
mtext(paste0(method_tested, "\n", format(Sys.Date(), "%Y-%m-%d")), 
      side=3, line=1, adj=1, cex=0.8)
axis(2, at=axTicks(2), labels=sprintf("%.0f", axTicks(2)), las=2)
# Signal and show mean values
segments(1:length(exec_means) - 0.25, exec_means,
         1:length(exec_means) + 0.25, exec_means, lwd=2, col="black")
text(1:length(exec_means) + 0.3, exec_means, sprintf("%.4f", exec_means),
     adj=0, cex=0.8)
axis(1, at=seq_along(levels(results$group)), labels=levels(results$group))
# Add legend to top plot with bold title and parameter values
legend("topleft", legend=c(
        expression(paste("Easy (", delta, " = 2.8, ", theta, " = ", pi/4, ", ", tau, " = 0.1)")),
        expression(paste("Hard (", delta, " = 0.7)"))), 
        fill=c(transparent_blue, transparent_green), title="Parameter Set",
        cex=1.2, bty="n")

# Plot circumference precision
plot(jitter(x_positions, amount=0.2), results$circumference_precision,
     col=point_colors, pch=19, main="Circumference Precision Distribution",
     xlab="Number of Trials", ylab="Average Distance from Boundary", xaxt="n",
     yaxt="n", xlim=c(0.5, 6.5))
axis(2, at=y_values, labels=sprintf("%.2e", y_values), las=2)
segments(1:length(circ_means) - 0.25, circ_means,
         1:length(circ_means) + 0.25, circ_means,
         lwd=2, col="black")
text(1:length(circ_means) + 0.3, circ_means, sprintf("%.4f", circ_means),
     adj=0, cex=0.8)
axis(1, at=seq_along(levels(results$group)), labels=levels(results$group))

# Close pdf
dev.off()

# Print summary
print(summary_stats)