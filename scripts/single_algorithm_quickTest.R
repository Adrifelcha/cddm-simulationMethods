#########################################################################
# S I N G L E     A L G O R I T H M     Q U I C K     T E S T ##########
# This script is used to quickly test any individual sampling algorithm
# contained in this repo.
# The *quick test* uses few trial sizes and only two parameter sets:
# one with strong drift and one with weak drift.
# This script is designed to be run from the command line with a 
# single argument specifying the algorithm to test.
#########################################################################
method_tested <- "RandomWalk"
# Possible methods:
# 1) "Metropolis"
# 2) "RandomWalk"
# 3) "inverseCDF"
# 4) "Rejection"
cat("-----------------------------------------------------------\n")
cat("\n\nSimulation algorithm to be tested:", method_tested, "\n\n")
cat("-----------------------------------------------------------\n")

#############################################################
# Load libraries and custom functions
#############################################################
cat("Loading R libraries...\n")
library("here")
library("circular")

cat("\nLoading custom function scripts from /code/cddm...\n\n")
source(here("code", "cddm", "sim_randomWalk.R"))        
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
cat("Setting parameter sets to test...\n")
param_sets <- list( fast = list(par = list(
                                 mu1 = 2.0,  mu2 = 2.0,  # (strong drift)
                                 boundary = 5, tzero = 0.1
                               )),
                    slow = list(par = list(
                                 mu1 = 0.5,  mu2 = 0.5,  # (weak drift)
                                 boundary = 5, tzero = 0.1
                               )))
cat("Easy drift length:", sqrt(sum(rep(param_sets$easy$par$mu1, 2)^2)), "\n")
cat("Hard drift length:", sqrt(sum(rep(param_sets$hard$par$mu1, 2)^2)), "\n")
cat("Boundary:", param_sets$easy$par$boundary, "\n")
cat("Tzero:", param_sets$easy$par$tzero, "\n")
cat("Theta:", atan2(param_sets$easy$par$mu2, param_sets$easy$par$mu1), "\n\n")
# Different trial sizes to test
cat("Setting trial sizes to test...\n")
trial_sizes <- c(100, 500, 1000)
cat("Trial sizes:", trial_sizes, "\n\n")
# Number of replications
n_reps <- 10
cat("Setting number of replications:", n_reps, "\n\n")

#############################################################
#### R U N N I N G     T E S T S ############################
#############################################################
cat("Running tests...\n")
results <- data.frame()
# Progress counter
total_iterations <- length(names(param_sets)) * length(trial_sizes) * n_reps
current_iteration <- 0
for(param_name in names(param_sets)) {
    for(n_trials in trial_sizes) {
        for(rep in 1:n_reps) {                       
            # Progress indicator
            current_iteration <- current_iteration + 1
            cat(sprintf("\rProgress: %d/%d (%.1f%%) - Running %s, trials=%d, rep=%d",
                current_iteration, total_iterations,
                100 * current_iteration/total_iterations,
                param_name, n_trials, rep))
            
            # Run test
            bench <- single_algorithm_test(param_sets[[param_name]], n_trials, method_tested)
            
            # Process results            
            output <- data.frame(param_set = param_name,
                                 n_trials = n_trials,
                                 replication = rep,
                                 execution_time = bench$execution_time,
                                 completion = bench$completion,
                                 prop_negative_rt = bench$prop_negative_rt,                                        
                                 mean_rt = bench$mean_rt,
                                 mean_angle = bench$mean_angle,
                                 angular_error = bench$angular_error)            
            # Add circumference precision if method is RandomWalk
            if(method_tested == "RandomWalk") {
                output$circumference_precision <- bench$circumference_precision
            }            
            results <- rbind(results, output)
        }
    }
}

#############################################################
#### G E T     R E S U L T S ################################
#############################################################
# Analyze results
summary_stats <- aggregate(
    cbind(execution_time, angular_error, prop_negative_rt, completion) 
    ~ param_set + n_trials, 
    data = results,
    FUN = function(x) c(mean = mean(x))
)
colnames(summary_stats) <- c("param_set", "n_trials", "mean_exec_time", "rad_difference", "prop_neg_rt", "completion_rate")

# Fix the order in which cells will be displayed
plot_order <- paste(rep(trial_sizes, each=2), 
                    rep(c("fast", "slow"), length(trial_sizes)))           

# Prepare results
results$group <- factor(paste(results$n_trials, results$param_set),
                       levels = plot_order)
# Save results! 
filename <- sprintf(here("results", "quickTest_%s_%s.RData"), method_tested, format(Sys.Date(), "%Y%m%d"))
save(results,file = filename)
cat("Saving results to:", filename, "\n\n")

#############################################################
#### P L O T T I N G     R E S U L T S ######################
#############################################################

# Plot algorithm performance per trial size and parameter set
# (1) execution time, (2) mean angle - theta error, (3) mean RT
figname_results <- sprintf(here("results", "quickTest_%s_%s.pdf"), method_tested, format(Sys.Date(), "%Y%m%d"))
plot_algorithm_performance(results, param_sets, trial_sizes, n_reps, method_tested, filename = figname_results)    


if(method_tested == "RandomWalk") {
    figname_circPrecision <- sprintf(here("results", "quickTest_%s_%s_circPrecision.pdf"), method_tested, format(Sys.Date(), "%Y%m%d"))
    pdf(figname_circPrecision, width=10, height=8)
   # Plot circumference precision
        y_range <- range(results$circumference_precision)
        y_mid <- mean(y_range)
        y_values <- c(y_range[1], y_mid, y_range[2])
        circ_means <- tapply(results$circumference_precision, results$group, mean)
        # Set up plotting parameters
        par(mfrow=c(1,1), mar=c(5,5,3,2),
        cex.main=1.8, cex.lab=1.2, cex.axis=1.1)

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
}

# Print summary
cat("\nSummary of results:\n")
print(summary_stats)

cat("\n Figures created:", figname_results, "\n\n")

cat("\n\nDone!\n\n")