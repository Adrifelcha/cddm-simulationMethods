#########################################################################
# S I N G L E     A L G O R I T H M     Q U I C K     T E S T ##########
# This script is used to quickly test any individual sampling algorithm
# contained in this repo.
# The *quick test* uses few trial sizes and only two parameter sets:
# one with strong drift and one with weak drift.
# This script is designed to be run from the command line with a 
# single argument specifying the algorithm to test.
#########################################################################
method_tested <- "Rejection"
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
source(here("code", "general_functions", "eCDF.R"))

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
cat("Fast drift length:", sqrt(sum(rep(param_sets$fast$par$mu1, 2)^2)), "\n")
cat("Slow drift length:", sqrt(sum(rep(param_sets$slow$par$mu1, 2)^2)), "\n")
cat("Boundary:", param_sets$fast$par$boundary, "\n")
cat("Tzero:", param_sets$fast$par$tzero, "\n")
cat("Theta:", atan2(param_sets$fast$par$mu2, param_sets$fast$par$mu1), "\n\n")
# Different trial sizes to test
cat("Setting trial sizes to test...\n")
trial_sizes <- c(80, 150, 300, 500, 1000)
cat("Trial sizes:", trial_sizes, "\n\n")
# Number of replications
n_reps <- 50
cat("Setting number of replications:", n_reps, "\n\n")

#############################################################
#### R U N N I N G     T E S T S ############################
#############################################################
cat("Running tests...\n")
results <- data.frame()

# Create 3D arrays for each parameter set to store bivariate data
data_arrays <- list()
for(param_name in names(param_sets)) {
    data_arrays[[param_name]] <- list()
    for(n_trial in trial_sizes) {
        data_arrays[[param_name]][[as.character(n_trial)]] <- array(
            NA, 
            dim = c(n_trial, 2, n_reps),
            dimnames = list(
                NULL,
                c("Choice", "RT"),
                paste0("rep", 1:n_reps)
            )
        )
    }
}

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
            
            # Store bivariate data in 3D array
            data_arrays[[param_name]][[as.character(n_trials)]][,1,rep] <- bench$data$Choice
            data_arrays[[param_name]][[as.character(n_trials)]][,2,rep] <- bench$data$RT
            
            # Process summary results            
            output <- data.frame(param_set = param_name,
                               n_trials = n_trials,
                               replication = rep,
                               execution_time = bench$execution_time,
                               completion = bench$completion,
                               prop_negative_rt = bench$prop_negative_rt,                                        
                               mean_rt = bench$mean_rt,
                               mean_angle = bench$mean_angle,
                               angular_error = bench$angular_error,
                               stringsAsFactors = FALSE)           
            # Add circumference precision if method is RandomWalk
            if(method_tested == "RandomWalk") {
                output$circumference_precision <- bench$circumference_precision
            }            
            results <- rbind(results, output)            
        }
    }
}

#############################################################
#### P R O C E S S     R E S U L T S  #######################
#############################################################
# Add CDFs to data arrays
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
data_arrays <- add_cdfs_to_arrays(data_arrays, param_sets)
results_cdfs <- calculate_cdf_metrics(data_arrays)

# Prepare for plotting
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_order <- paste(rep(trial_sizes, each=2), 
                    rep(c("fast", "slow"), length(trial_sizes)))           
results$group <- factor(paste(results$n_trials, results$param_set),
                       levels = plot_order)
results_cdfs$group <- factor(paste(results_cdfs$trial_size, results_cdfs$param_set),
                       levels = plot_order)

# Get summary statistics and metrics
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
summary_stats <- aggregate(
    cbind(execution_time, angular_error, prop_negative_rt, completion) 
    ~ param_set + n_trials, 
    data = results,
    FUN = function(x) c(mean = mean(x))
)
colnames(summary_stats) <- c("param_set", "n_trials", "mean_exec_time", "rad_difference", "prop_neg_rt", "completion_rate")


# Save both summary results and data arrays
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
filename <- sprintf(here("results", "quickTest_%s_%s.RData"), method_tested, format(Sys.Date(), "%Y%m%d"))
save(results, data_arrays, file = filename)
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
    plot_circumference_precision(results, param_sets, trial_sizes, method_tested, filename = figname_circPrecision)
}

# Print summary
cat("\nSummary of results:\n")
print(summary_stats)

cat("\n Figures created:", figname_results, "\n\n")

cat("\n\nDone!\n\n")


# Usage example:
#figure_cdf <- sprintf(here("results", "quickTest_%s_%s_cdfs.pdf"), method_tested, format(Sys.Date(), "%Y%m%d"))
#pdf(figure_cdf, width=10, height=8)
# plot_cdfs(data_arrays)
#dev.off()


# Generate the plots
filename_prefix <- here("results", 
                       sprintf("quickTest_%s_%s_metrics", 
                              method_tested, 
                              format(Sys.Date(), "%Y%m%d")))
plot_metrics_by_paramset(results_cdfs, results, param_sets, trial_sizes, n_reps, 
                        filename_prefix)