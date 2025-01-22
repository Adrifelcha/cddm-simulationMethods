#############################################################
method_tested <- "RandomWalk"
# Possible methods:
# 1) "Metropolis"
# 2) "RandomWalk"
# 3) "inverseCDF"
# 4) "Rejection"
#############################################################
superCalled <- TRUE
source("code/cddm/sim_randomWalk.R")
# Load all R files from code/cddm folder
r_files <- list.files(path = "code/cddm", 
                      pattern = "\\.R$", 
                      full.names = TRUE)
for(file in r_files) {
    source(file)
}

#############################################################
#### S E T T I N G S ########################################
#############################################################

# Define comprehensive parameter sets to test
param_sets <- list(
    very_easy = list(par = list(
        mu1 = 3.0, mu2 = 3.0,     # very strong drift
        boundary = 5, tzero = 0.1
    )),
    easy = list(par = list(
        mu1 = 2.0, mu2 = 2.0,     # strong drift
        boundary = 5, tzero = 0.1
    )),
    medium = list(par = list(
        mu1 = 1.0, mu2 = 1.0,     # moderate drift
        boundary = 5, tzero = 0.1
    )),
    hard = list(par = list(
        mu1 = 0.5, mu2 = 0.5,     # weak drift
        boundary = 5, tzero = 0.1
    )),
    very_hard = list(par = list(
        mu1 = 0.2, mu2 = 0.2,     # very weak drift
        boundary = 5, tzero = 0.1
    ))
)

# More comprehensive trial sizes
trial_sizes <- c(50, 150, 300, 500)

# Set to 1 for initial testing
n_reps <- 1

#############################################################
#### T E S T I N G     F U N C T I O N S ####################
#############################################################

# Function to run one benchmark test
run_single_test <- function(params, n_trials, method_tested) {
    # Create a new list with n and the existing parameters
    full_params <- list(
        n = n_trials,
        par = params$par
    )
    
    start_time <- Sys.time()
    
    # Select the appropriate sampling function based on method_tested
    result <- switch(method_tested,
        "RandomWalk" = do.call(sample.RW.cddm, full_params),
        "Metropolis" = do.call(sample.Metropolis.cddm, full_params),
        "inverseCDF" = do.call(sample.invCDF.cddm, full_params),
        "Rejection" = do.call(sample.Reject.cddm, full_params),
        stop(paste("Unknown method:", method_tested))
    )
    
    end_time <- Sys.time()
    
    # Get final states
    final_coords <- getFinalState(result$random.walk)
    
    # Calculate metrics
    execution_time <- as.numeric(difftime(end_time, start_time, units="secs"))
    completion <- mean(!is.na(result$bivariate.data$RT))
    
    # Calculate distance from boundary for each trial
    distances_from_boundary <- sqrt(rowSums(final_coords^2))
    circumference_precision <- mean(abs(distances_from_boundary - params$par$boundary))
    
    # Additional metrics for full test
    return(list(
        execution_time = execution_time,
        completion = completion,
        circumference_precision = circumference_precision,
        mean_rt = mean(result$bivariate.data$RT, na.rm=TRUE),
        sd_rt = sd(result$bivariate.data$RT, na.rm=TRUE),
        min_rt = min(result$bivariate.data$RT, na.rm=TRUE),
        max_rt = max(result$bivariate.data$RT, na.rm=TRUE),
        median_rt = median(result$bivariate.data$RT, na.rm=TRUE),
        na_count = sum(is.na(result$bivariate.data$RT))
    ))
}

#############################################################
#### R U N N I N G     T E S T S ############################
#############################################################

# Initialize results dataframe with additional columns
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
            
            bench <- run_single_test(param_sets[[param_name]], n_trials, method_tested)
            
            results <- rbind(results, data.frame(
                param_set = param_name,
                n_trials = n_trials,
                replication = rep,
                execution_time = bench$execution_time,
                completion = bench$completion,
                circumference_precision = bench$circumference_precision,
                mean_rt = bench$mean_rt,
                sd_rt = bench$sd_rt,
                min_rt = bench$min_rt,
                max_rt = bench$max_rt,
                median_rt = bench$median_rt,
                na_count = bench$na_count
            ))
        }
    }
}
cat("\nTesting completed!\n")

# Create ordered groups
plot_order <- paste(rep(trial_sizes, each=length(param_sets)), 
                   rep(names(param_sets), length(trial_sizes)))
results$group <- factor(paste(results$n_trials, results$param_set),
                       levels = plot_order)

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

# Create ordered groups that will be used across all plots
plot_order <- paste(rep(trial_sizes, each=2), 
                   rep(c("easy", "hard"), length(trial_sizes)))


#############################################################
#### P L O T T I N G     R E S U L T S ######################
#############################################################

# Create PDF with date in filename
filename <- sprintf("tests/test_%s_%s.pdf", method_tested, format(Sys.Date(), "%Y%m%d"))
pdf(filename, width=12, height=10)

# Set up 2x2 plotting layout
par(mfrow=c(2,2), mar=c(5,5,3,2),
    cex.main=1.5, cex.lab=1.2, cex.axis=1.1)

# Calculate means for plotting
exec_means <- tapply(results$execution_time, results$group, mean)
circ_means <- tapply(results$circumference_precision, results$group, mean)
rt_means <- tapply(results$mean_rt, results$group, mean)
compl_means <- tapply(results$completion, results$group, mean)

# Create color vector based on parameter set
point_colors <- rainbow(length(param_sets))[as.numeric(factor(results$param_set))]

# Function to add mean lines and labels
add_means <- function(means) {
    segments(1:length(means) - 0.25, means,
            1:length(means) + 0.25, means,
            lwd=2, col="black")
    text(1:length(means) + 0.3, means, 
         sprintf("%.4f", means), adj=0, cex=0.7)
}

# 1. Execution Time Plot
plot(jitter(as.numeric(results$group), amount=0.2), results$execution_time,
     col=point_colors, pch=19, main="Execution Time Distribution",
     xlab="Number of Trials", ylab="Time (seconds)", 
     xaxt="n", yaxt="n", xlim=c(0.5, length(levels(results$group))+0.5))
mtext(paste0(method_tested, "\n", format(Sys.Date(), "%Y-%m-%d")), 
      side=3, line=1, adj=1, cex=0.8)
axis(2, las=2)
axis(1, at=1:length(levels(results$group)), labels=levels(results$group), las=2)
add_means(exec_means)

# 2. Circumference Precision Plot
plot(jitter(as.numeric(results$group), amount=0.2), results$circumference_precision,
     col=point_colors, pch=19, main="Circumference Precision",
     xlab="Number of Trials", ylab="Average Distance from Boundary",
     xaxt="n", yaxt="n", xlim=c(0.5, length(levels(results$group))+0.5))
mtext(paste0(method_tested, "\n", format(Sys.Date(), "%Y-%m-%d")), 
      side=3, line=1, adj=1, cex=0.8)
axis(2, las=2)
axis(1, at=1:length(levels(results$group)), labels=levels(results$group), las=2)
add_means(circ_means)

# 3. Mean RT Plot
plot(jitter(as.numeric(results$group), amount=0.2), results$mean_rt,
     col=point_colors, pch=19, main="Mean Response Time",
     xlab="Number of Trials", ylab="RT (seconds)",
     xaxt="n", yaxt="n", xlim=c(0.5, length(levels(results$group))+0.5))
mtext(paste0(method_tested, "\n", format(Sys.Date(), "%Y-%m-%d")), 
      side=3, line=1, adj=1, cex=0.8)
axis(2, las=2)
axis(1, at=1:length(levels(results$group)), labels=levels(results$group), las=2)
add_means(rt_means)

# 4. Completion Rate Plot
plot(jitter(as.numeric(results$group), amount=0.2), results$completion,
     col=point_colors, pch=19, main="Completion Rate",
     xlab="Number of Trials", ylab="Completion Rate",
     xaxt="n", yaxt="n", xlim=c(0.5, length(levels(results$group))+0.5))
mtext(paste0(method_tested, "\n", format(Sys.Date(), "%Y-%m-%d")), 
      side=3, line=1, adj=1, cex=0.8)
axis(2, las=2)
axis(1, at=1:length(levels(results$group)), labels=levels(results$group), las=2)
add_means(compl_means)

# Add legend to the first plot
legend("topleft", 
       legend=names(param_sets),
       fill=rainbow(length(param_sets)),
       title="Parameter Set",
       title.font=2,
       cex=0.8,
       bty="n")

dev.off()

# Save results to RData file
save(results, file=sprintf("tests/results_%s_%s.RData", 
                          method_tested, format(Sys.Date(), "%Y%m%d")))

# Print summary statistics
summary_stats <- aggregate(
    cbind(execution_time, completion, circumference_precision, mean_rt) 
    ~ param_set + n_trials, 
    data = results,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
)
print(summary_stats)