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
n_reps <- 200


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
    if (method_tested == "RandomWalk") {
        result <- do.call(sample.RW.cddm, full_params)
    } else if (method_tested == "Metropolis") {
        result <- do.call(sample.Metropolis.cddm, full_params)
    } else if (method_tested == "inverseCDF") {
        result <- do.call(sample.invCDF.cddm, full_params)
    } else if (method_tested == "Rejection") {
        result <- do.call(sample.Reject.cddm, full_params)
    } else {
        stop(paste("Unknown method:", method_tested))
    }
    
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


#############################################################
#### R U N N I N G     T E S T S ############################
#############################################################

# Run full benchmark
results <- data.frame()

for(param_name in names(param_sets)) {
    for(n_trials in trial_sizes) {
        for(rep in 1:n_reps) {
            cat(sprintf("\nRunning %s, trials=%d, rep=%d", param_name, n_trials, rep))
            
            bench <- run_single_test(param_sets[[param_name]], n_trials, method_tested)
            
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
filename <- sprintf("tests/test_%s_%s.pdf", method_tested, format(Sys.Date(), "%Y%m%d"))
pdf(filename, width=10, height=8)
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

# Create color vector based on parameter set
point_colors <- ifelse(results$param_set == "easy", "lightblue", "lightgreen")

# Plot execution time
plot(jitter(x_positions, amount=0.2), results$execution_time,
    col=point_colors, pch=19, main="Execution Time Distribution",
    xlab="Number of Trials", ylab="Time (seconds)", xaxt="n",
    yaxt="n", xlim=c(0.5, 6.5))
mtext(paste0(method_tested, "\n", format(Sys.Date(), "%Y-%m-%d")), 
      side=3, line=1, adj=1, cex=0.8)
axis(2, at=axTicks(2), labels=sprintf("%.0f", axTicks(2)), las=2)
# Add mean lines and text
segments(1:length(exec_means) - 0.25, exec_means,
        1:length(exec_means) + 0.25, exec_means,
        lwd=2, col="black")
text(1:length(exec_means) + 0.3, exec_means, sprintf("%.4f", exec_means), 
    adj=0, cex=0.8)
axis(1, at=seq_along(levels(results$group)), labels=levels(results$group))
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
    col=point_colors, pch=19, main="Circumference Precision Distribution",
    xlab="Number of Trials", ylab="Average Distance from Boundary",
    xaxt="n", yaxt="n", xlim=c(0.5, 6.5))
axis(2, at=y_values, labels=sprintf("%.2e", y_values), las=2)
segments(1:length(circ_means) - 0.25, circ_means,
        1:length(circ_means) + 0.25, circ_means,
        lwd=2, col="black")
text(1:length(circ_means) + 0.3, circ_means, sprintf("%.4f", circ_means),  
    adj=0, cex=0.8)
axis(1, at=seq_along(levels(results$group)), labels=levels(results$group))
dev.off()


# Print summary
print(summary_stats)