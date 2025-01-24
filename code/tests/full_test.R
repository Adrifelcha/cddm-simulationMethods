#############################################################
method_tested <- "RandomWalk"
# Possible methods:
# 1) "Metropolis"
# 2) "RandomWalk"
# 3) "inverseCDF"
# 4) "Rejection"
#############################################################
library(circular)
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
        mu1 = 3.0, mu2 = 3.0,     # very strong drift (δ = 4.24)
        boundary = 5, tzero = 0.1
    )),
    easy = list(par = list(
        mu1 = 2.0, mu2 = 2.0,     # strong drift (δ = 2.83)
        boundary = 5, tzero = 0.1
    )),
    medium = list(par = list(
        mu1 = 1.0, mu2 = 1.0,     # moderate drift (δ = 1.41)
        boundary = 5, tzero = 0.1
    )),
    hard = list(par = list(
        mu1 = 0.5, mu2 = 0.5,     # weak drift (δ = 0.71)
        boundary = 5, tzero = 0.1
    )),
    very_hard = list(par = list(
        mu1 = 0.2, mu2 = 0.2,     # very weak drift (δ = 0.28)
        boundary = 5, tzero = 0.1
    ))
)

# More comprehensive trial sizes
trial_sizes <- c(50, 150, 300, 500)

# Number of replications
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
    
    # Convert angles to circular data type and calculate circular mean
    angles <- circular::circular(result$bivariate.data$Choice, type="angles", units="radians")
    mean_angle <- circular::mean.circular(angles)
    
    # Calculate theoretical angle (theta) from mu parameters
    theoretical_theta <- circular::circular(atan2(params$par$mu2, params$par$mu1))
    
    # Calculate angular error as the difference between mean_angle and theoretical_theta
    angular_error <- abs(as.numeric(mean_angle - theoretical_theta))
    
    # Calculate proportion of negative RTs (excluding NAs)
    prop_negative_rt <- mean(result$bivariate.data$RT < 0, na.rm=TRUE)
    
    # Calculate all metrics
    return(list(
        execution_time = execution_time,
        completion = completion,
        angular_error = angular_error,
        mean_rt = mean(result$bivariate.data$RT, na.rm=TRUE),
        sd_rt = sd(result$bivariate.data$RT, na.rm=TRUE),
        mean_angle = abs(as.numeric(mean_angle - theoretical_theta)),
        sd_angle = sd.circular(angles),
        prop_negative_rt = prop_negative_rt
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
                angular_error = bench$angular_error,
                mean_rt = bench$mean_rt,
                sd_rt = bench$sd_rt,
                mean_angle = bench$mean_angle,
                sd_angle = bench$sd_angle,
                prop_negative_rt = bench$prop_negative_rt
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
    cbind(execution_time, completion, angular_error) 
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
filename <- sprintf("tests/testFull_%s_%s.pdf", method_tested, format(Sys.Date(), "%Y%m%d"))
pdf(filename, width=12, height=10)

# Set up 2x2 plotting layout
par(mfrow=c(2,2), oma=c(2,2,2,2), mar=c(4,4,3,1), mgp=c(2,0.7,0),
    cex.main=1.5, cex.lab=1.2, cex.axis=1.1)

# Calculate means for plotting
exec_means <- tapply(results$execution_time, results$group, mean)
circ_means <- tapply(results$angular_error, results$group, mean)
rt_means <- tapply(results$mean_rt, results$group, mean)
compl_means <- tapply(results$completion, results$group, mean)

# Create color vector based on parameter set
point_colors <- rainbow(length(param_sets))[as.numeric(factor(results$param_set))]

# Function to add mean lines and labels
add_means <- function(means) {
    segments(1:length(means) - 0.25, means,
             1:length(means) + 0.25, means,
             lwd=3, col="black")
}

# Function to add background stripes with different colors
add_background_stripes <- function() {
    # Define four light colors for the stripes
    stripe_colors <- c("lavender", "azure", "honeydew", "cornsilk")
    
    # Add each colored stripe
    for(i in 1:4) {
        rect(xleft = seq(0.5 + (i-1)*5, max(as.numeric(results$group)) + 0.5, by=20) - 0.5,
             ybottom = par("usr")[3],
             xright = seq(5.5 + (i-1)*5, max(as.numeric(results$group)) + 0.5, by=20) - 0.5,
             ytop = par("usr")[4],
             col = stripe_colors[i],
             border = NA)
    }
}

# 1. Execution Time Plot
plot(1, type="n",
     xlim=c(0.5, length(levels(results$group))+0.5),
     ylim=range(results$execution_time),
     xlab="", ylab="Time (seconds)",
     main="Execution Time Distribution",
     xaxt="n", yaxt="n")
mtext(paste0(method_tested, "\n", format(Sys.Date(), "%Y-%m-%d")), 
      side=3, line=1, adj=0, cex=0.8)
add_background_stripes()
points(jitter(as.numeric(results$group), amount=0.2), results$execution_time,
       col=point_colors, pch=19)
axis(2, las=2)
axis(1, at=1:length(levels(results$group)), labels=levels(results$group), las=2)
add_means(exec_means)

# 2. Angular Error Plot
plot(1, type="n",
     xlim=c(0.5, length(levels(results$group))+0.5),
     ylim=range(results$angular_error),
     xlab="", ylab="Angular Distance (radians)",
     main=expression(bold(paste("Distance between mean angle and ", theta))),
     xaxt="n", yaxt="n")
add_background_stripes()
points(jitter(as.numeric(results$group), amount=0.2), results$angular_error,
       col=point_colors, pch=19)
axis(2, las=2)
axis(1, at=1:length(levels(results$group)), labels=levels(results$group), las=2)
add_means(circ_means)

# 3. Mean RT Plot
plot(1, type="n",
     xlim=c(0.5, length(levels(results$group))+0.5),
     ylim=range(results$mean_rt),
     xlab="", ylab="Time (seconds)",
     main="Mean Response Time",
     xaxt="n", yaxt="n")
add_background_stripes()
points(jitter(as.numeric(results$group), amount=0.2), results$mean_rt,
       col=point_colors, pch=19)
axis(2, las=2)
axis(1, at=1:length(levels(results$group)), labels=levels(results$group), las=2)
add_means(rt_means)

# 4. Legend Panel with both parameter sets and trial size legends
plot(1, type="n", xlab="", ylab="", main="", axes=FALSE)

# First legend for parameter sets (moved up)
legend("top", 
       legend=c(
           expression(paste("Very Easy (", mu[1], " = ", mu[2], " = 3.0)")),
           expression(paste("Easy (", mu[1], " = ", mu[2], " = 2.0)")),
           expression(paste("Medium (", mu[1], " = ", mu[2], " = 1.0)")),
           expression(paste("Hard (", mu[1], " = ", mu[2], " = 0.5)")),
           expression(paste("Very Hard (", mu[1], " = ", mu[2], " = 0.2)"))
       ),
       fill=point_colors,
       title=expression(bold("Parameter Sets")),
       cex=1.5,
       bty="n",
       inset=c(0, 0.075))

# Second legend for trial sizes
legend("bottom", 
       legend=paste(trial_sizes, "trials"),
       fill=c("lavender", "azure", "honeydew", "cornsilk"),
       title=expression(bold("Trial Sizes")),
       cex=1.5,
       bty="n",
       ncol=2,
       inset=c(0, 0.2))

# Second legend for trial sizes
legend("bottom", 
       legend=bquote(bold("Repetitions: ") * .(n_reps)),
       cex=1.5,
       bty="n",       
       inset=c(0, 0))

dev.off()

# Save results to RData file
save(results, file=sprintf("tests/results_%s_%s.RData", 
                          method_tested, format(Sys.Date(), "%Y%m%d")))

# Print summary statistics
summary_stats <- aggregate(
    cbind(execution_time, completion, angular_error, mean_rt) 
    ~ param_set + n_trials, 
    data = results,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
)
print(summary_stats)