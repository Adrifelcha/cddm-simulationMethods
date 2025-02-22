#############################################################
method_tested <- "RandomWalk"
# Possible methods:
# 1) "Metropolis"
# 2) "RandomWalk"
# 3) "inverseCDF"
# 4) "Rejection"

forceRun <- TRUE # Ignore existing results and run again
#############################################################
library(circular)
library(foreach)
library(doParallel)
library(here)
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
    very_fast = list(par = list(
        mu1 = 3.0, mu2 = 3.0,     # very fast drift (δ = 4.24)
        boundary = 5, tzero = 0.1
    )),
    fast = list(par = list(
        mu1 = 2.0, mu2 = 2.0,     # fast drift (δ = 2.83)
        boundary = 5, tzero = 0.1
    )),
    medium = list(par = list(
        mu1 = 1.0, mu2 = 1.0,     # moderate drift (δ = 1.41)
        boundary = 5, tzero = 0.1
    )),
    slow = list(par = list(
        mu1 = 0.5, mu2 = 0.5,     # slow drift (δ = 0.71)
        boundary = 5, tzero = 0.1
    )),
    very_slow = list(par = list(
        mu1 = 0.2, mu2 = 0.2,     # very slow drift (δ = 0.28)
        boundary = 5, tzero = 0.1
    ))
)

# More comprehensive trial sizes
trial_sizes <- c(50, 150, 300, 500)

# Number of replications
n_reps <- 10

# Define plot and results filenames
plotname <- sprintf("tests/testParallel_%s_%s.pdf", method_tested, format(Sys.Date(), "%Y%m%d"))
filename <- sprintf("tests/resultsParallel_%s_%s_n%s.RData", method_tested, format(Sys.Date(), "%Y%m%d"),n_reps)

#############################################################
#### R U N N I N G     T E S T S ############################
#############################################################

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Determine if we should run the tests
if(!forceRun) {
    if(file.exists(filename)) {
        results <- readRDS(filename)
        cat("Results already exist, skipping tests")
    }
} else {

# Initialize results dataframe with additional columns
results <- data.frame()

# Progress counter
total_iterations <- length(names(param_sets)) * length(trial_sizes) * n_reps
current_iteration <- 0

# Setup parallel processing
cores <- detectCores()
my.cluster <- makeCluster(cores[1]-2)  # Leave 2 cores free
registerDoParallel(cl = my.cluster)

# Run parallel tests
results <- foreach(seed = 1:n_reps,
                  .errorhandling = "pass",
                  .combine = 'rbind',
                  .packages = c('circular')) %dopar% {    
    # Set seed for reproducibility within each parallel process
    set.seed(seed)
   
    current_iteration <- current_iteration + 1   
    # Progress indicator
    cat(sprintf("\rProgress: %d/%d (%.1f%%) - Running seed=%d",
        current_iteration, total_iterations,
        100 * current_iteration/total_iterations, seed))

    # Container for this seed's results
    seed_results <- data.frame()    
    # Run tests for each parameter set and trial size
    for(param_name in names(param_sets)) {
        for(n_trials in trial_sizes) {
            # Run the test
            test_result <- run_full_test(
                params = param_sets[[param_name]],
                n_trials = n_trials,
                method_tested = method_tested
            )
            
            # Store results with metadata
            result_row <- data.frame(
                seed = seed,
                param_set = param_name,
                n_trials = n_trials,
                execution_time = test_result$execution_time,
                angular_error = test_result$angular_error,
                mean_rt = test_result$mean_rt,
                mean_angle = test_result$mean_angle
            )
            
            seed_results <- rbind(seed_results, result_row)
        }
    }
    seed_results    
}

stopCluster(cl = my.cluster)


# Create ordered groups
plot_order <- paste(rep(trial_sizes, each=length(param_sets)), 
                   rep(names(param_sets), length(trial_sizes)))
results$group <- factor(paste(results$n_trials, results$param_set),
                       levels = plot_order)
}
#############################################################
#### G E T     R E S U L T S ################################
#############################################################
# Save results
saveRDS(results, file = filename)

# Basic analysis of results
summary_stats <- aggregate(
    cbind(execution_time, angular_error, mean_rt) ~ param_set + n_trials,
    data = results,
    FUN = function(x) c(mean = mean(x), sd = sd(x))
)


#############################################################
#### P L O T T I N G     R E S U L T S ######################
#############################################################
pdf(plotname, width=12, height=10)

# Set up 2x2 plotting layout
par(mfrow=c(2,2), oma=c(2,2,2,2), mar=c(4,4,3,1), mgp=c(2,0.7,0),
    cex.main=1.5, cex.lab=1.2, cex.axis=1.1)

# Calculate means for plotting
exec_means <- tapply(results$execution_time, results$group, mean)
circ_means <- tapply(results$angular_error, results$group, mean)
rt_means <- tapply(results$mean_rt, results$group, mean)
compl_means <- tapply(results$completion, results$group, mean)
cex.axis = 0.9
xlim <- c(1, length(levels(results$group)))

# Create color vector based on parameter set
point_colors <- rainbow(length(param_sets))[as.numeric(factor(results$param_set))]
# Add transparency to point colors
point_colors_transparent <- adjustcolor(point_colors, alpha=0.3)  # 0.3 for 70% transparency

# Function to add mean lines and labels
add_means <- function(means) {
    segments(1:length(means) - 0.25, means,
             1:length(means) + 0.25, means,
             lwd=3, col="black")
}


# Define four light colors for the stripes
stripe_colors <- c("lavender", "azure", "honeydew", "cornsilk")

# Function to add background stripes with different colors
add_background_stripes <- function(stripe_colors) {    
    # Add each colored stripe
    for(i in 1:4) {
        rect(xleft = seq(0.5 + (i-1)*5, max(as.numeric(results$group)) + 0.5, by=20),
             ybottom = par("usr")[3],
             xright = seq(5.5 + (i-1)*5, max(as.numeric(results$group)) + 0.5, by=20),
             ytop = par("usr")[4],
             col = stripe_colors[i],
             border = NA)
        
        # Add trial size text
        text(x = 2.5 + (i-1)*5,  # Center of each stripe
             y = par("usr")[3] + (par("usr")[4] - par("usr")[3])*0.6,  # Middle of plot
             labels = paste(trial_sizes[i], "trials"),
             col = "gray50",  # Subtle gray color
             cex = 0.8,       # Smaller text
             srt = 0)        # Rotate text 90 degrees
    }
}

# Define speed labels
speed_labels <- c("Very Fast", "Fast", "Medium", "Slow", "Very Slow")

# 1. Execution Time Plot
plot(1, type="n", axes=FALSE,
     xlim=xlim,
     ylim=range(results$execution_time),
     xlab="", ylab="Time (seconds)",
     main="Execution Time Distribution",
     xaxt="n", yaxt="n")
mtext(paste0(method_tested, "\n", format(Sys.Date(), "%Y-%m-%d")), 
      side=3, line=0, adj=0, cex=0.8, outer=TRUE)
add_background_stripes(stripe_colors)
points(jitter(as.numeric(results$group), amount=0.2), results$execution_time,
       col=point_colors_transparent, pch=19)
axis(2, las=2, cex.axis=cex.axis, line=-0.3)
axis(1, at=1:length(levels(results$group)), labels=rep(speed_labels, length.out=length(levels(results$group))), las=2, cex.axis=cex.axis)
add_means(exec_means)

# 2. Angular Error Plot
plot(1, type="n", axes=FALSE,
     xlim=xlim,
     ylim=range(results$angular_error),
     xlab="", ylab="Angular Distance (radians)",
     main=expression(bold(paste("Distance between mean angle and ", theta))),
     xaxt="n", yaxt="n")
add_background_stripes(stripe_colors)
points(jitter(as.numeric(results$group), amount=0.2), results$angular_error,
       col=point_colors_transparent, pch=19)
axis(2, las=2, cex.axis=cex.axis, line=-0.3)
axis(1, at=1:length(levels(results$group)), labels=rep(speed_labels, length.out=length(levels(results$group))), las=2, cex.axis=cex.axis)
add_means(circ_means)

# 3. Mean RT Plot
plot(1, type="n", axes=FALSE,
     xlim=xlim,
     ylim=range(results$mean_rt),
     xlab="", ylab="Time (seconds)",
     main="Mean Response Time",
     xaxt="n", yaxt="n")
add_background_stripes(stripe_colors)
points(jitter(as.numeric(results$group), amount=0.2), results$mean_rt,
       col=point_colors_transparent, pch=19)
axis(2, las=2, cex.axis=cex.axis, line=-0.3)
axis(1, at=1:length(levels(results$group)), labels=rep(speed_labels, length.out=length(levels(results$group))), las=2, cex.axis=cex.axis)
add_means(rt_means)

# 4. Legend Panel with both parameter sets and trial size legends
plot(1, type="n", xlab="", ylab="", main="", axes=FALSE)

# First legend for parameter sets (moved up)
legend("top", 
       legend=c(
           expression(paste("Very Fast (", delta, " = 4.24)")),
           expression(paste("Fast (", delta, " = 2.83)")),
           expression(paste("Medium (", delta, " = 1.41)")),
           expression(paste("Slow (", delta, " = 0.71)")),
           expression(paste("Very Slow (", delta, " = 0.28)"))
       ),
       col=adjustcolor(unique(point_colors), alpha=0.5),
       pch=19,
       pt.cex=2,
       title=expression(bold("Parameter Sets")),
       cex=1.5,
       bty="n",
       inset=c(0, 0.075))

# Second legend for trial sizes
legend("bottom", 
       legend=paste(trial_sizes, "trials"),
       fill=stripe_colors,
       title=expression(bold("Trial Sizes")),
       cex=1.5,
       bty="n",
       ncol=2,
       inset=c(0, 0.15))

# Second legend for trial sizes
legend("bottom", 
       legend=bquote(bold("Repetitions: ") * .(n_reps)),
       cex=1.5,
       bty="n",       
       inset=c(0, 0.015))

dev.off()

# Save results to RData file
save(results, file=sprintf("tests/results_%s_%s.RData", 
                          method_tested, format(Sys.Date(), "%Y%m%d")))

# Print summary statistics
print(summary_stats)