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
param_sets <- list( easy = list(par = list(
                                 mu1 = 2.0,  mu2 = 2.0,  # (strong drift)
                                 boundary = 5, tzero = 0.1
                               )),
                    hard = list(par = list(
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

# Fix the order in which cells will be displayed
plot_order <- paste(rep(trial_sizes, each=2), 
                    rep(c("easy", "hard"), length(trial_sizes)))           

# Prepare results
results$group <- factor(paste(results$n_trials, results$param_set),
                       levels = plot_order)
# Save results! 
filename <- sprintf(here("results", "quickTest_%s_%s_execTime.RData"), method_tested, format(Sys.Date(), "%Y%m%d"))
save(results,file = filename)
cat("Saving results to:", filename, "\n\n")
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
figname_results <- sprintf(here("results", "quickTest_%s_%s.pdf"), method_tested, format(Sys.Date(), "%Y%m%d"))

# Create color vector based on parameter set
point_colors <- rainbow(length(param_sets))[as.numeric(factor(results$param_set))]
# Add transparency to point colors
point_colors_transparent <- adjustcolor(point_colors, alpha=0.3)  # 0.3 for 70% transparency
# Define plotting space variables   
x_positions <- as.numeric(results$group)
# Calculate means for each group
exec_means <- tapply(results$execution_time, results$group, mean)


pdf(figname_results, width=12, height=10)

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


cat("\nFigures saved to:", figname_execTime, "and", figname_circPrecision, "\n\n")

# Print summary
cat("\nSummary of results:\n")
print(summary_stats)

cat("\n\nDone!\n\n")