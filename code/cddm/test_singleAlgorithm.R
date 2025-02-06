#########################################################################
# S I N G L E     A L G O R I T H M     T E S T ########################
#########################################################################
# This source code contains a custom function to:
# 1) Run a single sampling algorithm on a fixed set of parameters.
# 2) Calculate summary statistics and performance metrics.
# 3) Return results as a list.
#########################################################################
single_algorithm_test <- function(params, n_trials, method_tested) {
    # Create a new list with n and the existing parameters
    full_params <- list(n = n_trials,
                        par = params$par)

    ##################################################################
    # Run the sampling function
    ##################################################################
    start_time <- Sys.time()
    if (method_tested == "RandomWalk") {
        result_raw <- do.call(sample.RW.cddm, full_params)
        result <- result_raw$bivariate.data
    } else if (method_tested == "Metropolis") {
        result <- do.call(sample.Metropolis.cddm, full_params)
    } else if (method_tested == "inverseCDF") {
        result <- do.call(rCDDM_inverse, full_params)
    } else if (method_tested == "Rejection") {
        result <- do.call(sample.Reject.cddm, full_params)
    } else {        stop(paste("Unknown method:", method_tested))       }
    end_time <- Sys.time()

    ##################################################################
    # Calculate summary statistics
    ##################################################################
    # Execution time
    execution_time <- as.numeric(difftime(end_time, start_time, units="secs"))
    # Completion rate (i.e. proportion of trials that did not return NA)
    completion <- mean(!is.na(result$RT))    
     # Difference between mean radian-choice and radian-theta
    angles <- circular::circular(result$Choice, type="angles", units="radians")
    mean_angle <- circular::mean.circular(angles)        
    theoretical_theta <- circular::circular(atan2(params$par$mu2, params$par$mu1))        
    angular_error <- as.numeric(mean_angle - theoretical_theta)    
    # Proportion of trials with negative RTs (excluding NAs)
    prop_negative_rt <- mean(result$RT < 0, na.rm=TRUE)

    ##################################################################
    # Return the results
    ##################################################################
    output <- list(
        execution_time = execution_time,
        completion = completion,
        mean_rt = mean(result$RT, na.rm=TRUE),        
        sd_rt = sd(result$RT, na.rm=TRUE),
        mean_angle = abs(as.numeric(mean_angle - theoretical_theta)),        
        sd_angle = sd.circular(angles),
        prop_negative_rt = prop_negative_rt,
        angular_error = angular_error,
        data = result)

    if(method_tested == "RandomWalk") {
        # Calculate distance from boundary for each trial
        final_coords <- getFinalState(result_raw$random.walk)
        distances_from_boundary <- sqrt(rowSums(final_coords^2))
        circumference_precision <- mean(abs(distances_from_boundary - params$par$boundary))
        output$circumference_precision <- circumference_precision
    }
    return(output)
}


#########################################################################
# A D D    E M P I R I C A L     A N D   T H E O R E T I C A L ##########
#     C U M U L A T I V E      D E N S I T Y     F U N C T I O N   ######
#########################################################################
# This source code contains a custom function to:
# 1) Add empirical and theoretical cumulative density functions 
#    to the arrays of bivariate data.
# 2) Return the updated arrays.
#########################################################################
add_cdfs_to_arrays <- function(data_arrays, param_sets) {    
    # Progress counter
    total_iterations <- length(data_arrays) * length(data_arrays[[1]]) * dim(data_arrays[[1]][[1]])[3]
    current_iteration <- 0
    
    # For each parameter set
    for(param_name in names(data_arrays)) {
        mu1 <- param_sets[[param_name]]$par$mu1
        mu2 <- param_sets[[param_name]]$par$mu2
        drift <- sqrt(mu1^2 + mu2^2)
        theta <- atan2(mu2, mu1)
        tzero <- param_sets[[param_name]]$par$tzero
        boundary <- param_sets[[param_name]]$par$boundary
        
        # Combine all bivariate data for this parameter set
        all_bivariate_data <- do.call(rbind, lapply(data_arrays[[param_name]], function(arr) {
            # Combine all repetitions for this trial size
            do.call(rbind, lapply(1:dim(arr)[3], function(rep) arr[,1:2,rep]))
        }))
        
        # Compute theoretical CDF for all data at once
        cat(sprintf("\nComputing theoretical CDFs for %s (δ = %.2f)....\n", param_name, drift))
        all_tcdf <- pCDDM(all_bivariate_data, drift, theta, tzero, boundary, method="monte_carlo")
        cat("...done!\n")
        
        # Current position in the all_tcdf vector
        current_pos <- 1
        
        # For each trial size
        for(n_trials_char in names(data_arrays[[param_name]])) {
            n_trials <- as.numeric(n_trials_char)  # Convert to numeric
            n_reps <- dim(data_arrays[[param_name]][[n_trials_char]])[3]
            current_array <- data_arrays[[param_name]][[n_trials_char]]
            
            # Create new array with extra columns
            new_array <- array(NA, 
                             dim = c(dim(current_array)[1], 4, dim(current_array)[3]),
                             dimnames = list(
                                 NULL,
                                 c("Choice", "RT", "eCDF", "tCDF"),
                                 dimnames(current_array)[[3]]
                             ))
            
            # Copy existing Choice and RT data (explicitly selecting first two columns)
            new_array[, 1:2, ] <- current_array[, 1:2, ]
            
            # For each repetition
            for(rep in 1:n_reps) {
                current_iteration <- current_iteration + 1
                cat(sprintf("\rComputing eCDFs: Progress %d/%d (%.1f%%) - trials=%d, rep=%d",
                    current_iteration, total_iterations,
                    100 * current_iteration/total_iterations,
                    n_trials, rep))
                
                # Extract bivariate data for this repetition
                bivariate_data <- current_array[, 1:2, rep]
                
                # Compute empirical CDF
                new_array[, 3, rep] <- myECDF(bivariate_data)
                
                # Extract corresponding theoretical CDF values
                new_array[, 4, rep] <- round(all_tcdf[current_pos:(current_pos + n_trials - 1)],4)
                current_pos <- current_pos + n_trials
            }
            
            # Store updated array
            data_arrays[[param_name]][[n_trials_char]] <- new_array
        }
        cat("\n")
    }
    
    return(data_arrays)
}


#########################################################################
# Helper functions #######################################################
#########################################################################
# Helper functions
    add_means <- function(means) {
        segments(1:length(means) - 0.25, means, 1:length(means) + 0.25, means,
                lwd=3, col="black")
    }
    
    add_background_stripes <- function(stripe_colors, trial_sizes, param_sets, results) {
        n_trials <- length(trial_sizes)
        for(i in 1:n_trials) {
            rect(xleft = seq(0.5 + (i-1)*length(param_sets), 
                           max(as.numeric(results$group)) + 0.5), 
                         ybottom = par("usr")[3],
                         xright = seq(length(param_sets) + 0.5 + (i-1)*length(param_sets), 
                                    max(as.numeric(results$group)) + 0.5, 
                                    by=n_trials*length(param_sets)),
                         ytop = par("usr")[4],
                         col = stripe_colors[i],
                         border = NA)
            
            text(x = (length(param_sets)/2) + 0.5 + (i-1)*length(param_sets),
                 y = par("usr")[3] + (par("usr")[4] - par("usr")[3])*0.8,
                 labels = paste(trial_sizes[i], "trials"),
                 col = "gray50",
                 cex = 0.8)
        }
    }


#########################################################################
# P L O T     A L G O R I T H M     P E R F O R M A N C E   #############
#########################################################################
# This source code contains a custom function to:
# 1) Plot the performance of a single sampling algorithm in terms of
#    a) Execution time, b) mean choice - theta, c) mean response time
# 2) Return the plots.
#########################################################################
plot_algorithm_performance <- function(results, param_sets, trial_sizes, n_reps, method_tested, filename = NA) {
    # At the start of the function
    if(nrow(results) == 0) stop("No results to plot")
    if(length(param_sets) == 0) stop("No parameter sets provided")
    if(length(trial_sizes) == 0) stop("No trial sizes provided")
    
    if(!is.na(filename)) {     pdf(filename, width=12, height=10)     }
    
    # Set up 2x2 plotting layout with adjusted margins
    par(mfrow=c(2,2), oma=c(2,2,2,0), mar=c(4,4,3,3), mgp=c(2,0.7,0),
        cex.main=1.5, cex.lab=1.2, cex.axis=1.1)    
    # Calculate means
    exec_means <- tapply(results$execution_time, results$group, mean)
    circ_means <- tapply(results$angular_error, results$group, mean)
    rt_means <- tapply(results$mean_rt, results$group, mean)    
    
    # Create color vector based on number of parameter sets and method tested    
    if (method_tested == "RandomWalk") {
        color_palette <- colorRampPalette(c("#0d2889", "#47c2f7"))  # Blue to light gray
    } else if (method_tested == "Metropolis") {
        color_palette <- colorRampPalette(c("#5b218b", "#ad67d0"))  # Dark green to light green
    } else if (method_tested == "inverseCDF") {
        color_palette <- colorRampPalette(c("#E04D4D", "#E5B4B4"))  # Red to light pink
    } else if (method_tested == "Rejection") {
        color_palette <- colorRampPalette(c("#0c6a2a", "#65d279"))  # Purple to light purple
    } else {        
        stop(paste("Unknown method:", method_tested))       
    }
    
    # Generate colors for each parameter set
    point_colors <- color_palette(length(param_sets))[as.numeric(factor(results$param_set))]
    point_colors_transparent <- adjustcolor(point_colors, alpha=0.3)
    
    # Generate stripe colors for trial sizes
    shade_change <- c((1:length(trial_sizes)*10))
    stripe_colors <- rgb((255-shade_change)/255,(255-shade_change)/255,(255-shade_change)/255,0.5)
    
    # Generate parameter set labels based on drift values
    param_labels <- names(param_sets)    
    # Calculate drift lengths for each parameter set
    drift_lengths <- sapply(param_sets, function(ps) {
        mu1 <- ps$par$mu1
        mu2 <- ps$par$mu2
        sqrt(mu1^2 + mu2^2)  # Magnitude of drift vector
    })
    
    # Create legend labels with expression for delta
    legend_labels <- sapply(seq_along(param_sets), function(i) {
        parse(text = sprintf('"%s" * " (" * delta * " = " * %.2f * ")"', 
                           names(param_sets)[i], 
                           drift_lengths[i]))
    })
    
    # Create plots
    plot_metrics <- list(
        list(data=results$execution_time, means=exec_means, 
             ylab="Time (seconds)", main="Execution Time Distribution",
             add_reference=FALSE),
        list(data=results$angular_error, means=circ_means,
             ylab="Angular Distance (radians)", 
             main=expression(bold(paste("Distance between mean angle and ", theta))),
             add_reference=TRUE),  # Added flag for reference line
        list(data=results$mean_rt, means=rt_means,
             ylab="Time (seconds)", main="Mean Response Time",
             add_reference=FALSE)
    )
    
    # Create the three metric plots
    for(plot_info in plot_metrics) {
        plot(1, type="n", axes=FALSE, xlim=c(1, length(levels(results$group))),
             ylim=range(plot_info$data), xlab="", ylab="", main=plot_info$main,
             xaxt="n", yaxt="n")
        add_background_stripes(stripe_colors, trial_sizes, param_sets, results)
        # Add reference line if specified
        if(plot_info$add_reference) {
            abline(h=0, lty=3, col="red", lwd=3)  # dotted horizontal line at y=0
        }
        mtext(plot_info$ylab, side=2, line=3)
        points(jitter(as.numeric(results$group), amount=0.2), plot_info$data,
               col=point_colors_transparent, pch=19)
        axis(2, las=2, cex.axis=0.9, line=-0)
        axis(1, at=1:length(levels(results$group)), 
             labels=rep(param_labels, length(trial_sizes)), 
             las=2, cex.axis=0.9)
        add_means(plot_info$means)        
    }
    
    # Create legend panel
    plot(1, type="n", xlab="", ylab="", main="", axes=FALSE)    
    # Parameter sets legend
    legend("top", legend=legend_labels,
           col=adjustcolor(unique(point_colors), alpha=0.5),
           pch=19, pt.cex=2, title=expression(bold("Parameter Sets")),
           cex=1.5, bty="n", inset=c(0, 0.075))    
    # Trial sizes legend
    legend("bottom", legend=paste(trial_sizes, "trials"),
           fill=stripe_colors, title=expression(bold("Trial Sizes")),
           cex=1.5, bty="n", ncol=min(2, length(trial_sizes)),
           inset=c(0, 0.6/length(param_sets)))    
    
    # Create simulation info with only labels in bold
    sim_info <- c(
        substitute(bold("Algorithm:") ~ x, list(x = method_tested)),
        substitute(bold("Repetitions:") ~ x, list(x = n_reps)),
        substitute(bold("Date:") ~ x, list(x = format(Sys.Date(), "%Y-%m-%d")))
    )
    
    # Add simulation information
    legend("bottom", legend = parse(text = sapply(sim_info, deparse)), 
           cex=1.5, bty="n", inset=c(0, 0))

    if(!is.na(filename)) {    dev.off()     }
}

######################################################################
# Function specific to the Random Walk algorithm #####################
######################################################################
# Show the distance between the walk end and the boundary 
######################################################################
plot_circumference_precision <- function(results, param_sets, trial_sizes, method_tested, filename = NA) {    
    if(!is.na(filename)){    
        pdf(filename, width=10, height=8)
    }
    
    # Calculate statistics
    y_range <- range(results$circumference_precision)
    y_mid <- mean(y_range)
    y_values <- c(y_range[1], y_mid, y_range[2])
    circ_means <- tapply(results$circumference_precision, results$group, mean)
    
    # Create color vector based on parameter sets
    point_colors <- rainbow(length(param_sets))[as.numeric(factor(results$param_set))]
    point_colors_transparent <- adjustcolor(point_colors, alpha=0.3)
    shade_change <- c((1:length(trial_sizes)*10))
    stripe_colors <- rgb((255-shade_change)/255,(255-shade_change)/255,(255-shade_change)/255,0.5)
    
    # Set up plotting parameters
    par(mfrow=c(1,1), 
        mar=c(5,7,3,2),    # Increased left margin from 5 to 7
        cex.main=1.8, 
        cex.lab=1.2, 
        cex.axis=1.1)
    
    # Create plot
    plot(jitter(as.numeric(results$group), amount=0.2), 
         results$circumference_precision,
         col=point_colors_transparent,  pch=19, 
         main="Circumference Precision Distribution",
         xlab="Parameter Set", ylab="", xaxt="n", yaxt="n", 
         xlim=c(0.5, length(levels(results$group)) + 0.5))
    
    # Add y-axis with scientific notation
    axis(2, at=y_values, labels=sprintf("%.2e", y_values), las=2)
    add_background_stripes(stripe_colors, trial_sizes, param_sets, results_cdfs) 
    mtext("Average Distance from Boundary", side=2, line=5.5)

    # Add mean lines
    segments(1:length(circ_means) - 0.25, circ_means,
             1:length(circ_means) + 0.25, circ_means,
             lwd=2, col="black")
    
    # Add mean values as text
    text(1:length(circ_means) + 0.3, 
         circ_means, 
         sprintf("%.4f", circ_means),
         adj=0, 
         cex=0.8)
    
    # Add x-axis labels
    axis(1, 
         at=seq_along(levels(results$group)), 
         labels=levels(results$group))
    
    # Calculate drift lengths for each parameter set
    drift_lengths <- sapply(param_sets, function(ps) {
        mu1 <- ps$par$mu1
        mu2 <- ps$par$mu2
        sqrt(mu1^2 + mu2^2)  # Magnitude of drift vector
    })
    
    # Create legend labels with expression for delta
    legend_labels <- sapply(seq_along(param_sets), function(i) {
        parse(text = sprintf('"%s" * " (" * delta * " = " * %.2f * ")"', 
                           names(param_sets)[i], 
                           drift_lengths[i]))
    })
    
    # Add legend with delta values
    legend("topright",
           legend=legend_labels,
           col=adjustcolor(unique(point_colors), alpha=0.5),
           pch=19,
           title="Parameter Sets",
           cex=1.2,
           bty="n")
    
    if(!is.na(filename)){
                dev.off()
    }
}

######################################################################
# Plot the empirical and theoretical CDFs ############################
######################################################################
plot_cdfs <- function(data_arrays) {
    # Check if scatterplot3d is installed
    if (!require("scatterplot3d")) {
        install.packages("scatterplot3d")
        library(scatterplot3d)
    }
    
    # Define pastel colors for each parameter set
    param_colors <- c("fast" = "#A7C7E7",  # Pastel blue
                     "slow" = "#FFB347")   # Pastel orange
    
    # Get dimensions for plot grid
    n_param_sets <- length(data_arrays)
    n_trial_sizes <- length(data_arrays[[1]])
    
    # Set up plotting grid
    par(mfrow=c(n_param_sets, n_trial_sizes), mar=c(1,1,2,0.5), 
                oma=c(2,2,2,1), mgp=c(1,0.5,0), tcl=-0.2, pty="s")
    
    # For each parameter set (rows)
    for(param_name in names(data_arrays)) {
        current_color <- param_colors[param_name]
        
        # For each trial size (columns)
        for(n_trials_char in names(data_arrays[[param_name]])) {
            n_trials <- as.numeric(n_trials_char)
            array_data <- data_arrays[[param_name]][[n_trials_char]]
            n_reps <- dim(array_data)[3]
            
            # Initialize empty 3D plot
            s3d <- scatterplot3d(x = array_data[, "Choice", 1],
                               y = array_data[, "RT", 1],
                               z = array_data[, "tCDF", 1],
                               angle = 30, type = "n",
                               xlab = "", ylab = "", zlab = "",
                               main = sprintf("%s\n%d trials", param_name, n_trials),
                               cex.main = 0.8, grid = TRUE, box = TRUE,
                               mar = c(1,1,2,0.5))   # Consistent with par margins
            
            # First, plot all empirical CDF points
            for(rep in 1:n_reps) {
                s3d$points3d(x = array_data[, "Choice", rep],
                            y = array_data[, "RT", rep],
                            z = array_data[, "eCDF", rep],
                            col = adjustcolor(current_color, alpha=0.8),
                            pch = 19, cex = 0.4)
            }
            
            # Then overlay theoretical CDF points
            for(rep in 1:n_reps) {
                s3d$points3d(x = array_data[, "Choice", rep],
                            y = array_data[, "RT", rep],
                            z = array_data[, "tCDF", rep],
                            col = "black", pch = 19, cex = 0.2)
            }
            
            # Add legend only to first plot
            if(param_name == names(data_arrays)[1] && 
               n_trials_char == names(data_arrays[[1]])[1]) {
                legend("topleft",
                       legend=c("Empirical CDFs", "Theoretical CDF"),
                       col=c(adjustcolor(current_color, alpha=0.1), "black"),
                       pch=c(19, 19), pt.cex=c(0.5, 0.3), bty="n", cex=0.6, inset=0.02)
            }
        }
    }
    
    # Add shared labels (closer to plots)
    mtext("Choice Angle (radians)", side=1, outer=TRUE, line=0.3, cex=0.8)
    mtext("Response Time", side=4, outer=TRUE, line=-0.1, cex=0.8)
    mtext("CDF", side=2, outer=TRUE, line=0.3, cex=0.8)
    mtext("Empirical vs Theoretical CDF", side=3, outer=TRUE, line=0.3, cex=0.9)
}


######################################################################
# Plot the differences between the empirical and theoretical CDFs ####
######################################################################
calculate_cdf_metrics <- function(data_arrays) {
    # Initialize empty results dataframe
    results <- data.frame()
    
    # Calculate metrics for each parameter set and trial size
    for(param_name in names(data_arrays)) {
        for(n_trials_char in names(data_arrays[[param_name]])) {
            array_data <- data_arrays[[param_name]][[n_trials_char]]
            n_trials <- as.numeric(n_trials_char)
            n_reps <- dim(array_data)[3]            
            
            for(rep in 1:n_reps) {
                # Calculate differences between empirical and theoretical CDFs
                differences <- array_data[, "eCDF", rep] - array_data[, "tCDF", rep]
                
                # Calculate summary statistics
                output <- data.frame(
                    mean_diff = mean(differences),
                    median_diff = median(differences),
                    ssd = sum(differences^2),  # Sum of squared differences
                    ks_stat = max(abs(differences)),  # Kolmogorov-Smirnov statistic                    
                    param_set = param_name,
                    trial_size = n_trials,
                    rep = rep
                )
                results <- rbind(results, output)
            }
        }
    }
    
    return(results)
}

plot_cdf_differences <- function(results_cdfs, data_arrays, param_sets, trial_sizes, n_reps, method_tested, filename = NA) {    
    if(!is.na(filename)) {    
        pdf(filename, width=12, height=6)  # Reduced height since we're removing one row
    }
    
    # Set up 2x2 plotting grid
    par(mfrow=c(2,2), oma=c(2,2,2,0), mar=c(4,4,3,3), mgp=c(2,0.7,0),
        cex.main=1.5, cex.lab=1.2, cex.axis=1.1)
    
    # Calculate drift information for parameter sets
    drift_info <- lapply(param_sets, function(ps) {
        magnitude <- sqrt(ps$par$mu1^2 + ps$par$mu2^2)
        angle_rad <- atan2(ps$par$mu2, ps$par$mu1)
        list(magnitude = magnitude, angle = angle_rad)
    })
    
    # Setup colors
    point_colors <- rainbow(length(param_sets))
    point_colors_transparent <- adjustcolor(point_colors, alpha=0.3)    
    indivPoint_colors_transparent <- adjustcolor(point_colors, alpha=0.1)
    shade_change <- c((1:length(trial_sizes)*10))
    stripe_colors <- rgb((255-shade_change)/255,(255-shade_change)/255,(255-shade_change)/255,0.5)

    # Collect all differences and group information
    all_diffs <- c()
    all_param_sets <- c()
    group_labels <- c()
    
    for(param_name in names(data_arrays)) {
        for(n_trials_char in names(data_arrays[[param_name]])) {
            array_data <- data_arrays[[param_name]][[n_trials_char]]
            n_trials <- as.numeric(n_trials_char)
            
            for(rep in 1:dim(array_data)[3]) {
                diffs <- array_data[, "eCDF", rep] - array_data[, "tCDF", rep]
                all_diffs <- c(all_diffs, diffs)
                all_param_sets <- c(all_param_sets, rep(param_name, length(diffs)))
                group_labels <- c(group_labels, rep(paste(n_trials, param_name), length(diffs)))
            }
        }
    }
    
    # Use the same group factor as results_cdfs for consistent ordering
    all_groups <- factor(group_labels, levels = levels(results_cdfs$group))
    
    # First plot: All differences
    means <- tapply(all_diffs, all_groups, mean)
    
    plot(1, type="n", axes=FALSE, 
         xlim=c(1, length(levels(all_groups))),
         ylim=range(all_diffs), 
         xlab="", ylab="", 
         main="CDF Differences Distribution",
         xaxt="n", yaxt="n")
    add_background_stripes(stripe_colors, trial_sizes, param_sets, results_cdfs)    
    # Add points and reference line
    mtext("CDF Difference (eCDF - tCDF)", side=2, line=3)
    points(jitter(as.numeric(all_groups), amount=0.2), all_diffs,
           col=indivPoint_colors_transparent[as.numeric(factor(all_param_sets))], 
           pch=19, cex=0.5)  # Smaller points due to larger number
    axis(2, las=2, cex.axis=0.9)
    axis(1, at=1:length(levels(all_groups)), 
         labels=levels(all_groups),
         las=2, cex.axis=0.9)
    
    # Add means
    segments(1:length(means) - 0.25, means,
            1:length(means) + 0.25, means,
            lwd=2, col="black")
    
    # Add reference line at zero
    abline(h=0, lty=2, col="gray")
    
    # Create remaining plots for summary metrics
    plot_metrics <- list(
        list(data = results_cdfs$ssd, 
             ylab = "Sum of Squared Differences", 
             main = "CDF Squared Error"),
        list(data = results_cdfs$ks_stat,
             ylab = "Maximum Absolute Difference", 
             main = "Kolmogorov-Smirnov Statistic")
    )
    
    # Create the remaining metric plots
    for(plot_info in plot_metrics) {
        means <- tapply(plot_info$data, results_cdfs$group, mean)
        
        plot(1, type="n", axes=FALSE, 
             xlim=c(1, length(levels(results_cdfs$group))),
             ylim=range(plot_info$data), 
             xlab="", ylab="", main=plot_info$main,
             xaxt="n", yaxt="n")        
        add_background_stripes(stripe_colors, trial_sizes, param_sets, results_cdfs)        
        # Add points and axes
        mtext(plot_info$ylab, side=2, line=3)
        points(jitter(as.numeric(results_cdfs$group), amount=0.2), 
               plot_info$data,
               col=point_colors_transparent[as.numeric(factor(results_cdfs$param_set))], 
               pch=19)
        axis(2, las=2, cex.axis=0.9)
        axis(1, at=1:length(levels(results_cdfs$group)), 
             labels=levels(results_cdfs$group),
             las=2, cex.axis=0.9)
        
        # Add means
        segments(1:length(means) - 0.25, means,
                1:length(means) + 0.25, means,
                lwd=2, col="black")
    }
    
    drift_lengths <- sapply(param_sets, function(ps) {
        mu1 <- ps$par$mu1
        mu2 <- ps$par$mu2
        sqrt(mu1^2 + mu2^2)  # Magnitude of drift vector
    })
    
    # Create legend labels with expression for delta
    legend_labels <- sapply(seq_along(param_sets), function(i) {
        parse(text = sprintf('"%s" * " (" * delta * " = " * %.2f * ")"', 
                           names(param_sets)[i], 
                           drift_lengths[i]))
    })

     # Create legend panel
    plot(1, type="n", xlab="", ylab="", main="", axes=FALSE)    
    # Parameter sets legend
    legend("top", legend=legend_labels,
           col=adjustcolor(unique(point_colors), alpha=0.5),
           pch=19, pt.cex=2, title=expression(bold("Parameter Sets")),
           cex=1.5, bty="n", inset=c(0, 0.075))    
    # Trial sizes legend
    legend("bottom", legend=paste(trial_sizes, "trials"),
           fill=stripe_colors, title=expression(bold("Trial Sizes")),
           cex=1.5, bty="n", ncol=min(2, length(trial_sizes)),
           inset=c(0, 0.6/length(param_sets)))    
    # Repetitions info
    legend("bottom",  legend=bquote(bold("Repetitions: ") * .(n_reps)),
           cex=1.5, bty="n", inset=c(0, 0.15))

    if(!is.na(filename)) {    
        dev.off()     
    }
}

# Usage example:
#figure_cdf_differences <- sprintf(here("results", "quickTest_%s_%s_cdfs_differences.pdf"), method_tested, format(Sys.Date(), "%Y%m%d"))
#plot_cdf_differences(results_cdfs, data_arrays, param_sets, trial_sizes, n_reps, filename = figure_cdf_differences)

plot_metrics_by_paramset <- function(results_cdfs, results, param_sets, trial_sizes, n_reps, method_tested, filename_prefix) {
    # Define metrics to plot
    metrics <- list(
        list(data_cdfs = results_cdfs$ks_stat, name = "ks",
             title = "Kolmogorov-Smirnov Statistic Distribution",
             ylab = "KS Statistic"),
        list(data_exec = results$execution_time, name = "exec",
             title = "Execution Time Distribution",
             ylab = "Time (seconds)")
    )
    
    # Start PDF device
    pdf(filename_prefix, width=12, height=7)
    
    # Calculate layout dimensions
    n_metrics <- length(metrics)
    n_param_sets <- length(param_sets)    
    layout_matrix <- matrix(1:(n_metrics * n_param_sets), 
                          nrow=n_metrics, 
                          ncol=n_param_sets, 
                          byrow=TRUE)
    layout(layout_matrix)    
    par(mar=c(1,3,1,1), oma=c(4,3,4,1), mgp=c(2.5,1,0))

    # Set colors based on method
    if (method_tested == "RandomWalk") {
        color_indiv <- adjustcolor("#4D9DE0", alpha=0.3); color_mean <- "#008cff"
    } else if (method_tested == "Metropolis") {
        color_indiv <- adjustcolor("#4D9DE0", alpha=0.3); color_mean <- "#008cff"
    } else if (method_tested == "inverseCDF") {
        color_indiv <- adjustcolor("#4D9DE0", alpha=0.3); color_mean <- "#008cff"
    } else if (method_tested == "Rejection") {
        color_indiv <- adjustcolor("#51E04D", alpha=0.3); color_mean <- "#00FF51"
    } else {        
        stop(paste("Unknown method:", method_tested))       
    }
    
    # Pre-calculate y-axis ranges for each metric
    y_ranges <- list()
    for(metric in metrics) {
        if(metric$name == "ks") {
            y_ranges[[metric$name]] <- range(metric$data_cdfs)
        } else {
            y_ranges[[metric$name]] <- range(metric$data_exec)
        }
    }
    
    # Create plots
    for(metric_idx in seq_along(metrics)) {
        metric <- metrics[[metric_idx]]
        for(param_idx in seq_along(names(param_sets))) {
            param_name <- names(param_sets)[param_idx]            
            
            # Get data for this parameter set
            if(metric$name == "ks") {
                param_data <- metric$data_cdfs[results_cdfs$param_set == param_name]
                trial_sizes_data <- results_cdfs$trial_size[results_cdfs$param_set == param_name]
            } else {
                param_data <- metric$data_exec[results$param_set == param_name]
                trial_sizes_data <- results$n_trials[results$param_set == param_name]
            }            
            
            # Calculate drift for title
            drift <- sqrt(param_sets[[param_name]]$par$mu1^2 + 
                        param_sets[[param_name]]$par$mu2^2)            
            
            # Create empty plot with consistent y-axis range for each metric
            plot(1, type="n", xlim=range(trial_sizes) + c(-50, 50), yaxt="n",
                 ylim=y_ranges[[metric$name]], xlab="", ylab="", main="", xaxt="n")  

            # Add y-axis with exactly 7 ticks
            y_ticks <- seq(y_ranges[[metric$name]][1], 
                          y_ranges[[metric$name]][2], 
                          length.out=7)
            axis(2, at=y_ticks, labels=round(y_ticks, 2), las=2, cex.axis=1.1)
            
            # Add x-axis with values only for bottom row plots
            if(metric_idx == length(metrics)) {
                axis(1, at=trial_sizes, labels=trial_sizes, cex.axis=1.5)
                mtext("Number of Trials", side=1, line=3.5)
            } else {
                axis(1, at=trial_sizes, labels=FALSE)
            }

            if(param_idx == 1){
                mtext(metric$ylab, side=2, line=3.5, cex=1, font=1)
            } 
            
            # Only show parameter set title on top row
            if(metric_idx == 1) {                
                mtext(bquote(.(sub("δ", "", param_name)) * " " * (delta * " = " * .(sprintf("%.2f", drift)))),
                      line=0.2, cex=1)
            } 

            # Store and plot mean values
            mean_values <- numeric(length(unique(trial_sizes_data)))
            trial_sizes_unique <- sort(unique(trial_sizes_data))            
            
            for(trial_idx in seq_along(trial_sizes_unique)) {
                n_trials <- trial_sizes_unique[trial_idx]
                data_points <- param_data[trial_sizes_data == n_trials]
                x_jittered <- jitter(rep(n_trials, length(data_points)), amount=20)
                points(x_jittered, data_points, col=color_indiv, pch=19)                
                abline(v=n_trials, lty=3, col="gray50")                
                mean_values[trial_idx] <- mean(data_points)                
                segments(n_trials-80, mean_values[trial_idx],
                        n_trials+80, mean_values[trial_idx],
                        col="black", lwd=3)
            }
            
            # Add connected mean points            
            lines(trial_sizes_unique, mean_values, type="b",
                  lwd=3, pch=21, cex=2.5,
                  bg=color_mean, col="black",
                  lty="dotted")            
            
            # Add grid
            grid(nx=NA, ny=NULL, col="gray", lty="dotted")
        }
    }
    
    # Add overall title
    mtext(sprintf("Performance Metrics by Parameter Set and Trial Size\n%s Algorithm (%s)", 
          method_tested, format(Sys.Date(), "%Y-%m-%d")), 
          outer=TRUE, line=1.2, cex=1.2, font=2)    
    
    dev.off()
}


# Generate the plots
#filename_prefix <- here("results", 
#                       sprintf("quickTest_%s_%s_metrics", 
#                              method_tested, 
#                              format(Sys.Date(), "%Y%m%d")))
#plot_metrics_by_paramset(results_cdfs, results, param_sets, trial_sizes, n_reps, 
#                        filename_prefix)