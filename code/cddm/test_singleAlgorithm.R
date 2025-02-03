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
        result <- do.call(sample.invCDF.cddm, full_params)
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
    # Create color vector based on number of parameter sets
    point_colors <- rainbow(length(param_sets))[as.numeric(factor(results$param_set))]
    point_colors_transparent <- adjustcolor(point_colors, alpha=0.3)    
    #stripe_colors <- colorRampPalette(c("lavender", "azure", "honeydew", "cornsilk"))(length(trial_sizes))
    shade_change <- c((1:length(trial_sizes)*10))
    stripe_colors <- rgb((255-shade_change)/255,(255-shade_change)/255,(255-shade_change)/255,0.5)
    
    # Helper functions
    add_means <- function(means) {
        segments(1:length(means) - 0.25, means, 1:length(means) + 0.25, means,
                lwd=3, col="black")
    }
    
    add_background_stripes <- function(stripe_colors) {
        n_trials <- length(trial_sizes)
        for(i in 1:n_trials) {
            rect(xleft = seq(0.5 + (i-1)*length(param_sets), 
                           max(as.numeric(results$group)) + 0.5, 
                           by=n_trials*length(param_sets)),
                 ybottom = par("usr")[3],
                 xright = seq(length(param_sets) + 0.5 + (i-1)*length(param_sets), 
                            max(as.numeric(results$group)) + 0.5, 
                            by=n_trials*length(param_sets)),
                 ytop = par("usr")[4],
                 col = stripe_colors[i],
                 border = NA)
            
            text(x = (length(param_sets)/2) + 0.5 + (i-1)*length(param_sets),
                 y = par("usr")[3] + (par("usr")[4] - par("usr")[3])*0.6,
                 labels = paste(trial_sizes[i], "trials"),
                 col = "gray50",
                 cex = 0.8)
        }
    }
    
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
             ylab="Time (seconds)", main="Execution Time Distribution"),
        list(data=results$angular_error, means=circ_means,
             ylab="Angular Distance (radians)", 
             main=expression(bold(paste("Distance between mean angle and ", theta)))),
        list(data=results$mean_rt, means=rt_means,
             ylab="Time (seconds)", main="Mean Response Time")
    )
    
    # Create the three metric plots
    for(plot_info in plot_metrics) {
        plot(1, type="n", axes=FALSE, xlim=c(1, length(levels(results$group))),
             ylim=range(plot_info$data), xlab="", ylab="", main=plot_info$main,
             xaxt="n", yaxt="n")
        add_background_stripes(stripe_colors)
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
    # Repetitions info
    legend("bottom",  legend=bquote(bold("Repetitions: ") * .(n_reps)),
           cex=1.5, bty="n", inset=c(0, 0.15))

    if(!is.na(filename)) {    dev.off()     }
}

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


add_cdfs_to_arrays <- function(data_arrays, param_sets) {    
    # Progress counter
    # Total iterations = number of parameter sets * number of trial sizes * number of repetitions
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
        all_tcdf <- pCDDM(all_bivariate_data, drift, theta, tzero, boundary, method="monte_carlo", n_points=100000)
        cat("...done!\n")
        
        # Current position in the all_tcdf vector
        current_pos <- 1
        
        # For each trial size
        for(n_trials_char in names(data_arrays[[param_name]])) {
            n_trials <- as.numeric(n_trials_char)  # Convert to numeric
            n_reps <- dim(data_arrays[[param_name]][[n_trials_char]])[3]
            current_array <- data_arrays[[param_name]][[n_trials_char]]
            
            # Create new array with extra columns for eCDF and tCDF
            new_array <- array(NA, 
                             dim = c(dim(current_array)[1], 4, dim(current_array)[3]),
                             dimnames = list(
                                 NULL,
                                 c("Choice", "RT", "eCDF", "tCDF"),
                                 dimnames(current_array)[[3]]
                             ))
            
            # Copy existing Choice and RT data
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
                new_array[, 4, rep] <- all_tcdf[current_pos:(current_pos + n_trials - 1)]
                current_pos <- current_pos + n_trials
            }
            
            # Store updated array
            data_arrays[[param_name]][[n_trials_char]] <- new_array
        }
        cat("  (Done!)\n")
    }
    
    return(data_arrays)
}

# Usage:
 x <- add_cdfs_to_arrays(data_arrays, param_sets)