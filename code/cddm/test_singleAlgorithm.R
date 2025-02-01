# Function to run one benchmark test
single_algorithm_test <- function(params, n_trials, method_tested) {
    # Create a new list with n and the existing parameters
    full_params <- list(n = n_trials,
                        par = params$par)
    start_time <- Sys.time()

    ##################################################################
    # Run the sampling function
    ##################################################################
    if (method_tested == "RandomWalk") {
        result_raw <- do.call(sample.RW.cddm, full_params)
        result <- result_raw$bivariate.data
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

    ##################################################################
    # Calculate summary statistics
    ##################################################################
    # Execution time
    execution_time <- as.numeric(difftime(end_time, start_time, units="secs"))
    # Completion rate
    completion <- mean(!is.na(result$RT))    
     # Convert angles to circular data type and calculate circular mean
    angles <- circular::circular(result$Choice, type="angles", units="radians")
    mean_angle <- circular::mean.circular(angles)    
    # Calculate theoretical angle (theta) from mu parameters
    theoretical_theta <- circular::circular(atan2(params$par$mu2, params$par$mu1))    
    # Calculate angular error (without absolute value)
    angular_error <- as.numeric(mean_angle - theoretical_theta)    
    # Calculate proportion of negative RTs (excluding NAs)
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
        angular_error = angular_error)

    if(method_tested == "RandomWalk") {
        # Calculate distance from boundary for each trial
        final_coords <- getFinalState(result_raw$random.walk)
        distances_from_boundary <- sqrt(rowSums(final_coords^2))
        circumference_precision <- mean(abs(distances_from_boundary - params$par$boundary))
        output$circumference_precision <- circumference_precision
    }
    return(output)
}


