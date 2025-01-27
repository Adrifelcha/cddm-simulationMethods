# Function to run one benchmark test
run_simple_test <- function(params, n_trials, method_tested) {
    # Create a new list with n and the existing parameters
    full_params <- list(n = n_trials,
                        par = params$par)
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
    # Execution time
    execution_time <- as.numeric(difftime(end_time, start_time, units="secs"))
    # Completion rate
    completion <- mean(!is.na(result$bivariate.data$RT))
    # Calculate distance from boundary for each trial
    final_coords <- getFinalState(result$random.walk)
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




# Function to run one benchmark test
run_full_test <- function(params, n_trials, method_tested) {
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
    
    # Calculate angular error (without absolute value)
    angular_error <- as.numeric(mean_angle - theoretical_theta)
    
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
