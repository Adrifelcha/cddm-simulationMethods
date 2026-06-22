#########################################################################
# N A M I N G    F U N C T I O N S ######################################
#########################################################################
name_algorithmType <- function(algorithm_type){
    if(algorithm_type == "Metropolis"){
        return("Metropolis sampling algorithm")
    } else if(algorithm_type == "RandomWalk"){
        return("RandomWalk emulation algorithm")
    } else if(algorithm_type == "Rejection"){
        return("Rejection sampling algorithm")
    } else if(algorithm_type == "inverseCDF"){
        return("Numerical approximation to the Probability transform")
    }
}

name_algorithm <- function(algorithm){
    if(algorithm == "Metropolis"){
        return("Proposal Distribution: Bivariate Normal")
    } else if(algorithm == "Metropolis_logRT"){
        return("Proposal Distribution: Bivariate Normal with log-transformed RTs")
    } else if(algorithm == "Rejection_exGvonM"){
        return("Proposal Distributions: exGaussian AND von Mises")
    } else if(algorithm == "Rejection_Uniform"){
        return("Proposal Distribution: Rectangular Uniform")
    } else if(algorithm == "Rejection_2DNormal"){
        return("Proposal Distribution: Bivariate Normal")
    }
}

########################

run_algorithm <- function(method_tested, full_params){
        #################################################################
    start_time <- Sys.time()
    if (method_tested == "Metropolis") {
        full_params$logRT <- FALSE
        result <- do.call(rCDDM_Metropolis, full_params)
    } else if (method_tested == "Metropolis_logRT") {
        full_params$logRT <- TRUE
        result <- do.call(rCDDM_Metropolis, full_params)
    } else if (method_tested == "inverseCDF") {
        result <- do.call(rCDDM_inverse, full_params)
    } else if (method_tested == "Rejection_Uniform") {
        full_params$type <- "Uniform"
        result <- do.call(rCDDM_Reject, full_params)
    } else if (method_tested == "Rejection_exGvonM") {
        full_params$type <- "exGvonM"
        result <- do.call(rCDDM_Reject, full_params)
    } else if (method_tested == "Rejection_2DNormal") {
        full_params$type <- "2DNormal"
        result <- do.call(rCDDM_Reject, full_params)
    } else {        stop(paste("Unknown method:", method_tested))       }
    end_time <- Sys.time()

    execution_time <- as.numeric(difftime(end_time, start_time, units="secs"))
    return(list(result = result, time = execution_time))
}