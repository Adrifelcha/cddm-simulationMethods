name_algorithmType <- function(algorithm_type){
    if(algorithm_type == "Metropolis"){
        return("Metropolis sampling algorithm")
    } else if(algorithm_type == "RandomWalk"){
        return("RandomWalk emulation algorithm")
    } else if(algorithm_type == "Rejection"){
        return("Rejection sampling algorithm")
    } else if(algorithm_type == "invCDF"){
        return("Numerical approximation to the Probability transform")
    }
}