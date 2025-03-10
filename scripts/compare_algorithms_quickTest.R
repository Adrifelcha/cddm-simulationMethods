#########################################################################
# C O M P A R E     A L G O R I T H M S     Q U I C K     T E S T #######
# This script is used to obtain a quick comparison of the different 
# sampling algorithms.
# The *quick test* generates 'n' samples using each algorithm and 
# compares the runtime and output of each algorithm, using a range of 
# drift values.
#########################################################################
forceRun <- TRUE

#############################################################
# Load libraries and custom functions
#############################################################
library(here)
source(here("scripts", "load_libraries.R"))

#############################################################
#### S E T T I N G S ########################################
#############################################################
cat("-----------------------------------------------------------\n")
cat("Quick comparison of the different sampling algorithms\n")
cat("-----------------------------------------------------------\n")
# Algorithms to compare
#algorithm_types <- c("Metropolis", "RandomWalk", "Rejection", "inverseCDF")
algorithm_types <- c("RandomWalk", "Rejection")
#algorithms <- list("Metropolis" = c("Metropolis", "Metropolis_logRT"),
#                   "Rejection" = c("Rejection_exGvonM", "Rejection_Uniform", "Rejection_2DNormal"))
algorithms <- list("RandomWalk" = c("RandomWalk"),
                   "Rejection" = c("Rejection_Uniform", "Rejection_exGvonM", "Rejection_2DNormal"))

# Model parameters
theta = pi
drift = c(0.5, 1.25, 3)
tzero = 0.1
boundary = 4
# Trial size
n = 100    

cat("Sample parameter set:\n")
cat("theta = ", theta, "\n")
cat("tzero = ", tzero, "\n")
cat("boundary = ", boundary, "\n")
cat(length(drift), " drift values: ", paste(drift, collapse = ", "), "\n")
cat("-----------------------------------------------------------\n")
cat("We will generate ", n, " samples using each algorithm.\n")
cat("-----------------------------------------------------------\n\n\n\n")

#############################################################
#### R U N N I N G   T H E   A L G O R I T H M S ############
#############################################################
# Create a nested list structure to store samples by algorithm and drift
samples <- list()

# Prepare parameters to pass to the algorithms
par <- list(n = n, theta = theta, tzero = tzero, boundary = boundary)

for(a in 1:length(algorithm_types)){
    # Get the name of the algorithm type    
    name_type <- name_algorithmType(algorithm_types[a])
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    cat("Category:", name_type, "\n")
    cat("~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\n")
    
    # Initialize list for this algorithm type
    samples[[algorithm_types[a]]] <- list()

    if(algorithm_types[a] == "RandomWalk"){
        # Initialize array for RandomWalk algorithm with proper dimensions
        # n rows (trials), 2 columns (choice, rt), length(drift) pages (drift values)
        samples[[algorithm_types[a]]][["RandomWalk"]] <- array(NA, dim = c(n, 2, length(drift)))
        dimnames(samples[[algorithm_types[a]]][["RandomWalk"]]) <- list(NULL, c("choice", "rt"), paste0("drift_", 1:length(drift)))
        
        for(d in 1:length(drift)){
            cat("Drift value:", drift[d], "\n")    
            par$drift <- drift[d]
            Mu <- polarToRect(drift[d], par$theta)
            par$mu1 <- Mu$x
            par$mu2 <- Mu$y
            x <- rCDDM_RandomWalk(n = n, par = par)
            
            # Store the bivariate data in the corresponding slice of the array
            samples[[algorithm_types[a]]][["RandomWalk"]][,,d] <- as.matrix(x$bivariate.data)
        }
    } else {
        # For other algorithm types
        for(b in 1:length(algorithms[[a]])){
            alg_name <- algorithms[[a]][b]
            method_tested <- name_algorithm(alg_name)
            cat(method_tested, "\n")           
            
            # Initialize array for this specific algorithm
            samples[[algorithm_types[a]]][[alg_name]] <- array(dim = c(n, 2, length(drift)))
            dimnames(samples[[algorithm_types[a]]][[alg_name]]) <- list(NULL, c("choice", "rt"), paste0("drift_", drift))
            
            for(d in 1:length(drift)){
                cat("Drift value:", drift[d], "\n")
                par$drift <- drift[d]
                Mu <- polarToRect(drift[d], par$theta)
                par$mu1 <- Mu$x
                par$mu2 <- Mu$y
                
                full_params <- list(n = n, par = par)
                x <- run_algorithm(alg_name, full_params)                
                samples[[algorithm_types[a]]][[alg_name]][,,d] <- as.matrix(x$result)
                    
            }
        }
    }
}