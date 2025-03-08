#########################################################################
# S I N G L E     A L G O R I T H M     Q U I C K     T E S T ##########
# This script is used to quickly test any individual sampling algorithm
# contained in this repo.
# The *quick test* uses few trial sizes and only two parameter sets:
# one with strong drift and one with weak drift.
# This script is designed to be run from the command line with a 
# single argument specifying the algorithm to test.
#########################################################################
forceRun <- TRUE

#############################################################
# Load libraries and custom functions
#############################################################
cat("Loading R libraries...\n")
library("here")
library("circular")

cat("\nLoading custom function scripts from /code/cddm...\n\n")
source(here("code", "cddm", "sim_randomWalk.R"))        
r_files <- list.files(path = here("code", "cddm"), 
                      pattern = "\\.R$", 
                      full.names = TRUE)
for(file in r_files) {
    source(file)
}
source(here("code", "general_functions", "eCDF.R"))

#############################################################
#### S E T T I N G S ########################################
#############################################################
cat("-----------------------------------------------------------\n")
cat("\n\nQuick comparison of the different sampling algorithms\n")
cat("-----------------------------------------------------------\n")

theta = pi
drift = c(0.5, 1.25, 3)
tzero = 0.1
boundary = 4

n = 5000    

cat("\n\nSample parameter set:\n")
cat("theta = ", theta, "\n")
cat("tzero = ", tzero, "\n")
cat("boundary = ", boundary, "\n")
cat(length(drift), " drift values: ", paste(drift, collapse = ", "), "\n")
cat("-----------------------------------------------------------\n")

cat("\n\nWe'll generate ", n, " samples using each algorithm.\n")
cat("-----------------------------------------------------------\n\n\n\n")



algorithm_types <- c("Metropolis", "RandomWalk", "Rejection")
algorithms <- list("Metropolis" = c("2DNormal", "2DNormal_2"),
                   "Rejection" = c("Rejection_exGvonM", "Rejection_Uniform", "Rejection_2DNormal"))

samples <- array(dim = c(n, 2, length(unlist(algorithms))))

for(a in length(algorithm_types)){

    name_type <- name_algorithmType(algorithm_types[a])
    cat("-----------------------------------------------------------\n")
    cat("\n\nCategory:", name_type, "\n\n")
    cat("-----------------------------------------------------------\n")

    if(algorithm_types[a] == "RandomWalk"){
        for(d in length(drift)){
            cat("Drift value:", drift[d], "\n")                        
        }

    }else{
            for(b in length(algorithms[[a]])){
                for(c in length(drift)){
                    
                }
            }
    }
}