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
# method_tested is a global variable that can be set by the user
# before running the script. Alternatively, the user can set it here.
if(!exists("method_tested")){
    method_tested <- "Rejection_exGvonM"
}
# Possible methods:
# 1) "Metropolis"
# 2) "RandomWalk"
# 3) "inverseCDF"
# 4) "Rejection_Uniform"
# 5) "Rejection_exGvonM"
# 6) "Rejection_2DNormal"
cat("-----------------------------------------------------------\n")
cat("\n\nSimulation algorithm to be tested:", method_tested, "\n\n")
cat("-----------------------------------------------------------\n")
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
theta = pi
drift = c(0.5, 1.25, 3)
tzero = 0.1
boundary = 4

n = 5000    



algorithm_types <- c("Metropolis", "RandomWalk", "Rejection")
algorithms <- list("Metropolis" = c("2DNormal", "2DNormal_2"),
                   "Rejection_exGvonM" = rejection_exGvonM_algorithm)


for(a in length(algorithm_types)){
    for(b in length(algorithms[[a]])){
        for(c in length(drift)){
            
        }
    }
}