#############################################################
# Load libraries and custom functions
#############################################################
cat("Loading R libraries...\n")
library("here")
library("circular")
library("mvtnorm")
library("scatterplot3d") 

cat("\nLoading custom function scripts from /code/cddm...\n\n")
source(here("src-code", "cddm", "sim_randomWalk.R"))        
r_files <- list.files(path = here("src-code", "cddm"), 
                      pattern = "\\.R$", 
                      full.names = TRUE)
for(file in r_files) {
    source(file)
}
source(here("src-code", "general_functions", "eCDF.R"))