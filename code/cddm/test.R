###############################################################################
###############################################################################
#####      A script to call, test and compare every simulation method
###############################################################################
########################################################   by Adriana F. Ch?vez 
# 1: Load relevant scripts and libraries ~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
test <- FALSE   # Turn-off built-in tests
source("./sim_MCMC.R")
source("./sim_randomWalk.R")

# 2: Define some parameter values and sample size ~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# Some parameter values
par <- list("drift" = 1, 
            "theta" = pi,
            "tzero" = 0.1,
            "boundary" = 7)
n <- 5000 # No. samples

# 3: ~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
X.RW <- sample.RW.cddm(n,par)
max.RT <- max(X.RW$bivariate.data[,2])
X.MCMC <- sample.MCMC.cddm(n,par,max.RT)

# 4: Plot and compare the data sampled ~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

# 5: Run accuracy tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


