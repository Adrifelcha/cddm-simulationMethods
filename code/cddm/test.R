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
source("./plottingFunctions.R")

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
X.invCDF <- sample.invCDF.cddm(n,par,max.RT)

# 4: Plot and compare the data sampled ~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
par(pty="m")          
par(mfrow=c(2,1),mar = c(3, 3, 3, 0)) 
# Choices
hist(X.RW$bivariate.data[,1], main="Choices on Random Walk",
     xlim=c(0,2*pi), col="goldenrod4")
hist(X.MCMC[,1], main= "Choices on MCMC", xlim=c(0,2*pi), col="cyan4")
# Response Times
hist(X.RW$bivariate.data[,2], main = "RT on Random Walk", col="goldenrod3")
hist(X.MCMC[,2], main = "RT on MCMC", col="cyan4")
# Decisions on a circle
plot.CDDM(X.MCMC, par, choice.col.RGB = c(0.15,.29,.80))
plot.CDDM(X.RW, par, choice.col.RGB = c(0.4,.9,.3))

# 5: Run accuracy tests ~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #


