library(here)

# Running the random walk emulator
# requires the following src-code files:
source(here("src-code", "cddm", "sim_randomWalk.R"))

###################################################
# 
###################################################
trials <- 1000
mu1 <- 1
mu2 <- 1
boundary <- 1
tzero <- 0.1

par <- list(mu1 = mu1, mu2 = mu2, boundary = boundary, tzero = tzero)

random_walk <- cddm.randomWalk(trials = trials, mu1 = mu1, mu2 = mu2, boundary = boundary)  

final_state <- getFinalState(randomWalk.states)


X <- rCDDM_RandomWalk(n = trials, par = par)