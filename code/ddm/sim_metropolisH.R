###############################################################################
###############################################################################
#####   A set of functions made to plot the data generated from simulations
###############################################################################
########################################################   by Adriana F. Chavez   
############## Load custom functions
library(RWiener)
###############################################################################
#####   1.2  Unidimensional Normal data: MH-MCMC sampling algorithm
###############################################################################
n <- 5000
par <- list("drift" = 1,
            "ndt" = 0.1,
            "boundary" = 3,
            "bias" = 0.5)

ddm_sim.MCMC <- function(n, par, plot=FALSE){
  # Load variables
  drift <- par$drift
  ndt <- par$ndt
  boundary <- par$boundary
  bias <- par$bias
  max.RT <- 10
  # Define maximum rejection value (RejVal) with respect to max. density point
  test.RT <- seq(ndt, max.RT, length.out = 500)
  test.responses <- c("upper", "lower")
  test.densities <- dwiener(rep(test.RT, 2), boundary, ndt, bias, drift, 
                            rep(test.responses, each=length(test.RT)))
  max.D <- max(test.densities)  # Density at mean (highest density point)
  height <- max.D*1.05          # RejVal ~ Uniform(0, 5% more than density at mean)
  # Define the range o candidate values to generate
  base.RT <- c(ndt,max.RT)
  ############################# !!!! ! ! ! ! ! ! ! ! !!!!!!!! !!!!! !! ! ! ! !
  # If plot is requested, draw the base Normal distribution
  if(plot){
  }
  # Start MH-MCMC sampling algorithm
  # This algorithm generates `n` samples from the specified Normal
  n.keep <- 0       # Start with 0 candidates stored
  samples <- NA     #   "    "  no     "        "
  n.try <- n        # Start by generating 'n' candidates
  while(n.keep < n){    # Until we reach a sample of size 'n'
    ######### First: Sample random candidates and rejectio values
    cand <- runif(n.try,base[1],base[2])    # We generate n.try candidate values
    eval <- dnorm(cand,mean,sd)           #    and compute their pdf
    rej.crit <- runif(n.try,0,height)     # We generate n.try Rej. Vals.
    # Keep the candidates whose density is greater than the corresponding Rej Val.
    keep <- (eval >= rej.crit)  
    ######### Second: Update
    n.pass <- sum(keep)      # Count number of accepted candidates in this run
    n.try <- n.try-n.pass    # Update number of samples missing to reach 'n'
    samples <- c(samples, cand[keep])   # Store candidate values approved
    n.keep <- length(samples)           # Count total no. of samples to test while() 
    ######### Third (Optional): Make plot
    if(plot){
      points(cand[!keep],rej.crit[!keep],col="red", pch=16, cex=0.3)
      points(cand[keep],rej.crit[keep],col="blue3", pch=16, cex=0.3)
    }
  }
  return(samples)
}

#~~~~~~~~~~~~~~~#
# Test/Example  #
#~~~~~~~~~~~~~~~#
if(!exists("test")){  test <- TRUE     }
if(test){
  par <- list("mean" = 10, "sd" = 1)
  n <- 5000
  sample.MCMC.normal(n,par, plot=TRUE)
}
