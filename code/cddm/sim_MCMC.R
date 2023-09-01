###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                   MCMC NUMERIC INTEGRATION ALGORITHM
###############################################################################
########################################################   by Adriana F. Ch?vez 
source("./dCDDM.R")
library(scatterplot3d) 

# Some data
par <- list("drift" = 1, 
            "theta" = pi,
            "tzero" = 0.1,
            "boundary" = 7)
n <- 5000

# Write MCMC algorithm for the CDDM pdf
sample.MCMC.cddm <- function(n, par, max.RT = 10, plot=FALSE){
  no.Dim <- 2
  drift <- par$drift
  theta <- par$theta
  tzero <- par$tzero
  boundary <- par$boundary

  test.RT <- seq(tzero, max.RT, length.out=10)
  test.densities <- NA
  for(i in 1:length(test.RT)){
      test.densities[i] <- dCDDM(c(theta,test.RT[i]),drift, theta, tzero, boundary)
  }
  max.Density <- max(test.densities)
  height <- max.Density*1.2
  base.C <- c(0, 2*pi)
  base.RT <- c(0, max.RT)
  
  if(plot){
    nSupp <- 100
    support.C <- seq(base.C[1],base.C[2],length.out=nSupp)
    support.RT <- seq(base.RT[1],base.RT[2],length.out=nSupp)
    
    z <- dCDDM(cbind(support.C,support.RT),drift, theta, tzero, boundary)
    a <- scatterplot3d(support.C, support.RT, z, 
                       xlim=base.C, ylim=base.RT, zlim=c(0,height),
                       color="blue", type="l", lwd=2)
  }
  
  n.keep <- 0
  n.try <- n
  samples <- matrix(NA, nrow=1, ncol=no.Dim)
  
  while(n.keep < n){
    cand <- matrix(NA, nrow=n.try, ncol=no.Dim)
    cand[,1] <- runif(n.try,base.C[1],base.C[2])
    cand[,2] <- runif(n.try,base.RT[1],base.RT[2])
    
    eval <- dCDDM(cand,drift, theta, tzero, boundary)
    rej.crit <- runif(n.try,0,height)  
    keep <- (eval >= rej.crit)
    
    n.keep <- sum(keep)
    n.try <- n.try-n.keep
    
    if(plot){
      a$points3d(cand[!keep,1], cand[!keep,2], rej.crit[!keep],
                 col = "red", pch = 16, cex = 0.2)
      a$points3d(cand[keep,1], cand[keep,2], rej.crit[keep],
                 col = "blue3", pch = 16, cex = 0.2)
    }
    
    samples <- rbind(samples, cand[keep,])
    n.keep <- nrow(samples)-1
  }
  
  samples <- samples[-1,]
  return(samples)
}

# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){              sample.MCMC.cddm(1000,par, plot=TRUE)  }