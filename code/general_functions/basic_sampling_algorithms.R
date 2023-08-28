###############################################################################
###############################################################################
#####   A set of simplified NORMAL EXAMPLES for the sampling algorithms 
#####     we use as...
#####     1) Numeric integration of the cdf of both the DDM and CDDM models
#####     2) Rejection algorithms to sample observations for DDM and CDDM
###############################################################################
########################################################   by Adriana F. Chavez   
############## Load custom functions
source("../general_functions/customFunctions.R")


###############################################################################
#####   Unidimensional Normal data
###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To approximate CDFs, we use a Trapezoid Numerical Integration algorithm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some data
par <- list("mean" = 10, "sd" = 1)
lower.bound <- 1
upper.bound <- 10

# Write Trapezoid N.I. algorithm
numInt.tpz.normal <- function(lower.bound, upper.bound, par,
                              kappa=300, plot=FALSE){
    mean <- par$mean
    sd <- par$sd
    
    if(plot){
      width <- (3*sd)
      ticks <- seq(mean-width,mean+width,length.out=10)
      support <- seq(mean-width,mean+width,length.out=9999)
      plot(support,dnorm(support,mean,sd),
           type="l", ann=F, axes=F, col=1)
      axis(1,ticks,ticks)
      legend("topright", paste("No. bins =",kappa), 
             cex = 0.75, bty ="n")
    }
    
    bins <- seq(lower.bound,upper.bound,length.out=kappa)
    bin.area <- rep(NA,kappa-1)
    for(b in 2:kappa){
        low.x <- bins[b-1]
        up.x  <- bins[b]
        low.y <- dnorm(low.x,par$mean,par$sd)
        up.y  <- dnorm(up.x,par$mean,par$sd)
        height <- (low.y+up.y)/2
        bin.area[b-1] <- height*(up.x-low.x)
        
        if(plot){
          polygon(x=c(low.x,up.x,up.x,low.x),
                  y=c(0,0,low.y,up.y),
                  col = (b %% 2)+1)
        }
    }
    total.area <- sum(bin.area)
    return(total.area)
}

# Test function
numInt.tpz.normal(lower.bound,upper.bound,par, plot=TRUE)

# Use Trapezoid Numeric integration to compute CDF
normal.cdf <- function(x,par){
  lower.bound <- par$mean-(par$sd*100)
  area <- numInt.tpz.normal(lower.bound,x,par)
  return(area)
}

# Test function
normal.cdf(10,par)


###############################################################################
#####   Bivariate Normal data
#####   This section is merely illustrative, it requires R packages
###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To approximate CDFs, we use a Trapezoid Numerical Integration algorithm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some data
library(mnormt)

par <- list("mean.vector" = c(10,10),
            "varcov.mat" = matrix(c(2, -1, -1, 2), nrow=2))
lower.X    <- 1
upper.X     <- 10
lower.Y    <- 1
upper.Y     <- 10

# Support function
f <- function(x, y) dmnorm(cbind(x, y), Mean, Sigma)

# Write Trapezoid N.I. algorithm
numInt.tpz.normal <- function(lower.X, upper.X,
                              lower.Y, upper.Y,
                              par, kappa=300, plot=FALSE){
  Mean <- par$mean.vector
  Sigma <- par$varcov.mat
  
  if(plot){
    width <- NA
    ticks <- matrix(NA,nrow=2,ncol=10)
    support.def <- 70
    support <- matrix(NA,nrow=2,ncol=support.def)
    for(i in 1:2){
        width[i] <- (3*Sigma[i,i])
        ticks[i,] <- seq(Mean[i]-width[i],
                         Mean[i]+width[i],length.out=10)
        support[i,] <- seq(Mean[i]-width[i],
                           Mean[i]+width[i],length.out=support.def)
    }
    z <- outer(support[1,],support[2,], f)
    persp(support[1,], support[2,], z, col="gray99",
          theta=-30, phi=25, expand=0.6, ticktype='detailed')
    legend("topright", paste("No. bins =",kappa), 
           cex = 0.75, bty ="n")
  }
  
  bin.X <- seq(lower.X,upper.X,length.out=kappa)
  bin.Y <- seq(lower.Y,upper.Y,length.out=kappa)
  bin.area <- rep(NA,kappa-1)
  for(b in 2:kappa){
    a <- bin.X[b-1]
    b  <- bin.X[b]
    c <- bin.Y[b-1]
    d  <- bin.Y[b]
    low.y <- dnorm(low.x,par$mean,par$sd)
    up.y  <- dnorm(up.x,par$mean,par$sd)
    height <- (low.y+up.y)/2
    bin.area[b-1] <- height*(up.x-low.x)
    
    if(plot){
      polygon(x=c(low.x,up.x,up.x,low.x),
              y=c(0,0,low.y,up.y),
              col = (b %% 2)+1)
    }
  }
  total.area <- sum(bin.area)
  return(total.area)
}

# Test function
numInt.tpz.normal(lower.bound,upper.bound,par, plot=TRUE)

# Use Trapezoid Numeric integration to compute CDF
normal.cdf <- function(x,par){
  lower.bound <- par$mean-(par$sd*100)
  area <- numInt.tpz.normal(lower.bound,x,par)
  return(area)
}

# Test function
normal.cdf(10,par)