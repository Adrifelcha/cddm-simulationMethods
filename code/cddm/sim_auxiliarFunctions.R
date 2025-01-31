###############################################################################
###############################################################################
#####      A script containing auxiliary functions that are used repeatedly
#####               across different sampling algorithms
###############################################################################
########################################################   by Adriana F. Ch?vez 
#if(!exists("superCalled")){superCalled <- FALSE}
#if(!superCalled){ source("./dCDDM.R") }
library(mvtnorm)


### Use CDDM density function to find key values
keyDensityPoints <- function(par, cutoff = 0.00009){
  drift <- par$drift;   theta <- par$theta
  tzero <- par$tzero;   boundary <- par$boundary
  
  ### Determine a reasonable minimum RT testing the density at theta
  density <- 0
  min.RT <- tzero+0.1
  while(density < cutoff){
    density <- dCDDM(c(theta,min.RT),drift,theta,tzero,boundary)
    min.RT <- min.RT+0.01
  }
  ### Find the maximum density and most likely RT testing density at theta
  density_increase <- TRUE; d <- 0
  pred.RT <- min.RT
  while(density_increase){
    pred.RT <- pred.RT+0.01
    density <- dCDDM(c(theta,pred.RT),drift,theta,tzero,boundary)
    density_increase <- density > d
    d <- density
  }
  max.Density <- density
  ### Determine a reasonable maximum RT testing the density at theta
  max.RT <- pred.RT
  while(density > cutoff){
    density <- dCDDM(c(theta,max.RT),drift,theta,tzero,boundary)
    max.RT <- max.RT+0.01
  }
return(list("min.RT" = min.RT,
            "pred.RT" = pred.RT,
            "max.Density" = max.Density,
            "max.RT" = max.RT))
}