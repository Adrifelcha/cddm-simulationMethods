###############################################################################
###############################################################################
#####      A script containing auxiliary functions that are used repeatedly
#####               across different sampling algorithms
###############################################################################
########################################################   by Adriana F. Ch?vez 
#if(!exists("superCalled")){superCalled <- FALSE}
#if(!superCalled){ source("./dCDDM.R") }
library(mvtnorm)



###############################################################################
# Transformation functions: ###################################################
###############################################################################
# Switch between Cardinal and Rectangular Coordinates 
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
rectToPolar <- function(x,y){
  n <- length(x)
  driftAngle <- atan2(y,x)
  driftLength <- sqrt((x^2)+(y^2))
  output <- as.data.frame(cbind(driftAngle,driftLength))
  colnames(output) <- c("dAngle","dLength")
  return(output)
}

polarToRect <- function(vectorAngle,vectorLength){
  x <- vectorLength*cos(vectorAngle)
  y <- vectorLength*sin(vectorAngle)
  X <-  as.data.frame(cbind(x,y))
  colnames(X) <-  c("x","y")
  return(X)
}


# Switch between degrees and radians
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
degToRad <- function(theta.deg){  
  theta <-  theta.deg * pi /180  #Transform to radians
  return(theta)
}

radToDeg <- function(theta.rad){
  theta <- theta.rad * (180/pi)
  return(theta)
}


###############################################################################
# Sampling space functions: ###################################################
###############################################################################

### Use CDDM density function to find key values for setting up the sampling space
keyDensityPoints <- function(par, cutoff = 0.00009){
  # Extract model parameters from input list
  drift <- par$drift;   theta <- par$theta      # Drift magnitude and direction
  tzero <- par$tzero;   boundary <- par$boundary # Non-decision time and boundary
  
  ### Find minimum response time (RT) where density exceeds cutoff
  # Start from non-decision time plus small increment
  density <- 0
  min.RT <- tzero+0.1
  while(density < cutoff){
    # Test density at drift direction (theta) until exceeding cutoff
    density <- dCDDM(c(theta,min.RT),drift,theta,tzero,boundary)
    min.RT <- min.RT+0.01
  }

  ### Find maximum density and most likely RT
  # Search for peak density by incrementing RT until density starts decreasing
  density_increase <- TRUE; d <- 0
  pred.RT <- min.RT
  while(density_increase){
    pred.RT <- pred.RT+0.01
    density <- dCDDM(c(theta,pred.RT),drift,theta,tzero,boundary)
    density_increase <- density > d
    d <- density
  }
  max.Density <- density

  ### Find maximum RT where density falls below cutoff
  # Start from predicted RT and increment until density drops below cutoff
  max.RT <- pred.RT
  while(density > cutoff){
    density <- dCDDM(c(theta,max.RT),drift,theta,tzero,boundary)
    max.RT <- max.RT+0.01
  }

  # Return key points needed for setting up rejection sampling
  return(list("min.RT" = min.RT,          # Minimum response time
              "pred.RT" = pred.RT,         # Most likely response time
              "max.Density" = max.Density, # Maximum density value
              "max.RT" = max.RT))         # Maximum response time
}

###############################################################################
# EZ CDDM functions: ##########################################################
###############################################################################

ezcddm_MRT <- function(drift, boundary, tzero) {
    # Calculate av product
    av <- drift * boundary
    
    # Calculate ratio of modified Bessel functions
    I1 <- besselI(av, nu = 1)  # First kind, order 1
    I0 <- besselI(av, nu = 0)  # First kind, order 0
    
    # Calculate expected RT
    MRT <- tzero + (boundary/drift) * (I1/I0)
    
    return(MRT)
}


ezcddm_VRT <- function(drift, boundary) {  
    # Calculate av product
    av <- drift * boundary
    
    # Calculate modified Bessel functions
    I1 <- besselI(av, nu = 1)  # First kind, order 1
    I0 <- besselI(av, nu = 0)  # First kind, order 0
    
    # Calculate ratio once to avoid repetition
    ratio <- I1/I0
    
    # Calculate variance using the full equation:
    # ((b^2/v^2)(I1^2/I0^2)) + ((2b/v^3)(I1/I0)) - (b^2/v^2)
    VRT <- ((boundary^2)/(drift^2)) * (ratio^2) + 
           ((2*boundary)/(drift^3)) * ratio - 
           ((boundary^2)/(drift^2))
    
    return(VRT)
}


ezcddm_VCA <- function(drift, boundary, tzero) {
    # Calculate av product
    av <- drift * boundary
    
    # Calculate ratio of modified Bessel functions
    I1 <- besselI(av, nu = 1)  # First kind, order 1
    I0 <- besselI(av, nu = 0)  # First kind, order 0
    
    # Calculate expected RT
    VCA <- 1 - (I1/I0)
    
    return(VCA)
}

###############################################################################
# Ex-Gaussian functions: ######################################################
###############################################################################

rexGAUS <- function(n, mu, sigma, tau) {
  return(rnorm(n, mu, sigma) + rexp(n, 1/tau))
}

dexGAUS <- function(x, mu, sigma, tau) {
  if(tau <= 0) return(dnorm(x, mu, sigma))
  
  arg1 <- sigma^2/(2*tau^2) + (mu-x)/tau
  arg2 <- (x-mu)/sigma - sigma/tau
  
  return(exp(arg1) * pnorm(arg2) / tau)
}

