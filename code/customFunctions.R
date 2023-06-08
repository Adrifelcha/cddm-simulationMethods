###############################################################################
###############################################################################
#####   A set of functions made for the simulation paper
###############################################################################
########################################################   by Adriana F. Chávez   


###############################################################################
###################   Functions related to the CDDM
###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Switch between Cardinal and Rectangular Coordinates 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Switch between degrees and radians
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
degToRad <- function(theta.deg){  
  theta <-  theta.deg * pi /180  #Transform to radians
  return(theta)
}

radToDeg <- function(theta.rad){
  theta <- theta.rad * (180/pi)
  return(theta)
}
