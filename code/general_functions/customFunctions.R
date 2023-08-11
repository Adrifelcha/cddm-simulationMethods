###############################################################################
###############################################################################
#####   A set of functions made for the simulation paper
###############################################################################
########################################################   by Adriana F. Ch?vez   



###############################################################################
###################   General Simulation Functions
###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Take full random walk coordinates and extract final states
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getFinalState <- function(randomWalk.states){
  randomWalk <- randomWalk.states
  dimensions <- dim(randomWalk)
  K <- nrow(randomWalk)
  
  if(length(dimensions)>2){
      I <- dimensions[3]
      coord <- matrix(NA, ncol=2,nrow=I)
      for(i in 1:I){
          for(k in 1:K){
              if(!is.na(randomWalk[k,1,i])){
                  a <- k
              }else{
                  break
              }
          }
          coord[i,] <- randomWalk[a,,i]
      }
      output <- coord
  }else{
      I <- dimensions[2]
      choice <- rep(NA, I)
      for(i in 1:I){
        for(k in 1:K){
          if(!is.na(randomWalk[k,i])){
            a <- k
          }else{
            break
          }
        }
        choice[i] <- randomWalk[a,i]
      }
      output <- choice
  }
  return(output)
}

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






