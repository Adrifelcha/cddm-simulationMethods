test = FALSE
source("./dCDDM.R")
library("scatterplot3d")
library("plot3D")

# Some data
cddm.par <- list("drift" = 1, 
            "theta" = pi,
            "tzero" = 0.1,
            "boundary" = 7)
rad <- 2*pi
rt <- 20
plot = TRUE

# Write Trapezoid N.I. algorithm
numInt.tpz.cddm <- function(rad,rt, cddm.par, plot=FALSE){
  ################
  # ~~ Set up ~~ #
  ################
  # Initial variables
  no.Dim <- 2
  max.RT <- 20
  # Load up CDDM parameters
  drift <- cddm.par$drift
  theta <- cddm.par$theta
  tzero <- cddm.par$tzero
  boundary <- cddm.par$boundary
  # Make sure radians are in 0-2pi scale
  rad <- rad %% (2*pi)
  # Since 2pi = 0pi, we use 2*pi for simplicity ****
  if(rad==0){ rad <- 2*pi}
  # Determine how far from lower rt we are
  range.RT <- rt-tzero
  # Determine number of 2D bins to draw based on
  # the radian and rt specified
  kappa <- max(rad,range.RT)*20
  # Define bins along the Radian dimension (2pi works better****)
  bin.C <- seq(0,rad,length.out=kappa)
  # Define bins along the RT dimension
  bin.RT <- seq(tzero,rt,length.out=kappa)
  
  # If requested, draw the base 3D density curve
  if(plot){
            nSupp <- 50 
            nLines <- 30
            base.C <- c(0, 2*pi)
            base.RT <- c(tzero,rt)
            support.C <- seq(base.C[1],base.C[2],length.out=nSupp)
            support.RT1 <- seq(base.RT[1],base.RT[2],length.out=nSupp)
            support.RT2 <- rev(support.RT1)
            support.theta <- rep(theta,nSupp)
            z.diag1 <- dCDDM(cbind(support.C,support.RT1),drift, theta, tzero, boundary)
            z.diag2 <- dCDDM(cbind(support.C,support.RT2),drift, theta, tzero, boundary)
            z.RT_at_theta <- dCDDM(cbind(support.theta,support.RT1),drift, theta, tzero, boundary)
            z.top <- max(c(z.diag1,z.diag2,z.RT_at_theta))
            a <- scatterplot3d(support.C, support.RT1, z.diag1, 
                               xlim=base.C, ylim=base.RT, zlim=c(0,z.top),
                               xlab="Choices", ylab="RT", zlab="Density",
                               color="blue", type="l")
            a$points3d(support.C, support.RT2, z.diag2, col = "blue", type="l")
            a$points3d(support.theta, support.RT1, z.RT_at_theta, col = "red", type="l")
            L <- round(nSupp/nLines,0)
            for(i in 1:nLines){
                  choose.RT <- rep(support.RT1[i*L],nSupp)
                  choose.C  <- rep(support.C[i*L],nSupp)
                  z.overRT <- dCDDM(cbind(choose.C,support.RT1),drift, theta, tzero, boundary)
                  z.overC <-  dCDDM(cbind(support.C,choose.RT),drift, theta, tzero, boundary)
                  a$points3d(support.C, choose.RT, z.overC, col = "blue", type="l")
                  a$points3d(choose.C, support.RT1, z.overRT, col = "blue", type="l")
            }
            legend("topright", c("p( RT | theta )"), col="red", cex=0.6, lwd=1)
  }
  
  ##############################
  # ~~ Base area of any bin ~~ #
  ##############################
  # All bins have the same length on both direction
  sidesLength.C <- bin.C[2]-bin.C[1]
  sidesLength.RT <- bin.RT[2]-bin.RT[1]
  # Compute the base area of all bins
  binBase.area <- sidesLength.C*sidesLength.RT
  
  ################################
  # ~~ Compute volume per bin ~~ #
  ################################
  # Initial values / Empty objects for storage
  bin.vol <- rep(NA,(kappa-1)^2)
  bin.count <- 1
  for(b.c in 2:kappa){  # Move along radian dimension
      for(b.rt in 2:kappa){  # Move along rt dimension
          C.space  <- c(bin.C[b.c-1],bin.C[b.c])
          RT.space <- c(bin.RT[b.rt-1],bin.RT[b.rt])
          v.height <- matrix(rep(NA, 4),nrow=2)
          for(t in 1:2){
              for(c in 1:2){
              vertix <- c(C.space[c],RT.space[t])
              v.density <- dCDDM(vertix,drift,theta,tzero,boundary)
              store <- max(c(0,v.density))
              v.height[c,t] <- store
                  if(plot){
                    a$points3d(vertix[1],vertix[2],store,
                               col="pink",type="p",pch = 16, cex = 0.2)
                  }
              }
          }
          height <- mean(v.height)
          bin.vol[bin.count] <- height*binBase.area
          # Alternative formula, exactly same result
          # VA = (1/6)*base.area*sum(as.vector(v.height)[-1])
          # VB = (1/6)*base.area*sum(as.vector(v.height)[-2])
          # VC = (1/6)*base.area*sum(as.vector(v.height)[-3])
          # VD = (1/6)*base.area*sum(as.vector(v.height)[-4])
          # bin.area[bin.count] <- (VA+VC+VB+VD)/2
          bin.count <- bin.count + 1
      }
  }
  total <- sum(bin.vol)
  if(plot){     
    mtext(paste("total =", round(total,3)),3, outer = TRUE, adj = 1)       
    }
  return(total.area)
}

# Test function
numInt.tpz.cddm(lower.C, upper.C,
                lower.RT, upper.RT,
                par = par, max.RT = 20, 
                kappa=300, plot=FALSE)

# Use Trapezoid Numeric integration to compute CDF
pCDDM <- function(data,drift, theta, tzero, boundary, plot=FALSE){
  lower.C <- 0
  lower.RT <- tzero
  upper.C <- data[1]
  upper.RT <- data[2]
  par <- list("drift" = drift, "theta" = theta, 
              "tzero" = tzero, "boundary" = boundary)
  area <- numInt.tpz.cddm(lower.C, upper.C,
                          lower.RT, upper.RT,
                          par, max.RT = 20, 
                          kappa=300, plot=FALSE)
  return(area)
} 

# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
    drift = 1
    theta = pi
    tzero = 0.1
    boundary = 7
    n <- 500
    C <- runif(1,0,2*pi)
    RT <- runif(1,0,15)
    data <- c(C,RT)
    pCDDM(data,drift, theta, tzero, boundary, plot=TRUE)
}