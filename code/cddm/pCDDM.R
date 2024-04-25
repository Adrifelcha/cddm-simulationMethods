test = FALSE
source("./dCDDM.R")
library("scatterplot3d")
library("plot3D")

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
  # Determine number of 2D bins to draw for current radian and rt
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
  
  #############################################
  # ~~ Compute the density at every vertix ~~ #
  #############################################
  # Part 1: Start empty storing objects
  density_mat_raw <- matrix(NA, nrow=kappa, ncol=kappa)
  # Part 2: Compute densities
  for(c in 1:kappa){
      for(t in 1:kappa){
          vertix <- c(bin.C[c],bin.RT[t])
          v.density <- dCDDM(vertix,drift,theta,tzero,boundary)
          density_mat_raw[c,t] <- v.density
      }
  }
  # Part 3: Densities computed from RT close to tzero are problematic 
  #         we locate them and remove them
  problem <- which(density_mat_raw < 0, arr.ind = T)
  bad.RT <- unique(problem[,2])
  density_matrix <- density_mat_raw[,-bad.RT]
  #         Update kappa upon removing the problematic RTs
  kappa.RT <- ncol(density_matrix)
  # Part 4: Plot
  if(plot){
      nLines2 <- nLines*2
      k.C  <- round(seq(1,kappa,length.out=nLines2),0)
      k.RT <- round(seq(1+length(bad.RT),kappa,length.out=nLines2),0)  
      for(i in k.C){  a$points3d(rep(bin.C[i],nLines2), bin.RT[k.RT], 
                                 density_mat_raw[i,k.RT], col="purple", type="l")
      }
      for(i in k.RT){ a$points3d(bin.C[k.C], rep(bin.RT[i],nLines2), 
                                 density_mat_raw[k.C,i], col="purple", type="l")
      }
  }
  
  ################################
  # ~~ Compute volume per bin ~~ #
  ################################
  # Part 1: Empty objects for storage
  bin.vol <- rep(NA,(kappa-1)*(kappa.RT-1))
  bin.area <- bin.vol
  bin.count <- 1
  for(b.c in 2:kappa){  # Move along radian dimension
      for(b.rt in 2:kappa.RT){  # Move along rt dimension
          bin.heights <- c(density_matrix[b.c,b.rt],
                           density_matrix[b.c-1,b.rt],
                           density_matrix[b.c,b.rt-1],
                           density_matrix[b.c-1,b.rt-1])
          mean_height <- mean(bin.heights)
          bin.vol[bin.count] <- mean_height*binBase.area
          # Alternative formula, exactly same result
          # VA = (1/6)*binBase.area*sum(bin.heights[-1])
          # VB = (1/6)*binBase.area*sum(bin.heights[-2])
          # VC = (1/6)*binBase.area*sum(bin.heights[-3])
          # VD = (1/6)*binBase.area*sum(bin.heights[-4])
          # bin.area[bin.count] <- (VA+VC+VB+VD)/2
          bin.count <- bin.count + 1
      }
  }
  total <- sum(bin.vol)
  if(plot){     
    mtext(paste("Total =", round(total,3)),3, adj = 1, cex=0.8)       
    }
  return(total)
}

# Test function
cddm.par <- list("drift" = 1, 
                 "theta" = pi,
                 "tzero" = 0.3,
                 "boundary" = )

numInt.tpz.cddm(rad=2*pi,rt=15, cddm.par, plot=TRUE)
numInt.tpz.cddm(rad=2*pi,rt=7, cddm.par, plot=TRUE)
numInt.tpz.cddm(rad=3*pi/2,rt=7, cddm.par, plot=TRUE)
numInt.tpz.cddm(rad=pi,rt=7, cddm.par, plot=TRUE)




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


