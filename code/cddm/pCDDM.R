###############################################################################
###############################################################################
#####   Approximate the theoretical CDF function using the total volume
###############################################################################
########################################################   by Adriana F. Chavez 
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){     source("./dCDDM.R")       }
library("scatterplot3d")
library("plot3D")

# Auxiliary function 1: A function to define the bins' ending points
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
define_bins <- function(tzero,obs.rt){
  ceiling(max(350, 20*(obs.rt-tzero)))
}

# Auxiliary function 2: Generate a plot of how the approximation works
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
embedded_plot <- function(bin.C, bin.RT, kappa, kappa.RT, total, rad,
                          rt, tzero, theta, drift, boundary,density_matrix, 
                          rgb1 = c(0.2,0.5,0.6), rgb2 = c(0.3,0.4,0.7)){
   par(pty="s")          
   par(mfrow=c(1,2),mar = c(0, 0, 0, 0)) 
   # Plot 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Base plot: Sketch the bivariate density curve
         # Define line positions over the X and Y dimensions
         nLines <- 30
         x.C <- seq(0,2*pi,length.out=nLines)            # X: Choices
         y.RT <- seq(tzero,max(rt,10),length.out=nLines) # Y: RT
         x.theta <- rep(theta,nLines)                    # Make sure to add theta
         # Compute the density at each intersection point
         z.Dens <- matrix(NA, nrow=nLines, ncol=nLines)
         for(c in 1:nLines){ for(t in 1:nLines){
           z.Dens[c,t] <- dCDDM(c(x.C[c],y.RT[t]),drift,theta,tzero,boundary)
         }}
         # Compute the density at Choice = theta
         theta.Dens <- dCDDM(cbind(x.theta,y.RT),drift,theta,tzero,boundary)
         # Draw bivariate density curve
         baseColor <- rgb(0,0,0,0.2)
         a <- scatterplot3d(x.C, y.RT, diag(z.Dens), zlab="Density", color = baseColor, type="l",
                            xlab="Choices", ylab="RT", zlim = c(0, max(z.Dens, theta.Dens)))
         a$points3d(x.C,rev(y.RT), diag(z.Dens[,c(nLines:1)]),type="l", col = baseColor)
         for(i in 1:nLines){  
           a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l", col = baseColor)
           a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i],  col = baseColor, type="l")
         }
         # Add a density line corresponding to Choice = Theta
         a$points3d(x.theta, y.RT, theta.Dens, col = "red", type="l")
         legend("topright", c("p( RT | theta )"), col="red", cex=0.6, lwd=1, bty = "n")
   #  plot: Color the density area 
         nLines2 <- nLines*2
         k.C  <- floor(seq(1,kappa,length.out=nLines2))
         k.RT <- floor(seq(1,kappa.RT,length.out=nLines2))  
         for(i in k.C){  a$points3d(rep(bin.C[i],nLines2), bin.RT[k.RT], type="l", 
                                    density_matrix[i,k.RT], col = rgb(rgb1[1],rgb1[2],rgb1[3],0.5))
         }
         for(i in k.RT){ a$points3d(bin.C[k.C], rep(bin.RT[i],nLines2), type="l",
                                    density_matrix[k.C,i], col = rgb(rgb1[1],rgb1[2],rgb1[3],0.5))
         }
         for(i in 1:length(rt)){
           row <- which.min(abs(bin.C-rad[i]))-1
           col <- which.min(abs(bin.RT-rt[i]))-1
           a$points3d(rad[i], rt[i], density_matrix[row,col], pch=16, cex=0.3,
                      col=rgb(rgb2[1],rgb2[2],rgb2[3],1))
         }
         mtext(paste("Total =", round(max(total),4)),3, adj = 1, cex=0.8)  
   # Plot 2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
   # Base plot: Sketch the bivariate density curve
         b <- scatterplot3d(rad, rt, total,  pch=16, xlab = "", ylab = "", zlab = "",
                            color = rgb(rgb2[1],rgb2[2],rgb2[3],0.5), cex.symbols = 0.5)
         mtext("Approximate CDF: Area under curve", f=2)
         mtext("Choices", side=1, line=0.5)
         mtext("RTs", side=4, line=-1, adj=0)
 }

# Base function: We write a 2D Trapezoid N.I. algorithm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
numInt.tpz.cddm <- function(rad,rt, cddm.par, nBins = NA, plot=FALSE){
    # ~~ Set up ~~ #
    ################
    # Load up CDDM parameters
    drift <- cddm.par$drift;     theta <- cddm.par$theta
    tzero <- cddm.par$tzero;     boundary <- cddm.par$boundary
    # Identify number of data points
    n <- length(rad)
    # Make sure radians are in 0-2pi scale
    rad <- rad %% (2*pi)
    # Since 2pi = 0pi, we use 2*pi for simplicity ****
    if(sum(rad==0)>0){ rad[which(rad==0)] <- 2*pi}
    # Define the bins' cut points
    if(is.na(nBins)){    
          kappa <- define_bins(tzero,rt)
    }else{
          kappa <- nBins
    }
    # Take largest Choice and RT and partition to form the bins
    bin.C <- seq(0,max(rad),length.out=kappa)
    bin.RT <- seq(tzero,max(rt),length.out=kappa)

    # ~~ Base area of any bin ~~ #
    ##############################
    # All bins have the same length on both direction
    sidesLength.C <- bin.C[2]-bin.C[1]
    sidesLength.RT <- bin.RT[2]-bin.RT[1]
    # Compute the base area of all bins
    binBase.area <- sidesLength.C*sidesLength.RT
    
    # ~~ Compute the density at every vertex ~~ #
    #############################################
    # Part 1: Start empty storing objects
    density_matrix <- matrix(NA, nrow=kappa, ncol=kappa)
    # Part 2: Compute densities
    for(c in 1:kappa){ for(t in 1:kappa){
            vertex <- c(bin.C[c],bin.RT[t])
            v.density <- dCDDM(vertex,drift,theta,tzero,boundary)
            density_matrix[c,t] <- v.density
    }}
    # Part 3: Densities computed from RT close to tzero are problematic 
    #         we locate them and remove them
    problem <- which(density_matrix < 0, arr.ind = T)
    if(length(problem)!=0){
       bad.RT <- as.numeric(unique(problem[,2]))
       density_matrix <- density_matrix[,-bad.RT]
       bin.RT <- bin.RT[-bad.RT]
    }
    kappa.RT <- ncol(density_matrix)

    # ~~ Compute the height of each bin ~~ #
    ########################################
    # Part 1: Empty objects for storage
    sumHeight_per_bin <- matrix(NA, nrow=kappa-1, ncol=kappa.RT-1)
    thisC <- 1
    # Part 2: Fill height_per_bin matrix
    for(b.c in 2:kappa){  # Move along radian dimension
        thisRT <- 1
        for(b.rt in 2:kappa.RT){  # Move along rt dimension
            sumHeight_per_bin[thisC,thisRT] <- sum(c(density_matrix[b.c,b.rt],
                                                     density_matrix[b.c-1,b.rt],
                                                     density_matrix[b.c,b.rt-1],
                                                     density_matrix[b.c-1,b.rt-1]))
            thisRT <- thisRT + 1
        }
        thisC <- thisC +1
    }
    # ~~ Compute the volume of each bin ~~ #
    ########################################
    volume_per_bin <- sumHeight_per_bin*binBase.area*0.25
    
    # ~~ Compute the volume under each data point ~~ #
    ##################################################
    total <- rep(NA, n)
    for(i in 1:n){
        row <- which.min(abs(bin.C-rad[i]))-1
        col <- which.min(abs(bin.RT-rt[i]))-1
        total[i] <- sum(volume_per_bin[1:row,1:col])
    }
    total <- total
    
    if(plot){   
        embedded_plot(bin.C, bin.RT, kappa, kappa.RT, total, rad,
                      rt, tzero, theta, drift, boundary,density_matrix, 
                      rgb1 = c(0.2,0.5,0.6), rgb2 = c(0.3,0.4,0.7))     
      }
return(total)
}

# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
    cddm.par <- list("drift" = 1, 
                     "theta" = pi,
                     "tzero" = 0.3,
                     "boundary" = 5)
    numInt.tpz.cddm(rad=2*pi,rt=20, cddm.par, plot=TRUE)
    numInt.tpz.cddm(rad=2*pi,rt=7, cddm.par, plot=TRUE)
    numInt.tpz.cddm(rad=2*pi,rt=2, cddm.par, plot=TRUE)
    numInt.tpz.cddm(rad=3*pi/2,rt=7, cddm.par, plot=TRUE)
    numInt.tpz.cddm(rad=pi,rt=7, cddm.par, plot=TRUE)
}

# Use Trapezoid Numeric integration to compute CDF
pCDDM <- function(data,drift, theta, tzero, boundary, plot=FALSE){
      cddm.par <- list("drift" = drift, "theta" = theta, 
                       "tzero" = tzero, "boundary" = boundary)
      if(is.vector(data)){
         volume <- numInt.tpz.cddm(rad = data[1], rt = data[2], cddm.par, plot=plot)
      }else{
         volume <- numInt.tpz.cddm(rad = data[,1],rt = data[,2], cddm.par, plot=plot)
      }
  return(volume)
} 

# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
    drift = 1
    theta = pi
    tzero = 0.1
    boundary = 7
    # Test pCDDM with a single pair of observations
    C <- runif(1,0,2*pi)
    RT <- runif(1,0,15)
    data <- c(C,RT)
    pCDDM(data,drift, theta, tzero, boundary, plot=TRUE)
    # Test pCDDM with n pairs of observations
    n <- 1000
    C <- runif(n,0,2*pi)
    RT <- runif(n,0,15)
    data <- cbind(C,RT)
    pCDDM(data,drift, theta, tzero, boundary, plot=TRUE)
}