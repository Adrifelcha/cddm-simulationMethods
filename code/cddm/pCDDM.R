if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){ source("./dCDDM.R") }
library("scatterplot3d")
library("plot3D")

# Auxiliary function 1: A function to define the bins' ending points
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
define_bins <- function(tzero,obs.rt){
  ceiling(max(300, 20*(obs.rt-tzero)))
}

# Auxiliary function 2: Generate a plot of how the approximation works
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
embedded_plot <- function(bin.C, bin.RT, kappa, kappa.RT, total, 
                          rt, tzero, theta, drift, boundary,
                          density_matrix){
   # Base plot: Just draw the density
         nLines <- 30
         x.C <- seq(0,2*pi,length.out=nLines)
         y.RT <- seq(tzero,max(rt,10),length.out=nLines)
         x.theta <- rep(theta,nLines)
         z.Dens <- matrix(NA, nrow=nLines, ncol=nLines)
         for(c in 1:nLines){ for(t in 1:nLines){
           z.Dens[c,t] <- dCDDM(c(x.C[c],y.RT[t]),drift,theta,tzero,boundary)
         }}
         density.theta <- dCDDM(cbind(x.theta,y.RT),drift,theta,tzero,boundary)
         high_density <- max(c(z.Dens,density.theta))
         baseColor <- rgb(0,0,0,0.2)
         a <- scatterplot3d(x.C, y.RT, diag(z.Dens), zlab="Density", color = baseColor,
                            xlim=range(x.C), ylim=range(y.RT), type="l",
                            zlim=c(0,high_density), xlab="Choices", ylab="RT")
         a$points3d(x.C,rev(y.RT), diag(z.Dens[,c(nLines:1)]),type="l", col = baseColor)
         for(i in 1:nLines){  
           a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l", col = baseColor)
           a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i],  col = baseColor, type="l")
         }
         a$points3d(x.theta, y.RT, density.theta, col = "red", type="l")
         legend("topright", c("p( RT | theta )"), col="red", cex=0.6, lwd=1)
         
         nLines2 <- nLines*2
         k.C  <- floor(seq(1,kappa,length.out=nLines2))
         k.RT <- floor(seq(1,kappa.RT,length.out=nLines2))  
         for(i in k.C){  a$points3d(rep(bin.C[i],nLines2), bin.RT[k.RT], 
                                    density_matrix[i,k.RT], col="purple", type="l")
         }
         for(i in k.RT){ a$points3d(bin.C[k.C], rep(bin.RT[i],nLines2), 
                                    density_matrix[k.C,i], col="purple", type="l")
         }
         mtext(paste("Total =", round(max(total),4)),3, adj = 1, cex=0.8)  
 }

# Base function: We write a 2D Trapezoid N.I. algorithm
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
numInt.tpz.cddm <- function(rad,rt, cddm.par, plot=FALSE){
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
    kappa <- define_bins(tzero,rt)
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
    total <- total / (2*pi)
    
    if(plot){   embedded_plot(bin.C, bin.RT, kappa, kappa.RT, total, rt, 
                              tzero, theta, drift, boundary, density_matrix)     }
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
         volume <- numInt.tpz.cddm(rad = data[1], rt = data[2], cddm.par, plot)
      }else{
         volume <- numInt.tpz.cddm(rad = data[,1],rt = data[,2], cddm.par, plot)
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
    n <- 10
    C <- runif(n,0,2*pi)
    RT <- runif(n,0,15)
    data <- cbind(C,RT)
    pCDDM(data,drift, theta, tzero, boundary, plot=FALSE)
}