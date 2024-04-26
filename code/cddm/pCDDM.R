if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){ source("./dCDDM.R") }
library("scatterplot3d")
library("plot3D")

# Write Trapezoid N.I. algorithm
numInt.tpz.cddm <- function(rad,rt, cddm.par, plot=FALSE){
    # ~~ Set up ~~ #
    ################
    # Load up CDDM parameters
    drift <- cddm.par$drift;     theta <- cddm.par$theta
    tzero <- cddm.par$tzero;     boundary <- cddm.par$boundary
    # Make sure radians are in 0-2pi scale
    rad <- rad %% (2*pi)
    # Since 2pi = 0pi, we use 2*pi for simplicity ****
    if(rad==0){ rad <- 2*pi}
    # Define number of 2D bins to draw based on choice and RT range
    kappa <- ceiling(max(rad,(rt-tzero))*20)
    # Define bin ends (i.e. vertices)
    bin.C <- seq(0,rad,length.out=kappa)
    bin.RT <- seq(tzero,rt,length.out=kappa)
    
    # ~~ If requested, plot base density ~~ #
    #########################################
    if(plot){
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
            a <- scatterplot3d(x.C, y.RT, diag(z.Dens), zlab="Density",
                               xlim=range(x.C), ylim=range(y.RT), type="l",
                               zlim=c(0,high_density), xlab="Choices", ylab="RT")
            a$points3d(x.C, rev(y.RT), diag(z.Dens[,c(nLines:1)]),type="l")
            for(i in 1:nLines){  
                a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l")
                a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i], type="l")
            }
            a$points3d(x.theta, y.RT, density.theta, col = "red", type="l")
            legend("topright", c("p( RT | theta )"), col="red", cex=0.6, lwd=1)
    }
    
    # ~~ Base area of any bin ~~ #
    ##############################
    # All bins have the same length on both direction
    sidesLength.C <- bin.C[2]-bin.C[1]
    sidesLength.RT <- bin.RT[2]-bin.RT[1]
    # Compute the base area of all bins
    binBase.area <- sidesLength.C*sidesLength.RT
    
    # ~~ Compute the density at every vertix ~~ #
    #############################################
    # Part 1: Start empty storing objects
    density_mat_raw <- matrix(NA, nrow=kappa, ncol=kappa)
    # Part 2: Compute densities
    for(c in 1:kappa){ for(t in 1:kappa){
            vertix <- c(bin.C[c],bin.RT[t])
            v.density <- dCDDM(vertix,drift,theta,tzero,boundary)
            density_mat_raw[c,t] <- v.density
    }}
    # Part 3: Densities computed from RT close to tzero are problematic 
    #         we locate them and remove them
    problem <- which(density_mat_raw < 0, arr.ind = T)
    if(length(problem)!=0){
       bad.RT <- unique(problem[,2])
       density_matrix <- density_mat_raw[,-bad.RT]
    }else{     density_matrix <- density_mat_raw       }
    kappa.RT <- ncol(density_matrix)
    # Part 4: Plot
    if(plot){
        nLines2 <- nLines*2
        k.C  <- floor(seq(1,kappa,length.out=nLines2))
        k.RT <- floor(seq(1+length(bad.RT),kappa.RT,length.out=nLines2))  
        for(i in k.C){  a$points3d(rep(bin.C[i],nLines2), bin.RT[k.RT], 
                                   density_mat_raw[i,k.RT], col="purple", type="l")
        }
        for(i in k.RT){ a$points3d(bin.C[k.C], rep(bin.RT[i],nLines2), 
                                   density_mat_raw[k.C,i], col="purple", type="l")
        }
    }
    
    # ~~ Compute volume per bin ~~ #
    ################################
    # Part 1: Empty objects for storage
    bin.vol <- rep(NA,(kappa-1)*(kappa.RT-1))
    bin.count <- 1
    for(b.c in 2:kappa){  # Move along radian dimension
        for(b.rt in 2:kappa.RT){  # Move along rt dimension
            bin.heights <- c(density_matrix[b.c,b.rt],
                             density_matrix[b.c-1,b.rt],
                             density_matrix[b.c,b.rt-1],
                             density_matrix[b.c-1,b.rt-1])
            mean_height <- mean(bin.heights)
            bin.vol[bin.count] <- mean_height*binBase.area
            bin.count <- bin.count + 1
        }
    }
    total <- sum(bin.vol) / (2*pi)
    if(plot){ mtext(paste("Total =", round(total,4)),3, adj = 1, cex=0.8)  }
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
  rad <- data[1]
  rt <- data[2]
  cddm.par <- list("drift" = drift, "theta" = theta, 
                   "tzero" = tzero, "boundary" = boundary)
  area <- numInt.tpz.cddm(rad,rt, cddm.par, plot)
  return(area)
} 

# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
    drift = 1
    theta = pi
    tzero = 0.1
    boundary = 7
    C <- runif(1,0,2*pi)
    RT <- runif(1,0,15)
    data <- c(C,RT)
    pCDDM(data,drift, theta, tzero, boundary, plot=TRUE)
}