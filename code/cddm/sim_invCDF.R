###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####   INVERSE PROBABILITY TRANSFORM ALGORITHM with a grid approximation
###############################################################################
########################################################   by Adriana F. Chavez 
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){ 
      source("./pCDDM.R") 
      source("./sim_auxiliarFunctions.R")}
library(scatterplot3d) 

sample.invCDF.cddm <- function(n, par, plot=FALSE, color = NA){
  # Set-up - Load important values and define important variables
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Load parameter values
    drift <- par$drift;   theta <- par$theta
    tzero <- par$tzero;   boundary <- par$boundary
    
    minDensity <- 0.0009
    test.Density <- keyDensityPoints(par,cutoff = minDensity)
    min.RT <- test.Density$min.RT
    max.RT <- test.Density$max.RT
    predRT <- test.Density$pred.RT
    
    space.C  <- seq(0,2*pi,0.01)
    predRT.D <- cbind(space.C,rep(predRT,length(space.C)))
    test.D <- dCDDM(predRT.D,drift,theta,tzero,boundary)
    keep.Angles <- (minDensity < test.D)
    space.C <- space.C[keep.Angles]
    nBins.C <- length(space.C)
    
    space.RT <- seq(min.RT,max.RT,0.005)
    nBins.RT <- length(space.RT)
    
    Histogram3D <-  binsVolume_3Dhist(space.C, space.RT,
                                      drift, theta, tzero, boundary)
    probs <- Histogram3D$volume_per_bin

    bin_at.C  <- space.C[-length(space.C)] + diff(space.C)/2
    bin_at.RT <- space.RT[-length(space.RT)] + diff(space.RT)/2
    cells <- expand.grid(bin_at.C,bin_at.RT)
    grid <- cbind(cells,c(probs),cumsum(c(probs)))
  
    if(plot){     
      if(!is.function(color)){ color <- function(alpha){rgb(0,204/255,0,alpha)} }
      par(pty="m") 
      plot(grid[,4], ann=F, axes=F)
      axis(2, seq(0,1,0.1), seq(0,1,0.1), las=2)
      mtext("RT x Choice pairings",1,line=2, f=2)
      axis(1,c(1,nrow(grid)),c("min RT", "max RT"),line=-1)
      axis(1,c(1,nrow(grid)),c(min.RT, max.RT), line=0, tick = FALSE)
    }
    
    u <- runif(n,0,1)
    data <- matrix(NA, nrow=n, ncol=2)
    for(i in 1:n){
        look.match <- abs(grid[,4]-u[i])
        better.match <- as.character(which(look.match==min(look.match)))
        found.at <- as.numeric(sample(better.match,1))
        data[i,] <- as.numeric(as.vector(grid[found.at,c(1,2)]))
        if(plot){
           lines(c(1,found.at),c(u[i],u[i]), col=color(0.1), lwd=0.1)
           lines(c(found.at,found.at),c(0,u[i]), col = color(0.1), lwd=0.5)
           points(found.at,u[i],col=color(1), cex=0.5, pch=16)
        }
    }
    colnames(data) <- c("Choice", "RT")
return(data)
}


# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
    par <- list("drift" = 3.5, 
                "theta" = pi,
                "tzero" = 0.2,
                "boundary" = 2)
    
    max.RT=NA
    n <- 5000
    data <- sample.invCDF.cddm(n,par)
    hist(data[,2], main="Histogram of RTs", xlab="RT")
    
    eCDF.RT <- myECDF(data[,2])
    plot(data[,2],eCDF.RT)

}