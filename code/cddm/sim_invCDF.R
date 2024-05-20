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

sample.invCDF.cddm <- function(n, par, plot=FALSE){
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
    
    cells <- expand.grid(space.C,space.RT)
    probs <-  numInt.tpz.cddm(cells[,1],cells[,2], par, plot=FALSE,
                              bin.RT = space.RT, bin.C = space.C)
  
    u <- runif(n,0,1)
    data <- matrix(NA, nrow=n, ncol=2)
    for(i in 1:n){
        look.match <- abs(probs-u[i])
        better.match <- as.character(which(look.match==min(look.match)))
        found.at <- as.numeric(sample(better.match,1))
        data[i,] <- as.numeric(as.vector(cells[found.at,]))
    }
    
    if(plot){
      par(mfrow=c(1,3),mar = c(0, 0, 0, 0)) 
      a <- scatterplot3d(cells[,1],cells[,2],probs, pch=16, zlim = c(0,1),
                         cex.symbols = 0.2, xlab = "", ylab = "", zlab = "",
                         xlim=range(cells[,1]), ylim=range(cells[,2]))
      a$points3d(data[,1], data[,2], u, col = "gray", pch=16, cex=0.5)
      mtext("Choices", side=1, line=0.5)
      mtext("RTs", side=4, line=-1, adj=0)
      plot(data[,1], u)
      plot(data[,2], u)
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