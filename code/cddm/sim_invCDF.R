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
    
    test.Density <- keyDensityPoints(par,cutoff = 0.0009)
    min.RT <- test.Density$min.RT
    max.RT <- test.Density$max.RT
    max.Density <- test.Density$max.Density
    predRT <- test.Density$pred.RT
    
    space.C  <- seq(0,2*pi,0.02)
    test <- cbind(space.C,rep(predRT,length(space.C)))
    test.D <- dCDDM(test,drift,theta,tzero,boundary)
    discard.Angles <- test.D<0.0009
    space.C <- space.C[!discard.Angles]
    nBins.C <- length(space.C)
    
    space.RT <- seq(min.RT,max.RT,0.009)
    nBins.RT <- length(space.RT)
    
    cells <- expand.grid(space.C,space.RT)
    probs <-  numInt.tpz.cddm(cells[,1],cells[,2], par, plot=FALSE,
                              nBins.RT = nBins.RT, nBins.C = nBins.C)
  
    if(plot){
      scatterplot3d(cells[,1],cells[,2],probs, pch=16, zlim = c(0,1),
                    cex.symbols = 0.1, xlab = "", ylab = "", zlab = "",
                    xlim=range(cells[,1]), ylim=range(cells[,2]))
      mtext("Choices", side=1, line=0.5)
      mtext("RTs", side=4, line=-1, adj=0)
    }
  
    u <- runif(n,0,1)
    data <- matrix(NA, nrow=n, ncol=2)
    for(i in 1:n){
        better.match <- min(abs(probs-u[i]))
        found.at <- sample(which(abs(probs-u[i])==better.match),1)
        data[i,] <- as.numeric(as.vector(cells[found.at,]))
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