###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####   INVERSE PROBABILITY TRANSFORM ALGORITHM with a grid approximation
###############################################################################
########################################################   by Adriana F. Chavez 
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){ source("./pCDDM.R") }
library(scatterplot3d) 

sample.invCDF.cddm <- function(n, par, max.RT=NA, plot=FALSE){
  # Set-up - Load important values and define important variables
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Load parameter values
    drift <- par$drift;   theta <- par$theta
    tzero <- par$tzero;   boundary <- par$boundary
    # Derive a sensical max.RT point, if it hasn't been provided
    if(is.na(max.RT)){  max.RT <- (boundary/drift)*4    }
  # 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    # Define the Choice and RT values for which we'll c
    nC <- 300
    space.C  <- seq(0,2*pi,length.out=nC)
    nRT <- define_bins(tzero, max.RT)-1
    space.RT <- seq(tzero,max.RT,length.out=nRT+1)
    space.RT <- space.RT[2:(nRT+1)]
    
    cells <- expand.grid(space.C,space.RT)
    start_time <- Sys.time()
    probs <-  pCDDM(cells,drift,theta,tzero,boundary)
    end_time <- Sys.time()
    probs_tabloid <- matrix(probs,nrow=nC,ncol=nRT,byrow=FALSE)
    end_time-start_time
  # Set-up
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~  
    if(plot){
      #show <- seq(1,nrow(cells),length.out = 1000)
      #scatterplot3d(cells[show,1],cells[show,2],probs[show], pch=3,
      #              cex.symbols = 0.3, #color = rgb(0.4), 
      #              xlab = "", ylab = "", zlab = "")
      scatterplot3d(cells[,1],cells[,2],probs, pch=16, zlim = c(0,1),
                    cex.symbols = 0.05, xlab = "", ylab = "", zlab = "",
                    xlim=range(cells[,1]), ylim=range(cells[,2]))
      #mtext("inverse-CDF algorithm", f=2)
      mtext("Choices", side=1, line=0.5)
      mtext("RTs", side=4, line=-1, adj=0)
      #mtext("inverse CDFs", outer=TRUE, line=-2, f=2, cex=1.2)
    }
  # Set-up
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    u <- runif(n,0,1)
    data <- matrix(NA, nrow=n, ncol=2)
    for(i in 1:n){
        search.for <- u[i]
        better.match <- min(abs(probs_tabloid-search.for))
        found.at <- which(abs(probs_tabloid-search.for)==better.match, arr.ind = TRUE)
        data[i,1] <- space.C[found.at[1]]
        data[i,2] <- space.RT[found.at[2]]
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