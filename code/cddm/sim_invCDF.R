###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####   INVERSE PROBABILITY TRANSFORM ALGORITHM with a grid approximation
###############################################################################
########################################################   by Adriana F. Chavez 
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){ source("./pCDDM.R") }
library(scatterplot3d) 

par <- list("drift" = 1, 
            "theta" = pi,
            "tzero" = 0.1,
            "boundary" = 7)
# par <- list("drift" = 3.5, 
#             "theta" = pi,
#             "tzero" = 0.2,
#             "boundary" = 2)

n <- 200
max.RT <- 20

sample.invCDF.cddm <- function(n, par, max.RT=20, plot=FALSE){
  
    drift <- par$drift;   theta <- par$theta
    tzero <- par$tzero;   boundary <- par$boundary
    
    n.space <- 300
    space.C  <- seq(0,2*pi,length.out=n.space)
    space.RT <- seq(tzero,max.RT,length.out=n.space)
    
    cells <- expand.grid(space.C,space.RT)
    start_time <- Sys.time()
    probs <-  pCDDM(cells,drift,theta,tzero,boundary)
    end_time <- Sys.time()
    probs_tabloid <- matrix(probs,nrow=n.space,ncol=n.space,byrow=FALSE)
    end_time-start_time
    
    if(plot){
      scatterplot3d(cells[,1],cells[,2],probs, pch=3,
                    cex.symbols = 0.3, color = col.Reject(0.4), 
                    xlab = "", ylab = "", zlab = "")
      #mtext("inverse-CDF algorithm", f=2)
      mtext("Choices", side=1, line=0.5)
      mtext("RTs", side=4, line=-1, adj=0)
      #mtext("inverse CDFs", outer=TRUE, line=-2, f=2, cex=1.2)
    }
    
    u <- runif(n,0,1)
    data <- matrix(NA, nrow=n, ncol=2)
    for(i in 1:n){
        search.for <- u[i]
        better.match <- min(abs(probs_tabloid-search.for))
        found.at <- which(abs(probs_tabloid-search.for)==better.match, arr.ind = TRUE)
        data[i,1] <- space.C[found.at[1]]
        data[i,2] <- space.C[found.at[2]]
    }
    colnames(data) <- c("Choice", "RT")
return(data)
}


# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
    par <- list("drift" = 1, 
                "theta" = pi,
                "tzero" = 0.1,
                "boundary" = 7)
    # par <- list("drift" = 3.5, 
    #             "theta" = pi,
    #             "tzero" = 0.2,
    #             "boundary" = 2)
    
    n <- 20
    max.RT <- 20
    sample.invCDF.cddm(n,par,max.RT = 20)
}