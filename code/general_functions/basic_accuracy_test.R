###############################################################################
###############################################################################
#####   A set of simplified functions to illustrate how Accuracy Test is 
#####     conducted in our study, for Normal and Multivariate Normal data
###############################################################################
########################################################   by Adriana F. Chavez   
############## Load custom functions
source("./customFunctions.R")
library(mnormt)
library("scatterplot3d") 

accuracyTest.normal <- function(data.vector, t.mean=0, t.sd=1){
  sorted.data <- sort(data.vector)
  tCDF <- pnorm(sorted.data,t.mean,t.sd)
  tCDF.color <- "darkblue"
  eCDF <- my_ecdf.1D(sorted.data)
  eCDF.color <- "forestgreen"
  
  plot(sorted.data, tCDF, type="l", 
       col=tCDF.color, lwd=3, ylim=c(0,1))
  points(sorted.data,eCDF, cex=0.2, col=eCDF.color)
  legend("topleft", c("Theoretical CDF", "Empirical CDF"),
         col=c(tCDF.color,eCDF.color), lwd=2, cex=0.65,
         lty=c(1,3))
  
  output <- getDifferences(eCDF,tCDF)
  return(output)
}

accuracyTest.mvnormal <- function(data, t.mean=c(0,0), 
                                  t.Sigma=matrix(c(2,-1,1,2),nrow=2,byrow=TRUE)){
  tCDF <- pmnorm(data,t.mean,t.Sigma)
  tCDF.color <- "darkblue"
  eCDF <- my_ecdf.MD(data)
  eCDF.color <- "forestgreen"
  
  a <- scatterplot3d(data[,1],data[,2],eCDF, pch=16,cex.symbols = 0.5,
                     color = tCDF.color)
  a$points3d(data[,1],data[,2],tCDF, pch=16, col=eCDF.color, cex=0.2)
  legend("topleft", c("Theoretical CDF", "Empirical CDF"),
         col=c(tCDF.color,eCDF.color), lwd=2, cex=0.65,
         lty=c(1,3))
  
  output <- getDifferences(eCDF,tCDF)
  return(output)
}