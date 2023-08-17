###############################################################################
###############################################################################
#####   A set of functions made to plot the data generated from simulations
###############################################################################
########################################################   by Adriana F. Chavez   
############## Load custom functions
source("../general_functions/customFunctions.R")

data.vector <- rexp(1000, 5)

my_ecdf <- function(data.vector){
  repeated.values <- !length(data.vector)==length(unique(data.vector))
  if(repeated.values){
    data.vector <- unique(data.vector)
  }
  
  n <- length(data.vector)
  percentiles <- seq(0,1,length.out=n)
  sorted.data <- order(data.vector)
  assign.percentiles <- percentiles[sorted.data]
  return(assign.percentiles)
}

accuracy.test <- function(data.vector){
  sorted.data <- sort(data.vector)
  tCDF <- pnorm(sorted.data,10,1)
  tCDF.color <- "darkblue"
  eCDF <- my_ecdf(sorted.data)
  eCDF.color <- "forestgreen"
  
  plot(sorted.data, tCDF, type="l", 
       col=tCDF.color, lwd=2, ylim=c(0,1))
  points(sorted.data,eCDF, cex=0.2, col=eCDF.color)
  legend("topleft", c("Theoretical CDF", "Empirical CDF"),
         col=c(tCDF.color,eCDF.color), lwd=2, cex=0.65,
         lty=c(1,3))
  
  difference <- tCDF - eCDF
  difference.sum <- sum(difference)
  abs.difference <- sum(abs(difference))
  sq.difference <- sum((difference)^2)
  
  output <- cbind(difference.sum,abs.difference,sq.difference)
  colnames(output) <- c("sumDiff","absDiff","SSDiff")
  return(output)
}
