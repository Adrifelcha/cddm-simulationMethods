###############################################################################
###############################################################################
#####   
###############################################################################
########################################################   by Adriana F. Chavez   
############## Load custom functions
source("../general_functions/customFunctions.R")

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
  Kldivergence <- max(abs(difference))
  sq.difference <- sum((difference)^2)
  
  output <- cbind(difference.sum,abs.difference,sq.difference)
  colnames(output) <- c("sumDiff","Kldivergence","SSDiff")
  return(output)
}
