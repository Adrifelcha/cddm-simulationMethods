###############################################################################
###############################################################################
#####   Accuracy test for the CDDM:
#####   A comparison between the bivariate empirical CDF and theoretical CDF
###############################################################################
########################################################   by Adriana F. Chavez   
test <- FALSE # Turn-off built-in examples.
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){       
      source("../cddm/pCDDM.R")}

accuracyTest.cddm <- function(data, par){
      drift <- par$drift
      theta <- par$theta
      tzero <- par$tzero
      boundary <- par$boundary
      
      tCDF <- pCDDM(data,par)
      tCDF.color <- "darkblue"
      eCDF <- my_ecdf.MD(data)
      eCDF.color <- "forestgreen"
      
      a <- scatterplot3d(data[,1],data[,2],eCDF, pch=16,cex.symbols = 0.5,
                         color = eCDF.color)
      a$points3d(data[,1],data[,2],tCDF, pch=16, col=eCDF.color, cex=0.2)
      legend("topleft", c("Theoretical CDF", "Empirical CDF"),
             col=c(tCDF.color,eCDF.color), lwd=2, cex=0.65,
             lty=c(1,3))
      
      output <- getDifferences(eCDF,tCDF)
      return(output)
}