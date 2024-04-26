###############################################################################
###############################################################################
#####       Customary functions to compute, plot and compare CDFs (eCDF)
###############################################################################
########################################################   by Adriana F. Chavez   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 1: Code to compute the empirical CDF from different types of data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myECDF <- function(data){
  # Check if data is contained in a vector
  if(is.vector(data)){
        repeated.values <- !length(data)==length(unique(data))
        if(repeated.values){
          data <- unique(data)
        }
        n <- length(data)
        cdf <- NA
        for(i in 1:n){
          cdf[i] <- mean(data[i] > data)
        }
  # Otherwise, check if data is bivariate
  }else{if(ncol(data)==2){
        data.matrix.NxD <- data
        no.Dim <- ncol(data.matrix.NxD)
        n <- nrow(data.matrix.NxD)
        count <- matrix(NA,nrow=n,ncol=no.Dim)
        cdf <- rep(NA,n)
        for(i in 1:nrow(data.matrix.NxD)){
          X <- rep(TRUE, n)
          for(j in 1:no.Dim){
            count[,j] <- data[i,j] > data[,j]
          }
          for(j in 2:no.Dim){
            Y <- count[,j-1] & count[,j]
            X <- X & Y
          }
          cdf[i] <- mean(X)
        }
  }}
  return(cdf)
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 2: Plot the empirical CDF
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myECDF.Plot <- function(data, color="forestgreen"){
  eCDF.color <- color
  if(is.vector(data)){
    sorted.data <- sort(data)
    lab.x <- "Data"
    eCDF <- myECDF(sorted.data)
    plot(sorted.data, eCDF, pch=16, cex=0.7,
         col=eCDF.color, lwd=2, ylim=c(0,1), xlab=lab.x)
    mtext("eCDF",2, line=2.1, f=2)
  }else{
    if(ncol(data)==2){
      eCDF <- myECDF(data)
      lab.x <- colnames(data)[1]
      lab.y <- colnames(data)[2]
      if(is.null(lab.x)){ lab.x <-"Dimension 1" }
      if(is.null(lab.y)){ lab.y <-"Dimension 2" }
      scatterplot3d(data[,1],data[,2],eCDF, pch=16,
                    cex.symbols = 0.5, color = eCDF.color,
                    xlab = lab.x, ylab = lab.y, zlab = "eCDF")
    }else{
      print("Cannot plot more than two dimensions")
    }
  }
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 3: Compare two CDF (ideally, an eCDF against a theoretical CDF)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getDifferences <- function(eCDF,tCDF){
  difference <- tCDF - eCDF
  difference.sum <- sum(difference)
  Kldivergence <- max(abs(difference))
  sq.difference <- sum((difference)^2)
  
  output <- round(cbind(difference.sum,Kldivergence,sq.difference),4)
  colnames(output) <- c("sumDiff","Kldivergence","SSDiff") 
  return(output)
}


#################
# Test/Examples
#################
# Test function
if(!exists("test")){  test <- TRUE     }
if(test){
    # Data comes from Normal
    true.mean <- 10
    true.sd <- 1
    x <- rnorm(1000,true.mean,true.sd)
    x.eCDF <- myECDF(x)
    myECDF.Plot(x)
    x.tCDF <- pnorm(x,true.mean,true.sd)
    getDifferences(x.eCDF,x.tCDF)
    # Data comes from bivariate normal
    true.means <- c(10,10)
    true.sds <- diag(length(true.means))
    y <- mvtnorm::rmvnorm(1000,true.means,true.sds)
    y.eCDF <- myECDF(y)
    myECDF.Plot(y)
    y.tCDF <- pnorm(y,true.mean,true.sd)
    getDifferences(y.eCDF,y.tCDF)
}