###############################################################################
###############################################################################
#####       Customary functions to compute, plot and compare CDFs (eCDF)
###############################################################################
########################################################   by Adriana F. Chavez   

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Part 1: Code to compute the empirical CDF from different types of data
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
myECDF <- function(data){    
    # Input validation
    if(!is.vector(data) && !is.matrix(data)) {
        stop("Input must be a vector or matrix")
    }
    if(is.matrix(data) && ncol(data) != 2) {
        stop("Matrix input must have exactly 2 columns")
    }
    
    # Case 1: Univariate data (vector input; P(X > x))
    if(is.vector(data)){
        n <- length(data)        
        cdf <- colMeans(outer(data, data, ">"))        
    # Case 2: Bivariate data (2-column matrix input; P(X1 > x1 AND X2 > x2))
    } else {
        n <- nrow(data)
        cdf <- numeric(n)  # Pre-allocate vector
        
        # For each observation, compute joint probability
        for(i in 1:n){
            # Vectorized comparison across both dimensions
            greater_than <- (data[i,1] > data[,1]) & (data[i,2] > data[,2])
            cdf[i] <- mean(greater_than)
        }
    }
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

#################
# Test/Examples
#################
# Test function
#if(!exists("test")){  test <- TRUE     }
#if(test){
    # Data comes from Normal
#    true.mean <- 10
#    true.sd <- 1
#    x <- rnorm(1000,true.mean,true.sd)
#    x.eCDF <- myECDF(x)
#    myECDF.Plot(x)
#    x.tCDF <- pnorm(x,true.mean,true.sd)
#    getDifferences(x.eCDF,x.tCDF)
    # Data comes from bivariate normal
#    true.means <- c(10,10)
#    true.sds <- diag(length(true.means))
#    y <- mvtnorm::rmvnorm(1000,true.means,true.sds)
#    y.eCDF <- myECDF(y)
#    myECDF.Plot(y)
    #y.tCDF <- mvtnorm::pmvnorm(lower = c(-Inf, -Inf), upper = y, mean = true.means, sigma = true.sds)
    #getDifferences(y.eCDF,y.tCDF)
#}