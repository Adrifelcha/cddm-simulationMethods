###############################################################################
###############################################################################
#####   A set of functions made for the simulation paper
###############################################################################
########################################################   by Adriana F. Ch?vez   
library("scatterplot3d") 

###############################################################################
###################   General Simulation Functions
###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Take full random walk coordinates and extract final states
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
getFinalState <- function(randomWalk.states){
  randomWalk <- randomWalk.states
  dimensions <- dim(randomWalk)
  K <- nrow(randomWalk)
  
  if(length(dimensions)>2){
      I <- dimensions[3]
      coord <- matrix(NA, ncol=2,nrow=I)
      for(i in 1:I){
          for(k in 1:K){
              if(!is.na(randomWalk[k,1,i])){
                  a <- k
              }else{
                  break
              }
          }
          coord[i,] <- randomWalk[a,,i]
      }
      output <- coord
  }else{
      I <- dimensions[2]
      choice <- rep(NA, I)
      for(i in 1:I){
        for(k in 1:K){
          if(!is.na(randomWalk[k,i])){
            a <- k
          }else{
            break
          }
        }
        choice[i] <- randomWalk[a,i]
      }
      output <- choice
  }
  return(output)
}


###############################################################################
###################   General Functions related to the Accuracy tests
###############################################################################
my_ecdf.1D <- function(data.vector){
    repeated.values <- !length(data.vector)==length(unique(data.vector))
    if(repeated.values){
       data.vector <- unique(data.vector)
    }
    n <- length(data.vector)
    cdf <- NA
    for(i in 1:n){
         cdf[i] <- mean(data.vector[i] > data.vector)
    }
    return(cdf)
}

my_ecdf.MD <- function(data.matrix.NxD){
  data <- data.matrix.NxD
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
  return(cdf)
}

my_ecdf.Plot <- function(data, color="forestgreen"){
  eCDF.color <- color
  if(is.vector(data)){
    sorted.data <- sort(data)
    eCDF <- my_ecdf.1D(sorted.data)
    plot(sorted.data, eCDF, pch=16, cex=0.7,
         col=eCDF.color, lwd=2, ylim=c(0,1))
  }else{
    if(ncol(data)==2){
      eCDF <- my_ecdf.MD(data)
      scatterplot3d(data[,1],data[,2],eCDF, pch=16,cex.symbols = 0.5)
    }else{
      print("Cannot plot more than two dimensions")
    }
  }
}

getDifferences <- function(eCDF,tCDF){
  difference <- tCDF - eCDF
  difference.sum <- sum(difference)
  Kldivergence <- max(abs(difference))
  sq.difference <- sum((difference)^2)
  
  output <- cbind(difference.sum,Kldivergence,sq.difference)
  colnames(output) <- c("sumDiff","Kldivergence","SSDiff") 
}