###############################################################################
###############################################################################
#####   A set of functions made for the simulation paper
###############################################################################
########################################################   by Adriana F. Ch?vez   
library("scatterplot3d") 
library("mnormt")

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


###############################################################################
###################   Auxiliary functions for Sampling algorithms
###############################################################################
testBins.bvnormal = function(bin.X,bin.Y,Mean,Sigma){
  dimm = 1:2
  n = min(length(bin.X),length(bin.Y))
  bin.vertices <- cbind(bin.X[1:n],bin.Y[1:n])
  d <- dmnorm(bin.vertices,Mean,Sigma)
  b <- which(d>0)[1]
  if(!is.na(b)){
      if(!b!=1){
          start.at <- NA
          for(i in dimm){
          fix.this.bin <- bin.vertices[,i]
          let.bin.vary <- bin.vertices[,-i]
          unique.densities <- 2
          a <- b
              while(unique.densities>1|a>0){
                fix = rep(fix.this.bin[a],a)
                bin.vertices <- matrix(NA,ncol=2,nrow=a)
                bin.vertices[,i] <- fix
                bin.vertices[,-i] <- let.bin.vary[a:1]
                d.try = dmnorm(bin.vertices,Mean,Sigma)
                unique.densities = length(unique(d.try))
                a = a-1 
              }
          start.at[i] <- a+1
          }
          bin.X <- bin.X[start.at[1]:n]
          bin.Y <- bin.Y[start.at[2]:n]
     }
  }else{
    bin.X <- NA
    bin.Y <- NA
  }
  output = list("bin.X" = bin.X, "bin.Y" = bin.Y)
  return(output)
}