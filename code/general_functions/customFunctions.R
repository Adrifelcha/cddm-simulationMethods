###############################################################################
###############################################################################
#####   A set of functions made for the simulation paper
###############################################################################
########################################################   by Adriana F. Ch?vez   



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
    assign.percentiles <- NA
    for(i in 1:n){
         cdf[i] <- mean(data.vector[i] > data.vector)
    }
    return(cdf)
}

n <- 10000
A <- rnorm(n,5,1)
B <- rnorm(n,5,1)
data.matrix.NxD <- cbind(A,B)

my_ecdf.MD <- function(data.matrix.NxD){
  data <- data.matrix.NxD
  no.Dim <- ncol(data.matrix.NxD)
  n <- nrow(data.matrix.NxD)
  count <- matrix(NA,nrow=n,ncol=no.Dim)
  X <- rep(TRUE, n)
  Z <- rep(NA,n)
  for(i in 1:nrow(data.matrix.NxD)){
      for(j in 1:no.Dim){
        count[,j] <- data[i,j] > data[,j]
      }
      for(j in 2:no.Dim){
        Y <- count[,j-1] & count[,j]
        X <- X & Y
      }
    Z[i] <- mean(X)
  }
  
  return(assign.percentiles)
}


