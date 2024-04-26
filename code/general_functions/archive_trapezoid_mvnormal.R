###############################################################################
#####   2.1  Bivariate Normal data: CDF approximation
###############################################################################
library(mnormt)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 2.1.a Auxiliary function: 
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
par <- list("mean" = c(10,10),
            "Sigma" = matrix(c(1, 0.5, 0.5, 1), nrow=2))
lower.bound.X <- -3
upper.bound.X <- 10
lower.bound.Y <- 1
upper.bound.Y <- 10
width.X <- upper.bound.X-lower.bound.X
width.Y <- upper.bound.Y-lower.bound.Y
#kappa <- max(width.X,width.Y)*10
kappa <- c(width.X,width.Y)*10
Mean  <- par$mean
Sigma <- par$Sigma

bin.X <- seq(lower.bound.X,upper.bound.X,length.out=kappa[1])
bin.Y <- seq(lower.bound.Y,upper.bound.Y,length.out=kappa[2])

testBins.bvnormal = function(bin.X,bin.Y,Mean,Sigma){
  
  n = min(length(bin.X),length(bin.Y))
  bin.vertices <- cbind(bin.X[1:n],bin.Y[1:n])
  d <- dmnorm(bin.vertices,Mean,Sigma)
  b <- which(d>0)[1]
  if(!is.na(b)){
    if(b==1){
      start.at <- NA
      for(i in 1:2){
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

numInt.tpz.bvnormal <- function(lower.bound.X, upper.bound.X,
                                lower.bound.Y, upper.bound.Y,
                                par, plot=FALSE){
  no.Dim <- 2
  Mean  <- par$mean
  Sigma <- par$Sigma
  
  width.X <- upper.bound.X-lower.bound.X
  width.Y <- upper.bound.Y-lower.bound.Y
  kappa <- c(width.X,width.Y)*10
  
  if(plot){
    width <- c(Sigma[1,1],Sigma[2,2])*5
    nSupp <- 100
    nLines <- 30
    base.X <- c(Mean[1]-width[1],Mean[1]+width[1])
    base.Y <- c(Mean[2]-width[2],Mean[2]+width[2])
    support.X <- seq(base.X[1],base.X[2],length.out=nSupp)
    support.Y <- seq(base.Y[1],base.Y[2],length.out=nSupp)
    support.muX <- rep(Mean[1],nSupp)
    support.muY <- rep(Mean[2],nSupp)
    z.diag1 <- dmnorm(cbind(support.X,support.Y),Mean,Sigma)
    z.diag2 <- dmnorm(cbind(support.X,rev(support.Y)),Mean,Sigma)
    z.Y_at_muX <- dmnorm(cbind(support.muX,support.Y),Mean,Sigma)
    z.X_at_muY <- dmnorm(cbind(support.X,support.muY),Mean,Sigma)
    z.top <- max(c(z.diag1,z.diag2,z.Y_at_muX,z.X_at_muY))
    a <- scatterplot3d(support.X, support.Y, z.diag1, 
                       xlim=base.X, ylim=base.Y, zlim=c(0,z.top),
                       color="blue", type="l", zlab="")
    a$points3d(support.X, rev(support.Y), z.diag2, col = "blue", type="l")
    a$points3d(support.X, support.muY, z.X_at_muY, col = "red", type="l")
    a$points3d(support.muX, support.Y, z.Y_at_muX, col = "red", type="l")
    L <- round(nSupp/nLines,0)
    for(i in 1:nLines){
      choose.X <- rep(support.X[i*L],nSupp)
      choose.Y <- rep(support.Y[i*L],nSupp)
      z.overX <- dmnorm(cbind(support.X,choose.Y),Mean,Sigma)
      z.overY <- dmnorm(cbind(choose.X,support.Y),Mean,Sigma)
      a$points3d(support.X, choose.Y, z.overX, col = "blue", type="l")
      a$points3d(choose.X, support.Y, z.overY, col = "blue", type="l")
    }
  }
  
  bin.X <- seq(lower.bound.X,upper.bound.X,length.out=kappa[1])
  bin.Y <- seq(lower.bound.Y,upper.bound.Y,length.out=kappa[2])
  
  test.bins = testBins.bvnormal(bin.X,bin.Y,Mean,Sigma)
  bin.X <- test.bins$bin.X
  bin.Y <- test.bins$bin.Y
  
  if(!unique(is.na(bin.X))){
    test.bins = testBins.bvnormal(rev(bin.X),rev(bin.Y),Mean,Sigma)
    bin.X <- sort(test.bins$bin.X)
    bin.Y <- sort(test.bins$bin.Y)
    b <- 2
  }
  nX <- length(bin.X)
  nY <- length(bin.Y)
  y <- 1
  y.up <- nY-1
  side.X <- bin.X[2]-bin.X[1]
  side.Y <- bin.Y[2]-bin.Y[1]
  base.area <- side.X*side.Y
  K <- length(bin.X)*length(bin.Y)
  bin.area <- rep(NA,K)
  if(unique(is.na(bin.X))){ b <- K+1 }
  
  X.megavector <- rep(bin.X,each=nY)
  Y.megavector <- rep(bin.Y, nX)
  d.megavector <- dmnorm(cbind(X.megavector,Y.megavector),Mean,Sigma)
  
  total.area = 0
  p <- c(1,rep(2:(nY-1),each=2),nY)
  while((total.area<1&b<=nX)==TRUE){
    D.from <- matrix(d.from[p],ncol=2,byrow=TRUE)
    D.to <- matrix(d.to[p],ncol=2,byrow=TRUE)
    D <- cbind(D.from,D.to)
    height <- apply(D,1,mean)
    bin.area[y:y.up] <- height*base.area
    y <- y.up+y
    y.up <- y.up+y.up
    total.area <- sum(bin.area,na.rm = TRUE)
    # if(plot){
    #   polygon(x=c(low.x,up.x,up.x,low.x),
    #           y=c(0,0,low.y,up.y),
    #           col = (b %% 2)+1)
    # }
    b <- b+1
  }
  total.area <- sum(bin.area)
  return(total.area)
}

#~~~~~~~~~~~~~~~#
# Test/Examples #
#~~~~~~~~~~~~~~~#
if(!exists("test")){  test <- TRUE     }
if(test){
  par <- list("mean" = c(10,10),
              "Sigma" = matrix(c(1, 0.5, 0.5, 1), nrow=2))
  lower.bound.X <- -3
  upper.bound.X <- 10
  lower.bound.Y <- 1
  upper.bound.Y <- 10
  plot = TRUE
  numInt.tpz.bvnormal(lower.bound.X, upper.bound.X, lower.bound.Y, upper.bound.Y, par)
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.1.b  User friendly cdf funciton
#        Compute the approximate normal cdf for a single value of x
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
mvnormal.cdf <- function(x,par){
  lower.bound <- par$mean-(par$sd*100)
  area <- numInt.tpz.bvnormal(lower.bound,x,par)
  return(area)
}


# Test function
mvnormal.cdf(10,par=par)