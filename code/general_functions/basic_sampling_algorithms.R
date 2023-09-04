###############################################################################
###############################################################################
#####   A set of simplified NORMAL EXAMPLES for the sampling algorithms 
#####     we use as...
#####     1) Numeric integration of the cdf of both the DDM and CDDM models
#####     2) Rejection algorithms to sample observations for DDM and CDDM
###############################################################################
########################################################   by Adriana F. Chavez   
############## Load custom functions
source("./customFunctions.R")
library(mnormt)
library(scatterplot3d) 

###############################################################################
#####   Unidimensional Normal data
###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To approximate CDFs, we use a Trapezoid Numerical Integration algorithm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some data
par <- list("mean" = 500, "sd" = 1)
lower.bound <- 1
upper.bound <- 5000
kappa=NA
plot=TRUE
# Write Trapezoid N.I. algorithm
numInt.tpz.normal <- function(lower.bound, upper.bound, par,
                              kappa=NA, plot=FALSE){
    mean <- par$mean
    sd <- par$sd
    if(is.na(kappa)){
          total.width <- upper.bound-lower.bound
          kappa <- total.width*20
    }
    
    if(plot){
      width <- (3*sd)
      ticks <- round(seq(mean-width,mean+width,length.out=10),2)
      support <- seq(mean-width,mean+width,length.out=9999)
      plot(support,dnorm(support,mean,sd),
           type="l", ann=F, axes=F, col=1)
      axis(1,ticks,ticks)
    }
    
    bins <- seq(lower.bound,upper.bound,length.out=kappa)
    bin.base = bins[2] - bins[1]
    d <- dnorm(bins,par$mean,par$sd)
    b <- which(d>0)[1]
    if(is.na(b)){
      b <- kappa+1
    }
    bin.area <- rep(NA,kappa)
    total.area <- 0
    while((total.area<1&b<=kappa)==TRUE){
        low.x = bins[b-1]
        up.x  = bins[b]
        density.x1 <- d[b-1]
        density.x2 <- d[b]
        height <- (density.x1+density.x2)/2
        bin.area[b-1] <- height*bin.base
        if(plot){
          polygon(x=c(low.x,up.x,up.x,low.x),
                  y=c(0,0,density.x1,density.x2),
                  col = (b %% 2)+1)
        }
        total.area <- sum(bin.area,na.rm=TRUE)
        b <- b+1
    }
    b <- b-1
    no.bins = sum(!is.na(bin.area))
    legend("topright", paste("No. bins =",no.bins), 
           cex = 0.75, bty ="n")
    return(total.area)
}

# Test function
numInt.tpz.normal(lower.bound,upper.bound,par, plot=TRUE)

# Use Trapezoid Numeric integration to compute CDF
normal.cdf <- function(x,par,plot=FALSE){
  lower.bound <- (par$mean-(par$sd*100))
  if(x<lower.bound){
   lower.bound  = x-lower.bound
  }
  area <- numInt.tpz.normal(lower.bound,x,par,plot=plot)
  return(area)
}

# Test function
normal.cdf(10,par,plot=TRUE)
normal.cdf(50,par,plot=TRUE)
normal.cdf(500,par,plot=TRUE)
normal.cdf(501,par,plot=TRUE)
normal.cdf(5000,par,plot=TRUE)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To sample data, we use a MCMC basic algorithm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some data
par <- list("mean" = 10, "sd" = 1)
n <- 1000
# Write Trapezoid N.I. algorithm
sample.MCMC.normal <- function(n, par, plot=FALSE){
  mean <- par$mean
  sd <- par$sd
  
  max.D <- dnorm(mean,mean,sd)
  height <- max.D*1.05
  width <- (3.5*sd)
  base.1 <- mean-width
  base.2 <- mean+width
  
  if(plot){
            plot.width <- (4*sd)
            x.low <- mean-plot.width
            x.top <- mean+plot.width
            ticks <- round(seq(x.low,x.top,length.out=10),2)
            support <- seq(x.low,x.top,length.out=9999)
            plot(support,dnorm(support,mean,sd), xlim=c(x.low,x.top),
                 type="l", ann=F, axes=F, col=1, ylim=c(0,height))
            axis(1,ticks,ticks)
            abline(h=height)
            abline(v=c(base.1,base.2))
  }
  
  n.keep <- 0
  n.try <- n
  samples <- NA
  while(n.keep < n){
        cand <- runif(n.try,base.1,base.2)
        eval <- dnorm(cand,mean,sd)
        rej.crit <- runif(n.try,0,height)  
        keep <- (eval >= rej.crit)
        
        n.keep <- sum(keep)
        n.try <- n.try-n.keep
        
        if(plot){
                  points(cand[!keep],rej.crit[!keep],col="red", pch=16, cex=0.3)
                  points(cand[keep],rej.crit[keep],col="blue3", pch=16, cex=0.3)
        }
        
        samples <- c(samples, cand[keep])
        n.keep <- length(samples)
  }
  
  return(samples)
}

# Test function
sample.MCMC.normal(5000,par, plot=TRUE)


###############################################################################
#####   Bivariate Normal data
#####   This section is merely illustrative, it requires R packages
###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To approximate CDFs, we use a Trapezoid Numerical Integration algorithm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some data
par <- list("mean" = c(10,10),
            "Sigma" = matrix(c(1, 0.5, 0.5, 1), nrow=2))
lower.bound.X <- 1
upper.bound.X <- 10
lower.bound.Y <- 1
upper.bound.Y <- 10
plot = TRUE


# Write Trapezoid N.I. algorithm
numInt.tpz.bvnormal <- function(lower.bound.X, upper.bound.X,
                                lower.bound.Y, upper.bound.Y,
                                par, plot=FALSE){
  no.Dim <- 2
  Mean  <- par$mean
  Sigma <- par$Sigma
  
  width.X <- upper.bound.X-lower.bound.X
  width.Y <- upper.bound.Y-lower.bound.Y
  kappa <- c(width.X,width.Y)*20
  
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
  
  side.X <- bin.X[2]-bin.X[1]
  side.Y <- bin.Y[2]-bin.Y[1]
  K <- length(bin.X)*length(bin.Y)
  bin.area <- rep(NA,K)
  if(unique(is.na(bin.X))){   
    b <- K+1}
  
  total.area = 0
  while((total.area<1&b<=K)==TRUE){
    X.from <- bin.X[b-1]
    X.to   <- bin.X[b]
    Y.from <- bin.Y[b-1]
    Y.to   <- bin.Y[b]
    low.y <- dnorm(low.x,par$mean,par$sd)
    up.y  <- dnorm(up.x,par$mean,par$sd)
    height <- (low.y+up.y)/2
    bin.area[b-1] <- height*(up.x-low.x)
    
    if(plot){
      polygon(x=c(low.x,up.x,up.x,low.x),
              y=c(0,0,low.y,up.y),
              col = (b %% 2)+1)
    }
  }
  total.area <- sum(bin.area)
  return(total.area)
}

# Test function
numInt.tpz.bvnormal(lower.bound,upper.bound,par, plot=TRUE)

# Use Trapezoid Numeric integration to compute CDF
mvnormal.cdf <- function(x,par){
  lower.bound <- par$mean-(par$sd*100)
  area <- numInt.tpz.bvnormal(lower.bound,x,par)
  return(area)
}

# Test function
mvnormal.cdf(10,par)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# To sample data, we use a MCMC basic algorithm
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Some data
par <- list("mean" = c(10,10), 
            "Sigma" = matrix(c(1,0.5,0.5,1.5),nrow=2,byrow=TRUE))
n <- 5000

# Write Trapezoid N.I. algorithm
sample.MCMC.mvnormal <- function(n, par, plot=FALSE){
      mean <- par$mean
      Sigma <- par$Sigma
      no.Dim <- length(mean)
      
      if(no.Dim>2 & plot==TRUE){
            plot <- FALSE
            print("Cannot make a plot with more than 2Dimensions")
      }
      
      max.Density <- dmnorm(mean,mean,Sigma)
      height <- max.Density*1.05
      
      if(plot){
          nSupp <- 100
          x.low <- NA
          x.top <- NA
          support <- matrix(NA, nrow=no.Dim, ncol=nSupp)
      }
      
      base.1 <- NA
      base.2 <- NA
      for(j in 1:no.Dim){
          sd <- sqrt(Sigma[j,j])
          width <- 3.5*sd
          base.1[j] <- mean[j]-width
          base.2[j] <- mean[j]+width
          if(plot){
              plot.width <- (4*sd)
              x.low[j] <- mean[j]-plot.width
              x.top[j] <- mean[j]+plot.width
              support[j,] <- seq(x.low[j],x.top[j],length.out=nSupp)
          }
      }
      
      if(plot){
        z <- dmnorm(t(support),mean,Sigma)
        a <- scatterplot3d(support[1,],support[2,],z, xlim=c(x.low[1],x.top[1]),
                           ylim=c(x.low[2],x.top[2]), color="blue",zlim=c(0,height),
                           type="l", lwd=2)
        a$points3d(c(base.1[1],base.2[1],base.2[1],base.1[1]), 
                   c(base.1[2],base.2[2],base.1[2],base.2[2]), rep(height,4),
                     col = "skyblue", type = "h", pch = 16)
        a$points3d(c(base.1[1],base.2[1],base.2[1],base.1[1],
                     base.1[1],base.2[1],base.2[1],base.1[1]), 
                   c(base.1[2],base.1[2],base.2[2],base.2[2],
                     base.1[2],base.1[2],base.2[2],base.2[2]), rep(height,8),
                   col = "skyblue", type = "l", pch = 16)
      }
    
      n.keep <- 0
      n.try <- n
      samples <- matrix(NA, nrow=1, ncol=no.Dim)
      
      while(n.keep < n){
        cand <- matrix(NA, nrow=n.try, ncol=no.Dim)
        for(j in 1:no.Dim){
            cand[,j] <- runif(n.try,base.1[j],base.2[j])
        }
        eval <- dmnorm(cand,mean,Sigma)
        rej.crit <- runif(n.try,0,height)  
        keep <- (eval >= rej.crit)
        
        n.keep <- sum(keep)
        n.try <- n.try-n.keep
        
         if(plot){
           a$points3d(cand[!keep,1], cand[!keep,2], rej.crit[!keep],
                      col = "red", pch = 16, cex = 0.2)
           a$points3d(cand[keep,1], cand[keep,2], rej.crit[keep],
                      col = "blue3", pch = 16, cex = 0.2)
         }
        
        samples <- rbind(samples, cand[keep,])
        n.keep <- nrow(samples)-1
      }
      
      samples <- samples[-1,]
      return(samples)
}

# Test function
sample.MCMC.mvnormal(1000,par, plot=TRUE)
