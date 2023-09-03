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
par <- list("mean" = 5000, "sd" = 1)
lower.bound <- 1
upper.bound <- 10
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
      ticks <- seq(mean-width,mean+width,length.out=10)
      support <- seq(mean-width,mean+width,length.out=9999)
      plot(support,dnorm(support,mean,sd),
           type="l", ann=F, axes=F, col=1)
      axis(1,ticks,ticks)
    }
    
    bins <- seq(lower.bound,upper.bound,length.out=kappa)
    d <- dnorm(bins,par$mean,par$sd)
    b <- which(d>0)[1]
    if(is.na(b)){
      bin.area <- rep(NA,kappa)
      b <- kappa+1
    }else{
      bin.area <- rep(NA,kappa-(b-1))}
    total.area <- 0
    while((total.area<1&b<=kappa)==TRUE){
        low.x <- bins[b-1]
        up.x  <- bins[b]
        low.y <- dnorm(low.x,par$mean,par$sd)
        up.y  <- dnorm(up.x,par$mean,par$sd)
        height <- (low.y+up.y)/2
        bin.area[b-1] <- height*(up.x-low.x)
        if(plot){
          polygon(x=c(low.x,up.x,up.x,low.x),
                  y=c(0,0,low.y,up.y),
                  col = (b %% 2)+1)
        }
        total.area <- sum(bin.area,na.rm=TRUE)
        b <- b+1
    }
    b <- b-1
    legend("topright", paste("No. bins =",b), 
           cex = 0.75, bty ="n")
    return(total.area)
}

# Test function
numInt.tpz.normal(lower.bound,upper.bound,par, plot=TRUE)

# Use Trapezoid Numeric integration to compute CDF
normal.cdf <- function(x,par,plot=FALSE){
  lower.bound <- x-(par$mean-(par$sd*100))
  area <- numInt.tpz.normal(lower.bound,x,par,plot=plot)
  return(area)
}

# Test function
normal.cdf(10,par,plot=TRUE)
normal.cdf(50,par,plot=TRUE)
normal.cdf(500,par,plot=TRUE)
normal.cdf(50000,par,plot=TRUE)

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
            "Sigma" = matrix(c(1, 1, 1, 1), nrow=2))
lower.bound.X <- 1
upper.bound.X <- 10
lower.bound.Y <- 1
upper.bound.Y <- 10
plot = TRUE


# Write Trapezoid N.I. algorithm
numInt.tpz.bvnormal <- function(lower.bound.X, upper.bound.X,
                              lower.bound.Y, upper.bound.Y,
                              par, plot=FALSE){
  Mean  <- par$mean
  Sigma <- par$Sigma
  
  width.X <- upper.bound.X-lower.bound.X
  kappa.X <- width.X*20
  width.Y <- upper.bound.Y-lower.bound.Y
  kappa.Y <- width.Y*20
  
  if(plot){
    width <- NA
    ticks <- matrix(NA,nrow=2,ncol=10)
    support.def <- 70
    support <- matrix(NA,nrow=2,ncol=support.def)
    for(i in 1:2){
        width[i] <- (3*Sigma[i,i])
        ticks[i,] <- seq(Mean[i]-width[i],
                         Mean[i]+width[i],length.out=10)
        support[i,] <- seq(Mean[i]-width[i],
                           Mean[i]+width[i],length.out=support.def)
    }
    z <- outer(support[1,],support[2,], f)
    persp(support[1,], support[2,], z, col="gray99",
          theta=-30, phi=25, expand=0.6, ticktype='detailed')
    legend("topright", paste("No. bins =",kappa), 
           cex = 0.75, bty ="n")
  }
  
  bin.X <- seq(lower.X,upper.X,length.out=kappa)
  bin.Y <- seq(lower.Y,upper.Y,length.out=kappa)
  bin.area <- rep(NA,kappa-1)
  for(b in 2:kappa){
    X.from <- bin.X[b-1]
    X.to   <- bin.X[b]
    Y.from <- bin.Y[b-1]
    Y.to   <- bin.Y[b]
    side.A <- X.from-X.to
    side.B <- Y.from-Y.to
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
numInt.tpz.normal(lower.bound,upper.bound,par, plot=TRUE)

# Use Trapezoid Numeric integration to compute CDF
mvnormal.cdf <- function(x,par){
  lower.bound <- par$mean-(par$sd*100)
  area <- numInt.tpz.normal(lower.bound,x,par)
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
