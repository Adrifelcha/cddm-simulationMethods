###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                     Metropolis Sampling Algorithm
###############################################################################
########################################################   by Adriana F. Ch?vez 
#if(!exists("superCalled")){superCalled <- FALSE}
#if(!superCalled){ 
#  source("./dCDDM.R") 
#  source("./sim_auxiliarFunctions.R")
#}
library(mvtnorm)

# Write a simple Rejection algorithm for the CDDM pdf
sample.Metropolis.cddm <- function(n, par, plot=FALSE){
  drift <- par$drift;     theta <- par$theta
  tzero <- par$tzero;     boundary <- par$boundary
  predChoice <- theta
  
  test.Density <- keyDensityPoints(par)
  min.RT <- test.Density$min.RT
  max.RT <- test.Density$max.RT
  max.Density <- test.Density$max.Density
  predRT <- test.Density$pred.RT
  
  if(plot){          
    par(mar = c(0, 0, 0, 0)) 
    # Base plot: Sketch the bivariate density curve
    nLines <- 40
    x.C <- seq(0,2*pi,length.out=nLines)            # X: Choices
    y.RT <- seq(min.RT,max.RT,length.out=nLines)    # Y: RT
    x.theta <- rep(theta,nLines)                    # Make sure to add theta
    # Compute the density at each intersection point
    z.Dens <- matrix(NA, nrow=nLines, ncol=nLines)
    for(c in 1:nLines){ for(t in 1:nLines){
      z.Dens[c,t] <- dCDDM(c(x.C[c],y.RT[t]),drift,theta,tzero,boundary)
    }}
    # Compute the density at Choice = theta
    theta.Dens <- dCDDM(cbind(x.theta,y.RT),drift,theta,tzero,boundary)
    # Draw bivariate density curve
    a <- scatterplot3d(x.C, y.RT, diag(z.Dens), zlab="Density", color = "white", type="l",
                       xlab="Choices", ylab="RT", zlim = c(0, max(z.Dens, theta.Dens)))
  }
  
  M <- 300
  ARate_des <- 0.4
  ARate_obs <- ARate_des
  Mu <- c(predChoice, predRT)
  Sigma <- diag(c((2*pi)^2,50))
  warm_up <- TRUE
  nLines2 <- 15
  x.C2 <- seq(0,2*pi,length.out=nLines2)            # X: Choices
  y.RT2 <- seq(min.RT,max.RT,length.out=nLines2)    # Y: RT
  
  while(warm_up){
    b <- runif(1,0,1);   e <- runif(1,0,1);   d <- runif(1,0,1)
    baseColor <- rgb(b,e,d,0.05);   meanColor <- rgb(b,e,d,0.4)
    newSigma <- Sigma*((ARate_obs/ARate_des)^(0.5))
    Sigma <- newSigma
    max.DensityMVN <- dmvnorm(Mu,Mu,Sigma)
    if(plot){
      z.Dens2 <- matrix(NA, nrow=nLines2, ncol=nLines2)
      for(c in 1:nLines2){ for(t in 1:nLines2){
        z.Dens2[c,t] <- dmvnorm(c(x.C2[c],y.RT2[t]),mean=Mu, sigma = Sigma)
      }}
      z.Dens2 <- (z.Dens2*max.Density)/max.DensityMVN
      for(i in 1:nLines2){  
        a$points3d(rep(x.C2[i],nLines2), y.RT2, z.Dens2[i,], type="l", col = baseColor)
        a$points3d(x.C2, rep(y.RT2[i],nLines2), z.Dens2[,i],  col = baseColor, type="l")
      }
      a$points3d(Mu[1],Mu[2],max.Density, col = meanColor, cex=0.5, pch=16)
    }
    u <- runif(M,0,1)
    change <- rep(NA,M)
    for(i in 1:M){
      cand <- rmvnorm(1,Mu,Sigma)
      cand[1] <- cand[1] %% (2*pi)
      ratio.num <- max(dCDDM(cand,drift, theta, tzero, boundary),0,na.rm = TRUE)
      ratio.den <- dCDDM(Mu,drift, theta, tzero, boundary)
      ratio <- ratio.num/ratio.den
      change[i] <- ratio>u[i]
      if(change[i]){
        Mu <- as.vector(cand)
        if(plot){          
          z.Dens2 <- matrix(NA, nrow=nLines2, ncol=nLines2)
          for(c in 1:nLines2){ for(t in 1:nLines2){
            z.Dens2[c,t] <- dmvnorm(c(x.C2[c],y.RT2[t]),mean=Mu, sigma = Sigma)
          }}
          z.Dens2 <- (z.Dens2*max.Density)/max.DensityMVN
          for(i in 1:nLines2){  
            a$points3d(rep(x.C2[i],nLines2), y.RT2, z.Dens2[i,], type="l", col = baseColor)
            a$points3d(x.C2, rep(y.RT2[i],nLines2), z.Dens2[,i],  col = baseColor, type="l")
          }
          a$points3d(Mu[1],Mu[2],max.Density, col = meanColor, cex=0.5, pch=16)
        }
      }
    }
    ARate_obs <- mean(change)
    warm_up <- ARate_obs<ARate_des
    if(ARate_obs==0){ ARate_obs <- ARate_des}
  }
  
  if(plot){
    baseColorCDDM <- rgb(0,0,0,0.5)
    for(i in 1:nLines){  
      a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l", col = baseColorCDDM)
      a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i],  col = baseColorCDDM, type="l")
    }
    a$points3d(x.theta, y.RT, theta.Dens, col = baseColorCDDM, type="l")
  }
  
  samples <- matrix(NA, nrow=n, ncol=2)
  u <- runif(n,0,1)
  for(i in 1:n){
    cand <- c(0,-1)
    while(cand[2]<=min.RT){  cand <- rmvnorm(1,Mu,Sigma)   }
    samples[i,1] <- cand[1] %% (2*pi)
    samples[i,2] <- cand[2]
    ratio.num <- max(dCDDM(cand,drift, theta, tzero, boundary),0,na.rm = TRUE)
    ratio.den <- dCDDM(Mu,drift, theta, tzero, boundary)
    ratio <- ratio.num/ratio.den
    if(ratio>u[i]){
      Mu <- as.vector(cand)
    }
  }
  
  if(plot){
    a$points3d(samples[,1],samples[,2],dCDDM(samples,drift,theta,tzero,boundary),
               col = rgb(b,e,d,0.5), pch = 16, cex = 0.5)
  }
  colnames(samples) <- c("Choice","RT")
  return(samples)
}

# Test function
#if(!exists("test")){    test <- TRUE                           }
#if(test){
#  par <- list("drift" = 1, 
#              "theta" = pi,
#              "tzero" = 0.1,
#              "boundary" = 7)
#  n <- 5000
#  sample.Reject.cddm(1000,par, plot=TRUE)  }