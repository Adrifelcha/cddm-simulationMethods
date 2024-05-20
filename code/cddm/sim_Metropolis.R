###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                     Metropolis Sampling Algorithm
###############################################################################
########################################################   by Adriana F. Ch?vez 
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){ 
      source("./dCDDM.R") 
      source("./sim_auxiliarFunctions.R")
  }
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
          baseColor <- rgb(0,0,0,0.5)
          a <- scatterplot3d(x.C, y.RT, diag(z.Dens), zlab="Density", color = baseColor, type="l",
                             xlab="Choices", ylab="RT", zlim = c(0, max(z.Dens, theta.Dens)))
          a$points3d(x.C,rev(y.RT), diag(z.Dens[,c(nLines:1)]),type="l", col = baseColor)
          for(i in 1:nLines){  
            a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l", col = baseColor)
            a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i],  col = baseColor, type="l")
          }
          # Add a density line corresponding to Choice = Theta
          a$points3d(x.theta, y.RT, theta.Dens, col = baseColor, type="l")
  }
  
  #baseColor <- rgb(0.75,0.75,0.75,0.3)
  baseColor <- rgb(0.75,0.75,0.75,1)
  M <- 500
  ARate_des <- 0.25
  ARate_obs <- ARate_des
  Mu <- c(predChoice, predRT)
  Sigma <- diag(c((2*pi)^2,10))
  warm_up <- TRUE
  
  while(warm_up){
    newSigma <- Sigma*((ARate_obs/ARate_des)^(0.5))
    Sigma <- newSigma
            if(plot){          
              z.Dens <- matrix(NA, nrow=nLines, ncol=nLines)
              for(c in 1:nLines){ for(t in 1:nLines){
                  z.Dens[c,t] <- dmvnorm(c(x.C[c],y.RT[t]),mean=Mu, sigma = Sigma)
              }}
              z.Dens <- (z.Dens*max.Density)/max(z.Dens)
              for(i in 1:nLines){  
                  a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l", col = baseColor)
                  a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i],  col = baseColor, type="l")
              }
            }
      u <- runif(M,0,1)
      valid.candidates <- 0
      test.cand <- c(NA,NA)
      while(valid.candidates<M){
        cand <- rmvnorm(M,Mu,Sigma)
        valid.RT <- which(cand[,2] > tzero)
        valid.candidates <- valid.candidates + length(valid.RT)
        test.cand <- rbind(test.cand, cand[valid.RT,])
      }
      cand <- test.cand[2:(M+1),]
      ratio.num <- dCDDM(cand,drift, theta, tzero, boundary)
      ratio.den <- (dmvnorm(cand,Mu,Sigma)*max.Density)/dmvnorm(Mu,Mu,Sigma)
      ratio <- ratio.num/ratio.den
      change <- ratio>u
      ARate_obs <- mean(change)
      warm_up <- ARate_obs<ARate_des
  }
  
  n.keep <- 0
  samples <- matrix(NA, nrow=1, ncol=2)
  while(n.keep < n){
        valid.candidates <- 0
        test.cand <- c(NA,NA)
        while(valid.candidates<n){
              cand <- rmvnorm(n,Mu,Sigma)
              valid.RT <- which(cand[,2] > tzero)
              valid.candidates <- valid.candidates + length(valid.RT)
              test.cand <- rbind(test.cand, cand[valid.RT,])
        }
        cand <- matrix(test.cand[2:(n+1),], ncol=2, byrow = FALSE)
        cand[,1] <- cand[,1] %% 2*pi
        
        # ratio.num <- dCDDM(cand,drift, theta, tzero, boundary)
        # ratio.den <- (dmvnorm(cand,Mu,Sigma)*max.Density)/dmvnorm(Mu,Mu,Sigma)
        # ratio <- ratio.num/ratio.den
        # rej.crit <- runif(n,0,1)  
        # keep <- (eval >= rej.crit)
        # if(plot){
        #   a$points3d(cand[!keep,1], cand[!keep,2], ratio.den[!keep],
        #              col = rgb(204/255,0,0,0.2), pch = 16, cex = 0.3)
        #   a$points3d(cand[keep,1], cand[keep,2], ratio.num[keep],
        #              col = rgb(0,204/255,0,0.5), pch = 16, cex = 0.3)
        # }
        
        eval <- dCDDM(cand,drift,theta,tzero,boundary)
        max.Density <- (dmvnorm(cand,Mu,Sigma)*max.Density)/dmvnorm(Mu,Mu,Sigma)
        rej.crit <- runif(n,0,max.Density)  
        keep <- (eval >= rej.crit)
        
        if(plot){
          a$points3d(cand[!keep,1], cand[!keep,2], rej.crit[!keep],
                     col = rgb(204/255,0,0,0.2), pch = 16, cex = 0.3)
          a$points3d(cand[keep,1], cand[keep,2], rej.crit[keep],
                     col = rgb(0,204/255,0,0.5), pch = 16, cex = 0.3)
        }
        
        n.keep <- n.keep + sum(keep)
        samples <- rbind(samples, cand[keep,])
  }
  colnames(samples) <- c("Choice","RT")
  samples <- samples[2:(n+1),]
  return(samples)
}

# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
          par <- list("drift" = 1, 
                      "theta" = pi,
                      "tzero" = 0.1,
                      "boundary" = 7)
          n <- 5000
          sample.Reject.cddm(1000,par, plot=TRUE)  }