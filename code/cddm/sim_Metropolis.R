###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####                     Metropolis Sampling Algorithm
###############################################################################
########################################################   by Adriana F. Ch?vez 
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){ source("./dCDDM.R") }
library(mvtnorm)


# Auxiliary function 2: Generate a plot of how the approximation works
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
embedded_plot <- function(bin.C, bin.RT, kappa, kappa.RT, total, 
                          rt, tzero, theta, drift, boundary,density_matrix, 
                          rgb1 = c(0.2,0.5,0.6), rgb2 = c(0.3,0.4,0.7)){
  par(pty="s")          
  par(mfrow=c(1,2),mar = c(0, 0, 0, 0)) 
  # Plot 1 ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Base plot: Sketch the bivariate density curve
  # Define line positions over the X and Y dimensions
  nLines <- 30
  x.C <- seq(0,2*pi,length.out=nLines)            # X: Choices
  y.RT <- seq(tzero,max.RT,length.out=nLines) # Y: RT
  x.theta <- rep(theta,nLines)                    # Make sure to add theta
  # Compute the density at each intersection point
  z.Dens <- matrix(NA, nrow=nLines, ncol=nLines)
  for(c in 1:nLines){ for(t in 1:nLines){
    z.Dens[c,t] <- dCDDM(c(x.C[c],y.RT[t]),drift,theta,tzero,boundary)
  }}
  # Compute the density at Choice = theta
  theta.Dens <- dCDDM(cbind(x.theta,y.RT),drift,theta,tzero,boundary)
  # Draw bivariate density curve
  baseColor <- rgb(0,0,0,0.2)
  a <- scatterplot3d(x.C, y.RT, diag(z.Dens), zlab="Density", color = baseColor, type="l",
                     xlab="Choices", ylab="RT", zlim = c(0, max(z.Dens, theta.Dens)))
  a$points3d(x.C,rev(y.RT), diag(z.Dens[,c(nLines:1)]),type="l", col = baseColor)
  for(i in 1:nLines){  
    a$points3d(rep(x.C[i],nLines), y.RT, z.Dens[i,], type="l", col = baseColor)
    a$points3d(x.C, rep(y.RT[i],nLines), z.Dens[,i],  col = baseColor, type="l")
  }
  # Add a density line corresponding to Choice = Theta
  a$points3d(x.theta, y.RT, theta.Dens, col = "red", type="l")
  legend("topright", c("p( RT | theta )"), col="red", cex=0.6, lwd=1, bty = "n")
  #  plot: Color the density area 
  nLines2 <- nLines*2
  k.C  <- floor(seq(1,kappa,length.out=nLines2))
  k.RT <- floor(seq(1,kappa.RT,length.out=nLines2))  
  for(i in k.C){  a$points3d(rep(bin.C[i],nLines2), bin.RT[k.RT], type="l", 
                             density_matrix[i,k.RT], col = rgb(rgb1[1],rgb1[2],rgb1[3],0.5))
  }
  for(i in k.RT){ a$points3d(bin.C[k.C], rep(bin.RT[i],nLines2), type="l",
                             density_matrix[k.C,i], col = rgb(rgb1[1],rgb1[2],rgb1[3],0.5))
  }
  for(i in 1:length(rt)){
    row <- which.min(abs(bin.C-rad[i]))-1
    col <- which.min(abs(bin.RT-rt[i]))-1
    a$points3d(rad[i], rt[i], density_matrix[row,col], pch=16, cex=0.3,
               col=rgb(rgb2[1],rgb2[2],rgb2[3],1))
  }
  mtext(paste("Total =", round(max(total),4)),3, adj = 1, cex=0.8)  
  # Plot 2  ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
  # Base plot: Sketch the bivariate density curve
  b <- scatterplot3d(rad, rt, total,  pch=16, xlab = "", ylab = "", zlab = "",
                     color = rgb(rgb2[1],rgb2[2],rgb2[3],0.5), cex.symbols = 0.5)
  mtext("Approximate CDF: Area under curve", f=2)
  mtext("Choices", side=1, line=0.5)
  mtext("RTs", side=4, line=-1, adj=0)
}


# Write a simple Rejection algorithm for the CDDM pdf
sample.Metropolis.cddm <- function(n, par, max.RT=NA, plot=FALSE){
  no.Dim <- 2
  drift <- par$drift;     theta <- par$theta
  tzero <- par$tzero;     boundary <- par$boundary
  
  predChoice <- theta
  predRT <- boundary/drift
  
  if(is.na(max.RT)){  max.RT <- predRT*4 }
  
  nTry.RT <- 30
  test.RT <- seq(tzero, max.RT, length.out=nTry.RT)
  test.data <- cbind(rep(theta,nTry.RT),test.RT)
  test.densities <- dCDDM(test.data,drift, theta, tzero, boundary)
  max.Density <- max(test.densities)
  height <- max.Density*1.1
  
  base.C <- c(0, 2*pi)
  base.RT <- c(tzero, max.RT)

  n.keep <- 0
  n.try <- n
  samples <- matrix(NA, nrow=1, ncol=no.Dim)
  Mu <- c(predChoice, predRT)
  Sigma <- diag(c(2,3))
  
  while(n.keep < n){
        valid.candidates <- 0
        test.cand <- c(NA,NA)
        while(valid.candidates<n.try){
              cand <- rmvnorm(n.try,Mu,Sigma)
              valid.RT <- which(cand[,2] > tzero)
              valid.candidates <- valid.candidates + length(valid.RT)
              test.cand <- rbind(test.cand, cand[valid.RT,])
        }
        cand <- cbind(test.cand[2:nrow(test.cand),1],
                      test.cand[2:nrow(test.cand),2])
        cand[,1] <- cand[,1] %% (2*pi)
        
        eval <- dCDDM(cand,drift, theta, tzero, boundary)
        rej.crit <- runif(nrow(cand),0,height)  
        keep <- (eval >= rej.crit)
        
        n.keep <- n.keep + sum(keep)
        samples <- rbind(samples, cand[keep,])
        
        # if(plot){
        #   a$points3d(cand[!keep,1], cand[!keep,2], rej.crit[!keep],
        #              col = "red", pch = 16, cex = 0.2)
        #   a$points3d(cand[keep,1], cand[keep,2], rej.crit[keep],
        #              col = "green", pch = 16, cex = 0.2)
        # }
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