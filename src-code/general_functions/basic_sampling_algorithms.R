###############################################################################
###############################################################################
#####   A set of simplified sampling algorithms for:
#####     1) Numeric integration of the cdf of both the DDM and CDDM models
#####     2) Rejection algorithms to sample observations for DDM and CDDM
#####   Functions for both Normal and Bivariate normal data are presented
###############################################################################
########################################################   by Adriana F. Chavez   
library(mnormt)
library(scatterplot3d) 

###############################################################################
#####   1.1  Unidimensional Normal data: CDF approximation
###############################################################################
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.1.a Trapezoid Numerical Integration for the Normal distribution
#       Calculate the area under the curve between two points
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
numInt.tpz.normal <- function(lower.bound, upper.bound, par,
                              kappa=NA, plot=FALSE){
    # Load variables
    mean <- par$mean
    sd <- par$sd
    if(is.na(kappa)){     total.width <- upper.bound-lower.bound
                          kappa <- total.width*20                      
    }
    # If plot is requested, draw the base Normal distribution
    if(plot){   width <- (3*sd)
                ticks <- round(seq(mean-width,mean+width,length.out=10),2)
                support <- seq(mean-width,mean+width,length.out=9999)
                plot(support,dnorm(support,mean,sd), type="l", ann=F, axes=F, col=1)
                axis(1,ticks,ticks)
    }
    
    # Define bins between the lower and upper bounds specified
    bins <- seq(lower.bound,upper.bound,length.out=kappa)
    bin.base = bins[2] - bins[1]   # Bin length
    # Compute the height of the normal distribution at each bin endpoint
    d <- dnorm(bins,par$mean,par$sd)
    # We want to locate the first bin d>0. (Discard really low bins where d=0)
    b <- which(d>0)[1]
    # Set up an initial value for the total area
    total.area <- 0
    # Update the total.area, only if the density is > 0 in the interval 
    if(!is.na(b)){   
                      bin.area <- rep(NA,kappa)  # Empty vector to store bin areas
                      # Add up the areas of each bin
                      # Stop at upper.bound or whenever area >= 1
                      while((total.area<1&b<=kappa)==TRUE){
                            # Get lower and upper bound of any bin
                            low.x = bins[b-1]; up.x  = bins[b]
                            # Get densities at each point
                            density.x1 <- d[b-1]; density.x2 <- d[b]
                            # Compute trapezoid area
                            height <- (density.x1+density.x2)/2
                            bin.area[b-1] <- height*bin.base
                            # Draw the trapezoid
                            if(plot){           polygon(x=c(low.x,up.x,up.x,low.x),
                                                        y=c(0,0,density.x1,density.x2),
                                                        col = (b %% 2)+1)
                            }
                            # Update total.area
                            total.area <- sum(bin.area,na.rm=TRUE)
                            b <- b+1   # Move to next bin
                      }
                      # Add legends to the plot
                      if(plot){ 
                        no.bins = sum(!is.na(bin.area))
                        legend("topright", paste("No. bins =",no.bins), 
                               cex = 0.75, bty ="n")
                      }
    }else{
                      if(plot){ 
                        legend("topright", "Region 'outside' of the curve", 
                               cex = 0.75, bty ="n")
                      }
    }

return(total.area)
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# 1.1.b  User friendly cdf funciton
#        Compute the approximate normal cdf for a single value of x
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
normal.cdf <- function(x,par,plot=FALSE){
  # Get a sufficiently low lower.bound
  lower.bound <- (par$mean-(par$sd*100))
  if(x<lower.bound){
   lower.bound  = x-lower.bound
  }
  # Compute the area using our trapezoid algorithm
  area <- numInt.tpz.normal(lower.bound,x,par,plot=plot)
  return(area)
}

#~~~~~~~~~~~~~~~#
# Test/Examples #
#~~~~~~~~~~~~~~~#
if(!exists("test")){  test <- TRUE     }
if(test){
          par <- list("mean" = 500, "sd" = 1)
          lower.bound <- 1
          upper.bound <- 5000
          kappa=NA
          plot=TRUE
          # Trapezoid Numerical Integration of whole support
          numInt.tpz.normal(lower.bound,upper.bound,par, plot=TRUE)
          # Compute various approximate cdf
          print(normal.cdf(10,par,plot=TRUE))
          print(normal.cdf(50,par,plot=TRUE))
          print(normal.cdf(500,par,plot=TRUE))
          print(normal.cdf(501,par,plot=TRUE))
          print(normal.cdf(502,par,plot=TRUE))
          print(normal.cdf(600,par,plot=TRUE))
}


###############################################################################
#####   1.2  Unidimensional Normal data: MH-MCMC sampling algorithm
###############################################################################
sample.MCMC.normal <- function(n, par, plot=FALSE){
      # Load variables
      mean <- par$mean
      sd <- par$sd
      # Define maximum rejection value (RejVal)
      max.D <- dnorm(mean,mean,sd)  # Density at mean (highest density point)
      height <- max.D*1.05          # RejVal ~ Uniform(0, 5% more than density at mean)
      # Define the range o candidate values to generate
      width <- (3.5*sd)
      base <- c(mean-width,mean+width)
      # If plot is requested, draw the base Normal distribution
      if(plot){
                plot.width <- (4*sd)    # Make plot slighly wider than candidate region
                x.low <- mean-plot.width
                x.top <- mean+plot.width
                ticks <- round(seq(x.low,x.top,length.out=10),2)
                support <- seq(x.low,x.top,length.out=9999)
                plot(support,dnorm(support,mean,sd), xlim=c(x.low,x.top),
                     type="l", ann=F, axes=F, col=1, ylim=c(0,height))
                axis(1,ticks,ticks)
                abline(h=height)
                abline(v=base)
      }
      # Start MH-MCMC sampling algorithm
      # This algorithm generates `n` samples from the specified Normal
      n.keep <- 0       # Start with 0 candidates stored
      samples <- NA     #   "    "  no     "        "
      n.try <- n        # Start by generating 'n' candidates
      while(n.keep < n){    # Until we reach a sample of size 'n'
            ######### First: Sample random candidates and rejectio values
            cand <- runif(n.try,base[1],base[2])    # We generate n.try candidate values
            eval <- dnorm(cand,mean,sd)           #    and compute their pdf
            rej.crit <- runif(n.try,0,height)     # We generate n.try Rej. Vals.
            # Keep the candidates whose density is greater than the corresponding Rej Val.
            keep <- (eval >= rej.crit)  
            ######### Second: Update
            n.pass <- sum(keep)      # Count number of accepted candidates in this run
            n.try <- n.try-n.pass    # Update number of samples missing to reach 'n'
            samples <- c(samples, cand[keep])   # Store candidate values approved
            n.keep <- length(samples)           # Count total no. of samples to test while() 
            ######### Third (Optional): Make plot
            if(plot){
                      points(cand[!keep],rej.crit[!keep],col="red", pch=16, cex=0.3)
                      points(cand[keep],rej.crit[keep],col="blue3", pch=16, cex=0.3)
            }
      }
  return(samples)
}

#~~~~~~~~~~~~~~~#
# Test/Example  #
#~~~~~~~~~~~~~~~#
if(!exists("test")){  test <- TRUE     }
if(test){
  par <- list("mean" = 10, "sd" = 1)
  n <- 5000
  sample.MCMC.normal(n,par, plot=TRUE)
}



###############################################################################
#####   2.2  Multivariate Normal data: MH-MCMC sampling algorithm
###############################################################################
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
