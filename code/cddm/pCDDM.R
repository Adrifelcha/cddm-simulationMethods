test = FALSE
source("./dCDDM.R")

# Some data
par <- list("drift" = 1, 
            "theta" = pi,
            "tzero" = 0.1,
            "boundary" = 3)
lower.C <- 0
upper.C <- 2*pi
lower.RT <- tzero
upper.RT <- 20
plot = TRUE
max.RT = 20
kappa = 200

# Write Trapezoid N.I. algorithm
numInt.tpz.cddm <- function(lower.C, upper.C,
                            lower.RT, upper.RT,
                            par, max.RT = 20, 
                            kappa=300, plot=FALSE){
  no.Dim <- 2
  drift <- par$drift
  theta <- par$theta
  tzero <- par$tzero
  boundary <- par$boundary

  if(plot){
    nSupp <- 100
    base.C <- c(0, 2*pi)
    base.RT <- c(tzero, max.RT)
    support.C <- seq(base.C[1],base.C[2],length.out=nSupp)
    support.RT <- seq(base.RT[1],base.RT[2],length.out=nSupp)
    z <- dCDDM(cbind(support.C,support.RT),drift, theta, tzero, boundary)
    z.top <- max(z)
    a <- scatterplot3d(support.C, support.RT, z, 
                       xlim=base.C, ylim=base.RT, zlim=c(0,z.top),
                       xlab="Choices", ylab="RT", zlab="Density",
                       color="blue", type="l", lwd=2)
    legend("topright", paste("No. bins =",kappa), 
           cex = 0.75)
  }
  
  bin.C <- seq(lower.C,upper.C,length.out=kappa)
  sidesLength.C <- bin.C[2]-bin.C[1]
  bin.RT <- seq(lower.RT,upper.RT,length.out=kappa)
  sidesLength.RT <- bin.RT[2]-bin.RT[1]
  base.area <- sidesLength.C*sidesLength.RT
  bin.area <- rep(NA,(kappa-1)^2)
  bin.count <- 1
  for(b.c in 2:kappa){
      for(b.rt in 2:kappa){
          C.space  <- c(bin.C[b.c-1],bin.C[b.c])
          RT.space <- c(bin.RT[b.rt-1],bin.RT[b.rt])
          v.height <- matrix(NA,nrow=2,ncol=2)
          for(rt in 1:2){
              for(c in 1:2){
              vertix <- c(C.space[c],RT.space[rt])
              v.density <- dCDDM(vertix,drift,theta,tzero,boundary)
              v.height[c,rt] <- max(c(0,v.density))
              }
          }
          height <- mean(v.height)
          bin.area[bin.count] <- height*base.area
          # VA = (1/6)*base.area*sum(as.vector(v.height)[-1])
          # VB = (1/6)*base.area*sum(as.vector(v.height)[-2])
          # VC = (1/6)*base.area*sum(as.vector(v.height)[-3])
          # VD = (1/6)*base.area*sum(as.vector(v.height)[-4])
          # bin.area[bin.count] <- (VA+VC+VB+VD)/2
          bin.count <- bin.count + 1
      }
  }
  total.area <- sum(bin.area)
  return(total.area)
}

# Test function
numInt.tpz.cddm(lower.C, upper.C,
                lower.RT, upper.RT,
                par = par, max.RT = 20, 
                kappa=300, plot=FALSE)

# Use Trapezoid Numeric integration to compute CDF
pCDDM <- function(data,par,plot=FALSE){
  lower.C <- 0
  lower.RT <- par$tzero
  upper.C <- data[1]
  upper.RT <- data[2]
  area <- numInt.tpz.cddm(lower.C, upper.C,
                          lower.RT, upper.RT,
                          par, max.RT = 20, 
                          kappa=300, plot=FALSE)
  return(area)
} 

# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
    par <- list("drift" = 1, 
                "theta" = pi,
                "tzero" = 0.1,
                "boundary" = 7)
    n <- 500
    C <- runif(1,0,2*pi)
    RT <- runif(1,0,15)
    data <- c(C,RT)
    pCDDM(data,par,plot=TRUE)
}