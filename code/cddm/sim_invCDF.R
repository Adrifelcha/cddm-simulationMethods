###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using a
#####   INVERSE PROBABILITY TRANSFORM ALGORITHM with a grid approximation
###############################################################################
########################################################   by Adriana F. Chavez 

par <- list("drift" = 1, 
            "theta" = pi,
            "tzero" = 0.1,
            "boundary" = 7)
n <- 100
max.RT <- 20

sample.invCDF.cddm <- function(n,par,max.RT){
  n.space <- 1000
  space.C  <- seq(0,2*pi,length.out=n.space)
  space.RT <- seq(0,max.RT,length.out=n.space)
  u <- runif(n,0,1)
}