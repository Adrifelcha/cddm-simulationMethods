library(brms)

cubrt <- function(x){   x^(1/3)   }

# Reference: https://statproofbook.github.io/P/exg-mome.html
exGauss_est <- function(rt){
  N <- length(rt)
  # Sample mean and variance
  s.mean <- mean(rt)
  s.var  <- var(rt)
  # Sample skewness
  num <- (sum((rt-s.mean)^3))*(1/N)
  denom <- ((sum((rt-s.mean)^2))*(1/N))^(3/2)
  s.skew <- num/denom
  #s.skew <- skew(rt)
  
  # Estimate Mu parameter using MOM
  # Note: This result is close to boundary/drift
  term1 <- (s.skew*(s.var^(3/2)))*0.5
  est.mu <- s.mean-cubrt(term1)
  # Estimate Sigma parameter using MOM
  # Note: This result is close to tzero/drift
  term2 <- 1-cubrt((s.skew^2)*0.25)
  est.sigma <- sqrt(s.var*term2)
  # Estimate Lambda parameter using MOM
  term3 <- 2/(s.skew*(s.var^(3/2)))
  est.lambda <- cubrt(term3)
  
return(list("mu"    = est.mu,
            "sigma" = est.sigma,
            "lambda"= est.lambda))
}

# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
        if(!exists("rt")){
            source("./sim_Reject.R")
            # Arbitrary set of parameter values
            par <- list("drift" = 3.5, 
                        "theta" = pi,
                        "tzero" = 0.2,
                        "boundary" = 2)
            n <- 2000 # No. samples
            rt <- sample.Reject.cddm(n,par,max.RT = 5)[,2]
        }
  g <- exGauss_est(rt)
  fake_rt <- rexgaussian(n,g$mu,g$sigma,1/g$lambda)
  par(pty="s", mfrow=c(1,2),mar = c(0, 1, 0, 0)) 
  hist(fake_rt, ann = F, axes = F)
  mtext("Replicated RT", 3, f = 2)
  axis(1, seq(0,2,0.2), seq(0,2,0.2))
  hist(rt, ann = F, axes = F)
  mtext("Observed RT", 3, f = 2)
  axis(1, seq(0,2,0.2), seq(0,2,0.2))
}
