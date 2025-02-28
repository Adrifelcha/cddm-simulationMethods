library(circular) 
source(here::here("code", "cddm", "sim_auxiliarFunctions.R"))

###############################################################################
# Select arbitrary parameter values
###############################################################################
mu1 <- 0.5
mu2 <- 0.5
boundary <- 5
tzero <- 0.1
par <- list(mu1 = mu1, mu2 = mu2,
            boundary = boundary,
            tzero = tzero)
drift <- sqrt(sum(par$mu1^2, par$mu2^2))
theta <- atan2(par$mu2, par$mu1)

###############################################################################
# Calculate expected variance of RT according to the EZ-CDDM
###############################################################################
ez <- ezcddm_VRT(drift = drift, boundary = par$boundary)

vM_Mu <- theta
vM_Kappa <- drift * boundary
###############################################################################
# Calculate variance of RT according to the von Mises distribution
###############################################################################
var_rt <- rep(NA, 1000)
mean_rt <- rep(NA, 1000)
for(i in 1:1000){
  x <- rvonmises(n = 1000, mu = circular(vM_Mu), kappa = vM_Kappa)
  var_rt[i] <- var(x)
  mean_rt[i] <- mean(x)
}

mean(var_rt)
mean(mean_rt)


