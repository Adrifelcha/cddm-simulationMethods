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
# Calculate expected Choice variance according to the EZ-CDDM
###############################################################################
ez_var <- ezcddm_VCA(drift = drift, boundary = par$boundary, tzero = par$tzero)

vM_Mu <- theta
vM_Kappa <- drift * boundary
###############################################################################
# Calculate variance of Choice according to the von Mises distribution
###############################################################################
var_choice <- rep(NA, 1000)
mean_choice <- rep(NA, 1000)
sd_choice <- rep(NA, 1000)
for(i in 1:1000){
  x <- rvonmises(n = 1000, mu = circular(vM_Mu), kappa = vM_Kappa)
  var_choice[i] <- var(x)
  mean_choice[i] <- mean(x)
  sd_choice[i] <- sd(x)
}

# Create PDF with appropriate dimensions for side-by-side plots
pdf(file = paste0(here::here("results", "varChoice_test.pdf")), 
    width = 8,    # Increased width to accommodate two plots
    height = 4)   # Height for a single row of plots

# Set up the plotting area for two side-by-side plots
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Plot 1: Distribution of means
plot(jitter(rep(1, length(mean_choice))), mean_choice,
     xlim = c(0.5, 1.5),
     ylim = range(c(mean_choice, theta)),
     xlab = "",
     ylab = "Mean angle",
     xaxt = "n")
abline(h = theta, col = "red", lwd = 2)
axis(1, at = 1, labels = "Simulated\nsamples")
title("Distribution of means")

# Plot 2: Distribution of variances
plot(jitter(rep(1, length(var_choice))), var_choice,
     xlim = c(0.5, 1.5),
     ylim = range(c(var_choice, ez_var)),
     xlab = "",
     ylab = "Variance",
     xaxt = "n")
abline(h = ez_var, col = "red", lwd = 2)
axis(1, at = 1, labels = "Simulated\nsamples")
title("Distribution of variances")

dev.off()