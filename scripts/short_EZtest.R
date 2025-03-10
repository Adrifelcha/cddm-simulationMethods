source(here::here("scripts", "load_libraries.R"))

###############################################################################
# Select arbitrary parameter values
###############################################################################
mu1 <- 0.5
mu2 <- 0.5
boundary <- 2.5
tzero <- 0.1
par <- list(mu1 = mu1, mu2 = mu2,
            boundary = boundary,
            tzero = tzero)
drift <- sqrt(sum(par$mu1^2, par$mu2^2))
par$drift <- drift
theta <- atan2(par$mu2, par$mu1)
par$theta <- theta

###############################################################################
# Calculate expected Choice variance according to the EZ-CDDM
###############################################################################
ez_var <- ezcddm_VCA(drift = drift, boundary = par$boundary)

vM_Mu <- theta
vM_Kappa <- drift * boundary
###############################################################################
# Calculate variance of Choice according to the von Mises distribution
###############################################################################
# Set number of iterations
m <- 200

var_choice_vonM <- rep(NA, m)
var_choice_vonM2 <- rep(NA, m)
mean_choice_vonM <- rep(NA, m)
sd_choice_vonM <- rep(NA, m)

var_choice_RW <- rep(NA, m)
var_choice_RW2 <- rep(NA, m)
var_choice_RW3 <- rep(NA, m)
mean_choice_RW <- rep(NA, m)
sd_choice_RW <- rep(NA, m)
for(i in 1:m){
  # Print progress message
  cat(sprintf("\rProcessing iteration %d of %d (%.1f%%)", i, m, i/m*100))
  
  # Von Mises samples
  x_vonM <- rvonmises(n = 300, mu = circular(vM_Mu), kappa = vM_Kappa)
  
  # Use circular variance for von Mises
  x_vonM_circular <- circular(x_vonM)
  var_choice_vonM[i] <- 1 - rho.circular(x_vonM_circular) # Circular variance
  var_choice_vonM2[i] <- var.circular(x_vonM_circular)
  mean_choice_vonM[i] <- mean.circular(x_vonM_circular)
  
  # Random walk samples
  x_RW <- rCDDM_RandomWalk(n = 300, par = par)
  x_RW <-  x_RW$bivariate.data
  
  # Convert random walk angles to circular
  x_RW_circular <- circular(x_RW[,1])
  var_choice_RW[i] <- 1 - rho.circular(x_RW_circular)  
  var_choice_RW2[i] <- var.circular(x_RW_circular)  
  var_choice_RW3[i] <- var(x_RW[,1])
  mean_choice_RW[i] <- mean.circular(x_RW_circular)
}
cat("\nSimulation complete!\n")


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