###############################################################################
# EXPLORING THE VARIANCE OF THE CHOICE SPACE
###############################################################################
# In this script, we compare:
# 1. The variance of the choice space theoretically prescribed by the EZ-CDDM
# 2. The variance of the choices generated using a von Mises distribution sampler
# 3. The variance of the choices generated using a Random Walk emulator

# Load libraries and source functions
source(here::here("scripts", "load_libraries.R"))
###############################################################################
###############################################################################
# Test 1: Fixed parameter set
# We use a FIXED PARAMETER SET and compare the variance of the choice space
# between the EZ-CDDM, the von Mises distribution, and the Random Walk emulator
###############################################################################
###############################################################################
# Set evidence accumulation parameters
boundary <- 2.5
tzero <- 0.1
# Drift vector
mu1 <- 0.5
mu2 <- 0.5
# Get corresponding polar coordinates
polar <- rectToPolar(mu1, mu2)
drift <- polar$dLength
theta <- polar$dAngle
# Store parameters in a list
par <- list(mu1 = mu1, mu2 = mu2,
            drift = drift,
            theta = theta,
            boundary = boundary,
            tzero = tzero)

###############################################################################
# THEORETICAL VARIANCE OF THE CHOICE SPACE (EZ-CDDM)
###############################################################################
ez_var <- ezcddm_VCA(drift = drift, boundary = par$boundary)


####################################################################################
# OBSERVED VARIANCE OF THE CHOICE SPACE (VON MISES SAMPLER AND RANDOM WALK EMULATOR)
####################################################################################
# Define simulation parameters
m <- 300  # Number of repetitions
n <- 500  # Number of trials sampled per repetition
# Define von Mises parameters
vM_Mu <- theta
vM_Kappa <- drift * boundary

# Empty vectors to store results
empty_vector <- rep(NA, m)
# Empty vectors to store results - von Mises
mean_vonM <- empty_vector
sd_vonM <- empty_vector
var_vonM <- empty_vector
var_vonM2 <- empty_vector
var_vonM3 <- empty_vector
# Empty vectors to store results - Random Walk
mean_RW <- empty_vector
sd_RW <- empty_vector
var_RW <- empty_vector
var_RW2 <- empty_vector
var_RW3 <- empty_vector

# Loop through repetitions
for(i in 1:m){  
          cat(sprintf("\rProcessing iteration %d of %d (%.1f%%)", i, m, i/m*100))          
          
          # Von Mises samples
          ######################
          x_vonM <- rvonmises(n = n, mu = circular(vM_Mu), kappa = vM_Kappa)    
          # Pass choices to circular space
          x_vonM_circular <- circular(x_vonM)
          mean_vonM[i] <- mean.circular(x_vonM_circular)
          sd_vonM[i] <- sd.circular(x_vonM_circular)
          var_vonM[i] <- var.circular(x_vonM_circular)
          var_vonM2[i] <- var(x_vonM_circular)
          var_vonM3[i] <- var(x_vonM)
          
          # Random walk samples
          ######################
          x_RW_full <- rCDDM_RandomWalk(n = n, par = par)
          x_RW <- x_RW_full$bivariate.data  
          # Pass choices to circular space
          x_RW_circular <- circular(x_RW[,1])
          mean_RW[i] <- mean.circular(x_RW_circular)
          sd_RW[i] <- sd.circular(x_RW_circular)
          var_RW[i] <- var.circular(x_RW_circular)
          var_RW2[i] <- var(x_RW_circular)
          var_RW3[i] <- var(x_RW[,1])  
}
cat("\nSimulation complete!\n")

range(var_vonM)
range(var_RW)
range(var_RW2)
range(var_RW3)

# Create PDF with appropriate dimensions for side-by-side plots
pdf(file = paste0(here::here("results", "varChoice_test1.pdf")), 
    width = 8, height = 4)   

# Set up the plotting area for two side-by-side plots
par(mfrow = c(1, 2), mar = c(4, 4, 2, 1))

# Plot 1: Distribution of means
plot(c(jitter(rep(1, length(mean_choice_vonM))), jitter(rep(2, length(mean_choice_RW)))),
     c(mean_choice_vonM, mean_choice_RW),
     xlim = c(0.5, 2.5),
     ylim = range(c(mean_choice_vonM, mean_choice_RW, theta)),
     xlab = "",
     ylab = "Mean angle",
     xaxt = "n",
     col = c(rep("blue", length(mean_choice_vonM)), rep("darkgreen", length(mean_choice_RW))),
     pch = 19)
abline(h = theta, col = "red", lwd = 2)
axis(1, at = 1:2, labels = c("von Mises", "Random Walk"))
title("Distribution of means")

# Plot 2: Distribution of variances
plot(c(jitter(rep(1, length(var_choice_vonM))), jitter(rep(2, length(var_choice_RW)))),
     c(var_choice_vonM, var_choice_RW),
     xlim = c(0.5, 2.5),
     ylim = range(c(var_choice_vonM, var_choice_RW, ez_var)),
     xlab = "",
     ylab = "Variance",
     xaxt = "n",
     col = c(rep("blue", length(var_choice_vonM)), rep("darkgreen", length(var_choice_RW))),
     pch = 19)
abline(h = ez_var, col = "red", lwd = 2)
axis(1, at = 1:2, labels = c("von Mises", "Random Walk"))
title("Distribution of variances")

# Add legend
legend("topright", 
       legend = c("von Mises", "Random Walk", "EZ-CDDM prediction"),
       col = c("blue", "darkgreen", "red"), 
       pch = c(19, 19, NA),
       lty = c(NA, NA, 1),
       lwd = c(NA, NA, 2))

dev.off()




###############################################################################
###############################################################################
# Test 2: Different parameter sets
# We generate n datasets, each using a different parameter set, and compare the 
# variance of the choice space between the EZ-CDDM, the von Mises distribution, 
# and the Random Walk emulator
###############################################################################
###############################################################################
# Define simulation parameters
m <- 300  # Number of repetitions
n <- 300  # Number of trials sampled per repetition

# Set evidence accumulation parameters
boundary <- 2.5
tzero <- 0.1
# Drift vector
mu1 <- runif(m, 0.5, 1.5)
mu2 <- runif(m, 0.5, 1.5)
# Get corresponding polar coordinates
polar <- rectToPolar(mu1, mu2)
drift <- polar$dLength
theta <- polar$dAngle
# Start parameter list
par <- list(boundary = boundary, tzero = tzero)
###############################################################################
# THEORETICAL VARIANCE OF THE CHOICE SPACE (EZ-CDDM)
###############################################################################
ez_var <- ezcddm_VCA(drift = drift, boundary = par$boundary)

####################################################################################
# OBSERVED VARIANCE OF THE CHOICE SPACE (VON MISES SAMPLER AND RANDOM WALK EMULATOR)
####################################################################################
# Define von Mises parameters
vM_Mu <- theta
vM_Kappa <- drift * boundary

# Empty vectors to store results
empty_vector <- rep(NA, m)
# Empty vectors to store results - von Mises
mean2_vonM <- empty_vector
sd2_vonM <- empty_vector
var2_vonM <- empty_vector
# Empty vectors to store results - Random Walk
mean2_RW <- empty_vector
sd2_RW <- empty_vector
var2_RW <- empty_vector

# Loop through repetitions
for(i in 1:m){  
          cat(sprintf("\rProcessing iteration %d of %d (%.1f%%)", i, m, i/m*100))          
                              # Von Mises samples
          ######################
          x_vonM <- rvonmises(n = n, mu = circular(vM_Mu[i]), kappa = vM_Kappa[i])    
          # Pass choices to circular space
          x_vonM_circular <- circular(x_vonM)
          mean2_vonM[i] <- mean.circular(x_vonM_circular)
          sd2_vonM[i] <- sd.circular(x_vonM_circular)
          var2_vonM[i] <- var.circular(x_vonM_circular)
          
          # Random walk samples
          ######################
          # Update parameter list
          par$mu1 <- mu1[i]
          par$mu2 <- mu2[i]
          par$drift <- drift[i]
          par$theta <- theta[i]
          # Generate random walk samples
          x_RW_full <- rCDDM_RandomWalk(n = n, par = par)
          x_RW <- x_RW_full$bivariate.data  
          # Pass choices to circular space
          x_RW_circular <- circular(x_RW[,1])
          mean2_RW[i] <- mean.circular(x_RW_circular)
          sd2_RW[i] <- sd.circular(x_RW_circular)
          var2_RW[i] <- var.circular(x_RW_circular)
}
cat("\nSimulation complete!\n")

# Create PDF with appropriate dimensions for 2x2 grid of plots
pdf(file = paste0(here::here("results", "varChoice_test2.pdf")), 
    width = 10, height = 8)   

# Set up the plotting area for 2x2 grid
par(mfrow = c(2, 2), mar = c(4, 4, 3, 2), oma = c(0, 0, 2, 0))

# Plot 1: von Mises variance vs EZ-CDDM variance (top left)
plot(var2_vonM, ez_var,
     ylab = "EZ-CDDM Variance",
     xlab = "von Mises Variance",
     main = "von Mises vs EZ-CDDM",
     pch = 19, col = "blue")
abline(0, 1, col = "red", lwd = 2)

# Plot 2: Random Walk variance vs EZ-CDDM variance (top right)
plot(var2_RW, ez_var,
     ylab = "EZ-CDDM Variance",
     xlab = "Random Walk Variance",
     main = "Random Walk vs EZ-CDDM",
     pch = 19, col = "darkgreen")
abline(0, 1, col = "red", lwd = 2)

# Plot 3: von Mises variance vs Random Walk variance (bottom left)
plot(var2_vonM, var2_RW,
     xlab = "von Mises Variance",
     ylab = "Random Walk Variance",
     main = "von Mises vs Random Walk",
     pch = 19, col = "purple")
abline(0, 1, col = "red", lwd = 2)

# Plot 4: Empty plot with parameter information (bottom right)
plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "", main = "Simulation Parameters")
text_x <- 0
text_y <- 0.5
text_step <- 0.1
params_text <- c(
  paste("Boundary:", boundary),
  paste("Non-decision time:", tzero),
  paste("Drift range (mu1):", round(min(mu1), 2), "to", round(max(mu1), 2)),
  paste("Drift range (mu2):", round(min(mu2), 2), "to", round(max(mu2), 2)),
  paste("Number of parameter sets (m):", m),
  paste("Trials per parameter set (n):", n)
)

for (i in 1:length(params_text)) {
  text(text_x, text_y - (i-1)*text_step, params_text[i], adj = c(0.5, 0.5))
}

# Add overall title
mtext("Exploring Choice Variance", outer = TRUE, line = 0, cex = 1.5)

dev.off()

cat("PDF created successfully at:", paste0(here::here("results", "varChoice_test2.pdf")), "\n")

