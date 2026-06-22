#############################################################
# Load libraries and custom functions
#############################################################
cat("Loading R libraries...\n")
library("here")
library("circular")

cat("\nLoading custom function scripts from /code/cddm...\n\n")
source(here("code", "cddm", "sim_randomWalk.R"))        
r_files <- list.files(path = here("code", "cddm"), 
                      pattern = "\\.R$", 
                      full.names = TRUE)
for(file in r_files) {
    source(file)
}
source(here("code", "general_functions", "eCDF.R"))

mu1 = -3
mu2 = 0
boundary = 3
tzero = 0.1

drift <- sqrt(mu1^2 + mu2^2)
theta <- atan2(mu2, mu1)

data <- rCDDM_Reject(n = 10000, par = list(mu1 = mu1, mu2 = mu2, boundary = boundary, tzero = tzero))


# 1. Order the data by RT and then Choice
data_ordered <- data[order(data$RT, data$Choice), ]

# Convert to matrix if needed
data_matrix <- as.matrix(data_ordered)

gold_eCDF = myECDF(data_matrix)

data_ordered$eCDF <- gold_eCDF


use <- 10000
x <- data_ordered[use,c(1,2)]


filename <- sprintf(here("results", "grid_pCDDM.pdf"))
pdf(filename)
p1b <- pCDDM(x, drift, theta, tzero, boundary, method="grid", show=TRUE)
dev.off()

filename2 <- sprintf(here("results", "monte_carlo_pCDDM.pdf"))
pdf(filename2)
p1b <- pCDDM(x, drift, theta, tzero, boundary, method="monte_carlo", show=TRUE)
dev.off()




#dCDDM(c(pi,6), drift, theta, tzero, boundary)

start_time <- Sys.time()
p1 <- pCDDM(x, drift, theta, tzero, boundary, method="grid")
end_time <- Sys.time()

start_time2 <- Sys.time()
p2 <- pCDDM(x, drift, theta, tzero, boundary, method="monte_carlo")
end_time2 <- Sys.time()
cat("Grid method result:", round(p1, 4), "\n")
cat("MCMC method result:", round(p2, 4), "\n")
cat("eCDF result:", round(data_ordered$eCDF[use], 4), "\n")
difftime(end_time, start_time, units = "secs")
difftime(end_time2, start_time2, units = "secs")
