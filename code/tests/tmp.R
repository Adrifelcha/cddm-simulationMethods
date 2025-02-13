mu1 = 1
mu2 = 1
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


use <- 20
x <- data_ordered[use,c(1,2)]


pdf("tmp.pdf")
p1b <- pCDDM(x, drift, theta, tzero, boundary, method="grid", show=TRUE)
dev.off()



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
