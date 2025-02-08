mu1 <- 1
mu2 <- 1
boundary <- 5
tzero <- 0.1

x <- sample.RW.cddm(n = 10000, par = list(mu1 = mu1, mu2 = mu2, boundary = boundary, tzero = tzero))
data <- x$bivariate.data

# 1. Order the data by RT and then Choice
data_ordered <- data[order(data$RT, data$Choice), ]

# Convert to matrix if needed
data_matrix <- as.matrix(data_ordered)

gold_eCDF = myECDF(data_matrix)

data_ordered$eCDF <- gold_eCDF


use <- 10
x <- data_ordered[use,c(1,2)]

start_time <- Sys.time()
p1 <- pCDDM(x, drift, theta, tzero, boundary, method="grid")
end_time <- Sys.time()
cat("Grid method result:", round(p1, 4), "\n")
cat("eCDF result:", round(data_ordered$eCDF[use], 4), "\n")
difftime(end_time, start_time, units = "secs")