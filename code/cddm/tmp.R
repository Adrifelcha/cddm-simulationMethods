set.seed(123)
N <- 1000
c <- rnorm(N,0,1)
d <- cbind(runif(N,0,2*pi), runif(N,0,2))



start1 <- Sys.time()
test1b <- pnorm(c,0,1)
end1 <- Sys.time()

start2 <- Sys.time()
test2d <- pCDDM(d, 1, pi/4, 0.1, 4, method="monte_carlo")
end2 <- Sys.time()

start3 <- Sys.time()
test3 <- pCDDM(d, 1, pi/4, 0.1, 4, method="monte_carlo")
end3 <- Sys.time()

print(end1 - start1) # Time difference of 0.03193784 secs
print(end2 - start2) # Time difference of 1.050278 mins
print(end3 - start3) # Time difference of 0.001000001 secs




