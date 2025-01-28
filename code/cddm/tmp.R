N <- 1000

c <- rnorm(N,0,1)

d <- cbind(runif(N,0,2*pi), runif(N,0,2))

start1 <- Sys.time()
test1 <- pnorm(c,0,1)
end1 <- Sys.time()

start2 <- Sys.time()
test2 <- pCDDM(d, 1, pi/4, 0.1, 4, method="monte_carlo")
end2 <- Sys.time()

print(end1 - start1)
print(end2 - start2)