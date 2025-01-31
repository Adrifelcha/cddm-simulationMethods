###############################################################################
###############################################################################
#####   Defining the bivariate probability density function for the CDDM
###############################################################################
########################################################   by Adriana F. Chavez   

# Compute the bivariate density and produce a suiting output for the data
dCDDM <- function(data, drift, theta, tzero, boundary) {
    # Input validation
    if(drift < 0) stop("drift must be non-negative")
    if(boundary <= 0) stop("boundary must be positive")
    if(tzero < 0) stop("tzero must be non-negative")
    
    if(is.vector(data)) {
        N <- 1
        # Ensure non-negative output
        pdf <- max(0, cddm.pdf.villarreal(data, drift, theta, tzero, boundary))
    } else {
        N <- nrow(data)
        pdf <- rep(NA, N)
        for(i in 1:N) {
            # Ensure non-negative output for each point
            pdf[i] <- max(0, cddm.pdf.villarreal(data[i,], drift, theta, tzero, boundary))
        }
    }
    return(pdf)
}

#################
# Test/Examples
#################
# Test function
if(!exists("test")){  test <- TRUE     }
if(test){
  n <- 10
  A <- runif(n, 0, 2*pi)
  B <- rexp(n,2)
  data <- cbind(A,B)
  drift <- 1
  theta <- pi
  tzero <- 0.1
  boundary <- 7
  dCDDM(data, drift, theta, tzero, boundary)
  dCDDM(data[1,], drift, theta, tzero, boundary)}