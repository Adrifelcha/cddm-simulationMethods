pCDDM <- function(data,par,plot=FALSE){
  
} 


# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
    par <- list("drift" = 1, 
                "theta" = pi,
                "tzero" = 0.1,
                "boundary" = 7)
    n <- 500
    C <- runif(n,0,2*pi)
    RT <- runif(n,0,15)
    data <- cbind(C,RT)
    pCDDM(data,par,plot=TRUE)
}