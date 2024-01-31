###############################################################################
###############################################################################
#####   A set of functions made to plot the data generated from simulations
###############################################################################
########################################################   by Adriana F. Chavez   
############## Load custom functions
source("../general_functions/customFunctions.R")


iter <- 100000 #Iterations

x<-NULL     #Empty array to start sampling
x[1] <- 0   #Initial value, arbitrary
change.counter <- 0    
for(i in 2:iter){
  x.cand <- x[i-1]+rnorm(1,0,5)      #Extract candidate from wide normal
  alpha <- dnorm(x.cand,0,1)/dnorm(x[i-1],0,1)   #Likelihood ratio
  keep <- min(alpha,1)   #Use minimun to make decisions
  if(keep==1){    #If ratio is greater than 1, keep candidate
    x[i] <- x.cand
    change.counter <- change.counter+1   #Update counter
  }else{        #If ratio is not greater than 1...
    update <- rbinom(1,1,alpha)     #Update x with probability alpha
    if(update==0){
      x[i] <- x[i-1]    #Either we keep the current value
    }else{
      x[i]<-x.cand      #Or we update to candidate
      change.counter<-change.counter+1}}
}