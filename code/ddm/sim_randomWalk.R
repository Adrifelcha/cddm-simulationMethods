###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the DDM using the
#####                   RANDOM-WALK EMULATION METHOD 
###############################################################################
########################################################   by Adriana F. Ch?vez   
source("../general_functions/customFunctions.R")

###############################################################################
# Variable dictionary: ########################################################
# mu - Drift rate
# boundary - Boundary separation
# ndt - Non decision time
# bias - Bias towards the upper boundary anser (0-1)
# drift.Coeff - Within-trial variability on the sampling process
# dt - Step size ("delta-t")
# state - evidence accumulated at a given point in time
###############################################################################
trials = 100
boundary = 1
ndt = 0.5
mu = 0.5
beta = 0.5
mean.Mu = NA
sd.Mu = NA

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate the full random walk across many trials (for each trial, 
# keeps the full chain of coordinates visited and response times)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Base Function: Simulate the random walk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
ddm.randomWalk <- function(trials, beta, boundary, mu=NA, mean.Mu = NA, 
                           sd.Mu =NA, ndt=0.1, drift.Coeff=1, dt=0.00015){
  
  if(is.na(mu)){
        mu <- rnorm(trials,mean.Mu,sd.Mu)
  }
  
  sqDT <- sqrt(dt)
  iter <- round(15/dt)  # Maximum number of iterations on the random walk 
  state <- matrix(NA, nrow=iter, ncol=trials)   # States are saved in a 2D matrix
  finalT <- rep(NA,trials) # Empty vector to store RT (a.k.a. total number of iterations)
  additional_steps_needed <- rep(0,trials) 
  
  # Arrays to be used in simulation
  random_deviations <- rnorm(trials*iter,0,1)*(drift.Coeff*sqDT)   # Deviations from Mu
  motion <- matrix(random_deviations,nrow=iter,ncol=trials)        # Store deviations in matrix
  steps <- motion+(mu*dt)
  
  # Set initial state for every trial
  state[1,] <- beta*boundary   # Set initial point for every random-walk on each trial
  
  for(a in 1:trials){   
      ### Random walk per trial
      for(t in 2:iter){
        this.step <- steps[t,a]
        this.state <- state[t-1,a]+this.step
        state[t,a] <- this.state
        
        # Stop random-walk if boundary is passed
        if(this.state >= boundary | this.state <= 0){
          finalT[a] <- t+(ndt/dt)   #Total no. of iterations required on each trial
          break
        }
      }
    
    # Test whether the random-walk reached the boundary, and re-sample if not.
    not.finished <- is.na(finalT[a])
    if(not.finished){ additional_steps_needed[a] <- 1 }
    
    whileLoopNo <- 1
    while(not.finished){
      last_state <- state[t,a]   # Store last state
      state[,a] <- NA            # Reset random-walk
      state[1,a] <- last_state   # Start at last state
      
      # Get a new list of random step sizes
      more_random_deviations <- rnorm(iter,0,1)*(drift.Coeff*sqDT)
      more_motion <- matrix(more_random_deviations,nrow=iter,ncol=2)
      more_steps <- more_motion+(mu*dt)
      
      for(t in 2:iter){
        this.step <- more_steps[t]
        this.state <- state[t-1,a]+this.step
        state[t,a] <- this.state
        
        if(this.state >= boundary | this.state <= 0){
          added_iterations <- iter*whileLoopNo
          finalT[a] <- (t+added_iterations)+(ndt/dt)   #Total no. of iterations required on each trial
          break
        }
      }
      
      not.finished <- is.na(finalT[a])  # Re-evaluate
      whileLoopNo <- whileLoopNo + 1    # Register while loop iteration
    }
    
    if(this.state > boundary){ # Once the boundary has been passed...
      state[t,a] <- boundary
    }else{
      state[t,a] <- 0
    }
  }
  
  finalT <- finalT*dt
  output <- list(state,finalT)
  names(output) <- c("state","RT")
  return(output)
}


# Final Function: Generate data from this method
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

ddm_sim.randomWalk <- function(trials, boundary, beta,
                               mu=NA, mean.Mu = NA, sd.Mu = NA, 
                               ndt=0.1, drift.Coeff=1, dt=0.0015){
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  #               Defensive Coding                                         #
        noMu <- is.na(mu)
        noDistribution <- is.na(mean.Mu) | is.na(sd.Mu)
        if(noDistribution&noMu){
               stop("Provide a drift rate (mu) OR define its distribution (mean.Mu and sd.Mu)", call. = FALSE)
          }
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
        
  # Get full Random walk using the function in the *customFunctions.R* file
  randomWalk <-  ddm.randomWalk(trials=trials, beta=beta, boundary=boundary, 
                                 mu=mu, mean.Mu = mean.Mu, sd.Mu = sd.Mu, 
                                 ndt= ndt, drift.Coeff= drift.Coeff, dt=dt)
  # Isolate important variables
  RT <- randomWalk$RT
  randomWalk <- randomWalk$state
  
  # Isolate coordinates for last choice
  choices <- getFinalState(randomWalk)
  
  data <- as.data.frame(cbind(choices,RT))
  colnames(data) <- c("Choice","RT")
  
  output <- list("random.walk" = randomWalk,
                 "bivariate.data" = data)
  return(output)
}