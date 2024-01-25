###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the DDM using the
#####                   RANDOM-WALK EMULATION METHOD 
###############################################################################
########################################################   by Adriana F. Ch?vez   
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


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate the full random walk across many trials (for each trial, 
# keeps the full chain of coordinates visited and response times)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

# Base Function: Simulate the random walk for a single trial
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
ddm.randomWalk <- function(drift, bias, boundary, ndt=0.1, 
                           drift.Coeff=1, dt=0.00015){
    # Define key variables
    sqDT <- sqrt(dt)
    iter <- round(15/dt)        # Maximum number of iterations on the random walk 
    state <- rep(NA, iter)      # Empty vector to store states
    # Arrays to be used in simulation
    motion <- rnorm(iter,0,1)*(drift.Coeff*sqDT)   # Deviations from the drift
    steps <- motion+(drift*dt)
    # Set initial state for every trial
    state[1] <- bias*boundary   # Set initial point for every random-walk on each trial
    
    for(t in 2:iter){
        this.step <- steps[t]
        this.state <- state[t-1]+this.step
        state[t] <- state[t-1]+this.step
        # Stop random-walk if boundary is passed
        if(this.state >= boundary | this.state <= 0){
            finalT <- t+(ndt/dt)   #Total no. of iterations required on each trial
            break
        }
    }
    # Test whether the random-walk reached the boundary, and re-sample if not.
    not.finished <- is.na(finalT)
    whileLoopNo <- 1
    while(not.finished){
          last_state <- state[t]     # Store last state
          state    <- rep(NA, iter)     # Reset random-walk
          state[1] <- last_state   # Start at last state
          # Get a new list of random step sizes
          more_motion <- rnorm(iter,0,1)*(drift.Coeff*sqDT)
          more_steps <- more_motion+(drift*dt)
  
          for(t in 2:iter){
            this.step <- more_steps[t]
            this.state <- state[t-1]+this.step
            state[t] <- this.state
            if(this.state >= boundary | this.state <= 0){
              added_iterations <- iter*whileLoopNo
              #Total no. of iterations required on each trial
              finalT <- (t+added_iterations)+(ndt/dt)   
              break
            }
          }
          not.finished <- is.na(finalT)  # Re-evaluate
          whileLoopNo <- whileLoopNo + 1    # Register while loop iteration
    }
    # Determine whether upper or lower boundary was passed
    if(this.state > boundary){ 
      state[t] <- boundary
    }else{  state[t] <- 0   }
    # Scale Response Time
    finalT <- finalT*dt
    output <- list("state" = state, "RT" = finalT)
  return(output)
}

# Auxiliary Function: Extract final state/choice
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
getFinalState <- function(states.vector){
    randomWalk <- states.vector
    I <- length(randomWalk)
    iters.not.used <- is.na(randomWalk)
    first.NA <- (which(iters.not.used==TRUE)[1])
    choice <- randomWalk[first.NA-1]
    return(choice)
}


# Final Function: Generate data from this method
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
ddm_sim.randomWalk <- function(trials, boundary, bias,
                               drift=NA, mean.drift = NA, sd.drift = NA, 
                               ndt=0.1, drift.Coeff=1, dt=0.0015){
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
    #               Defensive Coding                                         #
    if(is.na(drift)){
      if(is.na(mean.drift)){
        stop("Provide a drift rate (mu) OR define its distribution (mean.Mu and sd.Mu)", call. = FALSE)
      }else{
        drift <- rnorm(trials,mean.drift,sd.drift)
      }
    }
    #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  
    # Run first trial to use as reference point
    first.trial <- ddm.randomWalk(drift[1], bias, boundary)
    choice <- getFinalState(first.trial$state)
    # Define objects where we will store the rest of the trials  
    states <- matrix(NA, nrow = length(first.trial$state), ncol=trials)
    choices <- c(choice, rep(NA, trials-1))
    RT      <- c(first.trial$RT,rep(NA,trials-1))
    # Run the remaining trials
    for(t in 2:trials){
        x <- ddm.randomWalk(drift[t], bias, boundary)
        # Isolate important variables
        RT[t] <- x$RT
        states[,t] <- x$state
        # Isolate coordinates for last choice
        choices[t] <- getFinalState(x$state)
    }
    # Prepare output
    data <- data.frame("Choice" = choices, "RT" = RT)
    output <- list("random.walk" = states, "bivariate.data" = data)
  return(output)
}


#################
# Test/Examples
#################
# Test function
if(!exists("test")){  test <- TRUE     }
if(test){
      drift = 0.5
      boundary = 1
      bias = boundary/2
      mean.drift = 0.5
      sd.drift = 0.1
      x <- ddm.randomWalk(drift, bias, boundary)
      trials = 100
      data <- ddm_sim.randomWalk(trials, boundary, bias,
                                 mean.drift = mean.drift, sd.drift = sd.drift)
}