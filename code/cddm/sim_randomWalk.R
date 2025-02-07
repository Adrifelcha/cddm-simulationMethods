###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using the
#####                   RANDOM-WALK EMULATION METHOD 
###############################################################################
########################################################   by Adriana F. Ch?vez   



###############################################################################
# Variable dictionary: ########################################################
# mu1 and mu2 - Individual drift rates for the motion on the x and y axes
# drift.Angle - Direction of the drift vector
# drift.Length - Magnitude of the drift vector
# boundary - Boundary (radius)
# ndt - Non decision time
# drift.Coeff - Within-trial variability on the sampling process
# dt - Step size ("delta-t")
# state - rectangular coordinates recorded during the random walk
###############################################################################

###############################################################################
# Main functions: #############################################################
# Simulate the full random walk across many trials (for each trial, 
# keeps the full chain of coordinates visited and response times)
###############################################################################

# Base Function: Simulate the random walk
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
cddm.randomWalk <- function(trials, mu1, mu2, boundary, ndt=0.1, drift.Coeff=1, dt=0.00015){
  sqDT <- sqrt(dt)
  s.init <- c(0,0) 
  iter <- round(15/dt)  # Maximum number of iterations on the random walk 
  state <- array(NA, dim = c(iter, 2, trials))   # States are saved in a 3dimensional array
  finalT <- rep(NA,trials) # Empty vector to store RT (a.k.a. total number of iterations)
  additional_steps_needed <- rep(0,trials)
  
  # Arrays to be used in simulation
  random_deviations <- rnorm(trials*iter*2,0,1)*(drift.Coeff*sqDT)   # Deviations from step sizes mu1, mu2 (Noise)
  motion <- array(random_deviations,dim = c(iter,2,trials))          # Store deviations in array
  steps_d1 <- matrix(motion[,1,]+(mu1*dt), nrow=iter, ncol=trials)
  steps_d2 <- matrix(motion[,2,]+(mu2*dt), nrow=iter, ncol=trials)
  
  # Set initial state for every trial
  state[1,,] <- s.init # Set initial point for every random-walk on each trial
  
  for(a in 1:trials){   
    ### Random walk per trial
    for(t in 2:iter){
      d1 <- steps_d1[t,a]
      d2 <- steps_d2[t,a]
      state[t,,a] <- state[t-1,,a]+c(d1,d2)
      pass <- sqrt(sum(state[t,,a]^2))
      
      # Stop random-walk if boundary is passed
      if(pass >= boundary){
        finalT[a] <- t+(ndt/dt)   #Total no. of iterations required on each trial
        break
      }
    }
    
    # Test whether the random-walk reached the boundary, and re-sample if not.
    not.finished <- is.na(finalT[a])
    if(not.finished){ additional_steps_needed[a] <- 1 }
    
    whileLoopNo <- 1
    while(not.finished){
      last_state <- state[t,,a]   # Store last state
      state[,,a] <- NA            # Reset random-walk
      state[1,,a] <- last_state   # Start at last state
      
      # Get a new list of random step sizes
      more_random_deviations <- rnorm(iter*2,0,1)*(drift.Coeff*sqDT)
      more_motion <- array(more_random_deviations,dim = c(iter,2))
      more_steps_d1 <- more_motion[,1]+(mu1*dt)
      more_steps_d2 <- more_motion[,2]+(mu2*dt)
      
      for(t in 2:iter){
        d1 <- more_steps_d1[t]
        d2 <- more_steps_d2[t]
        state[t,,a] <- state[t-1,,a]+c(d1,d2)
        pass <- sqrt(sum(state[t,,a]^2))
        
        if(pass >= boundary){
          added_iterations <- iter*whileLoopNo
          finalT[a] <- (t+added_iterations)+(ndt/dt)   #Total no. of iterations required on each trial
          break
        }
      }
      
      not.finished <- is.na(finalT[a])  # Re-evaluate
      whileLoopNo <- whileLoopNo + 1    # Register while loop iteration
    }
    
    if(pass > boundary){ # Once the boundary has been passed...
      # Transform the rectangular coordinates of final state into polar coordinates
      get.Polar <- rectToPolar(state[t,1,a],state[t,2,a])
      # Isolate the radians
      get.Radians <- get.Polar[,"dAngle"] %% (2*pi)
      # Identify the exact point touching the circumference
      final.coord <- polarToRect(get.Radians,boundary)
      # Save these coordinate points on the circle
      final.x <- final.coord$x
      final.y <- final.coord$y
      state[t,,a] <- c(final.x,final.y)
    }
  }
  
  finalT <- finalT*dt
  output <- list(state,finalT)
  names(output) <- c("state","RT")
  return(output)
}

# Auxiliary Function: Extract final states
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
getFinalState <- function(randomWalk.states){
  randomWalk <- randomWalk.states
  dimensions <- dim(randomWalk)
  K <- nrow(randomWalk)
  
  if(length(dimensions)>2){
    I <- dimensions[3]
    coord <- matrix(NA, ncol=2,nrow=I)
    for(i in 1:I){
      for(k in 1:K){
        if(!is.na(randomWalk[k,1,i])){
          a <- k
        }else{
          break
        }
      }
      coord[i,] <- randomWalk[a,,i]
    }
    output <- coord
  }else{
    I <- dimensions[2]
    choice <- rep(NA, I)
    for(i in 1:I){
      for(k in 1:K){
        if(!is.na(randomWalk[k,i])){
          a <- k
        }else{
          break
        }
      }
      choice[i] <- randomWalk[a,i]
    }
    output <- choice
  }
  return(output)
}


# Final function: Generate bivariate data + random walk trace
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
sample.RW.cddm <- function(n, par, drift.Coeff=1, dt=0.0015){
  trials <- n
  boundary <- par$boundary
  drift.Angle <- par$theta
  drift.Length <- par$drift
  ndt <- par$tzero
  mu1 <- par$mu1
  mu2 <- par$mu2
  
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  #               Defensive Coding                                         #
  noPolar <- is.null(drift.Angle) & is.null(drift.Length)
  noRect <- is.null(mu1) & is.null(mu2)
  if(noRect){
    if(noPolar){
      stop("Provide Cartesian or Polar coordinates", call. = FALSE)
    }else{
      Mu <- polarToRect(drift.Angle,drift.Length)
      mu1 <- Mu$x
      mu2 <-Mu$y
    }
  }
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  # Get full Random walk using the function in the *customFunctions.R* file
  full.randomWalk <-  cddm.randomWalk(trials=trials,mu1=mu1,mu2=mu2,
                                      boundary=boundary,ndt=ndt,
                                      drift.Coeff=drift.Coeff,dt=dt)
  # Isolate important variables
  RT <- full.randomWalk$RT
  add.Iterations <- full.randomWalk$repeated.Walk
  randomWalk <- full.randomWalk$state
  
  # Isolate coordinates for last choice
  coord <- getFinalState(randomWalk)
  # Convert to radians
  polar <- rectToPolar(coord[,1],coord[,2])
  rad <- polar[,"dAngle"] %% (2*pi)
  radians <- round(rad,4)
  
  data <- as.data.frame(cbind(radians,RT))
  colnames(data) <- c("Choice","RT")
  
  output <- list("random.walk" = randomWalk,
                 "bivariate.data" = data)
  return(output)
}