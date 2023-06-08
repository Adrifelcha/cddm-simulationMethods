###############################################################################
###############################################################################
#####      A script to simulate bivariate data under the CDDM using the
#####                   RANDOM-WALK EMULATION METHOD 
###############################################################################
########################################################   by Adriana F. Chávez   
source("../customFunctions.R")

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

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Simulate the full random walk across many trials (for each trial, 
# keeps the full chain of coordinates visited and response times)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
cddm_sim.randomWalk <- function(trials, boundary, 
                                 drift.Angle=NA, drift.Length=NA, 
                                 mu1=NA,mu2=NA,
                                 ndt=0.1, drift.Coeff=1, dt=0.0015){
  
  #!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
  #               Defensive Coding                                         #
        noPolar <- is.na(drift.Angle) & is.na(drift.Length)
        noRect <- is.na(mu1) & is.na(mu2)
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
  randomWalk <-  cddm.randomWalk(trials=trials,mu1=mu1,mu2=mu2,
                                 boundary=boundary,ndt=ndt,
                                 drift.Coeff=drift.Coeff,dt=dt)
  # Isolate important variables
  RT <- randomWalk$RT
  add.Iterations <- randomWalk$repeated.Walk
  randomWalk <- randomWalk$state
  
  # Isolate coordinates for last choice
  coord <- getFinalState(randomWalk)
  # Convert to radians
  polar <- rectToPolar(coord[,1],coord[,2])
  rad <- polar[,"dAngle"] %% (2*pi)
  radians <- round(rad,4)
  
  data <- as.data.frame(cbind(radians,RT))
  colnames(data) <- c("Choice","RT")
  return(data)
}