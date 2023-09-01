###############################################################################
###############################################################################
#####  A set of functions made to plot the CDDM data generated from simulations
###############################################################################
########################################################   by Adriana F. Chavez   
############## Load custom functions
source("../general_functions/customFunctions.R")

plot.CDDM <- function(data, par=NA, 
                      choice.col.RGB = c(0.65,0.5,0.15)){
    randomWalk <- is.null(dim(data))
    params.available <- sum(is.na(par)) == 0
    if(randomWalk){
      state  <- data$random.walk
      finalState <- getFinalState(state)
      polar <- rectToPolar(finalState[1,1],finalState[1,2])
      boundary <- round(polar[,"dLength"],2)
      bivariate.data <- data$bivariate.data
      finalT <- bivariate.data$RT
    }else{
      choice  <- data[,1]
      finalT <- data[,2]
      if(!params.available){
        print("Please specify parameter values used in simulation")
      }else{
        boundary <- par$boundary
        finalState <- polarToRect(choice,boundary)
      }
    }
    
    trials <- length(finalT)
    
    # Formatting and plot settings
    cex.text <- 1
    par(pty="s")             # Square canvas
    par(mfrow=c(1,1),        # A single plot
        mar = c(0, 0, 0, 0)) # outer margins
    pm <- boundary + 0.2     # Set x/y lims
    pi.at <- boundary + 0.3     # Set position of pi indicators 
    
    # Create blank plotting space
    plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
         xlim=c(-pm,pm),ylim=c(-pm,pm))      
    # Draw base circle
    all.Angles <- seq(0,2*pi,0.001) 
    circle <- polarToRect(all.Angles,boundary)
    points(circle[,1],circle[,2], type="l")
    # Emphasize X and Y coordinates
    abline(h = 0, lty=4, col="gray60")  # X axis
    abline(v = 0, lty=4, col="gray60")  # Y axis
    # Add "pi" markers
    text(pi.at,0.15,"0", cex=cex.text, f=1, col="black")                # Pi markers
    text(pi.at,-0.15,expression(2*pi), cex=cex.text+0.1, f=1, col="black")
    text(-pi.at,0.15,expression(pi), cex=cex.text+0.1, f=1, col="black")
    text(-0.25,pi.at,expression(pi/2), cex=cex.text+0.1, f=1, col="black")
    text(-0.27,-pi.at,expression(3*pi/2), cex=cex.text+0.1, f=1, col="black")
    # Mark response observed
    z <- 40
    rgbCol = as.numeric(choice.col.RGB)
    if(trials>z){
      factor <- trials/40
      }else{
      factor <- 5}
    for(i in 1:trials){
      points(finalState[i,1],finalState[i,2], type = "p", pch =16, cex=2,
             col=rgb(rgbCol[1],rgbCol[2],rgbCol[3],1/factor))
    } 
    # Draw RW
    if(randomWalk){
      max.trials.plot = min(c(trials,200))
      color <- "#EEDB1C"
      for(i in 1:max.trials.plot){
        z = round(seq(1,sum(!is.na(state[,1,i])),length.out=75),0)
        points(state[z,,i], type = "l", lwd=2,
               col=rgb(rgbCol[1],rgbCol[2],rgbCol[3],50/max.trials.plot))
      }
    }
}

# Test function
if(!exists("test")){    test <- TRUE                           }
if(test){
  par <- list("drift" = 1, 
              "theta" = pi,
              "tzero" = 0.1,
              "boundary" = 7)
  n <- 5000
  C <- runif(n,0,2*pi)
  RT <- rexp(n,3)
  data <- cbind(C,RT)
  plot.CDDM(data, par, choice.col.RGB = c(0.15,.29,.80))       }

