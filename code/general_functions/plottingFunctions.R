###############################################################################
###############################################################################
#####   A set of functions made to plot the data generated from simulations
###############################################################################
########################################################   by Adriana F. Chavez   
############## Load custom functions
source("../general_functions/customFunctions.R")


###############################################################################
###################   Plotting whole random walk
###############################################################################
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function specific to CDDM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotRW.CDDM <- function(randomWalk.output){
  state  <- randomWalk.output$random.walk
  finalState <- getFinalState(state)
  
  bivariate.data <- randomWalk.output$bivariate.data
  finalT <- bivariate.data$RT
  
  trials <- length(finalT)
  polar <- rectToPolar(finalState[1,1],finalState[1,2])
  boundary <- round(polar[,"dLength"],2)
  
  # Formatting and plot settings
  cex.text <- 1
  par(pty="s")             # Square canvas
  par(mfrow=c(1,1),        # A single plot
      mar = c(0, 0, 0, 0)) # outer margins
  pm <- boundary + 0.2     # Set x/y lims
  
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
  text(3.2,0.15,"0", cex=cex.text, f=1, col="black")                # Pi markers
  text(3.25,-0.15,expression(2*pi), cex=cex.text+0.1, f=1, col="black")
  text(-3.15,0.15,expression(pi), cex=cex.text+0.1, f=1, col="black")
  text(-0.25,3.15,expression(pi/2), cex=cex.text+0.1, f=1, col="black")
  text(-0.27,-3.15,expression(3*pi/2), cex=cex.text+0.1, f=1, col="black")
  # Draw RW
  color <- "#EEDB1C"
  for(i in 1:trials){
    z = seq(1,sum(!is.na(state[,,i])),length.out=75)
    points(state[z,,i], type = "l", col=rgb(0.6,0.5,0.1,0.075), lwd=2)
  }
  # Mark response observed
  dot.color <- "#D5B31A"
  z <- 40
  if(trials>z){
    factor <- trials/40}
  else{
    factor <- 5}
  for(i in 1:trials){
    points(finalState[i,1],finalState[i,2], type = "p", pch =16, cex=2,
           col=rgb(0.65,0.5,0.15,1/factor))
  }  
}

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function specific to DDM
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Overall plotting function
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plotRW.output <- (randomWalk.output){
    bivariate.data <- randomWalk.output$bivariate.data
    choice <- bivariate.data$Choice
    binary <- length(table(choice))==2
    
    if(binary){
      
    }else{
      plotRW.CDDM(randomWalk.output)
    }
}


###############################################################################
###################   Plotting whole random walk
###############################################################################