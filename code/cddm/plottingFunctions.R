###############################################################################
#####                    Customized plotting functions
###############################################################################
########################################################   by Adriana F. Chavez   
if(!exists("superCalled")){superCalled <- FALSE}
if(!superCalled){ 
  source("./sim_randomWalk.R") 
  source("../general_functions/eCDF.R")
}
library(grid)
library(shape)
library(geostats)

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function 1: Plot choices contained in a bivariate data object, on a circle
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.CDDM_choiceData <- function(data,par=NA,choice.col.RGB = c(0.65,0.5,0.15)){
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
    par(pty="s",             # Square canvas
        mfrow=c(1,1),        # A single plot
        mar = c(0, 0, 0, 0)) # outer margins
    pm <- boundary + 0.2     # Set x/y lims
    pi.at <- boundary + 0.3  # Set position of pi indicators 
    
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
    
    if(trials>z){  factor <- trials/40    }else{    factor <- 5    }
    # Draw final choices
    for(i in 1:trials){
      points(finalState[i,1],finalState[i,2], type = "p", pch =16, cex=2,
             col=rgb(rgbCol[1],rgbCol[2],rgbCol[3],1/factor))
    } 
    # Draw full random walk
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



#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# Function 2: Make an illustrative Figure for the CDDM parameters
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot.CDDM_Fig1 <- function(trials=500, cddm.par=NA, return.RW = TRUE){
  #### If not specified, load some generic parameter values
  if(sum(is.na(cddm.par))>0){
    boundary = 3
    tzero = 0.1
    mu1 = 2.5      
    mu2 = -2.75    
    polar.coordinates = rectToPolar(mu1,mu2)
    theta = polar.coordinates[1]
    drift = polar.coordinates[2]
    cddm.par <- list("theta" = polar.coordinates[1],
                     "drift" = polar.coordinates[2],
                     "boundary" = 3, "tzero" = 0.1)
  }else{
    boundary <- cddm.par$boundary
    theta <- cddm.par$theta
    drift <- cddm.par$drift
    cart.coordinates = polarToRect(theta,drift)
    mu1 <- as.numeric(cart.coordinates[1])
    mu2 <- as.numeric(cart.coordinates[2])
  }
  ### Counterclockwise
  if(theta<=0){     theta <- (2*pi)+theta   }  
  ####   Generate data using random walk algorithm, so we can plot paths
  randomWalk = sample.RW.cddm(trials,cddm.par)
  
  #### Load relevant variables from the randomWalk output
  state  <- randomWalk$random.walk
  finalT <- randomWalk$bivariate.data[,2]
  choices <-   getFinalState(state)
  
  #######################
  ###  Make Figure 1  ###
  #######################
  ### Formatting details
  cex.text <- 1       # Size of text
  cex.greek <- 1.8    # Size of greek letters
  f = 1 
  ### Plot margins
    par(pty="s")             # Square canvas
  par(mfrow=c(1,1),        # A single plot
      mar = c(0, 0, 0, 0)) # outer margins
  pm <- boundary + 0.2     # Set x/y lims.
  ### Define points to draw circumference
  all.Angles <- seq(0,2*pi,0.001) 
  circle <- polarToRect(all.Angles,boundary)
  ### Create blank plotting space
  plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
       xlim=c(-pm,pm),ylim=c(-pm,pm))      
  ### Draw the circle
  points(circle[,1],circle[,2], type="l")
  ### Emphasize X and Y coordinates
  abline(h = 0, lty=4, col="gray60")  # X axis
  abline(v = 0, lty=4, col="gray60")  # Y axis
  ### Add "pi" markers
  text(3.2,0.15,"0", cex=cex.text, f=1, col="black")                # Pi markers
  text(3.25,-0.15,expression(2*pi), cex=cex.text+0.1, f=1, col="black")
  text(-3.15,0.15,expression(pi), cex=cex.text+0.1, f=1, col="black")
  text(-0.25,3.15,expression(pi/2), cex=cex.text+0.1, f=1, col="black")
  text(-0.27,-3.15,expression(3*pi/2), cex=cex.text+0.1, f=1, col="black")
  ### Select arbitrary number of paths to show
  show.trials <- 50
  color <- "#EEDB1C"
  for(i in 1:show.trials){
    # Plot 75 states per path, for a better resolution
    z = seq(1,sum(!is.na(state[,,i])),length.out=75)
    points(state[z,,i], type = "l", col=rgb(0.6,0.5,0.1,0.09), lwd=2)
  }
  ### Identify response observed
  dot.color <- "#D5B31A"
  text.color <- "#9B8D0D"
  for(i in 1:show.trials){
    points(choices[i,1],choices[i,2], type = "p", pch =16, cex=2,
           col=rgb(0.65,0.5,0.15,0.2))
  }
  ### Include a "Responses observed" label next to an arbitrary choice
  show.rw = 50
  text(choices[show.rw,1]+1,choices[show.rw,2]-0.05,"Responses",
       cex=cex.text, col=text.color, f=f)
  text(choices[show.rw,1]+1,choices[show.rw,2]-0.25,"observed",
       cex=cex.text, col=text.color, f=f)
  ### Plot / signal the boundary parameter
  arrow.color = "#C3DDC1"
  text.color = "#739070"
  Arrows(-0.05, 0.05, -2.65, 1.3, code = 2, arr.length = 0.2, arr.width = 0.2,
         arr.adj = 0.5, arr.type = "curved", segment = TRUE,
         col = arrow.color, lcol = arrow.color, lty = 2,
         arr.col = arrow.color, lwd = 1, arr.lwd = 2)
  Arrows(-2.65, 1.3, -0.05, 0.05, code = 2, arr.length = 0.2, arr.width = 0.2,
         arr.adj = 0.5, arr.type = "curved", segment = TRUE,
         col = arrow.color, lcol = arrow.color, lty = 2,
         arr.col = arrow.color, lwd = 1, arr.lwd = 2)
  text(-1.8,1.15,expression(eta), cex=cex.greek, col=text.color, f=f)
  ### Plot / signal the drift vector parameters
  rect <- as.numeric(polarToRect(theta,boundary))
  X <- c(0,rect[1])
  Y <- c(0,rect[2])
  draw.angle = polarToRect(seq(0,as.numeric(theta),0.01), 0.8)
  ### Draw drift vector path
  lines(X,Y, lwd=1, col="navy", lty=2)
  ### Drift angle
  color.line <- "#D091E6"
  color.text <- "#AE5BCB"
  points(draw.angle[,1],draw.angle[,2], type="l", 
         col=color.line, lwd=3, lty=3)
  text(.8,0.7,expression(theta), cex=cex.greek, col=color.text, f=f)
  ### MuX
  color1 <- "#4C8DED"
  color2 <- "#1964D3"
  lines(c(0,tail(draw.angle,1)[1]),
        c(tail(draw.angle,1)[2],tail(draw.angle,1)[2]),lty=3, col=color1, lwd=2)
  text(.35,tail(draw.angle,1)[2]-0.2,expression(paste(mu)), cex=cex.greek, col=color2, f=f)
  text(.56,tail(draw.angle,1)[2]-0.3,"x", cex=cex.greek-0.75, col=color2, f=f)
  ### MuY
  lines(c(tail(draw.angle,1)[1],tail(draw.angle,1)[1]),c(0,tail(draw.angle,1)[2]),
        lty=4, col=color1, lwd=2)
  text(0.75,-0.5,expression(paste(mu)), cex=cex.greek, col=color2, f=f)
  text(0.95,-0.62,"y", cex=cex.greek-0.75, col=color2, f=f)
  ### Drift length
  tmp <- tail(draw.angle,1)[1]
  X2 <- seq(0,as.numeric(tmp),0.01)
  Y2 <- (mu2/mu1)*X2
  color.line <- "#906BE0"
  color.text <- "#804BF1"
  lines(X2,Y2, lwd=3, col=color.line)
  text(0.33,-0.02,expression(delta), cex=cex.greek, col=color.text, f=f)
  ### MuX x MuY
  points(tail(draw.angle,1)[1],tail(draw.angle,1)[2], pch=16, col=color2)
  
  if(return.RW){
    return(randomWalk$bivariate.data)
  }
}

plot.CDDM_margECDF <- function(bivariate.data, color){
      rw.choices <- (sort(bivariate.data[,1]))
      rw.rt <- (sort(bivariate.data[,2]))  
      
      par(pty="m", mfrow=c(1,2), mar = c(3, 3, 3, 1)) 
      # Draw the marginal eCDF for Choices
      plot(rw.choices, myECDF(rw.choices), col="white", pch=16, ann=F, axes=F,
           xlim=c(0,2*pi))
      polygon(c(min(rw.choices),max(rw.choices),max(rw.choices),min(rw.choices)),
              c(-0.1,-0.1,1.1,1.1), col = "gray95", border = "gray80")
      points(rw.choices, myECDF(rw.choices), col=color, pch=16, cex=0.8)
      axis(2,seq(0,1,length.out=10),round(seq(0,1,length.out=10),1), las=2)
      axis(1, seq(0,2*pi,length.out=5),
           c("0",expression(paste(0.5,pi)),expression(paste(pi)),
             expression(paste(1.5,pi)), expression(paste(2,pi))))
      mtext("Choices - eCDF", 3, f=2)
      # Draw the marginal eCDF for the RTs
      plot(rw.rt, myECDF(rw.rt), col="white", pch=16, ann=F, axes=F, xlim=c(0,max(rw.rt)))
      polygon(c(min(rw.rt),max(rw.rt),max(rw.rt),min(rw.rt)),
              c(-0.1,-0.1,1.1,1.1), col = "gray95", border = "gray80")
      points(rw.rt, myECDF(rw.rt), col=color, pch=16, cex=0.8)
      axis(2,seq(0,1,length.out=10),round(seq(0,1,length.out=10),1), las=2)
      axis(1, seq(0,max(rw.rt),length.out=7),
           round(seq(0,max(rw.rt),length.out=7),1))
      mtext("Response Times - eCDF", 3, f=2)
}

myRThistogram <- function(rt.vector, color=NA, maxY=NA){
    if(!is.function(color)){
          color <- function(opacity){   rgb(0.2,0.6,0.9,opacity)    }
    }
    if(is.na(maxY)){
           maxY <- max(density(rt.vector)$y)
    }
    par(pty="m", mar = c(3, 3, 3, 0)) 
    hist(rt.vector, main = "Response Times", col = color(0.3), freq = FALSE, axes=F,
         ylim = c(0,maxY))
    lines(density(rt.vector), col = color(1), lwd=3)
    x.axis <- round(seq(min(rt.vector),max(rt.vector),length.out=10),1)
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
  plot.CDDM(data, par, choice.col.RGB = c(0.15,.29,.80))       
  plot.CDDM_Fig1()
}

