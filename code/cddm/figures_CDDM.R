############## Load custom functions
source("../general_functions/customFunctions.R")
source("../cddm/sim_randomWalk.R")
############## Set up environment
set.seed(123)
library(grid)
library(shape)
library(geostats)
############## Specify par. values to generate data
trials = 500
boundary = 3
mu1 = 2.5      #   Specify drift vector
mu2 = -2.75    #  (Cartesian coordinates)
polar.coordinates = rectToPolar(mu1,mu2)  # Derive polar
dangle <- polar.coordinates[1]  # Drift angle
if(dangle<=0){ dangle <- (2*pi)+dangle }  # Counterclockwise
dlength <- polar.coordinates[2] # Drift length
############## Generate data if necessary
if(!exists("randomWalk")){
   randomWalk = cddm.randomWalk(trials,mu1,mu2,boundary)
}
############## 
state  <- randomWalk$state
finalT <- randomWalk$RT
choices <- getFinalState(state)
#!!!!!!!Defensive coding !!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
trials <- length(finalT)
polar <- rectToPolar(choices[1,1],choices[1,2])
boundary <- round(polar[,"dLength"],2)
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#


# #############################################################
# Figure 1: Showing parameters in CDDM (single RW)
# #############################################################
# Set up format details
########################################
cex.text <- 1
cex.greek <- 1.8
f = 1
# Set up margins
########################################
par(pty="s")             # Square canvas
par(mfrow=c(1,1),        # A single plot
    mar = c(0, 0, 0, 0)) # outer margins
pm <- boundary + 0.2     # Set x/y lims
########################################
# Draw base circle
########################################
# Define points to draw circumference
all.Angles <- seq(0,2*pi,0.001) 
circle <- polarToRect(all.Angles,boundary)
# Create blank plotting space
plot(-10:10,-10:10,type="n", ann = FALSE, axes = FALSE,
     xlim=c(-pm,pm),ylim=c(-pm,pm))      
# Draw the circle
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
########################################
# Draw RW
########################################
# Select arbitrary number of RW to show
show.trials <- 50
# Draw RW
color <- "#EEDB1C"
for(i in 1:show.trials){
  z = seq(1,sum(!is.na(state[,,i])),length.out=75)
  points(state[z,,i], type = "l", col=rgb(0.6,0.5,0.1,0.09), lwd=2)
}
########################################
# Mark response observed
########################################
dot.color <- "#D5B31A"
text.color <- "#9B8D0D"
for(i in 1:show.trials){
  points(choices[i,1],choices[i,2], type = "p", pch =16, cex=2,
         col=rgb(0.65,0.5,0.15,0.2))
}
show.rw = 50
text(choices[show.rw,1]+1,choices[show.rw,2]-0.05,"Responses",
     cex=cex.text, col=text.color, f=f)
text(choices[show.rw,1]+1,choices[show.rw,2]-0.25,"observed",
     cex=cex.text, col=text.color, f=f)
########################################
# Boundary radius
########################################
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
########################################
# Drift vector
########################################
# Identify path prescribed by drift vector
rect <- as.numeric(polarToRect(dangle,boundary))
X <- c(0,rect[1])
Y <- c(0,rect[2])
# Sequence from 0 to driftangle
draw.angle = polarToRect(seq(0,as.numeric(dangle),0.01), 0.8)
# Draw drift vector path
lines(X,Y, lwd=1, col="navy", lty=2)
############################## Drift angle
color.line <- "#D091E6"
color.text <- "#AE5BCB"
points(draw.angle[,1],draw.angle[,2], type="l", 
       col=color.line, lwd=3, lty=3)
text(.8,0.7,expression(theta), cex=cex.greek, col=color.text, f=f)
############################## MuX
color1 <- "#4C8DED"
color2 <- "#1964D3"
lines(c(0,tail(draw.angle,1)[1]),
      c(tail(draw.angle,1)[2],tail(draw.angle,1)[2]),
      lty=3, col=color1, lwd=2)
text(.35,tail(draw.angle,1)[2]-0.2,expression(paste(mu)), cex=cex.greek, col=color2, f=f)
text(.56,tail(draw.angle,1)[2]-0.3,"x", cex=cex.greek-0.75, col=color2, f=f)
############################## MuY
lines(c(tail(draw.angle,1)[1],tail(draw.angle,1)[1]),c(0,tail(draw.angle,1)[2]),
      lty=4, col=color1, lwd=2)
text(0.75,-0.5,expression(paste(mu)), cex=cex.greek, col=color2, f=f)
text(0.95,-0.62,"y", cex=cex.greek-0.75, col=color2, f=f)
############################## Drift length
tmp <- tail(draw.angle,1)[1]
X2 <- seq(0,as.numeric(tmp),0.01)
Y2 <- (mu2/mu1)*X2
color.line <- "#906BE0"
color.text <- "#804BF1"
lines(X2,Y2, lwd=3, col=color.line)
text(0.33,-0.02,expression(delta), cex=cex.greek, col=color.text, f=f)
############################## MuX x MuY
points(tail(draw.angle,1)[1],tail(draw.angle,1)[2], 
       pch=16, col=color2)