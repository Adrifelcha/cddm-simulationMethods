###############################################################################
###############################################################################
#####       Customary functions to compute, plot and compare CDFs (eCDF)
###############################################################################
########################################################   by Adriana F. Chavez  
library("mnormt")


###############################################################################
###################   Auxiliary functions for Sampling algorithms
###############################################################################
testBins.bvnormal = function(bin.X,bin.Y,Mean,Sigma){
  dimm = 1:2
  n = min(length(bin.X),length(bin.Y))
  bin.vertices <- cbind(bin.X[1:n],bin.Y[1:n])
  d <- dmnorm(bin.vertices,Mean,Sigma)
  b <- which(d>0)[1]
  if(!is.na(b)){
      if(!b!=1){
          start.at <- NA
          for(i in dimm){
          fix.this.bin <- bin.vertices[,i]
          let.bin.vary <- bin.vertices[,-i]
          unique.densities <- 2
          a <- b
              while(unique.densities>1|a>0){
                fix = rep(fix.this.bin[a],a)
                bin.vertices <- matrix(NA,ncol=2,nrow=a)
                bin.vertices[,i] <- fix
                bin.vertices[,-i] <- let.bin.vary[a:1]
                d.try = dmnorm(bin.vertices,Mean,Sigma)
                unique.densities = length(unique(d.try))
                a = a-1 
              }
          start.at[i] <- a+1
          }
          bin.X <- bin.X[start.at[1]:n]
          bin.Y <- bin.Y[start.at[2]:n]
     }
  }else{
    bin.X <- NA
    bin.Y <- NA
  }
  output = list("bin.X" = bin.X, "bin.Y" = bin.Y)
  return(output)
}