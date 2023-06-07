trials = 500
boundary = 3
mu1 = 2.5      #   Specify drift vector
mu2 = -2.75    #  (Cartesian coordinates)
polar.coordinates = rectToPolar(mu1,mu2)  # Derive polar
dangle <- polar.coordinates[1]  # Drift angle
if(dangle<=0){ dangle <- (2*pi)+dangle }  # Counterclockwise
dlength <- polar.coordinates[2] # Drift length