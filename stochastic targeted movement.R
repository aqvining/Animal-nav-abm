seek = function (location, target, d = 1, sdAngle = (pi/8)) {
  #Moves an agent at location 'start' to location 'target', with step length d and angular deviations 
  #from directed movement drawn from a normal distrubution with mean 0 and sd = sdAngle
  #input: location = dataframe with 2 columns, x and y, with 1 row giving the starting coordinates
  #       target = vector with 2 values, x and y coodinates for target location
  #       d = float
  #       sdAngle = float, units in radians
  #ouput: dataframe giving x and y coorinates for each time step.
  
  if(sqrt(sum((location[nrow(location),] - target) ^ 2)) <= d) {   #end condition, if within one movement of target . . . (uses pyth thrm)
    return(rbind(location, target))
  } else {
    destDist = as.numeric(as.vector(target - location[nrow(location), ])) #vector with dist in x and dist in y
    destAngle = atan2(destDist[2], destDist[1])                             #angle to target from last location
    deviation = rnorm(1, mean = 0, sd = sdAngle)                          
    turn = destAngle + deviation
    newX = location$x[length(location$x)] + cos(turn) * d
    newY = location$y[length(location$x)] + sin(turn) * d
    location = rbind(location, data.frame(x = newX, y = newY))
    return(seek(location, target, d, sdAngle))                            #recursive call
  }
}

path2geom_segment = function(path) {
  #takes a path from seek function and formats for use with geom_segment plotting in ggplot2
  path$newX = NA
  path$newY = NA
  path$newX[1:(nrow(path) -1)] = path$x[2:nrow(path)]
  path$newY[1:(nrow(path) -1)] = path$y[2:nrow(path)]
  path[1:(nrow(path) -1),]
}

#######sample script################
library(ggplot2)
location = data.frame(x = 0, y = 0)
target = c(-20, -7)
path = seek(location, target, sdAngle = (pi/6))
path_ggplot = path2geom_segment(path)
#plot(x = path$x, y = path$y, col = heat.colors(nrow(path)))
ggplot(path_ggplot, aes(x, y)) + geom_segment(aes(xend = newX, yend = newY))
