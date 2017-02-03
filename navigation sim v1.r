##Version Updates
#Version 1: everything is new!
#Someday, I will describe the model here

###Notes for Alexander###
##Things that need fixing:
#~~noticed that detour sometimes goes to a resource further away than AG destinations. This should not be happening.

##Things to improve:
#~~Additive gravity rule (weight of j vs. G)
#~~Equation for when primates leave a patch
#~~Clustering and other spatial statistics of resources
#~~Magic numbers removed and added as parameters (extraction_rate)

##Things to add:
#~~resource birth/death
#~~new rules for starting locations
#~~limited LTM, LTM updating by sight
#~~increase functionality with multiple primates
#~~rank order behaviors and other group behavior rules

##Things to consider
#~~Detour will go to any resource closer than AG, consider requiring detour to move toward AG
#~~implementation of NNR could be informative, could also unnecessarily slow processing time

#####################################SET-UP##################################################
library(plyr)
library(calibrate)
library(graphics)

defaultEnergy = 0				#energy that primates start with by default. Can be overwritten during primate creation
fieldSize = 10				#maximum absolute value for x and y coordinates of a resource. Model currently only supports square fields
resourceAvail = 600			#average total resource availability. Expected mean resource size will be this value divided by resource num
resSizeRange = 7				#size by which any resource is able to deviate from the expected mean. Resource size is a random uniform draw.
resourceNum = 20				#number of resources generated
cluster = 0.2				#probability that a new resource will be created between 1 and 2 units away from the previous resource
depletionIndex = 0			#scales the diminishing return factor of resources. At 0, resource density has no impact on foraging return. At 1, foraging return is multiplied by the density
k = 1						#distance exponent in gravity model. The simple gravity score is Size/ Distance ^ k
primateSpeed = 1				#how many units a primate can move per step
turnDev = pi/4          #standard deviation of turn angle (mean is toward target dest)
groupSizes = c(1)				#each element of the vector represents a group, the value is the size of the group. Currently, all individuals are created on the same point, making more than one individual indistinguisable from increasing the foraging rate of a single individual. plotEnergy function does not handle multiple primates yet
moveCost = 2				#energy decrease each step a primate moves
metabCost = 2				#energy decrease each step. If primate moves, this is added to moveCost. Also currently used to determine if a primate should leave current resource (primate leaves when expected foraging return is less than this value)
rRegen = c(min = 0, max = 0)		#determines density increase of resources for each step. Resources are assigned a value uniformly drawn from this range during environment creation. To scale regen rate to resource value, manually change in fillRegen_rate function
types = c("AG")		#movement models to simulate across all created environments. Model can handle vector of any length, but currently only has methods for "AG" and "hybrid"

P = list(defEn = defaultEnergy, fieldSize = fieldSize, resAvail = resourceAvail, resSR = resSizeRange, resNum = resourceNum, clus = cluster, depI = depletionIndex, k = k, pSpeed = 1, turnDev = turnDev, grSizes = groupSizes, mCost = moveCost, metabCost = metabCost, rRegen = rRegen, types = types)  #saves all parameters to a list that can easily be passed to functions

###################################FUNCTIONS####################################################

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Environment Creation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
createEnviron = function(P) {
#input: List of parameters
#output: data frame with number of rows = resourceNum and capital letters for row names. 
#data frame created contains all resource information including coordinates, return_value, extract_rate, density (starts at 1), and regen_rate. If more than 26 resources, row-names will become NA which will affect some functions in model
	environ = data.frame(matrix(0, nrow = P$resNum, ncol = 5, dimnames = list(LETTERS[1:P$resNum], c("coordinates", "return_value", "extract_rate", "density", "regen_rate"))))			
	environ$coordinates = vector("list", length(environ$coordinates))																					#changes coordinates column to a list (allows mapply to fill each element with a vector)
	for (i in seq(nrow(environ))) environ$coordinates[[i]] = fillCoords(environ$coordinates, i, P)[[i]]															
	environ$return_value = sapply(environ$return_value, function(x) round(runif(n = 1, min = P$resAvail/P$resNum - P$resSR, max = P$resAvail/P$resNum + P$resSR), digits = 0))			
	environ$extract_rate = 1																											#If changed from one, add as a parameter
	environ$density = 1																												
	environ$regen_rate = mapply(fillRegen_rate, environ$regen_rate, MoreArgs = list(P = P))																	
	environ 																														
}

fillCoords = function(coords, i, P) {
#input: a list of coordinates (can be empty), index to fill, and list of parameters
#output: new list of coordinates
#only generates one set of coordinates, but takes whole list as input so previous coordinates can be referenced. element of list prior to i should be filled
#note* The output is a whole list of corrdinates, so be sure to index the result if you only want the coordinates that were generated.		
	if (i == 1) {coords[[i]] = c(x = round(runif(n = 1, min = -P$fieldSize, max = P$fieldSize), digits = 2), y = round(runif(n = 1, min = -P$fieldSize, max = P$fieldSize), digits = 2))						#First element of list is generated randomly. Adbsolute value of generated coordinates mus be less than field size
	 } else {
		while(! min(mapply(getDist, c1 = coords[-i], MoreArgs = list(c2 = coords[i])), na.rm = TRUE) > 1 || min(mapply(getDist, c1 = coords[-i], MoreArgs = list(c2 = coords[i])), na.rm = TRUE) == Inf) {  			#Prevents resources from being created too close to others. If the minumum distance of the current coordinate from any other coordinate is < 1(works if current coordinate has not been assigned) assign a new value to the current coordinate (Second condition allows for vector of NA's)
			if (runif(n = 1) >= P$clus) { coords[[i]] = c(x = round(runif(n = 1, min = -P$fieldSize, max = P$fieldSize), digits = 2), y = round(runif(n = 1, min = -P$fieldSize, max = P$fieldSize), digits = 2))		#**NEEDS ATTENTION**runs probability of clustering, and if TRUE new resource x AND y coordinates must be between 1 and 2 units from previous resource. Should limit total distance instead 
			 }else coords[[i]] = c(round(coords[[i-1]]["x"] + 1 + runif(n = 1, min = -1, max = 1), digits = 2), round(coords[[i-1]]["y"] + 1 + runif(n = 1, min = -1, max = 1), digits = 2))
		}	
	}
	coords
}


fillRegen_rate = function(x, P) {
#generates a single regenration rate using parameters. Can probably be integrated with createEnviron
	x = round(runif(n = 1, min = P$rRegen["min"], max = P$rRegen["max"]), digits = 4)
	x
}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Primate Creation~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
primateValidity = function(object) {
  #Checks that input values of a newly created primate are valid. Returns appropriate error messages in cases where they are not
  #NEEDS EXPANDING
  if (! length(object@location) == 2) {return("location contains more than 2 coordinates")
    if (! length(object@energy) == 1) return("energy is not a single value") }
  else TRUE
}

newPrimate = setClass("primate", slots = c(location = "numeric", LTM = "data.frame", energy = "numeric", forageEf = "numeric", speed = "numeric", turnDev = "numeric", type = "character"), validity = primateValidity)  #creates the primate class of objects
#*note there are currently no primate specific methods. This exists as a class mostly to allow expansion of subclasses by movement type or dominance in future versions



buildPrimate = function(location, environ, P, type) {
#input: location create a new primate, environment in which to create primate, model parameters, and type of movement primate should use
#output: an object of class primate
#most slots are filled from the parameters list. If forage efficiency is change, add as a parameter.
#LTM filled from environment, with columns added for gravity and j for location. *NOTE* Future version may implement limited memory size, which will need a more sophistacted LTM creation
	primate = newPrimate(location = location, LTM = environ, energy = P$defEn, forageEf = 0.2, speed = P$pSpeed, turnDev = P$turnDev, type = type)
	primate@LTM$gravity = mapply(calcGravity, coordinates = primate@LTM$coordinates, 
		return_value = primate@LTM$return_value, density = primate@LTM$density, MoreArgs= list(location = primate@location, k = P$k))
	primate@LTM$j = mapply(calc_j, i = 1:nrow(primate@LTM), 
		MoreArgs = list(coords = primate@LTM$coordinates, values = primate@LTM$return_value, dens = primate@LTM$density, k = P$k))
	primate@LTM = transform(primate@LTM, AG = gravity + j/nrow(primate@LTM))
	primate
}

createGroup = function(grSize, environ, P, type) {
#allows groups of primates to be created in seperate lists. This currently has no functionality, but allows for easier addition of group behavior rules in later versions
	group = replicate(grSize, buildPrimate(location = environ$coordinates[[which(environ$return_value == max(environ$return_value))[1]]], environ = environ, P = P, type = type)) #places all primates on most valueable resource
	group
}

	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Calculations and Tools~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

calcForage = function(primate, environ, P) {
	patch = which(environ$coordinates %in% list(primate@location))
	if (environ$density[patch] < primate@forageEf) { 0																				#disallows foraging if patch density is lower than primates foraging efficiency
	} else environ$return_value[patch] * primate@forageEf * environ$extract_rate[patch] * (environ$density[patch]/(P$depI + (1 - P$depI) * environ$density[patch]))
}

calcGravity = function(location, coordinates, return_value, density, k) {
#*note* uses actual coordinates for location rather than primate object because resource coodinates are also used as inout when calculating j.
	distance = getDist(location, coordinates)
	mass = return_value * density
	if (distance == 0) return(NA)
	if (distance < 1) {return(mass)				#moving to resource takes full step, won't divide by less than 1
	} else return(mass /(distance ** k))
}

calc_j = function(i, coords, values, dens, k) {
	sum(mapply(calcGravity, coordinates = coords, return_value = values, density = dens, MoreArgs = list(location = coords[[i]], k = k)), na.rm = TRUE)
}

calcDetour = function(primate, possDetours, AGPatch, P) {
	LTM = primate@LTM[possDetours,]
	target = sampleMod(which(LTM$dist == min(LTM$dist)))
	addRes = LTM$return_value[target] * LTM$density[target]
	addTrav = (P$mCost + P$metabCost) * (getDist(primate@location, LTM$coordinates[[target]]) + getDist(LTM$coordinates[[target]], LTM$coordinates[[which(row.names(LTM) == AGPatch)]]) - getDist(primate@location, LTM$coordinates[[which(row.names(LTM) == AGPatch)]]))
	if (addRes >= addTrav) {return(LTM$coordinates[[target]])
	 } else calcDetour(primate, possDetours[-which(possDetours == row.names(LTM[target,]))], AGPatch, P)				#Recursion! If closest patch not an acceptable detour, try next closest patch
}


getDist = function(c1, c2) {
#input: two sets of coordinates. x and y coordinates must be names
#output: Euclidean distance between input points
	c1 = unlist(c1)
	c2 = unlist(c2)
	if (! is.numeric(c1) || ! is.numeric(c2)) return(NA)
	if (identical(c1, c2)) return(0)
	unname(sqrt(((c1["x"] - c2["x"]) ** 2) + ((c1["y"] - c2["y"]) ** 2))) 		 #pythagorean theorem
}

sampleMod = function(x) {
#modifies the built in sample function so an input of a single value returns that value rather than sampling all integers between zero and input
	if (length(x) == 1 ) return(x) else return(sample(x, size = 1))
}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Step Actions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
runStep = function(prevStep) {
#fully updates data by one step
#*NOTE* returns NA if no valid action available (most commonly meaning all resources are depleted past the point that foraging outways metabolism cost)
	environ = prevStep$environ
	allPrimates = prevStep$allPrimates
	P = prevStep$P
	for (i in seq(length(allPrimates))) {													#execute step action for each primate
		for (j in seq(length(allPrimates[[i]]))) {
			primate = allPrimates[[i]][[j]]
			action = decideAction(primate, environ, P)
			if (action == "move") {primate = move(primate, environ, P)
				primate@energy = primate@energy - P$mCost
			} else { 
				if (action == "forage") {
					forageReturn = calcForage(primate, environ, P)
					primate@energy = primate@energy + forageReturn
					patch = which(environ$coordinates %in% list(primate@location))
					environ$density[patch] = round(environ$density[patch] - primate@forageEf, digits = 3)	#*NOTE* depI and extraction rate are calculated in forage return, but not here for desnity depletion. Consider updating
				} else stop("invalid action")	
			} 
			allPrimates[[i]][[j]] = primate
		}
	}
	for (i in seq(length(allPrimates))) {													#updates each primates LTM and energy after all primates have taken action
		for (j in seq(length(allPrimates[[i]]))) {
			allPrimates[[i]][[j]]@energy = allPrimates[[i]][[j]]@energy - P$metabCost					#this line could be moved to previous set of for loops
			allPrimates[[i]][[j]]@LTM = updateLTM(allPrimates[[i]][[j]], environ)
		}
	}
	environ = transform(environ, density = round(density + regen_rate, digits = 4)) 						#resource regeneration
	environ$density = sapply(environ$density, function(x) if (x > 1) 1 else x)  							#limits density to 1 *NOTE* may also need statment to not allow density below 0 if new rules are implemented for staying at a resource
	list(allPrimates = allPrimates, environ = environ, P = P)
}

decideAction = function(primate, environ, P) {
	if (! list(primate@location) %in% environ$coordinates) return("move")
	if (calcForage(primate, environ, P) < P$metabCost) { return("move") 
	} else return("forage")
}


move = function(primate, environ, P) {
	if (primate@type == "AG") {dest = getDestAG(primate)
	 } else { if (primate@type == "hybrid") {dest = getDestHybrid(primate, P)
	 	} else stop("invalid movement type")
	}
	if (getDist(primate@location, dest) <= 1) {primate@location = dest
	 } else primate@location = updateLoc(primate, dest)
	primate
}

getDestAG = function(primate) {
	dest = primate@LTM$coordinates[[sampleMod(which(primate@LTM$AG == max(primate@LTM$AG, na.rm = TRUE)))]]
}

getDestHybrid = function(primate, P) {
	AGPatch = row.names(primate@LTM[sampleMod(which(primate@LTM$AG == max(primate@LTM$AG, na.rm = TRUE))),])
	primate@LTM$dist = mapply(getDist, c1 = primate@LTM$coordinates, MoreArgs = list(c2 = primate@location))
	rCurLoc = primate@LTM[which(primate@LTM$dist > 0),]  #LTM without resource of current location
	rLowDens = rCurLoc[which(rCurLoc$density > primate@forageEf),]
	possDetours = row.names(rLowDens[which(rLowDens$dist <= rLowDens[AGPatch, "dist"]),])
	dest = calcDetour(primate, possDetours, AGPatch, P)
	dest
}

updateLoc = function(primate, dest) {
	targetAng = unname(atan2(x = dest["x"] - primate@location["x"], y = dest["y"] - primate@location["y"]))
	deviation = rnorm(1, mean = 0, sd = primate@turnDev)
	heading = targetAng + deviation
	moveVector = c(x = cos(heading) * primate@speed, y = sin(heading) * primate@speed)
	primate@location + moveVector
}

updateLTM = function(primate, environ) {
	primate@LTM$density = environ$density
	primate@LTM$gravity = mapply(calcGravity, coordinates = primate@LTM$coordinates, 
		return_value = primate@LTM$return_value, density = primate@LTM$density, MoreArgs= list(location = primate@location, k = P$k))
	primate@LTM$j = mapply(calc_j, i = 1:nrow(primate@LTM), 
		MoreArgs = list(coords = primate@LTM$coordinates, values = primate@LTM$return_value, dens = primate@LTM$density, k = P$k))
	primate@LTM = transform(primate@LTM, AG = gravity + j/nrow(primate@LTM))
	primate@LTM
}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Run Functions~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

initialize = function(P, environ, type) {
	allPrimates = as.list(P$grSizes)
	allPrimates = mapply(createGroup, allPrimates, MoreArgs = list(environ = environ, P = P, type = type), SIMPLIFY = FALSE)
	initialStep = list(allPrimates = allPrimates, environ = environ, P = P)
	initialStep
}

runSimulation = function(P, nsteps) {
	allTypes = vector("list", length(types))
	names(allTypes) = P$types
	initEnviron = createEnviron(P)
	for (type in P$types) {
		allTypes[[type]] = vector("list", nsteps)
		allTypes[[type]][[1]] = initialize(P, initEnviron, type)
		n = 2
		while (n <= nsteps && max(allTypes[[type]][[n - 1]]$environ$density) >= allTypes[[type]][[n - 1]]$allPrimates[[1]][[1]]@forageEf) {
			allTypes[[type]][[n]] = runStep(allTypes[[type]][[n - 1]])
			n = n + 1
		}
	}
	allTypes
}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~Plotting and Data Analysis~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
plotRun = function(run) {
#Plots an entire data set one step at a time. One plot in graphics window for each movement type simulated plus another for trace of energy in each type
#Stops if input is NA, generally indicating the primates has consumed all available resources in a time step.
	par(mfcol = c(1, length(run) + 1))
	for (i in seq(length(run[[1]]))) {
		for (j in seq(length(run))) {
			plotStep(run[[j]][[i]], i)
		}
		plotEnergy(run, i)
		Sys.sleep(1)
	}
}

plotStep = function(step, stepNum) {
	resX = vector("list", nrow(step$environ))
	resX = mapply(function(x, coords) coords["x"], x = resX, coords = step$environ$coordinates)
	resY = vector("list", nrow(step$environ))
	resY = mapply(function(y, coords) coords["y"], y = resY, coords = step$environ$coordinates)

	primX = vector("list", length(unlist(step$allPrimates)))
	primX = mapply(function(x, prim) prim@location["x"], x = primX, prim = unlist(step$allPrimates))
	primY = vector("list", length(unlist(step$allPrimates)))
	primY = mapply(function(y, prim) prim@location["y"], y = primY, prim = unlist(step$allPrimates))

	energies = vector("list", length(unlist(step$allPrimates)))
	energies = mapply(function(x, prim) prim@energy, x = energies, prim = unlist(step$allPrimates))


	plot(x = resX, y = resY, main = paste(step$allPrimates[[1]][[1]]@type, "step", stepNum, sep = " "))
	points(x = primX, y = primY, col = 2, pch = 4)
	textxy(primX, primY, energies, col = 2, cex = 1)
}

plotEnergy = function(run, steps = 0) {
#Does not work for more than one primate per step. Consider using sum(unlist(allPrimates)) to implement this functionality
	if (steps == 0) steps = length(run[[1]])
	if (steps == 1) {plot(x = seq(2), y = c(0, 0), type = "l")
	 } else {
		energies = vector("list", length(run))
		run = lapply(run, function(run) run[1:steps])
		energies = mapply(mapply, step = run, MoreArgs = list(FUN = function(step) step$allPrimates[[1]][[1]]@energy))
		plot(energies[,1], type = "l")
		for (i in 2:length(run)) lines(x = energies[,i], col = i)
	}
}

plotPath = function(run, primate = c(1,1)){
  #plots the path of a single primate from a run object. The "run" input should be only the subset of one "type" of run, eg. "AG".
  #The "primate" input should be a 2 value vector giving the group number and individual number of the target agent
  path = lapply(run, function(x, y) x$allPrimates[[y[1]]][[y[2]]]@location, y = primate) #get location coordinates for all times
  path = data.frame(matrix(unlist(path), nrow = length(path), ncol = 2, byrow = TRUE)) #convert location from list of cordinates to df with x and y columns
  colnames(path) = c("x", "y")
  pathPlot = ggplot(path2geom_segment(path), aes(x, y)) + geom_segment(aes(xend = newX, yend = newY, color = colorRampPalette(c("red", "blue"))(49))) + theme(legend.position = "none")
  envCoord = run[[1]]$environ$coordinates         #get patch coordinates from first time step
  envCoord = data.frame(matrix(unlist(envCoord), nrow = length(envCoord), ncol = 2, byrow = TRUE)) #convert patch locations from list of cordinates to df with x and y columns
  colnames(envCoord) = c("x", "y")
  pathPlot +  geom_point(data = envCoord, mapping = aes(x, y), shape = 1, size = 3)
}

path2geom_segment = function(path) {
  #takes a matrix with rows of x,y coords and and formats for use with geom_segment plotting in ggplot2
  path$newX = NA
  path$newY = NA
  path$newX[1:(nrow(path) -1)] = path$x[2:nrow(path)]
  path$newY[1:(nrow(path) -1)] = path$y[2:nrow(path)]
  path[1:(nrow(path) -1),]
}
	#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


####################################Run Scripts################################################

#~~~~~~Simple~~~~~~~~#
#uncomment following lines to run
nsteps = 50
run = runSimulation(P, nsteps)
plotRun(run)


