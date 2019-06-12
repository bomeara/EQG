########################################################################
#  Stochastic Simulations 
#  BM and OU simulations to help visualize the effects of the model 
#  parameters on the evolution of species in macroevolution.
#  In phylogenetic comparative methods, we use the BM and OU models 
#  to account for shared ancestry. Therefore, these simulations show a 
#  "cartoon" of the process of evolution that we are assuming
########################################################################
 
nsteps = 100               	# Number of evolutionary steps (generations, etc.)
x <- c(0:nsteps)			# initialize the evolutionary history, 100 gens
devs = rnorm(nsteps)		# deviations = amount of phenotypic change per step
							# simulated here with mean zero variance 1

# x[2] <- x[1] + devs[1]    # evolution = taking the ancestral value and adding  							# on a random amount of change

# More generally, x[i+1] <- x[i] + devs[i]

for (i in 1:nsteps)     # Put into a for loop to traverse the history
{
	x[i+1] <- x[i] + devs[i]
}
# Plot of phenotypic value of one lineage through time
plot(1:length(devs), devs, type = "l", col="red", ylim = c(-max(devs)*30, max(devs)*30), xlab="Time", ylab="Value", main="BM Simulation")

### Now make a nice simulation by ploting each evolutionary change 

plot(1:length(devs), devs, type = "n", col="red", ylim = c(-max(devs)*30, max(devs)*30), xlab="Time", ylab="Value", main="BM Simulation")

devs = rnorm(nsteps) 
for (i in 1:nsteps)     # Put into a for loop to traverse the history
{
  x[i+1] <- x[i] + devs[i]
  lines(i:(i+1), x[i:(i+1)], col="red")
}

# Now add in sigma parameter 
sigma= 20					# sigma parameter, arbitrarily set to 20

devs = rnorm(nsteps) 
for (i in 1:nsteps)     # Put into a for loop to traverse the history
{
  x[i+1] <- x[i] + sigma*devs[i]
  lines(i:(i+1), x[i:(i+1)], col="red")
}

# Wrap into a function:

bmsim <- function(sigma=1, nlineages=1, nsteps=100, yylim=c(-60,60) )  {
							# initalize the plot, but donʻt plot points
  plot(1:length(devs), devs, type = "n",   col="red", 	ylim = yylim, xlab="Time", ylab="Value", main="BM Simulation")

  for (j in 1:nlineages)    # loop over lineages
  { 
    devs = rnorm(nsteps) 		# draw the history of evolutionary changes
    for (i in 1:nsteps) {

      x[i+1] <- x[i] + sigma*devs[i]
      lines(i:(i+1), x[i:(i+1)], col="red")
	
    }
  }
}

### Save the final value for each lineage into a vector and return

bmsim <- function(sigma=1, nlineages=1, nsteps=100, yylim=c(-60,60) )  {

  final <- 1:nlineages # initialize the vector to hold final values
							# initialize the plot, but donʻt plot points
  plot(1:length(devs), devs, type = "n",   col="red", 	ylim = yylim, xlab="Time", ylab="Value", main="BM Simulation")

  for (j in 1:nlineages)    # loop over lineages
  { 
    devs = rnorm(nsteps) 		# draw the history of evolutionary changes
    for (i in 1:nsteps) {

      x[i+1] <- x[i] + sigma*devs[i]
      lines(i:(i+1), x[i:(i+1)], col="red")
	
    }
    final[j] <- x[nsteps+1] 
  }
  return(final)
}

dat <- bmsim(1,30)    # save the final values to "dat"
quartz()              # make a new plot window
hist(dat)             # plot a histogram of the phenotypic distribution


## The OU process is a generalization of BM, with a term for a deterministic 
## component. It is the "strength of selection" times 
##   "the distance from the optimum"
## Make all of the same modifications above to the OU code

ousim <- function(nsteps=100, sigma=1, alpha=.01, theta=20, yylim=c(-60,60))  {

  plot(1:length(devs), devs, type = "n",  col="red", ylim = yylim, xlab="Time", ylab="Value", main="OU Simulation")

  devs = rnorm(nsteps) 
  for (i in 1:nsteps) {

  	x[i+1] <- x[i] + alpha*(theta-x[i]) + sigma*devs[i]
  	lines(i:(i+1), x[i:(i+1)], col="red")
	
  }
}

# The algorithm for a cummulative sum for an OU process
#  	x[i+1] <- alpha*(theta-x[i]) + sigma*devs[i] + x[i]
