### R code from vignette source 'simulationplots.Rnw'

###################################################
### code chunk number 1: simulationplots.Rnw:3-4
###################################################
options(width=80)


###################################################
### code chunk number 2: simulationplots.Rnw:34-35 (eval = FALSE)
###################################################
##   x[i+1] <- x[i] + devs[i]


###################################################
### code chunk number 3: simulationplots.Rnw:39-41
###################################################
nsteps = 100       # number of steps in our simulation
devs =rnorm(nsteps)    # 100 draws from a normal distribution


###################################################
### code chunk number 4: devs
###################################################
x <- c(0:100)  # initialize x (create a vector of the appropriate length)

for (i in 1:nsteps)  # a loop which will run ngens times
{
   x[i+1] <- x[i] + devs[i]   # add random draw to old value of phenotype
}


###################################################
### code chunk number 5: simulationplots.Rnw:55-56
###################################################
x   # take a look at the values of x, this is how x changes each generation


###################################################
### code chunk number 6: simulationplots.Rnw:60-71
###################################################
plot(1:ngens, 1:ngens, type = "n", col="red", 
ylim = c(-max(devs)*30, max(devs)*30), 
xlab="Time", ylab="Value", main="BM Simulation")

x <- c(0:100)

for (i in 1:nsteps) 
{
   x[i+1] <- x[i] + devs[i]
   lines(i:(i+1), x[i:(i+1)], col="red")
}


###################################################
### code chunk number 7: simulationplots.Rnw:75-76 (eval = FALSE)
###################################################
## x[i+1] <- x[i] + sigma * devs[i]


###################################################
### code chunk number 8: simulationplots.Rnw:80-98
###################################################
bm.plot <- function( sigma=1, ngens=100, nlineages=100, yylim=c(-50, 50) )
{
    plot(1:ngens, 1:ngens, type = "n", col="red", 
    ylim = yylim, 
    xlab="Time", ylab="Value", main="BM Simulation")

    for (i in 1:nlineages)       # number of lineages to simulate
    {
        x <- c(0:100)
        devs =rnorm(ngens)    # draws from a normal distribution each gen
        for (i in 1:ngens) 
        {
           x[i+1] <- x[i] + sigma*devs[i]   # BM equation   
       	# step through time, increasing x a little bit each time
           lines(i:(i+1), x[i:(i+1)], col="red")   # plot line segment
        }
    }    
}


###################################################
### code chunk number 9: simulationplots.Rnw:103-104
###################################################
bm.plot()


###################################################
### code chunk number 10: simulationplots.Rnw:110-111
###################################################
bm.plot( sigma=1, ngens=300 ) 


###################################################
### code chunk number 11: simulationplots.Rnw:116-136
###################################################
bm.plot <- function( sigma=1, ngens=100, nlineages=100, yylim=c(-50, 50) ) 
{
    plot(1:ngens, 1:ngens, type = "n", col="red", ylim = yylim,
     xlab="Time", ylab="Value", main="BM Simulation")
     xfinal <- c(1:nlineages)  ### initialze final values
 
     for (j in 1:nlineages) # number of lineages to simulate
     {
         x <- c(0:ngens)
         devs =rnorm(ngens) # draws from a normal distribution, one per generation
	   for (i in 1:ngens) 
	   {
	      x[i+1] <- x[i] + sigma*devs[i]     # BM equation
	      # step through time, increasing x a little bit each time
	      lines(i:(i+1), x[i:(i+1)], col="red")  # plot line segment
	   }
	 xfinal[j] <- x[ngens+1]     ### collect the last value in each lineage
     }
    return(xfinal)     ### return the final values
} 


###################################################
### code chunk number 12: simulationplots.Rnw:144-146 (eval = FALSE)
###################################################
##     sigma=1
##     cumsum(rnorm(nsteps, sd=sigma)) 


###################################################
### code chunk number 13: simulationplots.Rnw:150-152 (eval = FALSE)
###################################################
##     y <- c(0, cumsum(rnorm(nsteps, sd=sigma)) )
##     lines(0:nsteps, y)


###################################################
### code chunk number 14: simulationplots.Rnw:156-158
###################################################
   sigma=1
   lines(0:nsteps, c(0, cumsum(rnorm(nsteps, sd=sigma))))


###################################################
### code chunk number 15: simulationplots.Rnw:162-172
###################################################
bm.plot <- function( sigma=1, nsteps=100, nlineages=100, yylim=c(-50, 50))
{
# Set up plotting environment
	plot(0, 0, type = "n", xlab = "Time", ylab = "Trait", 
		xlim=c(0, nsteps), ylim=yylim)

# Draw random deviates and plot
	lapply( 1:nlineages, function(x) 
	  lines(0:nsteps, c(0, cumsum(rnorm(nsteps, sd=sigma)))))
}


###################################################
### code chunk number 16: simulationplots.Rnw:210-211 (eval = FALSE)
###################################################
##          x[i+1] = alpha*(theta-x[i])+x[i] + sigma*devs[i]


###################################################
### code chunk number 17: simulationplots.Rnw:216-235
###################################################
ou.plot <- function( alpha=0.005, theta=0, sigma=1, ngens=100, nwalks=30,  yylim=c(-30, 30) )
{
     plot(1:ngens, 1:ngens, type = "n", col="red", ylim = yylim,
     xlab="Time", ylab="Value", main="OU Simulation")
     xfinal <- c(1:nwalks)

     for (j in 1:nwalks)   ### number of lineages to simulate
     {
         x <- c(0:ngens)
         devs =rnorm(ngens)
         for (i in 1:ngens) 
         {
             x[(i+1)] = alpha*(theta-x[i])+x[i] + sigma*devs[i]
             lines(i:(i+1), x[i:(i+1)], col="red") ### plot line segment 
         }
         xfinal[j] <- x[ngens+1]    ### collect last value at each lineage 
     }
     return(xfinal)    ### return the final values
}


###################################################
### code chunk number 18: simulationplots.Rnw:239-240
###################################################
ou.plot( theta=30 )


###################################################
### code chunk number 19: simulationplots.Rnw:268-269 (eval = FALSE)
###################################################
## source("OU.sim.branch.R")


###################################################
### code chunk number 20: simulationplots.Rnw:274-275 (eval = FALSE)
###################################################
## OU.sim.branch


###################################################
### code chunk number 21: simulationplots.Rnw:295-310 (eval = FALSE)
###################################################
## nsteps=100
## nlineages=30
## sigma=1
## sims <- sapply(1:nlineages, function(x) c(0, cumsum(rnorm(nsteps, sd=sigma))))
## yylim <- c(-30, 30)
## 
## png(filename="movies/Rplot%03d.png")  
##               # turn on png graphical device (write to file)
## for (i in 2:nlineages)
## {
## 	plot(0, 0, type = "n", xlab = "Time", ylab = "Trait", 
## 	xlim=c(0, nsteps), ylim=yylim)
## 	apply( sims[,1:i], 2 , function(x) lines(0:nsteps, x, col="red"))
## }
## dev.off()      # turn off png


