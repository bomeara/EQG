OU.sim.branch <- function( sigma=1, alpha=0, theta0=0, theta1=50, alpha2=0,
theta2=-50, ngens=200, nsplit=100, nsims=30, width=10, yylim=c(-60, 60))
{
  quartz( title=paste("OU simulation with sigma = ", sigma), width=9,height=6)   #open new X11 window
# Change quartz to x11 if on a PC or linux
#x11(display=paste("OU simulation with sigma = ", sigma),

  par(pty="s")                        # square plots
  par(mfrow=c(1,2))             # one row, two columns
  par(oma=c(0,0,2,0))   # outer margins, some space for overall title

  devs =rnorm(ngens)

  plot(1:length(devs), sigma*devs, type = "n", col="red", ylim = yylim, xlab="Time", ylab="Value") 
  mtext(c(paste("OU simulation with sigma =", sigma)),3,outer=T)

  x <- 1:ngens*0        # initialize vector for lineage 1
  x2 <- nsplit:ngens*0  # initialize vector for lineage 2 (just from split)
  xfin <- cbind(1:nsims*0,1:nsims*0)  # matrix to store final values for each simulation

  for(h in 1:nsims)      # BM simulation, plots simulations, saves end values
  {                      # h loops over the number of simulations

    devs =rnorm(ngens)
    devs2 =rnorm(ngens-nsplit)

    for(i in 1:(nsplit-1))           ## Before speciation (t=1 to t=99)
    {
      x[(i+1)] = alpha*(theta0-x[i])+x[i] + sigma*devs[i]
      lines(i:(i+1), x[i:(i+1)], col="red")    
    }

    x[(nsplit+1)] = x[nsplit] + sigma*devs[nsplit]     ##SPECIES 1 AT speciation (t=100 to t=101)
    x2[(nsplit+1-nsplit)] = x[nsplit] + sigma*devs2[1] ##SPECIES 2

    lines(nsplit:(nsplit+1), x[nsplit:(nsplit+1)], col="red")   ## plot SPECIES 1
    lines(nsplit:(nsplit+1), c(x[nsplit], x2[1]), col="red")    ## plot SPECIES 2

    for(i in (nsplit+1):(ngens-1))                  ## After speciation   (t=101 to 200)
    {
      x[(i+1)] = alpha*(theta1-x[i])+x[i] + sigma*devs[i]
      j = i-nsplit                ## because second "species" starts at 100 (after the split),
                                  ## but vector i starts at 1, this makes it less confusing
      x2[(j+1)] = alpha*(theta2-x2[j])+x2[j] + sigma*devs2[j]

      lines(i:(i+1), x[i:(i+1)], col="red")
      lines((i):(i+1), x2[j:(j+1)], col="blue")
    }

    xfin[h,1] = x[ngens]                 ## save the end results for x & x2 in matrix xfin
    xfin[h,2] = x2[ngens-nsplit]
  }

## plot the marginal distributions to the RHS
## using the values calculated from the histogram $density
## we could be more sophisticated and either plot these sideways with barplot
## or with the density function, I'll leave that as an excercise

  histdens<-hist(c(xfin[,1],xfin[,2]) , plot=F,nclass=ngens/20)    #the histogram output

## the rotated 'on their side' marginal plots
  plot(histdens$density,histdens$mids,type="l",ylim =yylim,
xlab="Density", ylab="Value")


}

######  The effect of varying sigma   #########

#print("Variance should be sigma^2 * t")
#OU.sim.branch(yylim=c(-100, 100))             # sigma=1
#OU.sim.branch(sigma=2, yylim=c(-100, 100))    # sigma=2
#OU.sim.branch(sigma=3, yylim=c(-100, 100))    # sigma=3 
#OU.sim.branch(sigma=2, alpha=.01, nsims=30, yylim=c(-75,75))