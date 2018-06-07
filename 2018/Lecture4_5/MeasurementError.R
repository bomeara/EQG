rm(list=ls()) # because by this point in the course your workspace is probably polluted
library(OUwie)
library(plyr)

# First set up the model parameters
data(sim.ex)
alpha=rep(2,2)
sigma.sq=rep(.45, 2)
theta0=10
theta=rep(10,2)

# Now make the data
sim.data<-OUwie.sim(tree,trait,simmap.tree=FALSE,scaleHeight=FALSE,
                    alpha=alpha,sigma.sq=sigma.sq,theta0=theta0,theta=theta)

# compute all combinations
mserr.vector <- seq(from=0, to=.25*max(sim.data$X), length.out=10)
models.vector <- c("BM1", "OU1")
mserr.argument.vector <- c("none", "known")

# combine all conditions
all.analysis.combinations <- expand.grid(model=models.vector, mserr = mserr.argument.vector, stringsAsFactors = FALSE)



# Assume no measurement error
sim.data$mserr <- 0

#' Do a single OUwie run for a model and summarize the results
#' @param model Model name, as for OUwie
#' @param mserr "known" or "none" as for OUwie
#' @param new.data OUwie-formatted data.frame
#' @param phy Phylo object with node labels
#' @return Data.frame summarizing this run
ouwie.single.run <- function(model, mserr, new.data, phy) {
  result <- OUwie(phy, new.data, model=model, mserr=mserr)
  result.df <- data.frame(mserr=mserr, alpha=result$solution['alpha',1], sigma.sq=result$solution['sigma.sq',1], theta=result$theta[1,1] , AICc=result$AICc, model=model, stringsAsFactors = FALSE)
  return(result.df)
}

#' Compute everything for a given mserr.value. Note we use the same dataset for each combo
#' @param mserr.value The value of measurement error
#' @param all.combos The data.frame of all analysis combinations
#' @param new.data The simulated data
#' @param phy Phylo object with node labels
#' @return A data.frame of results for each element in the combinations
compute.w.error <- function(mserr.value, all.combos = all.analysis.combinations, new.data = sim.data, phy=tree) {
  new.data$X <- rnorm(length(new.data$X), mean=new.data$X, sd=mserr.value)
  new.data$mserr <- mserr.value
  all.result.df <- data.frame()
  for (i in sequence(nrow(all.combos))) {
    # Note that this is slow: it has to keep adding to the object in memory. Look at plyr or dplyr for faster ways to do this
    all.result.df <- rbind(all.result.df, ouwie.single.run(model=all.combos$model[i], mserr=all.combos$mserr[i],  new.data=new.data, phy=phy))
  }
  rownames(all.result.df)<-NULL
  all.result.df$mserr.value = mserr.value
  all.result.df$delta.AICc = NA
  all.result.df$delta.AICc[which(all.result.df$mserr=="none")] <- all.result.df$AICc[which(all.result.df$mserr=="none")] - min(all.result.df$AICc[which(all.result.df$mserr=="none")])
  all.result.df$delta.AICc[which(all.result.df$mserr=="known")] <- all.result.df$AICc[which(all.result.df$mserr=="known")] - min(all.result.df$AICc[which(all.result.df$mserr=="known")]) # note the danger here: what if I'd forgotten to change a none to a known? A two element loop would probably be safer
  return(all.result.df)
}

# Yay, let's run everything. lapply will give a list of data.frames, then we'll use plyr to bind them. Remember, you can use mclapply in the parallel package instead of lapply to run things in parallel (though on Windows this is problematic at times)
all.results <- plyr::rbind.fill(lapply(mserr.vector, compute.w.error, all.combos = all.analysis.combinations, new.data = sim.data, phy=tree))
all.results$mserr.fraction <- all.results$mserr.value / max(sim.data$X)

# While that's running, talk to your neighbors about what you expect. 
#   How will having measurement error that you ignore affect OU vs BM fitting better?
#   How will having measurement error that you incorporate perfectly?
#   Are the values we're using for measurement error reasonable?
#   What do you expect about other parameter estimates?

# Once the run is done, let's plot

# First a plot of delta AICc for BM1 (which is not the generating model)
plot(x=range(all.results$mserr.fraction), y=range(subset(all.results, model=="BM1")$delta.AICc), xlab="Measurement error (as fraction of max trait value)", ylab="∆AICc of BM1 model", type="n", bty="n")
lines(subset(all.results, model=="BM1" & mserr=="known")$mserr.fraction, subset(all.results, model=="BM1" & mserr=="known")$delta.AICc)
lines(subset(all.results, model=="BM1" & mserr=="none")$mserr.fraction, subset(all.results, model=="BM1" & mserr=="none")$delta.AICc, col="red")
legend("topleft", legend=c("With known measurement error", "Assuming zero measurement error"), fill=c("black", "red"))

stop("Let's talk about what's happening here. Why does BM1 drop to zero ∆AICc as measurement error increases, even when we know the error? How does assuming no measurement error affect ∆AICc?")

# Next a plot of sigma-squared
plot(x=range(all.results$mserr.fraction), y=range(all.results$sigma.sq), xlab="Measurement error (as fraction of max trait value)", ylab="Sigma-squared estimate", type="n", bty="n")
lines(subset(all.results, model=="BM1" & mserr=="known")$mserr.fraction, subset(all.results, model=="BM1" & mserr=="known")$sigma.sq)
lines(subset(all.results, model=="OU1" & mserr=="known")$mserr.fraction, subset(all.results, model=="OU1" & mserr=="known")$sigma.sq, lty="dotted")
lines(subset(all.results, model=="BM1" & mserr=="none")$mserr.fraction, subset(all.results, model=="BM1" & mserr=="none")$sigma.sq, col="red")
lines(subset(all.results, model=="OU1" & mserr=="none")$mserr.fraction, subset(all.results, model=="OU1" & mserr=="none")$sigma.sq, lty="dotted", col="red")
abline(h=sigma.sq[1], col="blue", lty="dashed")
legend("topleft", legend=c("BM1 with known measurement error", "OU1 with known measurement error", "BM1 with assumed zero measurement error", "OU1 with assumed zero measurement error"), fill=c("black", "black", "red", "red"), lty=c("solid", "dotted"))

stop("What is affecting sigma-squared? Why?")

# A plot of alpha from OU1
plot(x=range(all.results$mserr.fraction), y=range(subset(all.results, model=="OU1")$alpha), xlab="Measurement error (as fraction of max trait value)", ylab="Alpha of OU1 model", type="n", bty="n", log="y")
lines(subset(all.results, model=="OU1" & mserr=="known")$mserr.fraction, subset(all.results, model=="OU1" & mserr=="known")$alpha)
lines(subset(all.results, model=="OU1" & mserr=="none")$mserr.fraction, subset(all.results, model=="OU1" & mserr=="none")$alpha, col="red")
legend("topleft", legend=c("With known measurement error", "Assuming zero measurement error"), fill=c("black", "red"))

stop("What do you see here about the estimates of alpha?")

