require(MASS)
require(rgl)


ou2drgl <- function(sigma=c(10,0,10), alpha=c(.005, 0, .005), theta0=c(-20, -20), theta1_1=c(-20,20), theta1_2=c(20,-20), regimetimes = 50, nn=100, sims=30, lwd=1) 
{
# Open RGL window, set up environment

# scale alpha and sigma by time, set up matrices, flags for later
  sigma <- sigma/sqrt(nn)   
  alpha <- alpha/nn
  sigma <- matrix(c(sigma[1], 0, sigma[2], sigma[3]), byrow=T, nrow=2)
  alpha <- matrix(c(alpha[1], alpha[2], alpha[2], alpha[3]), byrow=T, nrow=2)
  
  split = (regimetimes<nn)     # where does branch occur?
  optima= !(theta1_1[1]==theta0[1] & theta1_1[2]==theta0[2] & theta1_2[1]==theta0[1] & theta1_2[2]==theta0[2])            # are the optima equal to the root?
  nsplit_end <- nn-regimetimes+1
  
# Initialize arrays for storing simulation data: # time steps, two variables, # simulated lineages 
  x0_a <- array(, dim=c(regimetimes, 2, sims))
  x1_a <- array(, dim=c(nsplit_end, 2, sims))
  x2_a <- array(, dim=c(nsplit_end, 2, sims))
  theta1 <- rbind(theta1_1, theta1_2)

# Make simulations
  for (j in 1:sims)
  {
# Call random walks for each branch of phylogeny 
  	x0 <- randwalk( theta0, regimetimes, theta0, sigma, alpha )               # from start to split
    x1 <- randwalk( x0[regimetimes,], nsplit_end, theta1_1, sigma , alpha)    # lineage 1: split to end
    x2 <- randwalk( x0[regimetimes,], nsplit_end, theta1_2, sigma , alpha)    # lineage 2: split to end

    x0_a[,,j] <- x0      # save results for plotting later
    x1_a[,,j] <- x1
    x2_a[,,j] <- x2

   }
   
# Plot simulations   
  setuprgl()
  rgl.viewpoint(0,5, zoom=.75, fov=15); 
  spheres3d(0, theta0[1], theta0[2], col="black", radius=2)
  if (optima) spheres3d(rep((nn), length(theta1)), theta1[,1], theta1[,2], col="orange", radius=2.5)

   for (j in 1:sims)
   {
     for (k in 1:j)
     {
#      rgl.material(alpha=.05)
#       spheres3d(0, theta0[1], theta0[2], col="black", radius=2)    # plot starting point sphere
       if (optima) spheres3d(rep((nn), length(theta1)), theta1[,1], theta1[,2], col="orange", radius=2.5)
       lines3d(1:regimetimes, x0_a[,1, k], x0_a[,2, k], col="black", size=lwd)  # plot to split
       lines3d(regimetimes:(nn), x1_a[,1, k], x1_a[,2, k], col="red", size=lwd)  # plot after split
       if (split) lines3d(regimetimes:(nn),x2_a[,1, k], x2_a[,2, k], col="blue", size=lwd) # plot lineage2
#       rgl.viewpoint(0,5, zoom=.75, fov=15); 
     }
#    rgl.viewpoint(0,5, zoom=.75, fov=15); 
  }
  rgl.material(alpha=1)
  spheres3d(0, theta0[1], theta0[2], col="black", radius=2)
  spheres3d(rep((nn), sims), x1_a[nsplit_end,1, ], x1_a[nsplit_end,2, ], col="red", radius=2)
  spheres3d(rep((nn), sims), x2_a[nsplit_end,1, ], x2_a[nsplit_end,2, ], col="blue", radius=2)
  if (optima) spheres3d(rep((nn), length(theta1)), theta1[,1], theta1[,2], col="orange", radius=2.5)
  rgl.viewpoint(0,5, zoom=.75, fov=15); 
}

setuprgl <- function()
{
clear3d()#
xyz <- matrix(c(-30,30,0, 30, -30, 100), byrow=T, ncol=3)#
lines3d(xyz[,3], xyz[,1], xyz[,2], alpha=0)#
bbox3d(color=c("#333377","black"), emission="#333377", #
          specular="#3333FF", shininess=5, alpha=0.8)
rgl.viewpoint(0,5, zoom=.75, fov=15); 
rgl.viewpoint(0,5, zoom=.75, fov=15);
}

randwalk<- function(startx=c(0,0), nwalks=50, theta=c(0,0), sigma=matrix(c(1,0,0,1), ncol=2), alpha=matrix(0, ncol=2, nrow=2))
{ 
  x <- matrix(startx, ncol=2, nrow=nwalks)	
  db <- matrix(rnorm(2*nwalks), ncol=2)
  for (i in 1:(nwalks-1))
  {
    dx <- alpha%*%(theta-x[i,]) + sigma%*%db[i,]
    x[i+1,] <- x[i,]+t(dx)
  }
  return(x)
}
