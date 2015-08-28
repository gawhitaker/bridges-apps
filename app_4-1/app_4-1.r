##Acceptance probabilities for application 4.1 with option to plot
##Birth-death, fully observed
######################################################################
##Load helper functions and required libraries

#install.packages("deSolve") ##uncomment to install "deSolve" package if required 
library(deSolve)
source("helper_4-1.r")

######################################################################
##Options for Metropolis-Hastings samplers
##x0 = initial value of the process (x0) (left-endpoint) 
##x1 = observation at T (xT) (right-endpoint) 
##inter = T-0
##cvec = parameter vector (theta_1,theta_2)
##gamma = tuning parameter (LB only)
##m = discretisation to use
##it = number of iterations 
##plot = logical for plot
##xl = lower y limit for x in plot
##xu = upper y limit for x in plot

######################################################################
##Metropolis-Hastings samplers for median values for T=1, m=50

##MDB
MDBMH(x0=50,x1=24.615,inter=1,cvec=c(0.1,0.8),m=50,it=100,plot=TRUE,xl=20,xu=55)

##LB
LBMH(x0=50,x1=24.615,inter=1,cvec=c(0.1,0.8),m=50,gamma=0.1,it=100,plot=TRUE,xl=20,xu=55)

##RB 
RBMH(x0=50,x1=24.615,inter=1,cvec=c(0.1,0.8),m=50,it=100,plot=TRUE,xl=20,xu=55)

##RB^- 
RBminusMH(x0=50,x1=24.615,inter=1,cvec=c(0.1,0.8),m=50,it=100,plot=TRUE,xl=20,xu=55)

##GP-N 
GPNMH(x0=50,x1=24.615,inter=1,cvec=c(0.1,0.8),m=50,it=100,plot=TRUE,xl=20,xu=55)

##GP-S 
GPSMH(x0=50,x1=24.615,inter=1,cvec=c(0.1,0.8),m=50,it=100,plot=TRUE,xl=20,xu=55)

##GP 
GPMH(x0=50,x1=24.615,inter=1,cvec=c(0.1,0.8),m=50,it=100,plot=TRUE,xl=20,xu=55)

##GPMDB 
GPMDBMH(x0=50,x1=24.615,inter=1,cvec=c(0.1,0.8),m=50,it=100,plot=TRUE,xl=20,xu=55)
