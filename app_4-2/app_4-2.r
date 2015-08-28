##Acceptance probabilities for application 4.2 with option to plot
##Lotka-Volterra, fully observed
######################################################################
##Load helper functions and required libraries

#install.packages("deSolve") ##uncomment to install "deSolve" package if required 
library(deSolve)
source("helper_4-2.r")

######################################################################
##Options for Metropolis-Hastings samplers
##x0 = vector of initial value of the process (x1,x2) (left-endpoint) 
##x1 = observation at T (x1,x2) (right-endpoint) 
##inter = T-0
##cvec = parameter vector (theta_1,theta_2,theta_3)
##gamma = tuning parameter (LB only)
##m = discretisation to use
##it = number of iterations 
##plot = logical for plot
##x1l = lower y limit for x1 in plot
##x1u = upper y limit for x1 in plot
##x2l = lower y limit for x2 in plot
##x2u = upper y limit for x2 in plot

######################################################################
##Metropolis-Hastings samplers for median values for T=1

##MDB
MDBMH(x0=c(71,79),x1=c(96.816,71.926),inter=1,cvec=c(0.5,0.0025,0.3),m=50,it=100,plot=TRUE,x1l=65,x1u=105,x2l=65,x2u=85)

##LB
LBMH(x0=c(71,79),x1=c(96.816,71.926),inter=1,cvec=c(0.5,0.0025,0.3),m=50,gamma=0.01,it=100,plot=TRUE,x1l=65,x1u=105,x2l=65,x2u=85)

##RB 
RBMH(x0=c(71,79),x1=c(96.816,71.926),inter=1,cvec=c(0.5,0.0025,0.3),m=50,it=100,plot=TRUE,x1l=65,x1u=105,x2l=65,x2u=85)

##RB^-
RBminusMH(x0=c(71,79),x1=c(96.816,71.926),inter=1,cvec=c(0.5,0.0025,0.3),m=50,it=100,plot=TRUE,x1l=65,x1u=105,x2l=65,x2u=85)

##GP 
GPMH(x0=c(71,79),x1=c(96.816,71.926),inter=1,cvec=c(0.5,0.0025,0.3),m=50,it=100,plot=TRUE,x1l=65,x1u=105,x2l=65,x2u=85)

##GPMDB 
GPMDBMH(x0=c(71,79),x1=c(96.816,71.926),inter=1,cvec=c(0.5,0.0025,0.3),m=50,it=100,plot=TRUE,x1l=65,x1u=105,x2l=65,x2u=85)
