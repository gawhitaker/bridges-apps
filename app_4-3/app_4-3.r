##Acceptance probabilities for application 4.3 with option to plot
##Aphid growth, partially observed with error
######################################################################
##Load helper functions and required libraries

#install.packages("deSolve") ##uncomment to install "deSolve" package if required 
library(deSolve)
source("helper_4-3.r")

######################################################################
##Options for Metropolis-Hastings samplers
##x0 = vector of initial value of the process (N,C) (left-endpoint) 
##y = observation at T with error (N) (right-endpoint) 
##inter = T-0
##cvec = parameter vector (theta_1,theta_2)
##sigma = standard deviation of noise
##m = discretisation to use
##it = number of iterations 
##plot = logical for plot
##nl = lower y limit for N in plot
##nu = upper y limit for N in plot
##cl = lower y limit for C in plot
##cu = upper y limit for C in plot

######################################################################
##Metropolis-Hastings samplers for median values with sigma=5

##EM
EMMH(x0=c(347.548,398.937),y=786.089,inter=1.28,cvec=c(1.45,0.0009),sigma=5.0,m=50,it=100,plot=TRUE,nl=300,nu=900,cl=300,cu=1700)

##RB 
RBMH(x0=c(347.548,398.937),y=786.089,inter=1.28,cvec=c(1.45,0.0009),sigma=5.0,m=50,it=100,plot=TRUE,nl=300,nu=900,cl=300,cu=1700)

##RB^-
RBminusMH(x0=c(347.548,398.937),y=786.089,inter=1.28,cvec=c(1.45,0.0009),sigma=5.0,m=50,it=100,plot=TRUE,nl=300,nu=900,cl=300,cu=1700)

##GP 
GPMH(x0=c(347.548,398.937),y=786.089,inter=1.28,cvec=c(1.45,0.0009),sigma=5.0,m=50,it=100,plot=TRUE,nl=300,nu=900,cl=300,cu=1700)

##GPMDB 
GPMDBMH(x0=c(347.548,398.937),y=786.089,inter=1.28,cvec=c(1.45,0.0009),sigma=5.0,m=50,it=100,plot=TRUE,nl=300,nu=900,cl=300,cu=1700)
