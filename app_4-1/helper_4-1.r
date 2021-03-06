##Helper functions for application 4.1
##Birth-death, fully observed
######################################################################
#drift and diffusion

alpha=function(x,cvec)                  ##setting the drift (alpha)
{
  a=(cvec[1]-cvec[2])*x[1]
  a                                
}

beta=function(x,cvec)                   ##setting the diffusion (beta)
{
  b=(cvec[1]+cvec[2])*x[1]
  b                            
}

###################################################################### 
##suplimentary functions

matfun<-function(A,fun)   
{
  # applies a function fun to a matrix A, in a sensible way...
  s<-svd(A)
  ( s$u * fun(s$d) * s$v )
}

rmvn<-function(mean,V)
{
  eps<-rnorm(length(V)) 
  sq<-matfun(V,sqrt)
  mean + ( sq %*% eps )
}

lmvnpdf=function(x,mean,V) 
{
  detV=V
  Vinv=solve(V)
  -0.5*log(detV)-0.5*t(x-mean)%*%Vinv%*%(x-mean)
}

######################################################################
##LNA set up of the model and solving the system
##The ODE system here is can be solved analytically

BD_solved = function(t,parms,x_0) 
{
  eta = x_0*exp((parms[1]-parms[2])*t)
  G = exp((parms[1]-parms[2])*t)
  Psi = (((parms[1]+parms[2])*x_0)/(parms[1]-parms[2]))*(1 - exp(-(parms[1]-parms[2])*t))

  list(eta,Psi,G)
} 

######################################################################
##Targets and Proposals

##Euler likelihood (target)

llikE=function(xmat,cvec,dt)
{
  ll=0
  m=length(xmat[,1])-1
  for (i in 1:m) 
  {
    mean=xmat[i,]+alpha(xmat[i,],cvec)*dt
    ll=ll+lmvnpdf(xmat[i+1,],mean,beta(xmat[i,],cvec)*dt)
  }
  ll
}

##MDB likelihood (proposal)
llikMDB=function(xmat,cvec,dt)
{
  ll=0
  m=length(xmat[,1])-1
  x1=xmat[m+1,]
  inter=m*dt
  
  for (i in 1:(m-1)) 
  {
    t=i*dt
    mean=((inter-t)*xmat[i,]+x1*dt)/(inter-t+dt)
    var=((inter-t)/(inter-t+dt))*beta(xmat[i,],cvec)*dt             
    ll=ll+lmvnpdf(xmat[i+1,],mean,var)
  }
  ll
}

##LB likelihood (proposal)
llikLB=function(xmat,cvec,dt,gamma)
{
  ll = 0
  m = length(xmat[,1])-1
  x1 = xmat[m+1,]
  inter = m*dt
  
  for (i in 1:(m-1)) 
  {
    t = i*dt
    al = alpha(xmat[i,],cvec)
    bet = beta(xmat[i,],cvec)
    
    mean = xmat[i,] + al*dt + ((dt)/((inter-t+dt)+gamma*((inter-t)^2 )/(dt)))*(x1-xmat[i,]-al*(inter-t+dt))
    var = bet*dt - ((dt^2)/((inter-t+dt)+gamma*((inter-t)^2 )/(dt)))*bet
    
    ll = ll + lmvnpdf(xmat[i+1,],mean,var)
  }
  ll
}

##RB likelihood (proposal)
llikRB=function(xmat,cvec,dt,temp_x)
{
  ll=0
  m=length(xmat[,1])-1
  x1=xmat[m+1,]
  inter=m*dt
  
  for (i in 1:(m-1)) 
  {
    t=i*dt
    mean=((inter-t)*xmat[i,]+x1*dt)/(inter-t+dt)
    var=((inter-t)/(inter-t+dt))*beta(temp_x[i,],cvec)*dt             
    ll=ll+lmvnpdf(xmat[i+1,],mean,var)
  }
  ll
}

##RB^- likelihood (proposal)
llikRBminus=function(xmat,cvec,dt,temp_x)
{
  ll=0
  m=length(xmat[,1])-1
  x1=xmat[m+1,]
  inter=m*dt
  
  for (i in 1:(m-1)) 
  {
    t=i*dt
    mean=((inter-t)*xmat[i,]+x1*dt)/(inter-t+dt)
    var=((inter-t)/(inter-t+dt))*beta(temp_x[i,],cvec)*dt             
    ll=ll+lmvnpdf(xmat[i+1,],mean,var)
  }
  ll
}

#GP-N likelihood (proposal)
llikGPN=function(xmat,cvec,dt,eta)
{
  ll=0
  m=length(xmat[,1])-1
  x1=xmat[m+1,]
  inter=m*dt
  
  for (i in 1:(m-1)) 
  {
    t=i*dt

    A = exp((cvec[1]-cvec[2])*(inter-t+dt))
    mu_t = eta[m+1,1] + A*(xmat[i,1]-eta[i,1])
    out_end = BD_solved((inter-t+dt),cvec,eta[i,1])
    Psi_end = out_end[[2]]
    G_end = out_end[[3]]

    temp1 = A*(1/(G_end*Psi_end*G_end))*(x1-mu_t)

    mean = xmat[i,] + alpha(xmat[i,],cvec)*dt + beta(xmat[i,],cvec)*temp1*dt
    ll = ll + lmvnpdf(xmat[i+1,],mean,beta(xmat[i,],cvec)*dt)          
  }
  ll
}

#GP-S likelihood (proposal)
llikGPS=function(xmat,cvec,dt)
{
  ll=0
  m=length(xmat[,1])-1
  x1=xmat[m+1,]
  inter=m*dt
  
  bet_t = beta(x1,cvec)

  for (i in 1:(m-1)) 
  {
    t=i*dt

    eta_end = BD_solved((inter-t+dt),cvec,xmat[i,])[[1]]

    bet = beta(xmat[i,],cvec)
    al = alpha(xmat[i,],cvec)

    mean = xmat[i,] + (al + bet*(1/bet_t)*((x1-eta_end)/(inter-t+dt)))*dt	

    ll = ll + lmvnpdf(xmat[i+1,],mean,bet*dt)          
  }
  ll
}

#GP likelihood (proposal)
llikGP=function(xmat,cvec,dt)
{
  ll=0
  m=length(xmat[,1])-1
  x1=xmat[m+1,]
  inter=m*dt
  
  for (i in 1:(m-1)) 
  {
    t=i*dt

    bet = beta(xmat[i,],cvec)
    al = alpha(xmat[i,],cvec)

    A = exp((cvec[1]-cvec[2])*(inter-t+dt))
	
    out = BD_solved((inter-t+dt),cvec,xmat[i,1])

    mu_t = out[[1]] 

    Psi_end = out[[2]]
    G_end = out[[3]]
    
    temp1 = A*(1/(G_end*Psi_end*G_end))*(x1-mu_t)

    mean = xmat[i,] + (al + bet*temp1)*dt 
    var = bet*dt

    ll = ll + lmvnpdf(xmat[i+1,],mean,var)          
  }
  ll
}

#GP-MDB likelihood (proposal)
llikGPMDB=function(xmat,cvec,dt)
{
  ll=0
  m=length(xmat[,1])-1
  x1=xmat[m+1,]
  inter=m*dt
  
  for (i in 1:(m-1)) 
  {
    t=i*dt

    bet = beta(xmat[i,],cvec)
    al = alpha(xmat[i,],cvec)

    A = exp((cvec[1]-cvec[2])*(inter-t+dt))
	
    out = BD_solved((inter-t+dt),cvec,xmat[i,1])

    mu_t = out[[1]] 
    Psi_end = out[[2]]
    G_end = out[[3]]
    
    temp1 = A*(1/(G_end*Psi_end*G_end))*(x1-mu_t)

    mean = xmat[i,] + (al + bet*temp1)*dt 
    var = ((inter-t)/(inter-t+dt))*bet*dt

    ll = ll + lmvnpdf(xmat[i+1,],mean,var)          
  }
  ll
}

######################################################################
##Path proposals (single path)

##MDB
MDBpath=function(x0,x1,inter,cvec,m)
{
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)    						
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)										
  prop=cpath
  
  for (i in 1:(m-1)) 
  {    
    t = i*dt
    
    mean=((inter-t)*prop[i,] + x1*dt)/(inter-t+dt)						
    var=((inter-t)/(inter-t+dt))*beta(prop[i,],cvec)*dt 		
      
    prop[i+1,] = rmvn(mean,var)
  }
  prop
}

##LB
LBpath=function(x0,x1,inter,cvec,m,gamma)
{
  dt = inter/m
  
  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],x1[1],length.out=m+1)
  prop = cpath

  for (i in 1:(m-1)) 
  {           
    t = i*dt
    al = alpha(prop[i,],cvec)
    bet = beta(prop[i,],cvec)
    
    mean = prop[i,] + al*dt + ((dt)/((inter-t+dt)+gamma*((inter-t)^2 )/(dt)))*(x1-prop[i,]-al*(inter-t+dt))
    var = bet*dt - ((dt^2)/((inter-t+dt)+gamma*((inter-t)^2 )/(dt)))*bet
    
    prop[i+1,] = rmvn(mean,var)
  }
  prop
}

##RB
RBpath=function(x0,x1,inter,cvec,m,etaprop)
{
  dt = inter/m

  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],x1[1],length.out=m+1)
  prop = cpath
  
  mpath = matrix(0,ncol=length(x0),nrow=m+1)
  
  mpath[m+1,1] = x1[1] - etaprop[m+1,1]
  
  for(i in 1:(m-1))
  {           
    t = i*dt
    
    bet = beta(prop[i,],cvec)
    
    mean=((inter-t)*mpath[i,] + mpath[m+1,]*dt)/(inter-t+dt)    				
    var=((inter-t)/(inter-t+dt))*bet*dt 	
    
    mpath[i+1,] = rmvn(mean,var)
    
    prop[i+1,] = etaprop[i+1,] + mpath[i+1,]
  } 
  prop
}

##RB^-
RBminuspath=function(x0,x1,inter,cvec,m,mprop_c,etaprop)
{
  dt = inter/m

  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],x1[1],length.out=m+1)
  prop = cpath
  
  mpath = matrix(0,ncol=length(x0),nrow=m+1)
  
  mpath[m+1,1] = x1[1] - etaprop[m+1,1] - mprop_c[m+1,1]

  for(i in 1:(m-1))
  {           
    t = i*dt
    
    bet = beta(prop[i,],cvec)
    
    mean=((inter-t)*mpath[i,] + mpath[m+1,]*dt)/(inter-t+dt)    				
    var=((inter-t)/(inter-t+dt))*bet*dt 	
    
    mpath[i+1,] = rmvn(mean,var)
    
    prop[i+1,] = etaprop[i+1,] + mpath[i+1,] + mprop_c[i+1,]
  } 
  prop
}

##GP-N
GPNpath=function(x0,x1,inter,cvec,m,etaprop)
{
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)    						
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)										
  prop=cpath
  
  for (i in 1:(m-1)) 
  {    
    t = i*dt
  
    bet = beta(prop[i,],cvec)
    al = alpha(prop[i,],cvec)

    A = exp((cvec[1]-cvec[2])*(inter-t+dt))

    mu_t = etaprop[m+1,1] + A*(prop[i,1]-etaprop[i,1])
    out_end = BD_solved((inter-t+dt),cvec,etaprop[i,1])
    Psi_end = out_end[[2]]
    G_end = out_end[[3]]

    temp1 = A*(1/(G_end*Psi_end*G_end))*(x1-mu_t)

    mean = prop[i,] + (al + bet*temp1)*dt 
    var = bet*dt

    prop[i+1,] = rmvn(mean,var)
  }
  prop
}

##GP-S
GPSpath=function(x0,x1,inter,cvec,m)
{
  dt=inter/m

  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],x1[1],length.out=m+1)
  prop = cpath
  
  bet_t = beta(x1,cvec)

  for (i in 1:(m-1)) 
  {    
    t = i*dt

    bet = beta(prop[i,],cvec)
    al = alpha(prop[i,],cvec)

    eta_end = BD_solved((inter-t+dt),cvec,prop[i,])[[1]]

    mean = prop[i,] + (al + bet*(1/bet_t)*((x1-eta_end)/(inter-t+dt)))*dt
    var=bet*dt 

    prop[i+1,] = rmvn(mean,var)
  }
  prop
}

##GP
GPpath=function(x0,x1,inter,cvec,m)
{
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)    						
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)										
  prop=cpath
 
  for (i in 1:(m-1)) 
  {    
    t = i*dt
  
    bet = beta(prop[i,],cvec)
    al = alpha(prop[i,],cvec)

    A = exp((cvec[1]-cvec[2])*(inter-t+dt))
	
    out = BD_solved((inter-t+dt),cvec,prop[i,1])

    mu_t = out[[1]] 
    Psi_end = out[[2]]
    G_end = out[[3]]

    temp1 = A*(1/(G_end*Psi_end*G_end))*(x1-mu_t)

    mean = prop[i,] + (al + bet*temp1)*dt 
    var = bet*dt

    prop[i+1,] = rmvn(mean,var)	  
  }
  prop
}

##GP-MDB
GPMDBpath=function(x0,x1,inter,cvec,m)
{
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)    						
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)										
  prop=cpath
  
  for (i in 1:(m-1)) 
  {    
    t = i*dt
  
    bet = beta(prop[i,],cvec)
    al = alpha(prop[i,],cvec)

    A = exp((cvec[1]-cvec[2])*(inter-t+dt))
	
    out = BD_solved((inter-t+dt),cvec,prop[i,1])

    mu_t = out[[1]] 
    Psi_end = out[[2]]
    G_end = out[[3]]

    temp1 = A*(1/(G_end*Psi_end*G_end))*(x1-mu_t)

    mean = prop[i,] + (al + bet*temp1)*dt 
    var = ((inter-t)/(inter-t+dt))*bet*dt

    prop[i+1,] = rmvn(mean,var)	  
  }
  prop
}

######################################################################
##Metropolis-Hastings samplers

##MDB
MDBMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,xl,xu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)	
  
  for(j in 1:it)
  {
    prop=MDBpath(x0,x1,inter,cvec,m)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikMDB(cpath,cvec,dt)-llikMDB(prop,cvec,dt)} 

    if(log(runif(1))<lap)
    {
      cpath=prop
      count=count+1
			
      if(plot==T)
      {
        if(count==1)
        {
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(xl,xu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          lines(ts(prop[,1],start=0,deltat=dt),col=col[j])
        }
      }
    }
  }
  count/it
}

##LB
LBMH=function(x0,x1,inter,cvec,m,gamma,it,plot=FALSE,xl,xu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)		
  
  for(j in 1:it)
  {
    prop=LBpath(x0,x1,inter,cvec,m,gamma)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikLB(cpath,cvec,dt,gamma)-llikLB(prop,cvec,dt,gamma)} 

    if(log(runif(1))<lap)
    {
      cpath=prop
      count=count+1
			
      if(plot==T)
      {
        if(count==1)
        {
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(xl,xu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          lines(ts(prop[,1],start=0,deltat=dt),col=col[j])
        }
      }
    }
  }
  count/it
}

##RB
RBMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,xl,xu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],x1[1],length.out=m+1)	

  etapath = matrix(0,ncol=length(x0),nrow=m+1)  
  etapath[,1]=seq(x0[1],x1[1],length.out=m+1)

  for(i in 1:m)
  { 
    out = BD_solved(dt,cvec,etapath[i,])    
    etapath[i+1,] = out[[1]]
  }

  for(j in 1:it)
  {
    prop=RBpath(x0,x1,inter,cvec,m,etapath)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikRB(cpath-etapath,cvec,dt,cpath)-llikRB(prop-etapath,cvec,dt,prop)}

    if(log(runif(1))<lap)
    {
      cpath=prop
      count=count+1
			
      if(plot==T)
      {
        if(count==1)
        {
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(xl,xu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          lines(ts(prop[,1],start=0,deltat=dt),col=col[j])
        }
      }
    }
  }
  count/it
}

##RB^-
RBminusMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,xl,xu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],x1[1],length.out=m+1)	

  etapath = matrix(0,ncol=length(x0),nrow=m+1)  
  etapath[,1]=seq(x0[1],x1[1],length.out=m+1)

  mpath_c = matrix(0,ncol=length(x0),nrow=m+1)

  for(i in 1:(m))
  { 
    out = BD_solved(dt,cvec,etapath[i,])    
    etapath[i+1,] = out[[1]]
  }

  for(i in 1:m)
  {
    t = i*dt

    out_end = BD_solved(inter-t+dt,cvec,etapath[i,1])
    Psi_end = out_end[[2]]
    G_end = out_end[[3]]

    out_t = BD_solved(dt,cvec,etapath[i,1])
    Psi_t = out_t[[2]]
    G_t = out_t[[3]]

    stoch_mean = G_t%*%mpath_c[i,] + G_t%*%Psi_t%*%t(G_end)%*%solve(G_end%*%Psi_end%*%t(G_end))%*%((x1[1]-etapath[m+1,1])-G_end%*%mpath_c[i,1])

    mpath_c[i+1,] = stoch_mean
  } 

  for(j in 1:it)
  {
    prop=RBminuspath(x0,x1,inter,cvec,m,mpath_c,etapath)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikRBminus(cpath-etapath-mpath_c,cvec,dt,cpath)-llikRBminus(prop-etapath-mpath_c,cvec,dt,prop)}

    if(log(runif(1))<lap)
    {
      cpath=prop
      count=count+1
			
      if(plot==T)
      {
        if(count==1)
        {
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(xl,xu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          lines(ts(prop[,1],start=0,deltat=dt),col=col[j])
        }
      }
    }
  }
  count/it
}

##GP-N
GPNMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,xl,xu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)	
  etapath = matrix(0,ncol=length(x0),nrow=m+1)
  etapath[,1] = seq(x0,x1,length.out=m+1)

  for(i in 1:m)
  {           
    etapath[i+1,] = BD_solved(dt,cvec,etapath[i,1])[[1]][1]
  }

  for(j in 1:it)
  {
    prop=GPNpath(x0,x1,inter,cvec,m,etapath)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikGPN(cpath,cvec,dt,etapath)-llikGPN(prop,cvec,dt,etapath)} 
		
    if(log(runif(1))<lap)
    {
      cpath=prop
      count=count+1
			
      if(plot==T)
      {
        if(count==1)
        {
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(xl,xu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          lines(ts(prop[,1],start=0,deltat=dt),col=col[j])
        }
      }
    }
  }
  count/it
}

##GP-S
GPSMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,xl,xu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],x1[1],length.out=m+1)	

  for(j in 1:it)
  {
    prop=GPSpath(x0,x1,inter,cvec,m)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikGPS(cpath,cvec,dt)-llikGPS(prop,cvec,dt)}

    if(log(runif(1))<lap)
    {
      cpath=prop
      count=count+1
			
      if(plot==T)
      {
        if(count==1)
        {
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(xl,xu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          lines(ts(prop[,1],start=0,deltat=dt),col=col[j])
        }
      }
    }
  }
  count/it
}

##GP
GPMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,xl,xu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)	

  for(j in 1:it)
  {
    prop=GPpath(x0,x1,inter,cvec,m)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikGP(cpath,cvec,dt)-llikGP(prop,cvec,dt)}
 
    if(log(runif(1))<lap)
    {
      cpath=prop
      count=count+1
			
      if(plot==T)
      {
        if(count==1)
        {
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(xl,xu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          lines(ts(prop[,1],start=0,deltat=dt),col=col[j])
        }
      }
    }
  }
  count/it
}

##GP-MDB
GPMDBMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,xl,xu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)	

  for(j in 1:it)
  {
    prop=GPMDBpath(x0,x1,inter,cvec,m)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikGPMDB(cpath,cvec,dt)-llikGPMDB(prop,cvec,dt)}
 
    if(log(runif(1))<lap)
    {
      cpath=prop
      count=count+1
			
      if(plot==T)
      {
        if(count==1)
        {
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(xl,xu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          lines(ts(prop[,1],start=0,deltat=dt),col=col[j])
        }
      }
    } 
  }
  count/it
}
