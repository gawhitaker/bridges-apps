##Helper functions for application 4.2
##Lotka-Volterra, fully observed
######################################################################
##drift and diffusion

alpha=function(x,cvec)                  ##setting the drift (alpha)
{
  a=c(cvec[1]*x[1]-cvec[2]*x[1]*x[2],cvec[2]*x[1]*x[2]-cvec[3]*x[2])
  a                                  
}

beta=function(x,cvec)                   ##setting the diffusion (beta)
{
  b=matrix(c(cvec[1]*x[1]+cvec[2]*x[1]*x[2],-cvec[2]*x[1]*x[2],-cvec[2]*x[1]*x[2],cvec[3]*x[2]+cvec[2]*x[1]*x[2]),nrow=2,ncol=2)
  b                              
}

###################################################################### 
##suplimentary functions

matfun<-function(A,fun)
{
  ##applies a function fun to a matrix A, in a sensible way...
  s<-svd(A)
  ( s$u %*% diag( fun(s$d) ) %*% s$v )
}

rmvn<-function(mean,V)
{
  eps<-rnorm(dim(V)[1]) 
  sq<-matfun(V,sqrt)
  mean + ( sq %*% eps )
}

lmvnpdf=function(x,mean,V)
{
  detV=det(V)
  Vinv=solve(V)
  -0.5*log(detV)-0.5*t(x-mean)%*%Vinv%*%(x-mean)
}

######################################################################
##LNA set up of the model and solving the system
##The ODE system here is solved numerically
##The "deSolve" package is required

LV_solved = function(t,parms,x_0) 
{
  x_start = c(G=matrix(c(1,0,0,1),ncol=2),Psi=matrix(c(0,0,0,0),ncol=2),eta=matrix(c(x_0[1],x_0[2]),ncol=1))
  times = c(0,t) 
  out = lsoda(x_start,times,LVMod,parms)

  G = matrix(c(out[2,2],out[2,3],out[2,4],out[2,5]),ncol=2)
  Psi = matrix(c(out[2,6],out[2,7],out[2,8],out[2,9]),ncol=2)
  eta = c(out[2,10],out[2,11])

  list(eta,Psi,G)
}

LVMod = function(t,x,parms)  ##set-up the model 
{
  G = matrix(c(x[1],x[2],x[3],x[4]),ncol=2)
  Psi = matrix(c(x[5],x[6],x[7],x[8]),ncol=2)
  eta = matrix(c(x[9],x[10]),ncol=1)

  mat = matrix(c(parms[1]-parms[2]*eta[2],parms[2]*eta[2],-parms[2]*eta[1],parms[2]*eta[1]-parms[3]),nrow=2,ncol=2)
  
  dG = mat%*%G    
  dPsi = solve(G)%*%beta(eta,parms)%*%solve(t(G))
  deta = alpha(eta,parms)
 
  res = c(dG,dPsi,deta)
  list(res)
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

##GP likelihood (proposal)
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
	
    out = LV_solved((inter-t+dt),cvec,xmat[i,])

    mu_t = out[[1]] 
    Psi = out[[2]]
    G = out[[3]]

    mean = al + bet%*%t(G)%*%solve(G%*%Psi%*%t(G))%*%(x1-mu_t) 
    mean = xmat[i,] + mean*dt
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
	
    out = LV_solved((inter-t+dt),cvec,xmat[i,])

    mu_t = out[[1]] 
    Psi = out[[2]]
    G = out[[3]]

    mean = al + bet%*%t(G)%*%solve(G%*%Psi%*%t(G))%*%(x1-mu_t) 
    mean = xmat[i,] + mean*dt
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
  cpath[1:(m+1),2]=seq(x0[2],x1[2],length.out=m+1)											
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
  cpath[,2] = seq(x0[2],x1[2],length.out=m+1)
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
  cpath[,2] = seq(x0[2],x1[2],length.out=m+1)
  prop = cpath  
  mpath = matrix(0,ncol=length(x0),nrow=m+1)
  
  mpath[m+1,1] = x1[1] - etaprop[m+1,1]
  mpath[m+1,2] = x1[2] - etaprop[m+1,2]
  
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
  cpath[,2] = seq(x0[2],x1[2],length.out=m+1)
  prop = cpath
  mpath = matrix(0,ncol=length(x0),nrow=m+1) 
  
  mpath[m+1,1] = x1[1] - etaprop[m+1,1] - mprop_c[m+1,1]
  mpath[m+1,2] = x1[2] - etaprop[m+1,2] - mprop_c[m+1,2]

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

##GP
GPpath=function(x0,x1,inter,cvec,m)
{
  dt=inter/m

  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],x1[1],length.out=m+1)
  cpath[,2] = seq(x0[2],x1[2],length.out=m+1)
  prop = cpath
  
  for (i in 1:(m-1)) 
  {    
    t = i*dt
  
    bet = beta(prop[i,],cvec)
    al = alpha(prop[i,],cvec)

    out = LV_solved((inter-t+dt),cvec,prop[i,])

    mu_t = out[[1]] 
    Psi = out[[2]]
    G = out[[3]]

    mean = al + bet%*%t(G)%*%solve(G%*%Psi%*%t(G))%*%(x1-mu_t) 
    mean = prop[i,] + mean*dt
    var=bet*dt 

    prop[i+1,] = rmvn(mean,var)			 
  }
  prop
}

##GP-MDB
GPMDBpath=function(x0,x1,inter,cvec,m)
{
  dt=inter/m

  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],x1[1],length.out=m+1)
  cpath[,2] = seq(x0[2],x1[2],length.out=m+1)
  prop = cpath
  
  for (i in 1:(m-1)) 
  {    
    t = i*dt
  
    bet = beta(prop[i,],cvec)
    al = alpha(prop[i,],cvec)

    out = LV_solved((inter-t+dt),cvec,prop[i,])

    mu_t = out[[1]] 
    Psi = out[[2]]
    G = out[[3]]

    mean = al + bet%*%t(G)%*%solve(G%*%Psi%*%t(G))%*%(x1-mu_t) 
    mean = prop[i,] + mean*dt
	
    var=((inter-t)/(inter-t+dt))*bet*dt 

    prop[i+1,] = rmvn(mean,var)			 
  }
  prop
}

######################################################################
##Metropolis-Hastings samplers

##MDB
MDBMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,x1l,x1u,x2l,x2u)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)	
  cpath[1:(m+1),2]=seq(x0[2],x1[2],length.out=m+1)	
  
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
          split.screen(c(1, 2))
          screen(1)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(x1l,x1u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[1]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(x2l,x2u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[2]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(x1l,x1u),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(x2l,x2u),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}

##LB
LBMH=function(x0,x1,inter,cvec,m,gamma,it,plot=FALSE,x1l,x1u,x2l,x2u)
{
  count=0
  col=rainbow(it)
  dt=inter/m
  #Initialise
  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[1:(m+1),1]=seq(x0[1],x1[1],length.out=m+1)	
  cpath[1:(m+1),2]=seq(x0[2],x1[2],length.out=m+1)	
  
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
          split.screen(c(1, 2))
          screen(1)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(x1l,x1u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[1]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(x2l,x2u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[2]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(x1l,x1u),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(x2l,x2u),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}

##RB
RBMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,x1l,x1u,x2l,x2u)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],x1[1],length.out=m+1)	
  cpath[,2]=seq(x0[2],x1[2],length.out=m+1)	

  etapath = matrix(0,ncol=length(x0),nrow=m+1)  
  etapath[,1]=seq(x0[1],x1[1],length.out=m+1)	
  etapath[,2]=seq(x0[2],x1[2],length.out=m+1)	

  for(i in 1:(m))
  { 
    out = LV_solved(dt,cvec,etapath[i,])    
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
          split.screen(c(1, 2))
          screen(1)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(x1l,x1u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[1]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(x2l,x2u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[2]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(x1l,x1u),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(x2l,x2u),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}

##RB^-
RBminusMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,x1l,x1u,x2l,x2u)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],x1[1],length.out=m+1)	
  cpath[,2]=seq(x0[2],x1[2],length.out=m+1)	

  etapath = matrix(0,ncol=length(x0),nrow=m+1)  
  etapath[,1]=seq(x0[1],x1[1],length.out=m+1)	
  etapath[,2]=seq(x0[2],x1[2],length.out=m+1)	

  mpath_c = matrix(0,ncol=length(x0),nrow=m+1)

  for(i in 1:(m))
  { 
    out = LV_solved(dt,cvec,etapath[i,])    
    etapath[i+1,] = out[[1]]	
  }

  for(i in 1:(m))
  {
    t = i*dt

    out_end = LV_solved(inter-t+dt,cvec,etapath[i,])
    Psi_end = out_end[[2]]
    G_end = out_end[[3]]

    out_t = LV_solved(dt,cvec,etapath[i,])
    Psi_t = out_t[[2]]
    G_t = out_t[[3]]

    stoch_mean = G_t%*%mpath_c[i,] + G_t%*%Psi_t%*%t(G_end)%*%solve(G_end%*%Psi_end%*%t(G_end))%*%((x1-etapath[m+1,])-G_end%*%mpath_c[i,])

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
          split.screen(c(1, 2))
          screen(1)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(x1l,x1u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[1]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(x2l,x2u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[2]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(x1l,x1u),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(x2l,x2u),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}

##GP
GPMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,x1l,x1u,x2l,x2u)
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
          split.screen(c(1, 2))
          screen(1)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(x1l,x1u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[1]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(x2l,x2u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[2]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(x1l,x1u),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(x2l,x2u),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}

##GP-MDB
GPMDBMH=function(x0,x1,inter,cvec,m,it,plot=FALSE,x1l,x1u,x2l,x2u)
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
          split.screen(c(1, 2))
          screen(1)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(x1l,x1u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[1]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(x2l,x2u),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*X[2]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(x1l,x1u),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(x2l,x2u),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}
