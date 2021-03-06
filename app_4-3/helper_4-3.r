##Helper functions for application 4.3
##Aphid growth, partially observed with error
######################################################################
##drift and diffusion

alpha=function(x,cvec)                  ##setting the drift (alpha)
{
  a=c(cvec[1]*x[1]-cvec[2]*x[1]*x[2],cvec[1]*x[1])
  a                                
}

beta=function(x,cvec)                   ##setting the diffusion (beta)
{
  b=matrix(c(cvec[1]*x[1]+cvec[2]*x[1]*x[2],cvec[1]*x[1],cvec[1]*x[1],cvec[1]*x[1]),nrow=2,ncol=2)
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

Aphid_solved = function(t,parms,x_0) 
{
  x_start = c(G=matrix(c(1,0,0,1),ncol=2),Psi=matrix(c(0,0,0,0),ncol=2),eta=matrix(c(x_0[1],x_0[2]),ncol=1))
  times = c(0,t) 
  out = lsoda(x_start,times,AphidMod,parms)

  G = matrix(c(out[2,2],out[2,3],out[2,4],out[2,5]),ncol=2)
  Psi = matrix(c(out[2,6],out[2,7],out[2,8],out[2,9]),ncol=2)
  eta = c(out[2,10],out[2,11])

  list(eta,Psi,G)
}

AphidMod = function(t,x,parms)  ##set-up the model 
{
  G = matrix(c(x[1],x[2],x[3],x[4]),ncol=2)
  Psi = matrix(c(x[5],x[6],x[7],x[8]),ncol=2)
  eta = matrix(c(x[9],x[10]),ncol=1)

  mat = matrix(c(parms[1]-parms[2]*eta[2],parms[1],-parms[2]*eta[1],0),nrow=2,ncol=2)
  
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

##RB partial likelihood with noise (proposal)
llikRBpn=function(xmat,cvec,dt,temp_x,temp_eta,G,sigma,y)
{
  ll=0
  m=length(xmat[,1])-1
  inter=m*dt
  
  for (i in 1:m) 
  {
    t=i*dt

    al_x = alpha(temp_x[i,],cvec)
    al_eta = alpha(temp_eta[i,],cvec)
    bet = beta(temp_x[i,],cvec)
	
    mean = xmat[i,] + (al_x-al_eta)*dt + (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(y-temp_eta[m+1,1]-(t(G)%*%xmat[i,]+t(G)%*%(al_x-al_eta)*(inter-t+dt)))      				
    var =  bet*dt - (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(t(G)%*%bet*dt)    	
	     
    ll=ll+lmvnpdf(xmat[i+1,],mean,var)
  }
  ll
}

##RB^- partial likelihood with noise (proposal)
llikRBminuspn=function(xmat,cvec,dt,G,sigma,y,temp_x,temp_eta,temp_m)
{
  ll=0
  m=length(xmat[,1])-1
  inter=m*dt
  
  for (i in 1:m) 
  {
    t=i*dt

    al = alpha(temp_x[i,],cvec)
    bet = beta(temp_x[i,],cvec)

    mean = xmat[i,] + (al-((temp_eta[i+1,]-temp_eta[i,])/dt)-((temp_m[i+1,]-temp_m[i,])/dt))*dt + (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(y-temp_eta[m+1,1]-temp_m[m+1,1]-(t(G)%*%xmat[i,]+t(G)%*%(al-((temp_eta[i+1,]-temp_eta[i,])/dt)-((temp_m[i+1,]-temp_m[i,])/dt))*(inter-t+dt)))    
    var =  bet*dt - (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(t(G)%*%bet*dt)    	
	     
    ll=ll+lmvnpdf(xmat[i+1,],mean,var)
  }
  ll
}

##GP partial likelihood with noise (proposal)
llikGPpn=function(xmat,cvec,dt,G,sigma,y)
{
  ll=0
  m=length(xmat[,1])-1
  inter=m*dt
  
  for (i in 1:m) 
  {
    t=i*dt

    bet = beta(xmat[i,],cvec)
    al = alpha(xmat[i,],cvec)

    out = Aphid_solved((inter-t+dt),cvec,xmat[i,])

    mu_t = out[[1]] 
    Psi = out[[2]]
    Gmat = out[[3]]

    mean = al + bet%*%t(Gmat)%*%G%*%solve(t(G)%*%Gmat%*%Psi%*%t(Gmat)%*%G+sigma^2)%*%(y-t(G)%*%mu_t) 
    mean = xmat[i,] + mean*dt
    var = bet*dt

    ll = ll + lmvnpdf(xmat[i+1,],mean,var)          
  }
  ll
}

##GP-MDB partial likelihood with noise (proposal)
llikGPMDBpn=function(xmat,cvec,dt,G,sigma,y)
{
  ll=0
  m=length(xmat[,1])-1
  inter=m*dt
  
  for (i in 1:m) 
  {
    t=i*dt

    bet = beta(xmat[i,],cvec)
    al = alpha(xmat[i,],cvec)

    out = Aphid_solved((inter-t+dt),cvec,xmat[i,])

    mu_t = out[[1]] 
    Psi = out[[2]]
    Gmat = out[[3]]

    mean = al + bet%*%t(Gmat)%*%G%*%solve(t(G)%*%Gmat%*%Psi%*%t(Gmat)%*%G+sigma^2)%*%(y-t(G)%*%mu_t) 
    mean = xmat[i,] + mean*dt
    var = bet*dt - (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(t(G)%*%bet*dt) 

    ll = ll + lmvnpdf(xmat[i+1,],mean,var)          
  }
  ll
}

######################################################################
##Path proposals (single path)

##EM
EMpath=function(x0,y,inter,cvec,m)
{
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)    						
  cpath[,1]=seq(x0[1],y,length.out=m+1)	
  cpath[,2]=rep(x0[2],m+1)										
  prop=cpath
  
  for (i in 1:m) 
  {    
    t = i*dt
    
    mean = prop[i,] + alpha(prop[i,],cvec)*dt				
    var = beta(prop[i,],cvec)*dt 	
      
    prop[i+1,] = rmvn(mean,var)	  
  }
  prop
}

##RB
RBpath=function(x0,y,inter,cvec,sigma,m,G,etaprop)
{
  dt = inter/m

  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],y,length.out=m+1)
  cpath[,2] = rep(x0[2],m+1)
  prop = cpath
  
  mpath = matrix(0,ncol=length(x0),nrow=m+1) 
  mpath[m+1,1] = y - etaprop[m+1,1]  
  
  for(i in 1:m)
  {         
    t = i*dt

    al_x = alpha(prop[i,],cvec)
    al_eta=alpha(etaprop[i,],cvec)
    bet = beta(prop[i,],cvec)
  
    mean = mpath[i,] + (al_x-al_eta)*dt + (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(mpath[m+1,1]-(t(G)%*%mpath[i,]+t(G)%*%(al_x-al_eta)*(inter-t+dt)))      				
    var =  bet*dt - (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(t(G)%*%bet*dt)    	
    
    mpath[i+1,] = rmvn(mean,var)
    
    prop[i+1,] = etaprop[i+1,] + mpath[i+1,]
  } 
  prop
}

##RB^-
RBminuspath=function(x0,y,inter,cvec,sigma,m,G,mprop_c,etaprop)
{
  dt = inter/m

  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],y,length.out=m+1)
  cpath[,2] = rep(x0[2],m+1)
  prop = cpath
  
  mpath = matrix(0,ncol=length(x0),nrow=m+1)
  
  mpath[m+1,1] = y - etaprop[m+1,1] - mprop_c[m+1,1]

  for(i in 1:m)
  {           
    t = i*dt
    
    al = alpha(prop[i,],cvec)
    bet = beta(prop[i,],cvec)
  
    mean = mpath[i,] + (al-((etaprop[i+1,]-etaprop[i,])/dt)-((mprop_c[i+1,]-mprop_c[i,])/dt))*dt + (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(mpath[m+1,1]-(t(G)%*%mpath[i,]+t(G)%*%(al-((etaprop[i+1,]-etaprop[i,])/dt)-((mprop_c[i+1,]-mprop_c[i,])/dt))*(inter-t+dt)))      	
    var =  bet*dt - (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(t(G)%*%bet*dt)    		
    
    mpath[i+1,] = rmvn(mean,var)
    
    prop[i+1,] = etaprop[i+1,] + mpath[i+1,] + mprop_c[i+1,]
  } 
  prop
}

##GP
GPpath=function(x0,y,inter,cvec,sigma,m,G)
{
  dt=inter/m

  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],y,length.out=m+1)
  cpath[,2] = rep(x0[2],m+1)
  prop = cpath
  
  for (i in 1:m) 
  {    
    t = i*dt
  
    bet = beta(prop[i,],cvec)
    al = alpha(prop[i,],cvec)

    out = Aphid_solved((inter-t+dt),cvec,prop[i,])

    mu_t = out[[1]] 
    Psi = out[[2]]
    Gmat = out[[3]]

    mean = al + bet%*%t(Gmat)%*%G%*%solve(t(G)%*%Gmat%*%Psi%*%t(Gmat)%*%G+sigma^2)%*%(y-t(G)%*%mu_t) 
    mean = prop[i,] + mean*dt	
    var=bet*dt 

    prop[i+1,] = rmvn(mean,var)			 
  }
  prop
}

##GP-MDB
GPMDBpath=function(x0,y,inter,cvec,sigma,m,G)
{
  dt=inter/m
  cpath = matrix(0,ncol=length(x0),nrow=m+1)
  cpath[,1] = seq(x0[1],y,length.out=m+1)
  cpath[,2] = rep(x0[2],m+1)
  prop = cpath
  
  for (i in 1:m) 
  {    
    t = i*dt
  
    bet = beta(prop[i,],cvec)
    al = alpha(prop[i,],cvec)

    out = Aphid_solved((inter-t+dt),cvec,prop[i,])

    mu_t = out[[1]] 
    Psi = out[[2]]
    Gmat = out[[3]]

    mean = al + bet%*%t(Gmat)%*%G%*%solve(t(G)%*%Gmat%*%Psi%*%t(Gmat)%*%G+sigma^2)%*%(y-t(G)%*%mu_t) 
    mean = prop[i,] + mean*dt	
    var = bet*dt - (bet%*%G*dt)%*%solve((t(G)%*%bet%*%G)*(inter-t+dt)+sigma^2)%*%(t(G)%*%bet*dt) 

    prop[i+1,] = rmvn(mean,var)		 		 
  }
  prop
}

######################################################################
##Metropolis-Hastings samplers

##EM
EMMH=function(x0,y,inter,cvec,m,sigma,it,plot=FALSE,nl,nu,cl,cu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],y,length.out=m+1)	
  cpath[,2]=rep(x0[2],m+1)	
  
  for(j in 1:it)
  {
    prop=EMpath(x0,y,inter,cvec,m)

    if(j==1){lap=1}
    else{lap= dnorm(y,prop[m+1,1],sigma,log=T) - dnorm(y,cpath[m+1,1],sigma,log=T) }

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
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(nl,nu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*N[t]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(cl,cu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*C[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(nl,nu),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(cl,cu),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}

##RB
RBMH=function(x0,y,inter,cvec,sigma,m,G=matrix(c(1,0),ncol=1,nrow=2),it,plot=FALSE,nl,nu,cl,cu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],y,length.out=m+1)	
  cpath[,2]=rep(x0[2],m+1)

  etapath = matrix(0,ncol=length(x0),nrow=m+1)  
  etapath[,1]=seq(x0[1],y,length.out=m+1)	
  etapath[,2]=rep(x0[2],m+1)

  for(i in 1:(m))
  { 
    out = Aphid_solved(dt,cvec,etapath[i,])    
    etapath[i+1,] = out[[1]]		
  }

  for(j in 1:it)
  {
    prop=RBpath(x0,y,inter,cvec,sigma,m,G,etapath)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikRBpn(cpath-etapath,cvec,dt,cpath,etapath,G,sigma,y)-llikRBpn(prop-etapath,cvec,dt,prop,etapath,G,sigma,y)+dnorm(y,prop[m+1,1],sigma,log=T)-dnorm(y,cpath[m+1,1],sigma,log=T)}

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
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(nl,nu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*N[t]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(cl,cu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*C[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(nl,nu),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(cl,cu),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}

##RB^-
RBminusMH=function(x0,y,inter,cvec,sigma,m,G=matrix(c(1,0),ncol=1,nrow=2),it,plot=FALSE,nl,nu,cl,cu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],y,length.out=m+1)	
  cpath[,2]=rep(x0[2],m+1)

  etapath = matrix(0,ncol=length(x0),nrow=m+1)  
  etapath[,1]=seq(x0[1],y,length.out=m+1)	
  etapath[,2]=rep(x0[2],m+1)

  mpath_c = matrix(0,ncol=length(x0),nrow=m+1)

  for(i in 1:(m))
  { 
    out = Aphid_solved(dt,cvec,etapath[i,])    
    etapath[i+1,] = out[[1]]
  }

  out_end = Aphid_solved(inter,cvec,etapath[1,])
  Psi_end = out_end[[2]]
  G_end = out_end[[3]]

  for(i in 1:m)
  {
    t = i*dt

    out_t = Aphid_solved(t,cvec,etapath[1,])
    Psi_t = out_t[[2]]
    G_t = out_t[[3]]

    stoch_mean = G_t%*%Psi_t%*%t(G_end)%*%G%*%solve(t(G)%*%G_end%*%Psi_end%*%t(G_end)%*%G+sigma^2)%*%(y-t(G)%*%etapath[m+1,])
    mpath_c[i+1,] = stoch_mean
  }  

  for(j in 1:it)
  {
    prop=RBminuspath(x0,y,inter,cvec,sigma,m,G,mpath_c,etapath)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikRBminuspn(cpath-etapath-mpath_c,cvec,dt,G,sigma,y,cpath,etapath,mpath_c)-llikRBminuspn(prop-etapath-mpath_c,cvec,dt,G,sigma,y,prop,etapath,mpath_c)+dnorm(y,prop[m+1,1],sigma,log=T)-dnorm(y,cpath[m+1,1],sigma,log=T)}

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
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(nl,nu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*N[t]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(cl,cu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*C[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(nl,nu),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(cl,cu),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}

##GP
GPMH=function(x0,y,inter,cvec,sigma,m,G=matrix(c(1,0),ncol=1,nrow=2),it,plot=FALSE,nl,nu,cl,cu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],y,length.out=m+1)	
  cpath[,2]=rep(x0[2],m+1)

  for(j in 1:it)
  {
    prop=GPpath(x0,y,inter,cvec,sigma,m,G)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikGPpn(cpath,cvec,dt,G,sigma,y)-llikGPpn(prop,cvec,dt,G,sigma,y)+dnorm(y,prop[m+1,1],sigma,log=T)-dnorm(y,cpath[m+1,1],sigma,log=T)}
 
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
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(nl,nu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*N[t]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(cl,cu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*C[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(nl,nu),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(cl,cu),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}

##GP-MDB
GPMDBMH=function(x0,y,inter,cvec,sigma,m,G=matrix(c(1,0),ncol=1,nrow=2),it,plot=FALSE,nl,nu,cl,cu)
{
  count=0
  col=rainbow(it)
  dt=inter/m

  cpath=matrix(0,ncol=length(x0),nrow=m+1)  							
  cpath[,1]=seq(x0[1],y,length.out=m+1)	
  cpath[,2]=rep(x0[2],m+1)

  for(j in 1:it)
  {
    prop=GPMDBpath(x0,y,inter,cvec,sigma,m,G)

    if(j==1){lap=1}
    else{lap=llikE(prop,cvec,dt)-llikE(cpath,cvec,dt)+llikGPMDBpn(cpath,cvec,dt,G,sigma,y)-llikGPMDBpn(prop,cvec,dt,G,sigma,y)+dnorm(y,prop[m+1,1],sigma,log=T)-dnorm(y,cpath[m+1,1],sigma,log=T)}
 
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
          plot(ts(prop[,1],start=0,deltat=dt),ylim=c(nl,nu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*N[t]*""), side=2, line=3.1, cex=2,las=2)

          screen(2)
          par(mar=c(5.2,4.9,1.1,1.1))
          par(cex.axis=2)
          plot(ts(prop[,2],start=0,deltat=dt),ylim=c(cl,cu),ylab="",xlab="",main="",col=col[j])
          mtext("Time", side=1, line=3.5, cex=2)
          mtext(expression(""*C[t]*""), side=2, line=3.1, cex=2,las=2)
        }
        else
        {
          screen(1,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,1],start=0,deltat=dt),col=col[j],ylim=c(nl,nu),axes = FALSE,ylab="",xlab="")
				  
          screen(2,F)
          par(mar=c(5.2,4.9,1.1,1.1))
          plot(ts(prop[,2],start=0,deltat=dt),col=col[j],ylim=c(cl,cu),axes = FALSE,ylab="",xlab="")
        } 
      }
    }
  }
  close.screen(all = TRUE) 
  count/it
}
