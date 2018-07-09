VBGF<-function(Linf,k,t0,ages)
{
Lengths_exp<-Linf*(1-exp(-k*(ages-t0)))
return(Lengths_exp)
}

VBGF.fit<-function(p,obs,return.type=2)
{
    if(return.type==1)
    {
        exp.lts<-VBGF(p[1],p[2],p[3],obs[,1])
        return(exp.lts)
    }
    if(return.type==2)
    {
        exp.lts<-VBGF(p[1],p[2],p[3],obs[,1])
        sigma<-sqrt(sum((obs[,2]-exp.lts)^2)/length(obs[,2]))
        neglogsum<-sum(-log(dnorm(exp.lts,obs[,2],sigma)))
        return(neglogsum)
    }
}

VBGF.fit.fxn<-function(AG.dat,plot.it=F)
{
  vbgf.out<-nls(lts~Linf*(1-exp(-k*(ages-t0))),data=list(ages=AG.dat[,1],lts=AG.dat[,2]),start=list(Linf=max(AG.dat[,2]),k=0.1,t0=0),control = list(reltol=0.00000000001))
  vbgf.summ<-summary(vbgf.out)
  if(plot.it==T)
    {
      vbgf.fit.vals<-VBGF(vbgf.summ$par[,1][1],vbgf.summ$par[,1][2],vbgf.summ$par[,1][3],AG.dat[,1])
      plot(AG.dat[,1],vbgf.fit.vals, type="l",lwd=2,xlab="Age",ylab="Length")
      points(AG.dat[,1],AG.dat[,2],pch=21,col="black",bg="darkgreen")
    }
  return(vbgf.summ)
}
