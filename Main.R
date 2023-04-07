#######################################################################
#######################################################################
args <- commandArgs(trailingOnly=TRUE)

core <- as.numeric(args[1])
#######################################################################
#######################################################################

########### The Lambda Function #######################################
LAMBDA1=function(s){L=s;return(L)};LAMBDA2=function(s){L=3*sqrt(s);return(L)};
#######################################################################

################# Generate the Data Set ###############################
Gen.Data=function(n,Case,LAMBDA0){
  X1=rbinom(n,1,0.5);X2=runif(n);X=cbind(X1,X2);
  if(Case==1||Case==3){
    K=sample(1:6,n,TRUE);Obs=max(K);T=TT=matrix(NA,n,Obs);N=matrix(NA,n,Obs)
    for(i in 1:n){T[i,1:K[i]]=sort(runif(K[i]))*tau;TT[i,1:K[i]]=T[i,1:K[i]]*exp(rep(X[i,]%*%alpha0,K[i]))}
  }
  if(Case==2||Case==4){
    gamma=rgamma(n,4,4);index1=which(gamma>1);index0=which(gamma<=1);
    Temp.K=c(sample(1:8,length(index1),TRUE),sample(1:6,length(index0),TRUE));K=Temp.K;K[index1]=Temp.K[1:length(index1)];K[index0]=Temp.K[(1+length(index1)):n];
    Obs=max(K);T=TT=matrix(NA,n,Obs);N=matrix(NA,n,Obs);for(i in 1:n){T[i,1:K[i]]=sort(runif(K[i]))*tau;TT[i,1:K[i]]=T[i,1:K[i]]*exp(rep(X[i,]%*%alpha0,K[i]))}
  }
  if(Case==1){for(i in 1:n){
    N[i,1]=rpois(1,LAMBDA0(TT[i,1]));if(K[i]!=1){for(j in 2:K[i]){N[i,j]=rpois(1,(LAMBDA0(TT[i,j])-LAMBDA0(TT[i,j-1])))+N[i,j-1]}}
  }}
  if(Case==2){for(i in 1:n){
      N[i,1]=rpois(1,gamma[i]*LAMBDA0(TT[i,1]));if(K[i]!=1){for(j in 2:K[i]){N[i,j]=rpois(1,gamma[i]*(LAMBDA0(TT[i,j])-LAMBDA0(TT[i,j-1])))+N[i,j-1]}}
  }}
  if(Case==3){gamma=rgamma(n,4,4);for(i in 1:n){
    N[i,1]=rpois(1,gamma[i]*LAMBDA0(TT[i,1]));if(K[i]!=1){for(j in 2:K[i]){N[i,j]=rpois(1,gamma[i]*(LAMBDA0(TT[i,j])-LAMBDA0(TT[i,j-1])))+N[i,j-1]}}
  }}
  if(Case==4){for(i in 1:n){
    N[i,1]=rpois(1,gamma[i]*LAMBDA0(TT[i,1]));if(K[i]!=1){for(j in 2:K[i]){N[i,j]=rnbinom(1,gamma[i]*(LAMBDA0(TT[i,j])-LAMBDA0(TT[i,j-1])),0.5)+N[i,j-1]}}
  }}
  return(list(N=N,T=T,K=K,X=X,TT=TT))
}
#######################################################################

#################### Find the knots of splines ########################
find.knots=function(Data,theta){
  alphaknots=theta[1:p];TT=matrix(NA,length(Data$K),max(Data$K))
  for(i in 1:length(Data$K)){TT[i,1:Data$K[i]]=Data$T[i,1:Data$K[i]]*exp(rep(Data$X[i,]%*%alphaknots,Data$K[i]))}
  knots=c(quantile(TT,seq(1/(nknots+1),nknots/(nknots+1),1/(nknots+1)),na.rm=TRUE,names=FALSE))
  if(max(knots)>b){
    write.table(id.seed,file=paste0("core",core,"seed",seed,"bad.txt"),sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
    TT=Data$T;knots=c(quantile(TT,seq(1/(nknots+1),nknots/(nknots+1),1/(nknots+1)),na.rm=TRUE,names=FALSE))
  }
  return(knots)
}
#######################################################################

################## Generate Bootstrap Data ############################
Genbots.Data=function(Data){
  nnn=length(Data$K);id=sample(1:nnn,nnn,rep=T);N=Data$N[id,];T=Data$T[id,];K=Data$K[id];X=Data$X[id,];
  return(list(N=N,T=T,K=K,X=X))
}
#######################################################################

################# Calculate the Loss function #########################
Loss=function(theta){
  N=Data$N;T=Data$T;K=Data$K;nnn=length(Data$K);Obs=max(K);alpha=theta[1:p];xi=theta[(p+1):length(theta)];TT=matrix(NA,length(K),Obs);
  for(i in 1:length(Data$K)){TT[i,1:Data$K[i]]=Data$T[i,1:Data$K[i]]*exp(rep(Data$X[i,]%*%alpha,Data$K[i]))};TT[which(TT>b)]=b;
  II=iSpline(c(0,c(na.omit(as.vector(t(TT)))),b),df=INN,knots=knot,degree=(INN-1),intercept=TRUE)
  loss=matrix(0,nnn,Obs);for(i in 1:nnn){for(j in 1:K[i]){loss[i,j]=(N[i,j]-t(II[(1+sum(c(0,K)[1:i])+j),])%*%xi)^2}}
  ln=sum(loss)/nnn;return(ln)
}
D.Loss=function(theta){
  N=Data$N;T=Data$T;K=Data$K;nnn=length(Data$K);Obs=max(K);alpha=theta[1:p];xi=theta[(p+1):length(theta)];TT=matrix(NA,length(Data$K),max(Data$K));
  for(i in 1:length(Data$K)){TT[i,1:Data$K[i]]=Data$T[i,1:Data$K[i]]*exp(rep(Data$X[i,]%*%alpha,Data$K[i]))};TT[which(TT>b)]=b;
  II=iSpline(c(0,c(na.omit(as.vector(t(TT)))),b),df=INN,knots=knot,degree=(INN-1),intercept=TRUE);
  dII=mSpline(c(0,c(na.omit(as.vector(t(TT)))),b),df=INN,knots=knot,degree=(INN-1),intercept=TRUE);
  d.loss=array(0,c(nnn,Obs,(p+INN+nknots)))
  for(i in 1:nnn){for(j in 1:K[i]){
    d.loss[i,j,1:p]=2*c((t(II[(1+sum(c(0,K)[1:i])+j),])%*%xi-N[i,j])*t(dII[(1+sum(c(0,K)[1:i])+j),])%*%xi*TT[i,j])*c(Data$X[i,])
    d.loss[i,j,(p+1):(p+INN+nknots)]=2*c(t(II[(1+sum(c(0,K)[1:i])+j),])%*%xi-N[i,j])*c(II[(1+sum(c(0,K)[1:i])+j),])
  }}
  dln=apply(d.loss,3,sum)/nnn;return(dln)
}
bots.Loss=function(theta){
  N=bots.Data$N;T=bots.Data$T;K=bots.Data$K;nnn=length(bots.Data$K);Obs=max(K)
  alpha=theta[1:p];xi=theta[(p+1):length(theta)];TT=matrix(NA,length(bots.Data$K),max(bots.Data$K))
  for(i in 1:length(bots.Data$K)){TT[i,1:bots.Data$K[i]]=bots.Data$T[i,1:bots.Data$K[i]]*exp(rep(bots.Data$X[i,]%*%alpha,bots.Data$K[i]))};TT[which(TT>b)]=b;
  II=iSpline(c(0,c(na.omit(as.vector(t(TT)))),b),df=INN,knots=knot,degree=(INN-1),intercept=TRUE)
  loss=matrix(0,nnn,Obs);for(i in 1:nnn){for(j in 1:K[i]){loss[i,j]=(N[i,j]-t(II[(1+sum(c(0,K)[1:i])+j),])%*%xi)^2}}
  ln=sum(loss)/nnn;return(ln)
}
D.bots.Loss=function(theta){
  N=bots.Data$N;T=bots.Data$T;K=bots.Data$K;nnn=length(bots.Data$K);Obs=max(K)
  alpha=theta[1:p];xi=theta[(p+1):length(theta)];TT=matrix(NA,length(bots.Data$K),max(bots.Data$K))
  for(i in 1:length(bots.Data$K)){TT[i,1:bots.Data$K[i]]=bots.Data$T[i,1:bots.Data$K[i]]*exp(rep(bots.Data$X[i,]%*%alpha,bots.Data$K[i]))};TT[which(TT>b)]=b;
  II=iSpline(c(0,c(na.omit(as.vector(t(TT)))),b),df=INN,knots=knot,degree=(INN-1),intercept=TRUE);
  dII=mSpline(c(0,c(na.omit(as.vector(t(TT)))),b),df=INN,knots=knot,degree=(INN-1),intercept=TRUE);
  d.loss=array(0,c(nnn,Obs,(p+INN+nknots)))
  for(i in 1:nnn){for(j in 1:K[i]){
    d.loss[i,j,1:p]=2*c((t(II[(1+sum(c(0,K)[1:i])+j),])%*%xi-N[i,j])*t(dII[(1+sum(c(0,K)[1:i])+j),])%*%xi*TT[i,j])*c(bots.Data$X[i,])
    d.loss[i,j,(p+1):(p+INN+nknots)]=2*c(t(II[(1+sum(c(0,K)[1:i])+j),])%*%xi-N[i,j])*c(II[(1+sum(c(0,K)[1:i])+j),])
  }}
  dln=apply(d.loss,3,sum)/nnn;return(dln)
}
#######################################################################

library(splines2);
Case=1;LAMBDA=LAMBDA1;method=1;n=100; #################### Data Settings Changed as Needed
tau=6;b=200;alpha0=c(1,-0.5);p=length(alpha0);INN=4;nknots=round(2*n^(1/7));no.sim=20;num.bots=100;
theta.hat=matrix(0,no.sim,p+INN+nknots);knots.hat=matrix(0,no.sim,nknots);bots.theta=matrix(0,no.sim*num.bots,p+INN+nknots);
alpha.AMM=matrix(0,no.sim,p);base=matrix(0,no.sim,500)
for(seed in 1:no.sim){
  id.seed=seed+(core-1)*no.sim;set.seed(id.seed);Data=Gen.Data(n,Case,LAMBDA);Initheta=c(rep(0,p),rep(1,INN+nknots));
  if(method==1){
    knot=(1:nknots)/(nknots+1)*max(Data$T,na.rm=T)
    sol=optim(Initheta,Loss,D.Loss,method="L-BFGS-B",lower=c(rep(-Inf,p),rep(0,INN+nknots)),control=list(maxit=1000));
    TS=Data$T;ACC=exp(Data$X%*%sol$par[1:p]);for(i in 1:n){TS[i,]=TS[i,]*ACC[i]};knot=(1:nknots)/(nknots+1)*min(max(TS,na.rm=T),b);Initheta=sol$par
    sol=optim(Initheta,Loss,D.Loss,method="L-BFGS-B",lower=c(rep(-Inf,p),rep(0,INN+nknots)),control=list(maxit=1000));
  }
  if(method==2){
    knot=find.knots(Data,Initheta);
    sol=optim(Initheta,Loss,D.Loss,method="L-BFGS-B",lower=c(rep(-Inf,p),rep(0,INN+nknots)),control=list(maxit=1000));
    knot=find.knots(Data,sol$par);Initheta=sol$par
    sol=optim(Initheta,Loss,D.Loss,method="L-BFGS-B",lower=c(rep(-Inf,p),rep(0,INN+nknots)),control=list(maxit=1000));
  }
  theta.hat[seed,]=sol$par;knots.hat[seed,]=knot;
  for(i.bots in 1:num.bots){
    set.seed(id.seed*num.bots+i.bots+1000);bots.Data=Genbots.Data(Data)
    bots.sol=optim(Initheta,bots.Loss,D.bots.Loss,method="L-BFGS-B",lower=c(rep(-Inf,p),rep(0,INN+nknots)),control=list(maxit=1000))
    bots.theta[((seed-1)*num.bots+i.bots),]=bots.sol$par
  }
  write.table(id.seed,file=paste0("",core,".txt"),sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
}
write.table(theta.hat,file=paste0("core",core,"thetahat.txt"),sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
write.table(knots.hat,file=paste0("core",core,"knots.txt"),sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)
write.table(bots.theta,file=paste0("core",core,"botstheta.txt"),sep="\t",append=TRUE,col.names=FALSE,row.names=FALSE)


