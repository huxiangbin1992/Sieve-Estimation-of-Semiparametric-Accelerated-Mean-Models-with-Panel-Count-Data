##############################################################
### This is the simulation program for 
### "Sieve Estimation of Semiparametric Accelerated Mean Models with Panel Count Data" 
### published in Electronic Journal of Statistics.
### The program contains four data generation cases, two mean functions, two knots selection methods.
### The default setting is sample size n=100, linear mean function, data generated in Case 1 and knots selected by Method 1.
### Readers need to make changes as marked in this program field and "Main.R" for other settings.
##############################################################

################ RUNNING SIMULATION IN PARALLEL ##############
library(parallel)

workerFunc <- function(core){
  
  system(paste("Rscript --vanilla Main.R", core))
}


combinX=list()
no_cores <- 25


for(i in 1:no_cores)
{
  combinX[[i]]<- i
}

c1 <-makeCluster(no_cores)

parLapply(c1,combinX,workerFunc)

stopCluster(c1)
##############################################################

############# Record the Simulation Results ##################
library(splines2)
n=100;LAMBDA=1; #################### Data Settings Changed as Needed
no.sim=20;num.bots=100;INN=4;nknots=round(2*n^(1/7));p=2;l=INN+nknots+p;b=200
theta.hat=matrix(0,no.sim*no_cores,l);knots.hat=matrix(0,no.sim*no_cores,nknots)
bots.theta=array(0,c(no.sim*no_cores,num.bots,l));
for(i_core in 1:no_cores)
{
  theta.hat[(1+(i_core-1)*no.sim):(no.sim+(i_core-1)*no.sim),]=as.matrix(read.table(file=paste0("core",i_core,"thetahat.txt")),no.sim,l)
  knots.hat[(1+(i_core-1)*no.sim):(no.sim+(i_core-1)*no.sim),]=as.matrix(read.table(file=paste0("core",i_core,"knots.txt")),no.sim,nknots)
  temp.bot=as.matrix(read.table(file=paste0("core",i_core,"botstheta.txt")),no.sim*num.bots,l)
  for(i in 1:no.sim){bots.theta[((i_core-1)*no.sim+i),,]=temp.bot[((i-1)*num.bots+1):(i*num.bots),]}
}
bots.sd=apply(bots.theta,c(1,2),sd)
x=c(0:3000/3000*16.31,b);y1=matrix(0,no.sim*no_cores,length(x));
if(LAMBDA==1){y=x};if(LAMBDA==2){y=3*sqrt(x)};
for(i in 1:(no.sim*no_cores)){
  II=iSpline(x,df=INN,knots=knots.hat[i,],degree=(INN-1),intercept=TRUE);
  y1[i,]=c(II%*%theta.hat[i,(p+1):(p+INN+nknots)])
}
save.image("Result.RData")
##############################################################

############# Analysis of the Simulation Results #############
par(mar=c(4.5,4.5,1,2));y1.up=y1.dn=rep(NA,length(x))
for(i in 1:length(x)){y1.up[i]=quantile(y1[,i],probs=0.975);y1.dn[i]=quantile(y1[,i],probs=0.025)}
library(grDevices);library(extrafont);
windowsFonts(Times = windowsFont("Times New Roman"));
par(mar=c(4,4,1,1),family="Times");
plot(x[1:(length(x)-1)],y[1:(length(x)-1)],type="l",col="black",xlab="t",ylim=c(0,max(y1.up[1:(length(x)-1)])),ylab=expression(Lambda),cex.axis=1.4,cex.lab=1.4)
lines(x[1:(length(x)-1)],apply(y1,2,mean)[1:(length(x)-1)],type="l",lty=2,col="red")
lines(x[1:(length(x)-1)],y1.up[1:(length(x)-1)],type="l",lty=2,col="red");lines(x[1:(length(x)-1)],y1.dn[1:(length(x)-1)],type="l",lty=2,col="red");
legend("topleft", lty=c(1,2), col=c("black", "red"),legend=c("True Function","Estimate"),cex=1.4)
bots.sd=matrix(NA,no.sim*no_cores,(p+INN+nknots))
for(i in 1:(no.sim*no_cores)){
  id.temp=unique(c(which(abs(bots.theta[i,,1])>5),which(abs(bots.theta[i,,2])>5)))
  if(length(id.temp)>0){bots.sd[i,]=apply(bots.theta[i,-id.temp,],2,sd)}else{bots.sd[i,]=apply(bots.theta[i,,],2,sd)}
}
id.temp=unique(c(which(abs(theta.hat[,1])>5),which(abs(theta.hat[,2])>5)))
if(length(id.temp)>0){
  apply(theta.hat[-id.temp,],2,mean);apply(theta.hat[-id.temp,],2,sd);apply(bots.sd[-id.temp,],2,mean);apply(theta.hat[-id.temp,],2,mean)[1:2]-c(1,-0.5)
  1-apply(abs(theta.hat[-id.temp,1:2]-matrix(c(1,-0.5),nr=(no.sim*no_cores-length(-id.temp)),nc=p,byrow=T))>1.96*bots.sd[-id.temp,1:2],2,mean)
}else{
  apply(theta.hat,2,mean);apply(theta.hat,2,sd);apply(bots.sd,2,mean);apply(theta.hat,2,mean)[1:2]-c(1,-0.5)
  1-apply(abs(theta.hat[,1:2]-matrix(c(1,-0.5),nr=no.sim*no_cores,nc=p,byrow=T))>1.96*bots.sd[,1:2],2,mean)
}
##############################################################

