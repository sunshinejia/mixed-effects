library(mvtnorm)
library(lme4)
library(nlme)
library(MASS)
library(MSBVAR)
mean<-c(2,2)
sigma<-diag(2)
sigma[1,1] <-1
sigma[1,2] <-sigma[2,1]<-0.6
sigma[2,2] <-1
fix.beta<-2
Fx<-1-pnorm((2-2-2*0.05-2*0.1)/sqrt(1+0.05^2+0.1*0.8))
POD<-pnorm((2+2*0.05+2*0.1-2)/sqrt(1+0.05^2+0.1*0.8+1))
set.seed(100)
beta<- rmvnorm(15, mean, sigma)
set.seed(100)
error<-rnorm(300,0,1)
rand.beta<-matrix(NA,nrow=300,ncol=2)
for (i in 1:2)rand.beta[,i]<-rep(beta[,i],each=20)
x1<- seq(0.05,by=0.05,length=20)
x2<-seq(0.1,0.8,0.05)
x1<-rep(x1,15)
x2<-rep(x2,each=20)
dat<-rand.beta[,1]+rand.beta[,2]*x1+fix.beta*x2+error
subject<-rep(seq(1,15,1),each=20)
list<-data.frame(cbind(subject,x1,x2,dat))
list$subject<-as.factor(list$subject)
res<-lme(dat~x1+x2, data=list, random=~x1|subject,method="ML")

#MCMC sample 
library(bayesm)
V<-3*matrix(c(1,0,0,0,1,0,0,0,0),byrow=T,ncol=3)
regdata[[1]]<-list(dat,cbind(rep(1,300),x1,x2))
Z=c(1,100)
data1<-list(regdata=regdata,Z=Z)
out<-rhierLinearModel(Data=data1,Mcmc=list(R=5000,keep=1000),Prior=list(V=V))
library(lmm)
xcol<-1:3
zcol<-1:2
fmcmc.result<-fastmcmc.lmm(y=dat,subj=subject,pred=cbind(rep(1,300),list$x1,list$x2),xcol,zcol,prior=list(a=3*var(dat),b=3,c=2,Dinv=2*diag(2)),seed=2345,iter=10000)
gibbs.result<-mgibbs.lmm(y=dat,subj=subject,pred=cbind(rep(1,300),list$x1,list$x2),xcol,zcol,prior=list(a=3*var(dat),b=3,c=2,Dinv=2*diag(2)),seed=2345,iter=10000)
plot(gibbs.result$sigma2.series)
plot(fmcmc.result$sigma2.series)

par(mfrow=c(4,1))
acf(log(gibbs.result$psi.series[1,1,]),lag.max=10, ylim=0:1)
acf(log(fmcmc.result$psi.series[1,1,]),lag.max=10, ylim=0:1)
acf(gibbs.result$psi.series[1,2,],lag.max=10)
acf(fmcmc.result$psi.series[1,2,],lag.max=10)

n<-300
m<-15
p<-2
q<-1
BETA<-NULL
list$newdat<-list$dat-list$x2*res$coefficients$fixed[3]
for (j in 1:m){ BETA<-rbind(BETA,lm(newdat[subject==j]~x1[subject==j],data=list)$coef)}
mu0<-apply(BETA,2,mean)
S0<-cov(BETA)
eta0<-p+2
iL0<-iSigma<-solve(S0)
nu0<-3
sigma0<-var(dat)
S<-5000
beta.post<-list(NA);sigma2.post<-NULL;Sigma.post<-list(NA)
set.seed(1)
for(s in 1:5000){
    #update theta
	Lm<-solve(iL0+m*iSigma)
	mum<-Lm%*%(iL0%*%mu0+iSigma%*%apply(BETA,2,sum))
	theta<-t(rmvnorm(1,mum,Lm))
	#updata Sigma
	mtheta<-matrix(theta,m,p,byrow=T)
	iSigma<-rwishart(1,eta0+m,solve(S0+t(BETA-mtheta)%*%(BETA-mtheta)))
	#update sigma
	rand.BETA<-matrix(NA,nrow=300,ncol=2)
    for (i in 1:2)rand.BETA[,i]<-rep(BETA[,i],each=20)
	SSR<-sum((dat-rand.BETA[1]-rand.BETA[2]*x1)^2)
	sigma2<-1/rgamma(1,(n+nu0)/2,(nu0*sigma0+SSR)/2)
    #updata beta
	Sigma<-solve(iSigma)
	BETA<-rmvnorm(m,theta,Sigma)
	beta.post[[s]]<-BETA
    Sigma.post[[s]]<-Sigma 
    sigma2.post[s]<-sigma2	
	 }

quantile(sigma2.post,probs=c(0.025,0.975))
mean(sigma2.post)

############################try to simulate a different sample
x1<- seq(from=-1,to=1,length=20)
x2<-seq(from=-1,to=1,length=15)
x1<-rep(x1,15)
x2<-rep(x2,each=20)
subject<-rep(seq(1,15,1),each=20)
dt<-matrix(0,nrow=5000,ncol=300);dataframe<-list(NA);result<-list(NA)
for (i in 1:5000){
     set.seed(i)
     beta<- rmvnorm(15, mean, sigma)
     set.seed(i)
     error<-rnorm(300,0,1)
     rand.beta<-matrix(NA,nrow=300,ncol=2)
     for (j in 1:2)rand.beta[,j]<-rep(beta[,j],each=20)
     dat<-rand.beta[,1]+rand.beta[,2]*x1+fix.beta*x2+error
	 list<-data.frame(cbind(subject,x1,x2,dat))
     list$subject<-as.factor(list$subject)
	 res<-lme(dat~x1+x2, data=list, random=~x1|subject,method="ML")
     #res<-lme(dat~x1+x2, data=list, random=~x1|subject,method="ML",control=list(niterEM=2000))
	 dt[i,]<-dat
	 dataframe[[i]]<-data
	 result[[i]]<-res
	 rho[i]<-res
	 }
###########################################################test the intercept and slope for each subject when lme does not converge
x1_t<- seq(from=-1,to=1,length=20)
x2_t<-seq(from=-1,to=1,length=15)
x1_t<-rep(x1_t,15)
x2_t<-rep(x2_t,each=20)
subject_t<-rep(seq(1,15,1),each=20)
set.seed(863)
beta_t<- rmvnorm(15, mean, sigma)
set.seed(863)
error_t<-rnorm(300,0,1)
rand.beta_t<-matrix(NA,nrow=300,ncol=2)
for (j in 1:2)rand.beta_t[,j]<-rep(beta_t[,j],each=20)
dat_t<-rand.beta_t[,1]+rand.beta_t[,2]*x1_t+fix.beta*x2_t+error_t
list_t<-data.frame(cbind(subject_t,x1_t,x2_t,dat_t))
list_t$subject_t<-as.factor(list_t$subject_t)
res_t<-lme(dat_t~x1_t+x2_t, data=list_t, random=~x1_t|subject_t,method="ML") #does not work
res_t<-lme(dat_t~x1_t+x2_t, data=list_t, random=~x1_t|subject_t,method="ML",control=list(niterEM=2000)) 

#compare data using plots
par(mfrow=c(1,2))
plot(subject_t,dat_t,xlab="subject",ylab="data",main="lme does not converge")
plot(subject,dat,main="lme converges")
library(ggplot2)
par(mfrow=c(1,2))
plot(x=x1,y=dat,main="lme converges")
plot(x1_t,dat_t,xlab="x1",ylab="dat",main="lme does not converge")
par(mfrow=c(1,2))
hist(dat,main="lme converges")
hist(dat_t,xlab="dat",main="lme does not converge")
qplot(x=x1,y=dat,size=x2,main="lme converges")
qplot(x=x1_t,y=dat_t,size=x2_t,main="lme does not converge")

res1<-lmer(dat_t~1+x1_t+x2_t+(1+x1_t|subject_t),list_t)
b0<-NULL;b1<-NULL
newdat<-dat-x2*result$coefficients$fixed[3]
for (i in 1:15){
   res<-lm(newdat[(1:20)*i]~x1[(1:20)*i])
   b0<-c(b0,res$coefficients [1])
   b1<-c(b1,res$coefficients [2])
              }
mub0<-mean(b0)
mub1<-mean(b1)
sig<-cor(b0,b1)

set.seed(49)
beta<- rmvnorm(15, mean, sigma)
set.seed(49)
error<-rnorm(300,0,1)
rand.beta<-matrix(NA,nrow=300,ncol=2)
for (j in 1:2)rand.beta[,j]<-rep(beta[,j],each=20)
dat<-rand.beta[,1]+rand.beta[,2]*x1+fix.beta*x2+error
list<-data.frame(cbind(subject,x1,x2,dat))
list$subject<-as.factor(list$subject)
B0<-NULL;B1<-NULL
for (i in 1:15){
   Res<-lm(dat[(1:20)*i]~x1[(1:20)*i])
   B0<-c(b0,res$coefficients [1])
   B1<-c(b1,res$coefficients [2])
              }
Mub0<-mean(B0)
Mub1<-mean(B1)
Sig<-cor(B0,B1)
####################################################

yy<-matrix(0,nrow=(20+2),ncol=15)
for(i in 1:15){
yy[1:20,i]<-dat[((i-1)*20+1): (i*20)]
}
Xmat<-cbind(rep(1,300),list$x1,list$x2)
Xmat1<-matrix(NA, nrow=20,ncol=45)
for (i in 1:15){
     Xmat1[,(3*i-2):(i*3)]<-Xmat[(20*i-19):(20*i),]
	}
xx<-list(NA)
for( i in 1:15){
    xx[[i]]<-rbind(Xmat1[,(3*i-2):(i*3)],matrix(rep(0,6),nrow=2,ncol=3))
	           }
Zmat<-cbind(rep(1,20),list$x1[1:20])


#replace one of the linear parameters
uf<-2
x1<-0.05
x2<-0.1
full_loglike3<-function(ps){
    #browser()
    error<-ps[6]
	k<-qnorm(1-ps[7])
	delta<-matrix(c(ps[3:4],0,ps[5]),nrow=2,ncol=2,byrow=T)
	sigma0<-error^2*(ps[4]^2+ps[5]^2)/(ps[3]*ps[5])
	sigma1<-error^2/ps[5]^2
	covar<-error^2*ps[4]^2/(ps[3]*ps[5]^2)
	beta0<-uf-ps[1]*x1-ps[2]*x2-k*sqrt(sigma0+sigma1*x1^2+2*x1*covar)
	beta<-c(beta0,ps[1],ps[2])
	zz<-rbind(Zmat,delta)
	zmat_qr<-qr(zz)
    Ri<-qr.R(zmat_qr,complete=TRUE)
    Qi<-qr.Q(zmat_qr,complete=TRUE)
    R11<-Ri[1:2,]
    tmp<-matrix(0,nrow=15*20,ncol=4)
	a<-list(NA);b<-NULL
    for(i in 1:15){
		a[[i]]<-t(Qi)%*%xx[[i]]
		tmp[((i-1)*20+1):(i*20),1:3]<-a[[i]][3:22,]
		b<-t(Qi)%*%yy[,i]
		b<-as.vector(b)
		tmp[((i-1)*20+1):(i*20),4]<-b[3:22]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
    res<- -15*10*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-res)
}
max_tran_linear<-nlminb(c(2,2,3,-1,1,1,0.7),full_loglike3,lower=c(rep(-Inf,5),0.001,0.001),upper=c(rep(Inf,6),0.999))
est<-max_tran_linear$par
delta<-matrix(c(est[3:4],0,est[5]),nrow=2,ncol=2,byrow=T)
Dinv<-t(delta)%*%delta
Sigmahat<-solve(Dinv)*est[6]^2
sdev<-sqrt(diag(Sigmahat))
tmp<-diag(1/sdev)
corrm<- tmp%*%Sigmahat%*%tmp
sigma0<-est[6]^2*(est[4]^2+est[5]^2)/(est[3]*est[5])
sigma1<-est[6]^2/est[5]^2
covar<-est[6]^2*est[4]^2/(est[3]*est[5]^2)
beta0<-uf-est[1]*x1-est[2]*x2-qnorm(1-est[7])*sqrt(sigma0+sigma1*x1^2+2*x1*covar)
#beta0<-uf-est[1]*x1-est[2]*x2-qnorm(1-est[7])*sqrt(sdev[1]^2+sdev[2]^2*x1^2+2*x1*corrm[1,2]*sdev[1]*sdev[2])

prof_loglike3<-function(theta,ps){
    #browser()
	error<-ps[6]
	k<-qnorm(1-theta)	
	delta<-matrix(c(ps[3:4],0,ps[5]),nrow=2,ncol=2,byrow=T)
	sigma0<-error^2*(ps[4]^2+ps[5]^2)/(ps[3]*ps[5])
	sigma1<-error^2/ps[5]^2
	covar<-error^2*ps[4]^2/(ps[3]*ps[5]^2)
	beta0<-uf-ps[1]*x1-ps[2]*x2-k*sqrt(sigma0+sigma1*x1^2+2*x1*covar)
	beta<-c(beta0,ps[1],ps[2])
	zz<-rbind(Zmat,delta)
	zmat_qr<-qr(zz)
    Ri<-qr.R(zmat_qr,complete=TRUE)
    Qi<-qr.Q(zmat_qr,complete=TRUE)
    R11<-Ri[1:2,]
    tmp<-matrix(0,nrow=15*20,ncol=4)
	a<-list(NA);b<-NULL
    for(i in 1:15){
		a[[i]]<-t(Qi)%*%xx[[i]]
		tmp[((i-1)*20+1):(i*20),1:3]<-a[[i]][3:22,]
		b<-t(Qi)%*%yy[,i]
		b<-as.vector(b)
		tmp[((i-1)*20+1):(i*20),4]<-b[3:22]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
    result<- -15*10*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-result)
}


star<-max_tran_linear$par
stars.l <- seq(from=star[7],to=0.001,length.out=10)
opt_in<-nlminb(start=star[-7],prof_loglike3,theta=stars.l[2],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
yhat.l<-c(1,exp(-opt_in$objective+max_tran_linear$objective))
pars<-opt_in$par
for (st in stars.l[c(-1,-2)]){
  res_opt<-nlminb(start=pars,prof_loglike3,theta=st,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
    obj <- res_opt$objective
  yhat.l <-c(yhat.l,exp(-obj+max_tran_linear$objective))
  pars <- res_opt$par
}

stars.r <- seq(from=star[7],to=0.999,length.out=11)
#opt_in<-optim(star[-7],fn=prof_loglike2,theta=stars.r[2],lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
opt_in<-nlminb(start=star[-7],prof_loglike3,theta=stars.r[2],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
yhat.r<-c(1,exp(-opt_in$objective+max_tran_linear$objective))
pars<-opt_in$par
lst <- list()
lst <- c(opt_in)
for (st in stars.r[c(-1,-2)]){
  #res_opt<-optim(pars[-7],fn=prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
  res_opt<-nlminb(start=pars,prof_loglike3,theta=st,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
  lst <- c(lst,res_opt)
  obj <- res_opt$objective
  #if(exp(-res+max_tran2$objective)<0.5) browser()
  yhat.r <-c(yhat.r,exp(-obj+max_tran_linear$objective))
  pars <- res_opt$par
}
plot(c(stars.l,stars.r[-1]),c(yhat.l,yhat.r[-1]),main="20 points in total",xlab="F(x)",ylab="norm likelihood")
abline(h=exp(-qchisq(0.95,1)/2))  
lines(c(rev(stars.l),stars.r[-1]),c(rev(yhat.l),yhat.r[-1]))

####test to see whether the one-time optim and continuous optim agree;result: agree!
stars.t <- seq(from=star[7],to=Fx,length.out=20)
opt_test<-nlminb(start=star[-7],prof_loglike3,theta=stars.t[2],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
yhat.t<-c(1,exp(-opt_test$objective+max_tran_linear$objective))
pars.t<-opt_test$par
lst.t <- list()
lst.t <- c(opt_test)
for (st in stars.t[c(-1,-2)]){
  res_test<-nlminb(start=pars.t,prof_loglike3,theta=st,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
  lst.t <- c(lst.t,res_test)
  obj.t <- res_test$objective
  yhat.t <-c(yhat.t,exp(-obj.t+max_tran_linear$objective))
  pars.t <- res_test$par
}
opt.t<-nlminb(start=star[-7],prof_loglike3,theta=Fx,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
yhat.tt<-exp(-opt.t$objective+max_tran_linear$objective)


###############################################################
#probability of detection
yth<-2
x1<-0.05
x2<-0.1
full_loglike<-function(ps){
    #browser()
    error<-ps[6]
	k<-qnorm(ps[7])
	delta<-matrix(c(ps[3:4],0,ps[5]),nrow=2,ncol=2,byrow=T)
	sigma0<-error^2*(ps[4]^2+ps[5]^2)/(ps[3]*ps[5])
	sigma1<-error^2/ps[5]^2
	rhosigma0<-error*ps[4]^2/(ps[3]^2*ps[5])
	beta0<-k*sqrt(sigma0+sigma1*x1^2+2*x1*rhosigma0+error^2)+yth-ps[1]*x1-ps[2]*x2
	beta<-c(beta0,ps[1],ps[2])
	zz<-rbind(Zmat,delta)
	zmat_qr<-qr(zz)
    Ri<-qr.R(zmat_qr,complete=TRUE)
    Qi<-qr.Q(zmat_qr,complete=TRUE)
    R11<-Ri[1:2,]
    tmp<-matrix(0,nrow=15*20,ncol=4)
	a<-list(NA);b<-NULL
    for(i in 1:15){
		a[[i]]<-t(Qi)%*%xx[[i]]
		tmp[((i-1)*20+1):(i*20),1:3]<-a[[i]][3:22,]
		b<-t(Qi)%*%yy[,i]
		b<-as.vector(b)
		tmp[((i-1)*20+1):(i*20),4]<-b[3:22]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
    res<- -15*10*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-res)
}
max_tran_linear<-nlminb(c(2,2,3,-1,1,1,0.6),full_loglike,lower=c(rep(-Inf,5),0.001,0.001),upper=c(rep(Inf,6),0.999))
est<-max_tran_linear$par
delta<-matrix(c(est[3:4],0,est[5]),nrow=2,ncol=2,byrow=T)
Dinv<-t(delta)%*%delta
Sigmahat<-solve(Dinv)*est[6]^2
sdev<-sqrt(diag(Sigmahat))
tmp<-diag(1/sdev)
corrm<- tmp%*%Sigmahat%*%tmp
sigma0<-est[6]^2*(est[4]^2+est[5]^2)/(est[3]*est[5])
sigma1<-est[6]^2/est[5]^2
covar<-est[6]^2*est[4]^2/(est[3]*est[5]^2)
beta0<-uf-est[1]*x1-est[2]*x2-qnorm(1-est[7])*sqrt(sigma0+sigma1*x1^2+2*x1*covar)
#beta0<-uf-est[1]*x1-est[2]*x2-qnorm(1-est[7])*sqrt(sdev[1]^2+sdev[2]^2*x1^2+2*x1*corrm[1,2]*sdev[1]*sdev[2])

prof_loglike<-function(theta,ps){
    #browser()
	error<-ps[6]
	k<-qnorm(theta)	
	delta<-matrix(c(ps[3:4],0,ps[5]),nrow=2,ncol=2,byrow=T)
	sigma0<-error^2*(ps[4]^2+ps[5]^2)/(ps[3]*ps[5])
	sigma1<-error^2/ps[5]^2
	rhosigma0<-error*ps[4]^2/(ps[3]^2*ps[5])
	beta0<-k*sqrt(sigma0+sigma1*x1^2+2*x1*rhosigma0+error^2)+yth-ps[1]*x1-ps[2]*x2
	beta<-c(beta0,ps[1],ps[2])
	zz<-rbind(Zmat,delta)
	zmat_qr<-qr(zz)
    Ri<-qr.R(zmat_qr,complete=TRUE)
    Qi<-qr.Q(zmat_qr,complete=TRUE)
    R11<-Ri[1:2,]
    tmp<-matrix(0,nrow=15*20,ncol=4)
	a<-list(NA);b<-NULL
    for(i in 1:15){
		a[[i]]<-t(Qi)%*%xx[[i]]
		tmp[((i-1)*20+1):(i*20),1:3]<-a[[i]][3:22,]
		b<-t(Qi)%*%yy[,i]
		b<-as.vector(b)
		tmp[((i-1)*20+1):(i*20),4]<-b[3:22]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
    result<- -15*10*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-result)
}


star<-max_tran_linear$par
stars.l <- seq(from=star[7],to=0.001,length.out=10)
opt_in<-nlminb(start=star[-7],prof_loglike,theta=stars.l[2],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
yhat.l<-c(1,exp(-opt_in$objective+max_tran_linear$objective))
pars<-opt_in$par
for (st in stars.l[c(-1,-2)]){
  res_opt<-nlminb(start=pars,prof_loglike,theta=st,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
    obj <- res_opt$objective
  yhat.l <-c(yhat.l,exp(-obj+max_tran_linear$objective))
  pars <- res_opt$par
}

stars.r <- seq(from=star[7],to=0.999,length.out=11)
#opt_in<-optim(star[-7],fn=prof_loglike2,theta=stars.r[2],lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
opt_in<-nlminb(start=star[-7],prof_loglike,theta=stars.r[2],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
yhat.r<-c(1,exp(-opt_in$objective+max_tran_linear$objective))
pars<-opt_in$par
lst <- list()
lst <- c(opt_in)
for (st in stars.r[c(-1,-2)]){
  #res_opt<-optim(pars[-7],fn=prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
  res_opt<-nlminb(start=pars,prof_loglike,theta=st,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
  lst <- c(lst,res_opt)
  obj <- res_opt$objective
  #if(exp(-res+max_tran2$objective)<0.5) browser()
  yhat.r <-c(yhat.r,exp(-obj+max_tran_linear$objective))
  pars <- res_opt$par
}
plot(c(stars.l,stars.r[-1]),c(yhat.l,yhat.r[-1]),main="20 points in total",xlab="F(x)",ylab="norm likelihood")
abline(h=exp(-qchisq(0.95,1)/2))  
lines(c(rev(stars.l),stars.r[-1]),c(rev(yhat.l),yhat.r[-1]))
