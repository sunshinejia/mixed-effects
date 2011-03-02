library(mvtnorm)
library(lme4)
library(nlme)
library(MASS)
library(Rsolnp)
mean<-c(2,2)
sigma<-diag(2)
sigma[1,1] <-1
sigma[1,2] <-sigma[2,1]<-0.8
sigma[2,2] <-1
Fx<-1-pnorm((2-2-2*0.05-2*0.1)/sqrt(1+0.05^2+0.1*0.8))
set.seed(100)
beta<- rmvnorm(15, mean, sigma)
set.seed(100)
error<-rnorm(300,0,1)
rand.beta<-matrix(NA,nrow=300,ncol=2)
for (i in 1:2)rand.beta[,i]<-rep(beta[,i],each=20)
fix.beta<-2
x1<- seq(0.05,by=0.05,length=20)
x2<-seq(0.1,0.8,0.05)
x1<-rep(x1,15)
x2<-rep(x2,each=20)
dat<-rand.beta[,1]+rand.beta[,2]*x1+fix.beta*x2+error
subject<-rep(seq(1,15,1),each=20)
list<-data.frame(cbind(subject,x1,x2,dat))
list$subject<-as.factor(list$subject)
res<-lme(dat~x1+x2, data=list, random=~x1|subject,method="ML")

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