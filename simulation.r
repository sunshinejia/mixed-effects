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
fix.beta<-2
x1<- seq(0.05,by=0.05,length=20)
x2<-seq(0.1,0.8,0.05)
x1<-rep(x1,15)
x2<-rep(x2,each=20)
subject<-rep(seq(1,15,1),each=20)
dt<-matrix(0,nrow=100,ncol=300);dataframe<-list(NA);result<-list(NA)
i<-0
j <- 1
while (i <100){
    set.seed(j)
    beta<- rmvnorm(15, mean, sigma)
    set.seed(j)
    error<-rnorm(300,0,1)
	rand.beta<-cbind(rep(beta[,1],each=20),rep(beta[,2],each=20))
    dat<-rand.beta[,1]+rand.beta[,2]*x1+fix.beta*x2+error
    data<-data.frame(cbind(subject,x1,x2,dat))
    data$subject<-as.factor(data$subject)
    res<-try(lme(dat~x1+x2, data=data, random=~x1|subject,method="ML",control=list(maxIter=10000)),TRUE)
	j <- j+1	
	if (inherits(res, "try-error") != TRUE)
	{
	i<-i+1
	dt[i,]<-dat
	dataframe[[i]]<-data
	result[[i]]<-res
	}
}
yy<-list(NA)
for (j in 1:100){
    yy[[j]]<-matrix(0,nrow=(20+2),ncol=15)
    for(i in 1:15){
	     yy[[j]][1:20,i]<-dt[j,((i-1)*20+1): (i*20)]
                  }
				 }
Xmat<-cbind(rep(1,300),x1,x2)
Xmat1<-matrix(NA, nrow=20,ncol=45)
for (i in 1:15){
     Xmat1[,(3*i-2):(i*3)]<-Xmat[(20*i-19):(20*i),]
	}
xx<-list(NA)
for( i in 1:15){
    xx[[i]]<-rbind(Xmat1[,(3*i-2):(i*3)],matrix(rep(0,6),nrow=2,ncol=3))
	           }
Zmat<-cbind(rep(1,20),x1[1:20])


#replace one of the linear parameters
uf<-2
x1_test<-0.05
x2_test<-0.1
full_loglike3<-function(y,ps){
    #browser()
    error<-ps[6]
	k<-qnorm(1-ps[7])
	delta<-matrix(c(ps[3:4],0,ps[5]),nrow=2,ncol=2,byrow=T)
	sigma0<-error^2*(ps[4]^2+ps[5]^2)/(ps[3]*ps[5])
	sigma1<-error^2/ps[5]^2
	covar<-error^2*ps[4]^2/(ps[3]*ps[5]^2)
	beta0<-uf-ps[1]*x1_test-ps[2]*x2_test-k*sqrt(sigma0+sigma1*x1_test^2+2*x1_test*covar)
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
		b<-t(Qi)%*%y[,i]
		b<-as.vector(b)
		tmp[((i-1)*20+1):(i*20),4]<-b[3:22]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
    res<- -15*10*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-res)
}
max_tran_linear<-list(NA)
for (i in 1:100) max_tran_linear[[i]]<-nlminb(c(2,2,3,-1,1,1,0.6),full_loglike3,y=yy[[i]],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001,0.001),upper=c(rep(Inf,6),0.999))

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

prof_loglike3<-function(theta,y,ps){
    #browser()
	error<-ps[6]
	k<-qnorm(1-theta)	
	delta<-matrix(c(ps[3:4],0,ps[5]),nrow=2,ncol=2,byrow=T)
	sigma0<-error^2*(ps[4]^2+ps[5]^2)/(ps[3]*ps[5])
	sigma1<-error^2/ps[5]^2
	covar<-error^2*ps[4]^2/(ps[3]*ps[5]^2)
	beta0<-uf-ps[1]*x1_test-ps[2]*x2_test-k*sqrt(sigma0+sigma1*x1_test^2+2*x1_test*covar)
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
		b<-t(Qi)%*%y[,i]
		b<-as.vector(b)
		tmp[((i-1)*20+1):(i*20),4]<-b[3:22]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
    result<- -15*10*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-result)
}

standard<-exp(-qchisq(0.95,1)/2)
norm_value.l<-NULL;opt_in.l<-list(NA);norm_value.r<-NULL;opt_in.r<-list(NA);star<-matrix(NA,ncol=7,nrow=100)
for (i in 1:100){
    star[i,]<-max_tran_linear[[i]]$par
	if (star[i,7]<=Fx){
        opt_in.r[[i]]<-nlminb(start=star[i,-7],prof_loglike3,y=yy[[i]],theta=Fx,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
	    norm_value.r[i]<-exp(-opt_in.r[[i]]$objective+max_tran_linear[[i]]$objective)
	            } else {
		opt_in.l[[i]]<-nlminb(start=star[i,-7],prof_loglike3,y=yy[[i]],theta=Fx,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
	    norm_value.l[i]<-exp(-opt_in.l[[i]]$objective+max_tran_linear[[i]]$objective)
                       }
				} 
num.l<-sum(norm_value.l<standard,na.rm=T) #the true value falls in the left tail; ans:1
num.r<-sum(norm_value.r<standard,na.rm=T) #the true value falls in the right tail; ans:3
			
####test to see whether the one-time optim and continuous optim agree
stars.t<-matrix(NA,nrow=100,ncol=20);opt_test<-list(NA);yhat.t<-list(NA);res_test<-list(NA);pars.t<-matrix(NA,nrow=100,ncol=6)
res_test<-list(NA);obj.t<-NULL
for (i in 1:100){
     stars.t[i,] <- seq(from=star[i,7],to=Fx,length.out=20)
     opt_test[[i]]<-nlminb(start=star[i,-7],prof_loglike3,y=yy[[i]],theta=stars.t[i,2],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
     yhat.t[[i]]<-c(1,exp(-opt_test[[i]]$objective+max_tran_linear[[i]]$objective))
	 pars.t[i,]<-opt_test[[i]]$par
     for (st in stars.t[i,c(-1,-2)]){
           res_test[[i]]<-nlminb(start=pars.t[i,],prof_loglike3,y=yy[[i]],theta=st,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
           obj.t[i] <- res_test[[i]]$objective
           yhat.t[[i]] <-c(yhat.t[[i]],exp(-obj.t[i]+max_tran_linear[[i]]$objective))
           pars.t[i,] <- res_test[[i]]$par
                                  }
 }
yhat_true<-NULL
for (i in 1:100) yhat_true[i]<-yhat.t[[i]][20]
