library(mvtnorm)
library(lme4)
library(nlme)
library(MASS)
library(Rsolnp)
library(bayesm)
mean<-c(2,2)
sigma<-diag(2)
sigma[1,1] <-1
sigma[1,2] <-sigma[2,1]<-0.6
sigma[2,2] <-1
Fx<-1-pnorm((2-2-2*0.05-2*0.1)/sqrt(1+0.05^2+0.1*0.8))
fix.beta<-2
x1<- seq(from=-1,to=1,length=20)
x2<-seq(from=-1,to=1,length=15)
x1<-rep(x1,15)
x2<-rep(x2,each=20)
subject<-rep(seq(1,15,1),each=20)
Z<-cbind(rep(1,20),x1[1:20])
n<-50
dt<-matrix(0,nrow=n,ncol=300);dataframe<-list(NA);result<-list(NA)
i<-0
j <- 1
lst <- list()
k <- 1
while (i <n){
    set.seed(j)
    beta<- rmvnorm(15, mean, sigma)
    set.seed(j)
    error<-rnorm(300,0,1)
	rand.beta<-cbind(rep(beta[,1],each=20),rep(beta[,2],each=20))
    dat<-rand.beta[,1]+rand.beta[,2]*x1+fix.beta*x2+error
    data<-data.frame(cbind(subject,x1,x2,dat))
    data$subject<-as.factor(data$subject)
    res<-try(lme(dat~x1+x2, data=data, random=~x1|subject,method="ML"),TRUE)
	j <- j+1	
	if (£¡inherits(res, "try-error") )
	{
	i<-i+1
	dt[i,]<-dat
	dataframe[[i]]<-data
	result[[i]]<-res
	}else{
	lst[[k]] <- list(seed=j,data=data)
	k <- k + 1
	}
}
#j=2001 when -5 to 5; j=10008; j=2006 when -3 to 3 
yy<-list(NA)
for (j in 1:n){
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
for (i in 1:n) max_tran_linear[[i]]<-nlminb(c(2,2,3,-1,1,1,0.6),full_loglike3,y=yy[[i]],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001,0.001),upper=c(rep(Inf,6),0.999))

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
for (i in 1:n){
    star[i,]<-max_tran_linear[[i]]$par
	if (star[i,7]<=Fx){
        opt_in.r[[i]]<-nlminb(start=star[i,-7],prof_loglike3,y=yy[[i]],theta=Fx,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
	    norm_value.r <- c(norm_value.r,exp(-opt_in.r[[i]]$objective+max_tran_linear[[i]]$objective))
	            } else {
		opt_in.l[[i]]<-nlminb(start=star[i,-7],prof_loglike3,y=yy[[i]],theta=Fx,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
	    norm_value.l <- c(norm_value.l,exp(-opt_in.l[[i]]$objective+max_tran_linear[[i]]$objective))
                       }
				} 
num.l<-sum(norm_value.l<standard,na.rm=T) #the true value falls in the left tail; ans:0
num.r<-sum(norm_value.r<standard,na.rm=T) #the true value falls in the right tail; ans:6
			
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
#################################################
####test to see whether the one-time optim and two times optim agree
stars.t<-matrix(NA,nrow=n,ncol=3);opt_test<-list(NA);yhat.t<-list(NA);res_test<-list(NA);pars.t<-matrix(NA,nrow=n,ncol=6)
res_test<-list(NA);obj.t<-NULL
for (i in 1:n){
     stars.t[i,] <- seq(from=star[i,7],to=Fx,length.out=3)
     opt_test[[i]]<-nlminb(start=star[i,-7],prof_loglike3,y=yy[[i]],theta=stars.t[i,2],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
	 pars.t[i,]<-opt_test[[i]]$par
     res_test[[i]]<-nlminb(start=pars.t[i,],prof_loglike3,y=yy[[i]],theta=stars.t[i,3],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
     obj.t[i] <- res_test[[i]]$objective
     yhat.t[[i]] <-exp(-obj.t[i]+max_tran_linear[[i]]$objective)
   }

#################################POD
yth<-2
x1<-0.05
x2<-0.1
full_loglike<-function(y,ps){
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
		b<-t(Qi)%*%y[,i]
		b<-as.vector(b)
		tmp[((i-1)*20+1):(i*20),4]<-b[3:22]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
    res<- -15*10*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-res)
}
max_pod<-list(NA)
for (i in 1:n) max_pod[[i]]<-nlminb(c(2,2,3,-1,1,1,0.6),full_loglike,y=yy[[i]],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001,0.001),upper=c(rep(Inf,6),0.999))

prof_loglike<-function(theta,y,ps){
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
		b<-t(Qi)%*%y[,i]
		b<-as.vector(b)
		tmp[((i-1)*20+1):(i*20),4]<-b[3:22]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
    result<- -15*10*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-result)
}

stan_pod<-exp(-qchisq(0.95,1)/2)
norm_pod.l<-NULL;opt_pod.l<-list(NA);norm_pod.r<-NULL;opt_pod.r<-list(NA);start_pod<-matrix(NA,ncol=7,nrow=n)
for (i in 1:n){
    start_pod[i,]<-max_pod[[i]]$par
	if (start_pod[i,7]<=POD){
        opt_pod.r[[i]]<-nlminb(start=start_pod[i,-7],prof_loglike,y=yy[[i]],theta=POD,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
	    norm_pod.r <- c(norm_pod.r,exp(-opt_pod.r[[i]]$objective+max_pod[[i]]$objective))
	            } else {
		opt_pod.l[[i]]<-nlminb(start=start_pod[i,-7],prof_loglike,y=yy[[i]],theta=POD,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
	    norm_pod.l <- c(norm_pod.l,exp(-opt_pod.l[[i]]$objective+max_pod[[i]]$objective))
                       }
				} 
num_pod.l<-sum(norm_pod.l<stan_pod) #the true value falls in the left tail; ans:0
num_pod.r<-sum(norm_pod.r<stan_pod) #the true value falls in the right tail; ans:6

#########################################################################################################################################33
#Bayesian
intb1<-matrix(NA,nrow=20,ncol=2);intb2<-matrix(NA,nrow=20,ncol=2);intb3<-matrix(NA,nrow=20,ncol=2)
intvar1<-matrix(NA,nrow=20,ncol=2);intvar2<-matrix(NA,nrow=20,ncol=2);intcov<-matrix(NA,nrow=20,ncol=2);intsig<-matrix(NA,nrow=20,ncol=2)
for (i in 1:20){
  num<-300
  m<-15
  p<-2
  #BETA<-NULL
  #dataframe[[i]]$newdat<-dt[i,]-x2*result[[i]]$coefficients$fixed[3]
  #for (j in 1:m){ BETA<-rbind(BETA,lm(newdat[subject==j]~x1[subject==j],data=dataframe[[i]])$coef)}
  #mu0<-apply(BETA,2,mean)
  #mu0[3]<-2
  beta<-mu0<-c(2,2,2)
  #S0_b<-cov(BETA)
  S0_b<-matrix(c(1,0,0,1),nrow=2)
  eta0<-p+2
  L0<-100*diag(3)
  iSigma0<-solve(S0_b)
  nu0<-2
  #sigma0<-var(dt[i,])
  sigma0<-1
  S<- 10000
  bj<-matrix(rep(0,30),nrow=15,ncol=2)
  beta.post<-list(NA);sigma2.post<-NULL;Sigma.post<-list(NA);b.post<-list(NA)
  set.seed(as.integer(runif(1,0,1000)))
  #Rprof()
    for(s in 1:S){
	#updata Sigma
	 Sigma<-rwishart(eta0+m,solve(S0_b+t(bj)%*%bj))$IW
	#update sigma
	rand.BETA<-matrix(NA,nrow=300,ncol=2)
    for (l in 1:2)rand.BETA[,l]<-rep((bj[,l]+beta[l]),each=20)
	SSR<-sum((dt[i,]-rand.BETA[,1]-rand.BETA[,2]*x1-beta[3]*x2)^2)
	sigma2<-1/rgamma(1,(n+nu0)/2,(nu0*sigma0+SSR)/2)
    #updata beta
	XV<-matrix(0,nrow=3,ncol=3);Xlist<-matrix(NA,nrow=20,ncol=3);XVy<-matrix(0,nrow=3,ncol=1)
	for (r in 1:m){ 
	       Xlist<-cbind(rep(1,20),dataframe[[i]]$x1[dataframe[[i]]$subject==r],dataframe[[i]]$x2[dataframe[[i]]$subject==r])
		   XV<-t(Xlist)%*%Xlist+XV
		   XVy<-t(Xlist)%*%(dat[dataframe[[i]]$subject==r]-Z%*%bj[r,])+XVy
		       }
    Sigma.m<-solve(solve(L0)+XV/sigma2)
	mu.m<-Sigma.m%*%(solve(L0)%*%mu0+XVy/sigma2)
	beta<-rmvnorm(1,mu.m,Sigma.m)
	#updata b
	Sigma.b<-solve(t(Z)%*%Z/sigma2+solve(Sigma))
	mub<-matrix(NA,nrow=m,ncol=p)
	for (r in 1:m){
		 mub[r,]<-Sigma.b%*%t(Z)%*%(dat[dataframe[[i]]$subject==r]-cbind(rep(1,20),dataframe[[i]]$x1[dataframe[[i]]$subject==r],dataframe[[i]]$x2[dataframe[[i]]$subject==r])%*%t(beta))/sigma2
		 bj[r,]<-rmvnorm(1,mub[r,],Sigma.b)
		          }
	beta.post[[s]]<-beta
    Sigma.post[[s]]<-Sigma 
    sigma2.post[s]<-sigma2
    b.post[[s]]<-bj	
	}
	beta1<-NULL;beta2<-NULL;beta3<-NULL;var1<-NULL;var2<-NULL;covar<-NULL;tempsigma<-matrix(NA,nrow=2,ncol=2);temp<-NULL
	for (j in 2001:10000){
        temp<-beta.post[[j]]
	    beta1<-c(temp[1],beta1)
	    beta2<-c(temp[2],beta2)
	    beta3<-c(temp[3],beta3)
	    tempsigma<-Sigma.post[[j]]
	    var1<-c(tempsigma[1,1],var1)
	    var2<-c(tempsigma[2,2],var2)
	    covar<-c(tempsigma[1,2],covar)
	    }
    #Rprof(NULL)
    intb1[i,]<-quantile(beta1,probs=c(0.025,0.975));intb2[i,]<-quantile(beta2,probs=c(0.025,0.975));intb3[i,]<-quantile(beta3,probs=c(0.025,0.975))
    intvar1[i,]<-quantile(var1,probs=c(0.025,0.975))
	intvar2[i,]<-quantile(var2,probs=c(0.025,0.975))
	intcov[i,]<-quantile(covar,probs=c(0.025,0.975))
	intsig[i,]<-quantile(sigma2.post[2001:10000],probs=c(0.025,0.975))
}