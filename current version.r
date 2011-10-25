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
#set.seed(156)
set.seed(1121)
beta<- rmvnorm(15, mean, sigma)
#set.seed(156)
set.seed(1121)
error<-rnorm(300,0,1)
rand.beta<-matrix(NA,nrow=300,ncol=2)
for (i in 1:2)rand.beta[,i]<-rep(beta[,i],each=20)
x1<- seq(0.05,by=0.05,length=20)
x2<-seq(0.1,0.8,0.05)
x1<-rep(x1,15)
x2<-rep(x2,each=20)
dat<-rand.beta[,1]+rand.beta[,2]*x1+fix.beta*x2+error
subject<-rep(seq(1,15,1),each=20)
list1<-data.frame(cbind(subject,x1,x2,dat))
list1$subject<-as.factor(list1$subject)
res<-lme(dat~x1+x2, data=list1, random=~x1|subject,method="ML")
write.csv(list1,"C:/Jia/research/data.csv")
#MCMC sample 
library(bayesm)
nreg<-15
#V<-3*matrix(c(1,0,0,0,1,0,0,0,0),byrow=T,ncol=3)
regdata<-list(NA)
for (i in 1:nreg){
    regdata[[i]]<-list(list1$dat[list1$subject==i],cbind(rep(1,20),list1$x1[list1$subject==i],list1$x2[list1$subject==i]))
	           }
Z=as.matrix(rep(1,nreg))
data1<-list(regdata=regdata,Z=Z)
out<-rhierLinearModel(Data=data1,Mcmc=list(R=5000,keep=1000))
library(lmm)
xcol<-1:3
zcol<-1:2
fmcmc.result<-fastmcmc.lmm(y=dat,subj=subject,pred=cbind(rep(1,300),list1$x1,list1$x2),xcol,zcol,prior=list(a=3*var(dat),b=3,c=2,Dinv=2*diag(2)),seed=2345,iter=10000)
gibbs.result<-mgibbs.lmm(y=dat,subj=subject,pred=cbind(rep(1,300),list1$x1,list1$x2),xcol,zcol,prior=list(a=3*var(dat),b=3,c=2,Dinv=2*diag(2)),seed=2345,iter=10000)
plot(gibbs.result$sigma2.series)
plot(fmcmc.result$sigma2.series)

par(mfrow=c(4,1))
acf(log(gibbs.result$psi.series[1,1,]),lag.max=10, ylim=0:1)
acf(log(fmcmc.result$psi.series[1,1,]),lag.max=10, ylim=0:1)
acf(gibbs.result$psi.series[1,2,],lag.max=10)
acf(fmcmc.result$psi.series[1,2,],lag.max=10)


Z<-cbind(rep(1,20),list1$x1[1:20])	           
n<-300
m<-15
p<-2
#BETA<-NULL
#list1$newdat<-list1$dat-list1$x2*res$coefficients$fixed[3]
#for (j in 1:m){ BETA<-rbind(BETA,lm(newdat[subject==j]~x1[subject==j],data=list1)$coef)}
#mu0<-apply(BETA,2,mean)
#mu0[3]<-2
#beta<-mu0
#beta<-mu0<-c(2,2,2)
beta<-mu0<-c(beta1_winbugs[1],beta2_winbugs[2],beta3_winbugs[3])
#S0_b<-cov(BETA)
#S0_b<-matrix(c(1,0.6,0.6,1),nrow=2)
S0_b<-matrix(c(var1_winbugs[1],covar_winbugs[1],covar_winbugs[1],var2_winbugs[1],nrow=2)
eta0<-p+2
L0<-100*diag(3)
iSigma0<-solve(S0_b)
#nu0<-2
nu0<-0.002
#sigma0<-var(dat)
#sigma0<-1
sigma1<-sigma_winbugs[1]
S<-10000
bj<-matrix(rep(0,30),nrow=15,ncol=2)
beta.post<-list(NA);sigma2.post<-NULL;Sigma.post<-list(NA);b.post<-list(NA)
#set.seed(2345)
set.seed(5678)
#Rprof()
for(s in 1:S){
	#updata Sigma
	 iSigma<-rwishart(1,eta0+m,solve(S0_b+t(bj)%*%bj))
	 Sigma<-solve(iSigma)
	#update sigma
	rand.BETA<-matrix(NA,nrow=300,ncol=2)
    for (i in 1:2)rand.BETA[,i]<-rep((bj[,i]+beta[i]),each=20)
	SSR<-sum((dat-rand.BETA[,1]-rand.BETA[,2]*list1$x1-beta[3]*list1$x2)^2)
	sigma2<-1/rgamma(1,(n+nu0)/2,(nu0*sigma0+SSR)/2)
    #updata beta
    #V<-sigma2*diag(20)+Z%*%Sigma%*%t(Z)
	XV<-matrix(0,nrow=3,ncol=3);Xlist<-matrix(NA,nrow=20,ncol=3);XVy<-matrix(0,nrow=3,ncol=1)
	for (i in 1:m){ 
	       Xlist<-cbind(rep(1,20),list1$x1[list1$subject==i],list1$x2[list1$subject==i])
		   #XV<-t(Xlist)%*%solve(V)%*%Xlist+XV
		   XV<-t(Xlist)%*%Xlist/sigma2+XV
		   #XVy<-t(Xlist)%*%solve(V)%*%dat[list1$subject==i]+XVy
		   XVy<-t(Xlist)%*%(dat[list1$subject==i]-Z%*%bj[i,])/sigma2+XVy
		       }
	Sigma.m<-solve(solve(L0)+XV)
	mu.m<-Sigma.m%*%(solve(L0)%*%mu0+XVy)
	beta<-rmvnorm(1,mu.m,Sigma.m)
	#updata b
	#Sigma.b<-solve(t(Z)%*%solve(sigma2*diag(20))%*%Z+solve(Sigma))
	Sigma.b<-solve(t(Z)%*%Z/sigma2+solve(Sigma))
	mub<-matrix(NA,nrow=m,ncol=p)
	for (i in 1:m){
	     #mub[i,]<-Sigma.b%*%t(Z)%*%solve(sigma2*diag(20))%*%(dat[list1$subject==i]-cbind(rep(1,20),list1$x1[list1$subject==i],list1$x2[list1$subject==i])%*%t(BETA))
		 mub[i,]<-Sigma.b%*%t(Z)%*%(dat[list1$subject==i]-cbind(rep(1,20),list1$x1[list1$subject==i],list1$x2[list1$subject==i])%*%t(beta))/sigma2
		 bj[i,]<-rmvnorm(1,mub[i,],Sigma.b)
		          }
	beta.post[[s]]<-beta
    Sigma.post[[s]]<-Sigma 
    sigma2.post[s]<-sigma2
    b.post[[s]]<-bj	
	}
#Rprof(NULL)
quantile(sigma2.post[2001:10000],probs=c(0.025,0.975))
mean(sigma2.post[2001:10000])
hist(sigma2.post[2000:10000])
beta1<-NULL;beta2<-NULL;beta3<-NULL;var1<-NULL;var2<-NULL;covar<-NULL;tempsigma<-matrix(NA,nrow=2,ncol=2);temp<-NULL
for (i in 2001:10000){
    temp<-beta.post[[i]]
	beta1<-c(temp[1],beta1)
	beta2<-c(temp[2],beta2)
	beta3<-c(temp[3],beta3)
	tempsigma<-Sigma.post[[i]]
	var1<-c(tempsigma[1,1],var1)
	var2<-c(tempsigma[2,2],var2)
	covar<-c(tempsigma[1,2],covar)
	    }
beta1<-as.mcmc(beta1);beta2<-as.mcmc(beta2);beta3<-as.mcmc(beta3);var1<-as.mcmc(var1);var2<-as.mcmc(var2);covar<-as.mcmc(covar)
hist(beta3);hist(beta1);hist(beta2);hist(var1);hist(var2);hist(covar)
mean(beta3);mean(beta1);mean(beta2);mean(var1);mean(var2);mean(covar)
quantile(beta1,probs=c(0.025,0.975));quantile(beta2,probs=c(0.025,0.975));quantile(beta3,probs=c(0.025,0.975))
quantile(var1,probs=c(0.025,0.975));quantile(var2,probs=c(0.025,0.975));quantile(covar,probs=c(0.025,0.975))
cummean<-function(x){
   x.mx<-matrix(x,ncol=100,byrow=T)
   csum<-cumsum(rowSums(x.mx))
   Ntotal<-seq(from=100,to=length(x),by=100)
   res<-csum/Ntotal
   return(res)
       }   
cummean(beta1);cummean(beta2);cummean(beta3);cummean(var1);cummean(var2);cummean(covar);cummean(sigma2.post[2001:10000])

beta_winbugs<-read.table("C:/Jia/research/beta_winbugs.txt")
beta1_winbugs<-beta_winbugs[1:8000,2]
beta2_winbugs<-beta_winbugs[8001:16000,2]
beta3_winbugs<-beta_winbugs[16001:24000,2]
cummean(beta1_winbugs);cummean(beta2_winbugs);cummean(beta3_winbugs)

cov_winbugs<-read.table("C:/Jia/research/cov_winbugs.txt")
var1_winbugs<-cov_winbugs[1:8000,2]
covar_winbugs<-cov_winbugs[8001:16000,2]
var2_winbugs<-cov_winbugs[24001:32000,2]
cummean(var1_winbugs);cummean(var2_winbugs);cummean(covar_winbugs)

sigma_winbugs<-read.table("C:/Jia/research/sigma_winbugs.txt")
sigma_winbugs<-sigma_winbugs[,2]
cummean(sigma_winbugs)

library("ggplot2")
cm1 <- cummean(beta1)
cm1b <- cummean(beta1_winbugs)
cm2 <- cummean(beta2)
cm2b <- cummean(beta2_winbugs)
cm3 <- cummean(beta3)
cm3b <- cummean(beta3_winbugs)
cm4 <- cummean(var1)
cm4b <- cummean(var1_winbugs)
cm5 <- cummean(covar)
cm5b <- cummean(covar_winbugs)
cm6 <- cummean(var2)
cm6b <- cummean(var2_winbugs)
cm7 <- cummean(sigma2.post[2001:10000])
cm7b <- cummean(sigma_winbugs)

dataframe <- data.frame(x = rep(1:80, 7), mean = c(cm1, cm1b,cm2, cm2b,cm3, cm3b,cm4, cm4b,cm5, cm5b,cm6, cm6b,cm7, cm7b), 
                 var = rep(c("beta0", "beta1", "beta2", "var1", "covar", "var2", "sigma"), each = 80*2),
                 group = rep(c(rep("R", 80), rep("winbugs", 80)), 7))
plot_compare <- qplot(data = dataframe, y = mean, x = x, group = group, geom = "line", color = group)
pdf("C:/Jia/research/Rwin2.pdf")
plot_compare + facet_wrap(~var,scale="free")
dev.off()

plot(beta3);plot(beta1);plot(beta2);plot(var1);plot(var2);plot(covar)
library(coda)
install.packages("mcgibbsit") 
library("mcgibbsit") 
autocorr.plot(beta1);autocorr.plot(beta2);autocorr.plot(beta3);autocorr.plot(var1);autocorr.plot(var2);autocorr.plot(covar)
coda::traceplot(beta1);coda::traceplot(beta2);coda::traceplot(beta3);coda::traceplot(var1);coda::traceplot(var2);coda::traceplot(covar)
raftery.diag(beta1);raftery.diag(beta2);raftery.diag(beta3);raftery.diag(var1);raftery.diag(var2);raftery.diag(covar)
HPDinterval(beta1);HPDinterval(beta2);HPDinterval(beta3);HPDinterval(var1);HPDinterval(var2);HPDinterval(covar)
cumuplot(beta1);cumuplot(beta2);cumuplot(beta3);cumuplot(var1);cumuplot(var2);cumuplot(covar)
pred_dat<-matrix(NA,nrow=300,ncol=8001);diffdat<-matrix(NA,nrow=300,ncol=8001);cumdiff<-NULL
for (i in 2000:10000){
	bjpost<-b.post[[i]]
	pred_dat[,i-1999]<-bjpost[,1]+beta1[i-1999]+x1*(bjpost[,2]+beta2[i-1999])+beta3[i-1999]*x2+sqrt(sigma2.post[i])
	diffdat[,i-1999]<-pred_dat[,i-1999]-dat
	cumdiff[i-1999]<-sum(diffdat[,i-1999]^2)
	}
pred_dat<-rand.beta[,1]+rand.beta[,2]*x1+fix.beta*x2+error


install.packages("arm") 
install.packages("R2WinBUGS") 
library("arm")
library("BRugs")
library("R2WinBUGS")
y<-matrix(0,nrow=15,ncol=20);x1<-matrix(0,nrow=15,ncol=20);x2<-matrix(0,nrow=15,ncol=20)
for(i in 1:15){
    y[i,1:20]<-list1$dat[(20*i-19):(20*i)]
	x1[i,1:20]<-list1$x1[(20*i-19):(20*i)]
	x2[i,1:20]<-list1$x2[(20*i-19):(20*i)]
               }
mu.b=c(2,2,2)
prec = structure(.Data = c(0.01, 0, 0,0,0.01,0,0,0,0.01), .Dim = c(3, 3))
Omega = structure(.Data = c(2.5,0.55,0.55,1),.Dim = c(2, 2))
mu.beta<-c(0,0)
data<-list("y","x1","x2","mu.b","prec","Omega","mu.beta")
inits<-function() {list (mu.b=beta,prec = structure(.Data = c(0.01, 0, 0,0,0.01,0,0,0,0.01), .Dim = c(3, 3)))}
inits<-function() {list (mu.b=c(2,2,2),prec = structure(.Data = c(0.01, 0, 0,0,0.01,0,0,0,0.01), .Dim = c(3, 3)), Omega = structure(.Data = c(2.5,0.55,0.55,1),.Dim = c(2, 2)))}
parameters <- c("b", "sigma.y", "Sigma.b")
dput(list1$dat,'')

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
stars.l <- seq(from=star[7],to=0.001,length.out=50)
opt_in<-nlminb(start=star[-7],prof_loglike3,theta=stars.l[2],lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
yhat.l<-c(1,exp(-opt_in$objective+max_tran_linear$objective))
pars<-opt_in$par
for (st in stars.l[c(-1,-2)]){
  res_opt<-nlminb(start=pars,prof_loglike3,theta=st,lower=c(rep(-Inf,2),0.001,-Inf,0.001,0.001))
    obj <- res_opt$objective
  yhat.l <-c(yhat.l,exp(-obj+max_tran_linear$objective))
  pars <- res_opt$par
}

stars.r <- seq(from=star[7],to=0.999,length.out=51)
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
plot(c(stars.l,stars.r[-1]),c(yhat.l,yhat.r[-1]),main="100 points in total",xlab="F(x)",ylab="norm likelihood")
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


#########################Hessian matrix####################################
yy<-matrix(0,nrow=(20+2),ncol=15)
for(i in 1:15){
yy[1:20,i]<-dat[((i-1)*20+1): (i*20)]
}
Xmat<-cbind(rep(1,300),list1$x1,list1$x2)
Xmat1<-matrix(NA, nrow=20,ncol=45)
for (i in 1:15){
     Xmat1[,(3*i-2):(i*3)]<-Xmat[(20*i-19):(20*i),]
	}
xx<-list(NA)
for( i in 1:15){
    xx[[i]]<-rbind(Xmat1[,(3*i-2):(i*3)],matrix(rep(0,6),nrow=2,ncol=3))
	           }
Zmat<-cbind(rep(1,20),list1$x1[1:20])

loglike<-function(ps){
    beta<-c(ps[1],ps[2],ps[3])
	delta<-matrix(c(ps[4:5],0,ps[6]),nrow=2,ncol=2,byrow=T)
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
	error<-ps[7]
    result<- -15*10*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-result)
}
maximum<-nlm(loglike,c(2,2,2,1.7,0,1,1))
#maximum1<-nlm(loglike,c(1,1,1,2,0,5,5))
	
est<-maximum$estimate
delta<-matrix(c(est[4:5],0,est[6]),nrow=2,byrow=T)
Dinv<-t(delta)%*%delta
Sigmahat<-solve(Dinv)*est[7]^2
sdev<-sqrt(diag(Sigmahat))
tmp<-diag(1/sdev)
corrm<- tmp%*%Sigmahat%*%tmp

library(magic)
matrix1<-matrix(c(2*sdev[1],corrm[1,2]*sdev[2],corrm[1,2]*sdev[2],0),nrow=2)
matrix2<-matrix(c(0,corrm[1,2]*sdev[1],corrm[1,2]*sdev[1],2*sdev[2]),nrow=2)
matrix3<-matrix(c(0,sdev[1]*sdev[2],sdev[1]*sdev[2],0),nrow=2)
dersigma1<-Zmat%*%matrix1%*%t(Zmat)
dersigma2<-Zmat%*%matrix2%*%t(Zmat)
dersigma3<-Zmat%*%matrix3%*%t(Zmat)
dersigma4<-2*est[7]*diag(20)
V<-est[7]^2*diag(20)+Zmat%*%Sigmahat%*%t(Zmat)
M11<-sum(diag(solve(V)%*%dersigma1%*%solve(V)%*%dersigma1))/2
M12<-M21<-sum(diag(solve(V)%*%dersigma1%*%solve(V)%*%dersigma2))/2
M13<-M31<-sum(diag(solve(V)%*%dersigma1%*%solve(V)%*%dersigma3))/2
M14<-M41<-sum(diag(solve(V)%*%dersigma1%*%solve(V)%*%dersigma4))/2
M22<-sum(diag(solve(V)%*%dersigma2%*%solve(V)%*%dersigma2))/2
M23<-M32<-sum(diag(solve(V)%*%dersigma2%*%solve(V)%*%dersigma3))/2
M24<-M42<-sum(diag(solve(V)%*%dersigma2%*%solve(V)%*%dersigma4))/2
M33<-sum(diag(solve(V)%*%dersigma3%*%solve(V)%*%dersigma3))/2
M34<-M43<-sum(diag(solve(V)%*%dersigma3%*%solve(V)%*%dersigma4))/2
M44<-sum(diag(solve(V)%*%dersigma4%*%solve(V)%*%dersigma4))/2
M<-matrix(c(M11,M12,M13,M14,M21,M22,M23,M24,M31,M32,M33,M34,M41,M42,M43,M44),nrow=4,byrow=T)
H11<-matrix(rep(0,9),nrow=3)
for (i in 1:15){
     tempH<-t(Xmat[(20*i-19):(20*i),])%*%solve(V)%*%Xmat[(20*i-19):(20*i),]
	 H11<-tempH+H11
	  }
I<-adiag(H11,15*M)	
	  
#use delta method to get the var-cov matrix for the function 
uf<-2
x1<-0.05
x2<-0.1
fun1<-sqrt(sdev[1]^2+(x1*sdev[2])^2+2*x1*corrm[1,2]*sdev[1]*sdev[2])
fun2<-uf-est[1]-est[2]*x1-est[3]*x2
derkb0<--1/fun1
derkb1<--x1/fun1
derkb2<--x2/fun1
derksig0<--(sdev[1]+x1*corrm[1,2]*sdev[2])*fun2/(fun1^3)
derksig1<--(x1^2*sdev[2]+x1*corrm[1,2]*sdev[1])*fun2/(fun1^3)
derkrho<--x1*sdev[1]*sdev[2]*fun2/(fun1^3)
k<-fun2/fun1
derk<-c(derkb0,derkb1,derkb2,derksig0,derksig1,derkrho,0)
Sigmak<-derk%*%solve(I)%*%derk
int_k<-c(k-1.96*sqrt(Sigmak),k+1.96*sqrt(Sigmak))
F<-1-pnorm(k)
int_F<-1-pnorm(int_k)

derF<--dnorm(k)*derk
SigmaF<-derF%*%solve(I)%*%derF
w<-exp(1.96*sqrt(SigmaF)/(F*(1-F)))
int_F2<-c(F/(F+(1-F)*w),F/(F+(1-F)/w))

