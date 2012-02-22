library(mvtnorm)
library(lme4)
library(nlme)
library(MASS)
library(MCMCpack)
library(magic)
meanb<-c(2,2)
sigma<-diag(2)
sigma[1,2] <-sigma[2,1]<-0.6
sigma[2,2] <-1
fix.beta<-2
Fx<-1-pnorm((2-2-2*0.05-2*0.1)/sqrt(1+0.05^2+0.1*0.8))
POD<-pnorm((2+2*0.05+2*0.1-2)/sqrt(1+0.05^2+0.1*0.8+1))
##set.seed(156)
set.seed(156)
beta<- rmvnorm(15, meanb, sigma)
set.seed(156)
error<-rnorm(300,0,1)
rand.beta<-matrix(NA,nrow=300,ncol=2)
for (i in 1:2)rand.beta[,i]<-rep(beta[,i],each=20)
X1<- seq(0.05,by=0.05,length=20)
X2<-seq(0.1,0.8,0.05)
X1<-rep(X1,15)
X2<-rep(X2,each=20)
dat<-rand.beta[,1]+rand.beta[,2]*X1+fix.beta*X2+error
subject<-rep(seq(1,15,1),each=20)
list1<-data.frame(cbind(subject,X1,X2,dat))
list1$subject<-as.factor(list1$subject)
res<-lme(dat~X1+X2, data=list1, random=~X1|subject,method="ML")
Xmat<-cbind(rep(1,300),list1$X1,list1$X2)
Zmat<-cbind(rep(1,20),list1$X1[1:20])
uf<-2
x1<-0.05
x2<-0.1

beta0_est<-res$coefficients$fixed[1]
beta1_est<-res$coefficients$fixed[2]
beta2_est<-res$coefficients$fixed[3]
sigma_est<-res$sigma
var0_est<- as.numeric(VarCorr(res)[1,1])
var1_est<- as.numeric(VarCorr(res)[2,1])
corr_est<- as.numeric(VarCorr(res)[2,3])
meanb_sim<-c(beta0_est,beta1_est)
sigma_sim<-diag(2)
sigma_sim[1,1] <-var0_est
sigma_sim[1,2] <-sigma_sim[2,1]<-corr_est*sqrt(var0_est*var1_est)
sigma_sim[2,2] <-var1_est
fix.beta_sim<-beta2_est

## function
se_Z<-function(a1,a2,a3,a4,a5,a6,a7){
  matrix1<-matrix(c(2*sqrt(a4),a6*sqrt(a5),a6*sqrt(a5),0),nrow=2)
  matrix2<-matrix(c(0,a6*sqrt(a4),a6*sqrt(a4),2*sqrt(a5)),nrow=2)
  matrix3<-matrix(c(0,sqrt(a4)*sqrt(a5),sqrt(a4)*sqrt(a5),0),nrow=2)
  dersigma1<-Zmat%*%matrix1%*%t(Zmat)
  dersigma2<-Zmat%*%matrix2%*%t(Zmat)
  dersigma3<-Zmat%*%matrix3%*%t(Zmat)
  dersigma4<-2*a7*diag(20)
  Sigmahat<-matrix(c(a4,a6*sqrt(a4*a5),a6*sqrt(a4*a5),a5),nrow=2)
  V<-a7^2*diag(20)+Zmat%*%Sigmahat%*%t(Zmat)
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
  L<-adiag(H11,15*M)	
  fun1<-sqrt(a4+(x1^2)*a5+2*x1*a6*sqrt(a4*a5))
  fun2<-uf-a1-a2*x1-a3*x2
  k<-fun2/fun1
  derkb0<--1/fun1
  derkb1<--x1/fun1
  derkb2<--x2/fun1
  derksig0<--(sqrt(a4)+x1*a6*sqrt(a5))*fun2/(fun1^3)
  derksig1<--(x1^2*sqrt(a5)+x1*a6*sqrt(a4))*fun2/(fun1^3)
  derkrho<--x1*sqrt(a4)*sqrt(a5)*fun2/(fun1^3)
  derk<-c(derkb0,derkb1,derkb2,derksig0,derksig1,derkrho,0)
  Sigmak<-derk%*%solve(L)%*%derk
  return(sqrt(Sigmak))
}
k_fun<-function(a1,a2,a3,a4,a5,a6){
  fun1<-sqrt(a4+(x1^2)*a5+2*x1*a6*sqrt(a4*a5))
  fun2<-uf-a1-a2*x1-a3*x2
  k<-fun2/fun1
  return(k)
}
k_hat<-k_fun(beta0_est,beta1_est,beta2_est,var0_est,var1_est,corr_est)
## k_true<-k_fun(meanb,fix.beta,sigma[1,1],sigma[2,2],sigma[1,2])
k_true <- k_fun(meanb[1],meanb[2],fix.beta,sigma[1,1],sigma[2,2],sigma[1,2])
set.seed(156)
N <- 1000
t <- 0
t_boot_res <- numeric()

## version 1
error_num1 <- 0
error_num2 <- 0
system.time(
            while(t < N){
              beta_sim <- rmvnorm(15 , meanb_sim, sigma_sim)
              error_sim<-rnorm(300,0,sigma_est)
              rand.beta_sim<-matrix(NA,nrow=300,ncol=2)
              for (i in 1:2)rand.beta_sim[,i]<-rep(beta_sim[,i],each=20)
              dat_sim<-rand.beta_sim[,1]+rand.beta_sim[,2]*X1+fix.beta_sim*X2+error_sim
              list1_sim<-data.frame(cbind(subject,X1,X2,dat_sim))
              list1_sim$subject<-as.factor(list1$subject)
              tryres <- try(res_sim<-lme(dat_sim~X1+X2,
                                         data=list1_sim,
                                         random=~X1|subject,method="ML"), silent = TRUE)
              if(!inherits(tryres, "try-error")){
			 
		
		beta0_boot<-res_sim$coefficients$fixed[1]
		beta1_boot<-res_sim$coefficients$fixed[2]
		beta2_boot<-res_sim$coefficients$fixed[3]
		sigma_boot<-res_sim$sigma
		var0_boot<- as.numeric(VarCorr(res_sim)[1,1])
		var1_boot<- as.numeric(VarCorr(res_sim)[2,1])
		corr_boot<- as.numeric(VarCorr(res_sim)[2,3])
                tryres  <- try(se_boot<-se_Z(beta0_boot,beta1_boot,
                                             beta2_boot,var0_boot,
                                             var1_boot,corr_boot,sigma_boot),
                               silent = TRUE)
		if(inherits(tryres, "try-error")){
		error_num2 <- error_num2 + 1
                  next()
		}
		tryres <- try(se_hat<-se_Z(beta0_est,beta1_est,
                                           beta2_est,var0_est,
                                           var1_est,corr_est,sigma_est), silent = TRUE)
		if(inherits(tryres, "try-error")){
                  next()
		}
		t <- t + 1
		k_boot<-k_fun(beta0_boot,beta1_boot,beta2_boot,
                              var0_boot,var1_boot,corr_boot)
		t_boot<-(k_boot-k_hat)/se_boot
		t_boot_res <- c(t_boot_res, t_boot)
              }else{
			   error_num1 <- error_num1 + 1
			  }
            }
            )
## version 2
library(doMC)
registerDoMC()
getDoParWorkers()
set.seed(156)

t_boot_res2 <- numeric()
system.time(
            t_boot_res2<- foreach(icount(N), .combine = c) %dopar% {
              flag <- TRUE
            while(flag){
              beta_sim <- rmvnorm(15 , meanb_sim, sigma_sim)
              error_sim<-rnorm(300,0,sigma_est)
              rand.beta_sim<-matrix(NA,nrow=300,ncol=2)
              for (i in 1:2)rand.beta_sim[,i]<-rep(beta_sim[,i],each=20)
              dat_sim<-rand.beta_sim[,1]+rand.beta_sim[,2]*X1+fix.beta_sim*X2+error_sim
              list1_sim<-data.frame(cbind(subject,X1,X2,dat_sim))
              list1_sim$subject<-as.factor(list1$subject)
              tryres <- try(res_sim<-lme(dat_sim~X1+X2,
                                         data=list1_sim,
                                         random=~X1|subject,method="ML"), silent = TRUE)
              if(!inherits(tryres, "try-error")){
		
		beta0_boot<-res_sim$coefficients$fixed[1]
		beta1_boot<-res_sim$coefficients$fixed[2]
		beta2_boot<-res_sim$coefficients$fixed[3]
		sigma_boot<-res_sim$sigma
		var0_boot<- as.numeric(VarCorr(res_sim)[1,1])
		var1_boot<- as.numeric(VarCorr(res_sim)[2,1])
		corr_boot<- as.numeric(VarCorr(res_sim)[2,3])
                tryres  <- try(se_boot<-se_Z(beta0_boot,beta1_boot,
                                             beta2_boot,var0_boot,
                                             var1_boot,corr_boot,sigma_boot),
                               silent = TRUE)
		if(inherits(tryres, "try-error")){
                  next()
		}
		tryres <- try(se_hat<-se_Z(beta0_est,beta1_est,
                                           beta2_est,var0_est,
                                           var1_est,corr_est,sigma_est), silent = TRUE)
		if(inherits(tryres, "try-error")){
                  next()
		}
		k_boot<-k_fun(beta0_boot,beta1_boot,beta2_boot,
                              var0_boot,var1_boot,corr_boot)
		t_boot<-(k_boot-k_hat)/se_boot
		## t_boot_res2 <<- c(t_boot_res2, t_boot)
                flag <- FALSE
              }
            }
              ## t_boot_res2 <- c(t_boot_res2, t_boot)
              t_boot
          }
            )

identical(sort(t_boot_res), sort(t_boot_res2))

tstat<-(k_hat-k_true)/se_hat
int_boot<-k_hat-se_hat*quantile(t_boot_res,c(.975,.025))
int_boot
int_boot<-k_hat-se_hat*quantile(t_boot_res2,c(.975,.025))
int_boot
## sum(quantile(t_boot_res,0.025)>tstat)
## sum(quantile(t_boot_res,0.975)<tstat)


## sum(quantile(t_boot_res2,0.025)>tstat)
## sum(quantile(t_boot_res2,0.975)<tstat)
