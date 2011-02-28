rand.val1 <- rnorm(mean=3,sd=1,n=10)
rand.val2<-rnorm(mean=0,sd=2,n=100)
error.each <- rnorm(mean=0,sd=0.3,n=10*10*10)
rand.exp1 <- rep(rand.val1,each=100)
rand.exp2<-rep(rand.val2,each=10)
data.exp <- rep(1,1000)+rand.exp1+rand.exp2+error.each
ID1<-seq(1,10,1)
ID2<-seq(1,100,1)
ID.exp1 <- factor(rep(ID1,each=100))
ID.exp2<-factor(rep(ID2,each=10))

fit.1 <- lm(data.exp~0+ID.exp1+ID.exp2)
anova(fit.1)
library(lme4)
fit.2 <- lmer(data.exp~ (1|ID.exp1)+(1|ID.exp2))


write.table(data.exp,file="rand.frame.dat",eol=",\n",row.names=F)

library(mvtnorm)
library(lme4)
x1<-seq(1,10,1)
mean<-rep(5,2)
sigma<-diag(2)
sigma[1,1] <-10
sigma[1,2] <-6
sigma[2,2] <-10
sigma[2,1] <-6
beta<- rmvnorm(10, mean, sigma)
error<-rnorm(100,0,10)
rand.beta<-matrix(NA,nrow=100,ncol=2)
rand.beta[,1]<-rep(beta[,1],each=10)
rand.beta[,2]<-rep(beta[,2],each=10)
time<-rep(x1,10)
y<-rand.beta[,1]+rand.beta[,2]*time+error
subject<-rep(seq(1,10,1),each=10)
list<-data.frame(cbind(subject,time,y))
list$subject<-as.factor(list$subject)
res<-lmer(y~1+time+(1+time|subject),list)
set.seed(101)
sample<-mcmcsamp(res,n=100)

res2<-lmer(y~1+time+(1|subject)+(0+time|subject),list)
library(MCMCpack)
res_bay<-MCMCregress(y~time,list)
summary(res_bay)

library(mvtnorm)
library(lme4)
library(nlme)
library(MASS)
#library(languageR)
#library(lmm)

mean<-c(0,0,0)

sigma<-diag(3)
sigma[1,1] <-1
sigma[1,2] <-sigma[2,1]<-0
sigma[2,2] <-1
sigma[3,3]<-1
sigma[1,3]<-sigma[3,1]<-0.9
sigma[2,3]<-sigma[3,2]<-0.1

set.seed(100)
beta<- rmvnorm(15, mean, sigma)
set.seed(100)
error<-rnorm(240,0,1)
rand.beta<-matrix(NA,nrow=240,ncol=3)
for (i in 1:3)rand.beta[,i]<-rep(beta[,i],each=16)
fix.beta<-2
             
x1<-seq(0.1,0.85,0.25)
x2<-seq(0.25,0.85,0.2)
x3<-seq(0.1,0.8,0.05)
time<-rep(rep(x1,4),15)
temp<-rep(rep(x2,each=4),15)
length<-rep(x3,each=16)


dat<-rand.beta[,1]+rand.beta[,2]*time+rand.beta[,3]*temp+fix.beta*length+error
subject<-rep(seq(1,15,1),each=16)
list<-data.frame(cbind(subject,time,temp,length,dat))
list$subject<-as.factor(list$subject)
res<-lme(dat~length+time+temp, data=list, random=~time+temp|subject,method="ML")
qqnorm(res,~resid(.)|subject)
res_con<-lme(dat~length+time+temp, data=list, random=~time+temp|subject,control=list(msVerbose=TRUE,maxIter=100,opt="optim"),method="ML")

intervals(res)

res2<-lmer(dat~length+time+temp+(time+temp|subject),REML=F,data=list)

add2Matrix <- function(mx,skip=16,num=0,row=3){
 if(!is.null(dim(mx))){
  k <- nrow(mx)%/%skip-1
  k2 <- nrow(mx)%%skip
  lst <- lapply(0:k,function(i){
    st <- i*skip+1
    do.call('rbind',c(list(mx[st:(st+skip-1),]),as.list(rep(num,row))))
  })
  nmx <- do.call('rbind',lst)
  if(k2!=0){
  nmx <- rbind(nmx,mx[(tail(k,1)+1)*skip+(1:k2),])
}
  nmx
}else{
 if(length(mx)==0) stop('Please give a valid data')
  k <- length(mx)%/%skip-1
  k2 <- length(mx)%%skip
  lst <- lapply(0:k,function(i){
    st <- i*skip+1
    do.call('c',c(list(mx[st:(st+skip-1)]),as.list(rep(num,row))))
  })
  nmx <- do.call('c',lst)
  if(k2!=0){
  nmx <- c(nmx,mx[(tail(k,1)+1)*skip+(1:k2)])
}
  nmx
}
}

Xmat<-cbind(rep(1,240),list$length,list$time,list$temp)
lst <- lapply(0:14,function(i){
 st <- i*16+1
 do.call('rbind',list(Xmat[st:(st+15),],0,0,0))
})
Xmatcomp<- do.call('rbind',lst)
ymat<-list$dat
ye<-add2Matrix(ymat)
Zmati<-cbind(rep(1,16),list$time[1:16],list$temp[1:16])
Xmat1<-matrix(NA, nrow=16,ncol=60)
for (i in 1:15){
     Xmat1[,(4*i-3):(i*4)]<-Xmat[(16*i-15):(16*i),]
	}
Rqr<-qr(R22)
Rqr_q<-qr.Q(Rqr,complete=TRUE)
Rqr_q1<-qr.Q(Rqr)
Rqr_r<-qr.R(Rqr,complete=TRUE)
Rqr_r1<-qr.R(Rqr)	
w<-t(Rqr_q)%*%C2
w1<-w[1:4]
w2<-w[5:60]
C<-c(C,w1)
R<-rbind(R12,Rqr_r1)
k<-norm(as.matrix(w2),"F")^2+sum(C3)

library(magic)
loglhd_qr<-function(ps)
# beta1,beta2,beta3,beta4,phi1,phi12,phi13,phi2,phi23,phi3
# phi are the variance components, not the variance-covariance matrix
{
     #browser()
     beta_fix<-c(ps[1],ps[2],ps[3],ps[4])
     M<-15
     N<-240
     n<-16
     L11<-ps[5]
	 #if (L11<=0) L11<-0.001
     L21<-ps[6]/L11
	 #if (ps[8]<0) ps[8]<-0.001
     l22<-ps[8]^2-L21^2
     #if(l22<=0) l22<-0.001
     L22<-sqrt(l22)
     L31<-ps[7]/L11
     L32<-(ps[9]-L31*L21)/L22
	 #if (ps[10]<0) ps[10]<-0.001
     l33<-ps[10]^2-L31^2-L32^2
     #if (l33<0) l33<-0.001
     L33<-sqrt(l33)
     L<-matrix(c(L11,0,0,L21,L22,0,L31,L32,L33),nrow=3,ncol=3,byrow=T)
     phi<-L%*%t(L)
	 #if (ps[11]<0) ps[11]<-0.001
	 error<-ps[11]
     R11<-list(NA)
	 Vi<-list(NA)
	 for (i in 1:15){
	     R11[[i]]<-XZqr_r[[i]][1:3,1:3]
         Vi[[i]]<-diag(3)+R11[[i]]%*%phi%*%t(R11[[i]])
		 }
     V<-do.call("adiag",Vi)
	 V<-adiag(V,diag(4))
	 if (det(V)<2.220446e-16) V<-V+diag(0.01,49,49) else
	     if (1/det(V)<2.220446e-16) V<- diag(0.01,49,49)
     res<--log(det(V))/2-log(error^2)*N/2-error^2*((t(C-R%*%beta_fix)%*%solve(V)%*%(C-R%*%beta_fix))+k)/2
     return(-res)
         }
 maximum_qr<-nlm(f=loglhd_qr,p=c(-0.2,2,0.5,0.2,0.25,0.1,0.05,0.7,0.3,0.8,1),gradtol = 1e-16,iterlim = 1000,hessian=T) 
max_qr<-nlminb(c(-0.2,2,0.5,0.2,0.25,0.1,0.05,0.7,0.3,0.8,1),loglhd_qr,lower=c(-1,1,rep(1e-16,9)),upper=c(1e-16,3,rep(1,8),2))
max_qr<-optim(c(-0.2,2,0.5,0.2,0.25,0.1,0.05,0.7,0.3,0.8,1),loglhd_qr,lower=c(-1,1,rep(1e-16,9)),upper=c(1e-16,3,rep(1,8),2),method="L-BFGS-B")


loglhd_psd<-function(ps)
{
     #browser()
     beta_fix<-c(ps[1],ps[2],ps[3],ps[4])
     M<-15
     N<-240
     n<-16
     delta<-matrix(c(ps[5],ps[6],ps[7],0,ps[8],ps[9],0,0,ps[10]),nrow=3,ncol=3,byrow=T)
	 det_delta<-det(delta)
	 #if (det_delta==0) delta<-delta+diag(0.001,3)
	 mat<-list(NA)
	 mat_qr<-list(NA)
	 mat_qr.r<-list(NA)
	 mat_qr.q<-list(NA)
	 R11<-list(NA);R10<-list(NA);R00<-matrix(NA,nrow=240,ncol=4);c1<-list(NA);c0<-NULL;deter<-NULL
	 for (i in 1:15){
	           mat[[i]]<-cbind(rbind(Zmati,delta),rbind(Xmat1[,(4*i-3):(i*4)],matrix(rep(0,12),nrow=3,ncol=4)),c(ymat[(16*i-15):(16*i)],c(0,0,0)))
		       mat_qr[[i]]<-qr(mat[[i]])
               mat_qr.r[[i]]<-qr.R(mat_qr[[i]],complete=TRUE)
               mat_qr.q[[i]]<-qr.Q(mat_qr[[i]],complete=TRUE)
               R11[[i]]<-mat_qr.r[[i]][1:3,1:3]
               deter[i]<-det(R11[[i]])
               #if (deter[i]==0) R11[[i]]<-R11[[i]]+diag(0.001,3)			   
               R10[[i]]<-mat_qr.r[[i]][1:3,4:7]
               R00[(16*i-15):(i*16),]<-mat_qr.r[[i]][4:19,4:7]
               c1[[i]]<-mat_qr.r[[i]][1:3,8]
               c0[(16*i-15):(i*16)]<-mat_qr.r[[i]][4:19,8]		  
					 }
	 mat1<-cbind(R00,c0)
	 mat1_qr<-qr(mat1)
	 mat1_qr.r<-qr.R(mat1_qr,complete=TRUE)
	 mat1_qr.q<-qr.Q(mat1_qr,complete=TRUE)
	 if (ps[11]<=0) {res<--1e10} else{ 
	     error<-ps[11]
	     res<--log(error^2)*N/2-norm(mat1_qr.r%*%c(-beta_fix,1),"F")^2/(2*error^2)+sum(log(abs(det_delta/deter)))
		 }
     return(-res)
         }
 D<-matrix(c(0.25,0.15,0.05,0.15,0.7,0.3,0.05,0.3,0.9),nrow=3,byrow=T)
 d<-solve(D)
 delta<-chol(d)
initial1<-0.375*norm(as.matrix(Zmati[,1]),"F")
initial2<-0.375*norm(as.matrix(Zmati[,2]),"F")
initial3<-0.375*norm(as.matrix(Zmati[,3]),"F")

maximum_psd<-nlminb(c(-0.2,2,0.5,0.2,1.5,0,0,0.8,0,0.9,0.9),loglhd_psd,lower=c(rep(-Inf,10),0.1))
maximum_psd<-nlm(loglhd_psd,c(-0.2,2,0.5,0.2,1.5,0,0,0.8,0,0.9,0.9))
maximum_psd<-nlm(loglhd_psd,c(-0.23,2.16,0.55,0.21,2.15,-0.5,0.04,1.3,-0.43,1,0.1))


loglhd_psd_revised<-function(ps)
{
     #browser()
     beta_fix<-c(ps[1],ps[2],ps[3],ps[4])
     M<-15
     N<-240
     n<-16
     delta<-matrix(c(ps[5],ps[6],ps[7],0,ps[8],ps[9],0,0,ps[10]),nrow=3,ncol=3,byrow=T)
	 det_delta<-det(delta)
	 #if (det_delta==0) delta<-delta+diag(0.001,3)
	 yy<-list(NA)
	 xx<-list(NA)
	 zz<-rbind(Zmati,delta)
	 mat_qr<-qr(zz)
     mat_qr.r<-qr.R(mat_qr,complete=TRUE)
     mat_qr.q<-qr.Q(mat_qr,complete=TRUE)
     R11<-mat_qr.r[1:3,]
     det_R<-det(R11)	 
	 R00<-matrix(NA,nrow=240,ncol=4);c0<-NULL;deter<-NULL
	 for (i in 1:15){
	           yy[[i]]<-c(ymat[(16*i-15):(16*i)],c(0,0,0))
			   xx[[i]]<-rbind(Xmat1[,(4*i-3):(i*4)],matrix(rep(0,12),nrow=3,ncol=4))			            	   
               R00[(16*i-15):(i*16),]<-(t(mat_qr.q)%*%xx[[i]])[4:19,]               
               c0[(16*i-15):(i*16)]<-(t(mat_qr.q)%*%yy[[i]])[4:19]		  
					 }
	 mat1<-cbind(R00,c0)
	 mat1_qr<-qr(mat1)
	 mat1_qr.r<-qr.R(mat1_qr,complete=TRUE)
	 mat1_qr.q<-qr.Q(mat1_qr,complete=TRUE)
	 #if (ps[11]<=0) {res<--1e10} else{ 
	     error<-ps[11]
	     res<--log(error^2)*N/2-norm(mat1_qr.r%*%c(-beta_fix,1),"F")^2/(2*error^2)+sum(log(abs(det_delta/det_R)))
		# }
     return(-res)
         }
maximum_psd_revised<-nlm(loglhd_psd_revised,c(1,5,0.5,0.2,1.5,0,0,0.8,0,0.9,0.9))		 
maximum_psd_revised<-nlm(loglhd_psd_revised,c(-0.23,2.16,0.55,0.21,2.15,-0.5,0.04,1.3,-0.43,1,0.9))
est<-maximum_psd_revised$estimate
delta<-matrix(c(est[5:7],0,est[8:9],0,0,est[10]),nrow=3,byrow=T)
Dinv<-t(delta)%*%delta
Sigmahat<-solve(Dinv)*est[11]^2
sdev<-sqrt(diag(Sigmahat))
tmp<-diag(1/sdev)
corrm<- tmp%*%Sigmahat%*%tmp

library(magic)
profile_theta<-function(ps)
{
     #browser()
     delta<-matrix(c(ps[1],ps[2],ps[3],0,ps[4],ps[5],0,0,ps[6]),nrow=3,ncol=3,byrow=T)
	 det_delta<-det(delta)
	 Zlist<-list(NA)
	 for (i in 1:15){
	    Zlist[[i]]<-rbind(Zmati,delta)
		             }
	 Xe<-cbind(do.call("adiag",Zlist),Xmatcomp)
	 mat<-list(NA)
	 mat_qr<-list(NA)
	 mat_qr.r<-list(NA)
	 mat_qr.q<-list(NA)
	 R11<-list(NA);R10<-list(NA);R00<-matrix(NA,nrow=240,ncol=4);c1<-list(NA);c0<-NULL;deter<-NULL
	 for (i in 1:15){
	           mat[[i]]<-cbind(rbind(Zmati,delta),rbind(Xmat1[,(4*i-3):(i*4)],matrix(rep(0,12),nrow=3,ncol=4)),c(ymat[(16*i-15):(16*i)],c(0,0,0)))
		       mat_qr[[i]]<-qr(mat[[i]])
               mat_qr.r[[i]]<-qr.R(mat_qr[[i]],complete=TRUE)
               mat_qr.q[[i]]<-qr.Q(mat_qr[[i]],complete=TRUE)
               R11[[i]]<-mat_qr.r[[i]][1:3,1:3]
               deter[i]<-det(R11[[i]])
               #if (deter[i]==0) R11[[i]]<-R11[[i]]+diag(0.001,3)			   
               R10[[i]]<-mat_qr.r[[i]][1:3,4:7]
               R00[(16*i-15):(i*16),]<-mat_qr.r[[i]][4:19,4:7]
               c1[[i]]<-mat_qr.r[[i]][1:3,8]
               c0[(16*i-15):(i*16)]<-mat_qr.r[[i]][4:19,8]		  
					 }
	 mat1<-cbind(R00,c0)
	 mat1_qr<-qr(mat1)
	 mat1_qr.r<-qr.R(mat1_qr,complete=TRUE)
	 mat1_qr.q<-qr.Q(mat1_qr,complete=TRUE)
	 c_neg<-mat1_qr.r[-c(1:4),-c(1:4)]
	 res<--240*log(norm(as.matrix(c_neg),"F"))+sum(log(abs(det_delta/deter)))
     return(-res)
         }
max_theta<-nlm(profile_theta,c(2.14,-0.48,0.04,1.3,-0.43,1.05))	
max_theta<-nlm(profile_theta,c(1.5,0,0,0.8,0,0.9)) 

loglhd_psd1<-function(ps)
# beta1,beta2,beta3,beta4,phi1,phi12,phi13,phi2,phi23,phi3
# phi are the variance components, not the variance-covariance matrix
{
     #browser()
     beta_fix<-c(ps[1],ps[2],ps[3],ps[4])
     M<-15
     N<-240
     n<-16
     L11<-ps[5]
	 L21<-ps[6]/L11
	 l22<-ps[8]^2-L21^2
     #if(l22<=0) l22<-1e10
	 if (l22<=0) l22<-NA
     L22<-sqrt(l22)
	 #if (L22==NaNs) L22<-NA
     L31<-ps[7]/L11
     L32<-(ps[9]-L31*L21)/L22
	 l33<-ps[10]^2-L31^2-L32^2
     #if (l33<=0) l33<-1e10
	 if (is.na(l33)==F&l33<=0) l33<-NA
     L33<-sqrt(l33)
	 #if (L33==NaNs) L33<-NA
     a1<-1/L11
	 a2<--L21/(L22*L11)
	 a3<-(-L31/L33+(L21*L32)/(L33*L22))/L11
	 a4<-a7<-a8<-0
	 a5<-1/L22
	 a9<-1/L33
	 a6<--L32/(L33*L22)
	 L<-matrix(c(L11,0,0,L21,L22,0,L31,L32,L33),nrow=3,ncol=3,byrow=T)
     phi<-L%*%t(L)
	 delta<-matrix(c(a1,a2,a3,a4,a5,a6,a7,a8,a9),nrow=3,ncol=3,byrow=T)
	 #if (det(L)<2.220446e-16) L<-L+diag(0.01,3,3) 
	 #delta<-solve(t(L))
	 if (any(is.na(L))==T) {res<--1e15} else
	     {
	       mat<-list(NA)
	       mat_qr<-list(NA)
	       mat_qr.r<-list(NA)
	       mat_qr.q<-list(NA)
	       R11<-list(NA);R10<-list(NA);R00<-matrix(NA,nrow=240,ncol=4);c1<-list(NA);c0<-NULL;deter<-NULL
	       for (i in 1:15){
	           mat[[i]]<-cbind(rbind(Zmati,delta),rbind(Xmat1[,(4*i-3):(i*4)],matrix(rep(0,12),nrow=3,ncol=4)),c(ymat[(16*i-15):(16*i)],c(0,0,0)))
		       mat_qr[[i]]<-qr(mat[[i]])
               mat_qr.r[[i]]<-qr.R(mat_qr[[i]],complete=TRUE)
               mat_qr.q[[i]]<-qr.Q(mat_qr[[i]],complete=TRUE)
               R11[[i]]<-mat_qr.r[[i]][1:3,1:3]
               deter[i]<-det(R11[[i]])		  
               R10[[i]]<-mat_qr.r[[i]][1:3,4:7]
               R00[(16*i-15):(i*16),]<-mat_qr.r[[i]][4:19,4:7]
               c1[[i]]<-mat_qr.r[[i]][1:3,8]
               c0[(16*i-15):(i*16)]<-mat_qr.r[[i]][4:19,8]		  
					    }
	      mat1<-cbind(R00,c0)
	      mat1_qr<-qr(mat1)
	      mat1_qr.r<-qr.R(mat1_qr,complete=TRUE)
	      mat1_qr.q<-qr.Q(mat1_qr,complete=TRUE)
	      #if (ps[11]==0) ps[11]<-1e-15
	      error<-abs(ps[11])
	      #error<- ifelse(error<1e-7,1e-7,error)	 
	      deter<-sapply(deter,function(x){
	      ifelse(x==0,1e-15,x)
	            })
	      detp<-det(phi)
	      detp<-ifelse(detp==0,1e-15,detp)
          res<--log(abs(detp))*M/2-log(error^2)*N/2-norm(mat1_qr.r%*%c(beta_fix,-1),"F")^2/(2*error^2)-sum(log(abs(deter)))
		  }
     return(-res)
         }
maximum_psd1<-nlm(f=loglhd_psd1,p=c(-0.2,2,0.5,0.2,0.3,0.1,0.1,0.7,0.3,0.9,1),hessian=T)  

pred<-cbind(rep(1,240),list$length,list$time,list$temp,kronecker(diag(1,15),cbind(rep(1,16),list$time[1:16],list$temp[1:16])))
xcol<-1:4
zcol<-5:49
prior<-list(a=3*100,b=9,c=9,Dinv=100)
res3<-mgibbs.lmm(list$dat,list$subject,pred,xcol,zcol,prior=prior,seed=1234,iter=5000)




#test by large sample size
mean_t<-c(0,0)
sigma_t<-diag(2)
sigma_t[1,1] <-1
sigma_t[1,2] <-sigma[2,1]<-0
sigma_t[2,2] <-1
x1_t<-seq(1,100,1)
x3_t<-seq(1,100,1)
time_t<-rep(x1_t,100)
length_t<-rep(x3_t,each=100)
error_t<-rnorm(10000,0,1)
beta_t<- rmvnorm(100, mean, sigma)
rand.beta_t<-matrix(NA,nrow=10000,ncol=2)
for (i in 1:2)rand.beta_t[,i]<-rep(beta_t[,i],each=100)

dat_t<-rand.beta_t[,1]+rand.beta_t[,2]*time_t+fix.beta*length_t+error_t
subject_t<-rep(seq(1,100,1),each=100)
list_t<-data.frame(cbind(subject_t,time_t,length_t,dat_t))
list_t$subject_t<-as.factor(list_t$subject_t)
res_t<-lme(dat_t~length_t+time_t, data=list_t, random=~time_t|subject_t,method="ML")


loglhd_t<-function(ps)
# beta1,beta2,beta3,phi1,phi12,phi2
{
     beta_fix<-c(ps[1],ps[2],ps[3])
     #phi<-matrix(c(ps[4]^2,ps[5],ps[5],ps[6]^2),2,2,byrow=T)
     L11<-ps[4]
     if (L11==0) L11<-L11+0.001
     L21<-ps[5]/L11
     l22<-ps[6]^2-L21^2
     if (l22<=0) l22<-0.001
     L22<-sqrt(l22)
     L<-matrix(c(L11,0,L21,L22),nrow=2,ncol=2,byrow=T)
     phi<-L%*%t(L)
     if (det(phi)==0) phi<=phi+diag(0.001,2,2)
	 error<-ps[7]
     X<-cbind(rep(1,10000),list_t$length_t,list_t$time_t)
     y<-list_t$dat_t
     R<-y-X%*%beta_fix
     Rm<-matrix(R,nrow=100,ncol=100,byrow=F)
     Zi<-cbind(rep(1,100),list_t$time_t[1:100])
     Vi<-diag(100)+Zi%*%phi%*%t(Zi)
     res1<-0
     for (i in 1:ncol(Rm)){
         res1<-res1+t(Rm[,i])%*%solve(Vi)%*%Rm[,i]
                           }
     sigsquare<-res1/10000
     #res<--log(abs(det(Vi)))*15-log(res1)*1500
	 res<--log(det(Vi)*error^2)*100/2-res1/(error^2*2)
     return(-res)
         }
maximum_t<-nlm(f=loglhd_t,p=c(0.2,2,0.3,1,-0.1,1,1),hessian=T)  
   


loglhd<-function(ps)
# beta1,beta2,beta3,beta4,phi1,phi12,phi13,phi2,phi23,phi3
# phi are the variance components, not the variance-covariance matrix
{    
    
     beta_fix<-c(ps[1],ps[2],ps[3],ps[4])
     M<-15
     N<-240
     n<-16
     #phi<-matrix(c(ps[5],ps[6],ps[7],ps[6],ps[8],ps[9],ps[7],ps[9],ps[10]),3,3,byrow=T)
     L11<-ps[5]
     L21<-ps[6]/L11
     l22<-ps[8]^2-L21^2
     if(l22<=0) l22<-0.001
     L22<-sqrt(l22)
     L31<-ps[7]/L11
     L32<-(ps[9]-L31*L21)/L22
     l33<-ps[10]^2-L31^2-L32^2
     if (l33<0) l33<-0
     L33<-sqrt(l33)
     L<-matrix(c(L11,0,0,L21,L22,0,L31,L32,L33),nrow=3,ncol=3,byrow=T)
     phi<-L%*%t(L)
	 error<-ps[11]
     X<-cbind(rep(1,N),list$length,list$time,list$temp)
     y<-list$dat
     R<-y-X%*%beta_fix
     Rm<-matrix(R,nrow=n,ncol=M,byrow=F)
     Zi<-cbind(rep(1,n),list$time[1:n],list$temp[1:n])
     Vi<-diag(n)+Zi%*%phi%*%t(Zi)
     if (det(Vi)==0) Vi<-Vi+diag(0.01,n,n)
     res1<-0
     for (i in 1:ncol(Rm)){
         res1<-res1+t(Rm[,i])%*%solve(Vi)%*%Rm[,i]
                           }
     sigsquare<-res1/N
     #Z<-kronecker(diag(1,M),Zi)
     #D1<-kronecker(diag(1,M),phi)
     #v<-diag(N)+Z%*%D1%*%t(Z)
     #res<--log(abs(det(Vi)))*M/2-log(res1)*N/2
	 deter<-det(Vi)*error^2
	 #if (deter<=0) deter<-0.001 
	 res<--log(deter)*M/2-res1/(2*error^2)
     return(-res)
         }
 
maximum<-nlm(f=loglhd,p=c(-0.2,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9,1),hessian=T)  
varcov<-solve(maximum$hessian)
ci<-matrix(NA,nrow=4,ncol=2)
for (i in 1:4){
      ci[i,1]<-maximum$estimate[i]+qt(0.025,res$fixDF$X[i])*sqrt(varcov[i,i])
	  ci[i,2]<-maximum$estimate[i]+qt(0.975,res$fixDF$X[i])*sqrt(varcov[i,i])
			  }
			  
profile_beta1<-function(theta,ps)
# beta2,beta3,beta4,phi1,phi12,phi13,phi2,phi23,phi3
{
     M<-15
     N<-240
     n<-16
     beta_fix<-c(theta,ps[1],ps[2],ps[3])
     L11<-ps[4]
     L21<-ps[5]/L11
     l22<-ps[7]^2-L21^2
     if (l22<=0) l22<-0.001
     L22<-sqrt(l22)
     L31<-ps[6]/L11
     L32<-(ps[8]-L31*L21)/L22
     l33<-ps[9]^2-L31^2-L32^2
     if (l33<0) l33<-0
     L33<-sqrt(l33)
     L<-matrix(c(L11,0,0,L21,L22,0,L31,L32,L33),nrow=3,ncol=3,byrow=T)
     phi<-L%*%t(L)
     X<-cbind(rep(1,N),list$length,list$time,list$temp)
     y<-list$dat
     R<-y-X%*%beta_fix
     Rm<-matrix(R,nrow=n,ncol=M,byrow=F)
     Zi<-cbind(rep(1,n),list$time[1:n],list$temp[1:n])
     Vi<-diag(n)+Zi%*%phi%*%t(Zi)
     tryRes <- try(det(Vi)==0)
     if(is(tryRes,'try-error')){
            print(Vi)
     }else{
     if (det(Vi)==0) Vi<-Vi+diag(0.01,n,n)
     res1<-0
     for (i in 1:ncol(Rm)){
         res1<-res1+t(Rm[,i])%*%solve(Vi)%*%Rm[,i]
                           }
     res<--log(abs(det(Vi)))*M/2-log(res1)*N/2
     return(-res)
      }
 }

prof_beta1<-function(beta1){
  like_beta1<-NULL
  for (i in 1:length(beta1)){
     like_beta1[i]<-nlm(f=profile_beta1,p=c(2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),theta=beta1[i],hessian=T)$minimum
  }
  norm_beta1<-exp(-like_beta1+maximum$minimum)
 return(norm_beta1-exp(-qchisq(0.95,1)/2))
}
uniroot(prof_beta1,c(0.4,0.6))
[1] 0.5663474

uniroot(prof_beta1,c(-2,-0.8))
[1] -1.077680

beta1<-seq(-1.5,1,0.01)
res0 <-matrix(NA,nrow=length(beta1),ncol=2)
 for (i in 1:length(beta1)) {
  max_profile1<-nlm(f=profile_beta1,p=c(2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),theta=beta1[i],hessian=T)
  res0[i,1] <- max_profile1$minimum
  res0[i,2] <- beta1[i]
 }
plot(res0[,2],exp(-res0[,1]+maximum$minimum),xlab="beta1",ylab="normed likelihood")
abline(h=exp(-qchisq(0.95,1)/2))

traceback()

profile_betaid<-function(ps,replace.idx=5,theta)
# beta2,beta3,beta4,phi1,phi12,phi13,phi2,phi23,phi3
{
     M <- 15
     N <- 240
     n <- 16
     ps[replace.idx]<- theta
     beta_fix<-ps[1:4]
     ##beta_fix[replace.idx] <- theta
     L11<-ps[5]
     if (L11==0) L11<-L11+0.01
     L21<-ps[6]/L11
     l22<-ps[8]^2-L21^2
     if (l22<=0) l22<-0.01
     L22<-sqrt(l22)
     L31<-ps[7]/L11
     L32<-(ps[9]-L31*L21)/L22
     l33<-ps[10]^2-L31^2-L32^2
     if (l33<=0) l33<-0.01
     L33<-sqrt(l33)
     L<-matrix(c(L11,0,0,L21,L22,0,L31,L32,L33),nrow=3,ncol=3,byrow=T)
     phi<-L%*%t(L)
     X<-cbind(rep(1,N),list$length,list$time,list$temp)
     y<-list$dat
     R<-y-X%*%beta_fix
     Rm<-matrix(R,nrow=n,ncol=M,byrow=F)
     Zi<-cbind(rep(1,n),list$time[1:n],list$temp[1:n])
     Vi<-diag(n)+Zi%*%phi%*%t(Zi)
     if (det(Vi)==0) Vi<-Vi+diag(0.01,n,n)
     res1<-0
     for (i in 1:ncol(Rm)){
         res1<-res1+t(Rm[,i])%*%solve(Vi)%*%Rm[,i]
                           }
     res<--log(abs(det(Vi)))*M/2-log(res1)*N/2
     return(-res)
      }

prof_betaid1<-function(betaid){
 ## like_betaid<-NULL
 ## cat(length(betaid),'\n')
 ## for (i in 1:length(betaid)){
 ##    like_betaid[i]<-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),theta=betaid[i],hessian=T)$minimum
 ## }
like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=1,theta=betaid,hessian=T)$minimum
norm_betaid<-exp(-like_betaid+maximum$minimum)
return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }
uniroot(prof_betaid1,c(0.4,0.6))
[1] 0.5662083

uniroot(prof_betaid1,c(-1.5,-0.9))
[1] -1.073223

prof_betaid2<-function(betaid){
   like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=2,theta=betaid,hessian=T)$minimum
   norm_betaid<-exp(-like_betaid+maximum$minimum)
   return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }
uniroot(prof_betaid2,c(0.5,0.7))
[1] 0.6199054

uniroot(prof_betaid2,c(2,4))
[1] 3.405206

prof_betaid3<-function(betaid){
   like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=3,theta=betaid,hessian=T)$minimum
   norm_betaid<-exp(-like_betaid+maximum$minimum)
   return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }
uniroot(prof_betaid3,c(0,0.05))
[1] 0.001695459

uniroot(prof_betaid3,c(0.5,2))
[1] 1.104511

prof_betaid4<-function(betaid){
   like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=4,theta=betaid,hessian=T)$minimum
   norm_betaid<-exp(-like_betaid+maximum$minimum)
   return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }
uniroot(prof_betaid4,c(-0.5,0.1))
[1] -0.4760996

uniroot(prof_betaid4,c(0.6,1))
[1] 0.9097114

prof_betaid5<-function(betaid){
   like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=5,theta=betaid,hessian=T)$minimum
   norm_betaid<-exp(-like_betaid+maximum$minimum)
   return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }


uniroot(prof_betaid5,c(-1,-0.5))
[1] -0.8024155

uniroot(prof_betaid5,c(0.5,1))
[1] 0.8024155

prof_betaid6<-function(betaid){
   like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=6,theta=betaid,hessian=T)$minimum
   norm_betaid<-exp(-like_betaid+maximum$minimum)
   return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }
uniroot(prof_betaid6,c(0,1))
[1] 0.6116615
uniroot(prof_betaid6,c(-0.5,0))
[1] -0.2432511

prof_betaid7<-function(betaid){
   like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=7,theta=betaid,hessian=T)$minimum
   norm_betaid<-exp(-like_betaid+maximum$minimum)
   return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }
uniroot(prof_betaid7,c(0,1))
[1] 0.5047433

uniroot(prof_betaid7,c(-1,0))
[1] -0.8677278

prof_betaid8<-function(betaid){
   like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=8,theta=betaid,hessian=T)$minimum
   norm_betaid<-exp(-like_betaid+maximum$minimum)
   return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }

beta8<-seq(0,2,0.01)
res8<-matrix(NA,nrow=length(beta8),ncol=2)
 for (i in 1:length(beta8)) {
  max_profile8<-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=8,theta=beta8[i],hessian=T)
  res8[i,1] <- max_profile8$minimum
  res8[i,2] <- beta8[i]
 }
plot(res8[,2],exp(-res8[,1]+maximum$minimum),xlab="beta8",ylab="normed likelihood")
abline(h=exp(-qchisq(0.95,1)/2))


uniroot(prof_betaid8,c(0,2))
[1] 1.424451

uniroot(prof_betaid8,c(-10,10))

prof_betaid9<-function(betaid){
   like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=9,theta=betaid,hessian=T)$minimum
   norm_betaid<-exp(-like_betaid+maximum$minimum)
   return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }
uniroot(prof_betaid9,c(-1,1))
[1] -0.5998066

uniroot(prof_betaid9,c(1,2))
[1] 1.499535

prof_betaid10<-function(betaid){
   like_betaid <-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=10,theta=betaid,hessian=T)$minimum
   norm_betaid<-exp(-like_betaid+maximum$minimum)
   return(norm_betaid-exp(-qchisq(0.95,1)/2))
         }
#uniroot(prof_betaid10,c(0,1.5))
$root
[1] 1.283342

uniroot(prof_betaid10,c(1.5,2))
[1] 1.595103

beta10<-seq(1,2,0.01)
res10<-matrix(NA,nrow=length(beta10),ncol=2)
 for (i in 1:length(beta10)) {
  max_profile10<-nlm(f=profile_betaid,p=c(1,2,0.5,0.2,0.2,0.1,0.05,0.7,0.3,0.9),replace.idx=10,theta=beta10[i],hessian=T)
  res10[i,1] <- max_profile10$minimum
  res10[i,2] <- beta10[i]
 }
plot(res10[,2],exp(-res10[,1]+maximum$minimum),xlab="beta10",ylab="normed likelihood")
abline(h=exp(-qchisq(0.95,1)/2))

loglhd_s<-function(ps)
# beta1,beta2,beta3,beta4,phi1,phi12,phi13,phi2,phi23,phi3
# phi are the variance components, not the variance-covariance matrix
{
     #browser()
     beta_fix<-c(ps[1],ps[2],ps[3],ps[4])
     M<-15
     N<-240
     n<-16
     L11<-ps[5]
	 if (L11<=0) L11<-0.001
     L21<-ps[6]/L11
	 if (ps[8]<0) ps[8]<-0.001
     l22<-ps[8]^2-L21^2
     if(l22<=0) l22<-0.001
     L22<-sqrt(l22)
     L31<-ps[7]/L11
     L32<-(ps[9]-L31*L21)/L22
	 if (ps[10]<0) ps[10]<-0.001
     l33<-ps[10]^2-L31^2-L32^2
     if (l33<0) l33<-0.001
     L33<-sqrt(l33)
     L<-matrix(c(L11,0,0,L21,L22,0,L31,L32,L33),nrow=3,ncol=3,byrow=T)
     phi<-L%*%t(L)
	 if (ps[11]<0) ps[11]<-0.001
	 error<-ps[11]
     X<-cbind(rep(1,N),list$length,list$time,list$temp)
     y<-list$dat
     R<-y-X%*%beta_fix
     Rm<-matrix(R,nrow=n,ncol=M,byrow=F)
     Zi<-cbind(rep(1,n),list$time[1:n],list$temp[1:n])
     Vi<-diag(n)+Zi%*%phi%*%t(Zi)
	 if (det(Vi)==0) Vi<-Vi+diag(0.001,n,n) else
	     if (1/det(Vi)<2.220446e-16) Vi<- diag(0.001,n,n)
     res1<-0
     for (i in 1:ncol(Rm)){
         res1<-res1+t(Rm[,i])%*%solve(Vi)%*%Rm[,i]
                           }
     #sigsquare<-res1/N
     #res<--log(abs(det(Vi)))*M/2-log(res1)*N/2
	 res<--log(det(Vi))*M/2-log(error^2)*N/2-res1/(2*error^2)
     return(-res)
         }
 
#maximum_s<-nlm(f=loglhd_s,p=c(-0.23,2.16,0.55,0.21,0.21,0.13,0.056,0.65,0.27,0.9,0.9),hessian=T)  
maximum_s<-nlm(f=loglhd_s,p=c(-0.2,2,0.5,0.2,0.25,0.1,0.05,0.7,0.3,0.8,1),gradtol = 1e-13,iterlim = 200,hessian=T)  
est<-maximum_s$estimate

maximum_s<-nlm(f=loglhd_s,p=c(-0.2,2,0.5,0.2,0.3,0.1,0.1,0.7,0.3,0.8,1),hessian=T)  


##------------------------------------------------------------------------------------------------------------------------------------------------

mean<-c(0,0,0)
sigma<-diag(3)
sigma[1,1] <-1
sigma[1,2] <-sigma[2,1]<-0
sigma[2,2] <-1
sigma[3,3]<-1
sigma[1,3]<-sigma[3,1]<-0.9
sigma[2,3]<-sigma[3,2]<-0.1
set.seed(1000)
beta<- rmvnorm(15, mean, sigma)
set.seed(1000)
error<-rnorm(240,0,1)
rand.beta<-matrix(NA,nrow=240,ncol=3)
for (i in 1:3)rand.beta[,i]<-rep(beta[,i],each=16)
fix.beta<-2
x1<-seq(0.1,0.85,0.25)
x2<-seq(0.25,0.85,0.2)
x3<-seq(0.1,0.8,0.05)
time<-rep(rep(x1,4),15)
temp<-rep(rep(x2,each=4),15)
length<-rep(x3,each=16)
dat<-1+rand.beta[,1]+rand.beta[,2]*time+rand.beta[,3]*temp+fix.beta*length+error
subject<-rep(seq(1,15,1),each=16)
list<-data.frame(cbind(subject,time,temp,length,dat))
list$subject<-as.factor(list$subject)
res<-lme(dat~length, data=list, random=~time+temp|subject,method="ML")

yy<-matrix(0,nrow=(16+3),ncol=15)
for(i in 1:15){
yy[1:16,i]<-dat[((i-1)*16+1): (i*16)]
}
xx<-matrix(0,nrow=(16+3),ncol=15)
for(i in 1:15){
xx[1:16,i]<-length[((i-1)*16+1): (i*16)]
}
Zmat<-cbind(rep(1,16),list$time[1:16],list$temp[1:16])

loglike<-function(ps){
    beta<-c(ps[1],ps[2])
	delta<-matrix(c(ps[3:5],0,ps[6],ps[7],0,0,ps[8]),nrow=3,ncol=3,byrow=T)
	zz<-rbind(Zmat,delta)
	zmat_qr<-qr(zz)
    Ri<-qr.R(zmat_qr,complete=TRUE)
    Qi<-qr.Q(zmat_qr,complete=TRUE)
    R11<-Ri[1:3,]
    tmp<-matrix(0,nrow=15*16,ncol=3)
    for(i in 1:15){
		a<-matrix(0,nrow=(16+3),ncol=2)
		a[1:16,1]<-1
		a[,2]<-xx[,i]
		a<-t(Qi)%*%a
		tmp[((i-1)*16+1):(i*16),1:2]<-a[4:19,]
		a<-t(Qi)%*%yy[,i]
		a<-as.vector(a)
		tmp[((i-1)*16+1):(i*16),3]<-a[4:19]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
	error<-ps[9]
    result<- -15*8*log(error^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-beta,1),"F")^2/(2*error^2)
	return(-result)
}
max_loglike<-nlm(loglike,c(-0.4,5,2,-0.2,-1.6,1.2,-0.009,1,1))	
est<-max_loglike$estimate
delta<-matrix(c(est[3:5],0,est[6:7],0,0,est[8]),nrow=3,byrow=T)
Dinv<-t(delta)%*%delta
Sigmahat<-solve(Dinv)*est[9]^2
sdev<-sqrt(diag(Sigmahat))
tmp<-diag(1/sdev)
corrm<- tmp%*%Sigmahat%*%tmp
    
loglike_profile<-function(ps){
   
    delta<-matrix(c(ps[1:3],0,ps[4],ps[5],0,0,ps[6]),nrow=3,ncol=3,byrow=T)
	zz<-rbind(Zmat,delta)
	zmat_qr<-qr(zz)
    Ri<-qr.R(zmat_qr,complete=TRUE)
    Qi<-qr.Q(zmat_qr,complete=TRUE)
    R11<-Ri[1:3,]
    tmp<-matrix(0,nrow=15*16,ncol=3)
    for(i in 1:15){
		a<-matrix(0,nrow=(16+3),ncol=2)
		a[1:16,1]<-1
		a[,2]<-xx[,i]
		a<-t(Qi)%*%a
		tmp[((i-1)*16+1):(i*16),1:2]<-a[4:19,]
		a<-t(Qi)%*%yy[,i]
		a<-as.vector(a)
		tmp[((i-1)*16+1):(i*16),3]<-a[4:19]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
	c.minus<-Rtmp[,3]
	m<-length(c.minus)
	c.minus<-c.minus[3:m]
    c.minus.norm<- sqrt(sum(c.minus*c.minus))
    beta.hat<-solve(Rtmp[1:2,1:2])%*%Rtmp[1:2,3]
    sigma.squared.hat<- c.minus.norm^2/15*16
    result<- -15*16*log(c.minus.norm) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))
	return(-result)
}
max_pro<-nlm(loglike_profile,c(2,-0.2,-1.6,1.2,-0.009,1))	
     
D<-matrix(c(0.9,0.07,0.8,0.07,0.74,0.008,0.8,0.008,1),nrow=3,byrow=T)
 d<-solve(D)
 del<-chol(d)





##------------------------------------------------------------------------------------------------------------------------------------------------
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
	
est<-maximum$estimate
delta<-matrix(c(est[4:5],0,est[6]),nrow=2,byrow=T)
Dinv<-t(delta)%*%delta
Sigmahat<-solve(Dinv)*est[7]^2
sdev<-sqrt(diag(Sigmahat))
tmp<-diag(1/sdev)
corrm<- tmp%*%Sigmahat%*%tmp
    
#uf<-qnorm(0.2)*sqrt(0.48^2+0.05^2*0.62^2+0.1*0.634*0.48*0.62)+1.9+2.1*0.05+0.1*2.25

full_loglike<-function(ps){
    #browser()
    uf<-2
	x1<-0.05
	x2<-0.1
    beta<-c(ps[1],ps[2],ps[3])
	error<-ps[6]
	e1<-qnorm(1-ps[7])
	l1<--x1*ps[4]+sqrt(x1^2*ps[4]^2-x1^2*(ps[4]^2+ps[5]^2)+((uf-ps[1]-ps[2]*x1-ps[3]*x2)/(e1*error))^2)
	L<-matrix(c(l1,ps[4],0,ps[5]),nrow=2,ncol=2,byrow=F)
	delta<-solve(t(L))
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
max_tran<-nlminb(c(2,2,2,1,0.5,1,0.6),full_loglike,lower=c(rep(0,4),0.01,0.01,0.01),upper=c(rep(Inf,6),0.99))

est<-max_tran$par
l1<--0.05*est[4]+sqrt(0.05^2*est[4]^2-0.05^2*(est[4]^2+est[5]^2)+((2-est[1]-est[2]*0.05-est[3]*0.1)/(est[6]*qnorm(1-est[7])))^2)
L<-matrix(c(l1,0,est[4],est[5]),nrow=2,byrow=T)
delta<-solve(t(L))
Dinv<-t(delta)%*%delta
Sigmahat<-solve(Dinv)*est[6]^2
sdev<-sqrt(diag(Sigmahat))
tmp<-diag(1/sdev)
corrm<- tmp%*%Sigmahat%*%tmp

prof_loglike<-function(theta,ps){
    #browser()
    uf<-2
	x1<-0.05
	x2<-0.1
    beta<-c(ps[1],ps[2],ps[3])
	error<-ps[6]
	e1<-qnorm(1-theta)
	l1<--x1*ps[4]+sqrt(x1^2*ps[4]^2-x1^2*(ps[4]^2+ps[5]^2)+((uf-ps[1]-ps[2]*x1-ps[3]*x2)/(error*e1))^2)
	L<-matrix(c(l1,ps[4],0,ps[5]),nrow=2,ncol=2,byrow=T)
	delta<-solve(t(L))
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
prof_loglike(theta=0.5001,ps=c(2,2,2,0.3,0.6,1))

prof_fx<-function(fx){
  #browser()
  like_fx<-NULL
  for (i in 1:length(fx)){
     #like_fx[i]<-nlm(f=prof_loglike,p=c(2,2,2,0.3,0.6,1),theta=fx[i],hessian=T)$minimum
	 #like_fx[i]<-nlminb(start=c(2,2,2,0.3,0.6,1),f=prof_loglike,theta=fx[i],lower=c(rep(0,3),rep(0.01,3)))$objective
	 like_fx[i]<-optim(c(2,2,2,0.3,0.6,1),fn=prof_loglike,theta=fx[i],lower=c(rep(0,3),rep(0.01,3)))$value
  }
  norm_fx<-exp(-like_fx+max_tran$objective)
 return(norm_fx-exp(-qchisq(0.95,1)/2))
}

fx<-seq(0.5002,0.51,0.001)
prof_fx(0.50002)

res0 <-matrix(NA,nrow=length(fx),ncol=2)
 for (i in 1:length(fx)) {
  #browser()
  max_profile1<-nlm(f=prof_loglike,p=c(2,2,2,0.3,0.6,1),theta=fx[i],hessian=T)
  res0[i,1] <- max_profile1$minimum
  res0[i,2] <- fx[i]
 }
plot(res0[,2],exp(-res0[,1]+max_tran$objective),xlab="F(x)",ylab="normed likelihood")
abline(h=exp(-qchisq(0.95,1)/2))

uniroot(prof_fx,c(0.8,0.99))
#0.9763
uniroot(prof_fx,c(0.5001,0.51))
#0.02368

uf<-2
x1<-0.05
x2<-0.1
full_loglike2<-function(ps){
    #browser()
    beta<-c(ps[1],ps[2],ps[3])
	error<-ps[6]
	k<-qnorm(1-ps[7])
	if (k==0) k<-0.00001
	const<-uf-beta%*%c(1,x1,x2)
	const2<-const^2*ps[4]^2-error^2*k^2
	if (const2<=0) res <- -1e15 else
	{
	a3<-(ps[4]*x1+ps[5])*error*k/sqrt(const2)
	delta<-matrix(c(ps[4:5],0,a3),nrow=2,ncol=2,byrow=T)
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
	
	}
	return(-res)
}

max_tran2<-nlminb(c(2,2,2,3,-1,1,0.7),full_loglike2,lower=c(rep(0,3),0,-Inf,0.001,0.001),upper=c(rep(Inf,5),Inf,0.999))
max_tran3<-optim(c(2,2,2,3,-1,1,0.3),fn=full_loglike2,lower=c(rep(0,3),0,-Inf,0.001,0.001),upper=c(rep(Inf,5),Inf,0.999),hessian=T,method="L-BFGS-B")
A<-uf-est[1]-est[2]*x1-est[3]*x2
k<-sqrt(A^2*est[4]^2/((((est[4]*x1+est[5])^2/est[6]^2)+1)*est[7]^2))
F<-1-pnorm(-k)
full_loglike2(c(est[1:5],est[7],F))
full_loglike2(c(2,2,2,3,-1,1,0.3))

ineq<-function(ps){
	A<-uf-ps[1]-ps[2]*x1-ps[3]*x2
	k<-qnorm(1-ps[7])
	eqn1<-A^2*ps[4]^2-ps[6]^2*k^2
	#eqn2<-ps[7]
	#eqn3<-ps[6]
	return(eqn1)
	#return(c(eqn1,eqn2,eqn3))
	}
result<-solnp(pars=c(2,2,2,3,0,1,0.8),fun=full_loglike2,ineqfun=ineq,ineqLB=c(0.0001,0.0001,0),ineqUB=c(20,0.9999,2))
ans<-gosolnp(pars=c(2,2,2,3,0,1,0.8),fun=full_loglike2,LB=c(0,0,0,0,-Inf,0,0.0001),UB=c(Inf,Inf,Inf,Inf,Inf,2,0.9999),ineqfun=ineq,ineqLB=c(0.0001,0.0001,0.0001),ineqUB=c(10,0.9999,2))

est<-max_tran2$par
a3<-(est[4]*x1+est[5])*est[6]*qnorm(1-est[7])/sqrt((uf-c(est[1:3])%*%c(1,x1,x2))^2*est[4]^2-est[6]^2*(qnorm(1-est[7]))^2)
delta<-matrix(c(est[4:5],0,-a3),nrow=2,ncol=2,byrow=T)
Dinv<-t(delta)%*%delta
Sigmahat<-solve(Dinv)*est[6]^2
sdev<-sqrt(diag(Sigmahat))
tmp<-diag(1/sdev)
corrm<- tmp%*%Sigmahat%*%tmp


prof_loglike2<-function(theta,ps){
    #browser()
    beta<-c(ps[1],ps[2],ps[3])
	error<-ps[6]
	k<-qnorm(1-theta)	
	if (k==0) k<-k+0.00001
	const<-uf-beta%*%c(1,x1,x2)
	const2<-const^2*ps[4]^2-error^2*k^2
	if (const2<=0) {result<--1e15} else
	{
	a3<-(ps[4]*x1+ps[5])*error*k/sqrt(const2)
	delta<-matrix(c(ps[4:5],0,a3),nrow=2,ncol=2,byrow=T)
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
	}
	return(-result)
}

prof_fx2<-function(fx){
  #browser()
  like_fx<-NULL
  for (i in 1:length(fx)){
     #like_fx[i]<-nlm(f=prof_loglike2,p=c(2,2,2,3,-1,1),theta=fx[i],hessian=T)$minimum
	 like_fx[i]<-nlminb(start=c(2,2,2,3,-1,1),prof_loglike2,theta=fx[i],lower=c(rep(0,4),-Inf,0.01))$objective
	 #like_fx[i]<-optim(c(2,2,2,3,-1,1),fn=prof_loglike2,theta=fx[i],lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")$value
  }
  norm_fx<-exp(-like_fx+max_tran2$objective)
 return(norm_fx-exp(-qchisq(0.95,1)/2))
}
fx<-seq(0.001,0.499,0.001)
res0 <-matrix(NA,nrow=length(fx),ncol=2)
 for (i in 1:length(fx)) {
  #browser()
  max_profile1<-optim(c(2,2,2,3,-1,1),fn=prof_loglike2,theta=fx[i],lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
  res0[i,1] <- max_profile1$value
  res0[i,2] <- fx[i]
 }
plot(res0[,2],exp(-res0[,1]+max_tran2$objective),xlab="F(x)",ylab="normed likelihood")
abline(h=exp(-qchisq(0.95,1)/2))


uniroot(prof_fx2,c(0.2,0.4))
[1] 0.27

uniroot(prof_fx2,c(0.3,0.49))
[1] 0.4790949


##########################
stars.l <- seq(from=star[7],to=0,length.out=10)
#opt_in<-optim(star[-7],fn=prof_loglike2,theta=stars.l[2],lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
opt_in<-nlminb(start=star[-7],prof_loglike2,theta=stars.l[2],lower=c(rep(0,4),-Inf,0.01))
yhat.l<-c(1,exp(-opt_in$objective+max_tran2$objective))
pars<-opt_in$par
for (st in stars.l[c(-1,-2)]){
  #res_opt<-optim(pars[-7],fn=prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
  res_opt<-nlminb(start=pars[-7],prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01))
    obj <- res_opt$objective
  yhat.l <-c(yhat.l,exp(-obj+max_tran2$objective))
  pars <- res_opt$par
}

stars.r <- seq(from=star[7],to=1,length.out=200)
#opt_in<-optim(star[-7],fn=prof_loglike2,theta=stars.r[2],lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
opt_in<-nlminb(start=star[-7],prof_loglike2,theta=stars.r[2],lower=c(rep(0,4),-Inf,0.01))
yhat.r<-c(1,exp(-opt_in$objective+max_tran2$objective))
pars<-opt_in$par
lst <- list()
lst <- c(opt_in)
for (st in stars.r[c(-1,-2)]){
  #res_opt<-optim(pars[-7],fn=prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
  res_opt<-nlminb(start=pars,prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01))
  lst <- c(lst,res_opt)
  obj <- res_opt$objective
  #if(exp(-res+max_tran2$objective)<0.5) browser()
  yhat.r <-c(yhat.r,exp(-obj+max_tran2$objective))
  pars <- res_opt$par
}
plot(c(stars.l,stars.r[-1]),c(yhat.l,yhat.r[-1]),main="20 points in total",xlab="F(x)",ylab="norm likelihood")
abline(h=exp(-qchisq(0.95,1)/2))  
lines(c(rev(stars.l),stars.r[-1]),c(rev(yhat.l),yhat.r[-1]))

nlminb(start=lst[[6*13+1]],prof_loglike2,theta=stars.r[16],lower=c(rep(0,4),-Inf,0.01))

#############################################################################
stars.l <- seq(from=1-star[7],to=0.5,length.out=10)
#opt_in<-optim(star[-7],fn=prof_loglike2,theta=stars.l[2],lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
opt_in<-nlminb(start=star[-7],prof_loglike2,theta=stars.l[2],lower=c(rep(0,4),-Inf,0.01))
yhat.l<-c(1,exp(-opt_in$objective+max_tran2$objective))
pars<-opt_in$par
for (st in stars.l[c(-1,-2)]){
  #res_opt<-optim(pars[-7],fn=prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
  res_opt<-nlminb(start=pars[-7],prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01))
    obj <- res_opt$objective
  yhat.l <-c(yhat.l,exp(-obj+max_tran2$objective))
  pars <- res_opt$par
}

stars.r <- seq(from=1-star[7],to=1,length.out=11)
#opt_in<-optim(star[-7],fn=prof_loglike2,theta=stars.r[2],lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
opt_in<-nlminb(start=star[-7],prof_loglike2,theta=stars.r[2],lower=c(rep(0,4),-Inf,0.01))
yhat.r<-c(1,exp(-opt_in$objective+max_tran2$objective))
pars<-opt_in$par
lst <- list()
lst <- c(opt_in)
for (st in stars.r[c(-1,-2)]){
  #res_opt<-optim(pars[-7],fn=prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01),method="L-BFGS-B")
  res_opt<-nlminb(start=pars,prof_loglike2,theta=st,lower=c(rep(0,4),-Inf,0.01))
  lst <- c(lst,res_opt)
  obj <- res_opt$objective
  #if(exp(-res+max_tran2$objective)<0.5) browser()
  yhat.r <-c(yhat.r,exp(-obj+max_tran2$objective))
  pars <- res_opt$par
}
plot(c(stars.l,stars.r[-1]),c(yhat.l,yhat.r[-1]),main="20 points in total",xlab="F(x)",ylab="norm likelihood")
abline(h=exp(-qchisq(0.95,1)/2))  
lines(c(rev(stars.l),stars.r[-1]),c(rev(yhat.l),yhat.r[-1]))
	
#concentrate out the linear parameters
full_loglike3<-function(ps){
    #browser()
    beta<-c(ps[1],ps[2],ps[3])
	error<-ps[6]
	k<-qnorm(1-ps[7])
	if (k==0) k<-0.00001
	const<-uf-beta%*%c(1,x1,x2)
	const2<-const^2*ps[4]^2-error^2*k^2
	if (const2<=0) res <- -1e15 else
	{
	a3<-(ps[4]*x1+ps[5])*error*k/sqrt(const2)
	delta<-matrix(c(ps[4:5],0,a3),nrow=2,ncol=2,byrow=T)
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
	
	}
	return(-res)
}

