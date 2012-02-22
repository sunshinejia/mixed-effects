library(mvtnorm)
library(lme4)
#library(nlme)
library(MASS)
library(MCMCpack)
library(magic) 
meanb<-c(2,2)
sigma<-diag(2)
sigma[1,1] <-1
sigma[1,2] <-sigma[2,1]<-0.6
sigma[2,2] <-1
fix.beta<-2
set.seed(254)#set.seed(254) 
X1<- seq(0.05,by=0.05,length=20)
#X1<-X1-mean(X1)
X2<-seq(0.1,0.8,0.05)
#X2<-X2-mean(X2)
X1<-rep(X1,15)
X2<-rep(X2,each=20)
uf<-2
x1<-0.05
x2<-0.1
se_Zob<-function(a1,a2,a3,a4,a5,a6,a7){
    a4.sqrt <- sqrt(a4)
    a5.sqrt <- sqrt(a5)
    a45.sqrt <- sqrt(a4*a5)
	matrix1<-matrix(c(2*a4.sqrt,a6*a5.sqrt,a6*a5.sqrt,0),nrow=2)
    matrix2<-matrix(c(0,a6*a4.sqrt,a6*a4.sqrt,2*a5.sqrt),nrow=2)
    matrix3<-matrix(c(0,a4.sqrt*a5.sqrt,a4.sqrt*a5.sqrt,0),nrow=2)
	dersigma1<-Zmat%*%matrix1%*%t(Zmat)
	dersigma2<-Zmat%*%matrix2%*%t(Zmat)
	dersigma3<-Zmat%*%matrix3%*%t(Zmat)
	Sigmahat<-matrix(c(a4,a6*a45.sqrt,a6*a45.sqrt,a5),nrow=2)
	V<-a7^2*diag(20)+Zmat%*%Sigmahat%*%t(Zmat)
	invV<-solve(V)
    H11<-matrix(rep(0,9),nrow=3);offH11<-offH12<-offH13<-offH21<-offH22<-offH23<-offH31<-offH32<-offH33<-0
	diagH11<-diagH12<-diagH13<-diagH21<-diagH22<-diagH23<-diagH31<-diagH32<-diagH33<-0
	for (i in 1:15){
		tempH<--t(Xmat[(20*i-19):(20*i),])%*%invV%*%Xmat[(20*i-19):(20*i),]
		H11<-tempH+H11
		tempV<-invV%*%residual[,i]
		temp11<--t(Xmat[(20*i-19):(20*i),1])%*%invV%*%dersigma1%*%tempV
		offH11<-offH11+temp11
		temp12<--t(Xmat[(20*i-19):(20*i),1])%*%invV%*%dersigma2%*%tempV
		offH12<-offH12+temp12
		temp13<--t(Xmat[(20*i-19):(20*i),1])%*%invV%*%dersigma3%*%tempV
		offH13<-offH11+temp13
		temp21<--t(Xmat[(20*i-19):(20*i),2])%*%invV%*%dersigma1%*%tempV
		offH21<-offH21+temp21
		temp22<--t(Xmat[(20*i-19):(20*i),2])%*%invV%*%dersigma2%*%tempV
		offH22<-offH22+temp22
		temp23<--t(Xmat[(20*i-19):(20*i),2])%*%invV%*%dersigma3%*%tempV
		offH23<-offH23+temp23
		temp31<--t(Xmat[(20*i-19):(20*i),3])%*%invV%*%dersigma1%*%tempV
		offH31<-offH31+temp31
		temp32<--t(Xmat[(20*i-19):(20*i),3])%*%invV%*%dersigma2%*%tempV
		offH32<-offH32+temp32
		temp33<--t(Xmat[(20*i-19):(20*i),3])%*%invV%*%dersigma3%*%tempV
		offH33<-offH33+temp33
		tempV2<-invV%*%(2*residual[,i]%*%t(residual[,i])-V)%*%invV
		tempV3<-invV%*%(residual[,i]%*%t(residual[,i])-V)%*%invV
		dersigma12<-Zmat%*%matrix(c(0,a6,a6,0),nrow=2)%*%t(Zmat)
		dersigma13<-Zmat%*%matrix(c(0,a5.sqrt,a5.sqrt,0),nrow=2)%*%t(Zmat)
		dersigma21<-Zmat%*%matrix(c(0,a6,a6,0),nrow=2)%*%t(Zmat)
		dersigma23<-Zmat%*%matrix(c(0,a4.sqrt,a4.sqrt,0),nrow=2)%*%t(Zmat)
		dersigma31<-Zmat%*%matrix(c(0,a5.sqrt,a5.sqrt,0),nrow=2)%*%t(Zmat)
		dersigma32<-Zmat%*%matrix(c(0,a4.sqrt,a4.sqrt,0),nrow=2)%*%t(Zmat)
		diagH11<-diagH11-sum(diag(invV%*%dersigma1%*%tempV2%*%dersigma1))/2
		diagH12<-diagH12-sum(diag(invV%*%dersigma1%*%tempV2%*%dersigma2))/2+sum(diag(tempV3%*%dersigma12))/2
		diagH13<-diagH13-sum(diag(invV%*%dersigma1%*%tempV2%*%dersigma3))/2+sum(diag(tempV3%*%dersigma13))/2
		diagH21<-diagH21-sum(diag(invV%*%dersigma2%*%tempV2%*%dersigma1))/2+sum(diag(tempV3%*%dersigma21))/2
		diagH22<-diagH22-sum(diag(invV%*%dersigma2%*%tempV2%*%dersigma2))/2
		diagH23<-diagH23-sum(diag(invV%*%dersigma2%*%tempV2%*%dersigma3))/2+sum(diag(tempV3%*%dersigma23))/2
		diagH31<-diagH31-sum(diag(invV%*%dersigma3%*%tempV2%*%dersigma1))/2+sum(diag(tempV3%*%dersigma31))/2
		diagH32<-diagH32-sum(diag(invV%*%dersigma3%*%tempV2%*%dersigma2))/2+sum(diag(tempV3%*%dersigma32))/2
		diagH33<-diagH33-sum(diag(invV%*%dersigma3%*%tempV2%*%dersigma3))/2
	}
	offH<-matrix(c(offH11,offH12,offH13,offH21,offH22,offH23,offH31,offH32,offH33),nrow=3,ncol=3,byrow=T)
	diagH<-matrix(c(diagH11,diagH12,diagH13,diagH21,diagH22,diagH23,diagH31,diagH32,diagH33),nrow=3,ncol=3,byrow=T)
	L<--rbind(cbind(H11,offH),cbind(t(offH),diagH))
	fun1<-sqrt(a4+(x1^2)*a5+2*x1*a6*a45.sqrt)
	fun2<-uf-a1-a2*x1-a3*x2
    k<-fun2/fun1
	derkb0<--1/fun1
	derkb1<--x1/fun1
	derkb2<--x2/fun1
	derksig0<--(a4.sqrt+x1*a6*a5.sqrt)*fun2/(fun1^3)
    derksig1<--(x1^2*a5.sqrt+x1*a6*a4.sqrt)*fun2/(fun1^3)
    derkrho<--x1*a4.sqrt*a5.sqrt*fun2/(fun1^3)
	derk<-c(derkb0,derkb1,derkb2,derksig0,derksig1,derkrho)
	Sigmak<-derk%*%solve(L)%*%derk
	return(sqrt(Sigmak))
	}
se_Zob2<-function(b1,b2,b3,b4,b5){
#b1 is sigmab1^2,b2 is sigma,b3 is beta0,b4 is beta1, b5 is beta2	
	dersigma1<-Zmat2%*%t(Zmat2)*2*sqrt(b1)
	dersigma11<-Zmat2%*%t(Zmat2)*2
	V<-b2^2*diag(20)+Zmat2%*%t(Zmat2)*b1
	invV<-solve(V)
	H11<-matrix(rep(0,9),nrow=3);offH11<-offH21<-offH31<-diagH<-0
	for (i in 1:15){
		tempH<--t(Xmat[(20*i-19):(20*i),])%*%solve(V)%*%Xmat[(20*i-19):(20*i),]
		H11<-tempH+H11
		tempV<-invV%*%residual[,i]
		temp11<--t(Xmat[(20*i-19):(20*i),1])%*%invV%*%dersigma1%*%tempV
		offH11<-offH11+temp11
		temp21<--t(Xmat[(20*i-19):(20*i),2])%*%invV%*%dersigma1%*%tempV
		offH21<-offH21+temp21
		temp31<--t(Xmat[(20*i-19):(20*i),3])%*%invV%*%dersigma1%*%tempV
		offH31<-offH31+temp31
		tempV2<-invV%*%(2*residual[,i]%*%t(residual[,i])-V)%*%invV
		tempV3<-invV%*%(residual[,i]%*%t(residual[,i])-V)%*%invV
		diagH<-diagH-sum(diag(invV%*%dersigma1%*%tempV2%*%dersigma1))/2+sum(diag(tempV3%*%dersigma11))/2
	}
	offH<-matrix(c(offH11,offH21,offH31),nrow=3)
	L<-rbind(cbind(H11,offH),c(t(offH),diagH))
	fun1<-x1*sqrt(b1)
	fun2<-uf-b3-b4*x1-b5*x2
    k<-fun2/fun1
	derkb0<--1/fun1
	derkb1<--x1/fun1
	derkb2<--x2/fun1
	derksig1<--fun2/(x1*b1)
	derk<-c(derkb0,derkb1,derkb2,derksig1)
	Sigmak<-derk%*%solve(L)%*%derk
	return(sqrt(Sigmak))
}	
k_fun<-function(a1,a2,a3,a4,a5,a6){
	fun1<-sqrt(a4+(x1^2)*a5+2*x1*a6*sqrt(a4*a5))
	fun2<-uf-a1-a2*x1-a3*x2
    k<-fun2/fun1
	return(k)
	}
	
SIMUL_N <- 2
int_boot <- matrix(NA,nrow=SIMUL_N,ncol=2);seed1<-list();seed2<-list()
for(inter in 1:SIMUL_N){
seed1<- .Random.seed
beta<- rmvnorm(15, meanb, sigma)
seed2<-.Random.seed
error<-rnorm(300,0,1)
rand.beta<-matrix(NA,nrow=300,ncol=2)
for (i in 1:2)rand.beta[,i]<-rep(beta[,i],each=20)
dat<-rand.beta[,1]+rand.beta[,2]*X1+fix.beta*X2+error
subject<-rep(seq(1,15,1),each=20)
list1<-data.frame(cbind(subject,X1,X2,dat))
list1$subject<-as.factor(list1$subject)
res<-lmer(dat~1+X1+X2+(1+X1|subject),list1,REML =F)
Xmat<-cbind(rep(1,300),list1$X1,list1$X2)
Zmat<-cbind(rep(1,20),list1$X1[1:20])
Zmat2<-as.matrix(list1$X1[1:20])
residual<-matrix(NA,nrow=20,ncol=15)
for(i in 1:15){
residual[,i]<-list1$dat[((i-1)*20+1): (i*20)]-Xmat[(20*i-19):(20*i),]%*%c(beta0_est,beta1_est,beta2_est)
}
beta0_est<-res@fixef[1]
beta1_est<-res@fixef[2]
beta2_est<-res@fixef[3]
temp_est<-lme4::VarCorr(res)
sigma_est<-attr(temp_est,"sc")
var0_est<-temp_est$subject[1,1]
var1_est<-temp_est$subject[2,2]
corr_est<-attr(temp_est$subject,"correlation")[1,2]
meanb_sim<-c(beta0_est,beta1_est)
sigma_sim<-diag(2)
sigma_sim[1,1] <-var0_est
sigma_sim[1,2] <-sigma_sim[2,1]<-temp_est$subject[1,2]
sigma_sim[2,2] <-var1_est
fix.beta_sim<-beta2_est

se_hat<-se_Zob(beta0_est,beta1_est,beta2_est,var0_est,var1_est,corr_est,sigma_est)
k_hat<-k_fun(beta0_est,beta1_est,beta2_est,var0_est,var1_est,corr_est)
k_true<-k_fun(meanb[1],meanb[2],fix.beta,sigma[1,1],sigma[2,2],sigma[1,2])

options(warn = 2)
N <- 3000
t_boot_res <- numeric()
for (s in 1:N){
	beta_sim <- rmvnorm(15 , meanb_sim, sigma_sim)
	error_sim<-rnorm(300,0,sigma_est)
	rand.beta_sim<-matrix(NA,nrow=300,ncol=2)
	for (i in 1:2)rand.beta_sim[,i]<-rep(beta_sim[,i],each=20)
	dat_sim<-rand.beta_sim[,1]+rand.beta_sim[,2]*X1+fix.beta_sim*X2+error_sim
	list1_sim<-data.frame(cbind(subject,X1,X2,dat_sim))
	list1_sim$subject<-as.factor(list1$subject)
	res_sim<-lmer(dat_sim~1+X1+X2+(1+X1|subject),list1_sim,REML =F)
	beta0_boot<-res_sim@fixef[1]
	beta1_boot<-res_sim@fixef[2]
	beta2_boot<-res_sim@fixef[3]
	temp_boot<-lme4::VarCorr(res_sim)
	sigma_boot<-attr(temp_boot,"sc")
	var0_boot<-temp_boot$subject[1,1]
	var1_boot<-temp_boot$subject[2,2]
	corr_boot<-temp_boot$subject[1,2]/sqrt(var0_boot*var1_boot)
	if (var0_boot<1e-14) {
	var0_boot<-0;corr_boot<-0;
	se_boot<-se_Zob2(var1_boot,sigma_boot,beta0_boot,beta1_boot,beta2_boot)
	}else {
	se_boot<-se_Zob(beta0_boot,beta1_boot,beta2_boot,var0_boot,var1_boot,corr_boot,sigma_boot)
	}
	k_boot<-k_fun(beta0_boot,beta1_boot,beta2_boot,var0_boot,var1_boot,corr_boot)
	t_boot<-(k_boot-k_hat)/se_boot
	t_boot_res <- c(t_boot_res, t_boot)
	print(s)
	 }

tstat<-(k_hat-k_true)/se_hat
int_boot[inter,]<-k_hat-se_hat*quantile(t_boot_res,c(.975,.025))
}
int_boot
#[1] -1.050047  1.151586
[1] -1.042959  1.170223
k_true
 -0.2910428
 
 
 
 ###########################plot bad data#########################3
 library(ggplot2)
dts <- read.csv("C:/Jia/research/baddata.csv", header = TRUE)
dts <- dts[,-1]
ggplot(dts, aes(x = X1, y = dat_sim2)) + stat_smooth(method = lm) + geom_point() + facet_wrap(~ subject)
coeflst <- by(dts, dts$subject, function(x){
  coef(lm(dat_sim2~X1, data = x))
})
res_coef<- do.call(rbind, coeflst)
var(res_coef[,1]);var(res_coef[,2]);cor(res_coef[,1],res_coef[,2])
lmer(dat_sim2~1+X1+(1+X1|subject),dts,REML =F)

dat_full<-dts[,3]+X2*fix.beta_sim
dts[,4]<-dat_full
write.csv(list1_sim,"C:/Jia/research/dat_corr2.csv")
lmer(dat_full~1+X1+X2+(1+X1|subject),dts,REML =F)
yy<-matrix(0,nrow=(20+2),ncol=15)
for(i in 1:15){
   yy[1:20,i]<-dts[,4][((i-1)*20+1): (i*20)]
          }
Xmat1<-matrix(NA, nrow=20,ncol=45)
for (i in 1:15){
     Xmat1[,(3*i-2):(i*3)]<-Xmat[(20*i-19):(20*i),]
	}
xx<-list(NA)
for( i in 1:15){
    xx[[i]]<-rbind(Xmat1[,(3*i-2):(i*3)],matrix(rep(0,6),nrow=2,ncol=3))
	           }
			   
loglike<-function(ps){
    betafix<-c(ps[1],ps[2],ps[3])
	delta<-matrix(c(ps[4:5],0,ps[6]),nrow=2,ncol=2,byrow=T)
	zz<-rbind(Zmat,delta)
	zmat_qr<-qr(zz)
    Ri<-qr.R(zmat_qr,complete=TRUE)
    Qi<-qr.Q(zmat_qr,complete=TRUE)
    R11<-Ri[1:2,]
    tmp<-matrix(0,nrow=15*20,ncol=4)
	a<-as.list(rep(0, 15))
	Qi.t <- t(Qi)
    for(i in 1:15){
		a[[i]]<-Qi.t%*%xx[[i]]
		tmp[((i-1)*20+1):(i*20),1:3]<-a[[i]][3:22,]
		b<- Qi.t %*%yy[,i]
		b <- b[,,drop = TRUE]
		## b <-as.vector(b)
		tmp[((i-1)*20+1):(i*20),4]<-b[3:22]
     }
	tmp<-qr(tmp)
	Rtmp<-qr.R(tmp,complete=TRUE)
	err<-ps[7]
    result<- -15*10*log(err^2) + 15*log(abs(det(delta)))-15*log(abs(det(R11)))-norm(Rtmp%*%c(-betafix,1),"F")^2/(2*err^2)
	return(-result)
}
Rprof()
system.time(maximum<-nlm(loglike,c(2,2,2,1.7,0,1,1)))
Rprof(NULL)
est<-maximum$estimate
delta2<-matrix(c(est[4:5],0,est[6]),nrow=2,byrow=T)
Dinv<-t(delta2)%*%delta2
Sigmahat2<-solve(Dinv)*est[7]^2
sdev<-sqrt(diag(Sigmahat2))
tmp2<-diag(1/sdev)
corrm<- tmp2%*%Sigmahat2%*%tmp2

max2<-nlm(loglike,c(1.47,1.96,3.51,100£¬100,100,1))

max2<-nlm(loglike,c(1.47,1.96,3.51,0.82,-0.16,0.78,1))
est2<-max2$estimate
delta2<-matrix(c(est2[4:5],0,est2[6]),nrow=2,byrow=T)
Dinv2<-t(delta2)%*%delta2
Sigmahat2<-solve(Dinv2)*est2[7]^2
sdev2<-sqrt(diag(Sigmahat2))
tmp2<-diag(1/sdev2)
corrm2<- tmp2%*%Sigmahat2%*%tmp2
## how to profile 
Rprof() # begin
x <- 1:10
y <- rnorm(100000)
z <- seq(100)
Rprof(NULL) # end
summaryRprof()
## time of certain steps
system.time({
x<- 1:10
y<- rnorm(1000)
})

list1_sim$dat_sim2<-list1_sim$dat_sim-fix.beta_sim*list1_sim$X2
write.csv(list1_sim,"C:/Jia/research/dat_corr1.csv")
dts <- read.csv("C:/Jia/research/dat_corr1.csv", header = TRUE)
dts <- dts[,-1]
ggplot(dts, aes(x = X1, y =dat_sim2 )) + stat_smooth(method = lm) + geom_point() + facet_wrap(~ subject)
coeflst <- by(dts, dts$subject, function(x){
  coef(lm(dat_sim2~X1, data = x))
})
res_coef<- do.call(rbind, coeflst)
var(res_coef[,1]);var(res_coef[,2]);cor(res_coef[,1],res_coef[,2])
Vmatrix<-matrix(c(var(res_coef[,1]), cor(res_coef[,1],res_coef[,2])*sd(res_coef[,1])*sd(res_coef[,2]), cor(res_coef[,1],res_coef[,2])*sd(res_coef[,1])*sd(res_coef[,2]),var(res_coef[,2])),nrow=2)
pdMatrix(pdMat(Vmatrix))

fit<-lmer(dat_sim~1+X1+X2+(1+X1|subject),dts,REML =F)
res_t<-lme(dat_sim~X1+X2, data=dts, random=~X1|subject,method="ML",control=list(niterEM=3000)) 
 pdMatrix(res_t$modelStruct$reStruct)
lapply(pdMatrix(res_t$modelStruct$reStruct), "*", res_t$sigma^2)
plot(lme4::ranef(fit)$subject[,1],lme4::ranef(fit)$subject[,2])

######################diagnostic collinearity###############################
vif.mer <- function (fit) {
    ## adapted from rms::vif
    
    v <- vcov(fit)
    nam <- names(lme4::fixef(fit))

    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0) {
        v <- v[-(1:ns), -(1:ns), drop = FALSE]
        nam <- nam[-(1:ns)]
    }
    
    d <- diag(v)^0.5
    v <- diag(solve(v/(d %o% d)))
    names(v) <- nam
    v
}

kappa.mer <- function (fit,
                       scale = TRUE, center = FALSE,
                       add.intercept = TRUE,
                       exact = FALSE) {
    X <- fit@X
    nam <- names(lme4::fixef(fit))
    
    ## exclude intercepts
    nrp <- sum(1 * (nam == "(Intercept)"))
    if (nrp > 0) {
        X <- X[, -(1:nrp), drop = FALSE]
        nam <- nam[-(1:nrp)]
    }

    if (add.intercept) {
        X <- cbind(rep(1), scale(X, scale = scale, center = center))
        kappa(X, exact = exact)
    } else {
        kappa(scale(X, scale = scale, center = scale), exact = exact)
    }
}

colldiag.mer <- function (fit,
                          scale = TRUE, center = FALSE,
                          add.intercept = TRUE) {
    ## adapted from perturb::colldiag, method in Belsley, Kuh, and
    ## Welsch (1980).  look for a high condition index (> 30) with
    ## more than one high variance propotion.  see ?colldiag for more
    ## tips.
    result <- NULL
    if (center) 
        add.intercept <- FALSE
    if (is.matrix(fit) || is.data.frame(fit)) {
        X <- as.matrix(fit)
        nms <- colnames(fit)
    }
    else if (class(fit) == "mer") {
        nms <- names(lme4::fixef(fit))
        X <- fit@X
        if (any(grepl("(Intercept)", nms))) {
            add.intercept <- FALSE
        }
    }
    X <- X[!is.na(apply(X, 1, all)), ]

    if (add.intercept) {
        X <- cbind(1, X)
        colnames(X)[1] <- "(Intercept)"
    }
    X <- scale(X, scale = scale, center = center)

    svdX <- svd(X)
    svdX$d
    condindx <- max(svdX$d)/svdX$d
    dim(condindx) <- c(length(condindx), 1)

    Phi = svdX$v %*% diag(1/svdX$d)
    Phi <- t(Phi^2)
    pi <- prop.table(Phi, 2)
    colnames(condindx) <- "cond.index"
    if (!is.null(nms)) {
        rownames(condindx) <- nms
        colnames(pi) <- nms
        rownames(pi) <- nms
    } else {
        rownames(condindx) <- 1:length(condindx)
        colnames(pi) <- 1:ncol(pi)
        rownames(pi) <- 1:nrow(pi)
    }         

    result <- data.frame(cbind(condindx, pi))
    zapsmall(result)
}

maxcorr.mer <- function (fit,
                         exclude.intercept = TRUE) {
    so <- summary(fit)
    corF <- so@vcov@factors$correlation
    nam <- names(lme4::fixef(fit))

    ## exclude intercepts
    ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
    if (ns > 0 & exclude.intercept) {
        corF <- corF[-(1:ns), -(1:ns), drop = FALSE]
        nam <- nam[-(1:ns)]
    }
    corF[!lower.tri(corF)] <- 0
    maxCor <- max(corF)
    minCor <- min(corF)
    if (abs(maxCor) > abs(minCor)) {
        zapsmall(maxCor)
    } else {
        zapsmall(minCor)
    }
}


####################################check the bad data############################################################
graphdat<-function(dat){
#coeflst1<- by(sample1, sample1$subject, function(x){
#coef(lm(dat_sim2~X1, data = x))
#})
#res_coef1<- do.call(rbind, coeflst1)
#var(res_coef1[,1]);var(res_coef1[,2]);cor(res_coef1[,1],res_coef1[,2])
#fit1<-lmer(dat_sim~1+X1+X2+(1+X1|subject),sample1,REML =F)
dat2<-dat-X2*fix.beta_sim
samp<-data.frame(cbind(subject,X1,dat2))
samp$subject<-as.factor(samp$subject)
fit<-lmer(dat2~1+X1+(1+X1|subject),samp,REML =F)
#res_t1<-lme(dat_sim~X1+X2, data=sample1, random=~X1|subject,method="ML",control=list(niterEM=3000)) 
int<-lme4::coef(fit)$subject[,1];slope<-lme4::coef(fit)$subject[,2]
coef_fit<-as.data.frame(cbind(1:15,int,slope))
colnames(coef_fit)[1]<-"subject"
ggplot(samp, aes(x = X1, y = dat2)) + stat_smooth(method = lm) +geom_abline(data = coef_fit, aes(intercept = int, slope = slope)) +geom_point() + facet_wrap(~ subject) 
}

sample1 <- read.csv("C:/Jia/research/dat_corr1.csv", header = TRUE)
sample1 <- sample1[,-1]
graphdat(sample1$dat_sim)
sample2 <- read.csv("C:/Jia/research/dat_sim2.csv", header = TRUE)
sample2 <- sample2[,-1]
graphdat(sample2$dat_sim)
sample3<- read.csv("C:/Jia/research/dat_corr2.csv", header = TRUE)
sample3 <- sample3[,-1]
graphdat(sample3$dat_sim)
set.seed(5698)
beta_sim <- rmvnorm(15 , meanb_sim, sigma_sim)
error_sim<-rnorm(300,0,sigma_est)
rand.beta_sim<-matrix(NA,nrow=300,ncol=2)
for (i in 1:2)rand.beta_sim[,i]<-rep(beta_sim[,i],each=20)
dat_sim<-rand.beta_sim[,1]+rand.beta_sim[,2]*X1+fix.beta_sim*X2+error_sim
graphdat(dat_sim)

dat_sim2<-dat_sim-X2*fix.beta_sim
coeflst3<- by(sample1, sample1$subject, function(x){
coef(lm(dat_sim2~X1, data = x))
})
res_coef3<- do.call(rbind, coeflst3)