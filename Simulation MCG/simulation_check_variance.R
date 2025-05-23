library(MASS) # needed to generate from multivariate normal
library(ks)
library(mvtnorm)        # needed for rmvnorm in function data_gen_MCG
library(robustbase)     # needed for covMcd in function Robust_lme
library(robustvarComp)  # needed for varComprob
source("function_data_gen_MCG_CCMind_VG.R")
source("Robust_lme.R")
source("function_MLESMMcTAUind_estimates_MCG_CCMind_VG.R")
source("scenarios_MCG_simulation_VG.R")
#########################################################################
# This function orthogonalizes the d-by-2 matrix x=[1 v]
# for a d-vector v
# It returns a vector w, such that X^TX=20*ID, for X=[1 w]

data_sphe = function(x, tol = 1e-06){
  # x is a d by 2 matrix of the form [1 v]
  # where v is a d-vector
  Nm = nrow(x) - 1
  MEANS = colMeans(x)
  x.c = as.matrix(sweep(x, 2, MEANS, "-"))
  SVD = svd(x.c, nu = 0)
  SV = SVD$d
  if (!is.null(tol)) {
    rank = sum(SVD$d > (SVD$d[1L] * tol))
    if (rank < ncol(x)) {
      SVD$v = SVD$v[, 1L:rank, drop = FALSE]
      SVD$d = SVD$d[1L:rank]
    }
  }
  SIGMAs = SVD$d/sqrt(Nm+1)
  TRANS = SVD$v *(1/SIGMAs)
  # TRANS <- SVD$v %*% diag(1/SIGMAs)
  RES = x.c %*% TRANS
  #attr(RES, "center") <- MEANS
  #attr(RES, "transform") <- TRANS
  #attr(RES, "backtransform") <- diag(SIGMAs) %*% t(SVD$v)
  # attr(RES, "SingularValues") <- SV
  RES
}


##############################################################################
# Preparing the values at entry
n=100 # number of observations
k=4   # dimension of each observation

# Setting breakdown point and cut-off constant for Tukey's biweight
r=0.5
c0=rhoconst(k,r,0.01,100)
b0=expecrho(k,c0)

# settings of variance parameter and fixed effects
sigmaeps=1  # variance of error
sigmabeta=1 # variance of random effects
alpha1=1    # fixed effect 1
alpha2=1    # fixed effect 2
fixeff=c(alpha1,alpha2)     # vector of fixed effects parameters
theta=c(sigmabeta,sigmaeps) # vector of covariance parameters

# random settings
set.seed(123)

# Diagonal variances of measurement covariance
#varerror=diag(rep(1,4))   # setting in Copt & Victoria-Feser (2006)

# THIS SETTING GIVES DIFFERENCES IN COVARIANCE OF THETAHAT
# BUT ITERATIONS RUN OUT OF BOUNDS

#varerror=diag(runif(k)) 
varerror=diag(c(1,4,9,16))       
#varerror=cbind(c(1,r,r^2,r^3),c(r,1,r,r^2),c(r^2,r,1,r),c(r^3,r^2,r,1))

# construction design matrix X=[1 x], such that X^TX=4*I
X=cbind(1,rnorm(k))
X=cbind(1,data_sphe(X)) # setting in Copt & Victoria-Feser (2006) 
#X=as.matrix(cbind(1:4,rnorm(4)))

# Setting for design matrix random effects
Z=rep(1,times=k)   # setting in Copt & Victoria-Feser (2006)
#Z=1:4               # THIS SETTING GIVES DIFFERENCES IN COVRIANCE OF BETAHAT

# matrices L in V=sigmabeta*L1+sigmaeps*L2
L1=Z%*%t(Z)       # covariance part due to random effect
L2=varerror       # covariance part due to measurement error
Vtheta=theta[1]*L1+theta[2]*L2
L=as.matrix(cbind(vec(L1),vec(L2)))

Llist=list()          
Llist[[1]]=L1
Llist[[2]]=L2


# Asymptotic variances betahat
# Expression from our theory
varbeta=constbetahat(k,c0)*
  solve(t(X)%*%solve(Vtheta)%*%X)
# Expression from Copt & Victoria-Feser (2006)
varbetaCopt=constbetahat(k,c0)*
  solve(t(X)%*%X) %*% t(X) %*% Vtheta %*% X %*% solve(t(X)%*%X)

# Asymptotic variances thetahat
# Expression from our theory
vartheta=(2*sigma1(k,c0)*
            solve(t(L)%*%(solve(Vtheta)%x%solve(Vtheta))%*%L)+
            sigma2(k,c0)*theta%*%t(theta))
# Expression from Copt & Victoria-Feser (2006)
varthetaCopt=2*sigma1(k,c0)*
  solve((t(L)%*%L))%*%t(L)%*%(Vtheta%x%Vtheta)%*%L%*%solve((t(L)%*%L))+
  sigma2(k,c0)*theta%*%t(theta)


###################################################################
# setting up iterations for the simulation

nrep=10000
betahat=matrix(0,nrow=nrep,ncol = length(fixeff))   # matrix for nrep rows of betahat's with biweight
thetahat=matrix(0,nrow=nrep,ncol = length(theta))   # matrix for nrep rows of thetahat's with biweight
MDscale=matrix(0,nrow=nrep,ncol = 1)                # vector for nrep values of scaling
MDscalenew=matrix(0,nrow=nrep,ncol = 1)                # vector for nrep values of scaling

betahatmatS=matrix(0,nrow=nrep,ncol=length(fixeff))
thetahatmatS=matrix(0,nrow=nrep,ncol=length(theta))
asympvarbetaS=list()
asympvarthetaS=list()
betahatmatSclaudio=matrix(0,nrow=nrep,ncol=length(fixeff))
thetahatmatSclaudio=matrix(0,nrow=nrep,ncol=length(theta))
asympvarbetaSclaudio=list()
asympvarthetaSclaudio=list()
lbeta = length(fixeff)
ltheta = length(theta)
##############################################################################
for (m in 1:nrep){

# Generating data
y=NULL
item=NULL
obs=NULL
Xdesign=NULL

for (i in 1:n){
  itemnew=matrix(i,nrow=k,ncol=1)  # recording item index
  obsnew=matrix(1:k,nrow=k,ncol=1) # vector of observation indices
  
  # generating y=X*alpha+beta*Z+epsilon
  ynew=X%*%fixeff+
    rnorm(1,sd=sqrt(sigmabeta))*Z+
    sqrt(sigmaeps)*mvrnorm(1,Sigma=varerror,mu=rep(0,times=k))

  # stacking new y, item, obs, Xdesign to previous ones
  y=rbind(y,ynew)
  item=rbind(item,itemnew)
  obs=rbind(obs,obsnew)
  Xdesign=rbind(Xdesign,X)
}

# dataframe with n observations (y,X) in dimension kx(k*q)
simdata=data.frame(item=item,obs=obs,y=y,x1=Xdesign[,1],x2=Xdesign[,2])
XX = list()   
for (i in (1:n)){
  XX[[i]] <- X
}
Y=list()
for (i in 1:n){
  Y[[i]]=simdata[,3][simdata[,1]==i]
}
Ymat=matrix(0,n,k)
for (i in (1:n)){Ymat[i,] = Y[[i]]}

Zmat=Z
summaryS = tryCatch(
  expr  = {est0 =Roblme(Ymat,XX,Zmat,E=NULL,L=Llist,
                        rho="biweight",rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL)}, error  =  function(cond) {
                          fixedeffectsS = matrix(rep(NA,lbeta*4),lbeta,4)
                          fixedeffectsMM = matrix(rep(NA,lbeta*4),lbeta,4)
                          summarythetaS = matrix(rep(NA,ltheta*ltheta),ltheta,ltheta)
                          varbetaShat = matrix(rep(NA,lbeta*lbeta),lbeta,lbeta)
                          varbetaMMhat = matrix(rep(NA,lbeta*lbeta),lbeta,lbeta)
                          varthetahat = matrix(rep(NA,ltheta*ltheta),ltheta,ltheta)
                          w= rep(NA,n)
                          dis=rep(NA,n)
                          iterS =NA
                          iterM =NA
                          list(fixedeffectsS=fixedeffectsS,fixedeffectsMM=fixedeffectsMM,summarythetaS=summarythetaS,varbetaShat=varbetaShat,
                               varbetaMMhat=varbetaMMhat,varthetahat=varthetahat,
                               w=w,dis=dis,iterS=iterS,iterM=iterM)})

n_suj=n #nb of cases
n_mes=k #nb maximum number of measurements per case
J=n_mes
time<-as.numeric(rep(X[,2],n_suj))
id<-sort(as.numeric(rep(1:n_suj,n_mes)))
n_obs = n_suj*n_mes #nb of rows
groups <- cbind(rep(1:J, each=n),rep((1:n), J)) # a numeric matrix with two columns used to group the observations according to participant.

# Build the argument "varcov" of the varComprob() function
z1 = rep(1, J) #Value for intercept (=1) for the J observations by clusters
#z2 = unique(time) # Value for the time variable
#z1=Z 
K <- list() # the "varcov" object
K[[1]] <- tcrossprod(z1,z1) # Matrix for intercept
#K[[2]] <- tcrossprod(z2,z2) # Matrix for time variable
#K[[3]] <- tcrossprod(z1,z2) + tcrossprod(z2,z1) # Matrix of interaction Intercept by time variable
#names(K) = c("sigma2_Intercept", "sigma2_Time", "Covariance")
names(K) = c("sigma2_Intercept")
K <<- K 
y=vec(t(Ymat))
Dataset=data.frame(y,time,groups)
summarySclaudio=tryCatch(
  expr  = {est0 =varComprob(y ~ 1 +  time, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf),method = "S", psi = "bisquare")) }, error  =  function(cond) {
    beta = rep(NA,lbeta)
    eta = rep(NA,3)
    vcov.beta = matrix(rep(NA,lbeta*lbeta),lbeta,lbeta)
    eta0 = NA
    list(beta=beta,eta=eta,vcov.beta=vcov.beta,eta0=eta0)})

betahatmatS[m,]=summaryS$fixedeffectsS[,1]
thetahatmatS[m,]=summaryS$summarythetaS[,1]
asympvarbetaS[[m]]=summaryS$varbetaShat
asympvarthetaS[[m]]=summaryS$varthetahat

betahatmatSclaudio[m,]=summarySclaudio$beta
thetahatmatSclaudio[m,]=c(summarySclaudio$eta[1],summarySclaudio$eta[3],summarySclaudio$eta[2],summarySclaudio$eta0)
asympvarbetaSclaudio[[m]]=summarySclaudio$vcov.beta
print(m)
}

betahat10000sim1=betahatmatS
thetahat10000sim1=thetahatmatS
write.csv(betahat10000sim1,file='betahat10000sim1.csv')
write.csv(thetahat10000sim1,file='thetahat10000sim1.csv')
varbetasim1=varbeta
varbetaCoptsim1=varbetaCopt


#betahat10000sim2=betahatmatSclaudio
#thetahat10000sim2=thetahatmatSclaudio
#write.csv(betahat10000sim2,file='betahat10000sim2.csv')
#write.csv(thetahat10000sim2,file='thetahat10000sim2.csv')
varbetasim2=varbeta
varbetaCoptsim2=varbetaCopt

###############################################
betahat=betahat10000sim1
varbeta=varbetasim1
varbetaCopt=varbetaCoptsim1

fixeffest=c(mean(betahat[,1]),mean(betahat[,2]))
varbetaest=cov.wt(betahat)$cov*n

# Constructing contours for beta covariance plot
beta1=seq(qnorm(0.001,mean=alpha1,sd=sqrt(varbetaCopt[1,1]/n)),
          qnorm(0.999,mean=alpha1,sd=sqrt(varbetaCopt[1,1]/n)),length=100)
beta2=seq(qnorm(0.001,mean=alpha2,sd=sqrt(varbetaCopt[2,2]/n)),
          qnorm(0.999,mean=alpha2,sd=sqrt(varbetaCopt[2,2]/n)),length=100)
zbeta=matrix(0,nrow=length(beta1),ncol=length(beta2))
zbetaCopt=matrix(0,nrow=length(beta1),ncol=length(beta2))
zbetaEmpirical=matrix(0,nrow=length(beta1),ncol=length(beta2))

for (i in 1:length(beta1)){
  for (j in 1:length(beta2)){
    betavec=c(beta1[i],beta2[j])
    zbeta[i,j]=mahalanobis(betavec,center = fixeff,cov = varbeta/n)
    zbetaCopt[i,j]=mahalanobis(betavec,center = fixeff,cov = varbetaCopt/n)
    zbetaEmpirical[i,j]=mahalanobis(betavec,center = fixeffest,cov = varbetaest/n)
  }
}

par(mfrow=c(2,2))
ymax=max(dnorm(alpha1,mean=alpha1,sd=sqrt(varbeta1/n)),
         dnorm(alpha1,mean=alpha1,sd=sqrt(varbetaCopt1/n)))
xmin=min(qnorm(0.001,mean=alpha1,sd=sqrt(varbeta1/n)),
         qnorm(0.001,mean=alpha1,sd=sqrt(varbetaCopt1/n)))
xmax=max(qnorm(0.999,mean=alpha1,sd=sqrt(varbeta1/n)),
         qnorm(0.999,mean=alpha1,sd=sqrt(varbetaCopt1/n)))
hist(betahat[,1],nclass=60,prob=TRUE,ylim=c(0,ymax+1),xlim=c(xmin,xmax),
     xlab="S-estimates for alpha1 with biweight",main="")
varbeta1=varbeta[1,1]
varbetaCopt1=varbetaCopt[1,1]
xas=seq(min(betahat[,1])-1,max(betahat[,1])+1,length=1000)
lines(xas,dnorm(xas,mean=alpha1,sd=sqrt(varbeta1/n)),
      col="darkgreen")
lines(xas,dnorm(xas,mean=alpha1,sd=sqrt(varbetaCopt1/n)),col=2)
legend("topleft",lty=1,legend=c("Our theory","Copt"),bty="n",
       col=c("darkgreen",2),lwd=2)

varbeta2=varbeta[2,2]
varbetaCopt2=varbetaCopt[2,2]
ymax=max(dnorm(alpha2,mean=alpha2,sd=sqrt(varbeta2/n)),
         dnorm(alpha2,mean=alpha2,sd=sqrt(varbetaCopt2/n)))
xmin=min(qnorm(0.001,mean=alpha2,sd=sqrt(varbeta2/n)),
         qnorm(0.001,mean=alpha2,sd=sqrt(varbetaCopt2/n)))
xmax=max(qnorm(0.999,mean=alpha2,sd=sqrt(varbeta2/n)),
         qnorm(0.999,mean=alpha2,sd=sqrt(varbetaCopt2/n)))
hist(betahat[,2],nclass=60,prob=TRUE,ylim=c(0,ymax+1),xlim=c(xmin,xmax),
     xlab="S-estimates for alpha2 with biweight",main="")

xas=seq(min(betahat[,2])-1,max(betahat[,2])+1,length=1000)
lines(xas,dnorm(xas,mean=alpha2,sd=sqrt(varbeta2/n)),col="darkgreen")
lines(xas,dnorm(xas,mean=alpha2,sd=sqrt(varbetaCopt2/n)),col=2)
legend("topleft",lty=1,legend=c("Our theory","Copt"),bty="n",
       col=c("darkgreen",2),lwd=2)

boxplot(betahat[,1],betahat[,2],names=c("alpha1","alpha2"),
        xlab="S-estimates for alpha1 and alpha2 with biweight")
abline(h=1,lty=2)

xmin=min(beta1,betahat[,1])
xmax=max(beta1,betahat[,1])
ymin=min(beta2,betahat[,2])
ymax=max(beta2,betahat[,2])
plot(betahat[,1],betahat[,2],
     xlab="alpha1",ylab="alpha2",
     main="S-estimates with biweight",col="darkgray",
     xlim=c(xmin,xmax),
     ylim=c(ymin,ymax))
contour(beta1,beta2,zbeta,levels = qchisq(0.95,df=4),add=TRUE,
        drawlabels = FALSE,col="darkgreen",xlab="beta1",ylab="beta2")
contour(beta1,beta2,zbetaCopt,levels = qchisq(0.95,df=4),add=TRUE,drawlabels = FALSE,col=2)
contour(beta1,beta2,zbetaEmpirical,levels = qchisq(0.95,df=4),add=TRUE,drawlabels = FALSE,col=1,lty = 3)
legend("topleft",lty=1,legend=c("Our theory","Copt"),bty="n",
       col=c("darkgreen",2),lwd=2)
par(mfrow=c(1,1))

###############################################
thetahat=thetahat10000sim2
vartheta=varthetasim2
varthetaCopt=varthetaCoptsim2

fixeffest=c(mean(betahat[,1]),mean(betahat[,2]))
varbetaest=cov.wt(betahat)$cov*n

mean(thetahat[,1])
mean(thetahat[,2])
theta

round(cov.wt(thetahat)$cov,10)
round(vartheta/n,10)
round(varthetaCopt/n,10)

# Constructing contours for theta covariance plot
theta1=seq(qnorm(0.001,mean=sigmabeta,sd=sqrt(varthetaCopt[1,1]/n)),
           qnorm(0.999,mean=sigmabeta,sd=sqrt(varthetaCopt[1,1]/n)),length=100)
theta2=seq(qnorm(0.001,mean=sigmaeps,sd=sqrt(varthetaCopt[2,2]/n)),
           qnorm(0.999,mean=sigmaeps,sd=sqrt(varthetaCopt[2,2]/n)),length=100)
ztheta=matrix(0,nrow=length(theta1),ncol=length(theta2))
zthetaCopt=matrix(0,nrow=length(theta1),ncol=length(theta2))
for (i in 1:length(theta1)){
  for (j in 1:length(theta2)){
    thetavec=c(theta1[i],theta2[j])
    ztheta[i,j]=mahalanobis(thetavec,center = theta,cov = vartheta/n)
    zthetaCopt[i,j]=mahalanobis(thetavec,center = theta,cov = varthetaCopt/n)
  }
}



par(mfrow=c(2,2))
ymax=max(dnorm(theta[1],mean=theta[1],sd=sqrt(vartheta1/n)),
         dnorm(theta[1],mean=theta[1],sd=sqrt(varthetaCopt1/n)))
xmin=min(qnorm(0.001,mean=theta[1],sd=sqrt(vartheta1/n)),
         qnorm(0.001,mean=theta[1],sd=sqrt(varthetaCopt1/n)))
xmax=max(qnorm(0.999,mean=theta[1],sd=sqrt(vartheta1/n)),
         qnorm(0.999,mean=theta[1],sd=sqrt(varthetaCopt1/n)))
hist(thetahat[,1],prob=TRUE,nclass = 60,
     xlab="S-estimates for theta1 with biweight",main="",xlim=c(xmin,xmax),ylim=c(0,ymax+1))
vartheta1=vartheta[1,1]
varthetaCopt1=varthetaCopt[1,1]
xas=seq(min(thetahat[,1])-1,max(thetahat[,1])+1,length=1000)
lines(xas,dnorm(xas,mean=theta[1],sd=sqrt(vartheta1/n)),
      col="darkgreen")
lines(xas,dnorm(xas,mean=theta[1],sd=sqrt(varthetaCopt1/n)),col=2)
legend("topleft",lty=1,legend=c("Our theory","Copt"),bty="n",
       col=c("darkgreen",2),lwd=2)

ymax=max(dnorm(theta[2],mean=theta[2],sd=sqrt(vartheta2/n)),
         dnorm(theta[2],mean=theta[2],sd=sqrt(varthetaCopt2/n)))
xmin=min(qnorm(0.001,mean=theta[1],sd=sqrt(vartheta2/n)),
         qnorm(0.001,mean=theta[1],sd=sqrt(varthetaCopt2/n)))
xmax=max(qnorm(0.999,mean=theta[1],sd=sqrt(vartheta2/n)),
         qnorm(0.999,mean=theta[1],sd=sqrt(varthetaCopt2/n)))
hist(thetahat[,2],prob=TRUE,nclass = 60,
     xlab="S-estimates for theta2 with biweight",main="",xlim=c(xmin,xmax),ylim=c(0,ymax+1))
vartheta2=vartheta[2,2]
varthetaCopt2=varthetaCopt[2,2]
xas=seq(min(thetahat[,2])-1,max(thetahat[,2])+1,length=1000)
lines(xas,dnorm(xas,mean=theta[2],sd=sqrt(vartheta2/n)),
      col="darkgreen")
lines(xas,dnorm(xas,mean=theta[2],sd=sqrt(varthetaCopt2/n)),
      col=2)
legend("topleft",lty=1,legend=c("Our theory","Copt"),bty="n",
       col=c("darkgreen",2),lwd=2)

boxplot(thetahat[,1],thetahat[,2],names=c("theta1","theta2"),
        xlab="S-estimates for alpha1 and alpha2 with biweight")
abline(h=1,lty=2)

xmin=min(theta1,thetahat[,1])
xmax=max(theta1,thetahat[,1])
ymin=min(theta2,thetahat[,2])
ymax=max(theta2,thetahat[,2])
plot(thetahat[,1],thetahat[,2],col="darkgray",
     xlab="theta1",ylab="theta2",main="S-estimates for with biweight",
     xlim=c(xmin,xmax),
     ylim=c(ymin,ymax))

contour(theta1,theta2,ztheta,levels = qchisq(0.90,df=4),add=TRUE,
        drawlabels = FALSE,col="darkgreen",xlab="theta1",ylab="theta2")
contour(theta1,theta2,zthetaCopt,levels = qchisq(0.90,df=4),add=TRUE,
        drawlabels = FALSE,col=2)
legend("topleft",lty=1,legend=c("Our theory","Copt"),bty="n",
       col=c("darkgreen",2),lwd=2)
par(mfrow=c(1,1))

betahat