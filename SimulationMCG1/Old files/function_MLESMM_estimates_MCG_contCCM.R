# This script contains the code for 
# the function MLESMM_estimates_MCG_contCCM. 

# This function generates datasets according to the model 
# in Mason, Cantoni & Ghisletta (2021) with contamination generated 
# according to the CCM (central contamination model) in 
# the random effects, in the measurement error, and in the design matrix
# of the fixed effects.

# For each generated dataset the function also computes the MLE, 
# S and MM estimates with the function Roblme and saves the estimates,
# the estimated asymptotic variances and the number of outliers 
# in a list MLESMM

library(mvtnorm)    # needed for rmvnorm in function data_gen_MCG
library(robustbase) # needed for covMcd in function Robust_lme

MLESMM_estimates_MCG_contCCM=function(nrep=1,n=200,k=4,
                              pe=0,pb=0,px=0,
                              mec=0,mbc2=0,alphac=1,
                              no=0){

lbeta=2
ltheta=k

betahatmatMLE=matrix(0,nrow=nrep,ncol=lbeta)
betahatmatS=matrix(0,nrow=nrep,ncol=lbeta)
betahatmatMM=matrix(0,nrow=nrep,ncol=lbeta)

thetahatmatMLE=matrix(0,nrow=nrep,ncol=ltheta)
thetahatmatS=matrix(0,nrow=nrep,ncol=ltheta)
thetahatmatMM=matrix(0,nrow=nrep,ncol=ltheta)

asympvarbetaMLE=list()
asympvarbetaS=list()
asympvarbetaMM=list()

asympvarthetaMLE=list()
asympvarthetaS=list()
asympvarthetaMM=list()

no_outliers=matrix(0,nrow=nrep,ncol=5)

for (m in 1:nrep){
  # generating dataset
  # settings for data_gen_MCG are taken from settings MLESMM_estimates_MCG
  dat = data_gen_MCG_contCCM(n=n,k=k,
                             pe=pe,pb=pb,px=px,
                             mec=mec,mbc2=mbc2,alphac = alphac,no=no)
  
  # Roblme for MLE
  summaryMLE = Roblme(dat$Y,dat$X,dat$Z,E=NULL,L=dat$L,
                         rho="MLE",rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL)
  
  # Roblme for S
  summaryS = Roblme(dat$Y,dat$X,dat$Z,E=NULL,L=dat$L,
                   rho="biweight",rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL)
  
  # Roblme for MM
  summaryMM = Roblme(dat$Y,dat$X,dat$Z,E=NULL,L=dat$L,
                    rho="biweight",rhoMM="biweight",eps=1e-5,maxiter=100,eff=0.95,V0=NULL)
  
  betahatmatMLE[m,]=summaryMLE$fixedeffectsS[,1]
  betahatmatS[m,]=summaryS$fixedeffectsS[,1]
  betahatmatMM[m,]=summaryMM$fixedeffectsMM[,1]
  
  thetahatmatMLE[m,]=summaryMLE$summarythetaS[,1]
  thetahatmatS[m,]=summaryS$summarythetaS[,1]
  thetahatmatMM[m,]=summaryMM$summarythetaS[,1]
  
  asympvarbetaMLE[[m]]=summaryMLE$varbetaShat
  asympvarbetaS[[m]]=summaryS$varbetaShat
  asympvarbetaMM[[m]]=summaryMM$varbetaMMhat
  
  asympvarthetaMLE[[m]]=summaryMLE$varthetahat
  asympvarthetaS[[m]]=summaryS$varthetahat
  asympvarthetaMM[[m]]=summaryMM$varthetahat
  
  no_outliers[m,]=c(dat$nobi,dat$noei,dat$noxi,dat$noe,dat$nox)
  
  print(m)
}

MLE=list()
MLE[[1]]=betahatmatMLE
MLE[[2]]=asympvarbetaMLE
MLE[[3]]=thetahatmatMLE
MLE[[4]]=asympvarthetaMLE
names(MLE)=c("beta","varbeta","theta","vartheta")

Sest=list()
Sest[[1]]=betahatmatS
Sest[[2]]=asympvarbetaS
Sest[[3]]=thetahatmatS
Sest[[4]]=asympvarthetaS
names(Sest)=c("beta","varbeta","theta","vartheta")

MM=list()
MM[[1]]=betahatmatMM
MM[[2]]=asympvarbetaMM
MM[[3]]=thetahatmatMM
MM[[4]]=asympvarthetaMM
names(MM)=c("beta","varbeta","theta","vartheta")

no_outliers=data.frame(no_outliers)
names(no_outliers)=c("nobi","noei","noxi","noe","nox")

MLESMM_CCM=list()
MLESMM_CCM[[1]]=MLE
MLESMM_CCM[[2]]=Sest
MLESMM_CCM[[3]]=MM
MLESMM_CCM[[4]]=no_outliers
names(MLESMM_CCM)=c("MLE","S","MM","no_outliers")

return(MLESMM_CCM)

}
