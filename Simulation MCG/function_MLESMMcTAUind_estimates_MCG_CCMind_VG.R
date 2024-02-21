# This script contains the code for the function 
#   - MLESMMCOMPind_estimates_MCG_CCMind. 

# This function generates datasets according to the model 
# in Mason, Cantoni & Ghisletta (2021) with contamination generated 
# - ICM (independent contamination model or cellwise contamination) 
# - CCM (central contamination model or casewise contamination) 
# in the measurement error and according to 
# - CCM (central contamination model) 
# in the random effects.

# In addition, the function also generates contamination in 
# the design matrix of the fixed effects according to 
# - ICM (independent contamination model or cellwise contamination) 
# - CCM (central contamination model or casewise contamination) 

# For each generated dataset the function also computes the 
# MLE, S, MM, and COMP estimates with the functions 
# Robmle and varComprob 
# and saves the estimates, the estimated asymptotic variances 
# and the number of outliers  in a list MLESMMCOMP

library(mvtnorm)        # needed for rmvnorm in function data_gen_MCG
library(robustbase)     # needed for covMcd in function Robust_lme
library(robustvarComp)  # needed for varComprob


MLESMMTAUCOMPSTAUind_estimates_MCG_CCMind=function(nrep=1,n=200,k=4,pe=0,pb=0,px=0,
                                  mec=0,mbc2=0,alphac=1,
                                  randcont=0,COMPind=FALSE,CCMind=FALSE,
                                  Xa=FALSE,mux=0,Sclaudioind=FALSE){

lbeta=2   # X has 2 columns
ltheta=k

betahatmatMLE=matrix(0,nrow=nrep,ncol=lbeta)
betahatmatS=matrix(0,nrow=nrep,ncol=lbeta)
betahatmatMM=matrix(0,nrow=nrep,ncol=lbeta)
betahatmatTau=matrix(0,nrow=nrep,ncol=lbeta)

thetahatmatMLE=matrix(0,nrow=nrep,ncol=ltheta)
thetahatmatS=matrix(0,nrow=nrep,ncol=ltheta)
thetahatmatMM=matrix(0,nrow=nrep,ncol=ltheta)
thetahatmatTau=matrix(0,nrow=nrep,ncol=ltheta)

asympvarbetaMLE=list()
asympvarbetaS=list()
asympvarbetaMM=list()
asympvarbetaTau=list()

asympvarthetaMLE=list()
asympvarthetaS=list()
asympvarthetaMM=list()
asympvarthetaTau=list()

if (COMPind==TRUE){
  # Data preparation for varComprob
  n_suj=n #nb of cases
  n_mes=k #nb maximum number of measurements per case
  J=n_mes
  time<-as.numeric(rep(c(0:(n_mes-1)),n_suj))
  id<-sort(as.numeric(rep(1:n_suj,n_mes)))
  n_obs = n_suj*n_mes #nb of rows
  groups <- cbind(rep(1:J, each=n),rep((1:n), J)) # a numeric matrix with two columns used to group the observations according to participant.
  
  # Build the argument "varcov" of the varComprob() function
  z1 = rep(1, J) #Value for intercept (=1) for the J observations by clusters
  z2 = unique(time) # Value for the time variable
  
  K <- list() # the "varcov" object
  K[[1]] <- tcrossprod(z1,z1) # Matrix for intercept
  K[[2]] <- tcrossprod(z2,z2) # Matrix for time variable
  K[[3]] <- tcrossprod(z1,z2) + tcrossprod(z2,z1) # Matrix of interaction Intercept by time variable
  names(K) = c("sigma2_Intercept", "sigma2_Time", "Covariance")
  K <<- K  # This is Z%*%t(Z)
  
  betahatmatCOMPS=matrix(0,nrow=nrep,ncol=lbeta)
  thetahatmatCOMPS=matrix(0,nrow=nrep,ncol=ltheta)
  asympvarbetaCOMPS=list()
  asympvarthetaCOMPS=list()
  betahatmatCOMPTau=matrix(0,nrow=nrep,ncol=lbeta)
  thetahatmatCOMPTau=matrix(0,nrow=nrep,ncol=ltheta)
  asympvarbetaCOMPTau=list()
  asympvarthetaCOMPTau=list()
  
  if (Sclaudioind == TRUE){
    betahatmatSclaudio=matrix(0,nrow=nrep,ncol=lbeta)
    thetahatmatSclaudio=matrix(0,nrow=nrep,ncol=ltheta)
    asympvarbetaSclaudio=list()
    asympvarthetaSclaudio=list()}
  } # END of data preparation for varComprob

no_outliers=matrix(0,nrow=nrep,ncol=5)

for (m in 1:nrep){
  # generating dataset
  # settings for data_gen_MCG are taken from settings MLESMM_estimates_MCG
  dat =data_gen_MCG_CCMind(n=n,k=k,
                     pe=pe,pb=pb,px=px,
                     mec=mec,mbc2=mbc2,alphac = alphac,
                     randcont=randcont,CCMind=CCMind,
                     Xa=Xa,mux=mux)
  
  # Roblme for MLE
  summaryMLE = tryCatch(
    expr  = {est0 =Roblme(dat$Y,dat$X,dat$Z,E=NULL,L=dat$L,
                         rho="MLE",rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL)}, error  =  function(cond) {
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
  
  # Roblme for S
  summaryS = tryCatch(
    expr  = {est0 =Roblme(dat$Y,dat$X,dat$Z,E=NULL,L=dat$L,
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
  
  # Roblme for MM
  summaryMM = tryCatch(
    expr  = {est0 =Roblme(dat$Y,dat$X,dat$Z,E=NULL,L=dat$L,
                    rho="biweight",rhoMM="biweight",eps=1e-5,maxiter=100,eff=0.95,V0=NULL)}, error  =  function(cond) {
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

  # varComprob for TAU
  summaryTau=tryCatch(
    expr  = {est0 =varComprob(y ~ 1 +  time, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf),method = "Tau", psi = "bisquare")) }, error  =  function(cond) {
      beta = rep(NA,lbeta)
      eta = rep(NA,3)
      vcov.beta = matrix(rep(NA,lbeta*lbeta),lbeta,lbeta)
      eta0 = NA
      list(beta=beta,eta=eta,vcov.beta=vcov.beta,eta0=eta0)})
  

  betahatmatMLE[m,]=summaryMLE$fixedeffectsS[,1]
  betahatmatS[m,]=summaryS$fixedeffectsS[,1]
  betahatmatMM[m,]=summaryMM$fixedeffectsMM[,1]
  betahatmatTau[m,]=summaryTau$beta
  
  thetahatmatMLE[m,]=summaryMLE$summarythetaS[,1]
  thetahatmatS[m,]=summaryS$summarythetaS[,1]
  thetahatmatMM[m,]=summaryMM$summarythetaS[,1]
  thetahatmatTau[m,]=c(summaryTau$eta[1],summaryTau$eta[3],summaryTau$eta[2],summaryTau$eta0)
  
  asympvarbetaMLE[[m]]=summaryMLE$varbetaShat
  asympvarbetaS[[m]]=summaryS$varbetaShat
  asympvarbetaMM[[m]]=summaryMM$varbetaMMhat
  asympvarbetaTau[[m]]=summaryTau$vcov.beta
  
  asympvarthetaMLE[[m]]=summaryMLE$varthetahat
  asympvarthetaS[[m]]=summaryS$varthetahat
  asympvarthetaMM[[m]]=summaryMM$varthetahat
  
  if (COMPind==TRUE){
    # varComprob for composite TAU
    y=vec(t(dat$Y))
    Dataset=data.frame(y,time,groups)
    summaryCOMPTau=tryCatch(
      expr  = {est0 =varComprob(y ~ 1 +  time, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf),method = "compositeTau", psi = "bisquare")) }, error  =  function(cond) {
        beta = rep(NA,lbeta)
        eta = rep(NA,3)
        vcov.beta = matrix(rep(NA,lbeta*lbeta),lbeta,lbeta)
        eta0 = NA
        list(beta=beta,eta=eta,vcov.beta=vcov.beta,eta0=eta0)})
    
    betahatmatCOMPTau[m,]=summaryCOMPTau$beta
    thetahatmatCOMPTau[m,]=c(summaryCOMPTau$eta[1],summaryCOMPTau$eta[3],summaryCOMPTau$eta[2],summaryCOMPTau$eta0)
    asympvarbetaCOMPTau[[m]]=summaryCOMPTau$vcov.beta
    
    # varComprob for composite S
    summaryCOMPS=tryCatch(
      expr  = {est0 =varComprob(y ~ 1 +  time, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf),method = "compositeS", psi = "bisquare")) }, error  =  function(cond) {
        beta = rep(NA,lbeta)
        eta = rep(NA,3)
        vcov.beta = matrix(rep(NA,lbeta*lbeta),lbeta,lbeta)
        eta0 = NA
        list(beta=beta,eta=eta,vcov.beta=vcov.beta,eta0=eta0)})
    
    betahatmatCOMPS[m,]=summaryCOMPS$beta
    thetahatmatCOMPS[m,]=c(summaryCOMPS$eta[1],summaryCOMPS$eta[3],summaryCOMPS$eta[2],summaryCOMPS$eta0)
    asympvarbetaCOMPS[[m]]=summaryCOMPS$vcov.beta

    if (Sclaudio == TRUE){
      # varComprob for S
      summarySclaudio=tryCatch(
      expr  = {est0 =varComprob(y ~ 1 +  time, groups = groups, data = Dataset, varcov = K, control = varComprob.control(lower = c(0, 0, -Inf),method = "S", psi = "bisquare")) }, error  =  function(cond) {
        beta = rep(NA,lbeta)
        eta = rep(NA,3)
        vcov.beta = matrix(rep(NA,lbeta*lbeta),lbeta,lbeta)
        eta0 = NA
        list(beta=beta,eta=eta,vcov.beta=vcov.beta,eta0=eta0)})
    betahatmatSclaudio[m,]=summarySclaudio$beta
    thetahatmatSclaudio[m,]=c(summarySclaudio$eta[1],summarySclaudio$eta[3],summarySclaudio$eta[2],summarySclaudio$eta0)
    asympvarbetaSclaudio[[m]]=summarySclaudio$vcov.beta
    }
  }
  
  no_outliers[m,]=c(dat$nobi,dat$noei,dat$noxi,dat$noe,dat$nox)
  
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

Tau=list()
Tau[[1]]=betahatmatTau
Tau[[2]]=asympvarbetaTau
Tau[[3]]=thetahatmatTau
Tau[[4]]=asympvarthetaTau
names(Tau)=c("beta","varbeta","theta","vartheta")

if (COMPind==TRUE){
COMPS=list()
COMPS[[1]]=betahatmatCOMPS
COMPS[[2]]=asympvarbetaCOMPS
COMPS[[3]]=thetahatmatCOMPS
COMPS[[4]]=matrix(0,nrow=k,ncol=k)
names(COMPS)=c("beta","varbeta","theta","vartheta")

COMPTau=list()
COMPTau[[1]]=betahatmatCOMPTau
COMPTau[[2]]=asympvarbetaCOMPTau
COMPTau[[3]]=thetahatmatCOMPTau
COMPTau[[4]]=matrix(0,nrow=k,ncol=k)
names(COMPTau)=c("beta","varbeta","theta","vartheta")

if (Sclaudioind == TRUE){
  Sclaudio=list()
  Sclaudio[[1]]=betahatmatSclaudio
  Sclaudio[[2]]=asympvarbetaSclaudio
  Sclaudio[[3]]=thetahatmatSclaudio
  Sclaudio[[4]]=matrix(0,nrow=k,ncol=k)
  names(Sclaudio)=c("beta","varbeta","theta","vartheta")
  }
}

no_outliers=data.frame(no_outliers)
names(no_outliers)=c("nobi","noei","noxi","noe","nox")

if (COMPind==TRUE){
  if (Sclaudioind == FALSE){
  MLESMMTAUCOMPSTAU=list()
  MLESMMTAUCOMPSTAU[[1]]=MLE
  MLESMMTAUCOMPSTAU[[2]]=Sest
  MLESMMTAUCOMPSTAU[[3]]=MM
  MLESMMTAUCOMPSTAU[[4]]=Tau
  MLESMMTAUCOMPSTAU[[5]]=COMPS
  MLESMMTAUCOMPSTAU[[6]]=COMPTau
  MLESMMTAUCOMPSTAU[[7]]=no_outliers
  names(MLESMMTAUCOMPSTAU)=c("MLE","S","MM","Tau","COMPS","COMPTau","no_outliers")
  return(MLESMMTAUCOMPSTAU)
  }else{
    MLESMMTAUCOMPSTAUSclaudio=list()
    MLESMMTAUCOMPSTAUSclaudio[[1]]=MLE
    MLESMMTAUCOMPSTAUSclaudio[[2]]=Sest
    MLESMMTAUCOMPSTAUSclaudio[[3]]=MM
    MLESMMTAUCOMPSTAUSclaudio[[4]]=Tau
    MLESMMTAUCOMPSTAUSclaudio[[5]]=COMPS
    MLESMMTAUCOMPSTAUSclaudio[[6]]=COMPTau
    MLESMMTAUCOMPSTAUSclaudio[[7]]=Sclaudio
    MLESMMTAUCOMPSTAUSclaudio[[8]]=no_outliers
    names(MLESMMTAUCOMPSTAUSclaudio)=c("MLE","S","MM","Tau","COMPS","COMPTau","Sclaudio","no_outliers")
    return(MLESMMTAUCOMPSTAUSclaudio)
    }}else{MLESMMTAU=list()
        MLESMMTAU[[1]]=MLE
        MLESMMTAU[[2]]=Sest
        MLESMMTAU[[3]]=MM
        MLESMMTAU[[4]]=Tau
        MLESMMTAU[[5]]=no_outliers
        names(MLESMMTAU)=c("MLE","S","MM","Tau","no_outliers")
        return(MLESMMTAU)
        }

} # End of function MLESMMCOMPind_estimates_MCG_CCMind
