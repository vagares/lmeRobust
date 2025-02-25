# This script uses the contamination scenarios_CI and results of
# the corresponding simulation performed by
# Simulation_setting_model_MCG.R

# The script extracts the information from the .RData files
# corresponding to the different setups in dataframe 
# scenarios_CI and produces graph of the coverage probabilities

# Sets scenarios_CI equal to scenarios_CI
# scenarios_CI is used below to extract the .RData files

source("scenarios_MCG_simulation.R")

scenarios_CI=scenarios[scenarios$pe==0.05,]

if (max(scenarios_CI$mec)!=0){
  pname="mec";
  perc=unique(scenarios_CI$pe);
  levelsvec=mecvec}
if (max(scenarios_CI$mbc2)!=0){
  pname="mbc2";
  perc=unique(scenarios_CI$pb);
  levelsvec=mbc2vec}
if (max(scenarios_CI$alphac)!=1){
  pname="alphac";
  perc=unique(scenarios_CI$px);
  levelsvec=alphacvec}

#######################################################
# Extracting the information for beta from .Rdata files

# Dataframes for betahat and thetahat
CI_MLESMMbeta=NULL
CICP_MLESMMbeta=NULL
CP_MLESMMbeta=NULL

CI_MLESMMtheta=NULL
CICP_MLESMMtheta=NULL
CP_MLESMMtheta=NULL

for (i in 1:nrow(scenarios_CI)){
  nrep=scenarios_CI[i,1]
  nsample=scenarios_CI[i,2]
  ksample=scenarios_CI[i,3]
  pesample=scenarios_CI[i,4]
  pbsample=scenarios_CI[i,5]
  pxsample=scenarios_CI[i,6]
  mecsample=scenarios_CI[i,7]
  mbc2sample=scenarios_CI[i,8]
  alphacsample=scenarios_CI[i,9]
  rcsample=scenarios_CI[i,10]
  cTAUindsample=as.logical(scenarios_CI[i,11])
  CCMindsample=as.logical(scenarios_CI[i,12])
  Xasample=as.logical(scenarios_CI[i,13])
  muxsample=scenarios_CI[i,14]
  Sclaudiosample=as.logical(scenarios_CI[i,15])
  
  if (CCMindsample==TRUE){stop("Contains CCM scenarios_CI")}
  
  # Setting the filename depending on yes/no cTAU and yes/no CCM
  if (COMPindsample==FALSE){
    if (CCMindsample==FALSE){
      if(Xasample==TRUE){flnameEst="MLESMMTAU_ICM_Xa" }else{flnameEst="MLESMMTAU_ICM_Xf"}}else{
        if(Xasample==TRUE){flnameEst="MLESMMTAU_CCM_Xa" }else{flnameEst="MLESMMTAU_CCM_Xf"}}
  }else{if (Sclaudio == FALSE){
    if (CCMindsample==FALSE){
      if(Xasample==TRUE){flnameEst="MLESMMTAUCOMPSTAU_ICM_Xa"}else{flnameEst="MLESMMTAUCOMPSTAU_ICM_Xf"}}else{
        if(Xasample==TRUE){flnameEst="MLESMMTAUCOMPSTAU_CCM_Xa"}else{flnameEst="MLESMMTAUCOMPSTAU_CCM_Xf"}}
  }else
  {if (CCMindsample==FALSE){
    if(Xasample==TRUE){flnameEst="MLESMMTAUCOMPSTAUSclaudio_ICM_Xa"}else{flnameEst="MLESMMTAUCOMPSTAUSclaudio_ICM_Xf"}}else{
      if(Xasample==TRUE){flnameEst="MLESMMTAUCOMPSTAUSclaudio_CCM_Xa"}else{flnameEst="MLESMMTAUCOMPSTAUSclaudio_CCM_Xf"}}
  }   }
  
  if ((scenarios_CI$pe[i]==0)&
      (scenarios_CI$pb[i]==0)&
      (scenarios_CI$px[i]==0)){
    flname=paste0("./Results_Uncontaminated/",flnameEst,"_",
                  "nrep=",scenarios_CI$nrep[i],"_",
                  "n=",scenarios_CI$n[i],"_",
                  "k=",scenarios_CI$k[i],"_",
                  "pe=",scenarios_CI$pe[i],"_",
                  "pb=",scenarios_CI$pb[i],"_",
                  "px=",scenarios_CI$px[i],"_",
                  "mec=",scenarios_CI$mec[i],"_",
                  "mbc2=",scenarios_CI$mbc2[i],"_",
                  "alphac=",scenarios_CI$alphac[i],"_",
                  "rc=",scenarios_CI$rc[i],".RData")}
  
  if (scenarios_CI$pe[i]>0){
    flname=paste0("./Results_Epsilon_contamination/",flnameEst,"_",
                  "nrep=",scenarios_CI$nrep[i],"_",
                  "n=",scenarios_CI$n[i],"_",
                  "k=",scenarios_CI$k[i],"_",
                  "pe=",scenarios_CI$pe[i],"_",
                  "pb=",scenarios_CI$pb[i],"_",
                  "px=",scenarios_CI$px[i],"_",
                  "mec=",scenarios_CI$mec[i],"_",
                  "mbc2=",scenarios_CI$mbc2[i],"_",
                  "alphac=",scenarios_CI$alphac[i],"_",
                  "rc=",scenarios_CI$rc[i],".RData")}
  
  if (scenarios_CI$pb[i]>0){
    flname=paste0("./Results_Random_Effect_contamination/",flnameEst,"_",
                  "nrep=",scenarios_CI$nrep[i],"_",
                  "n=",scenarios_CI$n[i],"_",
                  "k=",scenarios_CI$k[i],"_",
                  "pe=",scenarios_CI$pe[i],"_",
                  "pb=",scenarios_CI$pb[i],"_",
                  "px=",scenarios_CI$px[i],"_",
                  "mec=",scenarios_CI$mec[i],"_",
                  "mbc2=",scenarios_CI$mbc2[i],"_",
                  "alphac=",scenarios_CI$alphac[i],"_",
                  "rc=",scenarios_CI$rc[i],".RData")}
  
  if (scenarios_CI$px[i]>0){
    flname=paste0("./Results_X_contamination/",flnameEst,"_",
                  "nrep=",scenarios_CI$nrep[i],"_",
                  "n=",scenarios_CI$n[i],"_",
                  "k=",scenarios_CI$k[i],"_",
                  "pe=",scenarios_CI$pe[i],"_",
                  "pb=",scenarios_CI$pb[i],"_",
                  "px=",scenarios_CI$px[i],"_",
                  "mec=",scenarios_CI$mec[i],"_",
                  "mbc2=",scenarios_CI$mbc2[i],"_",
                  "alphac=",scenarios_CI$alphac[i],"_",
                  "rc=",scenarios_CI$rc[i],".RData")}
  
  load(flname)
  
  #######################################################################
  # combining the extracted information for coverage probabiities 
  # for estimators for beta in a dataframe

  # dataframes for MLE beta
  CI_MLEbeta=NULL
  CICP_MLEbeta=NULL
  CP_MLEbeta=NULL

  for (j in 1:nrep){
    # Determine CI for beta1 and beta2 with MLE
    lwbbetaMLE=MLESMMTAUCOMPSTAUind$MLE$beta[j,]-qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$MLE$varbeta[[j]]))
    upbetaMLE=MLESMMTAUCOMPSTAUind$MLE$beta[j,]+qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$MLE$varbeta[[j]]))
    CI_MLEbeta=rbind(CI_MLEbeta,c(1,scenarios_CI[[pname]][i],lwbbetaMLE,upbetaMLE))
   
    # Check whether 250 and 10 are in MLE-CI
    CP_MLEbeta1=ifelse(lwbbetaMLE[1]<250,1,0)*ifelse(250<upbetaMLE[1],1,0)
    CP_MLEbeta2=ifelse(lwbbetaMLE[2]<10,1,0)*ifelse(10<upbetaMLE[2],1,0)
    CICP_MLEbeta=rbind(CICP_MLEbeta,c(1,scenarios_CI[[pname]][i],CP_MLEbeta1,CP_MLEbeta2))  
  }
  # Determine CP for CI-beta1 and CI-beta2 with MLE
  CP_MLEbeta=rbind(CP_MLEbeta,
                   c(1,scenarios_CI[[pname]][i],
                     sum(CICP_MLEbeta[,3])/nrep,sum(CICP_MLEbeta[,4])/nrep))
  
  # dataframes for S beta
  CI_Sbeta=NULL
  CICP_Sbeta=NULL
  CP_Sbeta=NULL
  
  for (j in 1:nrep){
    # Determine CI for beta1 and beta2 with S
    lwbbetaS=MLESMMTAUCOMPSTAUind$S$beta[j,]-qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$S$varbeta[[j]]))
    upbetaS=MLESMMTAUCOMPSTAUind$S$beta[j,]+qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$S$varbeta[[j]]))
    CI_Sbeta=rbind(CI_Sbeta,c(2,scenarios_CI[[pname]][i],lwbbetaS,upbetaS))
  
    # Check whether 250 and 10 are in S-CI
    CP_Sbeta1=ifelse(lwbbetaS[1]<250,1,0)*ifelse(250<upbetaS[1],1,0)
    CP_Sbeta2=ifelse(lwbbetaS[2]<10,1,0)*ifelse(10<upbetaS[2],1,0)
    CICP_Sbeta=rbind(CICP_Sbeta,c(2,scenarios_CI[[pname]][i],CP_Sbeta1,CP_Sbeta2))  
  }
  # Determine CP for CI-beta1 and CI-beta2 with S
  CP_Sbeta=rbind(CP_Sbeta,
                   c(2,scenarios_CI[[pname]][i],
                     sum(CICP_Sbeta[,3])/nrep,sum(CICP_Sbeta[,4])/nrep))
  
  # dataframes for MM beta
  CI_MMbeta=NULL
  CICP_MMbeta=NULL
  CP_MMbeta=NULL
  
  for (j in 1:nrep){
    # Determine CI for beta1 and beta2 with MM
    lwbbetaMM=MLESMMTAUCOMPSTAUind$MM$beta[j,]-qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$MM$varbeta[[j]]))
    upbetaMM=MLESMMTAUCOMPSTAUind$MM$beta[j,]+qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$MM$varbeta[[j]]))
    CI_MMbeta=rbind(CI_MMbeta,c(3,scenarios_CI[[pname]][i],lwbbetaMM,upbetaMM))
  
    # Check whether 250 and 10 are in MM-CI
    CP_MMbeta1=ifelse(lwbbetaMM[1]<250,1,0)*ifelse(250<upbetaMM[1],1,0)
    CP_MMbeta2=ifelse(lwbbetaMM[2]<10,1,0)*ifelse(10<upbetaMM[2],1,0)
    CICP_MMbeta=rbind(CICP_MMbeta,c(3,scenarios_CI[[pname]][i],CP_MMbeta1,CP_MMbeta2))  
  }
  # Determine CP for CI-beta1 and CI-beta2 with MM
  CP_MMbeta=rbind(CP_MMbeta,
                 c(3,scenarios_CI[[pname]][i],
                   sum(CICP_MMbeta[,3])/nrep,sum(CICP_MMbeta[,4])/nrep))
  
  # Combining all dataframes for beta in one.
  # This dataframe can be uses as input for classical boxplots or ggplot2
  CI_MLESMMbeta=rbind(CI_MLESMMbeta,
                      CI_MLEbeta,CI_Sbeta,CI_MMbeta)
  
  CICP_MLESMMbeta=rbind(CICP_MLESMMbeta,
                        CICP_MLEbeta,CICP_Sbeta,CICP_MMbeta)
  
  CP_MLESMMbeta=rbind(CP_MLESMMbeta,
                      CP_MLEbeta,CP_Sbeta,CP_MMbeta)
  
  #######################################################################
  # combining the extracted information for coverage probabiities 
  # for estimators for theta in a dataframe

  # dataframes for MLE theta
  CI_MLEtheta=NULL
  CICP_MLEtheta=NULL
  CP_MLEtheta=NULL
  
  for (j in 1:nrep){
    # Determine CI for theta1,theta2,theta3, and theta4 with MLE
    lwbthetaMLE=MLESMMTAUCOMPSTAUind$MLE$theta[j,]-qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$MLE$vartheta[[j]]))
    upbthetaMLE=MLESMMTAUCOMPSTAUind$MLE$theta[j,]+qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$MLE$vartheta[[j]]))
    CI_MLEtheta=rbind(CI_MLEtheta,c(1,scenarios_CI[[pname]][i],lwbthetaMLE,upbthetaMLE))
    
    # Check whether 790,-8.5,40, and 400 are in MLE-CI
    CP_MLEtheta1=ifelse(lwbthetaMLE[1]<790,1,0)*ifelse(790<upbthetaMLE[1],1,0)
    CP_MLEtheta2=ifelse(lwbthetaMLE[2]<-8.5,1,0)*ifelse(-8.5<upbthetaMLE[2],1,0)
    CP_MLEtheta3=ifelse(lwbthetaMLE[3]<40,1,0)*ifelse(40<upbthetaMLE[3],1,0)
    CP_MLEtheta4=ifelse(lwbthetaMLE[4]<400,1,0)*ifelse(400<upbthetaMLE[4],1,0)
    
    CICP_MLEtheta=rbind(CICP_MLEtheta,
                        c(1,scenarios_CI[[pname]][i],
                          CP_MLEtheta1,CP_MLEtheta2,CP_MLEtheta3,CP_MLEtheta4))  
  }
  # Determine CP for CI-theta1, CI-theta2, CI-theta3, and CI-theta4 with MLE
  CP_MLEtheta=rbind(CP_MLEtheta,
                    c(1,scenarios_CI[[pname]][i],
                      sum(CICP_MLEtheta[,3])/nrep,
                      sum(CICP_MLEtheta[,4])/nrep,
                      sum(CICP_MLEtheta[,5])/nrep,
                      sum(CICP_MLEtheta[,6])/nrep))
  
  # dataframes for S theta
  CI_Stheta=NULL
  CICP_Stheta=NULL
  CP_Stheta=NULL
  
  for (j in 1:nrep){
    # Determine CI for theta1,theta2,theta3, and theta4 with S
    lwbthetaS=MLESMMTAUCOMPSTAUind$S$theta[j,]-qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$S$vartheta[[j]]))
    upbthetaS=MLESMMTAUCOMPSTAUind$S$theta[j,]+qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$S$vartheta[[j]]))
    CI_Stheta=rbind(CI_Stheta,c(2,scenarios_CI[[pname]][i],lwbthetaS,upbthetaS))
    
    # Check whether 790,-8.5,40, and 400 are in S-CI
    CP_Stheta1=ifelse(lwbthetaS[1]<790,1,0)*ifelse(790<upbthetaS[1],1,0)
    CP_Stheta2=ifelse(lwbthetaS[2]<-8.5,1,0)*ifelse(-8.5<upbthetaS[2],1,0)
    CP_Stheta3=ifelse(lwbthetaS[3]<40,1,0)*ifelse(40<upbthetaS[3],1,0)
    CP_Stheta4=ifelse(lwbthetaS[4]<400,1,0)*ifelse(400<upbthetaS[4],1,0)
    
    CICP_Stheta=rbind(CICP_Stheta,
                      c(2,scenarios_CI[[pname]][i],
                        CP_Stheta1,CP_Stheta2,CP_Stheta3,CP_Stheta4))  
  }
  
  # Determine CP for CI-theta1, CI-theta2, CI-theta3, and CI-theta4 with S
  CP_Stheta=rbind(CP_Stheta,
                  c(2,scenarios_CI[[pname]][i],
                    sum(CICP_Stheta[,3])/nrep,
                    sum(CICP_Stheta[,4])/nrep,
                    sum(CICP_Stheta[,5])/nrep,
                    sum(CICP_Stheta[,6])/nrep))
  
  # dataframes for MM theta
  CI_MMtheta=NULL
  CICP_MMtheta=NULL
  CP_MMtheta=NULL
  
  for (j in 1:nrep){
    # Determine CI for theta1,theta2,theta3, and theta4 with MM
    lwbthetaMM=MLESMMTAUCOMPSTAUind$MM$theta[j,]-qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$MM$vartheta[[j]]))
    upbthetaMM=MLESMMTAUCOMPSTAUind$MM$theta[j,]+qnorm(0.975)*sqrt(diag(MLESMMTAUCOMPSTAUind$MM$vartheta[[j]]))
    CI_MMtheta=rbind(CI_MMtheta,c(3,scenarios_CI[[pname]][i],lwbthetaMM,upbthetaMM))
    
    # Check whether 790,-8.5,40, and 400 are in MM-CI
    CP_MMtheta1=ifelse(lwbthetaMM[1]<790,1,0)*ifelse(790<upbthetaMM[1],1,0)
    CP_MMtheta2=ifelse(lwbthetaMM[2]<-8.5,1,0)*ifelse(-8.5<upbthetaMM[2],1,0)
    CP_MMtheta3=ifelse(lwbthetaMM[3]<40,1,0)*ifelse(40<upbthetaMM[3],1,0)
    CP_MMtheta4=ifelse(lwbthetaMM[4]<400,1,0)*ifelse(400<upbthetaMM[4],1,0)
    
    CICP_MMtheta=rbind(CICP_MMtheta,
                       c(3,scenarios_CI[[pname]][i],
                         CP_MMtheta1,CP_MMtheta2,CP_MMtheta3,CP_MMtheta4))  
  }
  # Determine CP for CI-theta1, CI-theta2, CI-theta3, and CI-theta4 with MM
  CP_MMtheta=rbind(CP_MMtheta,
                   c(3,scenarios_CI[[pname]][i],
                     sum(CICP_MMtheta[,3])/nrep,
                     sum(CICP_MMtheta[,4])/nrep,
                     sum(CICP_MMtheta[,5])/nrep,
                     sum(CICP_MMtheta[,6])/nrep))
  
  # Combining all dataframes for theta in one.
  # This dataframe can be uses as input for classical boxplots or ggplot2
  CI_MLESMMtheta=rbind(CI_MLESMMtheta,
                       CI_MLEtheta,CI_Stheta)
  #,CI_MMtheta)
  
  CICP_MLESMMtheta=rbind(CICP_MLESMMtheta,
                         CICP_MLEtheta,CICP_Stheta)
  #,CICP_MMtheta)
  
  CP_MLESMMtheta=rbind(CP_MLESMMtheta,
                       CP_MLEtheta,CP_Stheta)
  #,CP_MMtheta)
} # END of scenarios_CI loop

####################################################################
# Setting the first two columns of CI_MLESMMbeta to factors
colnames(CI_MLESMMbeta)=c("Estimator",pname,
                          "lwbbeta1","lwbbeta2","upbbeta1","upbbeta2")
CI_MLESMMbeta=data.frame(CI_MLESMMbeta)
CI_MLESMMbeta[,1]=factor(CI_MLESMMbeta[,1],
                         levels=1:3,labels=c("MLE","S","MM"))
CI_MLESMMbeta[,2]=factor(CI_MLESMMbeta[,2],levels=levelsvec)

# Setting the first two columns of CICP_MLESMMbeta to factors
colnames(CICP_MLESMMbeta)=c("Estimator",pname,"CPindbeta1","CPindbeta2")
CICP_MLESMMbeta=data.frame(CICP_MLESMMbeta)
CICP_MLESMMbeta[,1]=factor(CICP_MLESMMbeta[,1],
                         levels=1:3,labels=c("MLE","S","MM"))
CICP_MLESMMbeta[,2]=factor(CICP_MLESMMbeta[,2],levels=levelsvec)

# Setting the first two columns of CP_MLESMMbeta to factors
colnames(CP_MLESMMbeta)=c("Estimator",pname,"CovPrbeta1","CovPrbeta2")
CP_MLESMMbeta=data.frame(CP_MLESMMbeta)
CP_MLESMMbeta[,1]=factor(CP_MLESMMbeta[,1],
                           levels=1:3,labels=c("MLE","S","MM"))
CP_MLESMMbeta[,2]=factor(CP_MLESMMbeta[,2],levels=levelsvec)

####################################################################
# Setting the first two columns of CI_MLESMMtheta to factors
colnames(CI_MLESMMtheta)=c("Estimator",pname,
                          "lwbtheta1","lwbtheta2","lwbtheta3","lwbtheta4",
                          "upbtheta1","upbtheta2","upbtheta3","upbtheta4")

CI_MLESMMtheta=data.frame(CI_MLESMMtheta)
CI_MLESMMtheta[,1]=factor(CI_MLESMMtheta[,1],
                         levels=1:2,labels=c("MLE","S"))
CI_MLESMMtheta[,2]=factor(CI_MLESMMtheta[,2],levels=levelsvec)

# Setting the first two columns of CICP_MLESMMtheta to factors
colnames(CICP_MLESMMtheta)=c("Estimator",pname,"CPindtheta1","CPindtheta2","CPindtheta3","CPindtheta4")
CICP_MLESMMtheta=data.frame(CICP_MLESMMtheta)
CICP_MLESMMtheta[,1]=factor(CICP_MLESMMtheta[,1],
                           levels=1:2,labels=c("MLE","S"))
CICP_MLESMMtheta[,2]=factor(CICP_MLESMMtheta[,2],levels=levelsvec)

# Setting the first two columns of CP_MLESMMtheta to factors
colnames(CP_MLESMMtheta)=c("Estimator",pname,"CovPrtheta1","CovPrtheta2","CovPrtheta3","CovPrtheta4")
CP_MLESMMtheta=data.frame(CP_MLESMMtheta)
CP_MLESMMtheta[,1]=factor(CP_MLESMMtheta[,1],
                         levels=1:2,labels=c("MLE","S"))
CP_MLESMMtheta[,2]=factor(CP_MLESMMtheta[,2],levels=levelsvec)

################################################################
# Preparing graphs

# Coverage probabilities using ggplot2
library(ggplot2)
library(gridExtra)

titlename=paste0("CI ICM: n=",unique(scenarios_CI[,2]),
                 ", k=",unique(scenarios_CI[,3]),
                 ", ",pname, ", perc=",perc)


# Coverage probabilities for beta
plotCPbeta1=ggplot(CP_MLESMMbeta,
                   aes(x=get(pname),y=CovPrbeta1,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Size of contamination")+
  ylab(quote(beta[1]))+
  ylim(c(0,1))+
  ggtitle(titlename)+
  geom_hline(yintercept = 0.95,lty=3)

plotCPbeta2=ggplot(CP_MLESMMbeta,aes(x=get(pname),y=CovPrbeta2,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Size of contamination")+
  ylab(quote(beta[2]))+
  ggtitle(titlename)+
  ylim(c(0,1))+
  geom_hline(yintercept = 0.95,lty=3)

grid.arrange(plotCPbeta1,plotCPbeta2,nrow=1)

# Coverage probabilities for theta
plotCPtheta1=ggplot(CP_MLESMMtheta,aes(x=get(pname),y=CovPrtheta1,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Size of contamination")+
  ylab(quote(sigma[gamma[1]]^2))+
  ggtitle(titlename)+
  ylim(c(0,1))+
  geom_hline(yintercept = 0.95,lty=3)

plotCPtheta2=ggplot(CP_MLESMMtheta,aes(x=get(pname),y=CovPrtheta2,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Size of contamination")+
  ylab(quote(sigma[gamma[12]]))+
  ggtitle(titlename)+
  ylim(c(0,1))+
  geom_hline(yintercept = 0.95,lty=3)

plotCPtheta3=ggplot(CP_MLESMMtheta,aes(x=get(pname),y=CovPrtheta3,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Size of contamination")+
  ylab(quote(sigma[gamma[2]]^2))+
  ggtitle(titlename)+
  ylim(c(0,1))+
  geom_hline(yintercept = 0.95,lty=3)

plotCPtheta4=ggplot(CP_MLESMMtheta,aes(x=get(pname),y=CovPrtheta4,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Size of contamination")+
  ylab(quote(sigma[epsilon]^2))+
  ggtitle(titlename)+
  ylim(c(0,1))+
  geom_hline(yintercept = 0.95,lty=3)

grid.arrange(plotCPtheta1,plotCPtheta2,plotCPtheta3,plotCPtheta4,nrow=2)
