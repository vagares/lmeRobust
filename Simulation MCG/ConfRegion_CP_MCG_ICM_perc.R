# This script uses the contamination scenarios and results of
# the corresponding simulation performed by
# Simulation_setting_model_MCG.R

# The script extracts the information from the .RData files
# corresponding to the different setups in dataframe 
# scenarios_CP and produces graph of the coverage probabilities

# Sets scenarios_CP equal to scenarios_CP
# scenarios_CP is used below to extract the .RData files

source("scenarios_MCG_simulation.R")

scenarios_CP=scenarios[(scenarios$mec==-80)&(scenarios$CCMind==0),]

if (max(scenarios_CP$pe)!=0){
  pname="pe";
  size=unique(scenarios_CP$mec);
  levelsvec=pevec}
if (max(scenarios_CP$pb)!=0){
  pname="pb";
  size=unique(scenarios_CP$mbc2);
  levelsvec=pbvec}
if (max(scenarios_CP$px)!=0){
  pname="px";
  size=unique(scenarios_CP$alphac);
  levelsvec=pxvec}

#######################################################
# Extracting the information for beta from .Rdata files

# Dataframes for betahat and thetahat
CPCR_MLESMMbeta=NULL
CPCR_MLESMMtheta=NULL

for (i in 1:nrow(scenarios_CP)){
  nrep=scenarios_CP[i,1]
  nsample=scenarios_CP[i,2]
  ksample=scenarios_CP[i,3]
  pesample=scenarios_CP[i,4]
  pbsample=scenarios_CP[i,5]
  pxsample=scenarios_CP[i,6]
  mecsample=scenarios_CP[i,7]
  mbc2sample=scenarios_CP[i,8]
  alphacsample=scenarios_CP[i,9]
  rcsample=scenarios_CP[i,10]
  cTAUindsample=as.logical(scenarios_CP[i,11])
  CCMindsample=as.logical(scenarios_CP[i,12])
  
  if (CCMindsample==TRUE){stop("Contains CCM scenarios_CP")}
  
  # Setting the filename depending on yes/no cTAU and yes/no CCM
  if (cTAUindsample==FALSE){
    if (CCMindsample==FALSE){
      flnameEst="MLESMM_ICM"}else{
        flnameEst="MLESMM_CCM"}
  }else{
    if (CCMindsample==FALSE){
      flnameEst="MLESMMcTAU_ICM"}else{
        flnameEst="MLESMMcTAU_CCM"}
  }
  
  if ((scenarios_CP$pe[i]==0)&
      (scenarios_CP$pb[i]==0)&
      (scenarios_CP$px[i]==0)){
    flname=paste0("./Results_Uncontaminated/",flnameEst,"_",
                  "nrep=",scenarios_CP$nrep[i],"_",
                  "n=",scenarios_CP$n[i],"_",
                  "k=",scenarios_CP$k[i],"_",
                  "pe=",scenarios_CP$pe[i],"_",
                  "pb=",scenarios_CP$pb[i],"_",
                  "px=",scenarios_CP$px[i],"_",
                  "mec=",scenarios_CP$mec[i],"_",
                  "mbc2=",scenarios_CP$mbc2[i],"_",
                  "alphac=",scenarios_CP$alphac[i],"_",
                  "rc=",scenarios_CP$rc[i],".RData")}
  
  if (scenarios_CP$pe[i]>0){
    flname=paste0("./Results_Epsilon_contamination/",flnameEst,"_",
                  "nrep=",scenarios_CP$nrep[i],"_",
                  "n=",scenarios_CP$n[i],"_",
                  "k=",scenarios_CP$k[i],"_",
                  "pe=",scenarios_CP$pe[i],"_",
                  "pb=",scenarios_CP$pb[i],"_",
                  "px=",scenarios_CP$px[i],"_",
                  "mec=",scenarios_CP$mec[i],"_",
                  "mbc2=",scenarios_CP$mbc2[i],"_",
                  "alphac=",scenarios_CP$alphac[i],"_",
                  "rc=",scenarios_CP$rc[i],".RData")}
  
  if (scenarios_CP$pb[i]>0){
    flname=paste0("./Results_Random_Effect_contamination/",flnameEst,"_",
                  "nrep=",scenarios_CP$nrep[i],"_",
                  "n=",scenarios_CP$n[i],"_",
                  "k=",scenarios_CP$k[i],"_",
                  "pe=",scenarios_CP$pe[i],"_",
                  "pb=",scenarios_CP$pb[i],"_",
                  "px=",scenarios_CP$px[i],"_",
                  "mec=",scenarios_CP$mec[i],"_",
                  "mbc2=",scenarios_CP$mbc2[i],"_",
                  "alphac=",scenarios_CP$alphac[i],"_",
                  "rc=",scenarios_CP$rc[i],".RData")}
  
  if (scenarios_CP$px[i]>0){
    flname=paste0("./Results_X_contamination/",flnameEst,"_",
                  "nrep=",scenarios_CP$nrep[i],"_",
                  "n=",scenarios_CP$n[i],"_",
                  "k=",scenarios_CP$k[i],"_",
                  "pe=",scenarios_CP$pe[i],"_",
                  "pb=",scenarios_CP$pb[i],"_",
                  "px=",scenarios_CP$px[i],"_",
                  "mec=",scenarios_CP$mec[i],"_",
                  "mbc2=",scenarios_CP$mbc2[i],"_",
                  "alphac=",scenarios_CP$alphac[i],"_",
                  "rc=",scenarios_CP$rc[i],".RData")}
  
  load(flname)
  
  #######################################################################
  # combining the extracted information for coverage probabiities 
  # for estimators for beta in a dataframe
  
  beta0=c(250,10)
  theta0=c(790,-8.5,40,400)
  X2beta=qchisq(0.95,length(beta0)) 
  X2theta=qchisq(0.95,length(theta0)) 
  
  ########################################
  # dataframes for MLE beta
  CPCR_MLEbeta=NULL
  CPCR_MLEbeta_count=0
  CPCR_MLEtheta=NULL
  CPCR_MLEtheta_count=0
  
  for (j in 1:nrep){
    CP_MLEbetaMD2=ifelse(mahalanobis(beta0,center=MLESMMcTAUind$MLE$beta[j,],
                                     cov=MLESMMcTAUind$MLE$varbeta[[j]])<X2beta,1,0)
    CPCR_MLEbeta_count=CPCR_MLEbeta_count+CP_MLEbetaMD2  
    
    CP_MLEthetaMD2=ifelse(mahalanobis(theta0,center=MLESMMcTAUind$MLE$theta[j,],
                                     cov=MLESMMcTAUind$MLE$vartheta[[j]])<X2theta,1,0)
    CPCR_MLEtheta_count=CPCR_MLEtheta_count+CP_MLEthetaMD2  
  }
  # Determine CP for CR-beta and CR-theta with MLE
  CPCR_MLEbeta=rbind(CPCR_MLEbeta,c(1,scenarios_CP[[pname]][i],CPCR_MLEbeta_count/nrep))
  CPCR_MLEtheta=rbind(CPCR_MLEtheta,c(1,scenarios_CP[[pname]][i],CPCR_MLEtheta_count/nrep))
  
  ########################################
  # dataframes for S beta
  CPCR_Sbeta=NULL
  CPCR_Sbeta_count=0
  CPCR_Stheta=NULL
  CPCR_Stheta_count=0
  
  for (j in 1:nrep){
    CP_SbetaMD2=ifelse(mahalanobis(beta0,center=MLESMMcTAUind$S$beta[j,],
                                     cov=MLESMMcTAUind$S$varbeta[[j]])<X2beta,1,0)
    CPCR_Sbeta_count=CPCR_Sbeta_count+CP_SbetaMD2  

    CP_SthetaMD2=ifelse(mahalanobis(theta0,center=MLESMMcTAUind$S$theta[j,],
                                      cov=MLESMMcTAUind$S$vartheta[[j]])<X2theta,1,0)
    CPCR_Stheta_count=CPCR_Stheta_count+CP_SthetaMD2  
  }
  # Determine CP for CR-beta and CR-theta with S
  CPCR_Sbeta=rbind(CPCR_Sbeta,c(2,scenarios_CP[[pname]][i],CPCR_Sbeta_count/nrep))
  CPCR_Stheta=rbind(CPCR_Stheta,c(2,scenarios_CP[[pname]][i],CPCR_Stheta_count/nrep))
  
  ########################################
  CPCR_MMbeta=NULL
  CPCR_MMbeta_count=0

  for (j in 1:nrep){
    CP_MMbetaMD2=ifelse(mahalanobis(beta0,center=MLESMMcTAUind$MM$beta[j,],
                                     cov=MLESMMcTAUind$MM$varbeta[[j]])<X2beta,1,0)
    CPCR_MMbeta_count=CPCR_MMbeta_count+CP_MMbetaMD2  
  }
  # Determine CP for CR-beta with MM
  CPCR_MMbeta=rbind(CPCR_MMbeta,c(3,scenarios_CP[[pname]][i],CPCR_MMbeta_count/nrep))

  
  # Combining all dataframes for beta and theta in one.
  # This dataframe can be uses as input for classical boxplots or ggplot2
  CPCR_MLESMMbeta=rbind(CPCR_MLESMMbeta,
                        CPCR_MLEbeta,CPCR_Sbeta,CPCR_MMbeta)
  
  CPCR_MLESMMtheta=rbind(CPCR_MLESMMtheta,
                        CPCR_MLEtheta,CPCR_Stheta)
} # END of scenarios loop

####################################################################
# Setting the first two columns of CPCR_MLESMMbeta to factors
colnames(CPCR_MLESMMbeta)=c("Estimator",pname,"CPbeta")
CPCR_MLESMMbeta=data.frame(CPCR_MLESMMbeta)
CPCR_MLESMMbeta[,1]=factor(CPCR_MLESMMbeta[,1],
                         levels=1:3,labels=c("MLE","S","MM"))
CPCR_MLESMMbeta[,2]=factor(CPCR_MLESMMbeta[,2],levels=levelsvec)

# Setting the first two columns of CPCR_MLESMMtheta to factors
colnames(CPCR_MLESMMtheta)=c("Estimator",pname,"CPtheta")
CPCR_MLESMMtheta=data.frame(CPCR_MLESMMtheta)
CPCR_MLESMMtheta[,1]=factor(CPCR_MLESMMtheta[,1],
                           levels=1:2,labels=c("MLE","S"))
CPCR_MLESMMtheta[,2]=factor(CPCR_MLESMMtheta[,2],levels=levelsvec)

################################################################
# Preparing graphs

# Coverage probabilities using ggplot2
library(ggplot2)
library(gridExtra)
titlename=paste0("CP ICM: n=",unique(scenarios_CP[,2]),
                 ", k=",unique(scenarios_CP[,3]),
                 ", ",pname, ", size=",size)


# Coverage probabilities for beta
plotCPCRbeta=ggplot(CPCR_MLESMMbeta,aes(x=get(pname),y=CPbeta,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Size of contamination")+
  ylab(quote(beta))+
  ylim(c(0,1))+
  ggtitle(titlename)+
  geom_hline(yintercept = 0.95,lty=3)

# Coverage probabilities for theta
plotCPCRtheta=ggplot(CPCR_MLESMMtheta,aes(x=get(pname),y=CPtheta,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Size of contamination")+
  ylab(quote(theta))+
  ylim(c(0,1))+
  ggtitle(titlename)+
  geom_hline(yintercept = 0.95,lty=3)

grid.arrange(plotCPCRbeta,plotCPCRtheta,nrow=1)