# This script uses the contamination scenarios and results of
# the corresponding simulation performed by
# Simulation_setting_model_MCG.R

# The script extracts the information from the .RData files
# corresponding to the different setups in dataframe scenarios
# and produces coverage probabilities for the confidence REGIONS
# of the vectors of estimators


#######################################################
# Extracting the information for beta from .Rdata files

# Dataframes for betahat and thetahat
CPCR_MLESMMbeta=NULL
CPCR_MLESMMtheta=NULL

for (i in 1:nrow(scenarios)){
if (scenarios$pe[i]>0){
    flname=paste0("./Results_Epsilon_contamination/","MLESMM","_",
                "nrep=",scenarios$nrep[i],"_",
                "n=",scenarios$n[i],"_",
                "k=",scenarios$k[i],"_",
                "pe=",scenarios$pe[i],"_",
                "pb=",scenarios$pb[i],"_",
                "px=",scenarios$px[i],"_",
                "mec=",scenarios$mec[i],"_",
                "mbc2=",scenarios$mbc2[i],"_",
                "alphac=",scenarios$alphac[i],".RData")}else{
                  flname=paste0("./Results_Uncontaminated/","MLESMM","_",
                                "nrep=",scenarios$nrep[i],"_",
                                "n=",scenarios$n[i],"_",
                                "k=",scenarios$k[i],"_",
                                "pe=",scenarios$pe[i],"_",
                                "pb=",scenarios$pb[i],"_",
                                "px=",scenarios$px[i],"_",
                                "mec=",scenarios$mec[i],"_",
                                "mbc2=",scenarios$mbc2[i],"_",
                                "alphac=",scenarios$alphac[i],".RData")}

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
    CP_MLEbetaMD2=ifelse(mahalanobis(beta0,center=MLESMM$MLE$beta[j,],
                                     cov=MLESMM$MLE$varbeta[[j]])<X2beta,1,0)
    CPCR_MLEbeta_count=CPCR_MLEbeta_count+CP_MLEbetaMD2  
    
    CP_MLEthetaMD2=ifelse(mahalanobis(theta0,center=MLESMM$MLE$theta[j,],
                                     cov=MLESMM$MLE$vartheta[[j]])<X2theta,1,0)
    CPCR_MLEtheta_count=CPCR_MLEtheta_count+CP_MLEthetaMD2  
  }
  # Determine CP for CR-beta and CR-theta with MLE
  CPCR_MLEbeta=rbind(CPCR_MLEbeta,c(1,scenarios$pe[i],CPCR_MLEbeta_count/nrep))
  CPCR_MLEtheta=rbind(CPCR_MLEtheta,c(1,scenarios$pe[i],CPCR_MLEtheta_count/nrep))
  
  ########################################
  # dataframes for S beta
  CPCR_Sbeta=NULL
  CPCR_Sbeta_count=0
  CPCR_Stheta=NULL
  CPCR_Stheta_count=0
  
  for (j in 1:nrep){
    CP_SbetaMD2=ifelse(mahalanobis(beta0,center=MLESMM$S$beta[j,],
                                     cov=MLESMM$S$varbeta[[j]])<X2beta,1,0)
    CPCR_Sbeta_count=CPCR_Sbeta_count+CP_SbetaMD2  

    CP_SthetaMD2=ifelse(mahalanobis(theta0,center=MLESMM$S$theta[j,],
                                      cov=MLESMM$S$vartheta[[j]])<X2theta,1,0)
    CPCR_Stheta_count=CPCR_Stheta_count+CP_SthetaMD2  
  }
  # Determine CP for CR-beta and CR-theta with S
  CPCR_Sbeta=rbind(CPCR_Sbeta,c(2,scenarios$pe[i],CPCR_Sbeta_count/nrep))
  CPCR_Stheta=rbind(CPCR_Stheta,c(2,scenarios$pe[i],CPCR_Stheta_count/nrep))
  
  ########################################
  CPCR_MMbeta=NULL
  CPCR_MMbeta_count=0

  for (j in 1:nrep){
    CP_MMbetaMD2=ifelse(mahalanobis(beta0,center=MLESMM$MM$beta[j,],
                                     cov=MLESMM$MM$varbeta[[j]])<X2beta,1,0)
    CPCR_MMbeta_count=CPCR_MMbeta_count+CP_MMbetaMD2  
  }
  # Determine CP for CR-beta with MM
  CPCR_MMbeta=rbind(CPCR_MMbeta,c(3,scenarios$pe[i],CPCR_MMbeta_count/nrep))

  
  # Combining all dataframes for beta and theta in one.
  # This dataframe can be uses as input for classical boxplots or ggplot2
  CPCR_MLESMMbeta=rbind(CPCR_MLESMMbeta,
                        CPCR_MLEbeta,CPCR_Sbeta,CPCR_MMbeta)
  
  CPCR_MLESMMtheta=rbind(CPCR_MLESMMtheta,
                        CPCR_MLEtheta,CPCR_Stheta)
} # END of scenatios loop

####################################################################
# Setting the first two columns of CPCR_MLESMMbeta to factors
colnames(CPCR_MLESMMbeta)=c("Estimator","pe","CPbeta")
CPCR_MLESMMbeta=data.frame(CPCR_MLESMMbeta)
CPCR_MLESMMbeta[,1]=factor(CPCR_MLESMMbeta[,1],
                         levels=1:3,labels=c("MLE","S","MM"))
CPCR_MLESMMbeta[,2]=factor(CPCR_MLESMMbeta[,2],levels=pevec,
                             labels = c("00","01","02","03","04","05",
                                        "06","07","08","09","10"))

# Setting the first two columns of CPCR_MLESMMtheta to factors
colnames(CPCR_MLESMMtheta)=c("Estimator","pe","CPtheta")
CPCR_MLESMMtheta=data.frame(CPCR_MLESMMtheta)
CPCR_MLESMMtheta[,1]=factor(CPCR_MLESMMtheta[,1],
                           levels=1:2,labels=c("MLE","S"))
CPCR_MLESMMtheta[,2]=factor(CPCR_MLESMMtheta[,2],levels=pevec,
                           labels = c("00","01","02","03","04","05",
                                      "06","07","08","09","10"))

################################################################
# Preparing graphs

# Coverage probabilities using ggplot2
library(ggplot2)
library(gridExtra)

# Coverage probabilities for beta
plotCPCRbeta=ggplot(CPCR_MLESMMbeta,aes(x=pe,y=CPbeta,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Percentage of contamination")+
  ylab("Coverage probability of confidence region")+
  ylim(c(0,1))+
  ggtitle("beta")+
  geom_hline(yintercept = 0.95,lty=3)

# Coverage probabilities for theta
plotCPCRtheta=ggplot(CPCR_MLESMMtheta,aes(x=pe,y=CPtheta,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Percentage of contamination")+
  ylab("Coverage probability of confidence region")+
  ylim(c(0,1))+
  ggtitle("theta")+
  geom_hline(yintercept = 0.95,lty=3)

grid.arrange(plotCPCRbeta,plotCPCRtheta,nrow=1)