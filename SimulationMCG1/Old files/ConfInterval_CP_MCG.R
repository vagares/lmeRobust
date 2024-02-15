# This script uses the contamination scenarios and results of
# the corresponding simulation performed by
# Simulation_setting_model_MCG.R

# The script extracts the information from the .RData files
# corresponding to the different setups in dataframe 
# scenarios_CI and produces graph of the coverage probabilities

#######################################################################
# creating dataframe containing different contamination schemes as rows
# SCENARIO 1: - -80 fixed
#             - pe 0, 0.01, ... , 0.10
#             - 250 repetitions

pevec=seq(0,0.10,by=0.05)
nrep=50
scenarios_CI=NULL
for (i in 1:length(pevec)){
  scenarios_CI=rbind(scenarios_CI,c(nrep,200,4,pevec[i],0,0,-80,0,1))
}
colnames(scenarios_CI)=c("nrep","n","k","pe","pb","px","mec","mbc2","alphac")

# Sets scenarios equal to scenarios_CI
# scenarios is used below to extract the .RData files
scenarios=data.frame(scenarios_CI)

#######################################################
# Extracting the information for beta from .Rdata files

# Dataframes for betahat and thetahat
CI_MLESMMbeta=NULL
CICP_MLESMMbeta=NULL
CP_MLESMMbeta=NULL

CI_MLESMMtheta=NULL
CICP_MLESMMtheta=NULL
CP_MLESMMtheta=NULL

for (i in 1:nrow(scenarios)){
if (scenarios$pe[i]>0){
  flname=paste0("./Results_Epsilon_contamination/","MLESMM","_",
#  flname=paste0("./Results_Epsilon_contamination/EpsilonCCM/","MLESMM","_",
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

  # dataframes for MLE beta
  CI_MLEbeta=NULL
  CICP_MLEbeta=NULL
  CP_MLEbeta=NULL

  for (j in 1:nrep){
    # Determine CI for beta1 and beta2 with MLE
    lwbbetaMLE=MLESMM$MLE$beta[j,]-qnorm(0.975)*sqrt(diag(MLESMM$MLE$varbeta[[j]]))
    upbetaMLE=MLESMM$MLE$beta[j,]+qnorm(0.975)*sqrt(diag(MLESMM$MLE$varbeta[[j]]))
    CI_MLEbeta=rbind(CI_MLEbeta,c(1,scenarios$pe[i],lwbbetaMLE,upbetaMLE))

    # Check whether 250 and 10 are in MLE-CI
    CP_MLEbeta1=ifelse(lwbbetaMLE[1]<250,1,0)*ifelse(250<upbetaMLE[1],1,0)
    CP_MLEbeta2=ifelse(lwbbetaMLE[2]<10,1,0)*ifelse(10<upbetaMLE[2],1,0)
    CICP_MLEbeta=rbind(CICP_MLEbeta,c(1,scenarios$pe[i],CP_MLEbeta1,CP_MLEbeta2))  
  }
  # Determine CP for CI-beta1 and CI-beta2 with MLE
  CP_MLEbeta=rbind(CP_MLEbeta,
                   c(1,scenarios$pe[i],
                     sum(CICP_MLEbeta[,3])/nrep,sum(CICP_MLEbeta[,4])/nrep))
  
  # dataframes for S beta
  CI_Sbeta=NULL
  CICP_Sbeta=NULL
  CP_Sbeta=NULL
  
  for (j in 1:nrep){
    # Determine CI for beta1 and beta2 with S
    lwbbetaS=MLESMM$S$beta[j,]-qnorm(0.975)*sqrt(diag(MLESMM$S$varbeta[[j]]))
    upbetaS=MLESMM$S$beta[j,]+qnorm(0.975)*sqrt(diag(MLESMM$S$varbeta[[j]]))
    CI_Sbeta=rbind(CI_Sbeta,c(2,scenarios$pe[i],lwbbetaS,upbetaS))
  
    # Check whether 250 and 10 are in S-CI
    CP_Sbeta1=ifelse(lwbbetaS[1]<250,1,0)*ifelse(250<upbetaS[1],1,0)
    CP_Sbeta2=ifelse(lwbbetaS[2]<10,1,0)*ifelse(10<upbetaS[2],1,0)
    CICP_Sbeta=rbind(CICP_Sbeta,c(2,scenarios$pe[i],CP_Sbeta1,CP_Sbeta2))  
  }
  # Determine CP for CI-beta1 and CI-beta2 with S
  CP_Sbeta=rbind(CP_Sbeta,
                   c(2,scenarios$pe[i],
                     sum(CICP_Sbeta[,3])/nrep,sum(CICP_Sbeta[,4])/nrep))
  
  # dataframes for MM beta
  CI_MMbeta=NULL
  CICP_MMbeta=NULL
  CP_MMbeta=NULL
  
  for (j in 1:nrep){
    # Determine CI for beta1 and beta2 with MM
    lwbbetaMM=MLESMM$MM$beta[j,]-qnorm(0.975)*sqrt(diag(MLESMM$MM$varbeta[[j]]))
    upbetaMM=MLESMM$MM$beta[j,]+qnorm(0.975)*sqrt(diag(MLESMM$MM$varbeta[[j]]))
    CI_MMbeta=rbind(CI_MMbeta,c(3,scenarios$pe[i],lwbbetaMM,upbetaMM))
  
    # Check whether 250 and 10 are in MM-CI
    CP_MMbeta1=ifelse(lwbbetaMM[1]<250,1,0)*ifelse(250<upbetaMM[1],1,0)
    CP_MMbeta2=ifelse(lwbbetaMM[2]<10,1,0)*ifelse(10<upbetaMM[2],1,0)
    CICP_MMbeta=rbind(CICP_MMbeta,c(3,scenarios$pe[i],CP_MMbeta1,CP_MMbeta2))  
  }
  # Determine CP for CI-beta1 and CI-beta2 with MM
  CP_MMbeta=rbind(CP_MMbeta,
                 c(3,scenarios$pe[i],
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
    lwbthetaMLE=MLESMM$MLE$theta[j,]-qnorm(0.975)*sqrt(diag(MLESMM$MLE$vartheta[[j]]))
    upbthetaMLE=MLESMM$MLE$theta[j,]+qnorm(0.975)*sqrt(diag(MLESMM$MLE$vartheta[[j]]))
    CI_MLEtheta=rbind(CI_MLEtheta,c(1,scenarios$pe[i],lwbthetaMLE,upbthetaMLE))
    
    # Check whether 790,-8.5,40, and 400 are in MLE-CI
    CP_MLEtheta1=ifelse(lwbthetaMLE[1]<790,1,0)*ifelse(790<upbthetaMLE[1],1,0)
    CP_MLEtheta2=ifelse(lwbthetaMLE[2]<-8.5,1,0)*ifelse(-8.5<upbthetaMLE[2],1,0)
    CP_MLEtheta3=ifelse(lwbthetaMLE[3]<40,1,0)*ifelse(40<upbthetaMLE[3],1,0)
    CP_MLEtheta4=ifelse(lwbthetaMLE[4]<400,1,0)*ifelse(400<upbthetaMLE[4],1,0)
    
    CICP_MLEtheta=rbind(CICP_MLEtheta,
                        c(1,scenarios$pe[i],
                          CP_MLEtheta1,CP_MLEtheta2,CP_MLEtheta3,CP_MLEtheta4))  
  }
  # Determine CP for CI-theta1, CI-theta2, CI-theta3, and CI-theta4 with MLE
  CP_MLEtheta=rbind(CP_MLEtheta,
                    c(1,scenarios$pe[i],
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
    lwbthetaS=MLESMM$S$theta[j,]-qnorm(0.975)*sqrt(diag(MLESMM$S$vartheta[[j]]))
    upbthetaS=MLESMM$S$theta[j,]+qnorm(0.975)*sqrt(diag(MLESMM$S$vartheta[[j]]))
    CI_Stheta=rbind(CI_Stheta,c(2,scenarios$pe[i],lwbthetaS,upbthetaS))
    
    # Check whether 790,-8.5,40, and 400 are in S-CI
    CP_Stheta1=ifelse(lwbthetaS[1]<790,1,0)*ifelse(790<upbthetaS[1],1,0)
    CP_Stheta2=ifelse(lwbthetaS[2]<-8.5,1,0)*ifelse(-8.5<upbthetaS[2],1,0)
    CP_Stheta3=ifelse(lwbthetaS[3]<40,1,0)*ifelse(40<upbthetaS[3],1,0)
    CP_Stheta4=ifelse(lwbthetaS[4]<400,1,0)*ifelse(400<upbthetaS[4],1,0)
    
    CICP_Stheta=rbind(CICP_Stheta,
                      c(2,scenarios$pe[i],
                        CP_Stheta1,CP_Stheta2,CP_Stheta3,CP_Stheta4))  
  }
  
  # Determine CP for CI-theta1, CI-theta2, CI-theta3, and CI-theta4 with S
  CP_Stheta=rbind(CP_Stheta,
                  c(2,scenarios$pe[i],
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
    lwbthetaMM=MLESMM$MM$theta[j,]-qnorm(0.975)*sqrt(diag(MLESMM$MM$vartheta[[j]]))
    upbthetaMM=MLESMM$MM$theta[j,]+qnorm(0.975)*sqrt(diag(MLESMM$MM$vartheta[[j]]))
    CI_MMtheta=rbind(CI_MMtheta,c(3,scenarios$pe[i],lwbthetaMM,upbthetaMM))
    
    # Check whether 790,-8.5,40, and 400 are in MM-CI
    CP_MMtheta1=ifelse(lwbthetaMM[1]<790,1,0)*ifelse(790<upbthetaMM[1],1,0)
    CP_MMtheta2=ifelse(lwbthetaMM[2]<-8.5,1,0)*ifelse(-8.5<upbthetaMM[2],1,0)
    CP_MMtheta3=ifelse(lwbthetaMM[3]<40,1,0)*ifelse(40<upbthetaMM[3],1,0)
    CP_MMtheta4=ifelse(lwbthetaMM[4]<400,1,0)*ifelse(400<upbthetaMM[4],1,0)
    
    CICP_MMtheta=rbind(CICP_MMtheta,
                       c(3,scenarios$pe[i],
                         CP_MMtheta1,CP_MMtheta2,CP_MMtheta3,CP_MMtheta4))  
  }
  # Determine CP for CI-theta1, CI-theta2, CI-theta3, and CI-theta4 with MM
  CP_MMtheta=rbind(CP_MMtheta,
                   c(3,scenarios$pe[i],
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
} # END of scenarios loop

####################################################################
# Setting the first two columns of CI_MLESMMbeta to factors
colnames(CI_MLESMMbeta)=c("Estimator","pe",
                          "lwbbeta1","lwbbeta2","upbbeta1","upbbeta2")
CI_MLESMMbeta=data.frame(CI_MLESMMbeta)
CI_MLESMMbeta[,1]=factor(CI_MLESMMbeta[,1],
                         levels=1:3,labels=c("MLE","S","MM"))
CI_MLESMMbeta[,2]=factor(CI_MLESMMbeta[,2],levels=pevec,
                         labels = as.character(pevec*100))

# Setting the first two columns of CICP_MLESMMbeta to factors
colnames(CICP_MLESMMbeta)=c("Estimator","pe","CPindbeta1","CPindbeta2")
CICP_MLESMMbeta=data.frame(CICP_MLESMMbeta)
CICP_MLESMMbeta[,1]=factor(CICP_MLESMMbeta[,1],
                         levels=1:3,labels=c("MLE","S","MM"))
CICP_MLESMMbeta[,2]=factor(CICP_MLESMMbeta[,2],levels=pevec,
                           labels = as.character(pevec*100))

# Setting the first two columns of CP_MLESMMbeta to factors
colnames(CP_MLESMMbeta)=c("Estimator","pe","CovPrbeta1","CovPrbeta2")
CP_MLESMMbeta=data.frame(CP_MLESMMbeta)
CP_MLESMMbeta[,1]=factor(CP_MLESMMbeta[,1],
                           levels=1:3,labels=c("MLE","S","MM"))
CP_MLESMMbeta[,2]=factor(CP_MLESMMbeta[,2],levels=pevec,
                         labels = as.character(pevec*100))

nrow(CI_MLESMMbeta)
nrow(CICP_MLESMMbeta)
nrow(CP_MLESMMbeta)

####################################################################
# Setting the first two columns of CI_MLESMMtheta to factors
colnames(CI_MLESMMtheta)=c("Estimator","pe",
                          "lwbtheta1","lwbtheta2","lwbtheta3","lwbtheta4",
                          "upbtheta1","upbtheta2","upbtheta3","upbtheta4")

CI_MLESMMtheta=data.frame(CI_MLESMMtheta)
CI_MLESMMtheta[,1]=factor(CI_MLESMMtheta[,1],
                         levels=1:2,labels=c("MLE","S"))
CI_MLESMMtheta[,2]=factor(CI_MLESMMtheta[,2],levels=pevec,
                          labels = as.character(pevec*100))

# Setting the first two columns of CICP_MLESMMtheta to factors
colnames(CICP_MLESMMtheta)=c("Estimator","pe","CPindtheta1","CPindtheta2","CPindtheta3","CPindtheta4")
CICP_MLESMMtheta=data.frame(CICP_MLESMMtheta)
CICP_MLESMMtheta[,1]=factor(CICP_MLESMMtheta[,1],
                           levels=1:2,labels=c("MLE","S"))
CICP_MLESMMtheta[,2]=factor(CICP_MLESMMtheta[,2],levels=pevec,
                            labels = as.character(pevec*100))

# Setting the first two columns of CP_MLESMMtheta to factors
colnames(CP_MLESMMtheta)=c("Estimator","pe","CovPrtheta1","CovPrtheta2","CovPrtheta3","CovPrtheta4")
CP_MLESMMtheta=data.frame(CP_MLESMMtheta)
CP_MLESMMtheta[,1]=factor(CP_MLESMMtheta[,1],
                         levels=1:2,labels=c("MLE","S"))
CP_MLESMMtheta[,2]=factor(CP_MLESMMtheta[,2],levels=pevec,
                          labels = as.character(pevec*100))

nrow(CI_MLESMMtheta)
nrow(CICP_MLESMMtheta)
nrow(CP_MLESMMtheta)

################################################################
# Preparing graphs

# Coverage probabilities using ggplot2
library(ggplot2)
library(gridExtra)

# Coverage probabilities for beta
plotCPbeta1=ggplot(CP_MLESMMbeta,aes(x=pe,y=CovPrbeta1,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Percentage of contamination")+
  ylab("Coverage probability")+
  ylim(c(0,1))+
  ggtitle("beta1")+
  geom_hline(yintercept = 0.95,lty=3)

plotCPbeta2=ggplot(CP_MLESMMbeta,aes(x=pe,y=CovPrbeta2,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Percentage of contamination")+
  ylab("Coverage probability")+
  ylim(c(0,1))+
  ggtitle("beta2")+
  geom_hline(yintercept = 0.95,lty=3)

grid.arrange(plotCPbeta1,plotCPbeta2,nrow=1)

# Coverage probabilities for theta
plotCPtheta1=ggplot(CP_MLESMMtheta,aes(x=pe,y=CovPrtheta1,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Percentage of contamination")+
  ylab("Coverage probability")+
  ylim(c(0,1))+
  ggtitle("theta1")+
  geom_hline(yintercept = 0.95,lty=3)

plotCPtheta2=ggplot(CP_MLESMMtheta,aes(x=pe,y=CovPrtheta2,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Percentage of contamination")+
  ylab("Coverage probability")+
  ylim(c(0,1))+
  ggtitle("theta2")+
  geom_hline(yintercept = 0.95,lty=3)

plotCPtheta3=ggplot(CP_MLESMMtheta,aes(x=pe,y=CovPrtheta3,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Percentage of contamination")+
  ylab("Coverage probability")+
  ylim(c(0,1))+
  ggtitle("theta3")+
  geom_hline(yintercept = 0.95,lty=3)

plotCPtheta4=ggplot(CP_MLESMMtheta,aes(x=pe,y=CovPrtheta4,group=Estimator,color=Estimator))+
  geom_line()+
  xlab("Percentage of contamination")+
  ylab("Coverage probability")+
  ylim(c(0,1))+
  ggtitle("theta4")+
  geom_hline(yintercept = 0.95,lty=3)

grid.arrange(plotCPtheta1,plotCPtheta2,plotCPtheta3,plotCPtheta4,nrow=2)
