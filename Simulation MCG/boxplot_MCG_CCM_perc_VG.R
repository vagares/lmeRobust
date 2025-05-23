# This script uses the contamination scenarios and results of
# the corresponding simulation performed by
# Simulation_setting_model_MCG.R

# The script extracts the information from the .RData files
# corresponding to the different setups in dataframe 
# scenarios_boxplot and produces boxplots of the estimators

source("scenarios_MCG_simulation_VG.R")

scenarios_boxplot=scenarios[(scenarios$mec==-80)&(scenarios$CCMind==1),]
#scenarios_boxplot=scenarios[(scenarios$mbc2==-25)&(scenarios$CCMind==1),]
#scenarios_boxplot=scenarios[(scenarios$alphac==5)&(scenarios$CCMind==1),]
scenarios_boxplot=scenarios[(scenarios$mux==0.5)&(scenarios$CCMind==1),]
#scenarios_boxplot=scenarios[(scenarios$mec==0)&(scenarios$CCMind==1),]

if (max(scenarios_boxplot$pe)!=0){
  scenarios_boxplot = scenarios_boxplot[scenarios_boxplot$pb==0,]
  scenarios_boxplot = scenarios_boxplot[scenarios_boxplot$px==0,]
  pname="pe";
  noname="noe";
  noiname="noei";
  size=unique(scenarios_boxplot$mec);
  levelsvec=pevec}
if (max(scenarios_boxplot$pb)!=0){
  scenarios_boxplot = scenarios_boxplot[scenarios_boxplot$pe==0,]
  scenarios_boxplot = scenarios_boxplot[scenarios_boxplot$px==0,]
  pname="pb";
  noname="nobi";
  noiname="nobi";
  size=unique(scenarios_boxplot$mbc2);
  levelsvec=pbvec}
if (max(scenarios_boxplot$px)!=0){
  scenarios_boxplot = scenarios_boxplot[scenarios_boxplot$pb==0,]
  scenarios_boxplot = scenarios_boxplot[scenarios_boxplot$pe==0,]
  pname="px";
  noname="nox";
  noiname="noxi";
  if (max(scenarios_boxplot$Xa) == 0) {size=unique(scenarios_boxplot$alphac)}else{size=unique(scenarios_boxplot$mux)};
  levelsvec=pxvec}

#######################################################
# Extracting the information for beta from .Rdata files

# Dataframes for betahat, thetahat, and number of outliers 
boxplotBETA=NULL
boxplotTHETA=NULL
boxplotOUTLIER=NULL


for (i in 1:nrow(scenarios_boxplot)){
  nrep=scenarios_boxplot[i,1]
  nsample=scenarios_boxplot[i,2]
  ksample=scenarios_boxplot[i,3]
  pesample=scenarios_boxplot[i,4]
  pbsample=scenarios_boxplot[i,5]
  pxsample=scenarios_boxplot[i,6]
  mecsample=scenarios_boxplot[i,7]
  mbc2sample=scenarios_boxplot[i,8]
  alphacsample=scenarios_boxplot[i,9]
  rcsample=scenarios_boxplot[i,10]
  COMPindsample=as.logical(scenarios_boxplot[i,11])
  CCMindsample=as.logical(scenarios_boxplot[i,12])
  Xasample=as.logical(scenarios_boxplot[i,13])
  muxsample=scenarios_boxplot[i,14]
  Sclaudiosample=as.logical(scenarios_boxplot[i,15])
  
  if (CCMindsample==FALSE){stop("Contains ICM scenarios")}
  
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
  
  
  
  if ((scenarios_boxplot$pe[i]==0)&
      (scenarios_boxplot$pb[i]==0)&
      (scenarios_boxplot$px[i]==0)){
    flname=paste0("./Results_Uncontaminated/",flnameEst,"_",
                  "nrep=",scenarios_boxplot$nrep[i],"_",
                  "n=",scenarios_boxplot$n[i],"_",
                  "k=",scenarios_boxplot$k[i],"_",
                  "pe=",scenarios_boxplot$pe[i],"_",
                  "pb=",scenarios_boxplot$pb[i],"_",
                  "px=",scenarios_boxplot$px[i],"_",
                  "mec=",scenarios_boxplot$mec[i],"_",
                  "mbc2=",scenarios_boxplot$mbc2[i],"_",
                  "alphac=",scenarios_boxplot$alphac[i],"_",
                  "rc=",scenarios_boxplot$rc[i],".RData")}
  
  if (scenarios_boxplot$pe[i]>0){
    flname=paste0("./Results_Epsilon_contamination/",flnameEst,"_",
                  "nrep=",scenarios_boxplot$nrep[i],"_",
                  "n=",scenarios_boxplot$n[i],"_",
                  "k=",scenarios_boxplot$k[i],"_",
                  "pe=",scenarios_boxplot$pe[i],"_",
                  "pb=",scenarios_boxplot$pb[i],"_",
                  "px=",scenarios_boxplot$px[i],"_",
                  "mec=",scenarios_boxplot$mec[i],"_",
                  "mbc2=",scenarios_boxplot$mbc2[i],"_",
                  "alphac=",scenarios_boxplot$alphac[i],"_",
                  "rc=",scenarios_boxplot$rc[i],".RData")}
  
  if (scenarios_boxplot$pb[i]>0){
    flname=paste0("./Results_Random_Effect_contamination/",flnameEst,"_",
                  "nrep=",scenarios_boxplot$nrep[i],"_",
                  "n=",scenarios_boxplot$n[i],"_",
                  "k=",scenarios_boxplot$k[i],"_",
                  "pe=",scenarios_boxplot$pe[i],"_",
                  "pb=",scenarios_boxplot$pb[i],"_",
                  "px=",scenarios_boxplot$px[i],"_",
                  "mec=",scenarios_boxplot$mec[i],"_",
                  "mbc2=",scenarios_boxplot$mbc2[i],"_",
                  "alphac=",scenarios_boxplot$alphac[i],"_",
                  "rc=",scenarios_boxplot$rc[i],".RData")}
  
  if (scenarios_boxplot$px[i]>0 &  scenarios_boxplot$Xa[i] == 0){
    flname=paste0("./Results_X_contamination/",flnameEst,"_",
                  "nrep=",scenarios_boxplot$nrep[i],"_",
                  "n=",scenarios_boxplot$n[i],"_",
                  "k=",scenarios_boxplot$k[i],"_",
                  "pe=",scenarios_boxplot$pe[i],"_",
                  "pb=",scenarios_boxplot$pb[i],"_",
                  "px=",scenarios_boxplot$px[i],"_",
                  "mec=",scenarios_boxplot$mec[i],"_",
                  "mbc2=",scenarios_boxplot$mbc2[i],"_",
                  "alphac=",scenarios_boxplot$alphac[i],"_",
                  "rc=",scenarios_boxplot$rc[i],".RData")}
  
  if (scenarios_boxplot$px[i]>0 & scenarios_boxplot$Xa[i] == 1){
    flname=paste0("./Results_X_contamination/",flnameEst,"_",
                  "nrep=",scenarios_boxplot$nrep[i],"_",
                  "n=",scenarios_boxplot$n[i],"_",
                  "k=",scenarios_boxplot$k[i],"_",
                  "pe=",scenarios_boxplot$pe[i],"_",
                  "pb=",scenarios_boxplot$pb[i],"_",
                  "px=",scenarios_boxplot$px[i],"_",
                  "mec=",scenarios_boxplot$mec[i],"_",
                  "mbc2=",scenarios_boxplot$mbc2[i],"_",
                  "mux=",scenarios_boxplot$mux[i],"_",
                  "rc=",scenarios_boxplot$rc[i],".RData")}
  
  load(flname)
  
  # combining the extracted information for boxplots 
  # for estimators for beta in a dataframe
  
  # dataframe for MLE beta
  boxplotMLEbeta=cbind(rep(1,times=nrep),
                       rep(scenarios_boxplot[[pname]][i],times=nrep),
                       MLESMMTAUCOMPSTAUind$MLE$beta)
  
  # dataframe for S beta
  boxplotSbeta=cbind(rep(2,times=nrep),
                     rep(scenarios_boxplot[[pname]][i],times=nrep),
                     MLESMMTAUCOMPSTAUind$S$beta)
  
  # dataframe for MM beta
  boxplotMMbeta=cbind(rep(3,times=nrep),
                      rep(scenarios_boxplot[[pname]][i],times=nrep),
                      MLESMMTAUCOMPSTAUind$MM$beta)
  
  boxplotTaubeta=cbind(rep(4,times=nrep),
                       rep(scenarios_boxplot[[pname]][i],times=nrep),
                       MLESMMTAUCOMPSTAUind$Tau$beta)
  
  if (COMPindsample==TRUE){
    # dataframe for cTAU beta
    boxplotcTAUbeta=cbind(rep(6,times=nrep),
                          rep(scenarios_boxplot[[pname]][i],times=nrep),
                          MLESMMTAUCOMPSTAUind$COMPTau$beta)
    boxplotcSbeta=cbind(rep(5,times=nrep),
                        rep(scenarios_boxplot[[pname]][i],times=nrep),
                        MLESMMTAUCOMPSTAUind$COMPS$beta)
  }
  
  if (Sclaudiosample==TRUE){
    boxplotSclaudiobeta=cbind(rep(7,times=nrep),rep(scenarios_boxplot[[pname]][i],times=nrep),
                              MLESMMTAUCOMPSTAUind$Sclaudio$beta)}
  
  # Combining all dataframes for beta in one.
  # This dataframe can be uses as input for 
  # classical boxplots or ggplot2
  boxplotBETA=rbind(boxplotBETA,
                    boxplotMLEbeta,
                    boxplotSbeta,
                    boxplotMMbeta,
                    boxplotTaubeta)
  
  if (COMPindsample==TRUE){ boxplotBETA=rbind(boxplotBETA,boxplotcTAUbeta,boxplotcSbeta)}   
  if (Sclaudiosample==TRUE){ boxplotBETA=rbind(boxplotBETA,boxplotcTAUSclaudiobeta)}
  # dataframe for MLE theta
  boxplotMLEtheta=cbind(rep(1,times=nrep),
                        rep(scenarios_boxplot[[pname]][i],times=nrep),
                        MLESMMTAUCOMPSTAUind$MLE$theta)
  # dataframe for S theta
  boxplotStheta=cbind(rep(2,times=nrep),
                      rep(scenarios_boxplot[[pname]][i],times=nrep),
                      MLESMMTAUCOMPSTAUind$S$theta)
  boxplotTautheta=cbind(rep(4,times=nrep),
                        rep(scenarios_boxplot[[pname]][i],times=nrep),
                        MLESMMTAUCOMPSTAUind$Tau$theta)
  
  if (COMPindsample==TRUE){
    # dataframe for cTAU theta
    boxplotcTAUtheta=cbind(rep(6,times=nrep),
                           rep(scenarios_boxplot[[pname]][i],times=nrep),
                           MLESMMTAUCOMPSTAUind$COMPTau$theta)
    boxplotcStheta=cbind(rep(5,times=nrep),
                         rep(scenarios_boxplot[[pname]][i],times=nrep),
                         MLESMMTAUCOMPSTAUind$COMPS$theta)}
  if (Sclaudiosample==TRUE){
    boxplotSclaudiotheta=cbind(rep(7,times=nrep),
                               rep(scenarios_boxplot[[pname]][i],times=nrep),
                               MLESMMTAUCOMPSTAUind$Sclaudio$theta)
  }
  
  
  # Combining all dataframes for theta in one.
  # This dataframe can be uses as input for 
  # classical boxplots or ggplot2
  boxplotTHETA=rbind(boxplotTHETA,
                     boxplotMLEtheta,
                     boxplotStheta,
                     boxplotTautheta)
  
  if (COMPindsample==TRUE){
    boxplotTHETA=rbind(boxplotTHETA,boxplotcTAUtheta,boxplotcStheta)}
  if (Sclaudiosample==TRUE){
    boxplotTHETA=rbind(boxplotTHETA,boxplotSclaudiotheta)
  }
  
  # dataframe for number of outliers in i-th scenario
  boxplotOUTLIER_tmp=cbind(rep(scenarios_boxplot$n,nrep),
                           rep(scenarios_boxplot[[pname]][i],times=nrep),
                           MLESMMTAUCOMPSTAUind$no_outliers)
  
  # Combining with dataframes of previous scenarios
  boxplotOUTLIER=rbind(boxplotOUTLIER,
                       boxplotOUTLIER_tmp)
} # END of scenarios loop

###################################################################
# Setting the first two of boxplotBETA columns to factors

# Setting the first two of boxplotBETA columns to factors

colnames(boxplotBETA)=c("Estimator",pname,"beta1","beta2")
boxplotBETA=data.frame(boxplotBETA)
if (COMPind==TRUE){
  if (Sclaudio==TRUE){
    levBETA=1:7
    labBETA=c("MLE","S","MM","Tau","COMPS","COMPTau","Sclaudio")}else{
      levBETA=1:6
      labBETA=c("MLE","S","MM","Tau","COMPS","COMPTau") 
    }
}else{
  levBETA=1:4
  labBETA=c("MLE","S","MM","Tau")
}
boxplotBETA[,1]=factor(boxplotBETA[,1],
                       levels=levBETA,
                       labels=labBETA)
boxplotBETA[,2]=factor(boxplotBETA[,2],levels=levelsvec,
                       labels = as.character(levelsvec*100))

# Setting the first two of boxplotTHETA columns to factors
colnames(boxplotTHETA)=c("Estimator",pname,
                         "theta1","theta2","theta3","theta4")
boxplotTHETA=data.frame(boxplotTHETA)
if (COMPind==TRUE){
  if (Sclaudio==TRUE){
    levTHETA=c(1,2,4,5,6,7)
    labTHETA=c("MLE","S","Tau","COMPS","COMPTau","Sclaudio")}else{
      levTHETA=c(1,2,4,5,6)
      labTHETA=c("MLE","S","Tau","COMPS","COMPTau")
    }}else{
      levTHETA=c(1,2,4)
      labTHETA=c("MLE","S","Tau")
    }
boxplotTHETA[,1]=factor(boxplotTHETA[,1],
                              levels=levTHETA,
                        labels=labTHETA)
boxplotTHETA[,2]=factor(boxplotTHETA[,2],levels=levelsvec,
                              labels = as.character(levelsvec*100))

# Setting the first two of boxplotOUTLIER columns to factors
colnames(boxplotOUTLIER)=c("nsample",pname,
                                  "nobi","noei","noxi","noe","nox")
boxplotOUTLIER=data.frame(boxplotOUTLIER)
boxplotOUTLIER[,2]=factor(boxplotOUTLIER[,2],levels=levelsvec,
                                 labels = as.character(levelsvec*100))

################################################################
# Preparing graphs

# Boxplots using ggplot2
library(ggplot2)
library(gridExtra)
library(ggpubr)
titlename=paste0("CCM: rc=",unique(scenarios_boxplot[,10]),
                 ", n=",unique(scenarios_boxplot[,2]),
                 ", k=",unique(scenarios_boxplot[,3]),
                 ", ",pname, ", size=",size)

# Boxplots for beta
plotbeta1=ggplot(boxplotBETA,aes(x=get(pname),y=beta1,fill=Estimator))+
  geom_boxplot()+
  xlab("Fixed % of contaminated cases")+
  ylab(quote(hat(beta)[1]))+
  geom_hline(yintercept=250,lty=1,col="orange")#+
  #ggtitle(titlename)

plotbeta2=ggplot(boxplotBETA,aes(x=get(pname),y=beta2,fill=Estimator))+
  geom_boxplot()+
  xlab("Fixed % of contaminated cases")+
  ylab(quote(hat(beta)[2]))+
  geom_hline(yintercept=10,lty=1,col="orange")#+
  #ggtitle(titlename)

ggarrange(plotbeta1,plotbeta2,nrow=1,common.legend = TRUE, legend = "bottom")


#pdf("tmp.pdf")
#print(plotbeta1)
#print(plotbeta2)
#dev.off()

################################################################
# Boxplots for theta

plottheta1=ggplot(boxplotTHETA,aes(x=get(pname),y=theta1,fill=Estimator))+
  geom_boxplot()+
  xlab("Fixed % of contaminated cases")+
  ylab(quote(hat(sigma)[gamma[1]]^2))+
  geom_hline(yintercept=790,lty=1,col="orange")#+
  #ggtitle(titlename)

plottheta2=ggplot(boxplotTHETA,aes(x=get(pname),y=theta2,fill=Estimator))+
  geom_boxplot()+
  xlab("Fixed % of contaminated cases")+
  ylab(quote(hat(sigma)[gamma[12]]))+
  geom_hline(yintercept=-8.5,lty=1,col="orange")#+
  #ggtitle(titlename)

plottheta3=ggplot(boxplotTHETA,aes(x=get(pname),y=theta3,fill=Estimator))+
  geom_boxplot()+
  ylab(quote(hat(sigma)[gamma[2]]^2))+
  xlab("Fixed % of contaminated cases")+
  geom_hline(yintercept=40,lty=1,col="orange")#+
  #ggtitle(titlename)

plottheta4=ggplot(boxplotTHETA,aes(x=get(pname),y=theta4,fill=Estimator))+
  geom_boxplot()+
  ylab(quote(hat(sigma)[epsilon]^2))+
  xlab("Fixed % of contaminated cases")+
  geom_hline(yintercept=400,lty=1,col="orange")#+
  #ggtitle(titlename)

ggarrange(plottheta1,plottheta2,plottheta3,plottheta4,nrow=2,ncol=2,common.legend = TRUE, legend = "bottom")


#grid.arrange(plotbeta1,plottheta1,plottheta2,plotbeta2,plottheta3,plottheta4,nrow=2,ncol=3)

################################################################
# Boxplots for number of outliers

plotnoei=ggplot(boxplotOUTLIER,aes(x=get(pname),y=get(noiname)/nsample))+
  geom_boxplot()+
  xlab("Fixed % of contaminated cells")+
  ylab("Proportion of contaminated cases")+
  ylim(c(0,1))+
  ggtitle(titlename)

plotnoe=ggplot(boxplotOUTLIER,aes(x=get(pname),y=get(noname)/(nsample*ksample)))+
  geom_boxplot()+
  xlab("Fixed % of contaminated cells")+
  ylab("Proportion of contaminated cells")+
  ylim(c(0,1))+
  ggtitle(titlename)

grid.arrange(plotnoei,plotnoe,nrow=1)

