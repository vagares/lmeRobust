# This script uses the contamination scenarios and results of
# the corresponding simulation performed by
# Simulation_setting_model_MCG.R

# The script extracts the information from the .RData files
# corresponding to the different setups in dataframe 
# scenarios_boxplot and produces boxplots of the estimators

#######################################################################
# creating dataframe containing different contamination schemes as rows
# SCENARIO 1: - -80 fixed
#             - pe 0, 0.01, ... , 0.10
#             - 250 repetitions

pevec=seq(0,0.10,by=0.05)
nrep=100
scenarios_boxplot=NULL
for (i in 1:length(pevec)){
  scenarios_boxplot=rbind(scenarios_boxplot,c(nrep,200,4,pevec[i],0,0,-80,0,1))
}
colnames(scenarios_boxplot)=c("nrep","n","k","pe","pb","px","mec","mbc2","alphac")

# Sets scenarios equal to scenarios_boxplot
# scenarios is used below to extract the .RData files
scenarios=data.frame(scenarios_boxplot)


#######################################################
# Extracting the information for beta from .Rdata files

# Dataframes for betahat, thetahat, and number of outliers 
boxplotMLESMMcTAUbeta=NULL
boxplotMLESMMcTAUtheta=NULL
boxplotMLESMMcTAU_outlier=NULL

for (i in 1:nrow(scenarios)){
if (scenarios$pe[i]>0){
    flname=paste0("./Results_Epsilon_contamination/","MLESMMcTAU","_",
#   flname=paste0("./Results_Epsilon_contamination/EpsilonCCM/","MLESMMcTAU","_",
  "nrep=",scenarios$nrep[i],"_",
                "n=",scenarios$n[i],"_",
                "k=",scenarios$k[i],"_",
                "pe=",scenarios$pe[i],"_",
                "pb=",scenarios$pb[i],"_",
                "px=",scenarios$px[i],"_",
                "mec=",scenarios$mec[i],"_",
                "mbc2=",scenarios$mbc2[i],"_",
                "alphac=",scenarios$alphac[i],".RData")}else{
                  flname=paste0("./Results_Uncontaminated/","MLESMMcTAU","_",
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

# combining the extracted information for boxplots 
# for estimators for beta in a dataframe

# dataframe for MLE beta
boxplotMLEbeta=cbind(rep(1,times=nrep),
                     rep(scenarios$pe[i],times=nrep),
                     MLESMMcTAU$MLE$beta)


# dataframe for S beta
boxplotSbeta=cbind(rep(2,times=nrep),
                     rep(scenarios$pe[i],times=nrep),
                     MLESMMcTAU$S$beta)
# dataframe for MM beta
boxplotMMbeta=cbind(rep(3,times=nrep),
                     rep(scenarios$pe[i],times=nrep),
                     MLESMMcTAU$MM$beta)

# dataframe for cTAU beta
boxplotcTAUbeta=cbind(rep(4,times=nrep),
                    rep(scenarios$pe[i],times=nrep),
                    MLESMMcTAU$cTAU$beta)


# Combining all dataframes for beta in one.
# This dataframe can be uses as input for 
# classical boxplots or ggplot2
boxplotMLESMMcTAUbeta=rbind(boxplotMLESMMcTAUbeta,
                        boxplotMLEbeta,
                        boxplotSbeta,boxplotMMbeta,boxplotcTAUbeta)

# dataframe for MLE theta
boxplotMLEtheta=cbind(rep(1,times=nrep),
                      rep(scenarios$pe[i],times=nrep),
                      MLESMMcTAU$MLE$theta)
# dataframe for S theta
boxplotStheta=cbind(rep(2,times=nrep),
                    rep(scenarios$pe[i],times=nrep),
                    MLESMMcTAU$S$theta)

# dataframe for cTAU theta
boxplotcTAUtheta=cbind(rep(4,times=nrep),
                    rep(scenarios$pe[i],times=nrep),
                    MLESMMcTAU$cTAU$theta)


# Combining all dataframes for theta in one.
# This dataframe can be uses as input for 
# classical boxplots or ggplot2
boxplotMLESMMcTAUtheta=rbind(boxplotMLESMMcTAUtheta,
                         boxplotMLEtheta,
                         boxplotStheta,boxplotcTAUtheta)

# dataframe for number of outliers in i-th scenario
boxplotMLESMMcTAU_outlier_tmp=cbind(rep(n,nrep),
                            rep(scenarios$pe[i],times=nrep),
                            MLESMMcTAU$no_outliers)

# Combining with dataframes of previous scenarios
boxplotMLESMMcTAU_outlier=rbind(boxplotMLESMMcTAU_outlier,
                            boxplotMLESMMcTAU_outlier_tmp)
} # END of scenarios loop

###################################################################
# Setting the first two of boxplotMLESMMbeta columns to factors

colnames(boxplotMLESMMcTAUbeta)=c("Estimator","pe","beta1","beta2")
boxplotMLESMMcTAUbeta=data.frame(boxplotMLESMMcTAUbeta)
boxplotMLESMMcTAUbeta[,1]=factor(boxplotMLESMMcTAUbeta[,1],
                             levels=1:4,labels=c("MLE","S","MM","cTAU"))
boxplotMLESMMcTAUbeta[,2]=factor(boxplotMLESMMcTAUbeta[,2],levels=pevec,
                             labels = as.character(pevec*100))

# Setting the first two of boxplotMLESMMcTAUtheta columns to factors
colnames(boxplotMLESMMcTAUtheta)=c("Estimator","pe",
                               "theta1","theta2","theta3","theta4")
boxplotMLESMMcTAUtheta=data.frame(boxplotMLESMMcTAUtheta)
boxplotMLESMMcTAUtheta[,1]=factor(boxplotMLESMMcTAUtheta[,1],
                             levels=c(1:2,4),labels=c("MLE","S","cTAU"))
boxplotMLESMMcTAUtheta[,2]=factor(boxplotMLESMMcTAUtheta[,2],levels=pevec,
                              labels = as.character(pevec*100))

# Setting the first two of boxplotMLESMMcTAU_outlier columns to factors
colnames(boxplotMLESMMcTAU_outlier)=c("nsample","pe",
                                  "nobi","noei","noxi","noe","nox")
boxplotMLESMMcTAU_outlier=data.frame(boxplotMLESMMcTAU_outlier)
boxplotMLESMMcTAU_outlier[,2]=factor(boxplotMLESMMcTAU_outlier[,2],levels=pevec,
                                 labels = as.character(pevec*100))

################################################################
# Preparing graphs

# Boxplots using ggplot2
library(ggplot2)
library(gridExtra)

# Boxplots for beta
plotbeta1=ggplot(boxplotMLESMMcTAUbeta,aes(x=pe,y=beta1,fill=Estimator))+
  geom_boxplot()+
  xlab("Percentage of contamination")+
  geom_hline(yintercept=250,lty=1,col="orange")


plotbeta2=ggplot(boxplotMLESMMcTAUbeta,aes(x=pe,y=beta2,fill=Estimator))+
  geom_boxplot()+
  xlab("Percentage of contamination")+
  geom_hline(yintercept=10,lty=1,col="orange")

grid.arrange(plotbeta1,plotbeta2,nrow=1)

################################################################
# Boxplots for theta

plottheta1=ggplot(boxplotMLESMMcTAUtheta,aes(x=pe,y=theta1,fill=Estimator))+
  geom_boxplot()+
  xlab("Percentage of contamination")+
  geom_hline(yintercept=790,lty=1,col="orange")

plottheta2=ggplot(boxplotMLESMMcTAUtheta,aes(x=pe,y=theta2,fill=Estimator))+
  geom_boxplot()+
  xlab("Percentage of contamination")+
  geom_hline(yintercept=-8.5,lty=1,col="orange")

plottheta3=ggplot(boxplotMLESMMcTAUtheta,aes(x=pe,y=theta3,fill=Estimator))+
  geom_boxplot()+
  xlab("Percentage of contamination")+
  geom_hline(yintercept=40,lty=1,col="orange")

plottheta4=ggplot(boxplotMLESMMcTAUtheta,aes(x=pe,y=theta4,fill=Estimator))+
  geom_boxplot()+
  xlab("Percentage of contamination")+
  geom_hline(yintercept=400,lty=1,col="orange")

grid.arrange(plottheta1,plottheta2,plottheta3,plottheta4,
             nrow=2,ncol=2)


grid.arrange(plotbeta1,plottheta1,plottheta2,
             plotbeta2,plottheta3,plottheta4,
             nrow=2,ncol=3)

################################################################
# Boxplots for number of outliers

plotnoei=ggplot(boxplotMLESMMcTAU_outlier,aes(x=pe,y=noei/scenarios$n))+
  geom_boxplot()+
  xlab("Percentage of contamination")+
  ylab("Proportion of outlier individuals")+
  ylim(c(0,0.5))+
  ggtitle("n=200 and k=4")

plotnoe=ggplot(boxplotMLESMMcTAU_outlier,aes(x=pe,y=noe/(scenarios$n*scenarios$k)))+
  geom_boxplot()+
  xlab("Percentage of contamination")+
  ylab("Proportion of outlier observations")+
  ylim(c(0,0.5))+
  ggtitle("n=200 and k=4")

grid.arrange(plotnoei,plotnoe,nrow=1)