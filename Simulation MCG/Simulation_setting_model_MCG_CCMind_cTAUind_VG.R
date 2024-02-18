# The script runs the simulation for these scenarios using
# the functions 
#   - MLESMMCOMPind_estimates_MCG_CCMind 
#   - data_gen_MCG_CCMind

# It can run from an empty environment and sources the required R-scripts.

# After the simulation for a scenario has been completed
# this script saves the results in .RData files 
# containing a list by the name of
#   - MLESMM_ICM
#   - MLESMM_CCM
#   - MLESMMCOMP_ICM
#   - MLESMMCOMP_ICM
# depending on the chosen scenario

# The .Rdata files are stored in designated subfolders
# BEFORE running this script, first create subfolders
# - Results_Uncontaminated
# - Results_Epsilon_contamination
# - Results_Random_Effect_contamination
# - Results_X_contamination

# This script sources the file
#   - scenarios_MCG_simulation.R
# which creates the different scenarios and puts them in a dataframe
#   - scenarios
# BEFORE running this script, 
# check whether the dataframe has created properly. 

library(robustbase)
library(lava)

source("biweight_functions.R")
source("asympt_norm_constants.R")
source("function_data_gen_MCG_CCMind_VG.R")
source("Robust_lme.R")
source("function_MLESMMcTAUind_estimates_MCG_CCMind_VG.R")
source("scenarios_MCG_simulation_VG.R")


# Information for ourselves
# beta=c(250,10)
# beta1=gamma0

# creating a dataframe that contains different contamination schemes as rows:
# n number of individuals
# k number of observations per individual
# pe probability of having outlier in component epsilon_ij
# pb probability of having outlier in vector of random effects b
# px probability of having outlier in component x_ij in X
# mec shift in the mean of component epsilon_ij
# mbc2 shift in the mean of random effect b2
# alphac multiplication factor in component x_ij in X
# rc index to use random number (1) or fixed number (0) of outliers
# COMPind index to incorporate COMP (TRUE) or not (FALSE)
# CCMind index to use CCM (TRUE) or ICM (FALSE)


set.seed(11131957)

for (i in (1:nrow(scenarios))){
  nrep=scenarios[i,1]
  nsample=scenarios[i,2]
  ksample=scenarios[i,3]
  pesample=scenarios[i,4]
  pbsample=scenarios[i,5]
  pxsample=scenarios[i,6]
  mecsample=scenarios[i,7]
  mbc2sample=scenarios[i,8]
  alphacsample=scenarios[i,9]
  rcsample=scenarios[i,10]
  COMPindsample=as.logical(scenarios[i,11])
  CCMindsample=as.logical(scenarios[i,12])
  Xasample=as.logical(scenarios[i,13])
  muxsample=scenarios[i,14]
  Sclaudiosample=as.logical(scenarios[i,15])
  # Setting the filename depending on yes/no COMP and yes/no CCM
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
  MLESMMTAUCOMPSTAUind=MLESMMTAUCOMPSTAUind_estimates_MCG_CCMind(nrep=nrep,n=nsample,k=ksample,
                                      pe=pesample,pb=pbsample,px=pxsample,
                                      mec=mecsample,mbc2=mbc2sample,
                                      alphac=alphacsample,
                                      randcont=rcsample,
                                      COMPind=COMPindsample,
                                      CCMind=CCMindsample,
                                      Xa=Xasample,
                                      mux=muxsample,
                                      Sclaudioind=Sclaudiosample)
  
  if ((pesample==0)&
      (pbsample==0)&
      (pxsample==0)){
      flname=paste0("./Results_Uncontaminated/",flnameEst,"_",
                    "nrep=",nrep,"_",
                    "n=",nsample,"_",
                    "k=",ksample,"_",
                    "pe=",pesample,"_",
                    "pb=",pbsample,"_",
                    "px=",pxsample,"_",
                    "mec=",mecsample,"_",
                    "mbc2=",mbc2sample,"_",
                    "alphac=",alphacsample,"_",
                    "rc=",rcsample,
                    ".RData")}

  if (pesample>0){
    flname=paste0("./Results_Epsilon_contamination/",flnameEst,"_",
                  "nrep=",nrep,"_",
                  "n=",nsample,"_",
                  "k=",ksample,"_",
                  "pe=",pesample,"_",
                  "pb=",pbsample,"_",
                  "px=",pxsample,"_",
                  "mec=",mecsample,"_",
                  "mbc2=",mbc2sample,"_",
                  "alphac=",alphacsample,"_",
                  "rc=",rcsample,
                  ".RData")}

    if (pbsample>0){
    flname=paste0("./Results_Random_Effect_contamination/",flnameEst,"_",
                  "nrep=",nrep,"_",
                  "n=",nsample,"_",
                  "k=",ksample,"_",
                  "pe=",pesample,"_",
                  "pb=",pbsample,"_",
                  "px=",pxsample,"_",
                  "mec=",mecsample,"_",
                  "mbc2=",mbc2sample,"_",
                  "alphac=",alphacsample,"_",
                  "rc=",rcsample,
                  ".RData")}
  
  if (pxsample>0){
    if(Xasample==FALSE){flname=paste0("./Results_X_contamination/",flnameEst,"_",
                  "nrep=",nrep,"_",
                  "n=",nsample,"_",
                  "k=",ksample,"_",
                  "pe=",pesample,"_",
                  "pb=",pbsample,"_",
                  "px=",pxsample,"_",
                  "mec=",mecsample,"_",
                  "mbc2=",mbc2sample,"_",
                  "alphac=",alphacsample,"_",
                  "rc=",rcsample,
                  ".RData")}else
                  {
                    flname=paste0("./Results_X_contamination/",flnameEst,"_",
                                  "nrep=",nrep,"_",
                                  "n=",nsample,"_",
                                  "k=",ksample,"_",
                                  "pe=",pesample,"_",
                                  "pb=",pbsample,"_",
                                  "px=",pxsample,"_",
                                  "mec=",mecsample,"_",
                                  "mbc2=",mbc2sample,"_",
                                  "mux=",muxsample,"_",
                                  "rc=",rcsample,
                                  ".RData") 
                    
                  }}
  
save(MLESMMTAUCOMPSTAUind,file=flname)

print(i)  
}

