# This script prepares the setting for the simulation for 
# the MLE, S and MM estimators with ICM contamination
# It sources the required R-scripts.
# 
# The user first specifies the contamination scenarios.
# The script runs the simulation for these scenarios using
# the function MLESMM_estimates_MCG and data_gen_MCR
#
# After the simulation for a scenario has been completed
# this script saves the results in .RData files containing a list MLESMM

# The .Rdata files are stored in designated subfolders
# BEFORE running this script, first create subfolders
# - Results_Uncontaminated
# - Results_Epsilon_contamination
# - Results_Random_Effect_contamination
# - Results_X_contamination


source("biweight_functions.R")
source("asympt_norm_constants.R")
source("function_data_gen_MCG.R")
source("Robust_lme.R")
source("function_MLESMM_estimates_MCG.R")

library(robustbase)

# Information for ourselves
# beta=c(250,10)
# beta1=gamma0
# beta2=gamma1

# theta=c(790,-8.5,40,400)
# theta1=sigma0^2
# theta2=sigma10
# theta3=sigma1^2
# theta4=sigmaeps^2

pe=0.05
pb=0.05
px=0.05
nrep=1

# creating a dataframe that contains different contamination schemes as rows:
# n number of individuals
# k number of observations per individual
# pe probability of having outlier in component epsilon_ij
# pb probability of having outlier in vector of random effects b
# px probability of having outlier in component x_ij in X
# mec shift in the mean of component epsilon_ij
# mbc2 shift in the mean of random effect b2
# alphac multiplication factor in component x_ij in X

scenarios=NULL
scenarios=rbind(scenarios,c(nrep,200,4,pe,0,0,-80,0,1))
scenarios=rbind(scenarios,c(nrep,200,4,0,pb,0,0,-25,1))
scenarios=rbind(scenarios,c(nrep,200,4,0,0,px,0,0,10))

colnames(scenarios)=c("nrep","n","k","pe","pb","px","mec","mbc2","alphac")
scenarios=data.frame(scenarios)

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

  MLESMM=MLESMM_estimates_MCG(nrep=nrep,n=nsample,k=ksample,
                              pe=pesample,pb=pbsample,px=pxsample,
                              mec=mecsample,mbc2=mbc2sample,
                              alphac=alphacsample)

  if (pe>0){
  flname=paste0("./Results_Epsilon_contamination/","MLESMM","_",
                "nrep=",nrep,"_",
                "n=",nsample,"_",
                "k=",ksample,"_",
                "pe=",pesample,"_",
                "pb=",pbsample,"_",
                "px=",pxsample,"_",
                "mec=",mecsample,"_",
                "mbc2=",mbc2sample,"_",
                "alphac=",alphacsample,"_",
                ".RData")}else{
  flname=paste0("./Results_Uncontaminated/","MLESMM","_",
                "nrep=",nrep,"_",
                "n=",nsample,"_",
                "k=",ksample,"_",
                "pe=",pesample,"_",
                "pb=",pbsample,"_",
                "px=",pxsample,"_",
                "mec=",mecsample,"_",
                "mbc2=",mbc2sample,"_",
                "alphac=",alphacsample,"_",
                ".RData")}
  
  if (pb>0){
    flname=paste0("./Results_Random_Effect_contamination/","MLESMM","_",
                  "nrep=",nrep,"_",
                  "n=",nsample,"_",
                  "k=",ksample,"_",
                  "pe=",pesample,"_",
                  "pb=",pbsample,"_",
                  "px=",pxsample,"_",
                  "mec=",mecsample,"_",
                  "mbc2=",mbc2sample,"_",
                  "alphac=",alphacsample,"_",
                  ".RData")}
                    
  if (px>0){
    flname=paste0("./Results_X_contamination/","MLESMM","_",
                  "nrep=",nrep,"_",
                  "n=",nsample,"_",
                  "k=",ksample,"_",
                  "pe=",pesample,"_",
                  "pb=",pbsample,"_",
                  "px=",pxsample,"_",
                  "mec=",mecsample,"_",
                  "mbc2=",mbc2sample,"_",
                  "alphac=",alphacsample,"_",
                  ".RData")}
  
  save(MLESMM,file=flname)

  }


