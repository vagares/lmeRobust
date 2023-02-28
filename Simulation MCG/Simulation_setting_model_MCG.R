# This script prepares the setting for the simulation.
# It sources the required R-scripts.
# 
# This script also sets the contamination scenarios.
# It runs the simulation for that scenario.
# 
# After the simulation for a scenario has been completed
# this script saves the results in .RData files 
# containing a list MLESMM

source("biweight_functions.R")
source("asympt_norm_constants.R")
source("Sim_data_MCG.R")
source("Robust_lme.R")
source("Simulation_Estimates_model_MCG.R")

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

################################################
# SCENARIO: - shift fixed at -80
#           - pe 0, 0.01, .... ,0.10
################################################

pevec=seq(0,0.10,by=0.01)
nrep=250

# creating a dataframe that contains different contamination schemes
# as rows

scenarios=NULL
for (i in 1:length(pevec)){
  scenarios=rbind(scenarios,c(nrep,200,4,pevec[i],0,0,-80,0,1))
  #scenarios[i,]=c(nrep,200,4,pevec[i],0,0,-80,0,1)
}
colnames(scenarios)=c("nrep","n","k","pe","pb","px","mec","mbc2","alphac")
scenarios=data.frame(scenarios)

set.seed(11131957)

for (i in (1:nrow(scenarios))){
  nrep=scenarios[i,1]
  n=scenarios[i,2]
  k=scenarios[i,3]
  pe=scenarios[i,4]
  pb=scenarios[i,5]
  px=scenarios[i,6]
  mec=scenarios[i,7]
  mbc2=scenarios[i,8]
  alphac=scenarios[i,9]
  
  MLESMM=MLESMM_estimates_MCG(nrep=nrep,n=n,k=k,pe=pe,pb=pb,px=px,
                       mec=mec,mbc2=mbc2,alphac=alphac)
  if (pe>0){
     flname=paste0("./Results_Epsilon_contamination/","MLESMM","_",
                   "nrep=",nrep,"_",
                   "n=",n,"_",
                   "k=",k,"_",
                   "pe=",pe,"_",
                   "pb=",pb,"_",
                   "px=",px,"_",
                   "mec=",mec,"_",
                   "mbc2=",mbc2,"_",
                   "alphac=",alphac,".RData")}else{
      flname=paste0("./Results_Uncontaminated/","MLESMM","_",
                    "nrep=",nrep,"_",
                    "n=",n,"_",
                    "k=",k,"_",
                    "pe=",pe,"_",
                    "pb=",pb,"_",
                    "px=",px,"_",
                    "mec=",mec,"_",
                    "mbc2=",mbc2,"_",
                    "alphac=",alphac,".RData")}
  
  save(MLESMM,file=flname)
  
}


