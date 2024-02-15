# This script prepares the setting for the simulation for 
# the MLE, S, MM, and cTAU estimators
# 
# The user has to specify 
#   - whether or not to incorporate the cTAU estimator
#   - whether to use ICM or CCM contamination scenarios
#   - whether to have a random or fixed number of outliers

# This script file prepares the contamination scenarios for the 
# simulations of the MCG model. 
# The file creats a dataframe, that can be sourced in the file
# - Simulation_setting_model_MCG_CCMind_cTAUind.R

# Set options for 
# cTAU included yes (cTAU==TRUE) or no (cTAU==FALSE)

cTAUind=TRUE
Sclaudio = TRUE
#######################################################################
# creating dataframe containing different contamination schemes as rows

pevec=0  # proportion of contamination in error
pbvec=0  # proportion of contamination in random effect
pxvec=0  # proportion of contamination in X-matrix
k=4       # dimension of y
#k=8
n=200     # sample size
nrep=250  # number of repetitions

scenarios=NULL
for (cs in 6:6){
  # 3 different contamination scenarios 
  # ICM and CCM with randcont 0/1, Xa=TRUE/FALSE,Xshiftall=TRUE/FALSE
  
 if (cs==6) {CCMind=TRUE;Sclaudio=TRUE;randcont=1;Xa=FALSE;Xshiftall=TRUE}
    # ONLY pevec varying

 scenarios=rbind(scenarios,c(nrep,n,k,0,0,0,0,0,1,
                      randcont,as.integer(cTAUind),as.integer(CCMind),
                      as.integer(Xa),as.integer(Xshiftall),0,as.integer(Sclaudio)))


} # END cs-loop 

colnames(scenarios)=c("nrep","n","k","pe","pb","px","mec","mbc2",
                              "alphac","rc","cTAUind","CCMind","Xa","Xshiftall","mux","Sclaudio")

scenarios=data.frame(scenarios)

