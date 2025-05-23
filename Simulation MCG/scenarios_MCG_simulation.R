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

#######################################################################
# creating dataframe containing different contamination schemes as rows

#pevec=c(0,0.05,0.10,0.20,0.30)  # proportion of contamination in error
#pbvec=c(0,0.05,0.10,0.20,0.30)  # proportion of contamination in random effect
#pxvec=c(0,0.05,0.10,0.20,0.30)  # proportion of contamination in X-matrix
#mecvec=c(-40,-80,-160,-500,-800)          # size of contamination-shift in error
#mbc2vec=c(-12.5,-25,-50)        # size of contamination-shift in random effect
#alphacvec=c(2,5,10,50,100)             # size of contamination-factor in X
pevec=c(0,0.10)  # proportion of contamination in error
pbvec=c(0,0.10)  # proportion of contamination in random effect
pxvec=c(0,0.10)  # proportion of contamination in X-matrix
mecvec=c(-40,-160)          # size of contamination-shift in error
mbc2vec=c(-12.5,-25)        # size of contamination-shift in random effect
alphacvec=c(2,10)             # size of contamination-factor in X

k=2      # dimension of y
#k=8
n=200     # sample size
nrep=50    # number of repetitions

scenarios=NULL
for (cs in 1:4){
  # 3 different contamination scenarios 
  # ICM and CCM with randcont 0/1
  
  if (cs==1) {CCMind=FALSE;randcont=1} # value of randcont has no influence
  if (cs==2) {CCMind=TRUE;randcont=1}
  if (cs==3) {CCMind=TRUE;randcont=0}

    # ONLY pevec varying
    for (i in 1:length(pevec)){
      for (j in 1:length(mecvec)){
      scenarios=rbind(scenarios,c(nrep,n,k,pevec[i],0,0,mecvec[j],0,1,
                      randcont,as.integer(cTAUind),as.integer(CCMind)))
      } # END j-loop mecvec
    } # END i-loop pevec

    # ONLY pbvec varying
  if (max(pbvec)>0){
    for (i in 1:length(pbvec)){
      for (j in 1:length(mbc2vec)){
        scenarios=rbind(scenarios,c(nrep,n,k,0,pbvec[i],0,0,mbc2vec[j],1,
                       randcont,as.integer(cTAUind),as.integer(CCMind)))
      } # END j-loop mbc2vec
    } # END i-loop pbvec
  }

    # ONLY pxvec varying
  if (max(pxvec)>0){ 
   for (i in 1:length(pxvec)){
      for (j in 1:length(alphacvec)){
        scenarios=rbind(scenarios,c(nrep,n,k,0,0,pxvec[i],0,0,alphacvec[j],
                      randcont,as.integer(cTAUind),as.integer(CCMind)))
      } # END j-loop alphacvec
    } # END i-loop pxvec
  }
} # END cs-loop 

colnames(scenarios)=c("nrep","n","k","pe","pb","px","mec","mbc2",
                              "alphac","rc","cTAUind","CCMind")

scenarios=data.frame(scenarios)

