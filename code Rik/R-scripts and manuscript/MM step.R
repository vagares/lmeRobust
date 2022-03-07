eff=0.95
objectf=function(c){
    1/constbetahat(k,c)-eff 
  }
c1=uniroot(objectf,interval=c(0.001,100))$root

# SETTING CUT-VALUES FOR MM WITH TUKEYS BIWEIGHT
M=0
c0t=c1

betaMMold=betaCBSmetalRESCALE
theta=thetaCBSmetalRESCALE
sqrt(theta)
Vtheta=theta[1]*L1+theta[2]*L2+theta[3]*L3++theta[4]*L4
eigen(Vtheta)


###########################################################################
# RUNNING THE ITERATION TO OBTAIN S-ESTIMATES
###########################################################################

# stop criterion
tol=1
iter=0

while (tol>=1e-06){
  ##############################################################################
  # Iteration step
  #############################################
  
  # Re-scaling Mahalanobis distances to satisfy S-constraint
  # array for n Mahalanobis distances
  MD=numeric(n)
  for (i in 1:n){
    #X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
    X=ifelse(Type[i]==1,1,0)*X1+ifelse(Type[i]==2,1,0)*X2
    y=as.matrix(Ymat[i,])
    #MD[i]=mahalanobis(y-X%*%betaMMold,center=FALSE,cov=Vtheta) # Note MD=d^2!
    MD[i]=t(y-X%*%betaMMold)%*%solve(Vtheta)%*%(y-X%*%betaMMold)
  }
  
  # determining scaling constant for MD
  # for translated biweight
  objfuntranslated=function(s){
    mean(biweightrhotranslated(sqrt(MD)/s,M,c0t))-expecrhotranslated(k,M,c0t)
  }
  MDscaletranslated=uniroot(f=objfuntranslated,c(0.01,100))$root
  
  
  # OPTION 1: re-scaling MD to satisfy S-constraint
  #MDtildetranslated=MD/MDscaletranslated^2
  
  # OPTION 2: without adapting the Mahalanobis distance
  MDtildetranslated=MD
  
  # OPTION 3: Rocke's adaptation as in Copt
  # CORRECTION USED IN PHD COPT? SEE (4.14)
  qconst=round((n+k+1)/2)
  kconstsq=sort(MD)[qconst]/qchisq(qconst/(n+1),df=k)
  #MDtildetranslated=MD/kconstsq
  
  # prepare for updating regression estimate in step 5
  betavectermtranslated=matrix(0,nrow=q,ncol=1) # vector-term in step 5
  betamattermtranslated=matrix(0,nrow=q,ncol=q) # matrix-term in step 5

    for (i in 1:n){
    #X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
    X=ifelse(Type[i]==1,1,0)*X1+ifelse(Type[i]==2,1,0)*X2
    y=matrix(as.numeric(Ymat[i,]),ncol=1,nrow=k)
    betavectermtranslated=
      betavectermtranslated+
      biweightutranslated(sqrt(MDtildetranslated[i]),M,c0t)*
      t(X)%*%solve(Vtheta)%*%y
    
    betamattermtranslated=
      betamattermtranslated+
      biweightutranslated(sqrt(MDtildetranslated[i]),M,c0t)*t(X)%*%solve(Vtheta)%*%X 
    
  }
  
  # STEP5 updating regression estimate
  betaMMnew=solve(betamattermtranslated)%*%betavectermtranslated
  # End of iteration step
  
  iter=iter+1
  
  diffbeta=norm(betaMMnew-betaMMold,type="F")
  tol=diffbeta
  
  # setting starting values for the next iteration within while-loop
  betaMMold=betaMMnew
}
###########################################################################

iter
betaMMnew

# RESISTANCE MM
betaCBSMMresistanceRESCALE=betaMMnew

# SEMANTIC MM
betaCBSMMsemanticRik=betaMMnew
betaCBSMMsemanticRESCALE=betaMMnew

# GROWTH MM
betaCBSMMgrowthRik=betaMMnew
betaCBSMMgrowthRESCALE=betaMMnew

# METALLIC MM
betaCBSMMmetalRik=betaMMnew
betaCBSMMmetalRESCALE=betaMMnew




