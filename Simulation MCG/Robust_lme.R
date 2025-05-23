# This script contains the code for the function Roblme. 
# This function computes the MLE, the S- and MM-estimates 
# for a single dataset, such as the ones in the output of
# the functions data_gen_MCG or data_gen_MCG_contCCM.

library(robustbase) # needed for covMcd

source("biweight_functions.R")
source("asympt_norm_constants.R")

Roblme = function(Ymat,X,Z,E=NULL,L=NULL,rho ="t-biweight",r =0.5,arp=0.01,rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL){
  # Ymat: nxk-matrix of observations: n individuals, dimension k
  # X: a list X[[i]] i = 1...n design matrix fixed effects
  # Z: a list Z[[z]] design matrices random effects
  # E: variance of error (default is identity)
  # L: list of matrices in linear decomposition of V (including error variance)
  
  # rho: rho function for the S estimation (t-biweight, biweight or MLE)
  # r: breakdown point
  # arp: asymptotic rejection rate
  # rhoMM: rho function for the MM estimation (biweight)
  # eps: bound for the stopping criteria
  # maxiter: maximum number of iteration
  # eff: efficiency
  # V0: initial covariance matrix
  
  lZ = length(Z) 
  lX = length(X)
  k = ncol(Ymat)
  n = nrow(Ymat)

  # Add identity matrix for error to Z-list or use L-list
  if(is.null(L) ==TRUE) {
    if(is.null(E) ==TRUE) {Z[[lZ+1]] =  diag(rep(1,k)) 
    }else{Z[[lZ+1]] =  E}
    lZ = length(Z)}else{lZ = length(L)}

  # Setting BDP and cut-off constant for translated biweight
  # and setting sigma1 and sigma2 in limiting variance of thetahat
  if (rho =="t-biweight"){
    c0 = rhotranslatedconst(k,r,arp,0.0001,1000)
    m0 = sqrt(qchisq(1-arp,df=k))-c0
    b0 = expecrhotranslated(k,m0,c0)
    s1 = sigma1t(k,m0,c0)
    s2 = sigma2t(k,m0,c0) 
    }
  if (rho =="biweight"){
    m0 = 0
    c0 = rhoconst(k,r,0.0001,10000)
    b0 = expecrho(k,c0)
    s1 = sigma1(k,c0)
    s2 = sigma2(k,c0)
    }
  if (rho =="MLE"){
    m0 = 10000
    c0 = rhotranslatedconst(k,r,arp,0.0001,1000)
    b0 = expecrhotranslated(k,m0,c0)
    s1 = sigma1t(k,m0,c0)
    s2 = sigma2t(k,m0,c0)  
    }
  
  # Creating list of matrices in linear decomposition of V
  # either with variance random effects equal to Id
  # or from given L-list
  if(is.null(L) ==TRUE) {L = array(data = 0, dim = c(k,k,lZ))
  for(i in 1:lZ) {
    L[,  , i] = Z[[i]] %*% t(Z[[i]])
  }}else{
    LL = array(data = 0, dim = c(k,k,lZ))
    for(i in 1:lZ) {
      LL[,  , i] = L[[i]] 
  }
  L=LL}

  #Setting initial values beta0 and V0 using MCD
  Vstart=covMcd(Ymat)
  if (is.null(V0) == TRUE){
    V0 = Vstart$cov
    }
    w = Vstart$mcd.wt

  V = V0 # initial value for V
  
  ##############################################################################
  #  RUN 1 ITERATION STEP
  ##############################################################################
  q = ncol(X[[1]])  # length of beta
  
  betavecterm = matrix(0,nrow=q,ncol=1)
  betamatterm = matrix(0,nrow=q,ncol=q)
  for(i in 1:n){
    betavecterm = betavecterm+w[i]*t(X[[i]])%*%solve(V)%*%as.vector(Ymat[i,])
    betamatterm = betamatterm+w[i]*t(X[[i]])%*%solve(V)%*%X[[i]]
    }
  beta = solve(betamatterm)%*%betavecterm # initial value beta0
  

  ##############################################################################
  #  ITERATION TO COMPUTE BETAHAT AND THETAHAT
  ##############################################################################
  tol = 1
  iter = 0
  theta=numeric(lZ)
  
  while ((iter <= maxiter) & (tol > eps)) {
    wold = w
    betaold  = beta
    thetaold = theta
    Vold  = V
    ##############################################################################
    # Iteration step
    ##############################################################################
    
    # Re-scaling Mahalanobis distances to satisfy S-constraint
    # Computing array for n Mahalanobis distances
    MD=numeric(n)
    for (i in 1:n){
      y = Ymat[i,]
      mu = X[[i]]%*%betaold
      MD[i] = mahalanobis(y,center = mu,cov = Vold) # Note MD=d^2!
    }
    
    # determining scaling constant for MD
    # for translated biweight
    objfuntranslated=function(s){
            mean(biweightrhotranslated(sqrt(MD)/s,m0,c0))-expecrhotranslated(k,m0,c0)
            }
     MDscaletranslated = uniroot(f=objfuntranslated,c(0.01,100))$root
     MD = MD/MDscaletranslated^2

     # Computing weights with updated MD
     w = sapply(MD,function(MD){biweightutranslated(sqrt(MD),m0,c0)})
     v = sapply(MD,function(MD){biweightu2translated(sqrt(MD),m0,c0)})

 
  ##############################################################################
  # updating betahat
    
    betavecterm = matrix(0,nrow=q,ncol=1)
    betamatterm = matrix(0,nrow=q,ncol=q)
    for(i in 1:n){
      betavecterm = betavecterm + w[i]*t(X[[i]]) %*%solve(Vold)%*%as.vector(Ymat[i,  ])
      betamatterm = betamatterm + w[i]*t(X[[i]])%*%solve(Vold)%*%X[[i]]
    }
    beta = solve(betamatterm)%*%betavecterm   
  
  ##############################################################################
  # updating thetahat
    
    Q = matrix(0,nrow=lZ,ncol=lZ)
    for(i in 1:lZ) {
      for(j in 1:lZ) {
        Q[i, j] = sum(diag(solve(Vold) %*%L[,  , i] %*% solve(Vold) %*%L[ ,  , j]))
      }
    }

    Uterm = matrix(0, nrow = lZ,ncol = n)
    for(i in 1:n) {
        y = as.vector(Ymat[i,])
        for(j in 1:lZ) {
          Uterm[j,i] = Uterm[j,i] + k*w[i] * (t(y- X[[i]]%*%betaold) %*% (solve(Vold) %*% L[,  , j] %*% solve(Vold)) %*% (y- X[[i]]%*%betaold))
      }
    }
    
    theta =(1/sum(v))* (solve(Q) %*% apply(Uterm,1,sum))
    
  ##############################################################################
  # updating V(thetahat)
  
    V = 0
    for (j in 1:lZ) {
      V = V + (theta[j] * L[,  , j])
    }
    
  ##############################################################################
  # check relative difference
    
    diffbeta=max(abs(beta/betaold-1))
    difftheta=max(abs(theta/thetaold-1))
    
    tol = max(diffbeta, difftheta)
    iter = iter + 1
  
  } #END of iteration for betahat and thetahet
  
  # Scaling V one last time to satisfy S-constraint
  MD=numeric(n)
  for (i in 1:n){
    y = Ymat[i,]
    mu = X[[i]]%*%beta
    MD[i] = mahalanobis(y,center = mu,cov = V) # Note MD=d^2!
  }
  
  # determining scaling constant for MD
  objfuntranslated=function(s){
    mean(biweightrhotranslated(sqrt(MD)/s,m0,c0))-expecrhotranslated(k,m0,c0)
  }
  s = uniroot(f=objfuntranslated,c(0.01,100))$root

  # Correction of theta and V
  theta = theta*s^2
  
  V = 0
  for (j in 1:lZ) {
    V = V + (theta[j] * L[,  , j])
  }

  iterS=iter
  
  ##################################################
  # Computing limiting variance and SE of betahat 
  ##################################################

  # Computing estimate for (E[X'V^{-1}X])^{-1}
  termtXX = matrix(0,nrow=length(beta),ncol=length(beta))
  for(i in 1:n){
    termtXX = termtXX+t(X[[i]])%*%solve(V)%*%X[[i]]
  }
  # Limiting variance of sqrt(n)*(betahat-beta)
  varbetaS = (constbetahattranslated(k,m0,c0))*solve(termtXX/n)
  
  # SE betahat
  betaS=beta
  SEbetaS = sqrt(diag(varbetaS)/n)
  tvalS=betaS/SEbetaS
  pvalueS=2*(1-pnorm(abs(tvalS)))
  
  ##################################################
  # Computing limiting variance and SE of thetahat 
  ##################################################
  
  # Computing estimate for (L'(V^{-1}xV^{-1})L)^{-1}
  vecL = matrix(0,k*k,lZ)
  for (j in 1:lZ) {
    vecL[,j] = as.vector(L[,  , j])
  }
  thetaS=theta
  term = solve(t(vecL) %*%(solve(V)%x%solve(V))%*%vecL)

  # Limiting variance of sqrt(n)*(thetahat-theta)
  varthetaS = 2 * s1 *term  + s2 * thetaS%*%t(thetaS)
  
  # SE thetahat
  SEthetaS = sqrt(diag(varthetaS)/n)
  tvalthetaS=thetaS/SEthetaS
  pvaluethetaS=2*(1-pnorm(abs(tvalthetaS)))
  
  ####################################
  # Computation of MM-estimator
  ####################################

  iterM=NULL
  
  # Determine cut-off constant rhoMM for efficiency eff
  if(is.null(rhoMM) ==FALSE) {
    objectf=function(c){
      1/constbetahat(k,c)-eff 
    }
    c1=uniroot(objectf,interval=c(0.001,100))$root

  # Computing expectation of rhoMM
  if (rhoMM =="biweight"){
      m1=0
      b0=expecrho(k,c1)
    }
    
    
    tol = 100
    iter = 0
    beta = betaS  # starting from S-estimator betahat

  ##############################################################################
    while((iter <= maxiter) & (tol > eps)) {

      betaold  = beta

      MD=numeric(n)
      for (i in 1:n){
        y=Ymat[i,]
        mu = X[[i]]%*%beta
        MD[i]=mahalanobis(y,center=mu,cov=V) # Note MD=d^2!
      }
      objfuntranslated=function(s){
        mean(biweightrhotranslated(sqrt(MD)/s,m1,c1))-expecrhotranslated(k,m1,c1)
      }
      w = sapply(MD,function(MD){biweightutranslated(sqrt(MD),m1,c1)})
      v = sapply(MD,function(MD){biweightu2translated(sqrt(MD),m1,c1)})
      
      betavecterm = matrix(0,nrow=q,ncol=1)
      betamatterm = matrix(0,nrow=q,ncol=q)
      for(i in 1:n){
        betavecterm = betavecterm + w[i]*t(X[[i]]) %*%solve(V)%*%as.vector(Ymat[i,  ])
        betamatterm = betamatterm + w[i]*t(X[[i]])%*%solve(V)%*%X[[i]]
      }
      
      beta = solve(betamatterm)%*%betavecterm  # MM-estimator for beta
      
  ##############################################################################
  
   tol=max(abs(beta/betaold-1))
      iter = iter + 1
    }
    iterM=iter
    
    # Computation of limiting variance and SE for MM
    termtXX = matrix(0,nrow=length(beta),ncol=length(beta))
    for(i in 1:n){
      termtXX = termtXX+t(X[[i]])%*%solve(V)%*%X[[i]]
    }
    betaMM=beta
    varbetaMM = (constbetahattranslated(k,m1,c1))*solve(termtXX/n)
    SEbetaMM = sqrt(diag(varbetaMM)/n)
  }else{
    betaMM = NA
    varbetaMM = NA
    SEbetaMM = NA
  }
  
  tvalMM=betaMM/SEbetaMM
  pvalueMM=2*(1-pnorm(abs(tvalMM)))
  
  
  # Summaries for fixed effects with estimates SE t-val p-value
  fixedeffectsS=cbind(betaS,SEbetaS,tvalS,pvalueS)
  fixedeffectsMM=cbind(betaMM,SEbetaMM,tvalMM,pvalueMM)
  
  summarythetaS = cbind(thetaS,SEthetaS,tvalthetaS,pvaluethetaS)
  
  fixedeffectsS = data.frame(fixedeffectsS)
  fixedeffectsMM = data.frame(fixedeffectsMM)
  colnames(fixedeffectsS) = c("beta","SEbeta","tval","p-value")
  colnames(fixedeffectsMM) = c("beta","SEbeta","tval","p-value")
  
  return(list(fixedeffectsS=fixedeffectsS,
              fixedeffectsMM=fixedeffectsMM,
              summarythetaS=summarythetaS,
              varbetaShat=varbetaS/n,
              varbetaMMhat=varbetaMM/n,
              varthetahat=varthetaS/n,
              w=w,dis=MD,iterS=iterS,iterM=iterM))
}