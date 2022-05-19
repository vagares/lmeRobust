Roblme = function(Ymat,X,Z,rho ="t-biweight",r =0.5,arp=0.01,rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,cov0="MCD"){
  # Ymat: outcome 
  # X: X in a list X[[i]] i = 1...n
  # Z: Z in a list Z[[z]]
  # rho: rho function for the S estimation (t-biweight, biweight or MLE)
  # r: breakdown point
  # arp: asymptotic rejection rate
  # rhoMM: rho function for the MM estimation (NULL, biweight)
  # eps: bound for the stopping criteria
  # maxiter: maximum number of iteration
  # eff: efficiency
  
  # cov0: initial covariance and weights (MCD from the robustbase package)

  ##############################################################################
  #  INITIALIZATION  
  ##############################################################################
  # Dimensions
  lZ = length(Z)
  lX = length(X)
  k = ncol(Ymat)
  n = nrow(Ymat)

  # Setting the breakdown point and cut-off constant for translated biweight
  if (rho =="t-biweight"){
    c0 = rhotranslatedconst(k,r,arp,0.01,10)
    m0 = sqrt(qchisq(1-arp,df=k))-c0
    b0 = expecrhotranslated(k,m0,c0)
        }
  if (rho =="biweight"){
    m0 = 0
    c0 = rhoconst(k,r,0.01,100)
    b0 = expecrho(k,c0)
    }
  if (rho =="MLE"){
    m0 = 10000
    c0 = rhotranslatedconst(k,r,arp,0.01,10)
    b0 = expecrhotranslated(k,m0,c0)
    }
  
  # L matrices
  L = array(data = 0, dim = c(k,k,lZ))
  for(i in 1:lZ) {
    L[,  , i] = Z[[i]] %*% t(Z[[i]])
  }
  
  
  # Initial values for the covariance and weights
  
  if (cov0 == "MCD"){
    Vstart=covMcd(Ymat)
    V0 = Vstart$cov
    
    w = Vstart$mcd.wt # Initial weights
  } else {
  print("Error: the input has to be MCD")
  }
  V = V0  # Initial covariance matrix
  ##############################################################################
  #  S-ESTIMATORS COMPUTATION
  ##############################################################################
  #####################################
  #First iteration step for S-estimators
  #####################################
  q = ncol(X[[1]])
  betavecterm = matrix(0,nrow=q,ncol=1)
  betamatterm = matrix(0,nrow=q,ncol=q)
  for(i in 1:n){
    betavecterm = betavecterm+w[i]*t(X[[i]])%*%solve(V)%*%as.vector(Ymat[i,])
    betamatterm = betamatterm+w[i]*t(X[[i]])%*%solve(V)%*%X[[i]]
  }
  beta = solve(betamatterm)%*%betavecterm # Initial Beta
  
  tol = 1
  iter = 0
  ##########################################################
    # Iteration steps for S-estimators    
  ##########################################################
    

  while ((iter <= maxiter) & (tol > eps)) {
    wold = w
    betaold  = beta
    Vold  = V
    
    ##########################################################
    # Re-scaling Mahalanobis distances to satisfy S-constraint
    ##########################################################
        MD=numeric(n) # array for n Mahalanobis distances    
    for (i in 1:n){
      y = Ymat[i,]
      mu = X[[i]]%*%betaold
      MD[i] = mahalanobis(y,center = mu,cov = Vold) # Note MD=d^2!
    }
    
    # Mahalanobis distance standardization     
    objfuntranslated=function(s){
            mean(biweightrhotranslated(sqrt(MD)/s,m0,c0))-expecrhotranslated(k,m0,c0)
            }
     MDscaletranslated = uniroot(f=objfuntranslated,c(0.01,100))$root # determining scaling constant for MDfor translated biweight    
     MD = MD/MDscaletranslated^2

     w = sapply(MD,function(MD){biweightutranslated(sqrt(MD),m0,c0)})# Computation of the translated biweight function u(s)s^2=psi(s)s at the Mahalanobis distance
     v = sapply(MD,function(MD){biweightu2translated(sqrt(MD),m0,c0)})#
    # Computation of the translated biweight function u(s)=psi(s)/s at the Mahalanobis distance
 
  
    ####################
    # Beta COMPUTATION
    ####################
    betavecterm = matrix(0,nrow=q,ncol=1)
    betamatterm = matrix(0,nrow=q,ncol=q)
    for(i in 1:n){
      betavecterm = betavecterm + w[i]*t(X[[i]]) %*%solve(Vold)%*%as.vector(Ymat[i,  ])
      betamatterm = betamatterm + w[i]*t(X[[i]])%*%solve(Vold)%*%X[[i]]
    }
    beta = solve(betamatterm)%*%betavecterm
    
    #########################
    # Theta and V COMPUTATION
    #########################
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
    V = 0
    for (j in 1:lZ) {
      V = V + (theta[j] * L[,  , j])
    }
    
    #############################
    # Checking stopping condition
    #############################
    diffbeta = norm(beta - betaold,type="F")
    difftheta = norm(V - Vold,type="F")
    tol = max(diffbeta, difftheta)
    iter = iter + 1
  } # end of the iterations for the S-estimators computation
  

  #########################################
  # Standard error of Beta for S-estimators
  #########################################
  termtXX = matrix(0,nrow=length(beta),ncol=length(beta))
  for(i in 1:n){
    termtXX = termtXX+t(X[[i]])%*%solve(V)%*%X[[i]]
  }
  varbeta = (constbetahattranslated(k,m0,c0))*solve(termtXX/n)
  SEbeta = sqrt(diag(varbeta)/n)
  
##### END of the COMPUTATIONS of the S-ESTIMATORS ####

##############################################################################
#  COMPUTATION of the MM-ESTIMATORS (if chosen)
##############################################################################
  if(is.null(rhoMM) ==FALSE) {
  if (rhoMM =="biweight"){    # Constants computation
    objectf=function(c){
      1/constbetahat(k,c)-eff 
    }
    c1=uniroot(objectf,interval=c(0.001,100))$root
    m1=0
    b0=expecrho(k,c1)
    
    tol = 100
    iter = 0
    ##########################################################    
    # Iteration steps for MM-estimators  
    ##########################################################                       
    while((iter <= maxiter) & (tol > eps)) {

      betaold  = beta

      # Computation of the squared Mahalanobis distance and the weights
      MD=numeric(n)
      for (i in 1:n){
        y=Ymat[i,]
        mu = X[[i]]%*%beta
        MD[i]=mahalanobis(y,center=mu,cov=V) # Note MD=d^2!
      }
      
      w = sapply(MD,function(MD){biweightutranslated(sqrt(MD),m1,c1)})
      v = sapply(MD,function(MD){biweightu2translated(sqrt(MD),m1,c1)})
      
            # Beta computation
      betavecterm = matrix(0,nrow=q,ncol=1)
      betamatterm = matrix(0,nrow=q,ncol=q)
      for(i in 1:n){
        betavecterm = betavecterm + w[i]*t(X[[i]]) %*%solve(V)%*%as.vector(Ymat[i,  ])
        betamatterm = betamatterm + w[i]*t(X[[i]])%*%solve(V)%*%X[[i]]
      }
      beta = solve(betamatterm)%*%betavecterm
      
      
      # Stopping criterion
      tol = norm(beta - betaold,type="F")
      iter = iter + 1
    }# End of the iterations for the MM computation

    ##########################################
    # Standard error of Beta for MM-estimators
    ##########################################
    termtXX = matrix(0,nrow=length(beta),ncol=length(beta))
    for(i in 1:n){
      termtXX = termtXX+t(X[[i]])%*%solve(V)%*%X[[i]]
    }
    varbeta = (constbetahattranslated(k,m1,c1))*solve(termtXX/n)
    SEbeta = sqrt(diag(varbeta)/n)
    }} # End of the MM computations
  
  
  tval=beta/SEbeta # Test statistics
  pvalue=2*(1-pnorm(abs(tval))) # p-value
  
  
  # Fixedeffects with estimates SE t-val p-value
  fixedeffects=cbind(beta,SEbeta,tval,pvalue)
  fixedeffects = data.frame(fixedeffects)
  colnames(fixedeffects) = c("beta","SEbeta","tval","p-value")
  return(list(fixedeffects=fixedeffects,theta=theta,w=w,dis=MD))
  }
