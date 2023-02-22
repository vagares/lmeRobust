filepath="C:\\Users\\vagares\\Documents\\lmeRobust\\"
source(paste(filepath,"biweight_functions.R",sep=""))
source(paste(filepath,"asympt_norm_constants.R",sep=""))

Roblme = function(Ymat,X,Z,E=NULL,L=NULL,rho ="t-biweight",r =0.5,arp=0.01,rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL){
  # Ymat: outcome 
  # X: X in a list X[[i]] i = 1...n
  # Z: Z in a list Z[[z]]
  # rho: rho function for the S estimation (t-biweight, biweight or MLE)
  # r: breakdown point
  # arp: asymptotic rejection rate
  # rhoMM: rho function for the MM estimation (biweight)
  # eps: bound for the stopping criteria
  # maxiter: maximum number of iteration
  # eff: efficiency
  # mu0: initial mean value
  # S0: initial covariance value
  lZ = length(Z)
  lX = length(X)
  k = ncol(Ymat)
  n = nrow(Ymat)
  if(is.null(L) ==TRUE) {
  if(is.null(E) ==TRUE) { Z[[lZ+1]] =  diag(rep(1,k)) # add identity matrix for error
  }else{Z[[lZ+1]] =  E}
  lZ = length(Z)}else{lZ = length(L)}
  # Setting the breakdown point and cut-off constant for translated biweight
  if (rho =="t-biweight"){
    c0 = rhotranslatedconst(k,r,arp,0.01,10)
    m0 = sqrt(qchisq(1-arp,df=k))-c0
    b0 = expecrhotranslated(k,m0,c0)
    s1 = sigma1t(k,m0,c0)
    s2 = sigma2t(k,m0,c0) 
    }
  if (rho =="biweight"){
    m0 = 0
    c0 = rhoconst(k,r,0.01,100)
    b0 = expecrho(k,c0)
    s1 = sigma1(k,c0)
    s2 = sigma2(k,c0)
    }
  if (rho =="MLE"){
    m0 = 10000
    c0 = rhotranslatedconst(k,r,arp,0.01,10)
    b0 = expecrhotranslated(k,m0,c0)
    s1 = sigma1t(k,m0,c0)
    s2 = sigma2t(k,m0,c0)  
    }
  

  if(is.null(L) ==TRUE) {L = array(data = 0, dim = c(k,k,lZ))
  for(i in 1:lZ) {
    L[,  , i] = Z[[i]] %*% t(Z[[i]])
  }}else{
    lZ=4
    LL = array(data = 0, dim = c(k,k,lZ))
    for(i in 1:lZ) {
      LL[,  , i] = L[[i]] 
    }
  L=LL}
  
  
  #INITIAL VALUES TO CHANGE mu0, S0 ==> MCD weight and first stape
  Vstart=covMcd(Ymat)
  if (is.null(V0) == TRUE){
    V0 = Vstart$cov
    }
    w = Vstart$mcd.wt

  V = V0
  ##############################################################################
  #  RUN 1 ITERATION STEP
  ##############################################################################
  q = ncol(X[[1]])
  betavecterm = matrix(0,nrow=q,ncol=1)
  betamatterm = matrix(0,nrow=q,ncol=q)
  for(i in 1:n){
    betavecterm = betavecterm+w[i]*t(X[[i]])%*%solve(V)%*%as.vector(Ymat[i,])
    betamatterm = betamatterm+w[i]*t(X[[i]])%*%solve(V)%*%X[[i]]
    }
  beta = solve(betamatterm)%*%betavecterm
  
  tol = 1
  iter = 0
  while ((iter <= maxiter) & (tol > eps)) {
    wold = w
    betaold  = beta
    Vold  = V
    ##############################################################################
    # Iteration step
    ##############################################################################
    
    # Re-scaling Mahalanobis distances to satisfy S-constraint
    # array for n Mahalanobis distances
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

     w = sapply(MD,function(MD){biweightutranslated(sqrt(MD),m0,c0)})
     v = sapply(MD,function(MD){biweightu2translated(sqrt(MD),m0,c0)})

 
  ##############################################################################
    betavecterm = matrix(0,nrow=q,ncol=1)
    betamatterm = matrix(0,nrow=q,ncol=q)
    for(i in 1:n){
      betavecterm = betavecterm + w[i]*t(X[[i]]) %*%solve(Vold)%*%as.vector(Ymat[i,  ])
      betamatterm = betamatterm + w[i]*t(X[[i]])%*%solve(Vold)%*%X[[i]]
    }
    beta = solve(betamatterm)%*%betavecterm
    ##############################################################################
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
    ##############################################################################
    diffbeta = norm(beta - betaold,type="F")
    difftheta = norm(V - Vold,type="F")
    tol = max(diffbeta, difftheta)
    iter = iter + 1
  }
  
  termtXX = matrix(0,nrow=length(beta),ncol=length(beta))
  for(i in 1:n){
    termtXX = termtXX+t(X[[i]])%*%solve(V)%*%X[[i]]
  }
  betaS=beta
  varbetaS = (constbetahattranslated(k,m0,c0))*solve(termtXX/n)
  
  vecL = matrix(0,k*k,lZ)
  for (j in 1:lZ) {
    vecL[,j] = as.vector(L[,  , j])
  }
  thetaS=theta
  term = solve(t(vecL) %*%(solve(V)%x%solve(V))%*%vecL)
  varthetaS = 2 * s1 *term  + s2 * thetaS%*%t(thetaS)
  SEthetaS = sqrt(diag(varthetaS)/n)
  tvalthetaS=thetaS/SEthetaS
  pvaluethetaS=2*(1-pnorm(abs(tvalthetaS)))
  
  SEbetaS = sqrt(diag(varbetaS)/n)
  tvalS=betaS/SEbetaS
  pvalueS=2*(1-pnorm(abs(tvalS)))
  
  if(is.null(rhoMM) ==FALSE) {
    objectf=function(c){
      1/constbetahat(k,c)-eff 
    }
    c1=uniroot(objectf,interval=c(0.001,100))$root
    

    #if (rhoMM =="t-biweight"){
    #  m1=sqrt(qchisq(1-arp,df=k))-c1
    #  b0=expecrhotranslated(k,m1,c1)
    #}
    if (rhoMM =="biweight"){
      m1=0
      b0=expecrho(k,c1)
    }
    
    
    tol = 100
    iter = 0
    beta = betaS
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
      beta = solve(betamatterm)%*%betavecterm
      
      ##############################################################################
      tol = norm(beta - betaold,type="F")
      iter = iter + 1
    }
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
  
  
  #fixedeffects with estimates SE t-val p-value
  fixedeffectsS=cbind(betaS,SEbetaS,tvalS,pvalueS)
  fixedeffectsMM=cbind(betaMM,SEbetaMM,tvalMM,pvalueMM)
  
  summarythetaS = cbind(theta,SEthetaS,tvalthetaS,pvaluethetaS)
  
  fixedeffectsS = data.frame(fixedeffectsS)
  fixedeffectsMM = data.frame(fixedeffectsMM)
  colnames(fixedeffectsS) = c("beta","SEbeta","tval","p-value")
  colnames(fixedeffectsMM) = c("beta","SEbeta","tval","p-value")
  return(list(fixedeffectsS=fixedeffectsS,fixedeffectsMM=fixedeffectsMM,summarythetaS=summarythetaS,w=w,dis=MD))
}