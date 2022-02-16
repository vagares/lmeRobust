filepath="C:\\Users\\vagares\\Documents\\lmeRobust\\"
source(paste(filepath,"utils_functions.R",sep=""))

Roblme = function(y,xfit,zfit,rho ="t-biweight",r =0.25,arp=0.01,rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,mu0=NULL,s0=NULL){
  # y: outcome 
  # xfit: X
  # zfit: Z
  # rho: rho function for the S estimation
  # r: breakdown point
  # arp: asymptotic rejection rate
  # rhoMM: rho function for the MM estimation
  # eps:
  # maxiter: maximum number of iteration
  # eff: efficiency
  # mu0: initial mean value
  # S0: initial covariance value
  zlength <- length(zfit)
  lxfit <- length(xfit)
  p <- ncol(y)
  n <- nrow(y)
  #TO CHANGE
  b0=2.312195
  c0 = Tbsc(r,p)
  if (rho =="t-biweight"){m0=0
          ##TO CHANGE
        }
  if (rho =="biweight"){m0=0}
  if (rho =="MLE"){m0=1000}
  

  qarray <- array(data = 0, dim = c(p,p,zlength))
  marray <- qarray
  for(i in 1:zlength) {
    qarray[,  , i] <- zfit[[i]] %*% t(zfit[[i]])
  }
  
  
  #INITIAL VALUES TO CHANGE mu0, S0 ==> MCD weight and first stape
  initial <- rogkmiss(y)
  if(is.null(mu0) || is.null(s0)) {
    if(is.null(mu0))
      mu0 <- initial$center
    if(is.null(s0))
      s0 <- initial$cov
  }
  w <- initial$w
  
  mats = s0
  inv <- solve(mats)
  ##############################################################################
  #  X matrices and XX' matrices
  ##############################################################################
  k <- ncol(xfit[[1]])
  sx <- matrix(0,k,1)
  tx <- matrix(0,k,k)
  for(i in 1:n){
    sx <- sx+w[i]*t(xfit[[i]])%*%inv%*%as.vector(y[i,])
    tx <- tx+w[i]*t(xfit[[i]])%*%inv%*%xfit[[i]]
  }
  alpha <- solve(tx)%*%sx
  
  crit <- 100
  iter <- 1
  w1d <- rep(1, n)
  d2 <- rep(0,n)
  w2d <- w1d
  vec <- rep(0,zlength)
  Q <- matrix(0,zlength,zlength)
  marray <- qarray
  
  while((iter <= maxiter) & (crit > eps)) {
    inv <- solve(mats)
    wt.old <- w1d
    alpha.old <- alpha
    mats.old <- mats
    for(i in 1:n){
      mu <- xfit[[i]]%*%alpha
      d2[i] <- mahalanobis(y[i,],mu,mats)
    }
    d <- sqrt(d2)
    #h <- floor((n + p + 1)/2)
    #quantile <- h/(n + 1)
    #d <- (d * sqrt(qchisq(quantile, p)))/(sort(d)[h])
    fun=function(kk){mean(sapply(d,function(d)biweightrhotranslated(d/kk,m0,c0)))-b0}
    kk = uniroot(f = fun, c(0.01, 10))$root
    d = d/kk
    w1d <- sapply(d,function(d){biweightu2translated(d,m0,c0)})
    w2d <- sapply(d,function(d){d*biweightpsitranslated(d,m0,c0)})
    ##############################################################################
    sx <- matrix(0,k,1)
    tx <- matrix(0,k,k)
    for(i in 1:n){
      sx <- sx + w1d[i]*t(xfit[[i]]) %*%inv%*%as.vector(y[i,  ])
      tx <- tx + w1d[i]*t(xfit[[i]])%*%inv%*%xfit[[i]]
    }
    alpha <- solve(tx)%*%sx
    ##############################################################################
    for(i in 1:zlength) {
      marray[,  , i] <- inv %*% qarray[,  , i]
    }
    for(i in 1:zlength) {
      for(j in 1:zlength) {
        Q[i, j] <- sum(diag(marray[,  , i] %*% marray[
          ,  , j]))
      }
    }
    Qinv <- solve(Q)
    vec <- vec * 0
    UT <- matrix(0, zlength, n)
    for(i in 1:zlength) {
      for(j in 1:n) {
        t1 <- xfit[[j]]%*%alpha
        xc <- as.vector(y[j,  ] - t1)
        vec[i] <- vec[i] + k*w1d[j] * (t(xc) %*% (inv %*%
                                                    qarray[,  , i] %*% inv) %*% xc)
        UT[i, j] <- t(xc) %*% (inv %*% qarray[, , i] %*% inv) %*% xc
      }
    }
    s <- (solve(Q) %*% vec)/sum(w2d)
    mats <- mats * 0
    for(i in 1:zlength) {
      mats <- mats + (s[i] * qarray[,  , i])
    }
    ##############################################################################
    crit1 <- max(abs(alpha - alpha.old))
    crit2 <- max(abs(mats - mats.old))
    crit <- max(crit1, crit2)
    iter <- iter + 1
  }
  if(is.null(rhoMM) ==FALSE) {
    if (rhoMM =="t-biweight"){m1=0;c1=Tbsc1(eff,p-2)}
    if (rhoMM =="biweight"){m1=0;c1=Tbsc1(eff,p-2)}
    crit <- 100
    iter <- 1
    w1d <- rep(1, n)
    d2 <- rep(0,n)
    w2d <- w1d
    vec <- rep(0,zlength)
    Q <- matrix(0,zlength,zlength)
    ##############################################################################
    while((iter <= maxiter) & (crit > eps)) {
      inv <- solve(mats)
      wt.old <- w1d
      alpha.old <- alpha
      mats.old <- mats
      for(i in 1:n){
        mu <- xfit[[i]]%*%alpha
        d2[i] <- mahalanobis(y[i,],mu,mats)
      }
      #d <- sqrt(d2)
      #h <- floor((n + p + 1)/2)
      #quantile <- h/(n + 1)
      #d <- (d * sqrt(qchisq(quantile, p)))/(sort(d)[h])
      w1d <- biweightu2translated(d,m1,c1)
      w2d <- d*biweightpsitranslated(d,m1,c1)
      ##############################################################################
      sx <- matrix(0,k,1)
      tx <- matrix(0,k,k)
      for(i in 1:n){
        sx <- sx + w1d[i]*t(xfit[[i]]) %*%inv%*%as.vector(y[i,  ])
        tx <- tx + w1d[i]*t(xfit[[i]])%*%inv%*%xfit[[i]]
      }
      alpha <- solve(tx)%*%sx
      mats <- mats
      ##############################################################################
      crit1 <- max(abs(alpha - alpha.old))
      crit2 <- max(abs(mats - mats.old))
      crit <- max(crit1, crit2)
      iter <- iter + 1
    }
  }
  
  e.SE=1.0521
  var.mu.R <- e.SE*mats/n
  S.xiT.xi <- matrix(0,length(alpha),length(alpha))
  S.xiT.Sig.xi <- matrix(0,length(alpha),length(alpha))
  for(i in 1:n){
    S.xiT.xi = S.xiT.xi+t(xfit[[i]])%*%xfit[[i]]
    S.xiT.Sig.xi = S.xiT.Sig.xi+n*t(xfit[[i]])%*%var.mu.R%*%xfit[[i]]
  }
  var.alpha.R = solve(S.xiT.xi)%*%S.xiT.Sig.xi%*%solve(S.xiT.xi)
  SE.alpha.R = sqrt(diag(var.alpha.R))

  
  
  #fixedeffects with estimates SE t-val p-value
  fixedeffects=cbind(alpha,SE.alpha.R)
  return(list(fixedeffects=fixedeffects,s=s,w=w1d,dis=d2))
}