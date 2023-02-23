library(mvtnorm) 

#########################################################################

data_sphe = function(x, tol = 1e-06){
  Nm = nrow(x) - 1
  MEANS = colMeans(x)
  x.c = as.matrix(sweep(x, 2, MEANS, "-"))
  SVD = svd(x.c, nu = 0)
  SV = SVD$d
  if (!is.null(tol)) {
    rank = sum(SVD$d > (SVD$d[1L] * tol))
    if (rank < ncol(x)) {
      SVD$v = SVD$v[, 1L:rank, drop = FALSE]
      SVD$d = SVD$d[1L:rank]
    }
  }
  SIGMAs = SVD$d/sqrt(Nm+1)
  TRANS = SVD$v *(1/SIGMAs)
  # TRANS <- SVD$v %*% diag(1/SIGMAs)
  RES = x.c %*% TRANS
  #attr(RES, "center") <- MEANS
  #attr(RES, "transform") <- TRANS
  #attr(RES, "backtransform") <- diag(SIGMAs) %*% t(SVD$v)
  # attr(RES, "SingularValues") <- SV
  RES
}



data_gen = function(n=1000,k=2,cont,p,mue=0,se=1,mug=0,sg=1){
  #cont = 0 no contamination
  #cont = 1 contamination on gamma with mug and sg
  #cont = 2 contamination on epsilon with mue and se
  #cont = 3 contamination on gamma and epsilon
  #cont = 4 contamination on X
  # proportion of contamination
  N = k * n # total number of observations
  no = rep(k, n)
  id = numeric(0)
  for (i in (1:n)){id = c(id, rep(i,k))}
  obs = rep(seq(1,k), n)
  y = numeric(N)
  X = numeric(N)
  X = rnorm(N,0,1)
  XX = cbind(rep(1,N), X)
  
  Z=NULL
  for (i in (1:n)){
    Xtemp = XX[((i-1)*k+1):((i-1)*k+k),]
    Ztemp = data_sphe(Xtemp)
    Z = rbind(Z,Ztemp)
  }
  
  X =   cbind(rep(1,N), Z)
  #t(X)%*%X #20 = 4*5
  for (i in (1:n)){
    if (cont %in% c(0,2) ){gamma = rnorm(1,0,1)} else{ if (cont %in% c(1,3) ) {gamma = ifelse(runif(1) <= (1-p),rnorm(1,0,1),rnorm(1,mug,sqrt(sg)))}}
    for (j in (1:no[i])){
      if (cont %in% c(0,1)){eps = rnorm(1,0,1)} else{if (cont %in% c(2,3) ) { eps = ifelse(runif(1) <= (1-p),rnorm(1,0,1),rnorm(1,mue,sqrt(se)))}}
      if (cont %in% c(0,1,2,3)){
        y[id == i & obs == j] =  1 + X[id == i & obs == j][2] + gamma + eps}
    }}
  data = data.frame(id = id,obs = obs,y = y,x = X[,2])
  colnames(data) = c("id", "obs", "y", "x")
  return(data)
}


data_gen = function(n=1000,k=2,pe,pb,mue=0,mub=0,px,alpha){
  id =rep(1,n)
  X = list()
  m = rep(0,k)
  s = diag(k)
  fc = rep(1,k)
  sc = seq(0,k-1,1)
  Z =  list(cbind(fc,sc))
  for (i in (1:n)){
    
    X[[i]] <- cbind(fc,sc)
  }
  Y = list()  

  mg0 = rep(0,2)
  sg0 = matrix(c(790,-8.5,-8.5,40),2,2)
  mgc = c(0,mg)
  me0 = 0
  se0 = 1
  nobi=0
  noe=0
  nox=0
  noei=0
  noxi=0
  for (i in (1:n)){
    noeitemp=0
    noxitemp=0
    if (runif(1) <= (1-pb)){b= rmvnorm(1, mean=mg0, sigma=sg0)}
    else{b=rmvnorm(1, mean=mgc, sigma=sg0);nobi=nobi+1}
    eps = numeric(k)
    for (kk in (1:k)){
      if (runif(1) <= (1-pe)){eps[kk] = rnorm(1, mean=me0, sd=se0)}
      else{eps[kk] = rnorm(1, mean=me, sd=0.5);noe=noe+1;noeitemp=1}
    }
    for (r in (1:k)){
      if(runif(1) <= (1-px)){ X[[i]][r,2]= X[[i]][r,2]}
      else{X[[i]][r,2]= alpha*X[[i]][r,2];nox=nox+1;noxitemp=1}
      }
  beta = c(250,10)
  Y[[i]] =  X[[i]]%*%beta + Z[[1]]%*%as.vector(b) + eps
  noei = noei + noeitemp
  noxi = noxi + noxitemp
  
  }
  Ymat=matrix(0,n,k)
  for (i in (1:n)){Ymat[i,] = Y[[i]]}
  
  return(list(Y=Ymat,X=X,Z=Z,nobi=nobi,noe=noe,nox=nox,noei=noei,noxi=noxi))
}

dat = data_gen(n=100,k=4,pe=0,pb=0,mue=0,mub=0,px=0.05,alpha=1)
#dat

summaryRoblme = Roblme(dat$Y,dat$X,dat$Z,E=NULL,L=NULL,rho ="biweight",r =0.5,arp=0.01,rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL)
  