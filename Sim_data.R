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
