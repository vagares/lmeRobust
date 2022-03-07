################################################################################
#Biweight and translated biweight functions
################################################################################

biweightrho=function(y,c0){
  # This is Tukey's biweight rho function
  ifelse(abs(y)<=c0,y^2/2-y^4/(2*c0^2)+y^6/(6*c0^4),c0^2/6)
}

biweightrhotranslated=function(d,M,c0){
  # This is Rocke's translated biweight rho function
  term1=d^2/2
  term2=M^2/2-
    M^2*(M^4-5*M^2*c0^2+15*c0^4)/(30*c0^4)+
    d^2*(1/2+M^4/(2*c0^4)-M^2/c0^2)+
    d^3*(4*M/(3*c0^2)-4*M^3/(3*c0^4))+
    d^4*(3*M^2/(2*c0^4)-1/(2*c0^2))-
    4*M*d^5/(5*c0^4)+
    d^6/(6*c0^4)
  term3=M^2/2+c0*(5*c0+16*M)/30
  ifelse(d<M,1,0)*term1+
    ifelse(d>=M,1,0)*ifelse(d<=M+c0,1,0)*term2+
    ifelse(d>M+c0,1,0)*term3
}

#############################################

biweightpsi=function(y,c0){
  # This is Tukey's biweight psi function
  ifelse(abs(y)<=c0,y-2*y^3/(c0^2)+y^5/(c0^4),0)
}

biweightpsitranslated=function(d,M,c0){
  # This is Rocke's translated biweight psi function
  ifelse(d<M,1,0)*d+
    ifelse(d>=M,1,0)*ifelse(d<=M+c0,1,0)*
    d*(1-((d-M)/c0)^2)^2
}

#############################################

biweightu=function(y,c0){
  # This is Tukey's biweight weight function u(s)=psi(s)/s
  ifelse(abs(y)<=c0,1-2*y^2/(c0^2)+y^4/(c0^4),0)
}

biweightutranslated=function(d,M,c0){
  # This is Rocke's translated biweight function u(s)=psi(s)/s
  ifelse(d<M,1,0)+
    ifelse(d>=M,1,0)*ifelse(d<=M+c0,1,0)*
    (1-((d-M)/c0)^2)^2
}

############################################################################


biweightu2=function(y,c0){
  # This is Tukey's biweight weight function u(s)s^2=psi(s)s
  ifelse(abs(y)<=c0,y^2-2*y^4/(c0^2)+y^6/(c0^4),0)
}

biweightu2translated=function(d,M,c0){
  # This is Rocke's translated biweight function u(s)s^2=psi(s)s
  ifelse(d<M,1,0)*d^2+
    ifelse(d>=M,1,0)*ifelse(d<=M+c0,1,0)*
    d^2*(1-((d-M)/c0)^2)^2
}

############################################################################


biweightv=function(y,c0,b0){
  # This is Tukey's biweight weight function v(s)=u(s)s^2-rho(s)+b0
  biweightu2(y,c0)-biweightrho(y,c0)+b0
}

biweightvtranslated=function(d,M,c0,b0){
  # This is Rocke's biweight weight function v(s)=u(s)s^2-rho(s)+b0
  biweightu2translated(d,M,c0)-biweightrhotranslated(d,M,c0)+b0
}

############################################################################
############################################################################

Tbsc <- function(bdp,p,maxit = 1e3,eps = 1e-8,diff = 1e6)
{
  # constant for sfirst Tukey Biweight rho-function for MM, for fixed
  # breakdown point
  talpha = sqrt(qchisq(1-bdp,p))
  ctest = talpha
  iter = 1
  while ((diff>eps)*(iter<maxit))
  {
    cold = ctest
    ctest = Tbsb(cold,p)/bdp
    diff = abs(cold-ctest)
    iter = iter+1
  }
  return(ctest)
}

Tbsb <- function(c,p)
{
  y1 = ksiint(c,1,p)*3/c-ksiint(c,2,p)*3/(c^3)+ksiint(c,3,p)/(c^5)
  y2 = c*(1-pchisq(c^2,p))
  return(y1+y2)
}
################################################################################

ksiint <- function(c,s,p) {(2^s)*gamma(s+p/2)*pgamma(c^2/2,s+p/2)/gamma(p/2)}

Tbsb <- function(c,p)
{
  y1 = ksiint(c,1,p)*3/c-ksiint(c,2,p)*3/(c^3)+ksiint(c,3,p)/(c^5)
  y2 = c*(1-pchisq(c^2,p))
  return(y1+y2)
}

################################################################################
Tbsc1 <- function(eff,p)
{# constant for second Tukey Biweight rho-function for MM, for fixed shape-efficiency
  #
  # Direct cals: sigma1, Tbsc
  
  maxit <- 1000
  eps <- 10^(-8)
  diff <- 10^6
  ctest <- Tbsc(0.5,p)
  iter <- 1
  while ((diff>eps) & (iter<maxit)){
    cold <- ctest
    ctest <- cold*eff*sigma1(cold,p)
    diff <- abs(cold-ctest)
    iter <- iter+1
    
  }
  return(ctest)
}

################################################################################
gint<-function(k,c,p)
{
  # this procedures computes the integral from zero to c of
  # the function r^k g(r^2), where g(||x||^2) is the density function
  # of a p-dimensional standardnormal distribution
  #
  # Direct calls: none
  e<-(k-p-1)/2
  numerator<-(2^e)*pgamma((c^2)/2,(k+1)/2)*gamma((k+1)/2)
  res<-(numerator/(pi^(p/2)))
  return(res)
}

################################################################################
#Constant computation
################################################################################
sigma1<-function(c,p)
  # Direct calls: gint
{
  Cp <- 2*(pi^(p/2))/gamma(p/2)
  gamma1.1 <- gint(p+1,c,p) - 6*gint(p+3,c,p)/(c^2) + 5*gint(p+5,c,p)/(c^4)
  gamma1.2 <- gint(p+1,c,p) - 2*gint(p+3,c,p)/(c^2) + gint(p+5,c,p)/(c^4)
  gamma1 <- Cp * ( gamma1.1 + (p+1)*gamma1.2 ) / (p+2)
  
  s1 <- Cp * ( gint(p+3,c,p) - 4*gint(p+5,c,p)/(c^2) + 6*gint(p+7,c,p)/(c^4) -
                 4*gint(p+9,c,p)/(c^6) + gint(p+11,c,p)/(c^8) ) /gamma1^2 * p/(p+2)
  return(s1)
}

"gnan0" = function(x, c1,c2)
{
  d <- apply(x, 2, mad)
  w <- sweep(x, 2, apply(x, 2, median), FUN = "-")
  w.c <- sweep(w, 2, d, FUN = "/")
  I <- ifelse(abs(w.c) <= c1, 1, 0)
  wc <- (1 - (w.c/c1)^2)^2 * I
  vec.mu <- apply(x * wc, 2, sum)/apply(wc, 2, sum)
  mat <- sweep(x, 2, vec.mu, FUN = "-")
  mat.c <- sweep(mat, 2, d, FUN = "/")
  rhoc <- ifelse(mat.c^2 <= c2^2, mat.c^2, c2^2)
  sig2 <- apply(rhoc, 2, mean) * d^2
  return(list(vec.mu=vec.mu,sig2=sig2))
}

################################################################################
#Initial values
################################################################################
"rogkmiss" = function(data, c1 = 4.5, c2 = 3){
  #
  # variables declaration
  #
  n <- nrow(data)
  p <- ncol(data)
  D <- matrix(0, p, p)
  M <- matrix(0, p, p)
  U <- matrix(0, p, p)
  k <- ((p - 1) * p)/2
  yp <- rep(0, n)
  ym <- rep(0, n)
  #
  # compute matrice D and standardized observations y
  #
  D.start <- gnan0(data, c1, c2)
  D.diag <- sqrt(D.start$sig2)
  y <- sweep(data, 2, D.diag, FUN = "/")
  diag(D) <- D.diag
  #
  # compute matrix U, eigenvectors and eigenvalues
  #
  i <- 1
  while(i <= p - 1) {
    j <- i + 1
    while(j <= p) {
      yp <- cbind(yp, y[, i] + y[, j], deparse.level = 1)
      ym <- cbind(ym, y[, i] - y[, j], deparse.level = 1)
      j <- j + 1
    }
    i = i + 1
  }
  y.plus <- yp[, 2:(k + 1)]
  y.moins <- ym[, 2:(k + 1)]
  d.plus <- apply(y.plus, 2, mad, na.rm = T)
  w <- sweep(y.plus, 2, apply(y.plus, 2, median, na.rm = T), FUN = "-")
  w.plus <- sweep(w, 2, d.plus, FUN = "/")
  I.plus <- ifelse(abs(w.plus) <= c1, 1, 0)
  wc.plus <- (1 - (w.plus/c1)^2)^2 * I.plus
  vec.mu.plus <- apply(y.plus * wc.plus, 2, sum, na.rm = T)/apply(wc.plus, 2, sum, na.rm = T)
  mat <- sweep(y.plus, 2, vec.mu.plus, FUN = "-")
  mat.plus <- sweep(mat, 2, d.plus, FUN = "/")
  rhoc.plus <- ifelse(mat.plus^2 <= c2^2, mat.plus^2, c2^2)
  sig2.plus <- apply(rhoc.plus, 2, mean, na.rm = T) * (d.plus^2)
  d.moins <- apply(y.moins, 2, mad, na.rm = T)
  w <- sweep(y.moins, 2, apply(y.moins, 2, median, na.rm = T), FUN = "-")
  w.moins <- sweep(w, 2, d.moins, FUN = "/")
  I.moins <- ifelse(abs(w.moins) <= c1, 1, 0)
  wc.moins <- (1 - (w.moins/c1)^2)^2 * I.moins
  vec.mu.moins <- apply(y.moins * wc.moins, 2, sum, na.rm = T)/apply(wc.moins, 2, sum, na.rm = T)
  mat <- sweep(y.moins, 2, vec.mu.moins, FUN = "-")
  mat.moins <- sweep(mat, 2, d.moins, FUN = "/")
  rhoc.moins <- ifelse(mat.moins^2 <= c2^2, mat.moins^2, c2^2)
  sig2.moins <- apply(rhoc.moins, 2, mean, na.rm = T) * (d.moins^2)
  vec.U <- 0.25 * (sig2.plus - sig2.moins)
  iter <- 0
  start <- 1
  vec.start <- p:1
  vec.end <- (p - 1):2
  for(i in 2:p) {
    iter <- iter + 1
    end <- sum((p - 1):vec.start[i])
    U[iter, i:p] <- vec.U[start:end]
    start <- start + vec.end[iter]
  }
  U <- U + t(U)
  diag(U) <- 1
  eig <- eigen(U)
  e <- eig$values
  E <- eig$vectors
  #
  # compute V(X) and t(X)
  #
  Z <- y %*% E
  A <- D %*% E
  Z.final <- gnan0(Z, c1, c2)
  diag(M) <- Z.final$sig2
  mu <- Z.final$vec.mu
  V <- A %*% M %*% t(A)
  t <- A %*% mu
  ####################################################3
  mal <- mahalanobis(data, t, V)
  qmal <- median(mal)
  qup <- qchisq(0.9, p)
  qdown <- qchisq(0.5, p)
  do <- (qup * qmal)/qdown
  w <- ifelse(mal <= do, 1, 0)
  #	rwdata <- sweep(data, 1, w, FUN = "*")
  #	tw <- apply(rwdata, 2, sum)/sum(w)
  #	mats <- matrix(0, p, p)
  #	for(i in 1:n) {
  #		xc <- as.vector(data[i,  ] - tw)
  #		mats <- mats + as.numeric(w[i]) * (xc %o% xc)
  #	}
  #	V <- mats/sum(w)
  return(list(center = t, cov = V,w=w))
}