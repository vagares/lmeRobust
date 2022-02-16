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
biwt.rho2<- function(d,c,M){
  hulp =M^2/2 -
    (M^2*(M^4-5*M^2*c^2+15*c^4))/(30*c^4)+
    d^2*(0.5+(M^4)/(2*c^4)-M^2/c^2) +
    d^3*((4*M)/(3*c^2)-(4* M^3)/(3*c^4))+
    d^4*((3*M^2)/(2*c^4)-1/(2*c^2))
  -(4*M*d^5)/(5*c^4)+ 
    d^6/(6*c^4)
  if(d<M){rho=d^2/2}else{if(d>=M & (d<=(M+c))){rho=hulp}else{rho=M^2/2 + c*(5*c+16*M)/30}}
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
  biweightu2translated(d,M,c0)-biweightrhotranslated(d,M,c0)+b0}

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