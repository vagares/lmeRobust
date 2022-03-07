# This script computes the constants involved in the limiting
# covariances of the S-estimators as decsribed in Corollary 5
# for the case y|X has a multivariate normal distribution

moment=function(m,k,c0){
  # computes the constrained expectation E|Z|^m{|z|<=c} 
  # where Z has a multivariate N(0,I_k) distribution
  # m=the power of the moment
  # k=dimension of the normal distribution
  # c0=cut-off value of the rho-function
  const=2^(m/2)*gamma((m+k)/2)/gamma(k/2)
  gamint=pgamma((c0^2)/2,shape=(m+k)/2,scale=1)
  return(const*gamint)
}

momenttranslated=function(d,k,a,b){
  # computes E[|Z|^d{a<|Z|<b}]
  # needed for expectations of translated biweight functions
  return(moment(d,k,b)-moment(d,k,a))
  }

expecrho=function(k,c0){
  # computes E[rho(|Z|)]
  # k=dimension multivariate normal
  a1=1/2
  a2=-1/(2*(c0^2))
  a3=1/(6*(c0^4))
  a4=c0^2/6
  return(a1*moment(2,k,c0)+
           a2*moment(4,k,c0)+
           a3*moment(6,k,c0)+
           a4*(1-moment(0,k,c0)))
}

expecrhotranslated=function(k,M,c0){
  # computes E[rho_M(|Z|)]
  # M=first parameter of translated biweight
  # c0=second parameter of translated biweight
  # k=dimension of multivariate normal 
  
  # compute E[|Z|^2/2{|Z|>=M}]
  term1=moment(2,k,M)/2
  
  # compute the expectation over [M,M+c0]
  term2=momenttranslated(0,k,M,M+c0)*M^2/2-
    momenttranslated(0,k,M,M+c0)*
      M^2*(M^4-5*M^2*c0^2+15*c0^4)/(30*c0^4)+
    momenttranslated(2,k,M,M+c0)*(1/2+M^4/(2*c0^4)-M^2/c0^2)+
    momenttranslated(3,k,M,M+c0)*(4*M/(3*c0^2)-4*M^3/(3*c0^4))+
    momenttranslated(4,k,M,M+c0)*(3*M^2/(2*c0^4)-1/(2*c0^2))-
    momenttranslated(5,k,M,M+c0)*4*M/(5*c0^4)+
    momenttranslated(6,k,M,M+c0)/(6*c0^4)
  
  # compute expectation over [M+c0,Inf]
  term3=(1-moment(0,k,M+c0))*(M^2/2+c0*(5*c0+16*M)/30)
  
  return(term1+term2+term3)
}

rhoconst=function(k,r,a,b){
  # computes the cut=off constant c0 for which
  # E[rho(|Z|)]=r*rho(c0),
  # where rho is the biweight
  # k=dimension of the multivariate normal
  # r=breakdown point
  # a,b are the boundaries of the interval over which we use uniroot
  objecfun=function(c0){
    expecrho(k,c0)-r*biweightrho(c0,c0)
  }
  return(uniroot(f=objecfun,c(a,b))$root)
}

rhotranslatedconst=function(k,r,p,a,b){
  # computes the cut=off constant c0 for which
  # E[rho_t(|Z|)]=r*rho_t(M+c0),
  # where rho_t is the translated biweight
  # M chosen such that p=pchisq((M+c0)^2,df=k)
  # k=dimension of the multivariate normal
  # r=breakdown point
  # p=asymptotic rejection rate
  # a,b are the boundaries of the interval over which we use uniroot
  
  # compute E[rho_t(|Z|)]-r*rho_t(M+c0)
  objectfun=function(c0){
    # M chosen such that p=pchisq((M+c0)^2,df=k)
    M=sqrt(qchisq(1-p,df=k))-c0
    expecrhotranslated(k,M,c0)-
      r*biweightrhotranslated(M+c0,M,c0)
  }
  return(uniroot(f=objectfun,c(a,b))$root)
}

alpha=function(k,c0){
  # computes the value of the constant alpha in (8.4)
  # k = dimension of the normal distribution
  # c0 = cut-off value of the rho-function
  
  # Compute E[rho'(|Z|)/|Z|]
  term1=moment(0,k,c0)-2*moment(2,k,c0)/(c0^2)+moment(4,k,c0)/(c0^4)
  
  # Compute E[rho''(|Z|)]
  term2=moment(0,k,c0)-6*moment(2,k,c0)/(c0^2)+5*moment(4,k,c0)/(c0^4)
  
  return((1-1/k)*term1+(1/k)*term2)
}

alphatranslated=function(k,M,c0){
  # computes the value of the constant alpha in (8.4)
  # for the translated biweight rho function
  # k = dimension of the normal distribution
  # M = first parameter of the translated rho-function
  # c0 = second parameter of the translated rho-function

  # compute E[rho_t'(d)/d]
  term1=moment(0,k,M)+
    momenttranslated(0,k,M,M+c0)*(1+M^4/(c0^4)-2*M^2/(c0^2))+
    momenttranslated(1,k,M,M+c0)*(4*M/(c0^2)-4*M^3/(c0^4))+
    momenttranslated(2,k,M,M+c0)*(6*M^2/(c0^4)-2/(c0^2))-
    momenttranslated(3,k,M,M+c0)*4*M/(c0^4)+
    momenttranslated(4,k,M,M+c0)/(c0^4)
  
  # compute E[rho_t''(d)]
  term2=moment(0,k,M)+
    momenttranslated(0,k,M,M+c0)*(1+M^4/(c0^4)-2*M^2/(c0^2))+
    momenttranslated(1,k,M,M+c0)*(8*M/(c0^2)-8*M^3/(c0^4))+
    momenttranslated(2,k,M,M+c0)*(18*M^2/(c0^4)-6/c0^2)-
    momenttranslated(3,k,M,M+c0)*16*M/(c0^4)+
    momenttranslated(4,k,M,M+c0)*5/c0^4
  
  return((1-1/k)*term1+(1/k)*term2)
}

momentcentered=function(l,m,k,M,a,b){
  # computes E[|Z|^l(|Z|-M)^m{a<|Z|<b}]
  # need for E[rho_t'(|Z|)^2] in alpha_t
  ss=0
  for (i in 0:m){
    ss=ss+choose(m,i)*((-M)^i)*momenttranslated(l+m-i,k,a,b)
  }
  return(ss)
}


psisquared=function(k,c0){
  # computes the numerator E[rho'(|Z|)^2] appearing in the 
  # limiting covariance of betahat
  # see Corollary 5
  expectation=moment(2,k,c0)-4*moment(4,k,c0)/(c0^2)+
    6*moment(6,k,c0)/(c0^4)-4*moment(8,k,c0)/(c0^6)+
    moment(10,k,c0)/(c0^8)
  return(expectation)
}

psitranslatedsquared=function(k,M,c0){
  # computes the numerator E[rho_t'(|Z|)^2] appearing in the 
  # limiting covariance of betahat
  # see Corollary 5

  # computes E[|Z|{|Z|<=M}]
  term1=moment(2,k,M)
  
  # computes expectation over [M,M+c0]
  term2=momenttranslated(2,k,M,M+c0)-
    (4/(c0^2))*momentcentered(2,2,k,M,M,M+c0)+
    (6/(c0^4))*momentcentered(2,4,k,M,M,M+c0)-
    (4/(c0^6))*momentcentered(2,6,k,M,M,M+c0)+
    (1/(c0^8))*momentcentered(2,8,k,M,M,M+c0)
  
  return(term1+term2)
}

constbetahat=function(k,c0){
  # determines the scalar in the limiting covariance of betahat
  # k=dimension of multivariate normal
  # c0=pre-determined cut-off value for given BDP r
  return(psisquared(k,c0)/(k*alpha(k,c0)^2))
}

constbetahattranslated=function(k,M,c0){
  # determines the scalar in the limiting covariance 
  # of betahat for the translated biweight
  # k=dimension of multivariate normal
  # M is chosen sqrt(qchisq(1-pi),df=k)-c0
  # where pi=asymptotic rejection rate 
  # c0=pre-determined cut-off value for given BDP r
  return(psitranslatedsquared(k,M,c0)/(k*alphatranslated(k,M,c0)^2))
}


sigma1=function(k,c0){
  # determines the constant sigma1 in the scalar in the 
  # limiting covariance of thetahat
  # k=dimension of multivariate normal
  # c0=pre-determined cut-off value for given BDP r
  
  # Compute E[u(d)^2d^4]
  term1=moment(4,k,c0)-
    4*moment(6,k,c0)/(c0^2)+
    6*moment(8,k,c0)/(c0^4)-
    4*moment(10,k,c0)/(c0^6)+
    moment(12,k,c0)/(c0^8)
  
  # Compute E[rho''(d)d^2]
  term2=moment(2,k,c0)-
    6*moment(4,k,c0)/(c0^2)+
    5*moment(6,k,c0)/(c0^4)
  
  # Compute E[rho'(d)d]
  term3=moment(2,k,c0)-
    2*moment(4,k,c0)/(c0^2)+
    moment(6,k,c0)/(c0^4)
  
  num=k*(k+2)*term1
  denum=(term2+(k+1)*term3)^2
  return(num/denum)
}

sigma1t=function(k,M,c0){
  # determines the constant sigma1 in the scalar in the 
  # limiting covariance of thetahat with translated biweight
  # k=dimension of multivariate normal
  # M= first parameter of translated biweight
  # c0= first parameter of translated biweight
  
  # Compute E[u_t(d)^2d^4]
  term1=moment(4,k,M)+
    momenttranslated(4,k,M,M+c0)-
    momentcentered(4,2,k,M,M,M+c0)*(4/(c0^2))+
    momentcentered(4,4,k,M,M,M+c0)*(6/(c0^4))-
    momentcentered(4,6,k,M,M,M+c0)*(4/(c0^6))+
    momentcentered(4,8,k,M,M,M+c0)/(c0^8)
  
  # Compute E[rho_t''(d)d^2]
  term2=moment(2,k,M)+
    momenttranslated(2,k,M,M+c0)*(1+M^4/(c0^4)-2*M^2/(c0^2))+
    momenttranslated(3,k,M,M+c0)*(8*M/(c0^2)-8*M^3/(c0^4))+
    momenttranslated(4,k,M,M+c0)*(18*M^2/c0^4-6/(c0^2))-
    momenttranslated(5,k,M,M+c0)*16*M/(c0^4)+
    momenttranslated(6,k,M,M+c0)*5/(c0^4)

  # Compute E[rho_t'(d)d]
  term3=moment(2,k,M)+
    momenttranslated(2,k,M,M+c0)*(1+M^4/(c0^4)-2*M^2/(c0^2))+
    momenttranslated(3,k,M,M+c0)*(4*M/(c0^2)-4*M^3/(c0^4))+
    momenttranslated(4,k,M,M+c0)*(6*M^2/c0^4-2/(c0^2))-
    momenttranslated(5,k,M,M+c0)*4*M/(c0^4)+
    momenttranslated(6,k,M,M+c0)/(c0^4)
  
  num=k*(k+2)*term1
  denum=(term2+(k+1)*term3)^2
  return(num/denum)
}

sigma2=function(k,c0){
  # determines the constant sigma2 in the scalar in the 
  # limiting covariance of thetahat
  # k=dimension of multivariate normal
  # c0=pre-determined cut-off value for given BDP r
  a1=1/2
  a2=-1/(2*(c0^2))
  a3=1/(6*(c0^4))
  a4=c0^2/6
  
  # Compute E[rho(d)^2]
  term1=(a1^2)*moment(4,k,c0)+
    (a2^2)*moment(8,k,c0)+
    (a3^2)*moment(12,k,c0)+
    2*a1*a2*moment(6,k,c0)+
    2*a1*a3*moment(8,k,c0)+
    2*a2*a3*moment(10,k,c0)+
    a4^2*(1-moment(0,k,c0))
  
  # Compute E[rho'(d)d]
  term2=moment(2,k,c0)-
    2*moment(4,k,c0)/(c0^2)+
    moment(6,k,c0)/(c0^4)
  
  # Compute E[rho(d)]
  term3=expecrho(k,c0)
  
  num=4*(term1-(term3)^2)
  denum=(term2)^2
  return(-(2/k)*sigma1(k,c0)+num/denum)
}

sigma2t=function(k,M,c0){
  # determines the constant sigma2 in the scalar in the 
  # limiting covariance of thetahat of the translated biweight
  # k=dimension of multivariate normal
  # M=first parameter of translated biweight
  # c0=second parameter of translated biweight
  
  a1t=M^2/2
  a2t=-M^2*(M^4-5*M^2*c0^2+15*c0^4)/(30*c0^4)
  a3t=1/2+M^4/(2*c0^4)-M^2/c0^2
  a4t=4*M/(3*c0^2)-4*M^3/(3*c0^4)
  a5t=3*M^2/(2*c0^4)-1/(2*c0^2)
  a6t=-4*M/(5*c0^4)
  a7t=1/(6*c0^4)
  a8t=M^2/2+c0*(5*c0+16*M)/30

  # Compute E[rho_t(d)^2]
  term1=moment(4,k,M)/4+
    a1t^2*momenttranslated(0,k,M,M+c0)+
    a2t^2*momenttranslated(0,k,M,M+c0)+
    a3t^2*momenttranslated(4,k,M,M+c0)+
    a4t^2*momenttranslated(6,k,M,M+c0)+
    a5t^2*momenttranslated(8,k,M,M+c0)+
    a6t^2*momenttranslated(10,k,M,M+c0)+
    a7t^2*momenttranslated(12,k,M,M+c0)+
    2*a1t*a2t*momenttranslated(0,k,M,M+c0)+
    2*a1t*a3t*momenttranslated(2,k,M,M+c0)+
    2*a1t*a4t*momenttranslated(3,k,M,M+c0)+
    2*a1t*a5t*momenttranslated(4,k,M,M+c0)+
    2*a1t*a6t*momenttranslated(5,k,M,M+c0)+
    2*a1t*a7t*momenttranslated(6,k,M,M+c0)+
    2*a2t*a3t*momenttranslated(2,k,M,M+c0)+
    2*a2t*a4t*momenttranslated(3,k,M,M+c0)+
    2*a2t*a5t*momenttranslated(4,k,M,M+c0)+
    2*a2t*a6t*momenttranslated(5,k,M,M+c0)+
    2*a2t*a7t*momenttranslated(6,k,M,M+c0)+
    2*a3t*a4t*momenttranslated(5,k,M,M+c0)+
    2*a3t*a5t*momenttranslated(6,k,M,M+c0)+
    2*a3t*a6t*momenttranslated(7,k,M,M+c0)+
    2*a3t*a7t*momenttranslated(8,k,M,M+c0)+
    2*a4t*a5t*momenttranslated(7,k,M,M+c0)+
    2*a4t*a6t*momenttranslated(8,k,M,M+c0)+
    2*a4t*a7t*momenttranslated(9,k,M,M+c0)+
    2*a5t*a6t*momenttranslated(9,k,M,M+c0)+
    2*a5t*a7t*momenttranslated(10,k,M,M+c0)+
    2*a6t*a7t*momenttranslated(11,k,M,M+c0)+
    a8t^2*(1-moment(0,k,M+c0))

  # Compute E[rho_t'(d)d]
  term2=moment(2,k,M)+
    momenttranslated(2,k,M,M+c0)*(1+M^4/(c0^4)-2*M^2/(c0^2))+
    momenttranslated(3,k,M,M+c0)*(4*M/(c0^2)-4*M^3/(c0^4))+
    momenttranslated(4,k,M,M+c0)*(6*M^2/c0^4-2/(c0^2))-
    momenttranslated(5,k,M,M+c0)*4*M/(c0^4)+
    momenttranslated(6,k,M,M+c0)/(c0^4)
  
  # Compute E[rho_t(d)]
  term3=expecrhotranslated(k,M,c0)
  
  num=4*(term1-(term3)^2)
  denum=(term2)^2
  return(-(2/k)*sigma1t(k,M,c0)+num/denum)
}

eta=function(k,c0){
  # computes the efficiency index of Tyler
  ifelse(k==1,1,0)*(2*sigma1(k,c0)+sigma2(k,c0))+
    ifelse(k>1,1,0)*sigma1(k,c0)
}

etat=function(k,M,c0){
  # computes the efficiency index of Tyler for translated biweight
  ifelse(k==1,1,0)*(2*sigma1t(k,M,c0)+sigma2t(k,M,c0))+
    ifelse(k>1,1,0)*sigma1t(k,M,c0)
}

