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
