library(MASS) # needed to generate from multivariate normal
library(ks)
library(robustbase)
library(nlme)


source("Biweight functions.R")
source("ASYMP NORM constants.R")

# From Potthoff and Roy (1964)
# Also in Demidenko (2013)
growth=read.csv("Growth.csv", header=TRUE,sep=";")
n=nrow(growth)

ymax=max(growth[,3:6])
ymin=min(growth[,3:6])

plot(1:4,growth[1,3:6],type="l",xaxt="n",
     ylim=c(ymin,ymax),xlim=c(1,4.5),
     main="Orthodontic Growth data",
     xlab="Age (years)",
     ylab="Bone density (g/cm2)")
axis(1,labels=8:14,at=seq(1,4,by=0.5))
for(i in 2:27){
  lines(1:4,growth[i,3:6],lty=i)
}

for (i in 12:27){
  lines(1:4,growth[i,3:6],lty=i,col=2)
}
legend("topright",legend=1:n,lty=1:n,
       col=c(rep(1,11),rep(2,16)),
       title="Subject",bty="n")
legend("topleft",legend=c("Girls","Boys"),col=1:2,lty=1,bty="n")


#############################################################################
# ESTIMATES REPORTED IN PHD COPT AND HERITIER
#############################################################################
# MLE, TABLE 8-7 IN PHD COPT
betaMLgrowthCopt=c(16.34,1.03,0.78,-0.30)
thetaMLgrowthCopt=c(2.336,0.007,1.894)

#############################################################################
# FROM LME FUNCTION
options(contrasts=c("contr.SAS","contr.poly"))
Dental=groupedData( distance ~ age | Subject, 
                    data = as.data.frame( Dental ),
                    FUN = mean,
                    outer = ~ Sex,
                    labels = list( x = "Age",
                                   y = "Distance from pituitary to pterygomaxillary fissure" ),
                    units = list( x = "(yr)", y = "(mm)") )


DentalML.lme <- lme(distance~Sex*age, data = Dental,random = list(Subject=pdDiag(~age)),method="ML")
summary(DentalML.lme)
betaMLgrowthLME=fixed.effects(Dental.lme)
thetaMLgrowthLME=c(1.499741, 0.08220457, 1.350634)^2

DentalREML.lme <- lme(distance~Sex*age, data = Dental,random = list(Subject=pdDiag(~age)),method="REML")
summary(DentalREML.lme)
betaREMLgrowthLME=fixed.effects(Dental.lme)
thetaREMLgrowthLME=c(1.554607,0.08801656, 1.365502)^2


###############################################################################
# CTBS TABLE 8-7 IN PHD COPT
betaCTBSgrowthCopt=c(17.24,0.14,0.70,-0.22)
thetaCTBSgrowthCopt=c(2.812,0.012,1.032)

###############################################################################
# CBS-MM TABLE 4.5 IN HERITIER
betaCBSMMgrowthHeritier=c(17.395,0.080,0.581,-0.110)
SQRTthetaCBSgrowthHeritier=c(1.584,0.115,1.04)

###############################################################################
# MODEL
t=c(8,10,12,14)

Ymat=growth[,3:6]

Ymat=as.matrix(Ymat)
n=nrow(Ymat)
k=ncol(Ymat)

# DESIGN MATRICES FOR BOYS AND GIRLS, SEE REMARK PAGE 91
Xb=cbind(c(1,1,1,1),c(0,0,0,0),c(8,10,12,14),c(0,0,0,0))
Xg=cbind(c(1,1,1,1),c(1,1,1,1),c(8,10,12,14),c(8,10,12,14))
q=4

Z1=matrix(rep(1,4),ncol=1)
Z2=t
varerror=diag(rep(1,4))

L1=Z1%*%t(Z1)
L2=Z2%*%t(Z2)
L3=varerror
l=3

###############################################################################
# COMPUTATION OF CTBS/CBS/ML ESTIMATES
###############################################################################
# Setting the breakdown point and cut-off constant for translated biweight
p=0.01
r=0.5
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)

# TO COMPUTE mle
#M=10000

# TO COMPUTE CTBS
c0=rhoconst(k,r,0.01,100)
b0=expecrho(k,c0)
#M=0
#c0t=c0

###############################################################################
# SETTING STARTING VALUES FOR MEAN AND COVARIANCE
###############################################################################

# OPTION 1: location-covariance from OGK
Vstart=covOGK(Ymat,sigmamu=s_mad)
mu0=Vstart$center
V0=Vstart$cov
w=Vstart$weights

# OPTION 2: location-covariance from weighted OGK
Vstart=covOGK(Ymat,sigmamu=s_mad)
mu0=Vstart$wcenter
V0=Vstart$wcov
w=Vstart$weights

# OPTION 3: location-covariance from  rogkmiss
Vstart=rogkmiss(as.matrix(Ymat))
mu0=Vstart$center
V0=Vstart$cov
w=Vstart$w

# OPTION 4: location-covariance from MCD
Vstart=covMcd(Ymat)
mu0=Vstart$center
V0=Vstart$cov
w=Vstart$mcd.wt

# OPTION 5: sample mean and covariance
Vstart=cov.wt(Ymat)
mu0=Vstart$center
V0=Vstart$cov
w=rep(1,times=n)


###########################################################################
# TRANSFER TO STARTING VALUES FOR BETA AND THETA
###########################################################################
# OPTION 1: RUN 1 ITERATION STEP
betavecterm=matrix(0,nrow=q,ncol=1) # vector-term in step 5
betamatterm=matrix(0,nrow=q,ncol=q) # matrix-term in step 5
for (i in 1:n){
  X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
  y=as.matrix(Ymat[i,],ncol=1,nrow=q)
  betavecterm=betavecterm+w[i]*t(X)%*%solve(V0)%*%y  
  betamatterm=betamatterm+w[i]*t(X)%*%solve(V0)%*%X 
}

betaoldtranslated=solve(betamatterm)%*%betavecterm
Voldtranslated=V0
thetaoldtranslated=numeric(l)

###########################################################################
# OPTION 2: SOLVING BETA
# setting starting values for beta and theta
Ymatb=Ymat[growth[,2]=="M",]
Ymatg=Ymat[growth[,2]=="F",]

# OPTION 3: Preparing starting values MCD
Vstartb=cov.mcd(Ymatb)
Vstartg=cov.mcd(Ymatg)
Vstart=cov.mcd(Ymat)

# STARTING VALUES FROM ROGKMISS
Vstartb=rogkmiss(as.matrix(Ymatb))
Vstartg=rogkmiss(as.matrix(Ymatg))
Vstart=rogkmiss(as.matrix(Ymat))

mu0b=Vstartb$center
mu0g=Vstartg$center

X0=cbind(c(1,1,1,1),t)
betaB=solve(t(X0)%*%X0)%*%t(X0)%*%mu0b
betaG=solve(t(X0)%*%X0)%*%t(X0)%*%mu0g-betaB
betaoldtranslated=c(betaB[1],betaG[1],betaB[2],betaG[2])


Voldtranslated=Vstart$cov

thetaoldtranslated=numeric(l)

###########################################################################
# RUNNING THE ITERATION TO OBTAIN S-ESTIMATES
###########################################################################

# stop criterion
tol=1
iter=0

while (tol>=1e-06){
  ##############################################################################
  # Iteration step
  #############################################
  
  # Re-scaling Mahalanobis distances to satisfy S-constraint
  # array for n Mahalanobis distances
  MD=numeric(n)
  for (i in 1:n){
    X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
    y=Ymat[i,]
    MD[i]=mahalanobis(y,center=X%*%betaoldtranslated,cov=Voldtranslated) # Note MD=d^2!
  }
  
  # determining scaling constant for MD
  # for translated biweight
  objfuntranslated=function(s){
    mean(biweightrhotranslated(sqrt(MD)/s,M,c0t))-expecrhotranslated(k,M,c0t)
  }
  MDscaletranslated=uniroot(f=objfuntranslated,c(0.01,100))$root
  
  
  # OPTION 1: re-scaling MD to satisfy S-constraint
  MDtildetranslated=MD/MDscaletranslated^2
  
  # OPTION 2: without adapting the Mahalanobis distance
  #MDtildetranslated=MD
  
  # OPTION 3: Rocke's adaptation as in Copt
  # CORRECTION USED IN PHD COPT? SEE (4.14)
  qconst=round((n+k+1)/2)
  kconstsq=sort(MD)[qconst]/qchisq(qconst/(n+1),df=k)
  #MDtildetranslated=MD/kconstsq
  
  # prepare for updating regression estimate in step 5
  betavectermtranslated=matrix(0,nrow=q,ncol=1) # vector-term in step 5
  betamattermtranslated=matrix(0,nrow=q,ncol=q) # matrix-term in step 5
  
  # prepare for updating covariance estimate in step 6
  Utranslated=matrix(0,nrow=l,ncol=1) # vector U in step 6
  Qtranslated=matrix(0,nrow=l,ncol=l) # matrix Q in step 6
  vtranslated=numeric(1)              # sum of v(sqrt(MDtilde)) with v(s)=u(s)s^2-rho(s)+b0
  
  # STEP 5-6: Updating steps
  # constructing summations in step 5
  Utermtranslated=matrix(0,nrow=l,ncol=n)
  vtermtranslated=numeric(n)
  for (i in 1:n){
    X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
    y=matrix(as.numeric(Ymat[i,]),ncol=1,nrow=k)
    betavectermtranslated=
      betavectermtranslated+
      biweightutranslated(sqrt(MDtildetranslated[i]),M,c0t)*
      t(X)%*%solve(Voldtranslated)%*%y
    
    betamattermtranslated=
      betamattermtranslated+
      biweightutranslated(sqrt(MDtildetranslated[i]),M,c0t)*t(X)%*%solve(Voldtranslated)%*%X 
    
    vtermtranslated[i]=biweightu2translated(sqrt(MDtildetranslated[i]),M,c0t)
    vtranslated=vtranslated+vtermtranslated[i]  
    
    # constructing summation components of U and matrix Q in step 6
    for (s in 1:l){
      Ls=(s==1)*L1+(s==2)*L2+(s==3)*L3
      Utermtranslated[s,i]=
        k*biweightutranslated(sqrt(MDtildetranslated[i]),M,c0t)*
        t(y-X%*%betaoldtranslated)%*%solve(Voldtranslated)%*%
        Ls%*%
        solve(Voldtranslated)%*%(y-X%*%betaoldtranslated)
      Utranslated[s]=Utranslated[s]+Utermtranslated[s,i]
      for (t in 1:3){
        Lt=(t==1)*L1+(t==2)*L2+(t==3)*L3
        Qtranslated[s,t]=
          sum(diag(solve(Voldtranslated)%*%Ls%*%
                     solve(Voldtranslated)%*%Lt))
      } 
    }
  }
  
  # STEP5 updating regression estimate
  betanewtranslated=solve(betamattermtranslated)%*%betavectermtranslated
  
  # STEP6 updating covariance estimate
  thetanewtranslated=(1/vtranslated)*solve(Qtranslated)%*%Utranslated
  Vnewtranslated=
    thetanewtranslated[1]*L1+thetanewtranslated[2]*L2+
    thetanewtranslated[3]*L3
  
  # End of iteration step
  
  iter=iter+1
  
  diffbeta=norm(betanewtranslated-betaoldtranslated,type="F")
  difftheta=norm(thetanewtranslated-thetaoldtranslated,type="F")
  tol=diffbeta+difftheta
  
  # setting starting values for the next iteration within while-loop
  betaoldtranslated=betanewtranslated
  thetaoldtranslated=thetanewtranslated
  Voldtranslated=Vnewtranslated
}
###########################################################################
iter
betanewtranslated
thetanewtranslated
sqrt(thetanewtranslated)

# SAVING THE DIFFERENT ESTIMATES AFTER RUNNING THE ITERATION

betaMLgrowthRik=betanewtranslated
thetaMLgrowthRik=thetanewtranslated

betaCTBSgrowthRik=betanewtranslated
thetaCTBSgrowthRik=thetanewtranslated

betaCBSgrowthRik=betanewtranslated
thetaCBSgrowthRik=thetanewtranslated

betaCTBSgrowthRESCALE=betanewtranslated
thetaCTBSgrowthRESCALE=thetanewtranslated

betaCBSgrowthRESCALE=betanewtranslated
thetaCBSgrowthRESCALE=thetanewtranslated




###########################################################################
# RESULTS MLE
data.frame(betaMLrik=betaMLgrowthRik,
           betaMLCopt=betaMLgrowthCopt,
           betaMLlme=betaMLgrowthLME,
           betaREMLlme=betaREMLgrowthLME)

data.frame(thetaMLrik=thetaMLgrowthRik,
           thetaMLCopt=thetaMLgrowthCopt,
           thetaMLlme=thetaMLgrowthLME,
           thetaREMLlme=thetaREMLgrowthLME)

###########################################################################
# RESULTS CTBS
data.frame(betaCTBSrik=betaCTBSgrowthRik,
           betaCTBSCopt=betaCTBSgrowthCopt)

data.frame(thetaCTBSrik=thetaCTBSgrowthRik,
           thetaCTBSCopt=thetaCTBSgrowthCopt)

###########################################################################
# RESULTS CBS-MM
data.frame(betaCBSrik=betaCBSgrowthRik,
           betaCBSMMrik=betaCBSMMgrowthRik,
           betaCBSMMheritier=betaCBSMMgrowthHeritier)

est=betaCBSMMgrowthHeritier
rbind(est[1]-est[2],
2*est[2],
est[3]-est[4],
2*est[4])
A=as.matrix(cbind(c(1,0,0,0),c(-1,2,0,0),c(0,0,1,0),c(0,0,-1,2)))
solve(A)


data.frame(thetaCBSrik=thetaCBSgrowthRik,
           thetaCBSheritier=SQRTthetaCBSgrowthHeritier^2)
data.frame(SQRTthetaCBSrik=sqrt(thetaCBSgrowthRik),
           SQRTthetaCBSheritier=SQRTthetaCBSgrowthHeritier)

###########################################################################
# RESULTS CTBS RE-SCALED
data.frame(betaCTBSrescale=betaCTBSgrowthRESCALE,
           betaCTBSCopt=betaCTBSgrowthCopt)

data.frame(thetaCTBSrescale=thetaCTBSgrowthRESCALE,
           thetaCTBSCopt=thetaCTBSgrowthCopt)

###########################################################################
# RESULTS CBSMM RE-SCALED
data.frame(betaCBSrescale=betaCBSgrowthRESCALE,
           betaCBSMMrescale=betaCBSMMgrowthRESCALE,
           betaCBSMMheritier=betaCBSMMgrowthHeritier)
A%*%betaCBSMMgrowthHeritier

data.frame(thetaCBSrescale=thetaCBSgrowthRESCALE,
           thetaCBSheritier=SQRTthetaCBSgrowthHeritier^2)

data.frame(SQRTthetaCBSrescale=sqrt(thetaCBSgrowthRESCALE),
           SQRTthetaCBSheritier=SQRTthetaCBSgrowthHeritier)


#############################################################################
# COMPARISION OF DETERMINANTS
#############################################################################
beta=betaCBSgrowthRik
theta=thetaCBSgrowthRik
Vtheta=theta[1]*L1+theta[2]*L2+theta[3]*L3

MD=numeric(n)
for (i in 1:n){
  X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
  y=Ymat[i,]
  MD[i]=mahalanobis(y,center=X%*%beta,cov=Vtheta) # Note MD=d^2!
}
objfuntranslated=function(s){
  mean(biweightrhotranslated(sqrt(MD)/s,M,c0t))-expecrhotranslated(k,M,c0t)
}
MDscaletranslated=uniroot(f=objfuntranslated,c(0.01,100))$root

detCTBSgrowthRESCALE=det(Vtheta*MDscaletranslated^2)
detCTBSgrowthCopt=det(Vtheta*MDscaletranslated^2)

detCBSgrowthRESCALE=det(Vtheta*MDscaletranslated^2)
detCBSgrowthHeritier=det(Vtheta*MDscaletranslated^2)

data.frame(DetCTBSrescale=detCTBSgrowthRESCALE,
           DetCTBSCopt=detCTBSgrowthCopt)

data.frame(DetCBSrescale=detCBSgrowthRESCALE,
           DetCBSHeritier=detCBSgrowthHeritier)

##########################################################################
# STANDARD ERRORS
#########################################################################
# Covariance beta alpha*(E(X^T*V^{-1}*X))^{-1}

#########################################################################
#MLE
M=10000
beta=betaMLgrowthRik
theta=thetaMLgrowthRik
VthetaMLE=theta[1]*L1+theta[2]*L2+theta[3]*L3

# Asymptotic variances betahat
# Expression from our theory

# ESTIMATE FOR X^T*v^{-1}*X
COVestMLErik=matrix(0,nrow=q,ncol=q)
for (i in 1:n){
  X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
  COVestMLErik=COVestMLErik+t(X)%*%solve(VthetaMLE)%*%X
}
varbetatranslatedMLE=constbetahattranslated(k,M,c0t)*solve(COVestMLErik/n)

SEbetaMLgrowthRik=sqrt(diag(varbetatranslatedMLE)/n)
SEbetaMLgrowthCopt=c(1.139,1.480,0.096,0.125)
SEbetaMLgrowthLME=c(0.9408690 ,1.4740585 ,0.0794421 ,0.1244618 )

data.frame(SEbetaMLErik=SEbetaMLgrowthRik,
           SEbetaMLEcopt=SEbetaMLgrowthCopt,
           SEbetaMLElme=SEbetaMLgrowthLME)

############################################################
# CTBS with re-scaling
p=0.01
r=0.5
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)

beta=betaCTBSgrowthRESCALE
theta=thetaCTBSgrowthRESCALE
VthetaCTBSrescale=theta[1]*L1+theta[2]*L2+theta[3]*L3


# ESTIMATE FOR X^T*v^{-1}*X
COVestCTBSrescale=matrix(0,nrow=q,ncol=q)
for (i in 1:n){
  X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
  COVestCTBSrescale=COVestCTBSrescale+t(X)%*%solve(VthetaCTBSrescale)%*%X
}
varbetatranslatedCTBSrescale=constbetahattranslated(k,M,c0t)*solve(COVestCTBSrescale/n)

SEbetaCTBSgrowthRESCALE=sqrt(diag(varbetatranslatedCTBSrescale)/n)
SEbetaCTBSgrowthCopt=c(1.037,1.347,0.086,0.112)

data.frame(SEbetaCTBSrescale=SEbetaCTBSgrowthRESCALE,
           SEbetaCTBSCopt=SEbetaCTBSgrowthCopt)

############################################################
# CBS-MM with re-scaling
p=0.01
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)
M=0
c0t=c1

beta=betaCBSMMgrowthRESCALE
theta=thetaCBSgrowthRESCALE
VthetaCBSrescale=theta[1]*L1+theta[2]*L2+theta[3]*L3

# ESTIMATE FOR X^T*v^{-1}*X
COVestCBSMMrescale=matrix(0,nrow=q,ncol=q)
for (i in 1:n){
  X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
  COVestCBSMMrescale=COVestCBSMMrescale+t(X)%*%solve(VthetaCBSrescale)%*%X
}
varbetatranslatedCBSMMrescale=constbetahattranslated(k,M,c0t)*solve(COVestCBSMMrescale/n)

SEbetaCBSMMgrowthRESCALE=sqrt(diag(varbetatranslatedCBSMMrescale)/n)
SEbetaCBSMMgrowthHeritier=c(0.613,0.613,0.052,0.052)

data.frame(SEbetaCBSMMrescale=SEbetaCBSMMgrowthRESCALE,
           SEbetaCBSMMgrowthHeritier=SEbetaCBSMMgrowthHeritier)

# SE for Heritier from out results
sqrt(diag(solve(A)%*%varbetatranslatedCBSMMrescale%*%t(solve(A)))/n)


############################################################################
# CTBS-MM with re-scaling
p=0.01
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)
M=0
c0t=c1

beta=betaCTBSMMgrowthRESCALE
theta=thetaCTBSgrowthRESCALE
VthetaCTBSrescale=theta[1]*L1+theta[2]*L2+theta[3]*L3

# ESTIMATE FOR X^T*v^{-1}*X
COVestCTBSMMrescale=matrix(0,nrow=q,ncol=q)
for (i in 1:n){
  X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
  COVestCTBSMMrescale=COVestCTBSMMrescale+t(X)%*%solve(VthetaCTBSrescale)%*%X
}
varbetatranslatedCTBSMMrescale=constbetahattranslated(k,M,c0t)*solve(COVestCTBSMMrescale/n)

SEbetaCTBSMMgrowthRESCALE=sqrt(diag(varbetatranslatedCTBSMMrescale)/n)
SEbetaCBSMMgrowthHeritier=c(0.613,0.613,0.052,0.052)

data.frame(SEbetaCTBSMMrescale=SEbetaCTBSMMgrowthRESCALE,
           SEbetaCBSMMrescale=SEbetaCBSMMgrowthRESCALE,
           SEbetaCBSMMgrowthHeritier=SEbetaCBSMMgrowthHeritier)


############################################################
# WALD TEST
###########################################################
# MLE
waldMLE=betaMLEgrowthRik/SEbetaMLEgrowthRik
pMLE=2*(1-pnorm(abs(waldMLE)))
2*(1-pt(abs(waldMLE),df=29))

waldMLEcopt=betaMLEgrowthCopt/SEbetaMLEgrowthCopt
pMLEcopt=2*(1-pnorm(abs(waldMLEcopt)))


data.frame(MLE=betaMLEgrowthRik,
           SE.MLE=round(SEbetaMLEgrowthRik,3),
           TMLE=waldMLE,
           PvalMLE=round(pMLE,4),
           PvalMLCopt=round(pMLEcopt,3))

# CTBS
waldCTBS=betaCTBSgrowthRESCALE/SEbetaCTBSgrowthRESCALE
pCTBS=2*(1-pnorm(abs(waldCTBS)))

waldCTBScopt=betaCTBSgrowthCopt/SEbetaCTBSgrowthCopt
pCTBScopt=2*(1-pnorm(abs(waldCTBScopt)))


data.frame(CTBS=betaCTBSgrowthRESCALE,
           SE.CTBS=round(SEbetaCTBSgrowthRESCALE,3),
           TCTBS=waldCTBS,
           PvalCTBS=round(pCTBS,4),
           PvalCTBScopt=round(pCTBScopt,4))

# CBS-MM
waldCBSMM=betaCBSMMgrowthRESCALE/SEbetaCBSMMgrowthRESCALE
pCBSMM=2*(1-pnorm(abs(waldCBSMM)))

waldCBSMMheritier=betaCBSMMgrowthHeritier/SEbetaCBSMMgrowthHeritier
pCBSMMheritier=2*(1-pnorm(abs(waldCBSMMheritier)))


data.frame(CBSMM=betaCBSMMgrowthRESCALE,
           SE.CBSMM=round(SEbetaCBSMMgrowthRESCALE,3),
           TCBSMM=waldCBSMM,
           PvalCBS=round(pCBSMM,4),
           PvalCBSheritier=round(pCBSMMheritier,4))

# CBS-MM comprable with Heritier
betaCTBSMMgrowthRESCALEadapted=solve(A)%*%betaCTBSMMgrowthRESCALE
SEbetaCTBSgrowthRESCALEadapted=sqrt(diag(solve(A)%*%varbetatranslatedCBSMMrescale%*%t(solve(A)))/n)

waldCTBSMMadapted=betaCTBSMMgrowthRESCALEadapted/SEbetaCTBSgrowthRESCALEadapted
pCTBSMMadapted=2*(1-pnorm(abs(waldCTBSMMadapted)))

data.frame(CTBSMM=betaCTBSMMgrowthRESCALEadapted,
           SE.CTBSMM=round(SEbetaCTBSgrowthRESCALEadapted,3),
           TCBSMM=waldCTBSMMadapted,
           PvalCTBS=round(pCTBSMMadapted,4))

############################################################################################
# RESISDUALS
############################################################################################
YMLEmat=matrix(0,nrow=n,ncol=k)
YCTBSmat=matrix(0,nrow=n,ncol=k)
YCBSMMmat=matrix(0,nrow=n,ncol=k)
resMLEgrowth=matrix(0,nrow=n,ncol=k)
resCTBSgrowth=matrix(0,nrow=n,ncol=k)
resCBSMMgrowth=matrix(0,nrow=n,ncol=k)

for (i in 1:n){
  X=ifelse(growth[i,2]=="F",1,0)*Xg+ifelse(growth[i,2]=="M",1,0)*Xb
  YMLEmat[i,]=X%*%betaMLEgrowthLME
  YCTBSmat[i,]=X%*%betaCTBSgrowthRik
  YCBSMMmat[i,]=X%*%betaCBSMMgrowthRik
}
resMLEgrowth=Ymat-YMLEmat
resCTBSgrowth=Ymat-YCTBSmat
resCBSMMgrowth=Ymat-YCBSMMmat

############################################################################################
# MAHALANOBIS DISTANCES OF RESIDUALS
a=0.025

MDerrorMLEgrowth=numeric(n)
MDerrorCTBSgrowth=numeric(n)
MDerrorCBSMMgrowth=numeric(n)


for (i in 1:n){
  MDerrorMLEgrowth[i]=mahalanobis(resMLEgrowth[i,],center=FALSE,cov=VthetaMLE) # Note MD=d^2!
  MDerrorCTBSgrowth[i]=mahalanobis(resCTBSgrowth[i,],center=FALSE,cov=VthetaCTBSrescale) # Note MD=d^2!
  MDerrorCBSMMgrowth[i]=mahalanobis(resCBSMMgrowth[i,],center=FALSE,cov=VthetaCBSrescale) # Note MD=d^2!
}

idMLE=(1:n)[MDerrorMLEgrowth>qchisq(0.975,df=k)]
idCTBS=(1:n)[MDerrorCTBSgrowth>qchisq(0.975,df=k)]
idCBSMM=(1:n)[MDerrorCBSMMgrowth>qchisq(0.975,df=k)]

ymax=sqrt(max(MDerrorMLEgrowth,MDerrorCTBSgrowth,MDerrorCBSMMgrowth))

par(mfrow=c(2,2))
plot(1:n,sqrt(MDerrorMLEgrowth),
     xlab="Index",
     ylab="MD's of MLE residuals",ylim=c(0,ymax))
points(idMLE,sqrt(MDerrorMLEgrowth)[idMLE],pch=16,col=2)
text(idMLE,sqrt(MDerrorMLEgrowth)[idMLE],label=growth[,1][idMLE],pos=1,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)

plot(1:n,sqrt(MDerrorCTBSgrowth),
     xlab="Index",
     ylab="MD's of CTBS residuals",ylim=c(0,ymax))
points(idCTBS,sqrt(MDerrorCTBSgrowth)[idCTBS],pch=16,col=2)
text(idCTBS,sqrt(MDerrorCTBSgrowth)[idCTBS],label=growth[,1][idCTBS],pos=1,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)

plot(1:n,sqrt(MDerrorCBSMMgrowth),
     xlab="Index",
     ylab="MD's of CBS-MM residuals",ylim=c(0,ymax))
points(idCBSMM,sqrt(MDerrorCBSMMgrowth)[idCBSMM],pch=16,col=2)
text(idCBSMM,sqrt(MDerrorCBSMMgrowth)[idCBSMM],label=growth[,1][idCBSMM],pos=1,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)

plot(sqrt(MDerrorMLEgrowth),sqrt(MDerrorCBSMMgrowth),
     xlab="MD's from MLE",
     ylab="MD's CBS-MM")
points(sqrt(MDerrorMLEgrowth)[idMLE],sqrt(MDerrorCBSMMgrowth)[idMLE],pch=16,col=2)
text(sqrt(MDerrorMLEgrowth)[idMLE],sqrt(MDerrorCBSMMgrowth)[idMLE],label=growth[,1][id],pos=1,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)
abline(v=sqrt(qchisq(0.975,df=k)),col=2)
par(mfrow=c(1,1))


#####################################################################
# Figure 4.3 in Demidenko

ymax=max(growth[,3:6])
ymin=min(growth[,3:6])
par(mfrow=c(1,2))
plot(1:4,growth[1,3:6],type="l",xaxt="n",
     ylim=c(ymin,ymax),
     main="Girls",xlab="Age (years)",ylab="Distance (mm)")
axis(1,labels=8:14,at=seq(1,4,by=0.5))
for(i in 2:11){
  lines(1:4,growth[i,3:6],lty=1)
}
lines(1:4,Xg%*%betaMLEgrowthRik,col=2,lwd=2)
lines(1:4,Xg%*%betaCTBSgrowthRESCALE,col=3,lwd=2)
lines(1:4,Xg%*%betaCBSMMgrowthRESCALE,col=4,lwd=2)
lines(1:4,Xg%*%betaCTBSMMgrowthRESCALE,col=6,lwd=2)

plot(1:4,growth[12,3:6],type="l",xaxt="n",
     ylim=c(ymin,ymax),
     main="Boys",xlab="Age (years)",ylab="Distance (mm)")
for(i in 12:27){
  lines(1:4,growth[i,3:6],lty=1)
}
axis(1,labels=8:14,at=seq(1,4,by=0.5))
lines(1:4,Xb%*%betaMLEgrowthRik,col=2,lwd=2)
lines(1:4,Xb%*%betaCTBSgrowthRESCALE,col=3,lwd=2)
lines(1:4,Xb%*%betaCBSMMgrowthRESCALE,col=4,lwd=2)
lines(1:4,Xb%*%betaCTBSMMgrowthRESCALE,col=6,lwd=2)

#lines(1:4,growth[20,3:6],col=2,lwd=2)
#lines(1:4,growth[24,3:6],col=3,lwd=2)
par(mfrow=c(1,1))


#####################################################################
