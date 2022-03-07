library(MASS) # needed to generate from multivariate normal
library(ks)
library(robustbase)
library(pracma)
library(nlme)

source("Biweight functions.R")
source("ASYMP NORM constants.R")

# EXAMPLE DISCUSSED IN phd Heritier and HERITIER, CANTONI, Heritier & VICTORIA-FESER
# Data are from Berry (1987)
resistance=read.csv("resistance.csv", header=TRUE,sep=";")

#############################################################################
# FIGURE 4.1 IN HERITIER
#############################################################################
ymax=max(resistance[,2:6])
ymin=0
plot(1:5,resistance[1,2:6],type="l",ylim=c(ymin,ymax),
     xlab="Electrode type",ylab="Resistance",
     xlim=c(1,6),xaxt="n")
for (i in 2:16){
  lines(1:5,resistance[i,2:6],lty=i)
}
axis(1,at=1:5,labels=c("E1","E2","E3","E4","E5"))
legend("topright",legend=1:16,lty=1:16,title="Subject",bty="n",cex=0.75)

#############################################################################
# ESTIMATES REPORTED IN COPT AND HERITIER
#############################################################################
# MLE IN COPT TABLE 8.1
betaMLEresistanceCopt=c(2.03,-0.21,0.84,0.55,-0.53)
thetaMLEresistanceCopt=c(1.329,2.098)

# Results from function lme with "sum to zero" contrasts
options(contrasts=c("contr.sum","cont.poly"))
electrode=Electrode
electrode[,1]=electrode[,1]/100
names(electrode)
resistance.lme=lme(fixed=Resis~Elec,
                   data=electrode,
                   random=list(Subject=pdIdent(~1)),
                   method="ML")
betaMLEresistanceLME=fixed.effects(resistance.lme)
summary(resistance.lme)
thetaMLEresistanceLME=c(1.152971, 1.448449)^2

#############################################################################
# REML IN TABLE 4.1 IN HERITIER
betaREMLresistanceHeritier=c(2.030,-0.213,0.842,0.549,-0.526)
SQRTthetaREMLresistanceHeritier=c(1.190,1.459)

options(contrasts=c("contr.sum","cont.poly"))
electrode=Electrode
electrode[,1]=electrode[,1]/100
names(electrode)
resistance.lme=lme(fixed=Resis~Elec,
                   data=electrode,
                   random=list(Subject=pdIdent(~1)),
                   method="REML")

betaREMLresistanceLME=fixed.effects(resistance.lme)
summary(resistance.lme)
SQRTthetaREMLresistanceLME=c(1.190784, 1.495952)

data.frame(betaREMLlme=betaREMLresistanceLME,
           betaREMLheritier=betaREMLresistanceHeritier)

data.frame(SQRTthetaREMLlme=SQRTthetaREMLresistanceLME,
           SQRTthetaREMLheritier=SQRTthetaMLEresistanceHeritier)

#############################################################################
# CTBS IN TABLE 8-1 IN PHD COPT
betaCTBSresistanceCopt=c(1.628,-0.068,0.512,0.142,-0.158)
thetaCTBSresistanceCopt=c(1.087,0.711)

#############################################################################
# CBS-MM IN TABLE 4.1 IN HERITIER
betaCBSMMresistanceHeritier=c(1.440, -0.161,  0.403,  0.243, -0.169)
SQRTthetaCBSresistanceHeritier=c(0.842,0.761)

###################################################################
# Model
Ymat=as.matrix(resistance[,2:6])/100 # scaling mentioned on page 98
n=nrow(Ymat)
k=ncol(Ymat)

# SUM TO ZERO CONTRAST
X=cbind(rep(1,5),contr.sum(5))
# TREATMENT CONTRAST
#X=cbind(rep(1,5),contr.treatment(5))
q=5

Z=rep(1,5)
varerror=diag(rep(1,k))
L1=Z%*%t(Z)
L2=varerror
l=2

#####################################################################
# COMPUTATION OF CTBS/CBS/MLE ESTIMATES
#####################################################################
# Setting the breakdown point and cut-off constant for translated biweight
r=0.5
p=0.01
c0t=rhotranslatedconst(k,r,p,0.01,10)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)

# TO COMPUTE MLE
#M=10000

# TO COMPUTE CBS
c0=rhoconst(k,r,0.01,100)
b0=expecrho(k,c0)
M=0
c0t=c0

#################################################################
# SETTING STARTING VALUES
###########################################################################

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

# OPTION 3: location-covariance from rogkmiss
Vstart=rogkmiss(Ymat)
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
betavectermtranslated=matrix(0,nrow=q,ncol=1) # vector-term in step 5
betamattermtranslated=matrix(0,nrow=q,ncol=q) # matrix-term in step 5
for (i in 1:n){
  y=matrix(Ymat[i,],ncol=1,nrow=5)
  betavectermtranslated=betavectermtranslated+w[i]*t(X)%*%solve(V0)%*%y  
  betamattermtranslated=betamattermtranslated+w[i]*t(X)%*%solve(V0)%*%X 
}

betaoldtranslated=solve(betamattermtranslated)%*%betavectermtranslated
Voldtranslated=V0
thetaoldtranslated=numeric(l)

#################################################
# OPTION 2 SIMPLE X MATRIX

betaoldtranslated=solve(t(X)%*%X)%*%t(X)%*%mu0

thetaoldtranslated=numeric(l)
thetaoldtranslated[1]=(sum(V0)-sum(diag(V0)))/(sum(Z%*%t(Z))-sum(diag(Z%*%t(Z))))
thetaoldtranslated[2]=(sum(diag(V0))-thetaoldtranslated[1]*sum(diag(Z%*%t(Z))))/sum(diag(varerror))
Voldtranslated=thetaoldtranslated[1]*L1+thetaoldtranslated[2]*L2

################################################################################
# RUNNING THE ITERATION TO OBTAIN S-ESTIMATES
################################################################################
# stop criterion
tol=1
iter=0
n=nrow(Ymat)

while (tol>=1e-06){
  ##############################################################################
  # Iteration step
  #############################################
  
  # Re-scaling Mahalanobis distances to satisfy S-constraint
  # array for n Mahalanobis distances
  MD=numeric(n)
  for (i in 1:n){
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
  
  # OPTION 3: Rocke's adaptation as in Heritier
  # CORRECTION USED IN PHD Heritier? SEE (4.14)
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
    y=matrix(Ymat[i,],ncol=1,nrow=5)
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
      Ls=(s==1)*L1+(s==2)*L2
      Utermtranslated[s,i]=
        k*biweightutranslated(sqrt(MDtildetranslated[i]),M,c0t)*
        t(y-X%*%betaoldtranslated)%*%solve(Voldtranslated)%*%
        Ls%*%
        solve(Voldtranslated)%*%(y-X%*%betaoldtranslated)
      Utranslated[s]=Utranslated[s]+Utermtranslated[s,i]
      for (t in 1:l){
        Lt=(t==1)*L1+(t==2)*L2
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
  Vnewtranslated=thetanewtranslated[1]*L1+thetanewtranslated[2]*L2
  
  # End of iteration step
  
  iter=iter+1
  
  diffbeta=norm(betanewtranslated-betaoldtranslated,type="F")
  difftheta=norm(thetanewtranslated-thetaoldtranslated,type="F")
  tol=max(diffbeta,difftheta)
  
  # setting starting values for the next iteration within while-loop
  betaoldtranslated=betanewtranslated
  thetaoldtranslated=thetanewtranslated
  Voldtranslated=Vnewtranslated
}
################################################################################

tol
iter
betanewtranslated
thetanewtranslated
sqrt(thetanewtranslated)

# SAVING THE DIFFERENT ESTIMATES AFTER RUNNING THE ITERATION

betaMLEresistanceRik=betanewtranslated
thetaMLEresistanceRik=thetanewtranslated

#betaCTBSresistanceRik=betanewtranslated
#thetaCTBSresistanceRik=thetanewtranslated

#betaCBSresistanceRik=betanewtranslated
#thetaCBSresistanceRik=thetanewtranslated

#betaCBSMMresistanceRik=betanewtranslated

#betaCTBSresistanceRESCALE=betanewtranslated
#thetaCTBSresistanceRESCALE=thetanewtranslated

betaCBSresistanceRESCALE=betanewtranslated
thetaCBSresistanceRESCALE=thetanewtranslated



# RESULTS MLE
data.frame(betaMLErik=betaMLEresistanceRik,
           betaRMLEHeritier=betaREMLresistanceLME,
           betaMLECopt=betaMLEresistanceLME)

data.frame(thetaMLErik=thetaMLEresistanceRik,
           thetaMLECopt=thetaMLEresistanceLME)

#RESULTS REML
data.frame(SQRTthetaMLErik=sqrt(thetaMLEresistanceRik),
           SQRTthetaRMLEHeritier=SQRTthetaREMLresistanceLME)

# RESULTS CTBS
data.frame(betaCTBSrik=betaCTBSresistanceRik,
           betaCTBScopt=betaCTBSresistanceCopt)

data.frame(thetaCTBSrik=thetaCTBSresistanceRik,
           thetaCTBSCopt=thetaCTBSresistanceCopt)

# RESULTS CBS-MM
data.frame(betaCBSrik=betaCBSresistanceRik,
           betaCBSMMrik=betaCBSMMresistanceRik,
           betaCBSMMHeritier=betaCBSMMresistanceHeritier)

data.frame(SQRTthetaCBSrik=sqrt(thetaCBSresistanceRik),
           SQRTthetaCBSHeritier=SQRTthetaCBSresistanceHeritier)

# RESULTS CTBS RE-SCALED
data.frame(betaCTBSrescale=betaCTBSresistanceRESCALE,
           betaCTBScopt=betaCTBSresistanceCopt)

data.frame(thetaCTBSrescale=thetaCTBSresistanceRESCALE,
           thetaCTBSCopt=thetaCTBSresistanceCopt)

# RESULTS CBS-MM RE-SCALED
data.frame(betaCBSrescale=betaCBSresistanceRESCALE,
           betaCBSMMrescale=betaCBSMMresistanceRESCALE,
           betaCBSMMHeritier=betaCBSMMresistanceHeritier)

data.frame(SQRTthetaCBSrescale=sqrt(thetaCBSresistanceRESCALE),
           SQRTthetaCBSHeritier=SQRTthetaCBSresistanceHeritier)

##########################################################################
# DTERMINANTS
#########################################################################
beta=betaCBSresistanceRESCALE
theta=thetaCBSresistanceRESCALE

beta=betaCBSresistanceRik
theta=SQRTthetaCBSresistanceHeritier^2

beta=betaCTBSresistanceCopt
theta=thetaCTBSresistanceCopt

MD=numeric(n)
for (i in 1:n){
  y=Ymat[i,]
  MD[i]=mahalanobis(y,center=X%*%beta,cov=theta[1]*L1+theta[2]*L2) # Note MD=d^2!
}

# determining scaling constant for MD
# for translated biweight
objfuntranslated=function(s){
  mean(biweightrhotranslated(sqrt(MD)/s,M,c0t))-expecrhotranslated(k,M,c0t)
}
MDscaletranslated=uniroot(f=objfuntranslated,c(0.01,100))$root
theta=theta*MDscaletranslated^2
#########################################################################

# CTBS-RE-SCALED
detCTBSresistanceRESCALE=det(theta[1]*L1+theta[2]*L2)
detCTBSresistanceCopt=det(theta[1]*L1+theta[2]*L2)

data.frame(detCTBSresistanceRESCALE=detCTBSresistanceRESCALE,
           detCTBSresistanceCopt=detCTBSresistanceCopt)

# CBS-MM RE-SCALED
detCBSresistanceRESCALE=det(theta[1]*L1+theta[2]*L2)
detCBSresistanceHeritier=det(theta[1]*L1+theta[2]*L2)

data.frame(detCBSresistanceRik=detCBSresistanceRESCALE,
           detCBSresistanceHeritier=detCBSresistanceHeritier)


##########################################################################
# STANDARD ERRORS
#########################################################################
# MLE
M=10000
beta=betaMLEresistanceRik
theta=thetaMLEresistanceRik
Vtheta=theta[1]*L1+theta[2]*L2
L=as.matrix(cbind(vec(L1),vec(L2)))

# Asymptotic variances betahat
# Expression from our theory
varbetatranslatedMLE=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%solve(Vtheta)%*%X)
# Expression from Heritier & Victoria-Feser (2006)
varbetaHeritiertranslatedMLE=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%X) %*% t(X) %*% Vtheta %*% X %*% solve(t(X)%*%X)

round(varbetatranslatedMLE,5)
round(varbetaHeritiertranslatedMLE,5)

SEbetaMLEresistanceRik=sqrt(diag(varbetatranslatedMLE)/n)
SEbetaRMLEresistanceHeritier=c(0.341,0.334,0.334,0.334,0.334)
SEbetaMLEresistanceCopt=c(0.3306,0.3239,0.3239,0.3239,0.3239)

data.frame(SEbetaMLErik=SEbetaMLEresistanceRik,
           SEbetaRMLEHeritier=SEbetaRMLEresistanceHeritier,
           SEbetaMLEcopt=SEbetaMLEresistanceCopt)

#########################################################################
# CTBS RE-SCALED
p=0.01
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)

beta=betaCTBSresistanceRESCALE
theta=thetaCTBSresistanceRESCALE
Vtheta=theta[1]*L1+theta[2]*L2
L=as.matrix(cbind(vec(L1),vec(L2)))

# Asymptotic variances betahat
# Expression from our theory
varbetatranslatedCTBS=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%solve(Vtheta)%*%X)
# Expression from Heritier & Victoria-Feser (2006)
varbetaHeritiertranslatedCTBS=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%X) %*% t(X) %*% Vtheta %*% X %*% solve(t(X)%*%X)

round(varbetatranslatedCTBS,5)
round(varbetaHeritiertranslatedCTBS,5)

SEbetaCTBSresistanceRESCALE=sqrt(diag(varbetatranslatedCTBS)/n)
SEbetaCTBSresistanceCopt=c(0.2638,0.2553,0.2553,0.2553,0.2553)

data.frame(SEbetCTBSrescale=SEbetaCTBSresistanceRESCALE,
           SEbetaCTBScopt=SEbetaCTBSresistanceCopt)

#########################################################################
# CBS-MM RE-SCALED
objectf=function(c){
  1/constbetahat(k,c)-0.95 
}
c1=uniroot(f=objectf,lower=1,upper=10)$root
p=0.01
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)
M=0
c0t=c1

beta=betaCBSMMresistanceRESCALE
theta=thetaCBSresistanceRESCALE
Vtheta=theta[1]*L1+theta[2]*L2
L=as.matrix(cbind(vec(L1),vec(L2)))

# Asymptotic variances betahat
# Expression from our theory
varbetatranslatedCBSMM=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%solve(Vtheta)%*%X)
# Expression from Heritier & Victoria-Feser (2006)
varbetaHeritiertranslatedCBSMM=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%X) %*% t(X) %*% Vtheta %*% X %*% solve(t(X)%*%X)

round(varbetatranslatedCBSMM,5)
round(varbetaHeritiertranslatedCBSMM,5)

SEbetaCBSMMresistanceRESCALE=sqrt(diag(varbetatranslatedCBSMM)/n)
SEbetaCBSMMresistanceHeritier=c(0.233,0.175,0.175,0.175,0.175)

data.frame(SEbetaCBSMMrescale=SEbetaCBSMMresistanceRESCALE,
           SEbetaCBSMMheritier=SEbetaCBSMMresistanceHeritier)


#########################################################################
# CTBSMM RESCALED
objectf=function(c){
  1/constbetahat(k,c)-0.95 
}
c1=uniroot(f=objectf,lower=1,upper=10)$root
p=0.01
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)
M=0
c0t=c1

beta=betanewtranslated
theta=thetanewtranslated
Vtheta=theta[1]*L1+theta[2]*L2
L=as.matrix(cbind(vec(L1),vec(L2)))

# Asymptotic variances betahat
# Expression from our theory
varbetatranslatedCTBSMM=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%solve(Vtheta)%*%X)
# Expression from Heritier & Victoria-Feser (2006)
varbetaHeritiertranslatedCTBSMM=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%X) %*% t(X) %*% Vtheta %*% X %*% solve(t(X)%*%X)

round(varbetatranslatedCTBSMM,5)
round(varbetaHeritiertranslatedCTBSMM,5)

SEbetaCTBSMMresistanceRESCALE=sqrt(diag(varbetatranslatedCTBSMM)/n)
SEbetaCBSMMresistanceHeritier=c(0.233,0.175,0.175,0.175,0.175)

data.frame(SEbetaCTBSMMrescale=SEbetaCTBSMMresistanceRESCALE,
           SEbetaCBSMMheritier=SEbetaCBSMMresistanceHeritier)


############################################################
# WALD TESTS
###########################################################
# MLE
waldMLE=betaMLEresistanceRik/SEbetaMLEresistanceRik
pMLE=2*(1-pnorm(abs(waldMLE)))
data.frame(MLE=betaMLEresistanceRik,
           SE.MLE=round(SEbetaMLEresistanceRik,3),
           TMLE=waldMLE,
           PvalMLE=round(pMLE,4))

# CTBS RE-SCALED
waldCTBS=betaCTBSresistanceRESCALE/SEbetaCTBSresistanceRESCALE
pCTBS=2*(1-pnorm(abs(waldCTBS)))
data.frame(CTBS=betaCTBSresistanceRESCALE,
           SE.CTBS=round(SEbetaCTBSresistanceRESCALE,3),
           TCTBS=waldCTBS,
           PvalCTBS=round(pCTBS,4))

# CBS-MM RESCALED
waldCBSMM=betaCBSMMresistanceRESCALE/SEbetaCBSMMresistanceRESCALE
pCBSMM=2*(1-pnorm(abs(waldCBSMM)))
data.frame(CBSMM=betaCBSMMresistanceRESCALE,
           SE.CBSMM=round(SEbetaCBSMMresistanceRESCALE,3),
           TCBSMM=waldCBSMM,
           PvalCBSMM=round(pCBSMM,4))


# CTBS-MM RESCALE
waldCTBSMM=betaCTBSMMresistanceRESCALE/SEbetaCTBSMMresistanceRESCALE
pCTBSMM=2*(1-pnorm(abs(waldCTBSMM)))
data.frame(CTBSMM=betaCTBSMMresistanceRESCALE,
           SE.CTBSMM=round(SEbetaCTBSMMresistanceRESCALE,3),
           TCTBSMM=waldCTBSMM,
           PvalCTBSMM=round(pCTBSMM,4))



############################################################################################
# RESISDUALS
############################################################################################
YMLEresistanceMat=matrix(0,nrow=n,ncol=k)
YCTBSresistanceMat=matrix(0,nrow=n,ncol=k)
YCBSMMresistanceMat=matrix(0,nrow=n,ncol=k)

for (i in 1:n){
  YMLEresistanceMat[i,]=X%*%betaMLEresistanceRik
  YCTBSresistanceMat[i,]=X%*%betaCTBSresistanceRESCALE
  YCBSMMresistanceMat[i,]=X%*%betaCBSMMresistanceRESCALE
}

resMLEresistance=Ymat-YMLEresistanceMat
resCTBSresistance=Ymat-YCTBSresistanceMat
resCBSMMresistance=Ymat-YCBSMMresistanceMat

# MAHALANOBIS DISTANCES OF RESIDUALS
a=0.025

MDerrorMLEresistance=numeric(n)
MDerrorCTBSresistance=numeric(n)
MDerrorCBSMMresistance=numeric(n)

VMLEresistance=thetaMLEresistanceRik[1]*L1+thetaMLEresistanceRik[2]*L2
VCTBSresistance=thetaCTBSresistanceRESCALE[1]*L1+thetaCTBSresistanceRESCALE[2]*L2
VCBSresistance=thetaCBSresistanceRESCALE[1]*L1+thetaCBSresistanceRESCALE[2]*L2

for (i in 1:n){
  MDerrorMLEresistance[i]=mahalanobis(resMLEresistance[i,],center=FALSE,cov=VMLEresistance) # Note MD=d^2!
  MDerrorCTBSresistance[i]=mahalanobis(resCTBSresistance[i,],center=FALSE,cov=VCTBSresistance) # Note MD=d^2!
  MDerrorCBSMMresistance[i]=mahalanobis(resCBSMMresistance[i,],center=FALSE,cov=VCBSresistance) # Note MD=d^2!
}

par(mfrow=c(2,2))
plot(sqrt(MDerrorMLEresistance),sqrt(MDerrorCTBSresistance),
     xlab="MD's of MLE residuals",
     ylab="MD's of CTBS residuals",
     main="Mahalanobis distances for Resistance data")
IDCTBSresistance=(1:n)[MDerrorCTBSresistance>qchisq(0.975,df=k)]

points(sqrt(MDerrorMLEresistance)[IDCTBSresistance],
       sqrt(MDerrorCTBSresistance)[IDCTBSresistance],pch=16,col=2)
text(sqrt(MDerrorMLEresistance)[IDCTBSresistance],
     sqrt(MDerrorCTBSresistance)[IDCTBSresistance],label=IDCTBSresistance,pos=2,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)
abline(v=sqrt(qchisq(0.975,df=k)),col=2)

plot(sqrt(MDerrorMLEresistance),sqrt(MDerrorCBSMMresistance),
     xlab="MD's of MLE residuals",
     ylab="MD's of CBS-MM residuals",
     main="Mahalanobis distances for Resistance data")
IDCBSMMresistance=(1:n)[MDerrorCBSMMresistance>qchisq(0.975,df=k)]

points(sqrt(MDerrorMLEresistance)[IDCBSMMresistance],
       sqrt(MDerrorCBSMMresistance)[IDCBSMMresistance],pch=16,col=2)
text(sqrt(MDerrorMLEresistance)[IDCBSMMresistance],
     sqrt(MDerrorCBSMMresistance)[IDCBSMMresistance],label=IDCBSMMresistance,pos=2,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)
abline(v=sqrt(qchisq(0.975,df=k)),col=2)

plot(sqrt(MDerrorCBSMMresistance),sqrt(MDerrorCTBSresistance),
     xlab="MD's of CBS-MM residuals",
     ylab="MD's of CTBS residuals",
     main="Mahalanobis distances for Resistance data")
abline(b=1,a=0)
par(mfrow=c(1,1))

# DATA REVISITED
ymax=max(resistance[,2:6])
ymin=0
plot(1:5,resistance[1,2:6],type="l",ylim=c(ymin,ymax),
     xlab="Electrode type",ylab="Resistance",
     main="Skin Resistance Data with identified outlying observations",
     xlim=c(1,6),xaxt="n")
for (i in 2:16){
  lines(1:5,resistance[i,2:6],lty=i)
}
axis(1,at=1:5,labels=c("E1","E2","E3","E4","E5"))
lines(1:5,resistance[15,2:6],lwd=2,col=2)
lines(1:5,resistance[2,2:6],lwd=2,col=3)
lines(1:5,resistance[1,2:6],lwd=2,col=4)
lines(1:5,resistance[5,2:6],lwd=2,col=5)
legend("topright",legend=1:16,
       lty=c(1,1,3,4,1,6:14,1,16),
       lwd=c(2,2,1,1,2,1,1,1,1,1,1,1,1,1,2,1),
       col=c(4,3,1,1,5,1,1,1,1,1,1,1,1,1,2,1),
       title="Subject",bty="n")

###############################################################
# SIGNIFICANCE
ymax=max(Ymat)
YmatCleaned=Ymat[-c(2,15),]
par(mfrow=c(1,2))
boxplot(Ymat[,1],Ymat[,2],Ymat[,3],Ymat[,4],Ymat[,5],
        ylim=c(0,ymax),
        names=c("E1","E2","E3","E4","E5"),
        xlab="Electrode type",
        ylab="Resistance",
        main="Resistance data")
abline(h=betaCBSMMresistanceRik[1])
boxplot(YmatCleaned[,1],YmatCleaned[,2],YmatCleaned[,3],
        YmatCleaned[,4],YmatCleaned[,5],
        ylim=c(0,ymax),
        names=c("E1","E2","E3","E4","E5"),
        xlab="Electrode type",
        ylab="Resistance",
        main="Resistance data without 2 and 15")
abline(h=betaCBSMMresistanceRik[1])
par(mfrow=c(1,1))

#####################################################################
# SIMPLE 1-WAY ANOVA
tmp=vec(t(resistance[-c(2,15),])[-1,])
fac=as.factor(names(tmp))
dataset=data.frame(type=fac,response=tmp)

summary(aov(response~type,data = dataset))

