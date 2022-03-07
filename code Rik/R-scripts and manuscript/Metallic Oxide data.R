library(MASS) # needed to generate from multivariate normal
library(ks)
library(robustbase)
library(nlme)

source("Biweight functions.R")
source("ASYMP NORM constants.R")

metallic=Metal
names(Metal)
n=nrow(metallic)/8
k=8

Ymat=matrix(0,nrow=n,ncol=k)
Type=matrix(0,nrow=n,ncol=1)
for (i in 1:n){
  Ymat[i,]=metallic[(8*(i-1)+1):(8*i),1]
  Type[i]=metallic[8*i,5]
}


#############################################################################
# ESTIMATES REPORTED IN PHD COPT, IN HERITIER AND BY LME
#############################################################################
# MLE, TABLE 8-6 IN PHD COPT
betaMLmetalCopt=c(3.85625,-0.79375)
thetaMLmetalCopt=c(0.565,0.043,0.032,0.043)

# CTBS, TABLE 8-6 IN PHD COPT
betaCTBSmetalCopt=c(3.9318973,-0.4017356)
thetaCTBSmetalCopt=c(0.097,0.012,0.040,0.036)

# CBS-MM TABLE 4.9 IN HERITIER
betaCBSMMmetalHeritier=c(3.726,0.184)
SQRTthetaCBSmetalHeritier=c(0.317,0.144,0.188,0.186)

# FROM THE FUNCTION LME
names(Metal)
options(contrasts=c("contr.treatment","cont.poly"))
metallicML.lme=lme(fixed=Response~Type,
                   data=Metal,
                   random=list(Lots=pdBlocked(list(pdIdent(~1),
                                                      pdIdent(~Sample-1),
                                                      pdIdent(~Chemist-1)))),
                   method="ML")
summary(metallicML.lme)

betaMLmetalLME=fixed.effects(metallicML.lme)
SQRTthetaMLmetalLME=c(0.7451026,0.2394382,0.1426704, 0.2256532)
thetaMLmetalLME=SQRTthetaMLmetalLME^2

########################################################
# MODEL
########################################################
X1=cbind(rep(1,8),rep(0,8))
X2=cbind(rep(1,8),rep(1,8))
q=2

Z1=cbind(rep(1,8))
Z2=cbind(c(rep(1,4),rep(0,4)),c(rep(0,4),rep(1,4)))
Z3=cbind(c(rep(1,2),rep(0,6)),
         c(rep(0,2),rep(1,2),rep(0,4)),
         c(rep(0,4),rep(1,2),rep(0,2)),
         c(rep(0,6),rep(1,2)))

varerror=diag(rep(1,8))

L1=Z1%*%t(Z1)
L2=Z2%*%t(Z2)
L3=Z3%*%t(Z3)
L4=varerror
l=4

#####################################################################
# COMPUTATION OF CTBS/CBS/ML ESTIMATES
#####################################################################
# Setting the breakdown point and cut-off constant for translated biweight
p=0.01
r=0.5
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)

# TO COMPUTE mle
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
  X=ifelse(Type[i]==1,1,0)*X1+ifelse(Type[i]==2,1,0)*X2
  y=as.matrix(Ymat[i,],ncol=1,nrow=5)
  betavectermtranslated=betavectermtranslated+w[i]*t(X)%*%solve(V0)%*%y  
  betamattermtranslated=betamattermtranslated+w[i]*t(X)%*%solve(V0)%*%X 
}

betaoldtranslated=solve(betamattermtranslated)%*%betavectermtranslated
Voldtranslated=V0
thetaoldtranslated=numeric(l)

###########################################################################
# OPTION 2: SOLVING BETA
# setting starting values for beta and theta
Ymat1=Ymat[Type==1,]
Ymat2=Ymat[Type==2,]

# STARTING VALUES FROM ROGKMISS
Vstart1=rogkmiss(as.matrix(Ymat1))
Vstart2=rogkmiss(as.matrix(Ymat2))
Vstart=rogkmiss(as.matrix(Ymat))

mu01=Vstart1$center
mu02=Vstart2$center

X0=rep(1,8)
beta1=solve(t(X0)%*%X0)%*%t(X0)%*%mu01
beta2=solve(t(X0)%*%X0)%*%t(X0)%*%mu02-beta1
betaoldtranslated=c(beta1,beta2)

Voldtranslated=Vstart$cov

thetaoldtranslated=numeric(l)

#################################################
# RUNNING THE ITERATION TO OBTAIN S-ESTIMATES
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
    X=ifelse(Type[i]==1,1,0)*X1+ifelse(Type[i]==2,1,0)*X2
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
    X=ifelse(Type[i]==1,1,0)*X1+ifelse(Type[i]==2,1,0)*X2
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
      Ls=(s==1)*L1+(s==2)*L2+(s==3)*L3+(s==4)*L4
      Utermtranslated[s,i]=
        k*biweightutranslated(sqrt(MDtildetranslated[i]),M,c0t)*
        t(y-X%*%betaoldtranslated)%*%solve(Voldtranslated)%*%
        Ls%*%
        solve(Voldtranslated)%*%(y-X%*%betaoldtranslated)
      Utranslated[s]=Utranslated[s]+Utermtranslated[s,i]
      for (t in 1:l){
        Lt=(t==1)*L1+(t==2)*L2+(t==3)*L3+(t==4)*L4
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
    thetanewtranslated[3]*L3+thetanewtranslated[4]*L4
  
  # End of iteration step
  
  iter=iter+1
  
  diffbeta=max(abs(betanewtranslated-betaoldtranslated))
  difftheta=max(abs(thetanewtranslated-thetaoldtranslated))
  tol=max(diffbeta,difftheta)
  
  # setting starting values for the next iteration within while-loop
  betaoldtranslated=betanewtranslated
  thetaoldtranslated=thetanewtranslated
  Voldtranslated=Vnewtranslated
}

iter
betanewtranslated
thetanewtranslated
sqrt(thetanewtranslated)

# SAVING THE DIFFERENT ESTIMATES AFTER RUNNING THE ITERATION

betaMLmetalRik=betanewtranslated
thetaMLmetalRik=thetanewtranslated

betaCTBSmetalRik=betanewtranslated
thetaCTBSmetalRik=thetanewtranslated

betaCBSmetalRik=betanewtranslated
thetaCBSmetalRik=thetanewtranslated

betaCTBSmetalRESCALE=betanewtranslated
thetaCTBSmetalRESCALE=thetanewtranslated

betaCBSmetalRESCALE=betanewtranslated
thetaCBSmetalRESCALE=thetanewtranslated


###########################################################################
# RESULTS MLE
data.frame(betaMLErik=betaMLmetalRik,
           betaMLEcopt=betaMLmetalCopt,
           betaMLElme=betaMLmetalLME)

data.frame(thetaMLErik=thetaMLmetalRik,
           thetaMLEcopt=thetaMLmetalCopt,
           thetaMLElme=thetaMLmetalLME)

###########################################################################
# CTBS
data.frame(betaCTBSrik=betaCTBSmetalRik,
           betaCTBSCopt=betaCTBSmetalCopt)
data.frame(thetaCTBSrik=thetaCTBSmetalRik,
           thetaCTBSCopt=thetaCTBSmetalCopt)

###########################################################################
# CBS-MM
data.frame(betaCBSrik=betaCBSmetalRik,
           betaCBSMMrik=betaCBSMMmetalRik,
           betaCBSMMheritier=betaCBSMMmetalHeritier)

data.frame(SQRTthetaCBSrik=sqrt(thetaCBSmetalRik),
           SQRTthetaCBSHeritier=SQRTthetaCBSmetalHeritier)

data.frame(thetaCBSrik=thetaCBSmetalRik,
           thetaCBSHeritier=SQRTthetaCBSmetalHeritier^2)

###########################################################################
# CTBS RESCALE
data.frame(betaCTBSrescale=betaCTBSmetalRESCALE,
           betaCTBSCopt=betaCTBSmetalCopt)

data.frame(thetaCTBSrescale=thetaCTBSmetalRESCALE,
           thetaCTBSCopt=thetaCTBSmetalCopt)

###########################################################################
# CBS-MM RESCALE
data.frame(betaCBSrescale=betaCBSmetalRESCALE,
           betaCBSMMrescale=betaCBSMMmetalRESCALE,
           betaCBSMMheritier=betaCBSMMmetalHeritier)

data.frame(SQRTthetaCBSrescale=sqrt(thetaCBSmetalRESCALE),
           SQRTthetaCBSHeritier=SQRTthetaCBSmetalHeritier)


##########################################################################
# STANDARD ERRORS
#########################################################################
# Covariance beta alpha*(E(X^T*V^{-1}*X))^{-1}

#MLE
M=10000
beta=betaMLmetalRik
theta=thetaMLmetalRik
VthetaMLE=theta[1]*L1+theta[2]*L2+theta[3]*L3++theta[4]*L4

# Asymptotic variances betahat
# Expression from our theory

# ESTIMATE FOR X^T*v^{-1}*X
COVestMLErik=matrix(0,nrow=q,ncol=q)
for (i in 1:n){
  X=ifelse(Type[i]==1,1,0)*X1+ifelse(Type[i]==2,1,0)*X2
  COVestMLErik=COVestMLErik+t(X)%*%solve(VthetaMLE)%*%X
}
varbetatranslatedMLE=constbetahattranslated(k,M,c0t)*solve(COVestMLErik/n)

SEbetaMLmetalRik=sqrt(diag(varbetatranslatedMLE)/n)
SEbetaMLmetalCopt=c(0.1826328,0.2820253)
SEbetaMLmetalLME=c(0.1833738,0.2820253)

data.frame(SEbetaMLErik=SEbetaMLmetalRik,
           SEbetaMLEcopt=SEbetaMLmetalCopt,
           SEbetaMLElme=SEbetaMLmetalLME)

#####################################################################
# CTBS with re-scaling
p=0.01
r=0.5
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)

beta=betaCTBSmetalRESCALE
theta=thetaCTBSmetalRESCALE
VthetaCTBSrescale=theta[1]*L1+theta[2]*L2+theta[3]*L3+theta[4]*L4
eigen(VthetaCTBSrescale)

# ESTIMATE FOR X^T*v^{-1}*X
COVestCTBSrescale=matrix(0,nrow=q,ncol=q)
for (i in 1:n){
  X=ifelse(Type[i]==1,1,0)*X1+ifelse(Type[i]==2,1,0)*X2
  COVestCTBSrescale=COVestCTBSrescale+t(X)%*%solve(VthetaCTBSrescale)%*%X
}
varbetatranslatedCTBSrescale=constbetahattranslated(k,M,c0t)*solve(COVestCTBSrescale/n)

SEbetaCTBSmetalRESCALE=sqrt(diag(varbetatranslatedCTBSrescale)/n)
SEbetaCTBSmetalCopt=c(0.09964693,0.15387677)

data.frame(SEbetaCTBSrescale=SEbetaCTBSmetalRESCALE,
           SEbetaCTBSCopt=SEbetaCTBSmetalCopt)


############################################################
# CBS-MM with re-scaling
eff=0.95
objectf=function(c){
  1/constbetahat(k,c)-eff 
}
c1=uniroot(objectf,interval=c(0.001,100))$root

p=0.01
c0t=rhotranslatedconst(k,r,p,0.01,5)
M=sqrt(qchisq(1-p,df=k))-c0t
b0t=expecrhotranslated(k,M,c0t)
M=0
c0t=c1

beta=betaCBSMMmetalRESCALE
theta=thetaCBSmetalRESCALE
VthetaCBSrescale=theta[1]*L1+theta[2]*L2+theta[3]*L3+theta[4]*L4

# ESTIMATE FOR X^T*v^{-1}*X
COVestCBSMMrescale=matrix(0,nrow=q,ncol=q)
for (i in 1:n){
  X=ifelse(Type[i]==1,1,0)*X1+ifelse(Type[i]==2,1,0)*X2
  COVestCBSMMrescale=COVestCBSMMrescale+t(X)%*%solve(VthetaCBSrescale)%*%X
}
varbetatranslatedCBSMMrescale=constbetahattranslated(k,M,c0t)*solve(COVestCBSMMrescale/n)

SEbetaCBSMMmetalRESCALE=sqrt(diag(varbetatranslatedCBSMMrescale)/n)
SEbetaCBSMMmetalHeritier=c(0.066,0.066)

data.frame(SEbetaCBSMMrescale=SEbetaCBSMMmetalRESCALE,
           SEbetaCBSMMHeritier=SEbetaCBSMMmetalHeritier)

# SE for Heritier from out results
sqrt(diag(solve(B)%*%varbetatranslatedCBSMMrescale%*%t(solve(B)))/n)


############################################################
# WALD TEST
###########################################################
# MLE
waldMLE=betaMLmetalRik/SEbetaMLmetalRik
pMLE=2*(1-pnorm(abs(waldMLE)))

waldMLEcopt=betaMLmetalCopt/SEbetaMLmetalCopt
pMLEcopt=2*(1-pnorm(abs(waldMLEcopt)))

data.frame(MLE=betaMLmetalRik,
           SE.MLE=round(SEbetaMLmetalRik,3),
           TMLE=waldMLE,
           PvalMLE=round(pMLE,4),
           PvalMLCopt=round(pMLEcopt,3))

# CTBS
waldCTBS=betaCTBSmetalRESCALE/SEbetaCTBSmetalRESCALE
pCTBS=2*(1-pnorm(abs(waldCTBS)))

waldCTBScopt=betaCTBSmetalCopt/SEbetaCTBSmetalCopt
pCTBScopt=2*(1-pnorm(abs(waldCTBScopt)))

data.frame(CTBS=betaCTBSmetalRESCALE,
           SE.CTBS=round(SEbetaCTBSmetalRESCALE,3),
           TCTBS=waldCTBS,
           PvalCTBS=round(pCTBS,4),
           PvalCTBScopt=round(pCTBScopt,4))

# CBS-MM
waldCBSMM=betaCBSMMmetalRESCALE/SEbetaCBSMMmetalRESCALE
pCBSMM=2*(1-pnorm(abs(waldCBSMM)))

waldCBSMMheritier=betaCBSMMmetalHeritier/SEbetaCBSMMmetalHeritier
pCBSMMheritier=2*(1-pnorm(abs(waldCBSMMheritier)))


data.frame(CBSMM=betaCBSMMmetalRESCALE,
           SE.CBSMM=round(SEbetaCBSMMmetalRESCALE,3),
           TCBSMM=waldCBSMM,
           PvalCBS=round(pCBSMM,4),
           PvalCBSheritier=round(pCBSMMheritier,4))

# CBS-MM comprable with Heritier
betaCBSMMmetalRESCALEadapted=solve(B)%*%betaCBSMMmetalRESCALE
SEbetaCBSmetalRESCALEadapted=sqrt(diag(solve(B)%*%varbetatranslatedCBSMMrescale%*%t(solve(B)))/n)

waldCBSMMadapted=betaCBSMMmetalRESCALEadapted/SEbetaCBSmetalRESCALEadapted
pCBSMMadapted=2*(1-pnorm(abs(waldCBSMMadapted)))

data.frame(CBSMM=betaCBSMMmetalRESCALEadapted,
           SE.CBSMM=round(SEbetaCBSmetalRESCALEadapted,3),
           TCBSMM=waldCBSMMadapted,
           PvalCBSMM=round(pCBSMMadapted,4))


############################################################################################
# RESISDUALS
############################################################################################
YMLmat=matrix(0,nrow=n,ncol=k)
YCTBSmat=matrix(0,nrow=n,ncol=k)
YCBSMMmat=matrix(0,nrow=n,ncol=k)

resMLE=matrix(0,nrow=n,ncol=k)
resCTBS=matrix(0,nrow=n,ncol=k)
resCBSMM=matrix(0,nrow=n,ncol=k)

for (i in 1:n){
  X=ifelse(Type[i]==1,1,0)*X1+ifelse(Type[i]==2,1,0)*X2
  YMLmat[i,]=X%*%betaMLmetalRik
  YCTBSmat[i,]=X%*%betaCTBSmetalRESCALE
  YCBSMMmat[i,]=X%*%betaCBSMMmetalRESCALE
}
resMLE=Ymat-YMLmat
resCTBS=Ymat-YCTBSmat
resCBSMM=Ymat-YCBSMMmat

cbind(resMLE[c(24,25,30),],
resCTBS[c(24,25,30),])
############################################################################################
# MAHALANOBIS DISTANCES OF RESIDUALS
a=0.025

MDerrorMLmetal=numeric(n)
MDerrorCTBSmetal=numeric(n)
MDerrorCBSMMmetal=numeric(n)

theta=thetaMLmetalRik
VthetaMLrescale=theta[1]*L1+theta[2]*L2+theta[3]*L3++theta[4]*L4

theta=thetaCTBSmetalRESCALE
VthetaCTBSrescale=theta[1]*L1+theta[2]*L2+theta[3]*L3++theta[4]*L4

theta=thetaCBSmetalRESCALE
VthetaCBSrescale=theta[1]*L1+theta[2]*L2+theta[3]*L3++theta[4]*L4


for (i in 1:n){
  MDerrorMLmetal[i]=mahalanobis(resMLE[i,],center=FALSE,cov=VthetaMLrescale) # Note MD=d^2!
  MDerrorCTBSmetal[i]=mahalanobis(resCTBS[i,],center=FALSE,cov=VthetaCTBSrescale) # Note MD=d^2!
  MDerrorCBSMMmetal[i]=mahalanobis(resCBSMM[i,],center=FALSE,cov=VthetaCBSrescale) # Note MD=d^2!
}

idCTBSmetal=(1:n)[MDerrorCTBSmetal>qchisq(0.975,df=k)]
idCBSMMmetal=(1:n)[MDerrorCBSMMmetal>qchisq(0.975,df=k)]

par(mfrow=c(1,2))
plot(sqrt(MDerrorMLmetal),sqrt(MDerrorCTBSmetal),
     xlab="MD's of ML residuals",
     ylab="MD's of CTBS residuals",
     main="CTBS residuals")
points(sqrt(MDerrorMLmetal)[idCTBSmetal],sqrt(MDerrorCTBSmetal)[idCTBSmetal],pch=16,col=2)
text(sqrt(MDerrorMLmetal)[idCTBSmetal],sqrt(MDerrorCTBSmetal)[idCTBSmetal],labels=idCTBSmetal,pos=1,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)
abline(v=sqrt(qchisq(0.975,df=k)),col=2)


plot(sqrt(MDerrorMLmetal),sqrt(MDerrorCBSMMmetal),
     xlab="MD's of ML residuals",
     ylab="MD's of CBS-MM residuals",
     main="CBS-MM residuals")
points(sqrt(MDerrorMLmetal)[idCBSMMmetal],sqrt(MDerrorCBSMMmetal)[idCBSMMmetal],pch=16,col=2)
text(sqrt(MDerrorMLmetal)[idCBSMMmetal],sqrt(MDerrorCBSMMmetal)[idCBSMMmetal],labels=idCTBSmetal,pos=1,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)
abline(v=sqrt(qchisq(0.975,df=k)),col=2)
par(mfrow=c(1,1))


ymax=max(sqrt(MDerrorCTBSmetal))
par(mfrow=c(1,2))
plot((1:n)[Type==1],sqrt(MDerrorCTBSmetal)[Type==1],
     ylim=c(0,ymax),
     xlab="Index",
     ylab="MD's of CTBS residuals",
     main="Type 1")
points(idCTBSmetal,sqrt(MDerrorCTBSmetal)[idCTBSmetal],pch=16,col=2)
text(idCTBSmetal,sqrt(MDerrorCTBSmetal)[idCTBSmetal],label=(1:n)[idCTBSmetal],pos=1,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)

plot((1:n)[Type==2],sqrt(MDerrorCTBSmetal)[Type==2],
     ylim=c(0,ymax),
     xlab="Index",
     ylab="MD's of CTBS residuals",
     main="Type 2")
points(idCTBSmetal,sqrt(MDerrorCTBSmetal)[idCTBSmetal],pch=16,col=2)
text(idCTBSmetal,sqrt(MDerrorCTBSmetal)[idCTBSmetal],label=(1:n)[idCTBSmetal],pos=1,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)
par(mfrow=c(1,1))

###############################################################
names(Metal)

typevec=ifelse(Metal[,5]=="type1",1,3)
colvec=rep(1,length(typevec))
colvec[c(193:200)]=2
colvec[c(185:192)]=3

coplot(Lots~Response|Sample+Chemist,data=Metal,pch=typevec,col=colvec,show.given = TRUE)
