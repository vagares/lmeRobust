library(MASS) # needed to generate from multivariate normal
library(ks)
library(robustbase)
library(pracma)
library(nlme)


source("Biweight functions.R")
source("ASYMP NORM constants.R")

semantic=Semantic
n=nlevels(semantic[,1])
k=6

# Renumbering subjects
subject=NULL
for (i in 1:n){
  subject=c(subject,rep(i,6))
}
semantic[,1]=subject

Ymat=matrix(0,nrow=n,ncol=k)
for (i in 1:n){
  Ymat[i,]=semantic[,4][semantic[,1]==i]
}
apply(Ymat,MARGIN=2,FUN=mean)
mean(Ymat)

# REARRANGING TO ORDER y11,y12,y13,y21,y22,y23
P=cbind(c(0,0,0,1,0,0),
        c(0,1,0,0,0,0),
        c(0,0,0,0,0,1),
        c(0,0,1,0,0,0),
        c(1,0,0,0,0,0),
        c(0,0,0,0,1,0))

###################################################################################
# DISPLAY OF DATA
###################################################################################

semanticSHORT=semantic[semantic[,2]=="short",]
semanticLONG=semantic[semantic[,2]=="long",]
semanticREL=semantic[semantic[,3]=="related",]
semanticNEU=semantic[semantic[,3]=="neutral",]
semanticUNR=semantic[semantic[,3]=="unrelated",]

Yshort=matrix(0,nrow=n,ncol=3)
Ylong=matrix(0,nrow=n,ncol=3)
Yrel=matrix(0,nrow=n,ncol=2)
Yneu=matrix(0,nrow=n,ncol=2)
Yunr=matrix(0,nrow=n,ncol=2)
for (i in 1:n){
  Yshort[i,]=semanticSHORT[semanticSHORT[,1]==i,4]
  Ylong[i,]=semanticLONG[semanticLONG[,1]==i,4]
  Yrel[i,]=semanticREL[semanticREL[,1]==i,4]
  Yneu[i,]=semanticNEU[semanticNEU[,1]==i,4]
  Yunr[i,]=semanticUNR[semanticUNR[,1]==i,4]
}


par(mfrow=c(2,2))
ymax=max((Yshort+Ylong)/2)
ymin=min((Yshort+Ylong)/2)
plot(1:3,(Yshort[1,]+Ylong[1,])/2,type="l",
     xlim=c(1,3.25),
     ylim=c(ymin,ymax),xaxt="n",
     xlab="Factor CONDITION",
     ylab="Mean Response",
     main="Mean response over CONDITION")
axis(1,at=1:3,labels = c("related","neutral","unrelated"))
for (i in 1:n){
  lines(1:3,(Yshort[i,]+Ylong[i,])/2,lty=i)
}
legend("topright",legend=1:n,lty=1:n,
       title="Subject",bty="n",cex=0.6)
lines(1:3,(Yshort[19,]+Ylong[19,])/2,col=2,lwd=2)
lines(1:3,(Yshort[12,]+Ylong[12,])/2,col=3,lwd=2)
lines(1:3,(Yshort[7,]+Ylong[7,])/2,col=4,lwd=2)
lines(1:3,(Yshort[14,]+Ylong[14,])/2,col=5,lwd=2)

ymax=max((Yrel+Yneu+Yunr)/3)
ymin=min((Yrel+Yneu+Yunr)/3)
plot(1:2,((Yrel+Yneu+Yunr)/3)[1,],type="l",
     xlim=c(1,2.125),
     ylim=c(ymin,ymax),xaxt="n",
     xlab="Factor DELAY",
     ylab="Mean Response",
     main="Mean response over DELAY")
axis(1,at=1:2,labels = c("short","long"))
for (i in 1:n){
  lines(1:2,((Yrel+Yneu+Yunr)/3)[i,],lty=i)
}
legend("topright",legend=1:n,lty=1:n,
       title="Subject",bty="n",cex=0.6)
lines(1:2,((Yrel+Yneu+Yunr)/3)[19,],col=2,lwd=2)
lines(1:2,((Yrel+Yneu+Yunr)/3)[12,],col=3,lwd=2)
lines(1:2,((Yrel+Yneu+Yunr)/3)[7,],col=4,lwd=2)
lines(1:2,((Yrel+Yneu+Yunr)/3)[14,],col=5,lwd=2)

ymax=max(Ymat)
ymin=min(Ymat)
plot(1:6,Ymat[1,],type="l",ylim=c(ymin,ymax),
     xlim=c(1,6.5),
     xlab="Factors",xaxt="n",
     ylab="Mean Response",
     main="Mean response for Semantic Priming Data")
axis(1,at=1:6,labels=c("short/rel","long/rel","short/neu","long/neu",
                       "short/unrel","long/unrel"))
for (i in 2:n){
  lines(1:6,Ymat[i,],lty=i)
}
legend("topright",legend=1:n,lty=1:n,
       title="Subject",bty="n",cex=0.75)
lines(1:6,Ymat[19,],col=2,lwd=2)
lines(1:6,Ymat[12,],col=3,lwd=2)
lines(1:6,Ymat[7,],col=4,lwd=2)
lines(1:6,Ymat[14,],col=4,lwd=2)

interaction.plot(factor(semantic$Condition,levels=c("related", "neutral","unrelated")),
                 factor(semantic$Delay,levels=c("short","long")),
                 semantic$Resp,
                 xlab = "Condition",
                 ylab= "Mean Response",
                 trace.label = "Delay")
par(mfrow=c(1,1))


#############################################################################
# ESTIMATES REPORTED IN HERITIER
#############################################################################
# REML, TABLE 4.5 IN HERITIER
betaREMLsemanticHeritier=c(633.436,-18.071,18.563,-51.22,-3.690,16.809)
SQRTthetaREMLsemanticHeritier=c(122.622,0.006,29.433,100.73)

# CBS-MM TABLE 4.5 IN HERITIER
betaCBSMMsemanticHeritier=
  c(586.420,-17.876,14.317,-56.994,12.706,8.844)
SQRTthetaCBSMMsemanticHeritier=
  c(77.991,NA,27.199,81.885)

# ML FROM THE FUNCTION LME
options(contrasts=c("contr.sum","cont.poly"))
semanticML.lme=lme(fixed=Resp~Delay*Condition,
                 data=Semantic,
                 random=list(Subject=pdBlocked(list(pdIdent(~1),
                                                    pdIdent(~Delay-1),
                                                    pdIdent(~Condition-1)))),
                 method="ML")

summary(semanticML.lme)

betaMLsemanticLME=fixed.effects(semanticML.lme)
SQRTthetaMLsemanticLME=c(119.6671,0.006590472,28.72451, 98.31039)

# REML FROM THE FUNCTION LME
options(contrasts=c("contr.sum","cont.poly"))
semanticREML.lme=lme(fixed=Resp~Delay*Condition,
                   data=Semantic,
                   random=list(Subject=pdBlocked(list(pdIdent(~1),
                                                      pdIdent(~Delay-1),
                                                      pdIdent(~Condition-1)))),
                   method="REML")

summary(semanticREML.lme)

betaREMLsemanticLME=fixed.effects(semanticREML.lme)
SQRTthetaREMLsemanticLME=c(122.6223,0.006854707,29.43386, 100.7382)




########################################################
# MODEL
########################################################
# DESIGN MATRIX FIXED EFFECTS SEE PAGE 90 IN HERITIER
X=cbind(rep(1,6),
        c(rep(1,3),rep(-1,3)),
        c(1,0,-1,1,0,-1),
        c(0,1,-1,0,1,-1),
        c(1,0,-1,-1,0,1),
        c(0,1,-1,0,-1,1))
X=P%*%X
q=6

# DESIGN MATRICES RANDOM FIXED EFFECTS
Z1=matrix(rep(1,6),nrow=6,ncol=1)
Z2=cbind(c(rep(1,3),rep(0,3)),c(rep(0,3),rep(1,3)))
Z3=rbind(diag(1,3),diag(1,3))
varerror=diag(rep(1,6))

Z1=P%*%Z1
Z2=P%*%Z2
Z3=P%*%Z3

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
#c0t=0.000001

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
  y=as.matrix(Ymat[i,],ncol=1,nrow=5)
  betavectermtranslated=betavectermtranslated+w[i]*t(X)%*%solve(V0)%*%y  
  betamattermtranslated=betamattermtranslated+w[i]*t(X)%*%solve(V0)%*%X 
}

betaoldtranslated=solve(betamattermtranslated)%*%betavectermtranslated
Voldtranslated=V0
thetaoldtranslated=numeric(l)

#################################################
# setting starting values for beta and theta
betaoldtranslated=solve(t(X)%*%X)%*%t(X)%*%mu0
Voldtranslated=V0

thetaoldtranslated=numeric(l)

tmp=L1*100+L2*10+L3
sum100=0
sum110=0
sum101=0
for (i in 1:6){
  for (j in 1:6){
    sum100=sum100+ifelse(tmp[i,j]==100,1,0)
    sum110=sum110+ifelse(tmp[i,j]==110,1,0)
    sum101=sum101+ifelse(tmp[i,j]==101,1,0)
  }
}
t1=sum(V0[tmp==100])/sum100
t3=sum(V0[tmp==101])/sum101-t1
t2=sum(V0[tmp==110])/sum110-t1
t4=sum(diag(V0))-t1-t2-t3
thetaoldtranslated=cbind(t1,t2,t3,t4)

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
    y=matrix(Ymat[i,],ncol=1,nrow=k)
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
  
  # AVOIDING NEGATIVE VARIANCES
  max0=function(x){max(x,0)}
  #thetanewtranslated=apply(thetanewtranslated,MARGIN = 1,FUN = max0)
  
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

betaMLEsemanticRik=betanewtranslated
thetaMLEsemanticRik=thetanewtranslated

betaCBSsemanticRik=betanewtranslated
thetaCBSsemanticRik=thetanewtranslated

betaCBSsemanticRESCALE=betanewtranslated
thetaCBSsemanticRESCALE=thetanewtranslated


##############################################################
# RESULTS MLE
data.frame(betaMLrik=betaMLsemanticRik,
           betaMLlme=betaMLsemanticLME)

data.frame(thetaMLrik=thetaMLsemanticRik,
           thetaMLlme=round(SQRTthetaMLsemanticLME^2,3))

data.frame(SQRTthetaMLrik=sqrt(thetaMLsemanticRik),
           SQRTthetaMLlme=SQRTthetaMLsemanticLME)

##############################################################
# RESULTS REML
data.frame(betaREMLlme=betaREMLsemanticLME,
           betaREMLheritier=betaREMLsemanticHeritier)

data.frame(thetaREMLlme=round(SQRTthetaREMLsemanticLME^2,6),
           thetaREMLheritier=SQRTthetaREMLsemanticHeritier^2)

data.frame(SQRTthetaREMLlme=round(SQRTthetaREMLsemanticLME,3),
           SQRTthetaREMLheritier=SQRTthetaREMLsemanticHeritier)

##############################################################
# RESULTS CBS/MM
data.frame(betaCBSrik=betaCBSsemanticRik,
           betaCBSMMsemanticRik=betaCBSMMsemanticRik,
           betaCBSMMHeritier=betaCBSMMsemanticHeritier)

data.frame(SQRTthetaCBSrik=sqrt(thetaCBSsemanticRik),
           SQRTthetaCBSHeritier=SQRTthetaCBSMMsemanticHeritier)


# RESULTS CBS/MM RESCALE
data.frame(betaCBSrescale=betaCBSsemanticRESCALE,
           betaCBSMMrescale=betaCBSMMsemanticRESCALE,
           betaCBSMMHeritier=betaCBSMMsemanticHeritier)

data.frame(SQRTthetaCBSrescale=sqrt(thetaCBSsemanticRESCALE),
           SQRTthetaCBSHeritier=SQRTthetaCBSMMsemanticHeritier)


##########################################################################
# STANDARD ERRORS
#########################################################################
# ML
M=10000
beta=betaMLsemanticRik
theta=thetaMLsemanticRik
#theta[2]=theta[3]=0
VthetaMLE=theta[1]*L1+theta[2]*L2+theta[3]*L3+theta[4]*L4

# Asymptotic variances betahat
# Expression from our theory
varbetatranslatedMLE=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%solve(VthetaMLE)%*%X)

beta=betaMLsemanticLME
theta=thetaMLsemanticLME
VthetaLME=theta[1]*L1+theta[2]*L2+theta[3]*L3+theta[4]*L4

# Asymptotic variances betahat
# Expression from our theory
varbetatranslatedLME=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%solve(VthetaLME)%*%X)

round(varbetatranslatedMLE,3)
round(varbetatranslatedLME,3)

SEbetaMLsemanticRik=sqrt(diag(varbetatranslatedMLE)/n)
SEbetaMLsemanticLME=sqrt(diag(varbetatranslatedLME)/n)
SEbetaREMLsemanticHeritier=c(28.465,8.974,13.732,13.732,12.691,12.691)

data.frame(SEbetaMLrik=SEbetaMLsemanticRik,
           SEbetaMLlme=SEbetaMLsemanticLME,
           SEbetaREMLheritier=SEbetaREMLsemanticHeritier)

##############################
# CBS-MM RESCALED
M=0
c0t=c1

beta=betaCBSMMsemanticRESCALE
theta=thetaCBSsemanticRESCALE
VthetaCBS=theta[1]*L1+theta[2]*L2+theta[3]*L3+theta[4]*L4
L=as.matrix(cbind(vec(L1),vec(L2),vec(L3),vec(L4)))

# Asymptotic variances betahat
# Expression from our theory
varbetatranslatedCBS=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%solve(VthetaCBS)%*%X)
# Expression from Copt & Victoria-Feser (2006)
varbetaCopttranslatedCBS=constbetahattranslated(k,M,c0t)*
  solve(t(X)%*%X) %*% t(X) %*% VthetaCBS %*% X %*% solve(t(X)%*%X)

round(varbetatranslatedCBS,3)
round(varbetaCopttranslatedCBS,3)

SEbetaCBSMMsemanticRESCALE=sqrt(diag(varbetatranslatedCBS)/n)
SEbetaCBSMMsemanticHeritier=c(18.817,6.082,11.691,11.691,10.582,10.582)

data.frame(SEbetaCBSMMrescale=SEbetaCBSMMsemanticRESCALE,
           SEbetaCBSMMheritier=SEbetaCBSMMsemanticHeritier)


############################################################
# WALD TEST
###########################################################
waldMLE=betaMLsemanticRik/SEbetaMLsemanticRik
pMLE=2*(1-pnorm(abs(waldMLE)))

data.frame(MLE=betaMLEsemanticRik,
           SE.MLE=round(SEbetaMLEsemanticRik,3),
           TMLE=waldMLE,
           PvalMLE=round(pMLE,4))

waldCBSMM=betaCBSMMsemanticRESCALE/SEbetaCBSMMsemanticRESCALE
pCBSMM=2*(1-pnorm(abs(waldCBSMM)))

data.frame(CBSMM=betaCBSMMsemanticRESCALE,
           SE.CBSMM=round(SEbetaCBSMMsemanticRESCALE,3),
           TCBS=waldCBSMM,
           PvalCBSMM=round(pCBSMM,4))


############################################################################################
# RESISDUALS
############################################################################################
YMLEsemanticMat=matrix(0,nrow=n,ncol=k)
YCBSMMsemanticMat=matrix(0,nrow=n,ncol=k)

for (i in 1:n){
  YMLEsemanticMat[i,]=X%*%betaMLEsemanticRik
  YCBSMMsemanticMat[i,]=X%*%betaCBSMMsemanticRESCALE
}

resMLEsemantic=Ymat-YMLEsemanticMat
resCBSMMsemantic=Ymat-YCBSMMsemanticMat


############################################################################################
# MAHALANOBIS DISTANCES OF RESIDUALS
a=0.025

MDerrorMLEsemantic=numeric(n)
MDerrorCBSMMsemantic=numeric(n)
VMLEsemantic=thetaMLEsemanticRik[1]*L1+thetaMLEsemanticRik[2]*L2+thetaMLEsemanticRik[3]*L3+thetaMLEsemanticRik[4]*L4
VCBSsemantic=thetaCBSsemanticRESCALE[1]*L1+thetaCBSsemanticRESCALE[2]*L2+thetaCBSsemanticRESCALE[3]*L3+thetaCBSsemanticRESCALE[4]*L4


for (i in 1:n){
  MDerrorMLEsemantic[i]=mahalanobis(resMLEsemantic[i,],center=FALSE,cov=VMLEsemantic) # Note MD=d^2!
  MDerrorCBSMMsemantic[i]=mahalanobis(resCBSMMsemantic[i,],center=FALSE,cov=VCBSsemantic) # Note MD=d^2!
}

plot(sqrt(MDerrorMLEsemantic),sqrt(MDerrorCBSMMsemantic),
     xlab="MD's of MLE residuals",
     ylab="MD's of CBS-MM residuals",
     main="Mahalanobis distances for semantic data")

IDsemantic=(1:n)[MDerrorCBSMMsemantic>qchisq(0.975,df=k)]

points(sqrt(MDerrorMLEsemantic)[IDsemantic],sqrt(MDerrorCBSMMsemantic)[IDsemantic],pch=16,col=2)
text(sqrt(MDerrorMLEsemantic)[IDsemantic],sqrt(MDerrorCBSMMsemantic)[IDsemantic],label=IDsemantic,pos=2,col=2)
abline(h=sqrt(qchisq(0.975,df=k)),col=2)
abline(v=sqrt(qchisq(0.975,df=k)),col=2)

#######################################################################################
# DATA REVISITED WITH MY OUTLIERS
ymax=max(Ymat)
ymin=min(Ymat)
plot(1:6,Ymat[1,],type="l",ylim=c(ymin,ymax),
     xlim=c(1,6.5),
     xlab="Factors",xaxt="n",
     ylab="Mean Response",
     main="Semantic Priming Data")
axis(1,at=1:6,labels=c("short/rel","long/rel","short/neu","long/neu",
                       "short/unrel","long/unrel"))
for (i in 2:n){
  lines(1:6,Ymat[i,],lty=i)
}
lines(1:6,Ymat[19,],col=2,lwd=2)
lines(1:6,Ymat[12,],col=3,lwd=2)
lines(1:6,Ymat[7,],col=4,lwd=2)
lines(1:6,Ymat[14,],col=5,lwd=2)

legend("topright",legend=1:n,
       lty=c(1:6,1,8:11,1,13,1,15:18,1,20:n),
       lwd=c(1,1,1,1,1,1,2,1,1,1,1,2,1,2,1,1,1,1,2,1,1),
       col=c(1,1,1,1,1,1,4,1,1,1,1,3,1,5,1,1,1,1,2,1,1),
       title="Subject",bty="n",cex=1)


# DATA REVISITED WITH HERITIER OUTLIERS
ymax=max(Ymat)
ymin=min(Ymat)
plot(1:6,Ymat[1,],type="l",ylim=c(ymin,ymax),
     xlim=c(1,6.5),
     xlab="Factors",xaxt="n",
     ylab="Mean Response",
     main="Semantic Priming Data")
axis(1,at=1:6,labels=c("short/rel","long/rel","short/neu","long/neu",
                       "short/unrel","long/unrel"))
for (i in 2:n){
  lines(1:6,Ymat[i,],lty=i)
}
lines(1:6,Ymat[3,],col=2,lwd=2)
lines(1:6,Ymat[16,],col=3,lwd=2)
lines(1:6,Ymat[8,],col=4,lwd=2)

legend("topright",legend=1:n,
       lty=c(1:2,1,4:7,1,9:15,1,17:n),
       lwd=c(1,1,2,1,1,1,1,2,1,1,1,1,1,1,1,2,1,1,1,1,1),
       col=c(1,1,2,1,1,1,1,4,1,1,1,1,1,1,1,3,1,1,1,1,1),
       title="Subject",bty="n",cex=1)


#########################################################
# INTERACTION PLOTS
# OVERALL EFFECT
boxplot(semantic[,4])

# MAIN EFFECTS
par(mfrow=c(1,2))
boxplot(semantic[semantic[,2]=="short",4],
        semantic[semantic[,2]=="long",4],
        names=c("short","long"),
        xlab="Delay",
        main="Response for Delay")

boxplot(semantic[semantic[,3]=="related",4],
        semantic[semantic[,3]=="neutral",4],
        semantic[semantic[,3]=="unrelated",4],
        names=c("related","neutral","unrelated"),
        xlab="Condition",
        main="Response for Condition")
par(mfrow=c(1,1))

names(semantic)
summary.aov(lm(Resp~Delay+Condition,data=semantic))

semanticSHORT=semantic[semantic[,2]=="short",]
semanticLONG=semantic[semantic[,2]=="long",]
semanticREL=semantic[semantic[,3]=="related",]
semanticNEU=semantic[semantic[,3]=="neutral",]
semanticUNR=semantic[semantic[,3]=="unrelated",]

Yshort=numeric(n)
Ylong=numeric(n)
Yrel=numeric(n)
Yneu=numeric(n)
Yunr=numeric(n)
for (i in 1:n){
  Yshort[i]=mean(semanticSHORT[semanticSHORT[,1]==i,4])
  Ylong[i]=mean(semanticLONG[semanticLONG[,1]==i,4])
  Yrel[i]=mean(semanticREL[semanticREL[,1]==i,4])
  Yneu[i]=mean(semanticNEU[semanticNEU[,1]==i,4])
  Yunr[i]=mean(semanticUNR[semanticUNR[,1]==i,4])
}

par(mfrow=c(2,2))
boxplot(apply(Ymat,MARGIN=1,FUN=mean),
        main="Mean",
        xlab="",ylab="Mean Response")

boxplot(Yshort,Ylong,names=c("Short","Long"),
        main="Mean over Delay",xlab="Delay",ylab="Mean Response Condition")

boxplot(Yrel,Yneu,Yunr,
        main="Mean over Condition",xlab="Condition",ylab="Mean Response Delay")

interaction.plot(factor(semantic$Delay,levels=c("short","long")),
                 factor(semantic$Condition,levels=c("related", "neutral","unrelated")),
                 semantic$Resp,
                 xlab = "Delay",
                 ylab= "Mean Response",
                 trace.label = "Condition")
par(mfrow=c(1,1))



ymax=max(Ymat)
ymin=min(Ymat)
Yunrelated=cbind(semanticUNR[semanticUNR[,2]=="short",4],
                 semanticUNR[semanticUNR[,2]=="long",4])
Yneutral=cbind(semanticNEU[semanticNEU[,2]=="short",4],
               semanticNEU[semanticNEU[,2]=="long",4])
Yrelated=cbind(semanticREL[semanticREL[,2]=="short",4],
               semanticREL[semanticREL[,2]=="long",4])

par(mfrow=c(2,2))
plot(1:2,Yunrelated[1,],type="l",ylim=c(ymin,ymax),
     xlim=c(0.5,2.5),xaxt="n",
     main="Unrelated",xlab="Delay",ylab="Mean repsonse for Unrelated")
axis(1,at=1:2,labels=c("short","long"))     
for (i in 2:n){
  lines(1:2,Yunrelated[i,],lty=i)
}
legend("topright",legend=1:n,lty=1:n,
       title="Subject",bty="n")
lines(1:2,Yunrelated[19,],lwd=2,col=2)
lines(1:2,Yunrelated[12,],lwd=2,col=3)

plot(1:2,Yneutral[1,],type="l",ylim=c(ymin,ymax),
     xlim=c(0.5,2.5),xaxt="n",
     main="Neutral",xlab="Delay",ylab="Mean repsonse for Neutral")
axis(1,at=1:2,labels=c("short","long"))     
for (i in 2:n){
  lines(1:2,Yneutral[i,],lty=i)
}
legend("topright",legend=1:n,lty=1:n,
       title="Subject",bty="n")
lines(1:2,Yneutral[19,],lwd=2,col=2)
lines(1:2,Yneutral[12,],lwd=2,col=3)

plot(1:2,Yrelated[1,],type="l",ylim=c(ymin,ymax),
     xlim=c(0.5,2.5),xaxt="n",
     main="Related",xlab="Delay",ylab="Mean repsonse for Related")
axis(1,at=1:2,labels=c("short","long"))     
for (i in 1:n){
  lines(1:2,Yrelated[i,],lty=i)
}

legend("topright",legend=1:n,lty=1:n,
       title="Subject",bty="n")
lines(1:2,Yrelated[19,],lwd=2,col=2)
lines(1:2,Yrelated[12,],lwd=2,col=3)

interaction.plot(factor(semantic$Delay,levels=c("short","long")),
                 factor(semantic$Condition,levels=c("related", "neutral","unrelated")),
                 semantic$Resp,
                 xlab = "Delay",
                 ylab= "Mean Response",
                 trace.label = "Condition")
par(mfrow=c(1,1))

par(mfrow=c(1,2))
interaction.plot(factor(semantic$Delay,levels=c("short","long")),
                 factor(semantic$Condition,levels=c("related", "neutral","unrelated")),
                 semantic$Resp,
                 xlab = "Delay",
                 ylab= "Mean Response",
                 trace.label = "Condition")

interaction.plot(factor(semantic$Condition,levels=c("related", "neutral","unrelated")),
  factor(semantic$Delay,levels=c("short","long")),
                 semantic$Resp,
                 xlab = "Condition",
                 ylab= "Mean Response",
                 trace.label = "Delay")
par(mfrow=c(1,1))