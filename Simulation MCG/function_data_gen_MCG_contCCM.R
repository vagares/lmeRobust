# This script contains the code for the function data_gen_MCG_contCCM.

# This function generates a single dataset according to the model in 
# Mason, Cantoni & Ghisletta (2021) with contamination generated 
# according to the CCM (central contamination model) 
# in the random effects, in the measurement error, 
# and in the design matrix of the fixed effects.


library(mvtnorm)    # needed for rmvnorm to generate 
# a vector of random effects 

data_gen_MCG_contCCM = function(n=200,k=4,pe=0,pb=0,px=0,mec=0,mbc2=0,alphac=1,no=0){
  # n number of individuals
  # k number of observations per individual
  # pe probability of having outlier in component epsilon_ij
  # pb probability of having outlier in vector of random effects b
  # px probability of having outlier in component x_ij in X
  # mec shift in the mean of component epsilon_ij
  # mbc2 shift in the mean of random effect b2
  # alphac multiplication factor in component x_ij in X
  # no indicator about a random number (no=0) of outliers or a fixed number (no=1)
  
  X = list()          # list of design matrices for each individual
  x1 = rep(1,k)       # first column of design matrix X
  x2 = seq(0,k-1,1)   # second column of design matrix X
  
  Zmat = cbind(x1,x2) # design matrix Z of vector of random effects
  Z=list(Zmat)        # Creating list for Roblme
  
  # preparing X-list for Roblme
  for (i in (1:n)){
    X[[i]] <- cbind(x1,x2)
  }
  
  # Setting up linear representation of V
  Sigma1=matrix(c(1,0,0,0),ncol=2)  # corresponds to var(b1)
  Sigma2=matrix(c(0,1,1,0),ncol=2)  # corresponds to covar(b1,b2)
  Sigma3=matrix(c(0,0,0,1),ncol=2)  # corresponds to var(b2)
  
  L1=Zmat%*%Sigma1%*%t(Zmat)          # term in V from var(b1)
  L2=Zmat%*%Sigma2%*%t(Zmat)          # term in V from covar(b1,b2)
  L3=Zmat%*%Sigma3%*%t(Zmat)          # term in V from var(b21)
  L4=diag(1,k)                        # term in V from var(epsilon)

  # preparing L-list for Roblme
  Llist=list()          
  Llist[[1]]=L1
  Llist[[2]]=L2
  Llist[[3]]=L3
  Llist[[4]]=L4
  
  Y = list()  # preparing Y-list containing individual observations
  
  sb0 = matrix(c(790,-8.5,-8.5,40),2,2)       # var of random effect b
  Se0 = 400*diag(1,k)                         # var of eps
  mbc = c(0,mbc2)                             # contamination in mean of b
  sbc =  matrix(c(7.9,-0.085,-0.085,0.4),2,2) # contamination in var of b
  sec = 0.5                                   # contamination sd of eps_ij
  
  # counters for number of contaminated observations/individuals
  nobi=0
  noe=0
  nox=0
  noei=0
  noxi=0

  beta0 = c(250,10)  # fixed effects
  
  # Contamination according to Central Contamination Model (CCM)
  for (i in (1:n)){
    noeitemp=0
    noxitemp=0
    
    if (no==0){
    # generating a random effect
      if (runif(1) <= (1-pb)){b= rmvnorm(1, mean=rep(0,2), sigma=sb0)}else{
        b=rmvnorm(1, mean=mbc, sigma=sbc);nobi=nobi+1}
    
    # generating measurement error
      if (runif(1) <= (1-pe)){eps= rmvnorm(1, mean=rep(0,k), sigma=Se0)}else{
        eps=rmvnorm(1, mean=c(mec,rep(0,(k-1))), sigma=sec^2*diag(1,k));nobi=nobi+1}

    # construct Y[[i]] according to MCG-model
      Y[[i]] =  X[[i]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
  
    # constructing x2 under contamination
      if(runif(1) <= (1-px)){ X[[i]][,2]= X[[i]][,2]}else{X[[i]][,2]= alphac*X[[i]][,2];nox=nox+1;noxitemp=1}

    # keeping track of number of outliers
      noei = noei + noeitemp
      noxi = noxi + noxitemp
    }else{
      b= rmvnorm(1, mean=rep(0,2), sigma=sb0)
      eps= rmvnorm(1, mean=rep(0,k), sigma=Se0)
      Y[[i]] =  X[[i]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
    }
  } # END of i-loop
  
  if (no!=0){   
    if (pe>0){
        no_eps=as.integer(floor(n*pe))
          for (j in 1:no_eps){
            b= rmvnorm(1, mean=rep(0,2), sigma=sb0)
            eps=rmvnorm(1, mean=c(mec,rep(0,(k-1))), sigma=sec^2*diag(1,k))
            Y[[j]] =  X[[j]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
          }
      }
  
      if (pb>0){
        no_b=as.integer(floor(n*pb))
          for (j in 1:no_b){
            b=rmvnorm(1, mean=mbc, sigma=sbc)
            eps= rmvnorm(1, mean=rep(0,k), sigma=Se0)
            Y[[j]] =  X[[j]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
          }
      }
  
      if (px>0){
        no_x=as.integer(floor(n*px))
        for (j in 1:no_x){
          b= rmvnorm(1, mean=rep(0,2), sigma=sb0)
          eps= rmvnorm(1, mean=rep(0,k), sigma=Se0)
          X[[j]][,2]= alphac*X[[j]][,2]
          Y[[j]] =  X[[j]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
        }
      }
    } # END of if (no!=0)

  # Setting up Y-matrix for Roblme
  Ymat=matrix(0,n,k)
  for (i in (1:n)){Ymat[i,] = Y[[i]]}
  
  
  return(list(Y=Ymat,X=X,Z=Z,L=Llist,
              pe=pe,pb=pb,px=px,
              mec=mec,mbc2=mbc2,alphac=alphac,
              nobi=nobi,noe=noe,nox=nox,noei=noei,noxi=noxi,
              no=no))
}

