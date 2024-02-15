# This script contains the code for the function data_gen_MCG_CCMind.

# CCMind=FALSE 
# Then this function generates a single dataset according to the model 
# in Mason, Cantoni & Ghisletta (2021) with contamination 
# generated according to the ICM (independent contamination model 
# or cellwise contamination) in the measurement error and 
# according to CCM (central contamination model) in the random effects.
# In addition, the function also generates contamination 
# in the design matrix of the fixed effects according to ICM.

# CCMind=TRUE 
# Then this function generates a single dataset according to the model 
# in Mason, Cantoni & Ghisletta (2021) with contamination generated 
# according to the CCM (central contamination model) 
# in the random effects, in the measurement error, 
# and in the design matrix of the fixed effects.

library(mvtnorm)    
# needed for rmvnorm to generate 
# a vector of random effects 

data_gen_MCG_CCMind = function(n=200,k=4,pe=0,pb=0,px=0,mec=0,mbc2=0,
                               alphac=1,randcont=0,CCMind=FALSE){
  # n number of individuals
  # k number of observations per individual
  # pe probability of having outlier in component epsilon_ij
  # pb probability of having outlier in vector of random effects b
  # px probability of having outlier in component x_ij in X
  # mec shift in the mean of component epsilon_ij
  # mbc2 shift in the mean of random effect b2
  # alphac multiplication factor in component x_ij in X
  # randcont indicator about a random number (randcont=1) of outliers or a fixed number (randcont=0)
  # CCMind index for casewise contamination (TRUE) or independent contamination (FALSE)
  
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
  sde = 20                                    # standard deviation of eps_ij  
  Se0 = sde^2*diag(1,k)                       # covariance matrix of eps
  mbc = c(0,mbc2)                             # contamination in mean of b
  sbc =  matrix(c(7.9,-0.085,-0.085,0.4),2,2) # contamination in var of b
  sec = 0.5                                   # contamination sd of eps_ij
  
  # counters for number of contaminated observations/individuals
  nobi=0  # number of contaminated random effect
  noe=0   # number of cell wise contaminated error 
  nox=0   # number of cell wise contaminated X-column
  noei=0  # number of case wise contaminated error
  noxi=0  # number of case wise contaminated X-column

  beta0 = c(250,10)  # fixed effects
  
  if (CCMind==TRUE){
  # Contamination according to Central Contamination Model (CCM)
  for (i in (1:n)){
    noeitemp=0
    noxitemp=0
    
    if (randcont!=0){
    # generating a random effect
      if (runif(1) <= (1-pb)){b= rmvnorm(1, mean=rep(0,2), sigma=sb0)}else{
        b=rmvnorm(1, mean=mbc, sigma=sbc);nobi=nobi+1}
    
    # generating measurement error
    # randomly selected being an outlier
      if (runif(1) <= (1-pe)){eps= rmvnorm(1, mean=rep(0,k), sigma=Se0)}else{
        eps=rmvnorm(1, mean=c(mec,rep(0,(k-1))), sigma=sec^2*diag(1,k));
        noeitemp=1}

    # construct Y[[i]] according to MCG-model
      Y[[i]] =  X[[i]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
  
    # constructing x2 under contamination AFTER Y to create leverage points
    # randomly selected being an outlier
      if(runif(1) <= (1-px)){ X[[i]][,2]= X[[i]][,2]}else{
        X[[i]][,2]= alphac*X[[i]][,2];
        noxitemp=1}

    # keeping track of number of case wise outliers
      noei = noei + noeitemp
      noxi = noxi + noxitemp
    # the numbers noe/nox of cell wise outliers being kept 0
    }else{
      # generate uncontaminated observation under randcont=0
      b= rmvnorm(1, mean=rep(0,2), sigma=sb0)
      eps= rmvnorm(1, mean=rep(0,k), sigma=Se0)
      Y[[i]] =  X[[i]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
    }
  } # END of i-loop
  
  if (randcont==0){  
    # we fix the number of case wise outliers 
    if (pe>0){
        no_eps=as.integer(floor(n*pe))
          for (j in 1:no_eps){
            b= rmvnorm(1, mean=rep(0,2), sigma=sb0)
            eps=rmvnorm(1, mean=c(mec,rep(0,(k-1))), sigma=sec^2*diag(1,k))
            Y[[j]] =  X[[j]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
          }
      noei=no_eps
      }
  
      if (pb>0){
        no_b=as.integer(floor(n*pb))
          for (j in 1:no_b){
            b=rmvnorm(1, mean=mbc, sigma=sbc)
            eps= rmvnorm(1, mean=rep(0,k), sigma=Se0)
            Y[[j]] =  X[[j]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
          }
      nobi=no_b
      }
  
      if (px>0){
        no_x=as.integer(floor(n*px))
        for (j in 1:no_x){
          b= rmvnorm(1, mean=rep(0,2), sigma=sb0)
          eps= rmvnorm(1, mean=rep(0,k), sigma=Se0)
          Y[[j]] =  X[[j]]%*%beta0 + Z[[1]]%*%as.vector(b) + as.vector(eps)
          
          # Contaminate X AFTER Y to create leverage points
          X[[j]][,2]= alphac*X[[j]][,2]
        }
      noxi=no_x
      }
    } # END of if (randcont==0)
}# END of if CCMind=TRUE
  
  if (CCMind==FALSE){
    # Generating cell wise contamination (ICM)
    for (i in (1:n)){
      noeitemp=0
      noxitemp=0
      
      # generating a random effect
      if (runif(1) <= (1-pb)){b= rmvnorm(1, mean=rep(0,2), sigma=sb0)}else{
        b=rmvnorm(1, mean=mbc, sigma=sbc);nobi=nobi+1}
      
      # generating measurement error
      eps = numeric(k)
      
      # generating individual components of measurement error
      for (kk in (1:k)){
        if (runif(1) <= (1-pe)){eps[kk] = rnorm(1, mean=0, sd=sde)}
        else{eps[kk] = rnorm(1, mean=mec, sd=sec);noe=noe+1;noeitemp=1}
      }
      
      beta0 = c(250,10)  # fixed effects
      
      # construct Y[[i]] according to MCG-model
      Y[[i]] =  X[[i]]%*%beta0 + Z[[1]]%*%as.vector(b) + eps

      # constructing components of x2 under contamination
      for (r in (1:k)){
        if(runif(1) <= (1-px)){ X[[i]][r,2]= X[[i]][r,2]}
        else{X[[i]][r,2]= alphac*X[[i]][r,2];nox=nox+1;noxitemp=1}
      }
      # The reason to contaminate X after the construction of Y is to investigate the effect
      # of a "leverage point" on the estimators


      # keeping track of number of outliers
      noei = noei + noeitemp
      noxi = noxi + noxitemp
    }
  } # END of if CCMind=FALSE
  
  
  # Setting up Y-matrix for Roblme
  Ymat=matrix(0,n,k)
  for (i in (1:n)){Ymat[i,] = Y[[i]]}
  
  
  return(list(Y=Ymat,X=X,Z=Z,L=Llist,
              pe=pe,pb=pb,px=px,
              mec=mec,mbc2=mbc2,alphac=alphac,
              noei=noei,noe=noe,nobi=nobi,noxi=noxi,nox=nox))
}

