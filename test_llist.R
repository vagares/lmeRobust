load("~/lmeRobust/ATT52969.RData")




est1 = Roblme(Ymat,Xlist,Zlist,E=diag(rep(1,4)),L=NULL,rho ="t-biweight",r =0.5,arp=0.01,rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL)


est2 = Roblme(Ymat,Xlist,Zlist,E=NULL,L=NULL,rho ="t-biweight",r =0.5,arp=0.01,rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL)

est1$fixedeffectsS
est2$fixedeffectsS

L=list()
lZ = length(Z)
k = ncol(Ymat)



s1 = matrix(c(1,0,0,0),2,2)

s2 = matrix(c(0,1,1,0),2,2)

s3 = matrix(c(0,0,0,1),2,2)

L[[1]] = Z[[1]] %*% s1%*% t(Z[[1]])
  
L[[2]] =  Z[[1]] %*% s2%*% t(Z[[1]])
  
L[[3]] =  Z[[1]] %*% s3%*% t(Z[[1]])
  
L[[4]] = diag(rep(1,k))
