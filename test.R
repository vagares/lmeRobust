library(robustvarComp)

data(Orthodont)

z1 <- rep(1, 4)
z2 <- c(8,10,12,14)
K <- list()
K[[1]] <- tcrossprod(z1,z1) ## Int
K[[2]] <- tcrossprod(z1,z2) + tcrossprod(z2,z1) ## Int:age
K[[3]] <- tcrossprod(z2,z2) ## age
names(K) <- c("Int", "Int:age", "age")
p <- 4
n <- 27
groups <- cbind(rep(1:p, each=n), rep((1:n), p))


# }
# NOT RUN {
## Composite S
OrthodontCompositeS <- varComprob(distance ~ age*Sex, groups = groups,
                                  data = Orthodont, varcov = K,
                                  control=varComprob.control(method="compositeS", lower=c(0,-Inf,0)))


est1 = Roblme(Ymat,Xlist,Zlist,E=diag(rep(1,4)),rho ="t-biweight",r =0.5,arp=0.01,rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL)


est2 = Roblme(Ymat,Xlist,Zlist,E=NULL,rho ="t-biweight",r =0.5,arp=0.01,rhoMM=NULL,eps=1e-5,maxiter=100,eff=0.95,V0=NULL)

est1$fixedeffectsS
est2$fixedeffectsS

