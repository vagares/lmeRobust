
#TO CHANGE take from lme

lmeRob=function(fixed, data, random,...){
  mc=match.call()
  return(mc)
}


lme.to.Robustlme=function(mc){
  #@ 17/12/06 (dlc)
  #@ - include new warnings and stop
  #@ - no problem with "data" anymore
  #@ - change slightly the output
  #@   - small matrices (for ctbs) now called y, x, and z (in yxz)
  #@   - new big matrices for the predictions and residuals (in YXZ)
  #@   - dim contains some dim of the fix, random, ...
  #
  # USEFUL FUNCTIONS
  #
  # (Bates/lme4)
  findbars <- function(term)
  {
    if (is.name(term) || !is.language(term)) return(NULL)
    if (term[[1]] == as.name("(")) return(findbars(term[[2]]))
    if (!is.call(term)) stop("term must be of class call")
    if (term[[1]] == as.name('|')) return(term)
    if (length(term) == 2) return(findbars(term[[2]]))
    c(findbars(term[[2]]), findbars(term[[3]]))
  }
  # (Bates/lme4)
  slashTerms <- function(x) {
    if (!("/" %in% all.names(x))) return(x)
    if (x[[1]] != as.name("/"))
      stop("unparseable formula for grouping factor")
    list(slashTerms(x[[2]]), slashTerms(x[[3]]))
  }
  # (Bates/lme4)
  makeInteraction <- function(x) {
    if (length(x) < 2) return(x)
    trm1 <- makeInteraction(x[[1]])
    trm11 <- if(is.list(trm1)) trm1[[1]] else trm1
    list(substitute(foo:bar, list(foo=x[[2]], bar = trm11)), trm1)
  }
  # (Bates/lme4)
  expandSlash <- function(bb) {
    if (!is.list(bb)) return(expandSlash(list(bb)))
    unlist(lapply(bb, function(x) {
      if (length(x) > 2 && is.list(trms <- slashTerms(x[[3]])))
        return(lapply(unlist(makeInteraction(trms)),
                      function(trm) substitute(foo|bar,
                                               list(foo = x[[2]],
                                                    bar = trm))))
      x
    }))
  }
  #
  # CHANGE THE MATCH.CALL (if different order given)
  #
  m <- match(c("fixed","data","random"),
             names(mc), 0)
  mc <- mc[c(1, m)]
  #
  # GET THE DATA / FIND THE GROUPING VARIABLE / ORDER THE DATASET / SIZE
  #
  data=get(as.character(mc[[3]]))
  re=mc[[4]]
  if(re[[1]]=="list"){stop("Stop:\n-> The random argument of lmeRob can't be a list\n")}
  re0= expandSlash(findbars(re))
  names(re0)=unlist(lapply(re0, function(x) deparse(x[[3]])))
  lre0=length(re0)
  gv=names(re0)[lre0]
  # the order of the dataset is paramount (as Z is the same for every unit)
  resp=all.vars(mc[[2]][[2]])
  allvar=dimnames(data)[[2]]
  allvar=allvar[allvar!=resp&allvar!=gv]
  #if (length(allvar)>=1){
  for(i in 1:length(allvar)){data=data[order(data[,allvar[i]]),]}
  #}
  data0=data0b=data[order(data[,gv]),]
  datan=split(data0,as.numeric(data0[,gv]))
  # n = number of indep unit / m = number of obs per unit / N = overall number of obs
  tgv=table(data0[,gv]);n=length(tgv)
  N=dim(data0)[1]
  m=N/n;if(sum(tgv==m)!=n){stop("Stop:\n-> Need balanced data\n")}
  #
  # FIXED EFFECTS
  #
  mf=mc[1]
  mf$formula=mc[[2]];mf$data=data0
  mf[[1]]=as.name("model.frame")
  mf=eval(mf,parent.frame(2))
  mt=attr(mf,"terms")
  # Build Y (matrix) and X (list of matrices)
  Y=matrix(model.response(mf, "numeric"),ncol=1)
  y=matrix(Y,ncol=m,nrow=n,byrow=T)
  X=model.matrix(mt,mf,contrasts)
  x=as.list(rep(NA,n))
  for(i in 1:n){x[[i]]=X[data0[,gv]==names(tgv[i]),]}
  # Change the contrast of X for the robust multivariate test:
  # (need only sum contrats for the X matrix)
  if(length(allvar)>1){factallvar = sapply(data0b[,allvar],class)=="factor"
  }else{factallvar = class(data0b[,allvar])=="factor"}
  factallvar=seq(1,length(factallvar))[factallvar]
  if (length(factallvar)>1){for(fav in 1:length(factallvar)){contrasts(data0b[,allvar[factallvar[fav]]])="contr.sum"}}
  #for(fav in 1:length(factallvar)){contrasts(data0b[,allvar[factallvar[fav]]])="contr.sum"}
  mf2=mc[1]
  mf2$formula=mc[[2]];mf2$data=data0b
  mf2[[1]]=as.name("model.frame")
  mf2=eval(mf2,parent.frame(2))
  mt2=attr(mf2,"terms")
  X2=model.matrix(mt2,mf2,contrasts)
  x2=as.list(rep(NA,n))
  for(i in 1:n){x2[[i]]=X2[data0[,gv]==names(tgv[i]),]}
  #
  # RANDOM EFFECTS
  #
  # case of nested structure...
  if(lre0>1){for(i in 1:lre0){names(re0)[i]=all.vars(re0[[i]])[1]}}
  re0=re0[lre0:1]
  # God save the Queen... for every level
  Zf=NULL
  for(l in 1:lre0){
    rel=as.list(re0[[l]])
    gel=as.character(pairlist(rel[[3]]))
    felv=all.vars(rel[[2]]);n.felv=length(felv)
    felft=terms(as.formula(paste("~",as.character(re0[[1]])[2])))
    # test if a fe is a between subject factor...
    if(n.felv>0){
      for(f in 1:n.felv){
        if(sum(unlist(lapply(datan,function(x,vect){sum(x[,felv[f]]==data0[1:m,felv[f]])==m})))!=n){
          cat(paste("Warning:\n->",felv[f],"seems to be a between factor: interaction with",gel,"not allowed\n"))
        }else{felv[f]=NA}
      }
      # delete interaction with between factors
      felv=felv[!is.na(felv)];n.felv=length(felv)
      fel=attr(felft,"factors")
      if(n.felv>0){for(f in 1:n.felv){fel=fel[,fel[felv[f],]==0,drop=F]}}
      fel=paste(dimnames(fel)[[2]],gel,sep=":")
    }else{fel=NULL}
    # add the wished formula in a list
    Zf=list(Zf,lapply(as.list(paste("~",c(gel,fel),"-1")),
                      as.formula))
  }
  # Build the Z and z matrices
  r=length(unlist(Zf))+1
  Z=as.list(rep(NA,r))
  Z[[r]]=diag(N);dimnames(Z[[r]])[[2]]=as.character(1:N)
  Z[-r]=lapply(unlist(Zf),
               function(x)model.matrix(eval(x),data0)
  )
  if(sum(sapply(Z[-r],function(x)dim(x)[2]==N)))stop("Stop:\n-> Too many random effects\n")
  z=lapply(Z,function(x)x[1:m,1:(dim(x)[2]/n),drop=F])
  dims=list(N=N,n=n,m=m,level=lre0,r=r,qv=sapply(Z,function(x)dim(x)[2]),gv=gv,
            idqv=lapply(Z,function(x)dimnames(x)[[2]]),
            idZ=c(sapply(unlist(Zf),function(x)as.character(pairlist(x[[2]][[2]]))),"Residuals"),
            idi=as.vector(data0[,gv])
  )
  names(dims$idqv)=dims$idZ
  # inf for scatterplot
  data1=data0[1:n,!(dimnames(data0)[[2]]==resp|dimnames(data0)[[2]]==gv)]
  if(is.data.frame(data1)){
    for(i in 1:(dim(data1)[2])){data1[,i]=as.character(data1[,i])}
    for(i in 1:n){data1[i,1]=paste(data1[i,],collapse="\n")}
    names.scat=as.vector(data1[,1])
  }else{names.scat=NA}
  # output
  invisible(list(yxz=list(y=y,x=x,z=z,x2=x2),YXZ=list(Y=Y,X=X,Z=Z,X2=X2),
                 mc=mc,mt=mt,mf=mf,names.scat=names.scat,dim=dims))
}