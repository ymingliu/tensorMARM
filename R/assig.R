assig <- function(n_args){
  cargs <- vector("list", length(n_args))
  for(i in 1:length(n_args)) cargs[[i]] <- 1:n_args[i]
  t(expand.grid(cargs))
}

##------------------------------------------------------------------------##
##--------------------produce the B-spline functions----------------------##
bsbasefun <- function(X,K,degr){
  n = dim(X)[1]
  p = dim(X)[2]
  nk = K - degr
  u.k = seq(0, 1, length=nk+2)[-c(1,nk+2)]
  BS = NULL
  for(j in 1:p){
    Knots = as.numeric(quantile(X[,j], u.k))  
    BS0 = bs(X[,j], knots=Knots, intercept=TRUE, degree=degr)
    BS = cbind(BS,BS0[,-1])
  }
  BS = scale(BS,center = T, scale = F)
  id = seq(1,p*K,K)
  Z = NULL
  for(j in 1:K){
    Z = cbind(Z,BS[,id+j-1])
  }
  return(Z)
}

##------------------------------------------------------------------------##
##-----------------generate scenario I data in IMMAM----------------------##
# e.g. mydata=immam3.sim.fbs(200,5,100,10,rep(1:4,each=25))
immam3.sim.fbs <- function(n,q,p,s,group=NULL,r10=2,r20=2,r30=2,isfixedR=0,D3=NULL,
                           K=6,degr=3,sigma2=NULL,seed_id=NULL,
                           r1_index=NULL,r2_index=NULL,r3_index=NULL,D0=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) group=rep(1,p)
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D3)) {
    set.seed(2)
    S3 <- matrix(runif(r10*r20*r30,5,10),nrow = r30)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r30),nrow = q)
    U3 <- qr.Q(qr(T1))
    D3 <- U3%*%S3%*%t(kronecker(U2,U1))
  }

  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  for(j in 1:ng){
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    f0 = f0 + bsbasefun(X1,K,degr)%*%t(D3)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f0 + eps*sigma2
  
  # initialize
  if(is.null(r1_index)) r1_index = 1:min(4,s)
  if(is.null(r2_index)) r2_index = 1:min(4,K)
  if(is.null(r3_index)) r3_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r30),q,r30)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30),r30,r10*r20)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r3_max),q,r3_max)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D3=D3,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,D0=D0))
}

##------------------------------------------------------------------------##
##----------------generate scenario II data in IMMAM----------------------##
# e.g. mydata=immam3.sim.fsin(200,5,100,10,rep(1:4,each=25))
immam3.sim.fsin <- function(n,q,p,s,group=NULL,r10=2,r20=2,r30=2,isfixedR=0,D2=NULL,
                            K=6,degr=3,sigma2=NULL,seed_id=NULL,
                            r1_index=NULL,r2_index=NULL,r3_index=NULL,D0=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) group=rep(1,p)
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D2)) {
    set.seed(2)
    S3 <- matrix(runif(r10*r20*r30,5,10),nrow = r30)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r30),nrow = q)
    U3 <- qr.Q(qr(T1))
    D3 <- U3%*%S3%*%t(kronecker(U2,U1))
    D2 <- TransferModalUnfoldingsT(D3,3,2,c(s,K,q))
  }
  
  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  for(j in 1:ng){
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    basefuns1 <- sin(2*pi*X1)
    basefuns2 <- cos(pi*X1)
    f0 <- f0 + basefuns1%*%matrix(D2[1,],nrow=s)+basefuns2%*%matrix(D2[2,],nrow=s)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f0 + sigma2*eps
  
  # initialize
  if(is.null(r1_index)) r1_index = 1:min(4,s)
  if(is.null(r2_index)) r2_index = 1:min(4,K)
  if(is.null(r3_index)) r3_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r30),q,r30)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30),r30,r10*r20)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r3_max),q,r3_max)
      C <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
      SABC = list(S=S,A=A,B=B,C=C)
      D0 = list()
      for(j in 1:ng) D0[[j]] = SABC
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D3=D3,D2=D2,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,D0=D0))
}

##------------------------------------------------------------------------##
##--------------generate scenario I data in structural IMMAM--------------##
# e.g. mydata=immam4.sim.fbs(200,5,100,10,rep(1:4,each=25))
immam4.sim.fbs <- function(n,q,p,s,group=NULL,r10=2,r20=2,r30=2,r40=2,isfixedR=0,D44=NULL,
                           K=6,degr=3,sigma2=NULL,seed_id=NULL,
                           r1_index=NULL,r2_index=NULL,r3_index=NULL,r4_index=NULL,D0=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) group=c(rep(1,as.integer(p/2)),rep(2,p-as.integer(p/2)))
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(ng<r30) stop("ng must be not smaller than r30")
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D44)) {
    set.seed(2)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(ng*r30),nrow = ng)
    U3 <- as.matrix(qr.Q(qr(T1)))
    T1 <- matrix(rnorm(q*r40),nrow = q)
    U4 <- qr.Q(qr(T1)) 
    S4 <- matrix(runif(r10*r20*r30*r40,5,10),nrow = r30)
    D44 <- U4%*%S4%*%t(kronecker(U3,kronecker(U2,U1)))
  }

  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  for(j in 1:ng){
    D3 = D44[,((j-1)*s*K+1):(j*s*K)]
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    f0 = f0 + bsbasefun(X1,K,degr)%*%t(D3)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f0 + eps*sigma2
  
  # initialize
  if(is.null(r1_index)) r1_index = 1:min(4,s)
  if(is.null(r2_index)) r2_index = 1:min(4,K)
  if(is.null(r3_index)) r3_index = 1:min(4,ng)
  if(is.null(r4_index)) r4_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r30),ng,r30)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r40),q,r40)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30*r40),r40,r10*r20*r30)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      r4_max = max(r4_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r3_max),ng,r3_max)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r4_max),q,r4_max)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max*r4_max),r4_max,r1_max*r2_max*r3_max)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D44=D44,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,r40=r40,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,r4_index=r4_index,D0=D0))
}

##------------------------------------------------------------------------##
##-------------generate scenario II data in structural IMMAM--------------##
# e.g. mydata=immam4.sim.fsin(200,5,100,10,rep(1:4,each=25))
immam4.sim.fsin <- function(n,q,p,s,group=NULL,r10=2,r20=2,r30=2,r40=2,isfixedR=0,D42=NULL,
                            K=6,degr=3,sigma2=NULL,seed_id=NULL,
                            r1_index=NULL,r2_index=NULL,r3_index=NULL,r4_index=NULL,D0=NULL){
  if(n<2) stop("n must be not smaller than 2")
  if(q<1) stop("q must be not smaller than 1")
  if(p<1) stop("p must be not smaller than 1")
  if(s<1) stop("s must be not smaller than 1")
  if(is.null(group)) group=c(rep(1,as.integer(p/2)),rep(2,p-as.integer(p/2)))
  gunique <- unique(group)
  ng = length(gunique)
  pg = p/ng
  if(ng<r30) stop("ng must be not smaller than r30")
  if(is.null(sigma2)) sigma2 = 0.1
  if(is.null(seed_id)) seed_id=1000
  
  if(is.null(D42)) {
    set.seed(2)
    T1 <- matrix(rnorm(s*r10),nrow = s)
    U1 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(K*r20),nrow = K)
    U2 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(ng*r30),nrow = ng)
    U3 <- qr.Q(qr(T1))
    T1 <- matrix(rnorm(q*r40),nrow = q)
    U4 <- qr.Q(qr(T1)) 
    S4 <- matrix(runif(r10*r20*r30*r40,5,10),nrow = r30)
    D44 <- U4%*%S4%*%t(kronecker(U3,kronecker(U2,U1)))
    D42 = TransferModalUnfoldingsT(D44,4,2,c(s,K,ng,q))
  }
  
  set.seed(seed_id)
  X <- matrix(runif(n*p), nrow = n)
  
  f0 = matrix(0,n,q)
  id = NULL
  for(j in 1:q) id = c(id, c(1:s)+(j-1)*s*ng)
  for(j in 1:ng){
    D2 = D42[,id+(j-1)*s]
    Xj <- X[,group==gunique[j]]
    X1 <- Xj[,1:s]
    basefuns1 <- sin(2*pi*X1)
    basefuns2 <- cos(pi*X1)
    f0 <- f0 + basefuns1%*%matrix(D2[1,],nrow=s)+basefuns2%*%matrix(D2[2,],nrow=s)
  }
  eps <- matrix(rnorm(n*q),n,q)
  Y <- f0 + sigma2*eps
  
  # initialize
  if(is.null(r1_index)) r1_index = 1:min(4,s)
  if(is.null(r2_index)) r2_index = 1:min(4,K)
  if(is.null(r3_index)) r3_index = 1:min(4,ng)
  if(is.null(r4_index)) r4_index = 1:min(4,q)
  if(is.null(D0)) {
    set.seed(10)
    if(isfixedR){
      T1 = matrix(rnorm(pg*r10),pg,r10)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K*r20),K,r20)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r30),ng,r30)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r40),q,r40)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r10*r20*r30*r40),r40,r10*r20*r30)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }else{
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      r4_max = max(r4_index) 
      K_max = K
      T1 = matrix(rnorm(pg*r1_max),pg,r1_max)
      A <- qr.Q(qr(T1))
      T1 = matrix(rnorm(K_max*r2_max),K_max,r2_max)
      B <- qr.Q(qr(T1))
      T1 = matrix(rnorm(ng*r3_max),ng,r3_max)
      C <- qr.Q(qr(T1))
      T1 = matrix(rnorm(q*r4_max),q,r4_max)
      D <- qr.Q(qr(T1))
      S = matrix(runif(r1_max*r2_max*r3_max*r4_max),r4_max,r1_max*r2_max*r3_max)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }
  }
  
  return(list(Y=Y,X=X,f0=f0,group=group,D44=D44,D42=D42,n=n,q=q,p=p,s=s,r10=r10,r20=r20,r30=r30,r40=r40,K=K,degr=degr,sigma2=sigma2,seed_id=seed_id,
              r1_index=r1_index,r2_index=r2_index,r3_index=r3_index,r4_index=r4_index,D0=D0))
}








