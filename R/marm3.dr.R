
##--------------Estimation with Penalty by CV----------------------##
marm3.dr <- 
  function(Y,X,group=NULL,K_index=NULL,r1_index=NULL,r2_index=NULL,r3_index=NULL,method="BIC",ncv=10,penalty="LASSO",lambda=NULL,
           D0=NULL,intercept=TRUE,nlam=50,degr=3,lam_min=0.01,eps=1e-4,max_step=20,eps1=1e-4,max_step1=20,gamma=2,dfmax=NULL,alpha=1){
    n <- nrow(Y)
    q <- ncol(Y)
    nx <- ncol(X)
    if(is.null(group)) group = rep(1,nx)
    gunique <- unique(group)
    G = length(gunique)
    p = rep(0,G)
    for(g in 1:G) p[g] = sum(group==gunique[g])
    K1 <- 6
    if(degr>min(6,K1-1)-1) stop("K must be larger than degree+1 !")
    if(is.null(K_index)) K_index = min(6,K1-1):max(8,K1+1)
    if(is.null(r1_index)) r1_index = 1:min(floor(log(n)),min(p))
    if(is.null(r2_index)) r2_index = 1:min(K_index)
    if(is.null(r3_index)) r3_index = 1:min(floor(log(n)),q)
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MCP penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = nx + 1
    # initial A,B,C,S
    if(is.null(D0)){
      set.seed(1)
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      K_max = max(K_index)
      
      B = rbind(diag(r2_max), matrix(0,K_max-r2_max,r2_max))
      C = rbind(diag(r3_max), matrix(0,q-r3_max,r3_max))
      S = matrix(rnorm(r1_max*r2_max*r3_max),r3_max,r1_max*r2_max)
      D0 = list(); 
      for(j in 1:G){
        A = rbind(diag(r1_max), matrix(0,p[j]-r1_max,r1_max))
        SABC = list(S=S,A=A,B=B,C=C)
        D0[[j]] = SABC
      }
    }
    opts = list(eps=eps,eps1=eps,max_step=max_step,max_step1=max_step,n=n,r1=2,r2=2,r3=2,p=p,q=q,degr=degr,K=max(K_index),G=G,nx=nx)
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1||is.null(nlam)) stop("nlambda must be at least 1")
      setlam = c(1,lam_min,alpha,nlam)
      Sinit = list()
      Ainit = list()
      Binit = list()
      Cinit = list()
      Z = list()
      for(i in 1:G){
        Sinit[[i]] = D0[[i]]$S
        Ainit[[i]] = D0[[i]]$A
        Binit[[i]] = D0[[i]]$B
        Cinit[[i]] = D0[[i]]$C
        Z[[i]] = bsbasefun(X[,group==gunique[i]],max(K_index),degr)
        Zbar1 = colMeans(Z[[i]])
        Z[[i]] = Z[[i]] - matrix(rep(Zbar1,each=n),n)
      }
      Ybar = colMeans(Y)
      Y1 = Y - matrix(rep(Ybar,each=n),n)
      lambda = setuplambda(Y1,Z,Sinit,Ainit,Binit,Cinit,nx,G,nlam,setlam)
    }
    else {
      is_setlam = 0
      nlam = length(lambda)
      setlam = c(1,lam_min,alpha,nlam)
    }
    opts_pen = list(gamma=gamma,dfmax=dfmax,pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,alpha=alpha,isPenColumn=1) 
    #---------------- The selection by CV  ---------------------#  
    if(method=="BIC") fit_dr = marm3.bic(Y,X,group,K_index,r1_index,r2_index,r3_index,lambda,D0,intercept,opts,opts_pen)
    if(method=="CV") fit_dr = marm3.cv(Y,X,group,ncv,K_index,r1_index,r2_index,r3_index,lambda,D0,intercept,opts,opts_pen)
    
    return(fit_dr)
  }