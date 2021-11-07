
##--------------Estimation with Penalty by CV----------------------##
marm4.dr <- 
  function(Y,X,group,K_index=NULL,r1_index=NULL,r2_index=NULL,r3_index=NULL,r4_index=NULL,method="BIC",ncv=10,penalty="LASSO",
           lambda=NULL,D0=NULL,intercept=TRUE,nlam=20,degr=3,lam_min=0.01,eps=1e-4,max_step=10,eps1=1e-4,max_step1=10,gamma=2,dfmax=NULL,alpha=1){
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    nx <- ncol(X)
    if(is.null(group)) group = rep(1,nx)
    gunique <- unique(group)
    G = length(gunique)
    p <- nx/G
    K1 <- 6
    if(degr>min(6,K1-1)-1) stop("K must be larger than degree+1 !")
    if(is.null(K_index)) K_index = min(6,K1-1):max(8,K1+1)
    if(is.null(r1_index)) r1_index = 1:min(floor(log(n)),p)
    if(is.null(r2_index)) r2_index = 1:min(K_index)
    if(is.null(r3_index)) r3_index = 1:min(floor(log(n)),G)
    if(is.null(r4_index)) r4_index = 1:min(floor(log(n)),q)
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MCP penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = p + 1
    opts = list(eps=eps,eps1=eps1,max_step=max_step,max_step1=max_step1,n=n,r1=max(r1_index),
                r2=max(r2_index),r3=max(r3_index),r4=max(r4_index),p=p,q=q,G=G,K=6,degr=degr) 
    # initial A,B,C,S
    if(is.null(D0)){
      set.seed(1)
      r1_max = max(r1_index) 
      r2_max = max(r2_index) 
      r3_max = max(r3_index) 
      r4_max = max(r4_index) 
      K_max = max(K_index)
      A = rbind(diag(r1_max), matrix(0,p-r1_max,r1_max))
      B = rbind(diag(r2_max), matrix(0,K_max-r2_max,r2_max))
      C = rbind(diag(r3_max), matrix(0,G-r3_max,r3_max))
      D = rbind(diag(r4_max), matrix(0,q-r4_max,r4_max))
      S = matrix(rnorm(r1_max*r2_max*r3_max*r4_max),r4_max,r1_max*r2_max*r3_max)
      D0 = list(S=S,A=A,B=B,C=C,D=D)
    }
    else{
      A = D0$A
      B = D0$B
      C = D0$C
      D = D0$D
      S = D0$S
    }
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1||is.null(nlam)) stop("nlambda must be at least 1")
      if (n<=p) lam_min = 1e-1
      setlam = c(1,lam_min,alpha,nlam)
      Z = NULL;    for(i in 1:G)  Z = cbind(Z, bsbasefun(X[,group==gunique[i]],max(K_index),degr))
      Zbar = colMeans(Z)
      Z = Z - matrix(rep(Zbar,each=n),n)
      
      Ybar = colMeans(Y)
      Y1 = Y - matrix(rep(Ybar,each=n),n)
      
      lambda = setuplambdaT4(Y1,Z,S,A,B,C,D,nlam,setlam)
    }
    else {
      is_setlam = 0
      nlam = length(lambda)
      setlam = c(1,lam_min,alpha,nlam)
    }
    #---------------- The selection by CV  ---------------------#  
    if((max(r1_index)>dim(A)[2])|(max(r2_index)>dim(B)[2])|(max(r3_index)>dim(C)[2]))
      stop("maximum number of index sequence of r1, r2, and r3 must not be larger than A, B, and C, respectively !")
    opts_pen = list(gamma=gamma,dfmax=dfmax,pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,alpha=alpha,isPenColumn=TRUE) 
    if(method=="BIC") fit_dr = marm4.bic(Y,X,group,K_index,r1_index,r2_index,r3_index,r4_index,lambda,D0,intercept,opts,opts_pen)
    if(method=="CV") fit_dr = marm4.cv(Y,X,group,ncv,K_index,r1_index,r2_index,r3_index,r4_index,lambda,D0,intercept,opts,opts_pen)
    
    return(fit_dr)
  }
