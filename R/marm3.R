marm3 <- 
  function(Y,X,group=NULL,K=6,r1=NULL,r2=NULL,r3=NULL,method="BIC",ncv=10,penalty="LASSO",lambda=NULL,D0=NULL,
           intercept=TRUE,degr=3,nlam=20,lam_min=0.01,eps=1e-4,max_step=20, eps1=1e-4,max_step1=20,gamma=2,dfmax=NULL,alpha=1){
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    nx <- ncol(X)
    if(is.null(group)) group = rep(1,nx)
    gunique <- unique(group)
    G = length(gunique)
    p = rep(0,G)
    for(g in 1:G) p[g] = sum(group==gunique[g])
    if(degr>K-1) stop("K must be larger than degree+1 !")
    if(is.null(r1)) r1 <- 2 
    if(is.null(r2)) r2 <- 2
    if(is.null(r3)) r3 <- 2
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MC penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = nx + 1
    # initial A,B,C,S
    if(is.null(D0)){
      set.seed(1)
      B = rbind(diag(r2), matrix(0,K-r2,r2))
      C = rbind(diag(r3), matrix(0,q-r3,r3))
      S = matrix(rnorm(r1*r2*r3),r3,r1*r2)
      D0 = list(); 
      for(j in 1:G){
        A = rbind(diag(r1), matrix(0,p[j]-r1,r1))
        SABC = list(S=S,A=A,B=B,C=C)
        D0[[j]] = SABC
      }
    }
    opts = list(eps=eps,eps1=eps,max_step=max_step,max_step1=max_step,n=n,r1=r1,r2=r2,r3=r3,p=p,q=q,degr=degr,K=K,G=G,nx=nx)
    Sinit = list()
    Ainit = list()
    Binit = list()
    Cinit = list()
    Z = list()
    Zbar = NULL
    for(i in 1:G){
      Sinit[[i]] = D0[[i]]$S[1:r3,1:(r1*r2)]
      Ainit[[i]] = D0[[i]]$A[,1:r1]
      Binit[[i]] = D0[[i]]$B[,1:r2]
      Cinit[[i]] = D0[[i]]$C[,1:r3]
      if(dim(Sinit[[i]])[1]!=dim(Cinit[[i]])[2]) Sinit[[i]] = t(Sinit[[i]])
      Z[[i]] = as.matrix(bsbasefun(X[,group==gunique[i]],K,degr))
      Zbar1 = colMeans(Z[[i]])
      Z[[i]] = Z[[i]] - matrix(rep(Zbar1,each=n),n)
      Zbar = cbind(Zbar, Zbar1)
    }
    Ybar = colMeans(Y)
    Y1 = Y - matrix(rep(Ybar,each=n),n)
    if (is.null(lambda)) {
      is_setlam = 1
      if (nlam < 1||is.null(nlam)) stop("nlambda must be at least 1")
      setlam = c(1,lam_min,alpha,nlam)
      lambda = setuplambda(Y1,Z,Sinit,Ainit,Binit,Cinit,nx,G,nlam,setlam)
    }
    else {
      is_setlam = 0
      nlam = length(lambda)
      setlam = c(1,lam_min,alpha,nlam)
    }
    opts_pen = list(gamma=gamma,dfmax=dfmax,pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,alpha=alpha,isPenColumn=1) 
    #---------------- The selection by BIC or CV  ---------------------# 
    if(method=="BIC"){
      fit = EstPenColumn(Y1,Z,Sinit,Ainit,Binit,Cinit,lambda,opts,opts_pen)
      df = 0; for(g in 1:G) df = df + r1*r2*r3 + fit$df[g,]*r1 + K*r2 + q*r3 - r1^2 - r2^2 - r3^2
      bic = log(fit$likhd/(n*q)) + log(n*q)*df/(n*q)
      
      selected = which.min(bic)
      lambda_opt = lambda[selected]

      opts_pen$nlam = length(lambda[1:selected])
      fit = EstPenColumn(Y1,Z,Sinit,Ainit,Binit,Cinit,lambda[1:selected],opts,opts_pen)

      activeX = fit$betapath[,selected]
      Snew = fit$Snew[[selected]]
      Anew = fit$Anew[[selected]]
      Bnew = fit$Bnew[[selected]]
      Cnew = fit$Cnew[[selected]]
      Dn = NULL
      for (g in 1:G) Dn = cbind(Dn, Cnew[[g]] %*% Snew[[g]] %*%t(kronecker(Bnew[[g]], Anew[[g]])))
      if(intercept)  mu = Ybar-Dn%*%as.vector(Zbar)
      else mu = rep(0,q)
    }
    if(method=="CV"&&nlam>1){
      len_cv = ceiling(n/ncv)
      RSS = rep(0,nlam)
      for(jj in 1:ncv){
        cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
        if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
        Ytrain = Y1[-cv.id,]
        Xtrain = X[-cv.id,]
        Ytest = Y1[cv.id,]
        Xtest = X[cv.id,]
        Ztrain = list()
        Ztest = list()
        for(i in 1:G){
          Ztrain[[i]] = bsbasefun(Xtrain[,group==gunique[i]],K,degr)
          Ztest[[i]] = bsbasefun(Xtest[,group==gunique[i]],K,degr)
        } 
        fit = EstPenColumnCV(Ytrain,Ztrain,Ytest,Ztest,Sinit,Ainit,Binit,Cinit,lambda,opts,opts_pen)
        RSS = RSS + fit$likhd
      } 
      selected = which.min(RSS)
      lambda_opt = lambda[selected]
      
      opts_pen$nlam = length(lambda[1:selected])
      fit_opt = EstPenColumn(Y1,Z,Sinit,Ainit,Binit,Cinit,lambda[1:selected],opts,opts_pen)
      
      activeX = fit_opt$betapath[,selected]
      Snew = fit_opt$Snew[[selected]]
      Anew = fit_opt$Anew[[selected]]
      Bnew = fit_opt$Bnew[[selected]]
      Cnew = fit_opt$Cnew[[selected]]
      Dn = NULL
      for (g in 1:G) Dn = cbind(Dn, Cnew[[g]] %*% Snew[[g]] %*%t(kronecker(Bnew[[g]], Anew[[g]])))
      if(intercept)  mu = Ybar-Dn%*%as.vector(Zbar)
      else mu = rep(0,q)
    }
    if(method=="CV"&&nlam==1){
      fit = EstPenColumn(Y1,Z,Sinit,Ainit,Binit,Cinit,lambda,opts,opts_pen)
      selected = 1
      lambda_opt = lambda
      activeX = fit$betapath
      Snew = fit$Snew[[selected]]
      Anew = fit$Anew[[selected]]
      Bnew = fit$Bnew[[selected]]
      Cnew = fit$Cnew[[selected]]
      Dn = NULL
      for (g in 1:G) Dn = cbind(Dn, Cnew[[g]] %*% Snew[[g]] %*%t(kronecker(Bnew[[g]], Anew[[g]])))
      if(intercept)  mu = Ybar-Dn%*%as.vector(Zbar)
      else mu = rep(0,q)
    }

    return(list(D          = Dn,
                mu         = mu,
                S.opt      = Snew,
                A.opt      = Anew,
                B.opt      = Bnew,
                C.opt      = Cnew,
                lambda.seq = lambda,
                lambda.opt = lambda_opt,
                rss        = fit$likhd[selected],
                df         = fit$df,
                activeX    = activeX,
                opts       = opts,
                opts_pen   = opts_pen))
  }