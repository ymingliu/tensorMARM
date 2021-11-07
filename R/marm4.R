marm4 <- 
  function(Y,X,group=NULL,K=6,r1=NULL,r2=NULL,r3=NULL,r4=NULL,method="BIC",ncv=10,penalty="LASSO",lambda=NULL,D0=NULL,
           intercept=TRUE,nlam=20,degr=3,lam_min=0.01, eps=1e-4,max_step=10,eps1=1e-4,max_step1=10,gamma=2,dfmax=NULL,alpha=1){
    n <- dim(Y)[1]
    q <- dim(Y)[2]
    nx <- ncol(X)
    if(is.null(group)) group = rep(1,nx)
    gunique <- unique(group)
    G = length(gunique)
    p <- nx/G
    if(degr>K-1) stop("K must be larger than degree+1 !")
    if(is.null(r1)) r1 <- 2 
    if(is.null(r2)) r2 <- 2
    if(is.null(r3)) r3 <- 2
    if(is.null(r4)) r4 <- 2
    if (penalty == "LASSO") pen <- 1
    if (penalty == "MCP")   pen <- 2 
    if (penalty=="SCAD"){    
      gamma <- 3
      pen <- 3
    }  
    if (gamma <= 1 & penalty=="MCP") stop("gamma must be greater than 1 for the MCP penalty")
    if (gamma <= 2 & penalty=="SCAD") stop("gamma must be greater than 2 for the SCAD penalty")
    
    if (is.null(dfmax)) dfmax = p + 1
    
    opts = list(eps=eps,eps1=eps1,max_step=max_step,max_step1=max_step1,
                n=n,r1=r1,r2=r2,r3=r3,r4=r4,p=p,q=q,G=G,K=K,degr=degr)  
    # initial A,B,C,S
    if(is.null(D0)){
      set.seed(1)
      A = rbind(diag(r1), matrix(0,nx-r1,r1))
      B = rbind(diag(r2), matrix(0,K-r2,r2))
      C = rbind(diag(r3), matrix(0,G-r3,r3))
      D = rbind(diag(r4), matrix(0,q-r4,r4))
      S = matrix(rnorm(r1*r2*r3*r4),r4,r1*r2*r3)
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
      Z = NULL;    for(i in 1:G)  Z = cbind(Z, bsbasefun(X[,group==gunique[i]],K,degr))
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
    opts_pen = list(gamma=gamma,dfmax=dfmax,pen=pen,nlam=nlam,lam_max=1,lam_min=lam_min,alpha=alpha,isPenColumn=TRUE) 
    #---------------- The selection by BIC or CV  ---------------------# 
    Z = NULL;  for(i in 1:G)  Z = cbind(Z, bsbasefun(X[,group==gunique[i]],K,degr))
    Zbar = colMeans(Z)
    Z = Z - matrix(rep(Zbar,each=n),n)
    Ybar = colMeans(Y)
    Y1 = Y - matrix(rep(Ybar,each=n),n)
    if(method=="BIC"){
      fit = EstPenColumnT4(Y1,Z,as.matrix(S),as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(D),lambda,opts,opts_pen)
      df = r1*r2*r3*r4 + fit$df*r1 + K*r2 + G*r3 + q*r4 - r1^2 - r2^2 - r3^2 - r4^2
      bic = log(fit$likhd/(n*q)) + log(n*q)*df/(n*q)
      
      selected = which.min(bic)
      lambda_opt = lambda[selected]
      
      opts_pen$nlam = length(lambda[1:selected])
      fit = EstPenColumnT4(Y1,Z,as.matrix(S),as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(D),lambda[1:selected],opts,opts_pen)
      
      activeX = fit$betapath[,selected]
      Anew=matrix(fit$Apath[,selected],nrow=p)
      Bnew=matrix(fit$Bpath[,selected],nrow=K)
      Cnew=matrix(fit$Cpath[,selected],nrow=G)
      Dnew=matrix(fit$Dpath[,selected],nrow=q)
      Snew=matrix(fit$Spath[,selected],nrow=r4)
      Dn = Dnew %*% Snew %*%t(kronecker(Cnew, kronecker(Bnew,Anew)))
      if(intercept)  mu = Ybar-Dn%*%Zbar
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
        Ztrain = NULL;  for(i in 1:G)  Ztrain = cbind(Ztrain, bsbasefun(Xtrain[,group==gunique[i]],K,degr))
        Ztest = NULL;  for(i in 1:G)  Ztest = cbind(Ztest, bsbasefun(Xtest[,group==gunique[i]],K,degr))
        fit = EstPenColumnT4CV(Ytrain,Ztrain,Ytest,Ztest,as.matrix(S),as.matrix(A),as.matrix(B),
                               as.matrix(C),as.matrix(D),lambda,opts,opts_pen)
        RSS = RSS + fit$likhd
      }
      selected = which.min(RSS)
      lambda_opt = lambda[selected]
      
      opts_pen$nlam = length(lambda[1:selected])
      fit = EstPenColumnT4(Y1,Z,as.matrix(S),as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(D),lambda[1:selected],opts,opts_pen)
      
      activeX = fit$betapath[,selected]
      Anew=matrix(fit$Apath[,selected],nrow=p)
      Bnew=matrix(fit$Bpath[,selected],nrow=K)
      Cnew=matrix(fit$Cpath[,selected],nrow=G)
      Dnew=matrix(fit$Dpath[,selected],nrow=q)
      Snew=matrix(fit$Spath[,selected],nrow=r4)
      Dn = Dnew %*% Snew %*%t(kronecker(Cnew, kronecker(Bnew,Anew)))
      if(intercept)  mu = Ybar-Dn%*%Zbar
      else mu = rep(0,q)
    }
    if(method=="CV"&&nlam==1){
      fit = EstPenColumnT4(Y1,Z,as.matrix(S),as.matrix(A),as.matrix(B),as.matrix(C),as.matrix(D),lambda_opt,opts,opts_pen)
      selected = 1
      lambda_opt = lambda
      activeX = fit$betapath
      Anew=matrix(fit$Apath,nrow=p)
      Bnew=matrix(fit$Bpath,nrow=K)
      Cnew=matrix(fit$Cpath,nrow=G)
      Dnew=matrix(fit$Dpath,nrow=q)
      Snew=matrix(fit$Spath,nrow=r4)
      Dn = Dnew %*% Snew %*%t(kronecker(Cnew, kronecker(Bnew,Anew)))
      if(intercept)  mu = Ybar-Dn%*%Zbar
      else mu = rep(0,q)
    }

    return(list(D          = Dn,
                mu         = mu,
                S.opt      = Snew,
                A.opt      = Anew,
                B.opt      = Bnew,
                C.opt      = Cnew,
                D.opt      = Dnew,
                lambda.seq = lambda,
                lambda.opt = lambda_opt,
                rss        = fit$likhd[selected],
                df         = fit$df,
                activeX    = activeX,
                opts       = opts,
                opts_pen   = opts_pen))
  }