
##--------------Estimation with Penalty by CV----------------------##
marm3.cv <- 
  function(Y,X,group,ncv,K_index,r1_index,r2_index,r3_index,lambda,D0,intercept,opts,opts_pen){
    n = opts$n
    p = opts$p
    q = opts$q
    G = opts$G
    nlam = length(lambda)
    degr = opts$degr
    gunique <- unique(group)
    Ybar = colMeans(Y)
    Y1 = Y - matrix(rep(Ybar,each=n),n)
    Sinit = list()
    Ainit = list()
    Binit = list()
    Cinit = list()
    
    len_cv = ceiling(n/ncv)
    RSS = matrix(0,nlam,length(r1_index)*length(r2_index)*length(r3_index)*length(K_index))
    for(jj in 1:ncv){ # start CV
      cv.id = ((jj-1)*len_cv+1):(jj*len_cv)
      if(jj==ncv) cv.id = ((jj-1)*len_cv+1):n
      Ytrain = Y1[-cv.id,]
      Xtrain = X[-cv.id,]
      Ytest = Y1[cv.id,]
      Xtest = X[cv.id,]
      
      RSS0 = NULL
      for(K in K_index){
        opts$K = K
        Ztrain = list()
        Ztest = list()
        for(i in 1:G){
          Ztrain[[i]] = bsbasefun(Xtrain[,group==gunique[i]],K,degr)
          Ztest[[i]] = bsbasefun(Xtest[,group==gunique[i]],K,degr)
        } 
        for(r3 in r3_index){
          opts$r3 = r3
          for(r2 in r2_index){
            opts$r2 = r2
            for(r1 in r1_index){
              opts$r1 = r1
              for(i in 1:G){
                Sinit[[i]] = as.matrix(D0[[i]]$S[1:r3,1:(r1*r2)])
                Ainit[[i]] = as.matrix(D0[[i]]$A[,1:r1])
                Binit[[i]] = as.matrix(D0[[i]]$B[,1:r2])
                Cinit[[i]] = as.matrix(D0[[i]]$C[,1:r3])
                if(dim(Sinit[[i]])[1]!=dim(Cinit[[i]])[2]) Sinit[[i]] = t(Sinit[[i]])
              }
              fit = EstPenColumnCV(Ytrain,Ztrain,Ytest,Ztest,Sinit,Ainit,Binit,Cinit,lambda,opts,opts_pen)
              RSS0 = cbind(RSS0,fit$likhd)
         
            }
          }
        }
      }
      RSS = RSS + RSS0
    } # end CV
    selected = which.min(RSS)
    qj = ceiling(selected/nlam)
    qj1 = selected%%nlam
    if(qj1==0) qj1 = nlam
    
    lambda_opt = lambda[qj1]
    opt = assig(c(length(r1_index),length(r2_index),length(r3_index),length(K_index)))[,qj]
    r1_opt = r1_index[opt[1]]
    r2_opt = r2_index[opt[2]]
    r3_opt = r3_index[opt[3]]
    K_opt = K_index[opt[4]]
    
    #---------------- The estimation after selection ---------------------#
    opts$r1 = r1_opt
    opts$r2 = r2_opt
    opts$r3 = r3_opt
    opts$K = K_opt
    Z = list()
    Zbar = NULL
    for(i in 1:G){
      Sinit[[i]] = as.matrix(D0[[i]]$S[1:r3_opt,1:(r1_opt*r2_opt)])
      Ainit[[i]] = as.matrix(D0[[i]]$A[,1:r1_opt])
      Binit[[i]] = as.matrix(D0[[i]]$B[,1:r2_opt])
      Cinit[[i]] = as.matrix(D0[[i]]$C[,1:r3_opt])
      if(dim(Sinit[[i]])[1]!=dim(Cinit[[i]])[2]) Sinit[[i]] = t(Sinit[[i]])
      Z[[i]] = bsbasefun(X[,group==gunique[i]],K_opt,degr)
      Zbar1 = colMeans(Z[[i]])
      Z[[i]] = Z[[i]] - matrix(rep(Zbar1,each=n),n)
      Zbar = cbind(Zbar, Zbar1)
    }
    opts_pen$nlam = length(lambda[1:qj1])
    opts_pen$lam_max = max(lambda[1:qj1])
    opts_pen$lam_min = min(lambda[1:qj1])
    fit_opt = EstPenColumn(Y1,Z,Sinit,Ainit,Binit,Cinit,lambda[1:qj1],opts,opts_pen) 
    activeX = fit_opt$betapath[,qj1]
    Snew = fit_opt$Snew[[qj1]]
    Anew = fit_opt$Anew[[qj1]]
    Bnew = fit_opt$Bnew[[qj1]]
    Cnew = fit_opt$Cnew[[qj1]]
    Dn = NULL
    for (g in 1:G) Dn = cbind(Dn, Cnew[[g]] %*% Snew[[g]] %*%t(kronecker(Bnew[[g]], Anew[[g]])))
    if(intercept)  mu = Ybar-Dn%*%as.vector(Zbar)
    else mu = rep(0,q)
    return(list(D          = Dn,
                mu         = mu,
                S.opt      = Snew,
                A.opt      = Anew,
                B.opt      = Bnew,
                C.opt      = Cnew,
                rk.opt     = c(r1_opt,r2_opt,r3_opt,K_opt),
                lambda.seq = lambda,
                lambda.opt = lambda_opt,
                rss        = fit_opt$likhd[qj1],
                df         = fit_opt$df,
                activeX    = activeX,
                opts       = opts,
                opts_pen   = opts_pen))
}
