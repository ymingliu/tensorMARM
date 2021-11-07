
##--------------Estimation with Penalty by BIC----------------------##
marm4.bic <- 
  function(Y,X,group,K_index,r1_index,r2_index,r3_index,r4_index,lambda,D0,intercept,opts,opts_pen){
  n <- opts$n
  q <- opts$q
  p <- opts$p
  G = opts$G
  degr = opts$degr
  gunique <- unique(group)
  nlam = opts_pen$nlam
  Ybar = colMeans(Y)
  Y1 = Y - matrix(rep(Ybar,each=n),n)
  RSS = NULL
  for(K in K_index){
    Z = NULL;  for(i in 1:G)  Z = cbind(Z, bsbasefun(X[,group==gunique[i]],K,degr))
    Zbar = colMeans(Z)
    Z = Z - matrix(rep(Zbar,each=n),n)
    
    opts$K = K
    for(r4 in r4_index){
      opts$r4 = r4
      for(r3 in r3_index){
        opts$r3 = r3
        for(r2 in r2_index){
          opts$r2 = r2
          for(r1 in r1_index){
            opts$r1 = r1
            S = as.matrix(D0$S[1:r4,1:(r1*r2*r3)])
            A = as.matrix(D0$A[,1:r1])
            B = as.matrix(D0$B[,1:r2])
            C = as.matrix(D0$C[,1:r3])
            D = as.matrix(D0$D[,1:r4])
            if(dim(S)[1]!=dim(D)[2]) S=t(S)
            fit = EstPenColumnT4(Y1,Z,S,A,B,C,D,lambda,opts,opts_pen)
            df = r1*r2*r3*r4 + fit$df*r1 + K*r2 + G*r3 + q*r4 - r1^2 - r2^2 - r3^2 - r4^2
            RSS = c(RSS,log(fit$likhd/(n*q)) + log(n*q)*df/(n*q))
          }
        }
      }
    }
  }
  selected = which.min(RSS)
  qj = ceiling(selected/nlam)
  qj1 = selected%%nlam
  if(qj1==0) qj1=nlam
  
  lambda_opt = lambda[qj1]
  opt = assig(c(length(r1_index),length(r2_index),length(r3_index),length(r4_index),length(K_index)))[,qj]
  r1_opt = r1_index[opt[1]]
  r2_opt = r2_index[opt[2]]
  r3_opt = r3_index[opt[3]]
  r4_opt = r4_index[opt[4]]
  K_opt = K_index[opt[5]]
  
  #---------------- The estimation after selection ---------------------#
  Z = NULL;  for(i in 1:G)  Z = cbind(Z, bsbasefun(X[,group==gunique[i]],K_opt,degr))
  Zbar = colMeans(Z)
  Z = Z - matrix(rep(Zbar,each=n),n)
  
  S = as.matrix(D0$S[1:r4_opt,1:(r1_opt*r2_opt*r3_opt)])
  A = as.matrix(D0$A[,1:r1_opt])
  B = as.matrix(D0$B[,1:r2_opt])
  C = as.matrix(D0$C[,1:r3_opt])
  D = as.matrix(D0$D[,1:r4_opt])
  if(dim(S)[1]!=dim(D)[2]) S = t(S)
  
  opts$K = K_opt
  opts$r1 = r1_opt
  opts$r2 = r2_opt
  opts$r3 = r3_opt
  opts$r4 = r4_opt
  
  opts_pen$nlam = length(lambda[1:qj1])
  fit_opt = EstPenColumnT4(Y1,Z,S,A,B,C,D,lambda[1:qj1],opts,opts_pen)
  
  activeX = fit_opt$betapath[,qj1]
  Anew = matrix(fit_opt$Apath[,qj1],nrow=p)
  Bnew = matrix(fit_opt$Bpath[,qj1],nrow=K_opt)
  Cnew = matrix(fit_opt$Cpath[,qj1],nrow=G)
  Dnew = matrix(fit_opt$Dpath[,qj1],nrow=q)
  Snew = matrix(fit_opt$Spath[,qj1],nrow=r4_opt)
  Dn = Dnew %*% Snew %*%t(kronecker(Cnew, kronecker(Bnew,Anew)))
  if(intercept)  mu = Ybar-Dn%*%Zbar
  else mu = rep(0,q)
   
  return(list(D          = Dn,
              mu         = mu,
              S.opt      = Snew,
              A.opt      = Anew,
              B.opt      = Bnew,
              C.opt      = Cnew,
              D.opt      = Dnew,
              rk.opt     = c(r1_opt,r2_opt,r3_opt,r4_opt,K_opt),
              lambda.seq = lambda,
              lambda.opt = lambda_opt,
              rss        = fit_opt$likhd[qj1],
              df         = fit_opt$df,
              activeX    = activeX,
              opts       = opts,
              opts_pen   = opts_pen))
}
