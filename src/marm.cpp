//[[Rcpp::depends(RcppEigen)]]

#include "marm.h"

//***------------------------------------------------------------------------------------**
//***--------------------------------updateS in T3---------------------------------------**
MatrixXd updateS(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C)
{
    int r1 = A.cols();
    int r2 = B.cols();
    int r3 = C.cols();
	int d = r1 * r2*r3, k,k1,j,j1;
	MatrixXd ztilde = Z * kroneckerProduct(B, A), vectorS;
	VectorXd U;
	U.setZero(d);
	MatrixXd  V = MatrixXd::Constant(d, d, 0);

	for (k = 0; k < r1*r2; k++) {
		for (j = 0; j < r3; j++) {
			U[k*r3 + j] = ztilde.col(k).transpose()*Y*C.col(j);
			for (k1 = 0; k1 < r1*r2; k1++) {
				for (j1 = 0; j1 < r3; j1++) {
					V(k*r3 + j, k1*r3 + j1) = kroneckerProduct(
						ztilde.col(k1).array()*ztilde.col(k).array(),
						(C.col(j1).array()*C.col(j).array()).transpose()).sum();
				}
			}
		}
	}
	vectorS = V.colPivHouseholderQr().solve(U);
	vectorS.resize(r3, r1*r2);
	return vectorS;
}

//***------------------------------------------------------------------------------------**
//***---------------------------------updateC in T3--------------------------------------**
MatrixXd updateC(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	int q = opts.q,j,kp;
	MatrixXd ztilde = Z * kroneckerProduct(B, A);
	MatrixXd StZ = ztilde * S.transpose();
	MatrixXd Cnew = MatrixXd::Constant(q, opts.r3, 0);
	HouseholderQR<MatrixXd> qr;
	qr.compute(StZ);
	MatrixXd R = qr.matrixQR().triangularView<Upper>();
	MatrixXd Q = qr.householderQ();

	
	kp = StZ.cols();
	MatrixXd temp, IDEN = MatrixXd::Identity(kp, kp);
	if (pow(condition_numberQRSym(R),2) > 1e10 || isnan(UpTriangularInv(R).sum())){
        temp = tRbyR(R).block(0, 0, kp, kp) + (IDEN.array()*1e-4).matrix();
	  for (j = 0; j < q; j++) Cnew.row(j) = (temp.colPivHouseholderQr().solve(StZ.transpose()*Y.col(j))).transpose();
	}
	else
	  for (j = 0; j < q; j++) Cnew.row(j) = (QbyR(Q.transpose(), UpTriangularInv(R), 0)*Y.col(j)).transpose();
	
	return Cnew;
}

//***------------------------------------------------------------------------------------**
//***----------------------------------updateA in T3-------------------------------------**
MatrixXd updateA(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	int r1 = opts.r1, r2 = opts.r2, K = opts.K, p = A.rows();
	int d = r1 * p,t1,t2,t3,t4,j;
	MatrixXd W = C * S, Wt1, Zt2, Wt3, Zt4, zbw1, zbw2, vectorA;
	VectorXd tU;
	tU.setZero(d);
	MatrixXd tV = MatrixXd::Constant(d, d, 0);

	for (t2 = 0; t2<p; t2++) {
		Zt2 = Z.col(t2);
		for (j = 1; j<K; j++)	Zt2 = cbind_rcpp(Zt2, Z.col(j*p + t2));
		for (t1 = 0; t1<r1; t1++) {
			Wt1 = W.col(t1);
			for (j = 1; j<r2; j++) Wt1 = cbind_rcpp(Wt1, W.col(j*r1 + t1));			
			zbw1 = Zt2 * B*(Wt1.transpose());
			tU[t2*r1 + t1] = (Y.array()*zbw1.array()).sum();

			for (t4 = 0; t4<p; t4++) {
				Zt4 = Z.col(t4);
				for (j = 1; j<K; j++) Zt4 = cbind_rcpp(Zt4, Z.col(j*p + t4));
				for (t3 = 0; t3<r1; t3++) {
					Wt3 = W.col(t3);
					for (j = 1; j<r2; j++)	Wt3 = cbind_rcpp(Wt3, W.col(j*r1 + t3));
					zbw2 = Zt4 * B*(Wt3.transpose());
					tV(t2*r1 + t1, t4*r1 + t3) = (zbw1.array()*zbw2.array()).sum();
				}
			}
		}
	}
	vectorA = tV.colPivHouseholderQr().solve(tU);
	vectorA.resize(r1, p);
	return vectorA.transpose();
}

//***------------------------------------------------------------------------------------**
//***--------------------------------updateB in T3---------------------------------------**
MatrixXd updateB(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S)
{
	int r1 = opts.r1, r2 = opts.r2, K = opts.K, q = opts.q, n = opts.n, p = A.rows();
	int d = r2 * K, t1, t2, t3, t4;
	MatrixXd W = C * S, Wt1, Zt2, Wt3, Zt4, zaw1, zaw2, vectorB;
	VectorXd tU;
	tU.setZero(d);
	MatrixXd tV = MatrixXd::Constant(d, d, 0);
	for (t2 = 0; t2<K; t2++) {
		Zt2 = Z.block(0, t2*p, n, p); 
		for (t1 = 0; t1<r2; t1++) {
			Wt1 = W.block(0, t1*r1, q, r1);			 
			zaw1 = Zt2 * A*(Wt1.transpose());
			tU[t2*r2 + t1] = (Y.array()*zaw1.array()).sum();
			for (t4 = 0; t4<K; t4++) {
				Zt4 = Z.block(0, t4*p, n, p); 
				for (t3 = 0; t3<r2; t3++) {
					Wt3 = W.block(0, t3*r1, q, r1);
					zaw2 = Zt4 * A*(Wt3.transpose()); 
					tV(t2*r2 + t1, t4*r2 + t3) = (zaw1.array()*zaw2.array()).sum();
				}
			}
		}
	}
    vectorB = tV.colPivHouseholderQr().solve(tU);
	vectorB.resize(r2, K);
	return vectorB.transpose();
}

//***------------------------------------------------------------------------------------**
//***------------------------------Estimation without penalty----------------------------**
List Estimation_Dj(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, VectorXi &convergence1, int k, double &likhd1)
{	
    double likhd0 = likhd1;
	MatrixXd Dnew, Anew, Bnew, Cnew, Snew;
	Snew = updateS(Y, Z, A, B, C);
	Dnew = C * Snew*kroneckerProduct(B.transpose(), A.transpose());
	likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
	if (likhd1<likhd0) {
		S = Snew;
		likhd0 = likhd1;
	}
	else convergence1[4*k]=0;	
	Cnew = updateC(Y, Z, A, B, C, S);
	Dnew = Cnew * S*kroneckerProduct(B.transpose(), A.transpose());
	likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
	if (likhd1<likhd0) {
		C = Cnew;
		likhd0 = likhd1;
	}
	else convergence1[4*k+1]=0;
	Bnew = updateB(Y, Z, A, B, C, S);
	Dnew = C * S*kroneckerProduct(Bnew.transpose(), A.transpose());
	likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
	if (likhd1<likhd0) {
		B = Bnew;
		likhd0 = likhd1;		
	}
	else convergence1[4*k+2]=0;	
	Anew = updateA(Y, Z, A, B, C, S);
	Dnew = C * S*kroneckerProduct(B.transpose(), Anew.transpose());
	likhd1 = (Y - Z * Dnew.transpose()).squaredNorm();
	if (likhd1<likhd0) {
		A = Anew;
		likhd0 = likhd1;
	}
	else convergence1[4*k+3]=0;		
	likhd1 = likhd0;
	return List::create(Named("S") = S, Named("A") = A, Named("B") = B, Named("C") = C);
}

// [[Rcpp::export]]
List Estimation(MatrixXd Y, List Z0, List Sinit, List Ainit, List Binit, List Cinit, List optsList)
{
    opts.eps = as<double>(optsList["eps"]);
    opts.max_step = as<int>(optsList["max_step"]);
    opts.eps1 = as<double>(optsList["eps1"]);
    opts.max_step1 = as<int>(optsList["max_step1"]);
    opts.r1 = as<int>(optsList["r1"]);
    opts.r2 = as<int>(optsList["r2"]);
    opts.r3 = as<int>(optsList["r3"]);
    opts.q = as<int>(optsList["q"]);
    opts.degr = as<int>(optsList["degr"]);
    opts.K = as<int>(optsList["K"]);
    opts.G = as<int>(optsList["G"]);
    opts.nx = as<int>(optsList["nx"]);
    opts.n = Y.rows();
    
    int i, j, ng = opts.G;
    double  likhd0 = 2*pow(10, 6), likhd1 = likhd0;
    MatrixXd S, A, B, C, Dnew, D0, Ynew = Y, Zj;
    VectorXi convergence1 = VectorXi::Constant(4*ng, 1);
    List fit, S0(ng), A0(ng), B0(ng), C0(ng), Z(ng);
    for(i=0; i<ng; i++){
        S0[i] = as<MatrixXd>(Sinit[i]);
        A0[i] = as<MatrixXd>(Ainit[i]);
        B0[i] = as<MatrixXd>(Binit[i]);
        C0[i] = as<MatrixXd>(Cinit[i]);
        Z[i] = as<MatrixXd>(Z0[i]);
    }
    int step = 0;
    for (i=0; i<ng; i++) {
        Zj = Z[i];
        S = S0[i]; A = A0[i]; B = B0[i]; C = C0[i];
        D0 = C * S * kroneckerProduct(B, A).transpose();
        Ynew = Ynew - Zj*D0.transpose();
    }
    while(step < opts.max_step){
        for(j=0;j<4*ng;j++) convergence1[j] = 1;
        step++;
        for(i=0; i<ng; i++){
            Zj = Z[i];
            S = S0[i]; A = A0[i]; B = B0[i]; C = C0[i];
            D0 = C * S * kroneckerProduct(B, A).transpose();
            Ynew = Ynew + Zj*D0.transpose();
            fit = Estimation_Dj(Ynew, Zj, S, A, B, C, convergence1, i, likhd1);
            S = fit["S"]; A = fit["A"]; B = fit["B"]; C = fit["C"];
            Dnew = C * S * kroneckerProduct(B, A).transpose();
            Ynew = Ynew - Zj * Dnew.transpose();
            S0[i] = S; A0[i] = A; B0[i] = B; C0[i] = C;
        }
        if(fabs(likhd0-likhd1)/likhd0 < opts.eps) break;
        if((D0-Dnew).norm()/D0.norm() < opts.eps) break;
        if(convergence1.sum()==0) break;
    }
    
    return List::create(Named("likhd") = likhd1, Named("Snew") = S0, Named("Anew") = A0, Named("Bnew") = B0, Named("Cnew") = C0);
}

//***------------------------------------------------------------------------------------**
//***------------------------------setup tuning parameters-------------------------------**
// [[Rcpp::export]]
VectorXd setuplambda(MatrixXd Y, List Z0, List Sinit, List Ainit, List Binit, List Cinit, int nx, int ng, int nlam, VectorXd setlam)
{
	int n = Y.rows(), q = Y.cols(), p, r1, r2, r3, K, j, jj, g,count=0;
	double lam_max, lam_min, alpha;
	VectorXd lambda, lambda1, tmp1, tmp = VectorXd::Constant(nx, 0);
	MatrixXd S,A,B,C,Z,S1, cbs, V, V_1, Y1 = Y, Gammaj, Gamma_sqrt, svdu, svdd, U;
	VectorXi id;
    VectorXi dims = VectorXi::Constant(3, 0);
	Y1.resize(n*q, 1);
	for(g=0;g<ng;g++){
		S = as<MatrixXd>(Sinit[g]);
		A = as<MatrixXd>(Ainit[g]);
		B = as<MatrixXd>(Binit[g]);
		C = as<MatrixXd>(Cinit[g]);
		Z = as<MatrixXd>(Z0[g]);
		p = A.rows();
		r1 = A.cols(); 
		r2 = B.cols(); 
		r3 = C.cols();
		K = B.rows();
        dims[0] = r1; dims[1] = r2; dims[2] = r3;
		S1 = TransferModalUnfoldingsT(S, 3, 1, dims);
		cbs = kroneckerProduct(C, B)*(S1.transpose());
		id = SEQ(1, p*K, p);		
		for (j = 0; j < p; j++)
		{
			V = submatrix_col(Z, id.array() + j)*(cbs.block(0, 0, K, r1));
			for (jj = 1; jj < q; jj++) {
				V_1 = submatrix_col(Z, id.array() + j)*(cbs.block(jj*K, 0, K, r1));
				V = rbind_rcpp(V, V_1);
			}
			Gammaj = ((V.transpose()*V).array() / n).matrix();		
			JacobiSVD<MatrixXd> svd(Gammaj, ComputeThinU | ComputeThinV);
			svdu = svd.matrixU();
			svdd = (svd.singularValues()).asDiagonal();
			Gamma_sqrt = svdu * ((1 / (svdd.diagonal().array().sqrt())).matrix().asDiagonal())*(svdu.transpose());	
			tmp1 = Y1.transpose()*V * Gamma_sqrt;
			tmp[count++]=tmp1.array().abs().sum();
		}
	}
	lam_max = setlam[0];
	lam_min = setlam[1];
	alpha = setlam[2];

	double max_tmp;
	max_tmp = (tmp.array()).maxCoeff()/sqrt(n*q*q);
	double max_lam;
	max_lam = lam_max * max_tmp / alpha;
	if (lam_min == 0) {
		lambda1.setLinSpaced(nlam - 1, log(max_lam), log(0.0001*max_lam));
		lambda.setLinSpaced(nlam, 0, 0);
		lambda.segment(0, nlam - 1) = lambda1.array().exp().matrix();
	}
	else {
		lambda1.setLinSpaced(nlam, log(max_lam), log(lam_min*max_lam));
		lambda = lambda1.array().exp();
	}
	return lambda;
}

//***------------------------------------------------------------------------------------**
//***---------------------update the jth row of matrix A with penalty--------------------**
MatrixXd updateA_penalty(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd S, VectorXi &activeA, double lambda1)
{
    /*
	Input:
	Y is n*q matrix
	Z is n*(K*p) matrix
	A is p*r1 matrix
	B is K*r2 matrix
	C is q*r3 matrix
	S is r3*(r1*r2) matrix
	penalty= 1(lasso),2(mcp), and 3(scad)
	lambda is preset tuning parameters, a L-vector
	alpha is tuning parameter being in (0,1) to control elasnet
	eps is error to control convergence
	nlam is the number of preset tuning parameters
	max_iter is the maxixum number of iteration
    */
    int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, K = opts.K, q = opts.q, n = opts.n, p = A.rows();
	int gamma = opts_pen.gamma, pen = opts_pen.pen;
	double alpha = opts_pen.alpha, eps1 = opts.eps1;
	int max_step1 = opts.max_step1, j,jj, active;
    int dfmax = opts_pen.dfmax;

    MatrixXd S1;
    VectorXd aj,ajnew,zj;
    VectorXi dims = VectorXi::Constant(3, 0);
    dims[0] = r1; dims[1] = r2; dims[2] = r3;
    S1 = TransferModalUnfoldingsT(S, 3, 1, dims);
    VectorXi id = SEQ(1, p*K, p);
    MatrixXd Vnew = MatrixXd::Constant(n*q, r1*p, 0);
    MatrixXd Gamma_sqrtn = MatrixXd::Constant(r1, r1*p, 0);
    MatrixXd cbs = kroneckerProduct(C, B)*(S1.transpose());
    MatrixXd A1 = A, V, Gammaj, Gamma_sqrt, V_1, D3,L;
    D3 = C * S*kroneckerProduct(B.transpose(), A.transpose());
    MatrixXd IDEN = MatrixXd::Identity(r1, r1);
    int count=0;
    
	for (j = 0; j < p; j++) {
		V = submatrix_col(Z, id.array() + j)*(cbs.block(0, 0, K, r1));
		for (jj = 1; jj < q; jj++) {
			V_1 = submatrix_col(Z, id.array() + j)*(cbs.block(jj*K, 0, K, r1));
			V = rbind_rcpp(V, V_1);
		}
		Gammaj = ((V.transpose()*V).array() / n).matrix();
		L = Gammaj.llt().matrixL();
        if (isnan(UpTriangularInv(L.transpose()).sum())) {
            count++;
            Gammaj = ((V.transpose()*V).array() / n).matrix() + (IDEN.array()*1e-4).matrix();
            L = Gammaj.llt().matrixL();
        }
		Gamma_sqrtn.block(0, j*r1, r1, r1) = UpTriangularInv(L.transpose());
		Vnew.block(0, j*r1, n*q, r1) = QbyR(V, UpTriangularInv(L.transpose()),1);
		A1.row(j) = (QbyR(A.row(j).transpose(), L.transpose(), 0)).transpose();
	}

	MatrixXd Anew = A1;
	MatrixXd r = Y - Z * (D3.transpose());
	r.resize(n*q, 1);
	VectorXd ajnorm_old, ajnorm;
	ajnorm_old = ajnorm = VectorXd::Constant(p, 0);
	double converged1;
	int step = 0;
	while (step<max_step1) 
	{
		step++;
		active = 0;
		for (j = 0; j < p; j++)
			if (ajnorm[j] != 0) active = active + 1;
        if (active>dfmax) {
            return A;
        }
		for (j = 0; j < p;j++) {
			aj = Anew.row(j).transpose();
			zj = Vnew.block(0, j*r1, n*q, r1).transpose()*r/n + aj;
            if (zj.norm()==0) {
                ajnew = aj;
            }else{
                ajnew = updateAj(zj, n, r1, lambda1, alpha, gamma, pen);
            }
            if (isnan(ajnew.norm())) {
                stop("error: ajnew is nan!");
            }
            r = r - Vnew.block(0, j*r1, n*q, r1)*(ajnew-aj);
			Anew.row(j) = ajnew.transpose();
			ajnorm[j] = ajnew.norm();
		}
		converged1 = 1;
		for (j = 0; j < p;j++) {
				if (ajnorm[j] != 0 && ajnorm_old[j] != 0) {
					if ((A1.row(j) - Anew.row(j)).norm() / ajnorm_old[j]>eps1) {
						converged1 = 0; break;
					}
				}
				else if (ajnorm[j] == 0 && ajnorm_old[j] != 0) {
					converged1 = 0; break;
				}
				else if (ajnorm[j] != 0 && ajnorm_old[j] == 0) {
					converged1 = 0; break;
				}
			}
		if (converged1) break;
		A1 = Anew;
		ajnorm_old = ajnorm;
    }//end while
	for (j = 0; j<p; j++) {
		Anew.row(j) = (QbyR(Anew.row(j).transpose(),Gamma_sqrtn.block(0, j*r1, r1, r1),0)).transpose();
		if (ajnorm[j]) activeA[j] = 1;
	}
	return Anew;
}

//----------------------------------------------------------------**
//***--------------------Estimation without penalty---------------**
List Estimation_pen_Dj(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C,
                       VectorXi &convergence1, int k, double &likhd1, double lambda1, MatrixXd Zm, MatrixXd Am, VectorXi &breakid)
{
    int j, p = A.rows(), r1 = A.cols(), K = opts.K;
    VectorXi activeA = VectorXi::Constant(p, 0);
    double likhd0 = likhd1;
    MatrixXd Dnew, Anew, Bnew, Cnew, Snew, Z1, A1, Z2, A2;
    Anew = A;
    A1 = Am;
    Z1 = Zm;
    
    Snew = updateS(Y, Z1, A1, B, C);
    Dnew = C * Snew * kroneckerProduct(B, A1).transpose();
    likhd1 = (Y - Z1 * Dnew.transpose()).squaredNorm();
    if (likhd1<likhd0) {
        S = Snew;
        likhd0 = likhd1;
    }
    else convergence1[4*k]=0;
    
    Cnew = updateC(Y, Z1, A1, B, C, S);
    Dnew = Cnew * S*kroneckerProduct(B, A1).transpose();
    likhd1 = (Y - Z1 * Dnew.transpose()).squaredNorm();
    if (likhd1<likhd0) {
        C = Cnew;
        likhd0 = likhd1;
    }
    else convergence1[4*k+1]=0;
    
    for(j=0;j<p;j++) activeA[j] = 0;
    Anew = updateA_penalty(Y, Z, A, B, C, S, activeA, lambda1);
    if(activeA.sum()<r1){
        A1 = A;
        Z1 = Z;
        breakid[k] = 1;
    }
    else{
        Z2 = extractColsZ(Z,p,K,activeA);
        A2 = extractRows(Anew, activeA);
        Dnew = C * S * kroneckerProduct(B.transpose(), A2.transpose());
        likhd1 = (Y - Z2 * Dnew.transpose()).squaredNorm();
        if (likhd1<likhd0) {
            A = Anew;
            Z1 = Z2;
            A1 = A2;
            likhd0 = likhd1;
        }
        else convergence1[4*k+3]=0;
    }
    
    if (breakid[k]==0) {
        Bnew = updateB(Y, Z1, A1, B, C, S);
        Dnew = C * S*kroneckerProduct(Bnew, A1).transpose();
        likhd1 = (Y - Z1 * Dnew.transpose()).squaredNorm();
        if (likhd1<likhd0) {
            B = Bnew;
            likhd0 = likhd1;
        }
        else convergence1[4*k+2]=0;
    }
    
    likhd1 = likhd0;
    return List::create(Named("S") = S, Named("A") = A, Named("B") = B, Named("C") = C, Named("Zm") = Z1, Named("Am") = A1);
}

//***------------------------------------------------------------------------------------**
//***--------------Estimation with penalizing functions in a whole column ---------------**
// [[Rcpp::export]]
List EstPenColumn(MatrixXd Y, List Z0, List Sinit, List Ainit, List Binit, List Cinit, VectorXd lambda, List optsList, List optsList_pen)
{
    /*
     Input:
     Y is n*q matrix
     Z is n*(K*p) matrix
     A is p*r1 matrix
     B is K*r2 matrix
     C is q*r3 matrix
     S is r3*(r1*r2) matrix
     lambda is preset tuning parameters, a L-vector
     alpha is tuning parameter being in (0,1) to control elasnet
     gamma is another tuning parameter for MCP and SCAD
     penalty= 1(lasso),2(mcp), and 3(scad)
     dfmax is the preset maximum digrees of freedom
     threshold is the error to control convergence for outer iteration
     eps is the error to control convergence for inner iteration
     max_step is the max step to control convergence for outer iteration
     max_iter is the max step to control convergence for inner iteration
     is_setlam is logical, 1 for set lambda by data; 0 for given lambda
     setlam is a vector saving preseting lam_max,lam_min,alpha,nlam
     nlam is the number of tuning parameters
     
     Output:
     Dnew is a estimator of D(3) = C%*%S%*%kronecker(t(B),t(A))
     */
    opts.eps = as<double>(optsList["eps"]);
    opts.max_step = as<int>(optsList["max_step"]);
    opts.eps1 = as<double>(optsList["eps1"]);
    opts.max_step1 = as<int>(optsList["max_step1"]);
    opts.r1 = as<int>(optsList["r1"]);
    opts.r2 = as<int>(optsList["r2"]);
    opts.r3 = as<int>(optsList["r3"]);
    opts.q = as<int>(optsList["q"]);
    opts.G = as<int>(optsList["G"]);
    opts.degr = as<int>(optsList["degr"]);
    opts.K = as<int>(optsList["K"]);
    opts.n = Y.rows();
    opts.nx = as<int>(optsList["nx"]);
    
    opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);
    opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
    opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
    opts_pen.gamma = as<double>(optsList_pen["gamma"]);
    opts_pen.alpha = as<double>(optsList_pen["alpha"]);
    opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
    opts_pen.pen = as<int>(optsList_pen["pen"]);
    opts_pen.nlam = as<int>(optsList_pen["nlam"]);
    
    int l,i,j,pg, step, max_step=opts.max_step, nlam=opts_pen.nlam, ng = opts.G, nx = opts.nx;
    double  likhd0 = 2*pow(10, 6), lambda1, likhd1=likhd0, eps = opts.eps;
    
    MatrixXd S, A, B, C, Dnew, D0, Ynew = Y, Zj, Zj1, A1;
    VectorXi activeAj, activeA = VectorXi::Constant(nx, 0), convergence1 = VectorXi::Constant(4*ng, 1), breakid = VectorXi::Constant(ng, 0);
    MatrixXi betapath = MatrixXi::Constant(nx, nlam, 0), df = MatrixXi::Constant(ng, nlam, 0);
    VectorXd likhd = VectorXd::Constant(nlam, 0);
    
    List fit, S0(ng), A0(ng), B0(ng), C0(ng), Z(ng), Zm(ng), Am(ng);
    List Shat(nlam), Ahat(nlam), Bhat(nlam), Chat(nlam);
    for(i=0; i<ng; i++){
        S0[i] = as<MatrixXd>(Sinit[i]);
        A0[i] = as<MatrixXd>(Ainit[i]);
        B0[i] = as<MatrixXd>(Binit[i]);
        C0[i] = as<MatrixXd>(Cinit[i]);
        Z[i] = as<MatrixXd>(Z0[i]);
        Zm[i] = as<MatrixXd>(Z0[i]);
        Am[i] = as<MatrixXd>(Ainit[i]);
    }
    for (l = 0; l < nlam; l++) {
        lambda1 = lambda[l];
        step = 0;
        
        Ynew = Y;
        for (i=0; i<ng; i++) {
            Zj = Z[i];
            S = S0[i]; A = A0[i]; B = B0[i]; C = C0[i];
            D0 = C * S * kroneckerProduct(B, A).transpose();
            Ynew = Ynew - Zj*D0.transpose();
        }

        while(step < max_step){
            for(j=0;j<4*ng;j++) convergence1[j] = 1;
            for(j=0;j<ng;j++) breakid[j] = 0;
            step++;
            for(i=0; i<ng; i++){
                Zj = Z[i]; Zj1 = Zm[i]; A1 = Am[i];
                S = S0[i]; A = A0[i]; B = B0[i]; C = C0[i];
                D0 = C * S * kroneckerProduct(B, A).transpose();
                Ynew = Ynew + Zj*D0.transpose();
                fit = Estimation_pen_Dj(Ynew, Zj, S, A, B, C, convergence1, i, likhd1, lambda1, Zj1, A1, breakid);
                S = fit["S"]; A = fit["A"]; B = fit["B"]; C = fit["C"]; Zj1 = fit["Zm"]; A1 = fit["Am"];
                Dnew = C * S * kroneckerProduct(B, A).transpose();
                Ynew = Ynew - Zj * Dnew.transpose();
                S0[i] = S; A0[i] = A; B0[i] = B; C0[i] = C; Zm[i] = Zj1; Am[i] = A1;
            }
            if (breakid.sum()==ng) break;
            if (likhd1<likhd0) {
                if((likhd0-likhd1)/likhd0 < eps) break;
                else likhd0 = likhd1;
            }else{
                if(convergence1.sum()==0) break;
            }
        } //end while
        
        activeA = VectorXi::Constant(nx, 0);
        for(i=0;i<ng;i++){
            A = A0[i]; Zj = Z0[i];
            pg = A.rows();
            activeAj = VectorXi::Constant(pg, 0);
            for(j=0;j<pg;j++) if(A.row(j).norm()){
                activeA[i*pg+j] = 1;
                activeAj[j] = 1;
            }
            df(i,l) = activeAj.sum();
        }
        betapath.col(l) = activeA;
        likhd[l] = likhd0;
        Shat[l] = S0;
        Ahat[l] = A0;
        Bhat[l] = B0;
        Chat[l] = C0;
        
    }// end for

    return List::create(Named("likhd") = likhd, Named("betapath") = betapath, Named("df") = df, Named("Snew") = Shat, Named("Anew") = Ahat, Named("Bnew") = Bhat, Named("Cnew") = Chat);
}

//***------------------------------------------------------------------------------------**
//***----------- Estimation with penalizing functions in a whole column by CV------------**
// [[Rcpp::export]]
List EstPenColumnCV(MatrixXd Y, List Z0, MatrixXd Ytest, List Ztest0, List Sinit, List Ainit, List Binit, List Cinit, VectorXd lambda, List optsList, List optsList_pen)
{
    /*
    Input:
    Y is n*q matrix
    Z is n*(K*p) matrix
    A is p*r1 matrix
    B is K*r2 matrix
    C is q*r3 matrix
    S is r3*(r1*r2) matrix
    lambda is preset tuning parameters, a L-vector
    alpha is tuning parameter being in (0,1) to control elasnet
    gamma is another tuning parameter for MCP and SCAD
    penalty= 1(lasso),2(mcp), and 3(scad)
    dfmax is the preset maximum digrees of freedom
    threshold is the error to control convergence for outer iteration
    eps is the error to control convergence for inner iteration
    max_step is the max step to control convergence for outer iteration
    max_iter is the max step to control convergence for inner iteration
    is_setlam is logical, 1 for set lambda by data; 0 for given lambda
    setlam is a vector saving preseting lam_max,lam_min,alpha,nlam
    nlam is the number of tuning parameters

    Output:
    Dnew is a estimator of D(3) = C%*%S%*%kronecker(t(B),t(A))
    */
    opts.eps = as<double>(optsList["eps"]);
    opts.max_step = as<double>(optsList["max_step"]);
    opts.eps1 = as<double>(optsList["eps1"]);
    opts.max_step1 = as<double>(optsList["max_step1"]);
    opts.r1 = as<int>(optsList["r1"]);
    opts.r2 = as<int>(optsList["r2"]);
    opts.r3 = as<int>(optsList["r3"]);
    opts.q = as<int>(optsList["q"]);
    opts.G = as<int>(optsList["G"]);
    opts.degr = as<int>(optsList["degr"]);
    opts.K = as<int>(optsList["K"]);
    opts.nx = as<int>(optsList["nx"]);
    opts.n = Y.rows();

    opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);
    opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
    opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
    opts_pen.gamma = as<double>(optsList_pen["gamma"]);
    opts_pen.alpha = as<double>(optsList_pen["alpha"]);
    opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
    opts_pen.pen = as<int>(optsList_pen["pen"]);
    opts_pen.nlam = as<int>(optsList_pen["nlam"]);


    int l,i,j, pg, step, max_step=opts.max_step, nlam=opts_pen.nlam, ng = opts.G, K = opts.K, nx = opts.nx;
    double  likhd0 = 2*pow(10, 6), lambda1, likhd1=likhd0, eps = opts.eps;

    MatrixXd S, A, B, C, Dnew, D0, Ynew = Y, Ztest, Zj, Z2, A2, Zj1, A1;
    VectorXi activeAj, activeA = VectorXi::Constant(nx, 0), convergence1 = VectorXi::Constant(4*ng, 1), breakid = VectorXi::Constant(ng, 0);
    MatrixXi betapath = MatrixXi::Constant(nx, nlam, 0), df = MatrixXi::Constant(ng, nlam, 0);
    VectorXd likhd = VectorXd::Constant(nlam, 0);

    List fit, S0(ng), A0(ng), B0(ng), C0(ng), Z(ng), Zm(ng), Am(ng);
    for(i=0; i<ng; i++){
        S0[i] = as<MatrixXd>(Sinit[i]);
        A0[i] = as<MatrixXd>(Ainit[i]);
        B0[i] = as<MatrixXd>(Binit[i]);
        C0[i] = as<MatrixXd>(Cinit[i]);
        Z[i] = as<MatrixXd>(Z0[i]);
        Zm[i] = as<MatrixXd>(Z0[i]);
        Am[i] = as<MatrixXd>(Ainit[i]);
    }
    
    for (l = 0; l < nlam; l++) {
        lambda1 = lambda[l];
        step = 0;
        
        Ynew = Y;
        for (i=0; i<ng; i++) {
            Zj = Z[i];
            S = S0[i]; A = A0[i]; B = B0[i]; C = C0[i];
            D0 = C * S * kroneckerProduct(B, A).transpose();
            Ynew = Ynew - Zj*D0.transpose();
        }
        
        while(step < max_step){
            for(j=0;j<4*ng;j++) convergence1[j] = 1;
            for(j=0;j<ng;j++) breakid[j] = 0;
            step = step + 1;
            for(i=0; i<ng; i++){
                Zj = Z[i]; Zj1 = Zm[i]; A1 = Am[i];
                S = S0[i]; A = A0[i]; B = B0[i]; C = C0[i];
                D0 = C * S * kroneckerProduct(B, A).transpose();
                Ynew = Ynew + Zj*D0.transpose();
                fit = Estimation_pen_Dj(Ynew, Zj, S, A, B, C, convergence1, i, likhd1, lambda1, Zj1, A1, breakid);
                S = fit["S"]; A = fit["A"]; B = fit["B"]; C = fit["C"]; Zj1 = fit["Zm"]; A1 = fit["Am"];
                Dnew = C * S * kroneckerProduct(B, A).transpose();
                Ynew = Ynew - Zj * Dnew.transpose();
                S0[i] = S; A0[i] = A; B0[i] = B; C0[i] = C; Zm[i] = Zj1; Am[i] = A1;
            }
            if (breakid.sum()==ng) break;
            if(fabs(likhd0-likhd1)/likhd0 < eps) break;
            if(convergence1.sum()==0) break;
        } //end while

        Ynew = Ytest;
        activeA = VectorXi::Constant(nx, 0);
        for(i=0;i<ng;i++){
            A = A0[i]; Ztest = Ztest0[i];
            pg = A.rows();
            activeAj = VectorXi::Constant(pg, 0);
            for(j=0;j<pg;j++) if(A.row(j).norm()){
                activeA[i*pg+j] = 1;
                activeAj[j] = 1;
            }
            df(i,l) = activeAj.sum();
            Z2 = extractColsZ(Ztest,pg,K,activeAj);
            A2 = extractRows(A, activeAj);
            C = C0[i]; S = S0[i]; B = B0[i];
            Dnew = C * S * kroneckerProduct(B, A2).transpose();
            Ynew = Ynew - Z2 * Dnew.transpose();
        }
        betapath.col(l) = activeA;
        likhd[l] = Ynew.squaredNorm();

    }// end for
    return List::create(Named("likhd") = likhd, Named("df") = df, Named("betapath")=betapath);
}



//***------------------------------------------------------------------------------------**
//***-----------------------------------updateS in T4------------------------------------**
MatrixXd updateT4S(MatrixXd Y, MatrixXd Z, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
    int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, r4 = opts.r4;
    int d = r1*r2*r3*r4, k,k1,j,j1;
    MatrixXd ztilde = Z * kroneckerProduct(C,kroneckerProduct(B, A)), vectorS;
    VectorXd U;
    U.setZero(d);
    MatrixXd  V = MatrixXd::Constant(d, d, 0);
    
    for (k = 0; k < r1*r2*r3; k++) {
        for (j = 0; j < r4; j++) {
            U[k*r4 + j] = ztilde.col(k).transpose()*Y*D.col(j);
            for (k1 = 0; k1 < r1*r2*r3; k1++) {
                for (j1 = 0; j1 < r4; j1++) {
                    V(k*r4 + j, k1*r4 + j1) = kroneckerProduct(
                                                               ztilde.col(k1).array()*ztilde.col(k).array(),
                                                               (D.col(j1).array()*D.col(j).array()).transpose()).sum();
                }
            }
        }
    }
    vectorS = V.colPivHouseholderQr().solve(U);
    vectorS.resize(r4, r1*r2*r3);
    return vectorS;
}
//***------------------------------------------------------------------------------------**
//***---------------------------------updateD in T4--------------------------------------**
MatrixXd updateT4D(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
    int q = opts.q,j,kp;
    MatrixXd ztilde = Z * kroneckerProduct(C,kroneckerProduct(B, A));
    MatrixXd StZ = ztilde * S.transpose();
    MatrixXd Cnew = MatrixXd::Constant(q, opts.r4, 0);
    HouseholderQR<MatrixXd> qr;
    qr.compute(StZ);
    MatrixXd R = qr.matrixQR().triangularView<Upper>();
    MatrixXd Q = qr.householderQ();
    
    kp = StZ.cols();
    MatrixXd temp, IDEN = MatrixXd::Identity(kp, kp);
    if (pow(condition_numberQRSym(R),2) > 1e10 || isnan(UpTriangularInv(R).sum())){
        temp = tRbyR(R).block(0, 0, kp, kp) + (IDEN.array()*1e-4).matrix();
        for (j = 0; j < q; j++) {
            Cnew.row(j) = (temp.colPivHouseholderQr().solve(StZ.transpose()*Y.col(j))).transpose();
        }
    }else{
        for (j = 0; j < q; j++) {
            Cnew.row(j) = (QbyR(Q.transpose(), UpTriangularInv(R), 0)*Y.col(j)).transpose();
        }
    }
    return Cnew;
}
//***------------------------------------------------------------------------------------**
//***----------------------------------updateA in T4-------------------------------------**
MatrixXd updateT4A(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
    int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, K = opts.K, G = opts.G, p = A.rows();
    int d = r1 * p, t1,t2,t3,t4,j;
    MatrixXd W = D * S, Wt1, Zt2, Wt3, Zt4, zbw1, zbw2, vectorA, Gamma = kroneckerProduct(C,B);
    VectorXd tU;
    tU.setZero(d);
    MatrixXd tV = MatrixXd::Constant(d, d, 0);
    
    for (t2 = 0; t2<p; t2++) {
        Zt2 = Z.col(t2);
        for (j = 1; j<K*G; j++)    Zt2 = cbind_rcpp(Zt2, Z.col(j*p + t2));
        for (t1 = 0; t1<r1; t1++) {
            Wt1 = W.col(t1);
            for (j = 1; j<r2*r3; j++) Wt1 = cbind_rcpp(Wt1, W.col(j*r1 + t1));
            zbw1 = Zt2 * Gamma *(Wt1.transpose());
            tU[t2*r1 + t1] = (Y.array()*zbw1.array()).sum();
            for (t4 = 0; t4<p; t4++) {
                Zt4 = Z.col(t4);
                for (j = 1; j<K*G; j++) Zt4 = cbind_rcpp(Zt4, Z.col(j*p + t4));
                for (t3 = 0; t3<r1; t3++) {
                    Wt3 = W.col(t3);
                    for (j = 1; j<r2*r3; j++)    Wt3 = cbind_rcpp(Wt3, W.col(j*r1 + t3));
                    zbw2 = Zt4 * Gamma *(Wt3.transpose());
                    tV(t2*r1 + t1, t4*r1 + t3) = (zbw1.array()*zbw2.array()).sum();
                }
            }
        }
    }
    vectorA = tV.colPivHouseholderQr().solve(tU);
    vectorA.resize(r1, p);
    return vectorA.transpose();
}
//***------------------------------------------------------------------------------------**
//***----------------------------------updateB in T4-------------------------------------**
MatrixXd updateT4B(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
    int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, K = opts.K, G = opts.G, q = opts.q, n = opts.n, p = A.rows();
    int d = r2 * K, t1, t2, t3, t4, j;
    MatrixXd W = D * S, Wt1, Zt2, Wt3, Zt4, zaw1, zaw2, vectorB, Gamma = kroneckerProduct(C,A);
    VectorXd tU;
    tU.setZero(d);
    MatrixXd tV = MatrixXd::Constant(d, d, 0);
    for (t2 = 0; t2<K; t2++) {
        Zt2 = Z.block(0, t2*p, n, p);
        for (j = 1; j<G; j++)    Zt2 = cbind_rcpp(Zt2, Z.block(0, j*p*K + t2*p, n, p));
        for (t1 = 0; t1<r2; t1++) {
            Wt1 = W.block(0, t1*r1, q, r1);
            for (j = 1; j<r3; j++) Wt1 = cbind_rcpp(Wt1, W.block(0, j*r1*r2 + t1*r1, q, r1));
            zaw1 = Zt2 * Gamma * (Wt1.transpose());
            tU[t2*r2 + t1] = (Y.array()*zaw1.array()).sum();
            for (t4 = 0; t4<K; t4++) {
                Zt4 = Z.block(0, t4*p, n, p);
                for (j = 1; j<G; j++)    Zt4 = cbind_rcpp(Zt4, Z.block(0, j*p*K + t4*p, n, p));
                for (t3 = 0; t3<r2; t3++) {
                    Wt3 = W.block(0, t3*r1, q, r1);
                    for (j = 1; j<r3; j++) Wt3 = cbind_rcpp(Wt3, W.block(0, j*r1*r2 + t3*r1, q, r1));
                    zaw2 = Zt4 * Gamma * (Wt3.transpose());
                    tV(t2*r2 + t1, t4*r2 + t3) = (zaw1.array()*zaw2.array()).sum();
                }
            }
        }
    }
    vectorB = tV.colPivHouseholderQr().solve(tU);
    vectorB.resize(r2, K);
    return vectorB.transpose();
}
//***------------------------------------------------------------------------------------**
//***-----------------------------------updateC in T4------------------------------------**
MatrixXd updateT4C(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D)
{
    int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, K = opts.K, G = opts.G, q = opts.q, n = opts.n, p = A.rows();
    int d = r3 * G, t1, t2, t3, t4;
    MatrixXd W = D * S, Wt1, Zt2, Wt3, Zt4, zaw1, zaw2, vectorC, Gamma = kroneckerProduct(B,A);
    VectorXd tU;
    tU.setZero(d);
    MatrixXd tV = MatrixXd::Constant(d, d, 0);
    for (t2 = 0; t2<G; t2++) {
        Zt2 = Z.block(0, t2*p*K, n, p*K);
        for (t1 = 0; t1<r3; t1++) {
            Wt1 = W.block(0, t1*r1*r2, q, r1*r2);
            zaw1 = Zt2 * Gamma *(Wt1.transpose());
            tU[t2*r3 + t1] = (Y.array()*zaw1.array()).sum();
            for (t4 = 0; t4<G; t4++) {
                Zt4 = Z.block(0, t4*p*K, n, p*K);
                for (t3 = 0; t3<r3; t3++) {
                    Wt3 = W.block(0, t3*r1*r2, q, r1*r2);
                    zaw2 = Zt4 * Gamma *(Wt3.transpose());
                    tV(t2*r3 + t1, t4*r3 + t3) = (zaw1.array()*zaw2.array()).sum();
                }
            }
        }
    }
    vectorC = tV.colPivHouseholderQr().solve(tU);
    vectorC.resize(r3, G);
    return vectorC.transpose();
}
//***------------------------------------------------------------------------------------**
//***---------------------------T4 Estimation without penalty----------------------------**
// [[Rcpp::export]]
List EstimationT4(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, List optsList)
{
    /*
     Input:
     Y is n*q matrix
     Z is n*(K*p*G) matrix
     A is p*r1 matrix
     B is K*r2 matrix
     C is G*r3 matrix
     D is q*r4 matrix
     S is r4*(r1*r2*r3) matrix
     lambda is preset tuning parameters, a L-vector
     alpha is tuning parameter being in (0,1) to control elasnet
     gamma is another tuning parameter for MCP and SCAD
     penalty= 1(lasso),2(mcp), and 3(scad)
     dfmax is the preset maximum digrees of freedom
     threshold is the error to control convergence for outer iteration
     eps is the error to control convergence for inner iteration
     max_step is the max step to control convergence for outer iteration
     max_iter is the max step to control convergence for inner iteration
     is_setlam is logical, 1 for set lambda by data; 0 for given lambda
     setlam is a vector saving preseting lam_max,lam_min,alpha,nlam
     nlam is the number of tuning parameters
     
     Output:
     Dnew is a estimator of D(4) = D%*%S%*%kronecker(t(C),kronecker(t(B),t(A)))
     */
    
    opts.eps = as<double>(optsList["eps"]);
    opts.max_step = as<int>(optsList["max_step"]);
    opts.eps1 = as<double>(optsList["eps1"]);
    opts.max_step1 = as<int>(optsList["max_step1"]);
    opts.r1 = as<int>(optsList["r1"]);
    opts.r2 = as<int>(optsList["r2"]);
    opts.r3 = as<int>(optsList["r3"]);
    opts.r4 = as<int>(optsList["r4"]);
    opts.p = as<int>(optsList["p"]);
    opts.q = as<int>(optsList["q"]);
    opts.G = as<int>(optsList["G"]);
    opts.degr = as<int>(optsList["degr"]);
    opts.K = as<int>(optsList["K"]);
    opts.n = Y.rows();
    
    int j, max_step=opts.max_step;
    double  likhd0 = 2*pow(10, 6), likhd1, eps = opts.eps;
    MatrixXd Snew, Anew, Bnew, Cnew, Dnew, Dn;
    VectorXi convergence1 = VectorXi::Constant(5, 1);
    
    int step = 0;
    while(step < max_step){
        for(j=0;j<5;j++) convergence1[j] = 1;
        step = step + 1;
        Snew = updateT4S(Y, Z, A, B, C, D);
        Dn = D * Snew*kroneckerProduct(C,kroneckerProduct(B, A)).transpose();
        likhd1 = (Y - Z * Dn.transpose()).squaredNorm();
        if (likhd1<likhd0) {
            S = Snew;
            likhd0 = likhd1;
        }
        else convergence1[0]=0;
        
        Dnew = updateT4D(Y, Z, S, A, B, C, D);
        Dn = Dnew * S*kroneckerProduct(C,kroneckerProduct(B, A)).transpose();
        likhd1 = (Y - Z * Dn.transpose()).squaredNorm();
        if (likhd1<likhd0) {
            D = Dnew;
            likhd0 = likhd1;
        }
        else convergence1[1]=0;
        
        Cnew = updateT4C(Y, Z, S, A, B, C, D);
        Dn = D * S* kroneckerProduct(Cnew,kroneckerProduct(B, A)).transpose();
        likhd1 = (Y - Z * Dn.transpose()).squaredNorm();
        if (likhd1<likhd0) {
            C = Cnew;
            likhd0 = likhd1;
        }
        else convergence1[2]=0;
        
        Bnew = updateT4B(Y, Z, S, A, B, C, D);
        Dn = D * S* kroneckerProduct(C,kroneckerProduct(Bnew, A)).transpose();
        likhd1 = (Y - Z * Dn.transpose()).squaredNorm();
        if (likhd1<likhd0) {
            B = Bnew;
            likhd0 = likhd1;
        }
        else convergence1[3]=0;
        
        Anew = updateT4A(Y, Z, S, A, B, C, D);
        Dn = D * S* kroneckerProduct(C,kroneckerProduct(B, Anew)).transpose();
        likhd1 = (Y - Z * Dn.transpose()).squaredNorm();
        if (likhd1<likhd0) {
            A = Anew;
            if(fabs(likhd0-likhd1)/fabs(likhd0+1) < eps) break;
            likhd0 = likhd1;
        }
        else convergence1[4]=0;
        if(convergence1.sum()==0) break;
    }
    
    return List::create(Named("likhd") = likhd1, Named("Snew") = S, Named("Anew") = A, Named("Bnew") = B, Named("Cnew") = C, Named("Dnew") = D);
}
//***------------------------------------------------------------------------------------**
//***------------------------------setup tuning parameters-------------------------------**
// [[Rcpp::export]]
VectorXd setuplambdaT4(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, int nlam, VectorXd setlam)
{
    int n = Y.rows(), q = Y.cols(), G = C.rows(), p = A.rows()/G, K = B.rows(), j, jj, KG=K*G;
    int r1 = A.cols(), r2 = B.cols(), r3 = C.cols(), r4 = D.cols();
    double lam_max=setlam[0], lam_min=setlam[1], alpha=setlam[2];
    VectorXd lambda, lambda1, tmp, tmp1;
    MatrixXd S1, dcbs, V, V_1, Y1 = Y, Gammaj, Gamma_sqrt, svdu, svdd, U, L;
    Y1.resize(n*q, 1);
    VectorXi dims = VectorXi::Constant(4, 0);
    dims[0] = r1; dims[1] = r2; dims[2] = r3; dims[3] = r4;
    S1 = TransferModalUnfoldingsT(S, 4, 1, dims);
    dcbs = kroneckerProduct(D,kroneckerProduct(C, B))*(S1.transpose());
    VectorXi id = SEQ(1, p*KG, p);
    tmp = VectorXd::Constant(p, 0);
    for (j = 0; j < p; j++)
    {
        V = submatrix_col(Z, id.array() + j)*(dcbs.block(0, 0, KG, r1));
        for (jj = 1; jj < q; jj++) {
            V_1 = submatrix_col(Z, id.array() + j)*(dcbs.block(jj*KG, 0, KG, r1));
            V = rbind_rcpp(V, V_1);
        }
        Gammaj = ((V.transpose()*V).array() / n).matrix();
        
        JacobiSVD<MatrixXd> svd(Gammaj, ComputeThinU | ComputeThinV);
        svdu = svd.matrixU();
        svdd = (svd.singularValues()).asDiagonal();
        Gamma_sqrt = svdu * ((1 / (svdd.diagonal().array().sqrt())).matrix().asDiagonal())*(svdu.transpose());
        tmp1 = Y1.transpose()* V * Gamma_sqrt;
        tmp[j]=tmp1.array().abs().sum();
    }
    
    double max_tmp;
    max_tmp = (tmp.array()).maxCoeff()/(sqrt(n)*q*G);
    double max_lam;
    max_lam = lam_max * max_tmp / alpha;
    if (lam_min == 0) {
        lambda1.setLinSpaced(nlam - 1, log(max_lam), log(0.0001*max_lam));
        lambda.setLinSpaced(nlam, 0, 0);
        lambda.segment(0, nlam - 1) = lambda1.array().exp().matrix();
    }
    else {
        lambda1.setLinSpaced(nlam, log(max_lam), log(lam_min*max_lam));
        lambda = lambda1.array().exp();
    }
    return lambda;
}
//***------------------------------------------------------------------------------------**
//***--------------------update the jth row of matrix A with penalty---------------------**
MatrixXd updateT4A_penalty(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, MatrixXi beta, VectorXi &activeA, double lambda1)
{
    /*
     Input:
     Y is n*q matrix
     Z is n*(K*p*G) matrix
     A is p*r1 matrix
     B is K*r2 matrix
     C is G*r3 matrix
     D is q*r4 matrix
     S is r4*(r1*r2*r3) matrix
     penalty= 1(lasso),2(mcp), and 3(scad)
     lambda is preset tuning parameters, a L-vector
     alpha is tuning parameter being in (0,1) to control elasnet
     eps is the error to control convergence for outer iteration
     eps1 is the error to control convergence for inner iteration
     max_step is the max step to control convergence for outer iteration
     max_step1 is the max step to control convergence for inner iteration
     nlam is the number of tuning parameters
     */
    int r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, r4 = opts.r4, K = opts.K, G = opts.G, q = opts.q, n = opts.n, p = A.rows();
    int nlam = opts_pen.nlam, gamma = opts_pen.gamma, pen = opts_pen.pen, dfmax = opts_pen.dfmax;
    double alpha = opts_pen.alpha, eps1 = opts.eps1;
    int max_step1 = opts.max_step1, j,jj, active, KG=K*G;
    
    MatrixXd S1,aj,ajnew,zj;
    VectorXi dims = VectorXi::Constant(4, 0);
    dims[0] = r1; dims[1] = r2; dims[2] = r3; dims[3] = r4;
    S1 = TransferModalUnfoldingsT(S, 4, 1, dims);
    VectorXi id = SEQ(1, p*KG, p);
    MatrixXd Vnew = MatrixXd::Constant(n*q, r1*p, 0);
    MatrixXd Gamma_sqrtn = MatrixXd::Constant(r1, r1*p, 0);
    MatrixXd dcbs = kroneckerProduct(D,kroneckerProduct(C, B))*(S1.transpose());
    MatrixXd A1 = A, V, Gammaj, Gamma_sqrt, V_1, D4,L;
    D4 = D*S*kroneckerProduct(C, kroneckerProduct(B, A)).transpose();
    MatrixXd IDEN = MatrixXd::Identity(r1, r1);
    int count=0;
    
    for (j = 0; j < p; j++) {
        V = submatrix_col(Z, id.array() + j)*(dcbs.block(0, 0, KG, r1));
        for (jj = 1; jj < q; jj++) {
            V_1 = submatrix_col(Z, id.array() + j)*(dcbs.block(jj*KG, 0, KG, r1));
            V = rbind_rcpp(V, V_1);
        }
        Gammaj = ((V.transpose()*V).array() / n).matrix();
        L = Gammaj.llt().matrixL();
        if (isnan(UpTriangularInv(L.transpose()).sum())) {
            count++;
            Gammaj = ((V.transpose()*V).array() / n).matrix() + (IDEN.array()*1e-4).matrix();
            L = Gammaj.llt().matrixL();
        }
        Gamma_sqrtn.block(0, j*r1, r1, r1) = UpTriangularInv(L.transpose());
        Vnew.block(0, j*r1, n*q, r1) = QbyR(V, UpTriangularInv(L.transpose()),1);
        A1.row(j) = (QbyR(A.row(j).transpose(), L.transpose(), 0)).transpose();
    }
    
    MatrixXd Anew = A1;
    MatrixXd r = Y - Z * (D4.transpose());
    r.resize(n*q, 1);
    VectorXd ajnorm_old, ajnorm;
    ajnorm_old = ajnorm = VectorXd::Constant(p, 0);
    int converged1, step = 0;
    while (step<max_step1)
    {
        step++;
        active = 0;
        for (j = 0; j < p; j++)
            if (ajnorm[j] != 0) active = active + 1;
        if (active>dfmax) {
            beta = MatrixXi::Constant(p*r1, nlam, -9);
            return A;
        }
        for (j = 0; j < p;j++) {
            aj = Anew.row(j).transpose();
            zj = Vnew.block(0, j*r1, n*q, r1).transpose()*r/n + aj;
            if (zj.norm()==0) {
                ajnew = aj;
            }else{
                ajnew = updateAj(zj, n, r1, lambda1, alpha, gamma, pen);
            }
            if (isnan(ajnew.norm())) {
                stop("error: ajnew is nan!");
            }
            r = r - Vnew.block(0, j*r1, n*q, r1)*(ajnew-aj);
            Anew.row(j) = ajnew.transpose();
            ajnorm[j] = ajnew.norm();
        }
        converged1 = 1;
        for (j = 0; j < p;j++) {
            if (ajnorm[j] != 0 && ajnorm_old[j] != 0) {
                if ((A1.row(j) - Anew.row(j)).norm() / ajnorm_old[j]>eps1) {
                    converged1 = 0; break;
                }
            }
            else if (ajnorm[j] == 0 && ajnorm_old[j] != 0) {
                converged1 = 0; break;
            }
            else if (ajnorm[j] != 0 && ajnorm_old[j] == 0) {
                converged1 = 0; break;
            }
        }
        if (converged1) break;
        A1 = Anew;
        ajnorm_old = ajnorm;
    }//end while
    for (j = 0; j<p; j++) {
        Anew.row(j) = (QbyR(Anew.row(j).transpose(),Gamma_sqrtn.block(0, j*r1, r1, r1),0)).transpose();
        if (ajnorm[j]) activeA[j] = 1;
    }
    return Anew;
}

//***------------------------------------------------------------------------------------**
//***-----Old main function: Estimation with penalizing functions in a whole column -----**
// [[Rcpp::export]]
List EstPenColumnT4(MatrixXd Y, MatrixXd Z, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D,
                    VectorXd lambda, List optsList, List optsList_pen)
{
    /*
     Input:
     Y is n*q matrix
     Z is n*(K*p*G) matrix
     A is p*r1 matrix
     B is K*r2 matrix
     C is G*r3 matrix
     D is q*r4 matrix
     S is r4*(r1*r2*r3) matrix
     lambda is preset tuning parameters, a L-vector
     alpha is tuning parameter being in (0,1) to control elasnet
     gamma is another tuning parameter for MCP and SCAD
     penalty= 1(lasso),2(mcp), and 3(scad)
     dfmax is the preset maximum digrees of freedom
     eps is the error to control convergence for outer iteration
     eps1 is the error to control convergence for inner iteration
     max_step is the max step to control convergence for outer iteration
     max_step1 is the max step to control convergence for inner iteration
     nlam is the number of tuning parameters
     
     Output:
     Dnew is a estimator of D(4) = D%*%S%*%kronecker(t(C),kronecker(t(B),t(A)))
     */
    opts.eps = as<double>(optsList["eps"]);
    opts.max_step = as<double>(optsList["max_step"]);
    opts.eps1 = as<double>(optsList["eps1"]);
    opts.max_step1 = as<double>(optsList["max_step1"]);
    opts.r1 = as<int>(optsList["r1"]);
    opts.r2 = as<int>(optsList["r2"]);
    opts.r3 = as<int>(optsList["r3"]);
    opts.r4 = as<int>(optsList["r4"]);
    opts.p = as<int>(optsList["p"]);
    opts.q = as<int>(optsList["q"]);
    opts.G = as<int>(optsList["G"]);
    opts.degr = as<int>(optsList["degr"]);
    opts.K = as<int>(optsList["K"]);
    opts.n = Y.rows();
    
    opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);
    opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
    opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
    opts_pen.gamma = as<double>(optsList_pen["gamma"]);
    opts_pen.alpha = as<double>(optsList_pen["alpha"]);
    opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
    opts_pen.pen = as<int>(optsList_pen["pen"]);
    opts_pen.nlam = lambda.size();
    
    
    int l,j, step, max_step=opts.max_step, nlam=opts_pen.nlam, G = opts.G, K = opts.K, KG = K*G, p = opts.p;
    int q = opts.q, r1 = opts.r1, r2 = opts.r2, r3 = opts.r3, r4 = opts.r4;
    double  likhd0 = 2*pow(10, 6), lambda1, likhd1, eps = opts.eps;
    MatrixXd Dnew, Anew, Bnew, Snew, Cnew, Dn, Z1, A1, Z2, A2;
    VectorXi activeA = VectorXi::Constant(p, 0), convergence1 = VectorXi::Constant(5, 1);
    VectorXi df = VectorXi::Constant(nlam, 0);
    MatrixXi betapath = MatrixXi::Constant(p, nlam, 0);
    VectorXd likhd = VectorXd::Constant(nlam, 0);
    MatrixXd Apath, Bpath, Cpath, Dpath, Spath, temp;
    Apath = MatrixXd::Constant(p*r1, nlam, 0);
    Bpath = MatrixXd::Constant(K*r2, nlam, 0);
    Cpath = MatrixXd::Constant(G*r3, nlam, 0);
    Dpath = MatrixXd::Constant(q*r4, nlam, 0);
    Spath = MatrixXd::Constant(r1*r2*r3*r4, nlam, 0);
    Anew = A;
    A1 = A;
    Z1 = Z;

    for (l = 0; l < nlam; l++) {
        lambda1 = lambda[l];
        step = 0;
        while (step<max_step) {
            for(j=0;j<5;j++) convergence1[j] = 1;
            step ++;
            Snew = updateT4S(Y, Z1, A1, B, C, D);
            Dn = D * Snew * kroneckerProduct(C,kroneckerProduct(B, A1)).transpose();
            likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            if (likhd1<likhd0) {
                S = Snew;
                likhd0 = likhd1;
            }
            else convergence1[0]=0;
            
            Dnew = updateT4D(Y, Z1, S, A1, B, C, D);
            Dn = Dnew * S*kroneckerProduct(C,kroneckerProduct(B, A1)).transpose();
            likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            if (likhd1<likhd0) {
                D = Dnew;
                likhd0 = likhd1;
            }
            else convergence1[1]=0;
            
            Cnew = updateT4C(Y, Z1, S, A1, B, C, D);
            Dn = D * S * kroneckerProduct(Cnew,kroneckerProduct(B, A1)).transpose();
            likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            if (likhd1<likhd0) {
                C = Cnew;
                likhd0 = likhd1;
            }
            else convergence1[2]=0;
            
            Bnew = updateT4B(Y, Z1, S, A1, B, C, D);
            Dn = D * S * kroneckerProduct(C,kroneckerProduct(Bnew, A1)).transpose();
            likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            if (likhd1<likhd0) {
                B = Bnew;
                likhd0 = likhd1;
            }
            else convergence1[3]=0;
            
            for(j=0;j<p;j++) activeA[j] = 0;
            Anew = updateT4A_penalty(Y, Z, S, A, B, C, D, betapath, activeA, lambda1);
            if(activeA.sum()<r1){
                A1 = A;
                Z1 = Z;
                break;
            }
            else{
                Z2 = extractColsZ(Z,p,KG,activeA);
                A2 = extractRows(Anew, activeA);
                Dn = D * S * kroneckerProduct(C,kroneckerProduct(B, A2)).transpose();
                likhd1 = (Y - Z2 * Dn.transpose()).squaredNorm();
                if (likhd1<likhd0) {
                    A = Anew;
                    Z1 = Z2;
                    A1 = A2;
                    if ((likhd0 - likhd1) / likhd0<eps) break;
                    else  likhd0 = likhd1;
                }
                else convergence1[4]=0;
            }
            if(convergence1.sum()==0) break;
        } //end while
        
        for(j=0;j<p;j++) activeA[j] = 0;
        for(j=0;j<p;j++) if(A.row(j).norm()) activeA[j] = 1;
        df[l] = activeA.sum();
        likhd[l] = likhd0;
        betapath.col(l) = activeA;
        
        temp = A; temp.resize(p*r1, 1);        Apath.col(l) = temp;
        temp = B; temp.resize(K*r2, 1);        Bpath.col(l) = temp;
        temp = C; temp.resize(G*r3, 1);        Cpath.col(l) = temp;
        temp = D; temp.resize(q*r4, 1);        Dpath.col(l) = temp;
        temp = S; temp.resize(r1*r2*r3*r4, 1);        Spath.col(l) = temp;
    }// end for
    
    return List::create(Named("likhd") = likhd, Named("betapath") = betapath, Named("df") = df, Named("lambda")=lambda, Named("Spath") = Spath, Named("Apath") = Apath, Named("Bpath") = Bpath, Named("Cpath") = Cpath, Named("Dpath") = Dpath);
}

//***------------------------------------------------------------------------------------**
//***---------- Estimation with penalizing functions in a whole column by CV-------------**
// [[Rcpp::export]]
List EstPenColumnT4CV(MatrixXd Y, MatrixXd Z, MatrixXd Ytest, MatrixXd Ztest, MatrixXd S, MatrixXd A, MatrixXd B, MatrixXd C, MatrixXd D, VectorXd lambda, List optsList, List optsList_pen)
{
     /*
     Input:
     Y is n*q matrix
     Z is n*(K*p*G) matrix
     A is p*r1 matrix
     B is K*r2 matrix
     C is G*r3 matrix
     D is q*r4 matrix
     S is r4*(r1*r2*r3) matrix
     lambda is preset tuning parameters, a L-vector
     alpha is tuning parameter being in (0,1) to control elasnet
     gamma is another tuning parameter for MCP and SCAD
     penalty= 1(lasso),2(mcp), and 3(scad)
     dfmax is the preset maximum digrees of freedom
     eps is the error to control convergence for outer iteration
     eps1 is the error to control convergence for inner iteration
     max_step is the max step to control convergence for outer iteration
     max_step1 is the max step to control convergence for inner iteration
     nlam is the number of tuning parameters
     
     Output:
     Dnew is a estimator of D(4) = D%*%S%*%kronecker(t(C),kronecker(t(B),t(A)))
     */
    opts.eps = as<double>(optsList["eps"]);
    opts.max_step = as<double>(optsList["max_step"]);
    opts.eps1 = as<double>(optsList["eps1"]);
    opts.max_step1 = as<double>(optsList["max_step1"]);
    opts.r1 = as<int>(optsList["r1"]);
    opts.r2 = as<int>(optsList["r2"]);
    opts.r3 = as<int>(optsList["r3"]);
    opts.r4 = as<int>(optsList["r4"]);
    opts.p = as<int>(optsList["p"]);
    opts.q = as<int>(optsList["q"]);
    opts.G = as<int>(optsList["G"]);
    opts.degr = as<int>(optsList["degr"]);
    opts.K = as<int>(optsList["K"]);
    opts.n = Y.rows();
    
    opts_pen.isPenColumn = as<int>(optsList_pen["isPenColumn"]);
    opts_pen.lam_max = as<double>(optsList_pen["lam_max"]);
    opts_pen.lam_min = as<double>(optsList_pen["lam_min"]);
    opts_pen.gamma = as<double>(optsList_pen["gamma"]);
    opts_pen.alpha = as<double>(optsList_pen["alpha"]);
    opts_pen.dfmax = as<int>(optsList_pen["dfmax"]);
    opts_pen.pen = as<int>(optsList_pen["pen"]);
    opts_pen.nlam = lambda.size();
    
    int l,j, step, max_step=opts.max_step, nlam=opts_pen.nlam, KG = opts.K*opts.G, p = opts.p, r1 = opts.r1;
    double  likhd0 = pow(10, 6), lambda1, likhd1, eps = opts.eps;
    MatrixXd Dnew, Anew, Bnew, Snew, Cnew, Dn, Z1, A1, Z2, A2;
    VectorXi activeA = VectorXi::Constant(p, 0), convergence1 = VectorXi::Constant(5, 1);
    VectorXi df = VectorXi::Constant(nlam, 0);
    MatrixXi betapath = MatrixXi::Constant(p, nlam, 0);
    VectorXd likhd = VectorXd::Constant(nlam, 0);
    Anew = A;
    A1 = A;
    Z1 = Z;
    
    for (l = 0; l < nlam; l++) {
        lambda1 = lambda[l];
        step = 0;
        while (step<max_step) {
            for(j=0;j<5;j++) convergence1[j] = 1;
            step ++;
            Snew = updateT4S(Y, Z1, A1, B, C, D);
            Dn = D * Snew * kroneckerProduct(C,kroneckerProduct(B, A1)).transpose();
            likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            if (likhd1<likhd0) {
                S = Snew;
                likhd0 = likhd1;
            }
            else convergence1[0]=0;
            
            Dnew = updateT4D(Y, Z1, S, A1, B, C, D);
            Dn = Dnew * S*kroneckerProduct(C,kroneckerProduct(B, A1)).transpose();
            likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            if (likhd1<likhd0) {
                D = Dnew;
                likhd0 = likhd1;
            }
            else convergence1[1]=0;
            
            Cnew = updateT4C(Y, Z1, S, A1, B, C, D);
            Dn = D * S * kroneckerProduct(Cnew,kroneckerProduct(B, A1)).transpose();
            likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            if (likhd1<likhd0) {
                C = Cnew;
                likhd0 = likhd1;
            }
            else convergence1[2]=0;
            
            Bnew = updateT4B(Y, Z1, S, A1, B, C, D);
            Dn = D * S * kroneckerProduct(C,kroneckerProduct(Bnew, A1)).transpose();
            likhd1 = (Y - Z1 * Dn.transpose()).squaredNorm();
            if (likhd1<likhd0) {
                B = Bnew;
                likhd0 = likhd1;
            }
            else convergence1[3]=0;
            
            for(j=0;j<p;j++) activeA[j] = 0;
            Anew = updateT4A_penalty(Y, Z, S, A, B, C, D, betapath, activeA, lambda1);
            if(activeA.sum()<r1){
                A1 = A;
                Z1 = Z;
                break;
            }
            else{
                Z2 = extractColsZ(Z,p,KG,activeA);
                A2 = extractRows(Anew, activeA);
                Dn = D * S * kroneckerProduct(C,kroneckerProduct(B, A2)).transpose();
                likhd1 = (Y - Z2 * Dn.transpose()).squaredNorm();
               if (likhd1<likhd0) {
                    A = Anew;
                    Z1 = Z2;
                    A1 = A2;
                    if ((likhd0 - likhd1) / likhd0<eps) break;
                    else  likhd0 = likhd1;
                }
                else convergence1[4]=0;
            }
            if(convergence1.sum()==0) break;
        } //end while
        
        for(j=0;j<p;j++) activeA[j] = 0;
        for(j=0;j<p;j++) if(A.row(j).norm()) activeA[j] = 1;
        df[l] = activeA.sum();
        Z2 = extractColsZ(Ztest,p,KG,activeA);
        A2 = extractRows(A, activeA);
        Dn = D * S * kroneckerProduct(C,kroneckerProduct(B, A2)).transpose();
        likhd[l] = (Ytest - Z2 * Dn.transpose()).squaredNorm();
        betapath.col(l) = activeA;
    }// end for
    
    return List::create(Named("likhd") = likhd, Named("df") = df, Named("betapath")=betapath);
}

