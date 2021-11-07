#include <Rcpp.h>
#include <RcppEigen.h>
#include <cmath>
#include <Eigen/Dense>
#include <iostream>
#include <Eigen/SVD>
#include <Eigen/Cholesky>
#include <Eigen/LU>
#include <list>
using namespace Rcpp;
using namespace Eigen;

//----------------------------------------------------------------**
//***----------------------parameters---------------**
struct Options
{
	int r1;
	int r2;
	int r3;
	int r4;
	int p;
	int q;
	int G;
	int n;
	int degr;
	int K;
	int nx;
	double eps;
	double eps1;
	int max_step;
    int max_step1;	
}opts;
//----------------------------------------------------------------**
//***----------------------parameters for penalization---------------**
struct Options_pen
{
	int pen; 
	int nlam;
	int dfmax;
	int isPenColumn;
	double lam_max;
	double lam_min;
	double alpha;
	double gamma;
}opts_pen;
//----------------------------------------------------------------**
//***----------------------------sequece--------------------------**
VectorXi SEQ(int min, int max, int by)
{
	int leng = static_cast<int>(floor((max - min) / by + pow(3, -12)) + 1);
	VectorXi C = VectorXi::Constant(leng, 0);
	for (int i = 0; i < leng; i++)
		C[i] = min + by * i;
	return C;
}
//----------------------------------------------------------------**
//***----------------------uppertriangleNorm1 of AA---------------**
double uppertriangleNorm1AA(MatrixXd A)
{
	int i, j, k, p, ii;
	double norm1 = 0, temp1, temp2;
	p = A.cols();

	for (k = 0; k < p; k++) {
		temp2 = 0;
		for (j = 0; j < p; j++) {
			temp1 = 0;
			ii = k < j ? k : j;
			for (i = 0; i <= ii; i++) temp1 += A(i, j) * A(i, k);
			temp2 += fabs(temp1);
		}
		if (temp2 > norm1) norm1 = temp2;
	}
	return norm1;
}
//----------------------------------------------------------------**
//***----------------------uppertriangleNorm1 of A----------------**
double uppertriangleNorm1(MatrixXd A)
{
	int j, k, p;
	double norm1 = 0, temp1;
	p = A.cols();

	for (k = 0; k < p; k++) {
		temp1 = 0;
		for (j = 0; j <= k; j++) {
			temp1 += fabs(A(j, k));
		}
		if (temp1 > norm1) norm1 = temp1;
	}
	return norm1;
}

//----------------------------------------------------------------**
//***----------------------UpTriangularInv------------------------**
MatrixXd UpTriangularInv(MatrixXd A)
{
	int i, j, k,n=A.cols();
	MatrixXd B = MatrixXd::Constant(n, n, 0);
	for (i = 0; i < n; i++) B(i, i) = 1;
	for (i = n - 1; i >= 0; i--)//rows
	{
		if (A(i, i) != 1)
			for (j = i; j<n; j++)
				B(i, j) = B(i, j) / A(i, i);
		if (i>0)
		{
			for (j = i; j<n; j++)// columns
				for (k = 0; k<i; k++)// rows
					B(k, j) = B(k, j) - A(k, i) * B(i, j);
		}
	}
	return B;
}
//----------------------------------------------------------------**
//***---------------------- Norm1 of a Matrix---------------------**
double MatrixNorm1(MatrixXd A)
{
  int j, k, p;
  double norm1 = 0, temp1;
  p = A.cols();
  
  for (k = 0; k < p; k++) {
    temp1 = 0;
    for (j = 0; j < p; j++) temp1 += fabs(A(j,k));
    if (temp1 > norm1) norm1 = temp1;
  }
  return norm1;
}
//----------------------------------------------------------------**
//***---------------------- Q*R of qr decomposition --------------**
MatrixXd QbyR(MatrixXd Q, MatrixXd R, int isQR)
{
	//isQR=1 denotes Q*R; otherwise R*Q
	int i, j, k, p = R.cols(),n;
	double temp1;
	MatrixXd A = Q;
	if(isQR){
		n = Q.rows();
		for (k = 0; k < p; k++)
			for (j = 0; j < n; j++) {
				temp1 = 0;
				for (i = 0; i <= k; i++) temp1 += Q(j, i)*R(i, k);
				A(j,k) = temp1;
			}
	}
	else{
		n = Q.cols();
		for (k = 0; k < n; k++)
			for (j = 0; j < p; j++) {
				temp1 = 0;
				for (i = j; i < p; i++) temp1 += R(j, i) * Q(i, k);
				A(j,k) = temp1;
			}
	}
	return A;
}
//----------------------------------------------------------------**
//***---------------------- R*R of Upper triangle ----------------**
MatrixXd tRbyR(MatrixXd R)
{
  int i, j, k, ii, p;
  double temp1;
  MatrixXd A = R;
  p = R.cols();
  
  for (k = 0; k < p; k++) {
    for (j = 0; j < p; j++) {
      temp1 = 0;
	  ii = k < j ? k : j;
      for (i = 0; i <= ii; i++) temp1 += R(i, j) * R(i, k);
      A(j,k) = temp1;
    }
  }
  return A;
}
//----------------------------------------------------------------**
//***----------------------cbind----------------------------------**
MatrixXd cbind_rcpp(MatrixXd A, MatrixXd B)
{
	int n = A.rows();
	int p1 = A.cols();
	int p2 = B.cols();
	MatrixXd C = MatrixXd::Constant(n, p1 + p2, 0);
	C.block(0, 0, n, p1) = A;
	C.block(0, p1, n, p2) = B;
	return C;
}
//----------------------------------------------------------------**
//***----------------------rbind----------------------------------**
MatrixXd rbind_rcpp(MatrixXd A, MatrixXd B)
{
	int n1 = A.rows();
	int n2 = B.rows();
	int p = A.cols();
	MatrixXd C = MatrixXd::Constant(n1 + n2, p, 0);
	C.block(0, 0, n1, p) = A;
	C.block(n1, 0, n2, p) = B;
	return C;
}

//----------------------------------------------------------------**
//***----------------------extract columns------------------------**
MatrixXd extractRows(MatrixXd A, VectorXi b)
{
	int p1 = A.rows(), n = A.cols(), p2=b.sum(), j, count = 0;
	if(p2>p1) stop("The length of index b must be not great than the number of columns of A!");
	MatrixXd C = MatrixXd::Constant(p2, n, 0);
	for(j=0; j< p1; j++)  if(b[j])   C.row(count++) = A.row(j);
	return C;
}
//----------------------------------------------------------------**
//***----------------------extract submatrix----------------------**
MatrixXd extractColsZ(MatrixXd Z, int p, int K, VectorXi b)
{
	int n = Z.rows(), p2=b.sum(), i, j, count;
	if(p2>p) stop("The length of index b must be not great than the number of columns of A!");
	MatrixXd C = MatrixXd::Constant(n, p2*K, 0);
	for(i=0;i<K;i++){
	  count = 0;
	  for(j=0; j< p; j++) if(b[j]) C.col(i*p2 + (count++)) = Z.col(i*p+j);
	}
	return C;
}
//----------------------------------------------------------------**
//***----------------------extract submatrix----------------------**
MatrixXd extractCols_by(MatrixXd Z, int n, int min, int K, int by)
{
	MatrixXd C = MatrixXd::Constant(n, K, 0);
	for(int i=0;i<K;i++) C.col(i) = Z.col(min+by*i);
	return C;
}
//----------------------------------------------------------------**
//***--------QRcondition_number for design matrix-----------------**
double condition_numberQR(MatrixXd R)
{
	return uppertriangleNorm1AA(R) * uppertriangleNorm1AA(UpTriangularInv(R));
}
//----------------------------------------------------------------**
//***------QRcondition_number for symetric matrix-----------------**
double condition_numberQRSym(MatrixXd R)
{
  return uppertriangleNorm1(R) * uppertriangleNorm1(UpTriangularInv(R));
}
//----------------------------------------------------------------**
//***----------------solve linear system by QR--------------------**
VectorXd solveEquationQR(MatrixXd A, VectorXd b)
{
  //solve linear system Ax = b, A = A.transpose();
  int p = b.size();
  HouseholderQR<MatrixXd> qr;
  qr.compute(A);
  MatrixXd R, Q;
  Q = qr.householderQ();
  R = qr.matrixQR().triangularView<Upper>();
  VectorXd X;
  
  if (condition_numberQRSym(R) > 1e10){
    MatrixXd temp, IDEN = MatrixXd::Identity(p, p);
    temp = A + (IDEN.array()*1e-4).matrix();
    X = temp.colPivHouseholderQr().solve(b);
  }
  else X = QbyR(Q.transpose(),UpTriangularInv(R),0)*b;
  return X;
}
//----------------------------------------------------------------**
//***--------------------transfer modal of unfoldings-------------**
//// [[Rcpp::export]]
MatrixXd TransferModalUnfoldings(MatrixXd S, int d1, int d2, int r1, int r2, int r3)
{
  //From S_(d1) to S_(d2)
  int j;
  MatrixXd S1,S_3;
  if (d1 == 3) {
    if (d2 == 1){
      S1 = S.row(0).transpose();
      S1.resize(r1, r2); 
      for (j = 1; j < r3;j++) {
        S_3 = S.row(j).transpose();
        S_3.resize(r1, r2);
        S1 = cbind_rcpp(S1, S_3);// S3 is r3 *(r2r1) matrix
      }
    }	
    if (d2 == 2) {
      S1 = S.block(0, 0, r3, r1).transpose();
      S1.resize(1,r1*r3);
      for (j = 1; j < r2; j++) {
        S_3 = S.block(0, j*r1, r3, r1).transpose();
        S_3.resize(1,r1*r3);
        S1 = rbind_rcpp(S1, S_3);
      }
    }
  }
  
  if (d1 == 2) {
    if (d2 == 1) {
      S1 = S.block(0, 0, r2, r1).transpose();
      for (j = 1; j < r3; j++) {
        S1 = cbind_rcpp(S1,S.block(0,j*r1,r2,r1).transpose());
      }
    }
    if (d2 == 3) {
      S1 = S.block(0, 0, r2, r1).transpose();
      S1.resize(1, r2*r1);
      for (j = 1; j < r3; j++) {
        S_3 = S.block(0, j*r1, r2, r1).transpose();
        S_3.resize(1, r2*r1);
        S1 = rbind_rcpp(S1, S_3);
      }
    }
    
  }
  
  if (d1 == 1) {
    if (d2 == 2) {
      S1 = S.block(0, 0, r1, r2).transpose();
      for (j = 1; j < r3; j++) {
        S1 = cbind_rcpp(S1, S.block(0, j*r2, r1, r2).transpose());
      }
    }
    if (d2 == 3) {
      S1 = S.block(0, 0, r1, r2);
      S1.resize(1, r1*r2);
      for (j = 1; j < r3; j++) {
        S_3 = S.block(0, j*r2, r1, r2);
        S_3.resize(1, r1*r2);
        S1 = rbind_rcpp(S1, S_3);
      }
    }
  }
  return S1;
}
//----------------------------------------------------------------**
//***--------------------transfer modal of unfoldings-------------**
//// [[Rcpp::export]]
MatrixXd TransferModalUnfoldingsT4(MatrixXd S, int d1, int d2, int r1, int r2, int r3, int r4)
{
    //From S_(d1) to S_(d2)
    int l, k;
    MatrixXd S1, S2, S_3;
    if (d1 == 1) {
        if (d2 == 2) {
            S1 = S.block(0, 0, r1, r2).transpose();
            for (l = 1; l < r3*r4; l++) {
                S1 = cbind_rcpp(S1, S.block(0, l*r2, r1, r2).transpose());
            }
        }
        if (d2 == 3) {
            S1 = S.block(0, 0, r1, r2);
            S1.resize(1, r1*r2);
            for (k = 1; k < r3; k++) {
                S_3 = S.block(0, k*r2, r1, r2);
                S_3.resize(1, r1*r2);
                S1 = rbind_rcpp(S1, S_3);
            }
            for (l = 1; l < r4; l++) {
                S2 = S.block(0, l*r2*r3, r1, r2);
                S2.resize(1, r1*r2);
                for (k = 1; k < r3; k++) {
                    S_3 = S.block(0, l*r2*r3+k*r2, r1, r2);
                    S_3.resize(1, r1*r2);
                    S2 = rbind_rcpp(S2, S_3);
                }
                S1 = cbind_rcpp(S1, S2);
            }
        }
        if (d2 == 4) {
            S1 = S.block(0, 0, r1, r2);
            S1.resize(1, r1*r2);
            for (l = 1; l < r4; l++) {
                S_3 = S.block(0, l*r2*r3, r1, r2);
                S_3.resize(1, r1*r2);
                S1 = rbind_rcpp(S1, S_3);
            }
            for (k = 1; k < r3; k++) {
                S2 = S.block(0, k*r2, r1, r2);
                S2.resize(1, r1*r2);
                for (l = 1; l < r4; l++) {
                    S_3 = S.block(0, l*r2*r3+k*r2, r1, r2);
                    S_3.resize(1, r1*r2);
                    S2 = rbind_rcpp(S2, S_3);
                }
                S1 = cbind_rcpp(S1, S2);
            }
        }
    }
    if (d1 == 2) {
        if (d2 == 1) {
            S1 = S.block(0, 0, r2, r1).transpose();
            for (l = 1; l < r3*r4; l++) {
                S1 = cbind_rcpp(S1, S.block(0,l*r1,r2,r1).transpose());
            }
        }
        if (d2 == 3) {
            S1 = S.block(0, 0, r2, r1).transpose();
            S1.resize(1, r2*r1);
            for (k = 1; k < r3; k++) {
                S_3 = S.block(0, k*r1, r2, r1).transpose();
                S_3.resize(1, r2*r1);
                S1 = rbind_rcpp(S1, S_3);
            }
            for (l = 1; l < r4; l++) {
                S2 = S.block(0, l*r1*r3, r2, r1).transpose();
                S2.resize(1, r2*r1);
                for (k = 1; k < r3; k++) {
                    S_3 = S.block(0, l*r1*r3+k*r1, r2, r1).transpose();
                    S_3.resize(1, r2*r1);
                    S2 = rbind_rcpp(S2, S_3);
                }
                S1 = cbind_rcpp(S1, S2);
            }
        }
        if (d2 == 4) {
            S1 = S.block(0, 0, r2, r1).transpose();
            S1.resize(1, r2*r1);
            for (l = 1; l < r4; l++) {
                S_3 = S.block(0, l*r1*r3, r2, r1).transpose();
                S_3.resize(1, r2*r1);
                S1 = rbind_rcpp(S1, S_3);
            }
            for (k = 1; k < r3; k++) {
                S2 = S.block(0, k*r1, r2, r1).transpose();
                S2.resize(1, r2*r1);
                for (l = 1; l < r4; l++) {
                    S_3 = S.block(0, l*r1*r3+k*r1, r2, r1).transpose();
                    S_3.resize(1, r2*r1);
                    S2 = rbind_rcpp(S2, S_3);
                }
                S1 = cbind_rcpp(S1, S2);
            }
        }
    }
    if (d1 == 3) {
        if (d2 == 1){
            S1 = S.block(0, 0, 1, r1*r2).transpose();
            S1.resize(r1, r2);
            for (k = 1; k < r3; k++) {
                S_3 = S.block(k, 0, 1, r1*r2).transpose();
                S_3.resize(r1, r2);
                S1 = cbind_rcpp(S1, S_3);
            }
            for (l = 1; l < r4; l++) {
                S2 = S.block(0, l*r1*r2, 1, r1*r2).transpose();
                S2.resize(r1, r2);
                for (k = 1; k < r3; k++) {
                    S_3 = S.block(k, l*r1*r2, 1, r1*r2).transpose();
                    S_3.resize(r1, r2);
                    S2 = cbind_rcpp(S2, S_3);
                }
                S1 = cbind_rcpp(S1, S2);
            }
        }
        if (d2 == 2) {
            S1 = S.block(0, 0, 1, r1*r2).transpose();
            S1.resize(r1, r2); S1 = S1.transpose();
            for (k = 1; k < r3; k++) {
                S_3 = S.block(k, 0, 1, r1*r2).transpose();
                S_3.resize(r1, r2); S_3 = S_3.transpose();
                S1 = cbind_rcpp(S1, S_3);
            }
            for (l = 1; l < r4; l++) {
                S2 = S.block(0, l*r1*r2, 1, r1*r2).transpose();
                S2.resize(r1, r2); S2 = S2.transpose();
                for (k = 1; k < r3; k++) {
                    S_3 = S.block(k, l*r1*r2, 1, r1*r2).transpose();
                    S_3.resize(r1, r2); S_3 = S_3.transpose();
                    S2 = cbind_rcpp(S2, S_3);
                }
                S1 = cbind_rcpp(S1, S2);
            }
        }
        if (d2 == 4) {
            S1 = S.row(0).transpose();
            S1.resize(r1*r2, r4); S1 = S1.transpose();
            for (k = 1; k < r3; k++) {
                S_3 = S.row(k).transpose();
                S_3.resize(r1*r2, r4); S_3 = S_3.transpose();
                S1 = cbind_rcpp(S1, S_3);
            }
        }
    }
    if (d1 == 4) {
        if (d2 == 1) {
            S1 = S.row(0).transpose();
            S1.resize(r1, r2*r3);
            for (l = 1; l < r4; l++) {
                S_3 = S.row(l).transpose();
                S_3.resize(r1, r2*r3);
                S1 = cbind_rcpp(S1, S_3);
            }
        }
        if (d2 == 2) {
            S1 = S.block(0, 0, 1, r1*r2).transpose();
            S1.resize(r1, r2); S1 = S1.transpose();
            for (k = 1; k < r3; k++) {
                S_3 = S.block(0, k*r1*r2, 1, r1*r2).transpose();
                S_3.resize(r1, r2); S_3 = S_3.transpose();
                S1 = cbind_rcpp(S1, S_3);
            }
            for (l = 1; l < r4; l++) {
                S2 = S.block(l, 0, 1, r1*r2).transpose();
                S2.resize(r1, r2); S2 = S2.transpose();
                for (k = 1; k < r3; k++) {
                    S_3 = S.block(l, k*r1*r2, 1, r1*r2).transpose();
                    S_3.resize(r1, r2); S_3 = S_3.transpose();
                    S2 = cbind_rcpp(S2, S_3);
                }
                S1 = cbind_rcpp(S1, S2);
            }
        }
        if (d2 == 3) {
            S1 = S.row(0).transpose();
            S1.resize(r1*r2, r3); S1 = S1.transpose();
            for (l = 1; l < r4; l++) {
                S_3 = S.row(l).transpose();
                S_3.resize(r1*r2, r3); S_3 = S_3.transpose();
                S1 = cbind_rcpp(S1, S_3);
            }
        }
    }
    return S1;
}
//----------------------------------------------------------------**
//***--------------------transfer modal 1 to 2 -------------------**
MatrixXd TransferModalUnfoldingsT12(MatrixXd S, VectorXi dim)
{
    int k,d=1, order = dim.size(), r1=dim[0], r2=dim[1];
	for(k=2; k<order; k++) d*=dim[k];
	MatrixXd S1 = S.block(0, 0, r1, r2).transpose();
	for(k=1; k<d; k++) S1 = cbind_rcpp(S1, S.block(0, k*r2, r1, r2).transpose());
	return S1;	
}
//----------------------------------------------------------------**
//***--------------------transfer modal 2 to 1 -------------------**
MatrixXd TransferModalUnfoldingsT21(MatrixXd S, VectorXi dim)
{
    int k,d=1, order = dim.size(), r1=dim[0], r2=dim[1];
	for(k=2; k<order; k++) d*=dim[k];
	MatrixXd S1 = S.block(0, 0, r2, r1).transpose();
	for(k=1; k<d; k++) S1 = cbind_rcpp(S1, S.block(0, k*r1, r2, r1).transpose());
	return S1;	
}
//----------------------------------------------------------------**
//***--------------------transfer modal 1 to d -------------------**
//// [[Rcpp::export]]
MatrixXd TransferModalUnfoldingsT1d(MatrixXd S, int d, VectorXi dim)
{
	if(d==1) return S;
	if(d==2) return TransferModalUnfoldingsT12(S, dim);
	else{
		d=d-1;
		int i,ii,j,jd,k,d1=1, d2=1, order = dim.size(),r1=dim[0], rd=dim[d];
		for(k=1; k<d; k++) d1*=dim[k];
		for(k=d+1; k<order; k++) d2*=dim[k];		
		MatrixXd S1 ,C = MatrixXd::Constant(r1, rd, 0), C1;
		
		for(ii=0;ii<rd;ii++) C.col(ii) = S.col(d1*ii);
		C1 = C.transpose();
		for(i=1;i<d1;i++){
			for(ii=0;ii<rd;ii++) C.col(ii) = S.col(d1*ii+i);
			C1 = cbind_rcpp(C1,C.transpose());
		}
		S1 = C1;
		for(j=1;j<d2;j++){
			jd = j*d1*rd;
			for(ii=0;ii<rd;ii++) C.col(ii) = S.col(jd+d1*ii);
			C1 = C.transpose();
			for(i=1;i<d1;i++){
				for(ii=0;ii<rd;ii++) C.col(ii) = S.col(jd+d1*ii+i);
				C1 = cbind_rcpp(C1,C.transpose());
			}
			S1 = cbind_rcpp(S1, C1);
			
		}
		return S1;
	}
}
//----------------------------------------------------------------**
//***--------------------transfer modal d to 1 -------------------**
//// [[Rcpp::export]]
MatrixXd TransferModalUnfoldingsTd1(MatrixXd S, int d, VectorXi dim)
{
	if(d==1) return S;
	if(d==2) return TransferModalUnfoldingsT21(S, dim);
	else{
		d=d-1;
		int i,ii,j,jd,k,d1=1, d2=1, order = dim.size(),r1=dim[0], rd=dim[d];
		for(k=1; k<d; k++) d1*=dim[k];
		for(k=d+1; k<order; k++) d2*=dim[k];		
		MatrixXd S1, C = MatrixXd::Constant(r1, d1*rd, 0), C0;
		for(i=0; i<d1; i++){
			C0 = S.block(0, i*r1, rd, r1).transpose();
			for(ii=0;ii<rd;ii++) C.col(d1*ii+i) = C0.col(ii);				
		}
		S1 = C;
		for(j=1;j<d2;j++){
			jd = j*d1*r1;
			for(i=0; i<d1; i++){
				C0 = S.block(0, jd+i*r1, rd, r1).transpose();
				for(ii=0;ii<rd;ii++) C.col(d1*ii+i) = C0.col(ii);				
			}
			S1 = cbind_rcpp(S1, C);
		}
		return S1;
	}
}
//----------------------------------------------------------------**
//***--------------------transfer modal d1 to d2 -----------------**
// [[Rcpp::export]]
MatrixXd TransferModalUnfoldingsT(MatrixXd S, int d1, int d2, VectorXi dim)
{
	if(dim.size()<3) stop("S must be greater than 3!");
	if(d1==1)
		return TransferModalUnfoldingsT1d(S, d2, dim);
	else{
		if(d2==1) return TransferModalUnfoldingsTd1(S, d1, dim);
		else{
			MatrixXd S1 = TransferModalUnfoldingsTd1(S, d1, dim);
			return TransferModalUnfoldingsT1d(S1, d2, dim);
		}
	}
}
//----------------------------------------------------------------**
//***----------------------reassign columns of a matrix-----------**
MatrixXd submatrix_col(MatrixXd A, VectorXi b)
{
	int n = A.rows();
	int p = b.size();
	int j;
	MatrixXd C = MatrixXd::Constant(n, p, 0);
	for (j = 0; j < p; j++)
		C.col(j) = A.col(b[j]-1);
	return C;
}
//----------------------------------------------------------------**
//***----------------------reassign rows of a matrix--------------**
MatrixXd submatrix_row(MatrixXd A, VectorXi b)
{
	int p = A.cols();
	int n = b.size();
	int i;
	MatrixXd C = MatrixXd::Constant(n, p, 0);
	for (i = 0; i < n; i++)
		C.row(i) = A.row(b[i] - 1);
	return C;
}
//----------------------------------------------------------------**
//***--------------------penalty----------------------------------**
double penalties(double z, double v, double lambda, double alpha, double gamma, double penalty) {
	double beta=0,l1,l2;
	l1 = lambda*alpha; 
	l2 = lambda*(1-alpha);
	if (penalty==1)
	{
        if (z > l1) beta = (z-l1)/(v*(1+l2));
		if (z < -l1) beta = (z+l1)/(v*(1+l2));
	}
	if (penalty==2)
	{
		double s = 0;
		if (z > 0) s = 1;
		else if (z < 0) s = -1;
		if (fabs(z) <= l1) beta = 0;
		else if (fabs(z) <= gamma*l1*(1+l2)) beta = s*(fabs(z)-l1)/(v*(1+l2-1/gamma));
		else beta = z/(v*(1+l2));
	}
	if (penalty==3)
	{
		double s = 0;
		if (z > 0) s = 1;
		else if (z < 0) s = -1;
        if (fabs(z) <= l1) {
            beta = 0;
        }
        else if (fabs(z) <= (l1*(1+l2)+l1)) {
            beta = s*(fabs(z)-l1)/(v*(1+l2));
        }
        else if (fabs(z) <= gamma*l1*(1+l2)) {
            beta = s*(fabs(z)-gamma*l1/(gamma-1))/(v*(1-1/(gamma-1)+l2));
        }
        else {
            beta = z/(v*(1+l2));
        }
	}
	return(beta);
}
//----------------------------------------------------------------**
//***----update the jth row of matrix A with penalty--------------**
VectorXd updateAj(VectorXd z, int n, int r1, double lambda, double alpha, double gamma, int penalty)
{
	double znorm = 0;
	int j;
	VectorXd b = VectorXd::Constant(r1, 0);
	znorm = z.norm();
    if(znorm==0){
        b = z;
        return b;
    }else if(isinf(znorm)){
        znorm = fabs(gamma*(lambda*alpha)*(1+lambda*(1-alpha)))+1;   // penalty=3
        znorm = penalties(znorm, 1, lambda, alpha, gamma, penalty) / znorm;
        for (j = 0; j<r1; j++) b[j] = znorm * z[j];
        //    for (j = 0; j<r1; j++) b[j] = znorm * z[j]/n;
        return b;
    }else{
        znorm = penalties(znorm, 1, lambda, alpha, gamma, penalty) / znorm;
        for (j = 0; j<r1; j++) b[j] = znorm * z[j];
        //    for (j = 0; j<r1; j++) b[j] = znorm * z[j]/n;
        return b;
    }
}

