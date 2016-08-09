#ifndef BASIC_MATH_H
#define BASIC_MATH_H

#include <iostream>
#include <vector>
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include <set>
#include <Eigen/Dense>
using namespace Eigen;
using namespace std;

typedef vector<vector<double> > doubleM;
typedef vector<double> doubleV;
typedef std::vector<double>::iterator doubleitV;
typedef std::vector<vector<double> >::iterator doubleitM;


void transpose(doubleM a, doubleM& at);
double Determinant(doubleM a);
double Determinant(double **a,int n);
void CoFactor(doubleM a, doubleM& b);
void invMat(doubleM a, doubleM& b); // Matrix inverse
void matMul(doubleM a, doubleM b, doubleM& c); // Matrix multiply
void matMulc(doubleM a, double b, doubleM& c); // Matrix multiply constant
double vecMul(doubleV a, doubleV b); // vector multiply a' * b
void matSum(doubleM a, doubleM b, doubleM& c); // Summation of mat
void matSub(doubleM a, doubleM b, doubleM& c); // Subtraction of mat
void diag(doubleV a, doubleM& b); // vector to diagonal matrix
void diag(doubleV a, MatrixXd &b);
void vectoMat(doubleV a, doubleM& b); // vertical vector to matrix
void mattoVec(doubleM a, doubleV& b); // matrix (vertical matrix) to vector
void mattoVecR(MatrixXd a, VectorXd& b, int row); // matrix to vector row
void mattoVecC(MatrixXd a, VectorXd& b, int col); // matrix to vector column
void vecSub(doubleV a, doubleV b, doubleV& c); // vector subtraction
double normcdf(double x0, double mu, double sigma);
double lognormpdf(double x0, double mu, double sigma);
void matSum(doubleM &a,int dimension, doubleV &b); // implement sum function in matlab
void matSum(doubleM &a, set<int> idx1, set<int> idx2, int dimension, doubleV &b); // implement sum with idx
void matSum(MatrixXd &a, set<int> idx1, set<int> idx2, int dimension, doubleV &b);
void matSum(MatrixXd &a, VectorXd &b, int dim);
int maxidxV(doubleV &b);

void transpose(doubleM a, doubleM& at)
{
	int N = a.size();
	assert(N > 0);
	int M = a[0].size();
	at.resize(M, doubleV(N));
	doubleitM itat = at.begin();
	for (int i = 0; i < M; i++)
	{
		doubleitV itatj = (*itat).begin();
		for (int j = 0; j < N; j++)
		{
			*itatj = a[j][i];
			std::advance(itatj,1);
		}
		std::advance(itat,1);
	}
}

double Determinant(doubleM a)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;
   int n = a.size();

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = new double*[n-1]; //malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = new double[n-1]; //malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}

double Determinant(double **a,int n)
{
   int i,j,j1,j2;
   double det = 0;
   double **m = NULL;

   if (n < 1) { /* Error */

   } else if (n == 1) { /* Shouldn't get used */
      det = a[0][0];
   } else if (n == 2) {
      det = a[0][0] * a[1][1] - a[1][0] * a[0][1];
   } else {
      det = 0;
      for (j1=0;j1<n;j1++) {
         m = new double*[n-1]; // malloc((n-1)*sizeof(double *));
         for (i=0;i<n-1;i++)
            m[i] = new double[n-1]; // malloc((n-1)*sizeof(double));
         for (i=1;i<n;i++) {
            j2 = 0;
            for (j=0;j<n;j++) {
               if (j == j1)
                  continue;
               m[i-1][j2] = a[i][j];
               j2++;
            }
         }
         det += pow(-1.0,j1+2.0) * a[0][j1] * Determinant(m,n-1);
         for (i=0;i<n-1;i++)
            free(m[i]);
         free(m);
      }
   }
   return(det);
}


void CoFactor(doubleM a, doubleM& b)
{
	int N = a.size();
	assert(N > 0);
	int M = a[0].size();
	b.resize(N, doubleV(M));
	int i,j,ii,jj,i1,j1;
	double det;
	int n = a.size();
	double **c;

	c = new double*[n-1];//malloc((n-1)*sizeof(double *));
	for (i=0;i<n-1;i++)
	 c[i] = new double[n-1];//malloc((n-1)*sizeof(double));

	for (j=0;j<n;j++) {
		cout << "j: " << j << endl;
	  for (i=0;i<n;i++) {

		 /* Form the adjoint a_ij */
		 i1 = 0;
		 for (ii=0;ii<n;ii++) 
		 {
		    if (ii == i)
		       continue;
		    j1 = 0;
		    for (jj=0;jj<n;jj++) 
		    {
		       if (jj == j)
		          continue;
		       c[i1][j1] = a[ii][jj];
		       j1++;
		    }
		    i1++;
		 }

		 /* Calculate the determinate */
		 det = Determinant(c,n-1);

		 /* Fill in the elements of the cofactor */
		 /*doubleitM itb = b.begin();
		 std::advance(itb,i);
		 doubleitV itbj = (*itb).begin();
		 std::advance(itbj,j);
		 (*itbj) = pow(-1.0,i+j+2.0) * det;*/
		 b[i][j] = pow(-1.0,i+j+2.0) * det;
	  }
	}
	for (i=0;i<n-1;i++)
	  free(c[i]);
	free(c);
}

void invMat(doubleM a, doubleM& b)
{
	int N = a.size();
	assert(N > 0);
	int M = a[0].size();
	assert(M > 0);
	assert(N == M);
	doubleM b1;
	cout << "Start CoFactor" << endl;
	CoFactor(a,b1);
	cout << "Finish CoFactor" << endl;
	double det = Determinant(a);
	cout << "Finish Determinant" << endl;
	transpose(b1,b);
	cout << "Finish transpose" << endl;
	doubleitM itb = b.begin();
	for (int i = 0; i < N; i++)
	{
		doubleitV itbj = (*itb).begin();
		for (int j = 0; j < M; j++)
		{	
			(*itbj) = (*itbj) / det;
			std::advance(itbj,1);
		}
		std::advance(itb,1);
	}
	cout << "Finish invMat" << endl;
}

void matMul(doubleM a, doubleM b, doubleM& c)
{
	int Na = a.size();
	assert(Na > 0);
	int Ma = a[0].size();
	assert(Ma > 0);
	int Nb = b.size();
	assert(Nb > 0);
	int Mb = b[0].size();
	assert(Mb > 0);
	assert(Ma == Nb);
	doubleM b1;
	transpose(b,b1);
	c.resize(Na, doubleV(Mb));
	for (int i = 0; i < c.size(); i++)
	{
		for (int j = 0; j < c[0].size(); j++)
		{
			doubleitM itc = c.begin();
			std::advance(itc,i);
			doubleitV itcj = (*itc).begin();
			std::advance(itcj,j);
			(*itcj) = vecMul(a[i],b1[j]);
			//c[i][j] = vecMul(a[i],b1[j]);
		}
	}
}

double vecMul(doubleV a, doubleV b)
{
	double c = 0;
	assert(a.size() == b.size());
	for (int i = 0; i < a.size(); i++)
		c += a[i] * b[i];
	return c;
}

void matSum(doubleM a, doubleM b, doubleM& c)
{
	int Na = a.size();
	assert(Na > 0);
	int Ma = a[0].size();
	assert(Ma > 0);
	int Nb = b.size();
	assert(Nb > 0);
	int Mb = b[0].size();
	assert(Mb > 0);
	assert(Na == Nb && Ma == Mb);
	c.resize(Na, doubleV(Ma));
	for (int i = 0; i < Na; i++)
	{
		for (int j = 0; j < Ma; j++)
		{
			doubleitM itc = c.begin();
			std::advance(itc,i);
			doubleitV itcj = (*itc).begin();
			std::advance(itcj,j);
			(*itcj) = a[i][j] + b[i][j];
			//c[i][j] = a[i][j] + b[i][j];
		}
	}
		
}


void diag(doubleV a, doubleM& b)
{
	int N = a.size();
	assert(N > 0);
	b.resize(N, doubleV(N,0));
	for (int i = 0; i < N; i++)
	{
		doubleitM itb = b.begin();
		std::advance(itb,i);
		doubleitV itbj = (*itb).begin();
		std::advance(itbj,i);
		(*itbj) = a[i];
		//b[i][i] = a[i];
	}
}

void diag(doubleV a, MatrixXd& b)
{
	int N = a.size();
	assert(N > 0);
	b = MatrixXd::Zero(N,N);
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j < N; j++)
		{
			if (i == j)
				b(i,i) = a[i];
			else
				b(i,j) = 0;
		}
	}
}

void vectoMat(doubleV a, doubleM& b)
{
	int N = a.size();
	assert(N > 0);
	b.resize(N, doubleV(1,0));
	for (int i = 0; i < N; i++)
	{
		doubleitM itb = b.begin();
		std::advance(itb,i);
		doubleitV itbj = (*itb).begin();
		(*itbj) = a[i];
		//b[i][0] = a[i];
	}
}

void mattoVec(doubleM a, doubleV& b)
{
	int N = a.size();
	assert(N > 0);
	b.resize(N,0);
	for (int i = 0; i < N; i++)
	{
		doubleitV itb  = b.begin();
		std::advance(itb,i);
		(*itb) = a[i][0];
		//b[i] = a[i][0];
	}
}

void matSub(doubleM a, doubleM b, doubleM& c)
{
	int Na = a.size();
	assert(Na > 0);
	int Ma = a[0].size();
	assert(Ma > 0);
	int Nb = b.size();
	assert(Nb > 0);
	int Mb = b[0].size();
	assert(Mb > 0);
	assert(Na == Nb && Ma == Mb);
	c.resize(Na, doubleV(Ma));
	for (int i = 0; i < Na; i++)
	{
		for (int j = 0; j < Ma; j++)
		{
			doubleitM itc = c.begin();
			std::advance(itc,i);
			doubleitV itcj = (*itc).begin();
			std::advance(itcj,j);
			(*itcj) = a[i][j] - b[i][j];
			//c[i][j] = a[i][j] - b[i][j];
		}
	}

}

void matMulc(doubleM a, double b, doubleM& c)
{
	int Na = a.size();
	assert(Na > 0);
	int Ma = a[0].size();
	assert(Ma > 0);
	c.resize(Na, doubleV(Ma));
	for (int i = 0; i < Na; i++)
	{
		for (int j = 0; j < Ma; j++)
		{
			doubleitM itc = c.begin();
			std::advance(itc,i);
			doubleitV itcj = (*itc).begin();
			std::advance(itcj,j);
			(*itcj) = a[i][j] * b;
			//c[i][j] = a[i][j] * b;
		}
	}

}

void vecSub(doubleV a, doubleV b, doubleV& c)
{
	int N = a.size();
	assert(N>0&&a.size()==b.size());
	c.resize(N,0);
	for (int i = 0; i < N; i++)
	{
		doubleitV itc = c.begin();
		std::advance(itc,i);
		(*itc) = a[i] - b[i];
		//c[i] = a[i] - b[i];
	}
}

double normcdf(double x0, double mu, double sigma)
{
	double x = (x0-mu)/sigma;
    // constants
    double a1 =  0.254829592;
    double a2 = -0.284496736;
    double a3 =  1.421413741;
    double a4 = -1.453152027;
    double a5 =  1.061405429;
    double p  =  0.3275911;
 
    // Save the sign of x
    int sign = 1;
    if (x < 0)
        sign = -1;
    x = fabs(x)/sqrt(2.0);
 
    // A&S formula 7.1.26
    double t = 1.0/(1.0 + p*x);
    double y = 1.0 - (((((a5*t + a4)*t) + a3)*t + a2)*t + a1)*t*exp(-x*x);
 
    return 0.5*(1.0 + sign*y);
}

double lognormpdf(double x0, double mu, double sigma)
{
	//double x = (x0-mu)/sigma;
	//return(-0.5*log(2*3.1415926*pow(sigma,2))-pow(x0,2)/2/pow(sigma,2));
	return(-0.5*log(2*3.1415926*pow(sigma,2))-pow(x0,2)/2/pow(sigma,2));
	//return(-0.5*log(2*3.1415926*1)-pow(x,2)/2/1);
}

void matSum(doubleM &a,int dimension, doubleV &b)
{
	if (dimension == 1)
	{
		b.clear();
		b.resize(a[0].size(),0);
		for (int i = 0; i < b.size(); i++)
		{
			for (int ii = 0; ii < a.size(); ii++)
				b[i] += a[ii][i];
		}
	}
	else if (dimension == 2)
	{
		b.clear();
		b.resize(a.size(),0);
		for (int i = 0; i < b.size(); i++)
		{
			for (int ii = 0; ii < a[0].size(); ii++)
				b[i] += a[i][ii];
		}
	}
}

void matSum(doubleM &a, set<int> idx1, set<int> idx2, int dimension, doubleV &b)
{
	if (dimension == 1)
	{
		b.clear();
		b.resize(idx2.size(),0);

		for (set<int>::iterator iti = idx2.begin(); iti != idx2.end(); iti++)
		{
			for (set<int>::iterator itii = idx1.begin(); itii != idx1.end(); itii++)
				b[*iti] += a[*itii][*iti];
		}
	}
	else if (dimension == 2)
	{
		b.clear();
		b.resize(idx1.size(),0);

		for (set<int>::iterator iti = idx1.begin(); iti != idx1.end(); iti++)
		{
			for (set<int>::iterator itii = idx2.begin(); itii != idx2.end(); itii++)
				b[*iti] += a[*iti][*itii];
		}
	}
}

void matSum(MatrixXd &a, set<int> idx1, set<int> idx2, int dimension, doubleV &b)
{
	if (dimension == 1)
	{
		b.clear();
		b.resize(idx2.size(),0);

		for (set<int>::iterator iti = idx2.begin(); iti != idx2.end(); iti++)
		{
			for (set<int>::iterator itii = idx1.begin(); itii != idx1.end(); itii++)
				b[*iti] += a(*itii,*iti);
		}
	}
	else if (dimension == 2)
	{
		b.clear();
		b.resize(idx1.size(),0);

		for (set<int>::iterator iti = idx1.begin(); iti != idx1.end(); iti++)
		{
			for (set<int>::iterator itii = idx2.begin(); itii != idx2.end(); itii++)
				b[*iti] += a(*iti,*itii);
		}
	}
}

void matSum(MatrixXd &a, VectorXd &b, int dim)
{
	if (dim == 1)
	{
		b = VectorXd::Zero(a.cols());
		for (int i = 0; i < a.cols(); i++)
		{
			double sum = 0;
			for (int j = 0; j < a.rows(); j++)
				sum += a(j,i);
			b(i) = sum;
		}
	}
	else if(dim == 2)
	{
		b = VectorXd::Zero(a.rows());
		for (int i = 0; i < a.rows(); i++)
		{
			double sum = 0;
			for (int j = 0; j < a.cols(); j++)
				sum += a(i,j);
			b(i) = sum;
		}
	}
}

int maxidxV(doubleV &b)
{
	double maxv = -1;
	int maxidx = -1;
	for (int i = 0; i < b.size(); i++)
	{
		if (maxv < b[i])
		{
			maxv = b[i];
			maxidx = i;
		}
	}
	return maxidx;
}

void mattoVecR(MatrixXd a, VectorXd& b, int row) // matrix to vector row
{
	assert(row < a.rows());
	b = VectorXd::Zero(a.cols());
	for (int i = 0; i < a.cols(); i++)
		b(i) = a(row,i);
}
void mattoVecC(MatrixXd a, VectorXd& b, int col) // matrix to vector column
{
	assert(col < a.cols());
	b = VectorXd::Zero(a.rows());
	for (int i = 0; i < a.rows(); i++)
		b(i) = a(i,col);
}

#endif
