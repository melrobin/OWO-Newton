// kron.cpp : Defines the entry point for the console application.
#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <valarray>
#include <algorithm>
#include <cassert>
#ifdef _WIN32
  #include <mkl.h>
  #define F77_FUNC(func) (func)
#else
  #include "fortran.h"
  #define F77_FUNC(func) func##_
#endif

#include "matrix.h"
#include "matops.h"
using namespace std;

matrix reshape(const valarray<double> & x, size_t rows,size_t cols)
{
  if (x.size() != rows*cols)
    throw 1;

  matrix temp(cols,rows);
 
  *temp.data=x;
  temp=transpose(temp);
  return(temp);
}

valarray<double> unroll(const matrix &X)
{
  int N=X.data[0].size(),c=X.getCol(),r=X.getRow();
  valarray<double> temp(N),col(r);
  for (int i=0;i < c;i++)
    {
        col=X.data[0][slice(i,r,c)];
       temp[slice(i*r,r,1)]=col;
    }
  return(temp);
}

matrix element_mult(const matrix &A,const matrix &B)
{
  size_t M1=A.getRow(),M2=B.getRow(),N1=A.getCol(),N2=B.getCol();
  if ((M1 != M2) || (N1 != N2))
  {
      cout << "element_mult:  Matrices are of wrong dimensions.  Fix that! " << endl;
      exit(0);
  }
 
  matrix temp(M1,N1);
  (*temp.data) = (*A.data) * (*B.data);

  return(temp);
}

void zero(matrix &X)
{
  int N = X.data[0].size();
  for (int i=0;i<N;i++)
    X.data[0][i]=0.;
}

valarray<double> diag(const matrix &A)
{
   size_t M=A.getRow(),N=A.getCol();
   size_t i;
   valarray<double> temp(M);
   if (M != N)
   {
      cout << "diag:  Matrix must be square" << endl;
      exit(0);
   }
   for (i=0;i<M;i++)
    temp[i]=A[i][i];

   return(temp);

}

double rcond(const matrix &X)
{
   int M=X.getRow();
   int N=X.getCol();
   int i,j;
   int ipiv[M],info,iwork[M];
   matrix A(M,N);
   //const char NORM="1";
   double anorm,result,work[4*N];
   valarray<double> temprow(N);
   A=X;
   for (i=0;i<M;i++)
     temprow += abs(A.extract_row(i));
   anorm = temprow.max();
   A=transpose(A); 
   dgetrf_(&M,&N,&A[0][0],&N,ipiv,&info);
   assert(info==0);
   dgecon_("1",&M,&A[0][0],&M,&anorm,&result,work,iwork,&info);
   assert(info==0);
   return(result);
}
      
matrix vertcat(bool TransposeA, bool TransposeB, const matrix &A, const matrix &B)
{
  matrix tempA,tempB;
  int M1=A.getRow();
  int M2=B.getRow();
  int N1=A.getCol();
  int N2=B.getCol();
  size_t lengths1[2],lengths2[2],strides1[2],strides2[2];
  if (TransposeA)
   {
      lengths1[0]=N1;
      lengths1[1]=M1;
      strides1[0]=1;
      strides1[1]=N1;
      swap(M1,N1);
   }
  else
  { 
     lengths1[0]=M1;
     lengths1[1]=N1;
     strides1[0]=N1;
     strides1[1]=1;
  }

  if (TransposeB)
  {
    lengths2[0]=N2;
    lengths2[1]=M2;
    strides2[0]=1;
    strides2[1]=N2;
    swap(M2,N2);
  }
  else
  {
    lengths2[0]=M2;
    lengths2[1]=N2;
    strides2[0]=N2;
    strides2[1]=1;
  }
  if (N1 != N2)
  {
     cout << "vertcat: dimensions do not match" << endl;
     exit(0);
  }
  matrix temp(M1+M2,N2);
  gslice A_slice(0,valarray<size_t>(lengths1,2),valarray<size_t>(strides1,2));
  gslice B_slice(0,valarray<size_t>(lengths2,2),valarray<size_t>(strides2,2));
  slice A_part(0,M1*N1,1);
  slice B_part(M1*N1,M2*N2,1);
  temp.data[0][A_part]=A.data[0][A_slice];
  temp.data[0][B_part]=B.data[0][B_slice];
 return(temp);
}
matrix horzcat(const matrix &A, const matrix &B)
{
   return(transpose(vertcat(true,true,A,B)));
}
matrix kron(const matrix &A,const matrix & B)
{
  int rows, cols,*f,*g,*h,*e;
  int i,j;
  int ra,rb,ca,cb;
  ra = A.getRow();
  rb = B.getRow();
  ca = A.getCol();
  cb = B.getCol();
        rows = ra*rb;
        cols = ca * cb;
  matrix temp(rows,cols);
//  vector< valarray<double> > vals;
  f = new int[rows];
  g = new int[cols];
  e = new int[rows];
  h = new int[cols];

for (i=0;i<rb;i++)
  for (j=0;j<ra;j++)
        {
        f[i+rb*j]=j;
        e[i+rb*j]=i;
        }

for (i=0;i<cb;i++)
  for (j=0;j<ca;j++)
  {
        g[i+cb*j] = j;
        h[i+cb*j]=i;
  }

for (i=0;i<rows;i++)
  for (j=0;j<cols;j++)
     temp[i][j]=A[f[i]][g[j]]*B[e[i]][h[j]]; 
//for (i=0;i<ra;i++)
// for (j=0;j<ca;j++)
  // vals.push_back(A[i][j]*B.data[0]);

return(temp);
}
/*vector<int> generate_rand( int range_min, int range_max, int n )
{
   // Generate random numbers in the half-closed interval
   // [range_min, range_max). In other words,
   // range_min <= random number < range_max
   int i;
   vector<int> temp;
   for ( i = 0; i < n; i++ )
  
	   temp.push_back((double)rand() / (RAND_MAX + 1) * (range_max - range_min)
            + range_min);
   return(temp);
   
} */
//void get_correct(vector<int> & number,string guess)
//{
//	int i,n = number.size();
//	bool match=false;

//}
//void get_correct_position(vector<int> & number,string guess)
//{
//}
//void process(vector<int> & number, string guess)
//{
//	int n,i;
//	n=guess.length();

//	get_correct(number,guess);
//	get_correct_position(number,guess);
	
//}
vector<double> pack(const matrix &X)
{
  int i,j;
  vector<double> temp;

  for (i=1;i<=X.getRow();i++)
		for (j=0;j<i;j++)
			temp.push_back(X(i-1,j));
	return(temp);
}
matrix transpose(const matrix &X)
{
	matrix temp(X.getCol(),X.getRow());
	int i,j;

	for (i=0;i<X.getCol();i++)
		for (j=0;j<X.getRow();j++)
			temp(i,j)=X(j,i);
	return(temp);
}

matrix chol(const matrix &A)
{
  int c,d,info,i,j;
  matrix C;
	
  C=transpose(A);
			
  c=C.getRow();
  d=C.getCol();
  if (c != d)
    throw 1;
  for (j=1;j<=d;j++)
    for (i=c-1;i>=j;i--)
	C(i,j-1)=0.;
  F77_FUNC(dpotrf)("L",&c,&C(0,0),&c,&info);
  if (info != 0)
    cout << "Not positive definite...you better check yourself!  Minor " << info << " is not working for me" << endl;
  return(C);
}

matrix shrink_matrix(size_t rows, size_t cols,const vector<int> dependent_rows)
{
  size_t column_index = 0,nlin_index=0;
  matrix T(rows-dependent_rows.size(),cols); //Create transformation matrix
  for (size_t row_index=0;row_index< T.getRow();row_index++)
  {
     while (column_index == dependent_rows[nlin_index])
     {
        column_index++;
        nlin_index++;
     }
     T[row_index][column_index++]=1.;
  }
  return(T);
}
matrix qr(const matrix &A)
{
	matrix C=transpose(A);
	int lwork = -1,info,r,c,tau_size;
	vector<double> work(1);
	
	r=A.getRow();
	c=A.getCol();

	tau_size=min(r,c);
	vector<double> tau(tau_size);
	F77_FUNC(dgeqrf)(&r,&c,&C(0,0),&r,&tau[0],&work[0],&lwork,&info); //workspace query
	if (info == 0)
		lwork=(int)work[0];
	else
		cout << "Still some problems in QR factorization" << endl;
	work.resize(lwork);
	F77_FUNC(dgeqrf)(&r,&c,&C(0,0),&r,&tau[0],&work[0],&lwork,&info);
	C=transpose(C);
	return(C);
}

SVDRECORD svd(const matrix &A)
{
  SVDRECORD results;
  int M,N,K;
  int info,lwork=-1; 
  double wkopt;
  M=A.getRow();
  N=A.getCol();
  K=min(M,N);
  valarray<double> S(K);
  matrix temp(M,N),U(M,M),VT(N,N);
  temp=A;
  dgesvd_("A","A",&M,&N,&temp[0][0],&M,&S[0],&U[0][0],&M,&VT[0][0],&N,&wkopt,&lwork,&info);
  lwork = (int)wkopt;
  double work[lwork];
  dgesvd_("A","A",&M,&N,&temp[0][0],&M,&S[0],&U[0][0],&M,&VT[0][0],&N,work,&lwork,&info);
  results.U=matrix(M,M);
  results.V=matrix(N,N);
  results.S=valarray<double>(K);
  results.U=U;
  results.V=VT;
  results.S=S;
  return(results);
}
