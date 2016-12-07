#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include "matrix.h"
#include "fortran.h"
using namespace std;
// NAME: KRON
//PURPOSE:  IMPLEMENT THE KRONECKER TENSOR PRODUCT
//PARAMETERS: 
//             A       AN RA X CA DOUBLE PRECISION INPUT MATRIX
//             B       AN RB X CB DOUBLE PRECISION INPUT MATRIX
//             C       AN RA X RB BY CA X CB DOUBLE PRECISION INPUT MATRIX
//             ALPHA   A SCALAR MULTIPLIER
//OPERATION:
//             C = ALPHA * C + KRON(A,B)
//REFERENCES:
//            IMPLEMENTS THE KRON FUNCTION AND IS BASED ON THE OCTAVE SOURCE
//            CODE GIVEN BELOW
//            function c = mykron(a,b)
//                  [ra,ca]=size(a);
//                  [rb,cb]=size(b);
//  
//                  c = a((ones(rb,1)*[1:ra])(:),(ones(cb,1)*[1:ca])(:)).*...
//                      b(([1:rb]'*ones(1,ra))(:),([1:cb]'*ones(1,ca))(:));

matrix kron(const matrix &A,const matrix & B)
{
  unsigned int rows, cols,*f,*g,*h,*e;
  unsigned int i,j;
  unsigned int ra,rb,ca,cb;
  matrix *temp;
  ra = A.nrows;
  rb = B.nrows;
  ca = A.ncols;
  cb = B.ncols;
        rows = ra*rb;
        cols = ca * cb;
  temp = new matrix(rows,cols);
  f = new unsigned int[rows];
  g = new unsigned int[cols];
  e = new unsigned int[rows];
  h = new unsigned int[cols];

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
        temp->data[i*cols+j]=A(f[i],g[j])*B(e[i],h[j]);
       temp->nrows=rows;
     temp->ncols=cols;
       return(*temp);
}
vector<int> generate_rand( int range_min, int range_max, int n )
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
   
}
void get_correct(vector<int> & number,string guess)
{
	int i,n = number.size();
	bool match=false;

}
void get_correct_position(vector<int> & number,string guess)
{
}
void process(vector<int> & number, string guess)
{
	int n,i;
//	n=guess.length();

	get_correct(number,guess);
	get_correct_position(number,guess);
	
}
vector<double> pack(const matrix &X)
{
	unsigned int i,j;
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
	matrix C=transpose(A);
			
	c=C.getRow();
	d=C.getCol();
	if (c != d)
		throw 1;
	for (j=1;j<=d;j++)
		for (i=c-1;i>=j;i--)
			C(i,j-1)=0.;
	dpotrf_("L",&c,C.data,&c,&info);
	if (info != 0)
		cout << "Not positive definite...you better check yourself!" << endl;
	return(C);
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
	dgeqrf_(&r,&c,C.data,&r,&tau[0],&work[0],&lwork,&info); //workspace query
	if (info == 0)
		lwork=(int)work[0];
	else
		cout << "Still some problems in QR factorization" << endl;
	work.resize(lwork);
	dgeqrf_(&r,&c,C.data,&r,&tau[0],&work[0],&lwork,&info);
	C=transpose(C);
	return(C);
}


