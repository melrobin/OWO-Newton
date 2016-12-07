// kron.cpp : Defines the entry point for the console application.
#include <string>
#include <iostream>
#include <cstdlib>
#include <vector>
#include <mkl.h>
#include "matrix.h"
using namespace std;
//extern "C"
//{
  //  void dump_matrix(int ,int ,const char *,double **);
    //double  **FarAllocateMatrixMemory(int , int );
//}
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
  ra = A.getRow();
  rb = B.getRow();
  ca = A.getCol();
  cb = B.getCol();
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
		(*temp)(i,j)=A(f[i],g[j])*B(e[i],h[j]);

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
	matrix C;
	
	C=transpose(A);
			
	c=C.getRow();
	d=C.getCol();
	if (c != d)
		throw 1;
	for (j=1;j<=d;j++)
		for (i=c-1;i>=j;i--)
			C(i,j-1)=0.;
	dpotrf("L",&c,&C(0,0),&c,&info);
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
	dgeqrf(&r,&c,&C(0,0),&r,&tau[0],&work[0],&lwork,&info); //workspace query
	if (info == 0)
		lwork=(int)work[0];
	else
		cout << "Still some problems in QR factorization" << endl;
	work.resize(lwork);
	dgeqrf(&r,&c,&C(0,0),&r,&tau[0],&work[0],&lwork,&info);
	C=transpose(C);
	return(C);
}

int main()
{
  matrix A(3,2),B(3,3),C;
  vector<int> number,guess;
  vector<double> packedB;
  vector<string> guesses;
	  int i,c,info,Bsize;

  A(0,0)=1;
  A(0,1) = .5;
  A(1,0)=4;
  A(1,1)=-.2;
  A(2,0)=-.3;
  A(2,1)=17;
  B(0,0)=B(1,0)=B(0,1)=3;
  B(0,2)=B(2,0)=0;
  B(1,2)=B(2,1)=7;
  B(1,1)=10;
  B(2,2)=8;

  C=chol(B);
  cout << C << endl;
        C=kron(A,B);
       cout << C << endl; 
		C=qr(B);
		cout << C << endl;
		number=generate_rand(0,10,6);
		for (i=0;i<number.size();i++)
			cout << number[i];
		cout << endl;
		
		string s;
		getline(cin, s,'\n');
		cout << "read:" << s << endl;
		guesses.push_back(s); //store the guess
		process(guess,s);
		system("pause");
        return 0;
}


