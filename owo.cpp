#include <string>
#include <iostream>
#include <valarray>
#include <fstream>
#include <string>
#include <cstring>
#include <cassert>
#include "fortran.h"
#include "mlp.h"
#include "matrix.h"
#include "matops.h"
using namespace std;
extern double ols(int ,matrix &,matrix &,matrix &,int,const valarray<double> &,double *,matrix &);

extern "C"
{
	double CGrad(int,double *,double **,double *,double,double);
	double  **FarAllocateMatrixMemory(int,int);
	void newCGrad(int,int,double *,double **,double *);
}

double MLP::owo()

//INPUT:  R  autocorrelation matrix
//                C  crosscorrelation matrix
//                W  array from inputs to hidden units
//OUTPUT: W_oh array of output weights from hidden units
//                W_oi array of bypass weights

{
    double MSE;
    int i,j,Nu;
            
    Nu = N + Nh + 1;
    valarray<double> E(M);
    matrix A(Nu,Nu),Wa(M,Nu),tempR(Nu,Nu),tempC(M,Nu);
    int ipiv[Nu],ONE=1,info,lwork=-1,rank;
    valarray<double> S(Nu);
    calculate_corr();

    MSE=ols(Nu,Wa,R,C,M,Et,&E[0],A);
  // A=R;
  //tempC=transpose(C);
  //dgesv_(&Nu,&M,&A[0][0],&Nu,ipiv,&Wa[0][0],&Nu,&info);

 // assert(info==0);

        //Remap weights back to conventional arrays
    for (i=0;i<M;i++)
    {
       for (j=0;j<N;j++)
         Woi[i][j]=Wa[i][j];
       for (j=N;j<N+Nh;j++)
         Woh[i][j-N]=Wa[i][j];
                
       Woi[i][N]=Wa[i][N+Nh];
    }
       
   return(MSE);
}

