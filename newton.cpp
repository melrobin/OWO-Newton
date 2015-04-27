#include <iostream>
#include <valarray>
#include <cassert>
#include "matrix.h"
#ifdef _WIN32
  #include <mkl.h>
  #define F77_FUNC(func) (func)
#else
  #include "fortran.h"
  #define F77_FUNC(func) func##_
#endif

//#include "fortran.h"
using namespace std;
extern double ols(int ,matrix &,matrix &,matrix &,int,const valarray<double> &,double *,matrix &);
int ols_linear_solve(valarray<double> &,matrix &,valarray<double>,vector<int> & );

valarray<double> compute_newton_direction(const matrix &H,const valarray<double> &g, size_t& num_dependent_rows)
{
  int M=H.getRow(),N=H.getCol(),K=g.size(),ONE=1,info,lwork=-1,rank;
  double *work,rcond=1e-16,work_size,MINUSONE=-1.;
  valarray<int> jpvt(N);
  int iwork_size;
  vector<int> dependent_rows;
  valarray<double> temp(M),tempg(K),S(M);
  if ((K != N) || (M != N))
     throw 1;
  matrix tempH(M,M);
  tempH=H;
  tempg=g;
  //num_dependent_rows=0;
 // num_dependent_rows=ols_linear_solve(tempg,tempH,g,dependent_rows);
  F77_FUNC(dgelsd)(&M,&N,&ONE,&tempH[0][0],&N,&tempg[0],&N,&S[0],&MINUSONE,&rank,&work_size,&lwork,&iwork_size,&info);
 // F77_FUNC(dgelsy)(&M,&N,&ONE,&tempH[0][0],&N,&tempg[0],&N,&jpvt[0],&rcond,&rank,&work_size,&lwork,&info);
  assert(info==0);
  int iwork[iwork_size];
  lwork=int(work_size);
  work = new double[lwork];
 // F77_FUNC(dgelsy)(&M,&N,&ONE,&tempH[0][0],&N,&tempg[0],&N,&jpvt[0],&rcond,&rank,work,&lwork,&info);
  F77_FUNC(dgelsd)(&M,&N,&ONE,&tempH[0][0],&N,&tempg[0],&N,&S[0],&MINUSONE,&rank,work,&lwork,iwork,&info);

  assert(info==0); 
  num_dependent_rows=rank;
  delete work;
  return(tempg);
}

matrix compute_HWO_direction(const matrix& R,const matrix& G_bp, size_t& num_dependent_rows)
{
  int info;
  int M1=R.getRow(),N1=R.getCol(),M2=G_bp.getRow(),N2=G_bp.getCol();
  matrix tempR(M1,N1),tempG_bp(M2,N2);
  vector<int> dependent_rows;
  num_dependent_rows=0;
  if ((M1 != N1) || (N1 != N2))
     throw 1;
  tempR=R;
  tempG_bp=G_bp;
  dposv_("U",&N1,&M2,&tempR[0][0],&N1,&tempG_bp[0][0],&N2,&info);
  assert(info==0);
  return(tempG_bp);
}

