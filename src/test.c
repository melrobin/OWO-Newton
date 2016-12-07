#include <stdio.h>
#include <stdlib.h>
extern void elemmult_(double *,double *,double **, int * , int *);
int main()
{
  double *A,*B,**C;
  int i,j,TWO=2;

  A=(double *)calloc(4,sizeof(double));
  B=(double *)calloc(4,sizeof(double));
  C=(double **)calloc(4,sizeof(double));

  A[0]=4;
  A[1]=3;
  A[2]=-2;
  A[3]=2;
  B[0]=1;
  B[1]=0;
  B[2]=5;
  B[3]=.6;
  elemmult_(A,B,C,&TWO,&TWO);
  for (i=0;i<2;i++)
    for (j=0;j<2;j++)
       printf("%ld ",C[i][j]);  
  return(0);
}


