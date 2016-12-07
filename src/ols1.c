#include <stdio.h>
#include "allocmem.h"
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <mkl.h>
extern void dump_matrix(int ,int ,char *,double **);
extern void dump_vector(int N,char *fname,double *mat);

void ols_order_bases(int start, int end,double **A,double **R, double **C,
					 int *o,int *oc,int Noc,double **tempW,int N,
					 int Nh,int M,double *E)
{
	int i,j,k,m,n,ILin,Ibest,mbest,NLin=0,ILLin;
	double EEbest,EE,*b,*c,sum,g,*aa,*w,total_error;
	const double eps=1e-6;

	   c=(double *)calloc(N+Nh+1,sizeof(double));
   b=(double *)calloc(N+Nh+1,sizeof(double));
   aa=(double *)calloc(N+Nh+1,sizeof(double));
   w=(double *)calloc(M,sizeof(double));
	
	for (n=start;n<end;n++)
   {
	   EEbest = -1;
	   for (m=0;m<Noc;m++)
	   {
		   ILin = 0;
		   o[n]=oc[m];
		   for (j=0;j<n;j++)
		   {
			   
			   c[j]=0;
			   for (k=0;k<=j;k++)
			   {
				   c[j]+= A[j][k]*R[o[k]][o[n]];
				   
			   }
		   }
		   b[n]=1;
		   for (j=0;j<n;j++)
		   {
			   b[j]=0;
			   for (k=j;k<n;k++)
				   b[j]-=c[k]*A[k][j];
		   }
		   sum=0;

		   for (k=0;k<n;k++)
			   sum +=c[k]*c[k];
		   assert((R[o[n]][o[n]]-sum)>0);
		   g=sqrt(R[o[n]][o[n]]-sum);
		   if (g < eps)
		   {
			   for (k=0;k<=n;k++)
				   A[n][k]=0;
			   ILin = 1;
			   EE=0;
			   for (i=0;i<M;i++)
				   tempW[i][n]=0;
		   }
		   else
		   {
			   for (k=0;k<=n;k++)
				   A[n][k]=b[k]/g;

			   EE=0;
				for (i=0;i<M;i++)
				{
				   tempW[i][n]=0;
				   for (k=0;k<=n;k++)
					   tempW[i][n]+=A[n][k]*C[i][o[k]];
				   EE += tempW[i][n]*tempW[i][n];
				}
		   }
		   if (EE > EEbest)
		   {
			   ILLin=ILin;
			   EEbest=EE;
			   Ibest=o[n];
			mbest=m;
			for (i=0;i<M;i++)
				w[i]=tempW[i][n];

			for (k=0;k<=n;k++)
				aa[k]=A[n][k];
		   }
	   }//end loop over m
	   NLin +=ILLin;
	   // Best value of o[n] has now been found
	   for (i=0;i<M;i++)
		   tempW[i][n]=w[i];

	   for (k=0;k<=n;k++)
		   A[n][k]=aa[k];

	   //Update training errors for the nth basis function
	   total_error = 0;
	   for (i=0;i<M;i++)
		   {
			   E[i]-= tempW[i][n]*tempW[i][n];
			   total_error += E[i];
		   }
	   printf("%d\t%f\t%f\t%d\n",n+1,total_error,EEbest,oc[mbest]);

	   o[n]=Ibest;
	   if (mbest < Noc)
		   for (m=mbest;m<Noc;m++)
			   oc[m]=oc[m+1];
	   Noc--;


   } //end loop over n
}


double ols_train_error(int basis,int M,double **W,double *curr_E,double *energy)
{
   int i,ONE=1;
   double training_error=0;

   for (i=0;i<M;i++)
     {
        if (basis==0)
           curr_E[i]=energy[i];

        curr_E[i]-= W[i][basis]*W[i][basis];
        training_error += curr_E[i];
     }

   return(training_error);
}
int ols(int Nu,double **W,double **R,double **C,int M,double *Et,double *E,
	double **A)

	//Factors Matrix R as W*A
	//In general R can be an NxM matrix

{
   double g,*c,*b,sum;
   double **tempW,MSE;
   int i,j,k,NLin=0,n;
   const double eps=1e-6;
   const  double ONE=1.0,ZERO=0.0;
   tempW=FarAllocateMatrixMemory(M,Nu);
   c=(double *)calloc(Nu,sizeof(double));
   b=(double *)calloc(Nu,sizeof(double));

  
   // Handle the A[0][0] case
	g=sqrt(R[0][0]);
   if (g<eps)
     {
        A[0][0]=0;
        NLin++;
     }
   else
     A[0][0]=1./g;

   //let's start the 2nd basis function
   c[0]=A[0][0]*R[0][1];
   b[1]=1;
   b[0]=-c[0]*A[0][0];
   g=sqrt(R[1][1]-c[0]*c[0]);
   if (g<eps)
    {
       A[1][0]=0;
       A[1][1]=0;
       NLin++;
    }
    else
    {
       A[1][0]=b[0]/g;
       A[1][1]=b[1]/g;
    }      
  
   for (n=2;n < Nu;n++)
   {
      for (j=0;j< n;j++)
      {
        c[j]=0;
        for (k=0;k<=j;k++)
          c[j]+= A[j][k]*R[k][n];
      }
      b[n]=1;
 
	  for (j=0;j< n;j++)
      {
	  b[j]=0;
	  for (k=j;k< n;k++)
             b[j]-=c[k]*A[k][j];
      }
      sum=0;
      for (k=0;k<n;k++)
	sum +=c[k]*c[k];
      g=sqrt(R[n][n]-sum);
      if (g < eps)
      {
	for (k=0;k<n;k++)
	  A[n][k]=0;
	NLin++;
      }
      else
	  for (k=0;k<=n;k++)
		  A[n][k]=b[k]/g;
   }
   dgemm("T","N",&Nu,&M,&Nu,&ONE,&A[0][0],&Nu,&C[0][0],&Nu,&ZERO,
&tempW[0][0],&Nu);

   dgemm("N","N",&Nu,&M,&Nu,&ONE,&A[0][0],&Nu,&tempW[0][0],&Nu,
	&ZERO,&W[0][0],&Nu);             
  
  MSE = 0.0;
  for (i=0;i<M;i++)
     {
		 E[i]=Et[i];
		 for (k=0;k<Nu;k++)
       	    E[i]-= tempW[i][k]*tempW[i][k];
    	MSE+=E[i];
 	printf("%d\t%f\n",i,MSE);
     }
   printf("OLS Basis Error Report\n");
   for (i=0;i<Nu;i++)
     printf("%d\t%f\n",i,ols_train_error(i,M,tempW,E,Et));
   

   free(tempW);
   return(NLin);
}

int ols2(int Nu,double **W,double **R,double **C,int M,double *Et,double *E,
	double **A,int *o)

	//Factors Matrix R as W*A
	//In general R can be an NxM matrix
	//C is M by Nu
	//R is Nu by Nu
	//E and Et are of dimension M

{
   double g,*c,*b,sum;
   double **tempW,*tempColumn,total_error,*node_error;
   int i,j,k,NLin=0,n;
   const double eps=1e-6,ONE=1.0,ZERO=0.0;;
   tempW=FarAllocateMatrixMemory(M,Nu);
   c=(double *)calloc(Nu,sizeof(double));
   b=(double *)calloc(Nu,sizeof(double));
   tempColumn=(double *)calloc(M,sizeof(double));
   node_error = (double *)calloc(M,sizeof(double));
//   E=(double *)calloc(M,sizeof(double));
   
  
   // Handle the A[0][0] case
	g=sqrt(R[o[0]][o[0]]);
   if (g<eps)
     {
        A[0][0]=0;
        NLin++;
     }
   else
     A[0][0]=1./g;

   //*****New *****
   for (i=0;i<M;i++)
	   //tempColumn[i]=A[0][0]*C[o[0]][i];
	   tempColumn[i]=A[0][0]*C[i][o[0]];

   total_error = 0;
   for (i=0;i<M;i++)
   {
	   E[i] = Et[i]- tempColumn[i]*tempColumn[i];
	   total_error += E[i];
   }
   printf("Nu\tE\n");
   printf("%d\t%f\n",1,total_error);



   //let's start the 2nd basis function
   c[0]=A[0][0]*R[o[0]][o[1]];
   b[1]=1;
   b[0]=-c[0]*A[0][0];
   g=sqrt(R[o[1]][o[1]]-c[0]*c[0]);
   if (g<eps)
    {
       A[1][0]=0;
       A[1][1]=0;
       NLin++;
    }
    else
    {
       A[1][0]=b[0]/g;
       A[1][1]=b[1]/g;
    }      
  
   for (i=0;i<M;i++)
   {
	   tempColumn[i]=0;
	   for (j=0;j<2;j++)
		   tempColumn[i]+=A[1][j]*C[i][o[j]];
   }

   total_error = 0;
   for (i=0;i<M;i++)
   {
	   E[i]-= tempColumn[i]*tempColumn[i];
	   total_error += E[i];
   }
   printf("%d\t%f\n",2,total_error);
   for (n=2;n < Nu;n++)
   {
      for (j=0;j< n;j++)
      {
        c[j]=0;
        for (k=0;k<=j;k++)
          c[j]+= A[j][k]*R[o[k]][o[n]];
      }
      b[n]=1;
      for (j=0;j<=n-1;j++)
      {
	  b[j]=0;
	  for (k=j;k<=n-1;k++)
             b[j]-=c[k]*A[k][j];
      }
      sum=0;
      for (k=0;k<=n-1;k++)
		  sum +=c[k]*c[k];
      g=sqrt(R[o[n]][o[n]]-sum);
      if (g < eps)
      {
		  for (k=0;k<n;k++)
			  A[n][k]=0;
		  NLin++;
      }
      else
		  for (k=0;k<=n;k++)
			  A[n][k]=b[k]/g;

	  // New
	    
   for (i=0;i<M;i++)
   {
	   tempColumn[i]=0;
	   for (j=0;j<=n;j++)
		   tempColumn[i]+=A[n][j]*C[i][o[j]];
   }
 
   total_error = 0;
   for (i=0;i<M;i++)
   {
	   E[i]-= tempColumn[i]*tempColumn[i];
	   total_error += E[i];
   }
   printf("%d\t%f\n",n+1,total_error);
   }

   //dgemm("T","N",&Nu,&M,&Nu,&ONE,&A[0][0],&Nu,&C[0][0],&Nu,&ZERO,
//&tempW[0][0],&Nu);

  // dgemm("N","N",&Nu,&M,&Nu,&ONE,&A[0][0],&Nu,&tempW[0][0],&Nu,
	//&ZERO,&W[0][0],&Nu);             
  
 /* MSE = 0.0;
  for (i=0;i<M;i++)
     {
	E[i]=Et[i];
           for (k=0;k<Nu;k++)
       	    E[i]-= tempW[i][k]*tempW[i][k];
    	MSE+=E[i];
 	printf("%d\t%f\n",i,MSE);
     }
   printf("OLS Basis Error Report\n");
   for (i=0;i<Nu;i++)
     printf("%d\t%f\n",i,ols_train_error(i,M,tempW,E,Et));
   */

   free(tempW);
   return(NLin);
}
int ols3(int Nu,double **W,double **R,double **C,int M,double *Et,double *E,
	double **A)

	//Factors Matrix R as W*A
	//In general R can be an NxM matrix
	//C is M by Nu
	//R is Nu by Nu
	//E and Et are of dimension M

{
	double g,*c,*b,sum,EE,EEbest,*aa,*w;
   double **tempW,total_error;
   int i,j,k,NLin=0,Noc,n,*o,m,Ibest,mbest,Nocc,*occ,No,*oc;
   const double eps=1e-6,ONE=1.0,ZERO=0.0;

   tempW=FarAllocateMatrixMemory(M,Nu);
   c=(double *)calloc(Nu,sizeof(double));
   b=(double *)calloc(Nu,sizeof(double));
   o=(int *)calloc(Nu,sizeof(int));
   aa=(double *)calloc(Nu,sizeof(double));
   w=(double *)calloc(M,sizeof(double));
   occ=(int *)calloc(Nu,sizeof(int));
    oc=(int *)calloc(Nu,sizeof(int));


   for (i=0;i<Nu;i++)
	   oc[i]=i;

   o[0]=oc[0];
   Noc=Nu;
  
	g=sqrt(R[o[0]][o[0]]);
     A[0][0]=1./g;

   for (i=0;i<M;i++)
	   W[i][0]=A[0][0]*C[i][o[0]];

   total_error = 0;
   for (i=0;i<M;i++)
   {
	   E[i] = Et[i]- W[i][0]*W[i][0];
	   total_error += E[i];
   }
   printf("Nu\tE\n");
   printf("%d\t%f\n",1,total_error);

	for (i=0;i<Noc-1;i++)
		oc[i]=oc[i+1];
	Noc--;

	No=1;

	if (Noc==0)
		goto step6;

	EEbest=-10;

	for (m=0;m<Noc;m++)
	{
		o[1]=oc[m];
		c[0]=A[0][0]*R[o[0]][o[1]];
		b[1]=1;
		b[0]=-c[0]*A[0][0];
		g=sqrt(R[o[1]][o[1]]-c[0]*c[0]);
		if (g<eps)
		{
			oc[m]=-1;
			EE=-1;
			A[1][0]=0;
			A[1][1]=0;
			for (i=0;i<M;i++)
				tempW[i][1]=0;
		}
		else
		{
			A[1][0]=b[0]/g;
			A[1][1]=b[1]/g;
			EE=0;

			for (i=0;i<M;i++)
			{
				tempW[i][1]=0;
				for (k=0;k<2;k++)
					tempW[i][1]+=A[1][k]*C[i][o[k]];
				EE+= tempW[i][1]*tempW[i][1];
			}
		}

		if (EE > EEbest)
		{
			EEbest=EE;
			Ibest=o[1];
			mbest=m;
			for (i=0;i<M;i++)
				w[i]=tempW[i][1];
			aa[0]=A[1][0];
			aa[1]=A[1][1];
		}
	} //end of m loop

	if (oc[mbest]==-1)
		goto step6;

	for (i=0;i<M;i++)
		tempW[i][1]=w[i];
	A[1][0]=aa[0];
	A[1][1]=aa[1];
	total_error = 0;
   for (i=0;i<M;i++)
   {
	   E[i]-= tempW[i][1]*tempW[i][1];
	   total_error += E[i];
   }
   printf("%d\t%f\t%f\t%d\n",2,total_error,EEbest,oc[mbest]);
   o[1]=Ibest;
   No++;
   oc[mbest]=-1;
	Nocc=0;  //Number of good elements in the array occ
	for (m=0;m<Noc;m++)
		if (oc[m]!=-1)
			occ[Nocc++]=oc[m];

	Noc=Nocc;
	if (Noc==0)
		goto step6;

	for (m=0;m<Noc;m++)
		oc[m]=occ[m];

   for (n=2;n < Nu;n++)
   {
	   EEbest=-10;
	   for (m=0;m<Noc;m++)
	   {
		   o[n]=oc[m];
		   for (j=0;j<n;j++)
		   {
			   c[j]=0;
			   for (k=0;k<=j;k++)
				   c[j]+= A[j][k]*R[o[k]][o[n]];
		   }
		   b[n]=1;
		   for (j=0;j<n;j++)
		   {
			   b[j]=0;
			   for (k=j;k<n;k++)
				   b[j]-=c[k]*A[k][j];
		   }
		   sum=0;
		   for (k=0;k<n;k++)
			   sum +=c[k]*c[k];
		   assert((R[o[n]][o[n]]-sum)>0);
		   g=sqrt(R[o[n]][o[n]]-sum);
		   if (g < eps)
		   {
			   oc[m]=-1;
			   for (k=0;k<=n;k++)
				   A[n][k]=0;
			   EE = -1;
			   for (i=0;i<M;i++)
				   tempW[i][n]=0;
		   }
		   else
		   {
			   for (k=0;k<=n;k++)
				   A[n][k]=b[k]/g;

			   EE=0;
				for (i=0;i<M;i++)
				{
				   tempW[i][n]=0;
				   for (k=0;k<=n;k++)
					   tempW[i][n]+=A[n][k]*C[i][o[k]];
				   EE += tempW[i][n]*tempW[i][n];
				}
		   }
		   if (EE > EEbest)
		   {
			   EEbest=EE;
			   Ibest=o[n];
			mbest=m;
			for (i=0;i<M;i++)
				w[i]=tempW[i][n];

			for (k=0;k<=n;k++)
				aa[k]=A[n][k];
		   }
	   } //end loop over m
	   if (oc[mbest]==-1)
		   goto step6;

	   for (i=0;i<M;i++)
		   tempW[i][n]=w[i];

	   for (k=0;k<=n;k++)
		   A[n][k]=aa[k];

	   total_error = 0;
	   for (i=0;i<M;i++)
		   {
			   E[i]-= tempW[i][n]*tempW[i][n];
			   total_error += E[i];
		   }
	   printf("%d\t%f\t%f\t%d\n",n+1,total_error,EEbest,oc[mbest]);

	   o[n]=Ibest;
	   No++;
	   oc[mbest]=-1;
	   Nocc=0;
	   for (m=0;m<Noc;m++)
		   if (oc[m]!=-1)
			   occ[Nocc++]=oc[m];

	   Noc=Nocc;
	   if (Noc==0)
		goto step6;

	   for (m=0;m<Noc;m++)
		   oc[m]=occ[m];
   } //end of loop over n
step6:   printf("# good basis fcns= %d\n",No);
   NLin=Nu-No;
   printf("# dep basis fcns = %d\n",NLin);
   for (i=0;i<M;i++)
	   for (k=0;k<Nu;k++)
		   W[i][k]=0;

   for (i=0;i<M;i++)
	   for (k=0;k<No;k++)
		   for (m=k;m<No;m++)
			   W[i][o[k]] += A[m][k]*tempW[i][m];


   free(c);
   free(b);
   free(o);
   free(occ);
   free(oc);
   free(w);
   free(aa);
   free(tempW);
   return(NLin);
}
int ols4(int N,int Nh,double **W,double **R,double **C,int M,double *Et,double *E,
	double **A)

	//Factors Matrix R as W*A
	//In general R can be an NxM matrix
	//C is M by Nu
	//R is Nu by Nu
	//E and Et are of dimension M

{
	double g,*c,*b,sum,EE,EEbest,*aa,*w;
   double **tempW,total_error;
   int i,j,k,NLin=0,Noc,n,*o,m,Ibest,mbest,*oc;
   int Nu,ILLin,ILin=0;
   const double eps=1e-6,ONE=1.0,ZERO=0.0;

   tempW=FarAllocateMatrixMemory(M,N+Nh+1);
   c=(double *)calloc(N+Nh+1,sizeof(double));
   b=(double *)calloc(N+Nh+1,sizeof(double));
   o=(int *)calloc(N+Nh+1,sizeof(int));
   aa=(double *)calloc(N+Nh+1,sizeof(double));
   w=(double *)calloc(M,sizeof(double));
   oc=(int *)calloc(N+Nh+1,sizeof(int));

   Noc=N+1;
   Nu=Noc;

   for (i=0;i<Nu;i++)
	   oc[i]=i;

   o[0]=oc[0];
   
  
	g=sqrt(R[o[0]][o[0]]);
     A[0][0]=1./g;

   for (i=0;i<M;i++)
	   W[i][0]=A[0][0]*C[i][o[0]];

   total_error = 0;
   for (i=0;i<M;i++)
   {
	   E[i] = Et[i]- W[i][0]*W[i][0];
	   total_error += E[i];
   }
   printf("\nOLS4\nNu\tE\n");
   printf("%d\t%f\n",1,total_error);

	for (i=0;i<Noc-1;i++)
		oc[i]=oc[i+1];
	Noc--;
ols_order_bases(1,Nu,A,R,C,o,oc,Noc,tempW,N,Nh,M,E);
  /* for (n=1;n < Nu;n++)
   {
	   EEbest=-1;
	   for (m=0;m<Noc;m++)
	   {
		   ILin=0;
		   o[n]=oc[m];
		   for (j=0;j<n;j++)
		   {
			   c[j]=0;
			   for (k=0;k<=j;k++)
				   c[j]+= A[j][k]*R[o[k]][o[n]];
		   }
		   b[n]=1;
		   for (j=0;j<n;j++)
		   {
			   b[j]=0;
			   for (k=j;k<n;k++)
				   b[j]-=c[k]*A[k][j];
		   }
		   sum=0;
		   
		   
		   for (k=0;k<n;k++)
			   sum +=c[k]*c[k];
		   assert((R[o[n]][o[n]]-sum)>0);
		   g=sqrt(R[o[n]][o[n]]-sum);
		   if (g < eps)
		   {
			   for (k=0;k<=n;k++)
				   A[n][k]=0;
			   ILin=1;
			   EE=0;
			   for (i=0;i<M;i++)
				   tempW[i][n]=0;
		   }
		   else
		   {
			   for (k=0;k<=n;k++)
				   A[n][k]=b[k]/g;

			   EE=0;
				for (i=0;i<M;i++)
				{
				   tempW[i][n]=0;
				   for (k=0;k<=n;k++)
					   tempW[i][n]+=A[n][k]*C[i][o[k]];
				   EE += tempW[i][n]*tempW[i][n];
				}
		   }
		   if (EE > EEbest)
		   {
			   ILLin=ILin;
			   EEbest=EE;
			   Ibest=o[n];
			   mbest=m;
			   for (i=0;i<M;i++)
				   w[i]=tempW[i][n];

			for (k=0;k<=n;k++)
				aa[k]=A[n][k];
		   }
	   } //end loop over m

	   NLin +=ILLin;
	   // Best value of o[n] has now been found
	   for (i=0;i<M;i++)
		   tempW[i][n]=w[i];

	   for (k=0;k<=n;k++)
		   A[n][k]=aa[k];

	   //Update training errors for the nth basis function
	   total_error = 0;
	   for (i=0;i<M;i++)
		   {
			   E[i]-= tempW[i][n]*tempW[i][n];
			   total_error += E[i];
		   }
	   printf("%d\t%f\t%f\t%d\n",n+1,total_error,EEbest,oc[mbest]);

	   o[n]=Ibest;
	   if (mbest < Noc-1)
		   for (m=mbest;m<Noc-1;m++)
			   oc[m]=oc[m+1];
	   Noc--;

   } //end of loop over n */

// Order the hidden units
   for (n=0;n<Nh;n++)
	   oc[n] = N + 1 + n;

   Nu = N + Nh + 1;
   Noc = Nh;
ols_order_bases(N+1,Nu,A,R,C,o,oc,Noc,tempW,N,Nh,M,E);
  /* for (n=N+1;n<Nu;n++)
   {
	   EEbest = -1;
	   for (m=0;m<Noc;m++)
	   {
		   ILin = 0;
		   o[n]=oc[m];
		   for (j=0;j<n;j++)
		   {
			   
			   c[j]=0;
			   for (k=0;k<=j;k++)
			   {
				   c[j]+= A[j][k]*R[o[k]][o[n]];
				   
			   }
		   }
		   b[n]=1;
		   for (j=0;j<n;j++)
		   {
			   b[j]=0;
			   for (k=j;k<n;k++)
				   b[j]-=c[k]*A[k][j];
		   }
		   sum=0;

		   for (k=0;k<n;k++)
			   sum +=c[k]*c[k];
		   assert((R[o[n]][o[n]]-sum)>0);
		   g=sqrt(R[o[n]][o[n]]-sum);
		   if (g < eps)
		   {
			   for (k=0;k<=n;k++)
				   A[n][k]=0;
			   ILin = 1;
			   EE=0;
			   for (i=0;i<M;i++)
				   tempW[i][n]=0;
		   }
		   else
		   {
			   for (k=0;k<=n;k++)
				   A[n][k]=b[k]/g;

			   EE=0;
				for (i=0;i<M;i++)
				{
				   tempW[i][n]=0;
				   for (k=0;k<=n;k++)
					   tempW[i][n]+=A[n][k]*C[i][o[k]];
				   EE += tempW[i][n]*tempW[i][n];
				}
		   }
		   if (EE > EEbest)
		   {
			   ILLin=ILin;
			   EEbest=EE;
			   Ibest=o[n];
			mbest=m;
			for (i=0;i<M;i++)
				w[i]=tempW[i][n];

			for (k=0;k<=n;k++)
				aa[k]=A[n][k];
		   }
	   }//end loop over m
	   NLin +=ILLin;
	   // Best value of o[n] has now been found
	   for (i=0;i<M;i++)
		   tempW[i][n]=w[i];

	   for (k=0;k<=n;k++)
		   A[n][k]=aa[k];

	   //Update training errors for the nth basis function
	   total_error = 0;
	   for (i=0;i<M;i++)
		   {
			   E[i]-= tempW[i][n]*tempW[i][n];
			   total_error += E[i];
		   }
	   printf("%d\t%f\t%f\t%d\n",n+1,total_error,EEbest,oc[mbest]);

	   o[n]=Ibest;
	   if (mbest < Noc)
		   for (m=mbest;m<Noc;m++)
			   oc[m]=oc[m+1];
	   Noc--;


   } //end loop over n
   */
   printf("# dep basis fcns = %d\n",NLin);
   for (k=0;k<=N;k++)
	   for (i=0;i<M;i++)
	   {
	   
		   W[i][k]=0;
		   for (m=k;m<Nu;m++)
			   W[i][o[k]] += A[m][k]*tempW[i][m];
	   }

   for (k=N+1;k<Nu;k++)
	   for (i=0;i<M;i++)
	   {
		   W[i][k]=0;
		   for (m=k;m<Nu;m++)
			   W[i][o[k]] += A[m][k]*tempW[i][m];
	   }

   for (k=0;k<=N;k++)
	   for (i=0;i<M;i++)
	   {
		   W[i][o[k]]=0;
		   for (j=k;j<Nu;j++)
			   W[i][o[k]]+=A[j][k]*tempW[i][j];
	   }

   for (k=N+1;k<Nu;k++)
	   for (i=0;i<M;i++)
	   {
		   W[i][o[k]]=0;
		   for (j=k;j<Nu;j++)
			   W[i][o[k]]+=A[j][k]*tempW[i][j];
	   }


   free(c);
   free(b);
   free(o);
   free(oc);
   free(w);
   free(aa);
   free(tempW);
   return(NLin);
}
    

