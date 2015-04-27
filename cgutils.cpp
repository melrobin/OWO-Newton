#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <assert.h>
#include <valarray>
#include "matrix.h"
using namespace std;
double* createArray(int);
double** createMatrix(int, int);
void freeArray(double *);
void freeMatrix(double **, int);
double** Matr(int , int );
void ArrZero(valarray<double> &, int);
void MatZero(double **, int, int );
double Schmidt(double **, double **, double **, double *, int , int );
double* Arr(int size)
{
	double *ptr = (double *)malloc(sizeof(double) * size);
	if(!ptr)
	{
		printf("Memory allocation error!\n");
		exit(-1);
	}
	return ptr;
}
double **cal_opw(double **gm,double **hm,int M,int Nh)
{
	double e1=0.000001,e2,g,sum;
	double **am,*cm,*bm,**z1,**zs;
	int i,j,k;
	am=Matr(Nh,Nh);
	cm=Arr(Nh);
	bm=Arr(Nh);
	z1=Matr(M,Nh);
	zs=Matr(M,Nh);
	e2=e1*e1;
	g=hm[0][0];
	if(g>0)
	{
		g=sqrt(g);
	}
	for(i=0;i<Nh;i++)
	{  
		bm[i]=0.0;
		cm[i]=0.0;
		for(j=0;j<Nh;j++)
		{
			am[i][j]=0.0;
		}
	}
	am[0][0]=1/g;

	cm[0]=am[0][0]*hm[0][1];
	g=hm[1][1]-cm[0]*cm[0];
	g=sqrt(g);
	if(g<e2)
	{
		am[1][0]=0.0;
		am[1][1]=0.0;
	}
	else
	{
		am[1][0]=(-cm[0]*am[0][0])/g;
		am[1][1]=1/g;
	}
	for(i=2;i<Nh;i++)
	{
		for(j=0;j<i;j++)
		{	cm[j]=0.0;
		for(k=0;k<j+1;k++)
			cm[j]+=am[j][k]*hm[k][i];
		}
		bm[i]=1;
		for(j=0;j<i;j++)
		{
			bm[j]=0;
			for(k=j;k<i;k++)
				bm[j]-=cm[k]*am[k][j];
		}
		sum=0.0;
		for(k=0;k<i;k++)
			sum+=cm[k]*cm[k];
		g=hm[i][i]-sum;
		if(g>0)
		{
			g=sqrt(g);
		}
		if(g<e2)
		{
			for(k=0;k<i+1;k++)
				am[i][k]=0.0;
		}
		else
		{
			for(k=0;k<i+1;k++)
				am[i][k]=bm[k]/g;
		}
	}

	for(i=0;i<M;i++)
		
		for(k=0;k<Nh;k++)
		{
			z1[i][k]=0.0;
			zs[i][k]=0.0;
		}
	
	for(i=0;i<M;i++)
	
		for(k=0;k<Nh;k++)
		{	//z1[i][k]=0.0;
			for(j=0;j<k+1;j++)
				z1[i][k]+=am[k][j]*gm[i][j];
		}
	
	for(i=0;i<M;i++)
	{
		for(k=0;k<Nh;k++)
		{	//zs[i][k]=0.0;
			for(j=k;j<Nh;j++)
				zs[i][k]+=am[j][k]*z1[i][j];
		}
	}
	return zs;
	//fMatr(am,Nh);
	//fArr(cm);
	//fArr(bm);
	//fMatr(z1,Nh);

}

double *cal_molf(double *gm,double **hm,int Nh)
{
	double e1=0.000001,e2,g,sum;
	double **am,*cm,*bm,*z1,*zs;
	int i,j,k;
	am=Matr(Nh,Nh);
	cm=Arr(Nh);
	bm=Arr(Nh);
	z1=Arr(Nh);
	zs=Arr(Nh);
	e2=e1*e1;
	g=hm[0][0];
	if(g>0)
	{
		g=sqrt(g);
	}
	for(i=0;i<Nh;i++)
	{  
		bm[i]=0.0;
		cm[i]=0.0;
		for(j=0;j<Nh;j++)
		{
			am[i][j]=0.0;
		}
	}
	if(g<e2)
	{
		am[0][0]=0.0;
	}		
	else
	{
		am[0][0]=1/g;
	}
	if(Nh>1)
	{

		cm[0]=am[0][0]*hm[0][1];
		g=hm[1][1]-cm[0]*cm[0];
		if(g>0)
		{
			g=sqrt(g);
		}
		if(g<e2)
		{
			am[1][0]=0.0;
			am[1][1]=0.0;
		}
		else
		{
			am[1][0]=(-cm[0]*am[0][0])/g;
			am[1][1]=1/g;
		}
		for(i=2;i<Nh;i++)
		{
			for(j=0;j<i;j++)
			{	cm[j]=0.0;
			for(k=0;k<j+1;k++)
				cm[j]+=am[j][k]*hm[k][i];
			}
			bm[i]=1;
			for(j=0;j<i;j++)
			{
				bm[j]=0;
				for(k=j;k<i;k++)
					bm[j]-=cm[k]*am[k][j];
			}
			sum=0.0;
			for(k=0;k<i;k++)
				sum+=cm[k]*cm[k];
			g=hm[i][i]-sum;
			if(g>0)
			{
				g=sqrt(g);
			}
			if(g<e2)
			{
				for(k=0;k<i+1;k++)
					am[i][k]=0.0;
			}
			else
			{
				for(k=0;k<i+1;k++)
					am[i][k]=bm[k]/g;
			}
		}
	}

	for(k=0;k<Nh;k++)
	{	z1[k]=0.0;
	for(j=0;j<k+1;j++)
		z1[k]+=am[k][j]*gm[j];
	}   

	for(k=0;k<Nh;k++)
	{	zs[k]=0.0;
	for(j=k;j<Nh;j++)
		zs[k]+=am[j][k]*z1[j];
	}

	return zs;
	//fArr(z1);
	//fArr(bm);
	//fArr(cm);
	//fMatr(am,Nh);
}
double** Matr(int rows, int cols)
{
	int i;
	double **mat;

	mat = (double **)malloc(sizeof(double*) * rows);
	if(mat==NULL)
	{
		printf("Memory allocation error.\n");
		exit(-1);
	}
	for(i=0; i<rows; i++)
	{
		mat[i] = (double *)malloc(sizeof(double) * cols);
		if(mat[i]==NULL)
		{
			printf("Memory allocation error.\n");
			exit(-1);
		}
	}

	return mat;
}
void ArrZero(valarray<double> &AR, int size)
{
	int i;
	for(i=0;i<size;i++)
		AR[i]=0;
}

void MatZero(double **AR, int rows, int cols)
{
	for(int i=0;i<rows;i++)
		for(int j=0; j<cols; j++)
			AR[i][j]=0;
}
// This function reads a pattern from a text file
int read1(FILE *infile, int N, int M, valarray<double> &x, valarray<double> &ty)
{
	int readErr = 0;
	for(int n=0; n<N; n++) {
		if(fscanf(infile, "%lf", &x[n]) < 1)
			readErr = 1;
	}
	if(readErr)
		return -1;

	readErr = 0;
	for(int i=0; i<M; i++) {
		if(fscanf(infile, "%lf", &ty[i]) < 1)
			readErr = 1;
	}
	if(readErr)
		return -1;

	return 0;
}

// This macro rduces the number of parameters from the function call
//#define read(infile,x,t) read1(infile, N, M, x, t)

void net_calc(valarray<double> &net, valarray<double> &o, valarray<double> &x, matrix &wih, int Nh, int N)
{
	for(int k=0; k<Nh; k++)
	{	
		net[k] = 0.0;
		for(int n=0; n<N+1; n++)
			net[k] += wih[k][n]*x[n];
		o[k] = 1/(1+exp(-net[k]));
		x[k+N+1] = o[k];
	}
}

//#define NET_CALC net_calc(net, o, x, wih, Nh, N)

 void y_calc(valarray<double>& x, valarray<double>& y, matrix &w, int L, int M)
{
	for(int k=0; k<M; k++)
	{
		y[k]=0.0;
		for(int j=0;j<L;j++)
			y[k]+=w[k][j]*x[j];
	}
}

//#define Y_CALC y_calc(x, y, w, L, M)

 double mse_calc(double *t, double *y, int M)
{
	double ret = 0;
	for(int i=0; i<M; i++)
	{
		ret += (t[i] - y[i])*(t[i] - y[i]);
	}
	return ret;
}

double hardProd(double **a, double **b, int rows, int cols)
{
	double prod=0.0;
	int i, j;
	for(i=0; i<rows; i++)
	{
		for(j=0; j<cols; j++)
		{
			prod+=a[i][j]*b[i][j];
		}
	}
	return prod;
}

double **Multiply(double **A, double **B, int r1, int c1, int c2, double **Result)
{
	for(int i=0; i<r1; i++)
	{
		for(int j=0; j<c2; j++)
		{
			Result[i][j] = 0.0;
			for(int k=0; k<c1; k++)
			{
				Result[i][j] += A[i][k] * B[k][j];
			}
		}
	}

	return Result;
}

double OLS1(double **Wout, int N, int M, int L, double **R, double **C,std::valarray<double> &Et)
{
	int NLin=0 ;
	int n,m,j,i,k,Nv=0;
	double g,b1,b2,c1,g1=0;  
	double *E,*b,*c;
	double **Wt,**A;
	double e;

	A=createMatrix(L, L);
	Wt=createMatrix(M, L);
	b=createArray(L);
	c=createArray(L);
	E = createArray(M);

	//	Initializing the weights and thresholds, Auto and Cross correlation matrices to zero

	for (i = 0; i < L; i++)
	{
		for (j = 0; j < L; j++)
			A[i][j]=0.0;


		for (j = 0; j < M; j++)
			Wt[j][i] = 0.0;
	
	}

	for (i = 0; i < M; i++ )
		E[i] = 0.0;		

	Nv=0;	


	//OLS TYPE 1 CODE(OWO with no reordering).
	g=R[0][0];
	if(g < pow(10.0,-12.0))
	{
		A[0][0]=0;
		NLin=NLin+1;
	}

	else
	{
		g=sqrt(g);
		A[0][0]=(1/g);
	}

	c1=A[0][0]*(R[0][1]);
	b2=1;	
	b1=((-c1)*A[0][0]);
	g=R[1][1]-(c1*c1);

	if(g < pow(10.0,-12.0))
	{
		A[1][0]=0;
		A[1][1]=0;
		NLin=NLin+1;
	}
	else
	{
		g=sqrt(g);
		A[1][0]=(1/g)*b1;
		A[1][1]=(1/g)*b2;
	}

	for(n=2;n<L;n++)
	{
		for(j=0;j<=(n-1);j++)
		{
			c[j]=0;
			for(k=0;k<=j;k++)
				c[j]+=A[j][k]*R[k][n];
		}
		b[n]=1;

		for(j=0;j<=(n-1);j++)//change
		{
			b[j]=0;
			for(k=j;k<=(n-1);k++)
				b[j]=b[j]-c[k]*A[k][j];
		}


		g=R[n][n];
		g1=0;
		for(k=0;k<=n-1;k++)
			g1=g1+(c[k]*c[k]);
		g=g-g1;
                assert(g>0);
		if(g <pow(10.0,-12.0))
		{
			for(k=0;k<=n;k++)		// change n-1 to n : rohit 5/20/2011
				A[n][k]=0;
			NLin=NLin+1;
		}
		else
		{
			g=1/(sqrt(g));
			for(k=0;k<=n;k++)
				A[n][k]=g*b[k];
		}


	}

	//Finding Orthonormal System Weights

	for(i=0;i<M;i++)
	{
		for(m=0;m<L;m++)
		{
			Wt[i][m]=0;
			for(k=0;k<=m;k++)
				Wt[i][m]+= A[m][k]*C[i][k];

		}

	}

	//Finding Training Errors
	e=0;
	for(i=0;i<M;i++)
	{
		E[i]=Et[i];
		for(k=0;k<L;k++)
			E[i]=E[i]-(Wt[i][k]*Wt[i][k]);
		e=e +E[i];
	}
	
	//Find output Weights for original System 
	for(i=0;i<M;i++)
	{
		for(k=0;k<L;k++)

		{
			Wout[i][k]=0;
			for(m=k;m<L;m++)
				Wout[i][k]=Wout[i][k] + A[m][k]*Wt[i][m];
		}
	}

	freeMatrix(A, L);
	freeMatrix(Wt, M);
	freeArray(b);
	freeArray(c);
	freeArray(E);

	return(e);
}

double OLS1(const matrix &Wout, int N, int M, int L, double **R, double **C,std::valarray<double> &Et)
{
	int NLin=0 ;
	int n,m,j,i,k,Nv=0;
	double g,b1,b2,c1,g1=0;  
	double *E,*b,*c;
	double **Wt,**A;
	double e;

	A=createMatrix(L, L);
	Wt=createMatrix(M, L);
	b=createArray(L);
	c=createArray(L);
	E = createArray(M);

	//	Initializing the weights and thresholds, Auto and Cross correlation matrices to zero

	for (i = 0; i < L; i++)
	{
		for (j = 0; j < L; j++)
			A[i][j]=0.0;


		for (j = 0; j < M; j++)
			Wt[j][i] = 0.0;
	
	}

	for (i = 0; i < M; i++ )
		E[i] = 0.0;		

	Nv=0;	


	//OLS TYPE 1 CODE(OWO with no reordering).
	g=R[0][0];
	if(g < pow(10.0,-12.0))
	{
		A[0][0]=0;
		NLin=NLin+1;
	}

	else
	{
		g=sqrt(g);
		A[0][0]=(1/g);
	}

	c1=A[0][0]*(R[0][1]);
	b2=1;	
	b1=((-c1)*A[0][0]);
	g=R[1][1]-(c1*c1);

	if(g < pow(10.0,-12.0))
	{
		A[1][0]=0;
		A[1][1]=0;
		NLin=NLin+1;
	}
	else
	{
		g=sqrt(g);
		A[1][0]=(1/g)*b1;
		A[1][1]=(1/g)*b2;
	}

	for(n=2;n<L;n++)
	{
		for(j=0;j<=(n-1);j++)
		{
			c[j]=0;
			for(k=0;k<=j;k++)
				c[j]+=A[j][k]*R[k][n];
		}
		b[n]=1;

		for(j=0;j<=(n-1);j++)//change
		{
			b[j]=0;
			for(k=j;k<=(n-1);k++)
				b[j]=b[j]-c[k]*A[k][j];
		}


		g=R[n][n];
		g1=0;
		for(k=0;k<=n-1;k++)
			g1=g1+(c[k]*c[k]);
		g=g-g1;
                assert(g>0);
		if(g <pow(10.0,-12.0))
		{
			for(k=0;k<=n;k++)		// change n-1 to n : rohit 5/20/2011
				A[n][k]=0;
			NLin=NLin+1;
		}
		else
		{
			g=1/(sqrt(g));
			for(k=0;k<=n;k++)
				A[n][k]=g*b[k];
		}


	}

	//Finding Orthonormal System Weights

	for(i=0;i<M;i++)
	{
		for(m=0;m<L;m++)
		{
			Wt[i][m]=0;
			for(k=0;k<=m;k++)
				Wt[i][m]+= A[m][k]*C[i][k];

		}

	}

	//Finding Training Errors
	e=0;
	for(i=0;i<M;i++)
	{
		E[i]=Et[i];
		for(k=0;k<L;k++)
			E[i]=E[i]-(Wt[i][k]*Wt[i][k]);
		e=e +E[i];
	}
	
	//Find output Weights for original System 
	for(i=0;i<M;i++)
	{
		for(k=0;k<L;k++)

		{
			Wout[i][k]=0;
			for(m=k;m<L;m++)
				Wout[i][k]=Wout[i][k] + A[m][k]*Wt[i][m];
		}
	}

	freeMatrix(A, L);
	freeMatrix(Wt, M);
	freeArray(b);
	freeArray(c);
	freeArray(E);

	return(e);
}

double innerProd(double **a, double **b, int rows, int cols)
{
	double prod=0.0;
	int i, j;
	for(i=0; i<rows; i++)
	{
		for(j=0; j<cols; j++)
		{
			prod+=a[i][j]*b[i][j];
		}
	}
	return prod;
}

double* createArray(int size)
{
    double *ptr = (double *)malloc(sizeof(double) * size);
    if(!ptr)
    {
        printf("Memory allocation error!\n");
        exit(0);
    }
    return ptr;
}

double** createMatrix(int rows, int cols)
{
    int i;
    double **mat;
    
    mat = (double **)malloc(sizeof(double*) * rows);
    if(mat==NULL)
    {
        printf("Memory allocation error.\n");
        exit(0);
    }
    for(i=0; i<rows; i++)
    {
        mat[i] = (double *)malloc(sizeof(double) * cols);
        if(mat[i]==NULL)
        {
            printf("Memory allocation error.\n");
            exit(0);
        }
    }

	return mat;
}

void freeArray(double *ptr)
{
	assert(ptr!=NULL);
	free(ptr);
}

void freeMatrix(double **ptr, int rows)
{
	int i;
	assert(ptr!=NULL);
	for(i=0; i<rows; i++)
	{
		assert(ptr[i]!=NULL);
		free(ptr[i]);
	}
	free(ptr);
}

void safe_gets(char *str, int size)
{
	int i;
	fgets(str, size-1, stdin);
	i = strlen(str)-1;
	if(str[i] == '\n') 
		str[i] = '\0';
}

double Schmidt(double **R, double **C, double **W, valarray<double> &outputEnergies, int Nu, int M)
{
	int i, j, k, m, q;
	double g=0.0;
	double tempc = 0.0, eps = 1E-5, E, temp;

	double **schmidt_A = Matr(Nu, Nu);
	double **schmidt_Woth = Matr(M, Nu);
	double *schmidt_c = Arr(Nu);
	double *schmidt_b = Arr(Nu);

	schmidt_A[0][0] = 1/sqrt(R[0][0]);

	for(m = 1; m < Nu; m++)
	{
		for(i = 0; i < m; i++)
		{
			schmidt_c[i] = 0.0;
			for( q = 0; q <= i; q++)
				schmidt_c[i] += schmidt_A[i][q] * R[q][m];
		}

		for(k = 0; k < m; k++)
		{
			schmidt_b[k] = 0.0;
			for(i = k; i <= m-1; i++)
				schmidt_b[k] -= schmidt_c[i] * schmidt_A[i][k];
		}

		schmidt_b[m] = 1.0;

		tempc = 0.0;

		for(i = 0; i < m; i++)
			tempc += pow(schmidt_c[i], 2);

		g = R[m][m] - tempc;

		if(g < eps)
		{
			for(k = 0; k <= m; k++)
				schmidt_A[m][k] = 0.0;
		}
		else
		{
			for(k = 0; k <= m; k++)
				schmidt_A[m][k] = schmidt_b[k] / sqrt(g);
		}
	}

	for(i = 0; i < M; i++)
	{
		for(q = 0; q < Nu; q++)
		{
			schmidt_Woth[i][q] = 0.0;
			for(k = 0; k <= q; k++)
				schmidt_Woth[i][q] += schmidt_A[q][k] * C[i][k];
		}
	}

	for(i = 0; i < M; i++)
	{
		for(q = 0; q < Nu; q++)
		{
			W[i][q] = 0.0;
			for(k = q; k < Nu; k++)
				W[i][q] += schmidt_Woth[i][k] * schmidt_A[k][q];
		}
	}
	
	E = 0.0;
	for(k=0; k<M; k++)
	{
		for(i=0; i<Nu; i++)
		{
				temp = 0.0;
				for(j=0; j<Nu; j++)
				{
					temp += W[k][j]*R[j][i];
				}
				E += (-2.0*C[k][i]*W[k][i] + W[k][i]*temp);
		}

		E += outputEnergies[k];
	}

	return (E);
}

void variance(int N, int M, const char *infile, valarray<double> x, valarray<double> t, valarray<double>xm,valarray<double>xv)
{
    int i, k;
    FILE *ifs2;
    int Nv=0;
    ifs2=fopen(infile,"r");
	for(i=0;i<N;i++)
         {
              xm[i]=0.0;
			  xv[i]=0.0;
         }
    while(!feof(ifs2))
    {
         for(i=0;i<N;i++)
         {
              fscanf(ifs2,"%lf", &x[i]);
              xm[i]+=x[i];
         }         
         for(k=0;k<M;k++)
         fscanf(ifs2,"%lf", &t[k]);
         if(feof(ifs2)) break;
         Nv++;
    }
    rewind(ifs2);

    for(i=0;i<N;i++)
    {
         xm[i]/=Nv;
    }        
    Nv=0;

    while(!feof(ifs2))
    {
         for(i=0;i<N;i++)
         {
              fscanf(ifs2,"%lf", &x[i]);
              xv[i]+=(x[i]-xm[i])*(x[i]-xm[i]);

              x[i]=x[i]-xm[i];
         }         
         for(k=0;k<M;k++)
         fscanf(ifs2,"%lf", &t[k]);
         if(feof(ifs2)) break;
         Nv++;
    }
    rewind(ifs2);

    for(i=0;i<N;i++)
    {
         xv[i]/=Nv;
    }
    fclose(ifs2);
}

double hnet(int k, int N, int Nh, matrix wih,valarray<double> x)
	{
        int in;
		double h1=0.0;
		for(in=0;in<N;++in)
		{
			h1+=x[in]*wih[k][in];
		}
		h1=h1+wih[k][N];
		return(h1);
	}

double hact(double h1)
   {
         double h2=0.0;
         h2=1/(1+exp(-h1));
         return (h2);
   }

double onet(int i, int N, int Nh, int M, double **wio, matrix wih, double **who, double *th, double *to, double *h, valarray<double> x)
	{
        int in,hu;
		double o1=0;
		for(in=0;in<N;++in)
		{
			o1+=x[in]*wio[i][in];
		}
		for(hu=0;in<Nh;++hu)
		{
			h[hu]=hnet(hu,N, Nh, wih, x);
			o1+=h[hu]*who[i][hu];
		}
		o1=o1+to[i];
		return(o1);
	}

double CGrad(int Nh, int Nu, int N,int M, valarray<double> t1,double **w,matrix & who,matrix R,matrix C,matrix &woo)
{

//int d=0;
int u=0;
int NLin=0,i,m,j,n,k;
double e,e2,c1=0,C1=0,b2=0,b1=0,tt=0;
double **a,*c2,*b;
double *Error;
double g=0;
//long double **a;
//long double **c2;
double MSE=0.0;
a=(double **)malloc(Nu*sizeof(double *));
for(u=0;u<Nu;++u)
  a[u]=(double *)malloc(Nu*sizeof(double));
c2=(double *)malloc(Nu*sizeof(double));
b=(double *)malloc(Nu*sizeof(double));
//double tt=0;
//long double woo[7][45];
Error = (double *) malloc(sizeof(double)*(M));
e=0.000001;
e2=e*e;
g=R[0][0];

for(i=0;i<Nu;++i)
{  
for(j=0;j<Nu;++j)
  {

	  a[i][j]=0.0;
      
  }

}
for(j=0;j<Nu;++j)
{
  b[j]=0.0;
c2[j]=0.0;
}


	//b=(long double *)malloc(Nu*sizeof(long double *));
//for(j=0;j<1;++j)
		//b[j]=(double *)malloc(Nu*sizeof(double));


if(g<e2)
{
a[0][0]=0.0;
NLin=NLin+1;
}
else
{
g=sqrt(g);
a[0][0]=1/g;
}
tt=R[0][1];
c1=a[0][0]*R[0][1];
b2=1;
b1=-c1*a[0][0];
g=R[1][1]-(c1*c1);

if(g<e2)
{
a[1][0]=0.0;
a[1][1]=0.0;
NLin=NLin+1;
}
else
{
g=pow(g,0.5);
a[1][0]=(1/g)*b1;
a[1][1]=(1/g)*b2;
}


for(n=2;n<Nu;n++)
{
for(j=0;j<=n-1;j++)
{
c2[j]=0.0;
for(k=0;k<=j;k++)
{
tt=R[k][n];
	c2[j]=c2[j]+a[j][k]*R[k][n];
}
}
b[n]=1.0;
for(j=0;j<=n-1;j++)
{
b[j]=0.0;
for(k=j;k<=n-1;k++)
{
b[j]=b[j]-c2[k]*a[k][j];
}
}
C1=0.0;
for(k=0;k<=n-1;k++)
{
C1=C1+c2[k]*c2[k];
}
g=R[n][n]-C1;

if(g<e2)
{
for(k=0;k<=n;k++)
{
a[n][k]=0.0;
}
NLin=NLin+1;
}
else
{
g=1/sqrt(g);
for(k=0;k<=n;k++)
{
a[n][k]=g*b[k];
}
}   
}

for(i=0;i<M;i++)
{	
for(m=0;m<Nu;m++)
	{
		woo[i][m]=0.0;
		for(k=0;k<=m;k++)
		{
		tt=C[i][k];
			//woo[i][m]=woo[i][m]+C[i][k]*a[m][k];
        woo[i][m]=woo[i][m]+C[i][k]*a[m][k];
		}
	}

}

for (i=0;i<M;i++)
      {
		  Error[i]=0.0;
	  }

for(i=0;i<M;i++)
{
	Error[i]=t1[i];
}
for(i=0;i<M;i++)
{
	for(m=0;m<Nu;m++)
	{
		Error[i]=Error[i]-(woo[i][m]*woo[i][m]);
	}
}
	for (i=0;i<M;i++) 
    {
		MSE+=Error[i];
		//printf("\nError at node %d : %lf",(i+1),Error[i]/Nv);
	}

	//getch();
	//return(MSE);
for(i=0;i<M;i++)
{
for (k=0;k<Nu;k++)
{ 
   w[i][k]=0;
   for(m=k;m<Nu;m++)
	  { 
	   w[i][k]=w[i][k]+a[m][k]*woo[i][m];
   }
}
}

for(k=0;k<M;k++)
{
	for(j=0;j<Nh;j++)
	{
	who[k][j]=w[k][j+N+1];
}
}

   return(MSE);
}     
