/******************************************************************************
MODIFIED ON :	24 Nov 98 belongs to LM for Web project

	  major modification done on 2 oct 98

	  aim : out of the two iterations, return the solution from that
	  iteration that gives lowest error

	  complication : same function is called to solve for output weights
	  as well as hidden weights... there are different auto, cross
	  correlation matrices involved

	  knowing which ones to use complicates error comparison




PROGRAM		:	conjugate.C
AUTHOR		:	Hema Chandrasekaran
CONJUGATE.C	:      Belongs to \LM project
SYNOPSIS 	:

INPUTS		:

OUTPUTS		:
*******************************************************************************/
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include "conjugat.h"
extern struct  UserEntry   *userptr;
extern struct  TopologyInfo *tptr;
extern struct  MiscArrays  *arrayptr;
extern double **R,**C;
/*****************************************************************************/
void ConjugateGradient(int InputDim)
{
//	double		**R, *T;
double *T;
double	*g, *p,  *x, *x2;
 long double B1, B2,  b1numerator, b1denominator = 1.0,wts_cor = 0.0,
 p_cor, b2numerator, b2denominator,Error1 = 0.0, Error2 = 0.0,
 BestError = 0.0, eps = 1e-20;

 int i, j, k,  N, OuterIteration;

N =  InputDim;
g = arrayptr->g;
p = arrayptr->p;
x = arrayptr->ColumnWts;
x2 = arrayptr->ColumnWts2;

T = arrayptr->ColumnCrossCor ;


	/* 1 Oct 98	*/
	/* R is the Auto Correlation matrix */
	/* T is one column of Cross Correlation Matrix */


//R = arrayptr->GlobalAutoCor;

	/* initialize the direction vector and the weights	*/
	/* iteration is one more than the input vector dimension	*/

	/* not needed	at all */

	/* OuterIteration  is needed to improve accuracy of the solution;
	particularly valuable in numerically ill-conditioned problems	*/

	for(OuterIteration =1;OuterIteration <=2;OuterIteration++)
	{
		b1denominator = 1.0;


/* commented out on 27 sept 98 ..  solution   improved*/

		for(i=0; i < N; i++)
			p[i] = 0.0;

		for ( i = 0; i <= N; i++)
		{
			/* index i is the conjugate gradient iteration index ;
			loop terminates after (N+1) iterations	*/
			/*-----------------------------------------------------------------*/
						/*	CALCULATION OF B1	*/
			/*-----------------------------------------------------------------*/
			b1numerator = 0.0;
			for ( j = 0; j < N; j++)
			/* index j is the index of  g(j)		*/
			{
				/* reset wts_cor to 0.0 for the next g(j)	*/
				wts_cor = 0.0;
				for ( k = 0; k < N; k++)
				wts_cor += x[k] * R[k][j];
				g[j] =  wts_cor - T[j] ;

				b1numerator += g[j] * g[j];
			}

			B1 = b1numerator / (b1denominator+eps);

			/* for the next iteration denominator becomes the numerator	*/
			b1denominator = b1numerator;
			for ( j = 0; j < N; j++)
				p[j] = -g[j] + B1 * p[j];

			/*-----------------------------------------------------------------*/
						/*	CALCULATION OF B2	*/
			/*-----------------------------------------------------------------*/

			b2numerator = 0.0;
			b2denominator	= 0.0;

			for ( j = 0; j < N; j++)
			/* index j is the index of  p(j)		*/
			{
			p_cor	 = 0.0;
			for ( k = 0; k < N; k++)
			p_cor+= p[k] * R[j][k];
			b2numerator += p[j] * g[j];
			b2denominator += p[j] * p_cor;
			}
			B2 = -b2numerator / (b2denominator+eps);
			for ( k = 0; k < N; k++)
				x[k] = x[k] + B2 * p[k];
		}

		/* calculate error	and keep track of the best solution */

	/* code from the modular net program.... adapted to the need to pick
the best solution given by one of the two iterations of conjugate gradient
	algorithm
	for derivation refer to notes under modular net program as to how to
	calculate error efficiently */

		if (OuterIteration == 1)
		{
			for(j =0; j < N; j++)
			{
				Error1 -= 2.0 * x[j] * T[j];
				for(k =0; k < N; k++)
					Error1 += x[j]* x[k]* R[j][k];
			}
			BestError = Error1;

			/* save the column wts, since it will be overwritten */
			for(j =0; j < N; j++)
				x2[j] = x[j];
		}
		else
		{
			Error2 = 0.0;
			for(j =0; j < N; j++)
			{
				Error2 -= 2.0 * x[j] * T[j];
				for(k =0; k < N; k++)
					Error2 += x[j]* x[k]* R[j][k];
			}
			if (Error2 > Error1)
			{
				/* restore columnwts from first iteration */
				/* x[j] contains the current solution which is not optimal
					if Error2 > Error1	*/
				for(j =0; j < N; j++)
					x[j] = x2[j];
				BestError = Error1;
/*				printf(" Error increased in Conjugate gradient second iteration \n");


	printf("\nconjugate gradient \n");
	printf("\nError1 = %Lf, Error2 = %Lf\n", Error1, Error2);

*	fprintf(userptr->fplog,"\nconjugate gradient \n");
	fprintf(userptr->fplog,"\nError1 = %lf, Error2 = %lf\n", Error1, Error2);

*/

			}
			else
				BestError = Error2;
		}
	}


		/* modified on 25 Sept 98 */
	if (arrayptr->HwoOrOwo == 1)
		/* calculate error	*/
		arrayptr->ErrorCG = BestError;
}
/*****************************************************************************/

