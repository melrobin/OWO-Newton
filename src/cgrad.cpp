#include <stdlib.h>
#include <stdio.h>
#include "fortran.h"
/*******************************************************************************
Title    : Subroutine- Conjugate gradient
Purpose: Minimize a quadratic of the form f(w)=w^T*R*w-C^T*w-b
Comments : Follow the derivation done in the class for better understanding of
			  this routine
*******************************************************************************/

double CGrad(int Nu,double *W,double **R,double *C,double E)
{
	double XD, XN, Num, Den, B1, B2, Ei = 0.0;
	double *p, *g,*tempg1;
	int iter;
	int intone=1;
	double one=1.0,zero=0.0,two=2.0;
	double minusone=-1.0,minustwo=-2.0;
	FILE *outfile;

	if ((outfile=fopen("error.txt","a"))==NULL)
	   printf("Could not open up output file\n");
	XD = 1;
/* Initialize the direction and the gradient vectors to zero */
	p = (double *)calloc(Nu,sizeof(double));
	g = (double *)calloc(Nu,sizeof(double));
	tempg1=(double *)calloc(Nu,sizeof(double));

	for(iter = 0; iter < Nu; iter++)	  					/* start of Iteration loop */
	{
	    dgemv_("N",&Nu,&Nu,&two,&R[0][0],&Nu,W,&intone,&zero,g,&intone);
	    daxpy_(&Nu,&minustwo,C,&intone,g,&intone);

	    XN=ddot_(&Nu,g,&intone,g,&intone);

	    B1 = XN/XD;
	    XD = XN;

	    dscal_(&Nu,&B1,p,&intone);
	    daxpy_(&Nu,&minusone,g,&intone,p,&intone);

     	    Den = Num = 0.0;

	    dgemv_("N",&Nu,&Nu,&one,&R[0][0],&Nu,p,&intone,
		&zero,tempg1,&intone);
	    Den=ddot_(&Nu,tempg1,&intone,p,&intone);
	    Num=-ddot_(&Nu,p,&intone,g,&intone)/2.0;
			
	    B2 = Num/Den;

	    daxpy_(&Nu,&B2,p,&intone,W,&intone);

	   dgemv_("N",&Nu,&Nu,&one,&R[0][0],&Nu,W,&intone,&zero,tempg1,&intone);
	    Ei=E+ddot_(&Nu,tempg1,&intone,W,&intone)
		-2*ddot_(&Nu,C,&intone,W,&intone);
	    fprintf(outfile,"Iteration:  %d\tError %f\tGradient: %f\tB2: %f\n",
		iter+1,Ei,dnrm2_(&Nu,g,&intone),B2);

	}              							/* end of Iteration loop */
        fprintf(outfile,"\n");
/*Here we calculate the error by taking advantage of the fact that we
  already have solved for W and the error will simply be the residual:
  W^T*R*W-2*C */
	dgemv_("N",&Nu,&Nu,&one,&R[0][0],&Nu,W,&intone,&zero,tempg1,&intone);
	Ei=E+ddot_(&Nu,tempg1,&intone,W,&intone)-2*ddot_(&Nu,C,&intone,W,&intone);

	free(p);
	free(g);
	free(tempg1);
	fclose(outfile);
	return (Ei);
}      
