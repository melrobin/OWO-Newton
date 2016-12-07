#include <stdlib.h>
#include <stdio.h>
#include "allocmem.h"

double  **FarAllocateMatrixMemory(int , int );

/***************************************************************************/
FILE *OpenFile(char *filename, char *mode)
{
	FILE *fp;
	/* no user interaction here ! if unable to open , return NULL pointer
		announcing error; if successful, seek beginning of opened file
		indicating success	*/
	if((fp = fopen(filename,mode))== NULL)
	{
		printf("%s file could not be opened \n", filename);
		return(NULL);
	}
	/* if successful, position the file pointer at the beginning	*/
	fseek(fp, (long)0, 0);
	return(fp);
}
/***************************************************************************/
double *FarAllocateMemory(int size)
{
	double *dptr;

	dptr = (double *)calloc(size, sizeof(double));
	if(dptr == NULL)
	{
		printf("Not enough memory to allocate buffer. \n");
		printf("Aborting......\n");
		/* terminate program if out of memory 	*/
		exit(1);
	}
	return(dptr);
}
/*****************************************************************************/
double *AllocateMemory(int size)
{
	double  *dptr;

	dptr = (double *)calloc(size, sizeof(double));
	if(dptr == NULL)
	{
		printf("Not enough memory to allocate buffer. \n");
		printf("Aborting......\n");
		/* terminate program if out of memory 	*/
		exit(1);
	}
	return(dptr);
}

/*****************************************************************************/
double  **FarAllocateMatrixMemory(int m, int n)
{
	double  **dptr;
	int i;

	dptr = (double **)calloc(m, sizeof(double *));
	if(dptr == NULL)
	{
		printf("Not enough memory to allocate buffer. \n");
		printf("Aborting......\n");
		/* terminate program if out of memory 	*/
		exit(1);
	}
	dptr[0]=(double *)calloc(m*n,sizeof(double));
	for(i = 1; i < m ; i++)
	{
		dptr[i]=dptr[0]+i*n;
		if(dptr[i] == NULL)
		{
			printf("Not enough memory to allocate buffer. \n");
			printf("Aborting......\n");
			/* terminate program if out of memory 	*/
			exit(1);
		}
	}
	return(dptr);
}
/*****************************************************************************/

