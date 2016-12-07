/************************************************************************************************************
							  
								
							  
									 Author :  Jiang Li  

										EE Department, 
							Image Processing and Neural Network Lab
								University of Texas at Arlington
										  7/1/2002

		This program is the process program for the MLP Network, it uses the trained network to process 
	the data which may have inputs and outputs or only has inputs. If the data has outputs this program can
	compare the predicted outputs and the actual outputs, if the data has only inputs this program will provide
	the predicted outputs. Our Lab has a standard weights file format, the following is the format details of 
	the weights file:


	Weights File Format:
	

	1. Data File Name
	2. No. of Inputs
	3. No. of Hidden Units
	4. No. of Outputs
	5. Training Algorithm
	6. Output Weights, (Inputs to Outputs, Hidden Units to Outputs), from inputs first, then from hidden units
	7. Hidden weights, From Inputs to Hidden Units.


***************************************************************************************************************/

#include "math.h"
#include "stdio.h"
#include "stdlib.h"
#include "string.h"


#define WITHOUTPUT  1
#define NOOUTPUT    0

#define RESULT   "ProcessingResult.txt"



main()
{
	FILE *fpData, *fpWeights, *fpResult;
	double *x, *new_xa, **Wo, **W,*ty,*y, *Error, TotalError,*net;
	int N, Nout, Nh,Nv, flagoutput,pattern;
	int i, j,k;
	char str[30], *FileName, *WeightsFileName,
		TrainingAlgorithm[200];
	
	char *get_string(char *);                                    /* Input string */
	int get_int(char *,int, int);                                /* Input int varible */

	/* Get the processing file Name */
	do
	{
		FileName = get_string("Enter testing or processign file name : ");
		fpData = fopen(FileName,"r");
	}
	while(fpData == NULL);
	
	/* does the testing file has a disired output? */
	printf("\n Does the testing file has desired output? \n");
	printf("\n Choose (0) for NO \n");
	printf("\n        (1) for YES : ");
	flagoutput = get_int("",0,1);
	
	/* Get the Weights file Name */
	do
	{
		WeightsFileName = get_string("Enter weights file name : ");
		fpWeights = fopen(WeightsFileName,"r");
	}
	while(fpWeights == NULL);
	

	fscanf(fpWeights,"%s",FileName);
	fscanf(fpWeights,"%d",&N);
	fscanf(fpWeights,"%d",&Nout);
	fscanf(fpWeights,"%d",&Nh);
	fscanf(fpWeights,"%s",TrainingAlgorithm);
	//initilinize the weights matrixs
	x = (double *)malloc(sizeof(double)*(N+1));
	new_xa = (double *)malloc(sizeof(double)*(N+Nh+1));
	net = (double *)malloc(sizeof(double)*Nh);
	ty = (double *)malloc(sizeof(double)*Nout);//actual output
	y = (double *)malloc(sizeof(double)*Nout);//predict output
	Error = (double *)malloc(sizeof(double)*Nout);
	Wo = (double **)malloc(sizeof(double *)*Nout);
	W = (double**)malloc(sizeof(double *)*Nh);
	
	for(i = 0 ; i < Nout; i ++)
	{
		Wo[i] = (double *)malloc(sizeof(double)*(N+Nh+1));
	}
	for(i = 0; i < Nh; i ++)
	{
		W[i] = (double *)malloc(sizeof(double)*(N+1));
	}
	//read the outputs weights
	for(i = 0 ; i < N + Nh + 1; i ++)
	{
		for( j = 0 ; j < Nout; j ++)
		{
			fscanf(fpWeights,"%s",str);
			Wo[j][i] = atof(str);
		}
	}
	//read the hidden units weights

	for(i = 0; i < Nh; i ++)
	{
		for( j = 0; j < N + 1; j ++)
		{
			fscanf(fpWeights,"%s",str);
			W[i][j] = atof(str);
		}
	}

	
	if((fpResult = fopen(RESULT,"w")) == NULL)
	{
		printf("Can't open result.txt to write!\n");
		exit(1);
	}

	fprintf(fpResult,"\n\n\n\t\tMLP  Process  Program");
	fprintf(fpResult,"\n\n\tTo be processed File Name:%s\t\t\t",FileName);
	fprintf(fpResult,"\n\n\tWeights File Name:\t%s",WeightsFileName);
	fprintf(fpResult,"\n\n\tNo of Inputs:\t%d",N);
	fprintf(fpResult,"\n\n\tNo of Hidden Units: \t%d",Nh);
	fprintf(fpResult,"\n\n\tNo of Outputs:\t%d",Nout);
	fprintf(fpResult,"\n\n\t%s.",TrainingAlgorithm);
	if(flagoutput == WITHOUTPUT)
	{
		fprintf(fpResult,"\n\n\tThis data has desired output.\n");
	}
	else
	{
		fprintf(fpResult,"\n\n\tThis data do not has desired output.\n");
	}
	
	fprintf(fpResult,"\n\n\tThe following is the format of the result:\n\tPatterns, Input(s), Predict Output(s), Actual Output(s) if it has. \n\n");


	Nv = 0;
	for(i = 0; i < Nout; i ++)
	{
		Error[i] = 0.0;
	}
	TotalError = 0.0;
	while(!feof(fpData))
	{
		Nv ++;
		for(i = 0; i < N ; i++)
		{
			fscanf(fpData,"%s",str);
			x[i] = atof(str);
			new_xa[i] = x[i];
		}
		if(flagoutput == WITHOUTPUT)
		{
			for(i = 0; i < Nout; i ++)
			{
				fscanf(fpData,"%s",str);
				ty[i] = atof(str);
			}
		}
		new_xa[N] = 1.0;
		
	}
	Nv--;
	rewind(fpData);
	pattern = 0;
	for(k = 0; k < Nv ; k ++)
	{
		pattern++;
		
		for(i = 0; i < N ; i++)
		{
			fscanf(fpData,"%s",str);
			x[i] = atof(str);
			new_xa[i] = x[i];
		}
		if(flagoutput == WITHOUTPUT)
		{
			for(i = 0; i < Nout; i ++)
			{
				fscanf(fpData,"%s",str);
				ty[i] = atof(str);
			}
		}
		new_xa[N] = 1.0;

		
		for(i = 0 ; i < Nh; i ++)
		{
			net[i] = 0.0;
			for(j = 0; j < N+1; j ++)
			{
				net[i] += W[i][j]*new_xa[j];
			}
			new_xa[i+N+1] = 1/(1+exp(-net[i]));
		}

		for(i = 0; i < Nout; i ++)
		{
			y[i] = 0.0;
			for(j = 0;j < N+Nh+1; j ++)
			{
				y[i] += Wo[i][j]*new_xa[j];
			}
			Error[i] += (ty[i]-y[i])*(ty[i]-y[i]);
		}
		fprintf(fpResult,"\t%d",pattern);
		
		for(i = 0 ; i < N ; i ++)
		{
			fprintf(fpResult,"\t%f",x[i]);
		}

		for( i = 0; i < Nout; i ++)
		{
			fprintf(fpResult,"\t%f",y[i]);
		}
		if(flagoutput == WITHOUTPUT)
		{
			for(i = 0; i < Nout; i ++)
				fprintf(fpResult,"\t%f",ty[i]);
		}
		fprintf(fpResult,"\n");
		
		
	}

	

	for(i = 0; i < Nout; i ++)
	{
//		Error[i] /= Nv;
		TotalError += Error[i]/Nv;
	}
	//output to the file
	for(i = 0; i < Nout; i ++)
	{
		printf("Error at Node %d :\t%f\n",i,Error[i]/Nv);
	}
	printf("The Total Error is %f:\t",TotalError);
	
	
	for(i = 0; i < Nout; i ++)
	{
		fprintf(fpResult,"\n\tThe Error at Node[%d] is: \t%f",i,Error[i]/Nv);
	}
	fprintf(fpResult,"\n\n\tThe total Error is: \t%f",TotalError);

	
	fcloseall();

}

char* GetRidofQuotation(char* str)
{
	int i, j;
	char result[100];

	i = 0;
	j = 0;;

	
	do
	{
		if(str[i] == '"')
		{
			i ++;
			continue;
		}
		else
		{
			result[j] = str[i];
			i ++;
			j ++;
		}
	}while(str[i] != 0);
	result[j] = 0;
	str = result;
	return str;
	
	
}

/*************************************************************************/

int get_int(char *title_string,int low_limit, int up_limit)
{
	 int i,error_flag;
	 char *get_string();             /* get string routine */
	 char *cp,*endcp;                /* char pointer */
	 char *stemp;                    /* temp string */

/* check for limit error, low may equal high but not greater */
	 if(low_limit > up_limit) {
		  printf("\nLimit error, lower > upper\n");
		  exit(1);
	 }

/* make prompt string */
	 stemp = (char *) malloc(strlen(title_string) + 60);
	 if(!stemp) {
		  printf("\nString allocation error in get_int\n");
		  exit(1);
	 }
	 sprintf(stemp,"%s [%d...%d]",title_string,low_limit,up_limit);

/* get the string and make sure i is in range and valid */
	 do {
		  cp = get_string(stemp);
		  i = (int) strtol(cp,&endcp,10);
		  error_flag = (cp == endcp) || (*endcp != '\0'); /* detect errors */
		  free(cp);                                   /* free string space */
	 } while(i < low_limit || i > up_limit || error_flag);

/* free temp string and return result */
	 free(stemp);
	 return(i);
}

/*****************************************************************************/
char *get_string(char *title_string)
{
	 char *alpha;                            /* result string pointer */

	 alpha = (char *) malloc(80);
	 if(!alpha) {
		  printf("\nString allocation error in get_string\n");
		  exit(1);
	 }
	 printf(" %s ",title_string);
	 gets(alpha);

	 return(alpha);
}