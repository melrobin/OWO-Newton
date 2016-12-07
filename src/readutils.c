#include "utils.h"
#include <stdlib.h>
#include <math.h>
#include <stdio.h>
/***************************************************************************
Function ReadOnePattern
***************************************************************************/
int	ReadOnePattern(struct UserEntry *userptr,struct MiscArrays *arrayptr)
{
	char str[FILENAMELEN];
	int i,j, Nf, Nmax, Nout, retval;
	double  *Inputs, *Outputs,  **Hout;
	FILE *fp;


	Nf = userptr->Features;
	Nout = userptr->DesiredOutputs;
	Nmax = Nf+Nout;
	fp = userptr->fpfeat;
	Inputs = arrayptr->Inputs;
	Outputs = arrayptr->Outputs;
	Hout = arrayptr->LayerOutputs;

	if (!feof(fp))
	{
		retval = 1;
		for ( i = 0; i < Nmax; i++)
		{
			if (fscanf(fp, "%s",str ))
			{
				if (feof(fp))
				{
					retval = EOF;
					break;
				}
				if ( i < Nf)
				{
					Inputs[i] = atof(str);
					Hout[0][i] = Inputs[i];
				}
				else
				{
					j = i-Nf;
					Outputs[j] = atof(str);
				}
			}
			else
			{
				break;
			}
		}
	}
	else
	{
		retval = EOF;
	}
	return(retval);
}
/***************************************************************************/
void	 CalculateNet(struct TopologyInfo *tptr,struct MiscArrays *arrayptr)
{
	int i,j,k,layernum,  TotalLayers,  Lhid, *LUnits, OutLayer;
	double net,   **Nets, **Thresholds, **Hout;

	LUnits = tptr->LayerUnits;
	Lhid = tptr->HiddenLayers;
	OutLayer = Lhid+1;
	TotalLayers = Lhid+2;

	Nets = arrayptr->Nets;
	Hout = arrayptr->LayerOutputs;
	Thresholds = arrayptr->Thresholds;



	/* need to reinitialize Nets[i][j] = 0.0 for each pattern	*/

	for ( i= 1; i < TotalLayers; i++)
	{
		for ( j= 0; j < LUnits[i]; j++)
		{
			Nets[i][j] = 0.0;
		}
	}

	for ( i= 1; i < TotalLayers; i++)
	{
		for ( j= 0; j < LUnits[i]; j++)
		{
			Nets[i][j] += Thresholds[i][j];

			for ( layernum = 0; layernum < i; layernum++)
			{
				for ( k = 0; k < LUnits[layernum]; k++)
				{
					Nets[i][j] += arrayptr->Wts[i][layernum][j][k]*Hout[layernum][k];
				}
			}
			/* once net function for each unit has been calculated, f(net) = output
			for that unit should be calculated and written back to the array
			to be ready for the next iteration*/
			net =  Nets[i][j] ;
			/* output layer has  linear activation	*/
			if ( i == OutLayer)
			{
				Hout[i][j] = net;
			}
			else
			{
				/* to take care of underflow, overflow problems	*/
				if (net >= 25.0)
				{
					Hout[i][j] = 1.0;
				}
				else if (net <= -25.0)
				{
					Hout[i][j] = 0.0;
				}
				else
				{
					Hout[i][j] = (1.0 /(1.0 + exp(-net)));
				}
			}
		}
	}
}
/*****************************************************************************
Function CalculateError calculates error in the k th output
***************************************************************************/
double CalculateError(struct UserEntry *userptr,struct TopologyInfo *tptr,
							struct MiscArrays *arrayptr)
{
	int 		Nout,i,Lhid,	OutLayer;
	double   *Outputs,  **Hout, error, TotalMse = 0.0;



	Nout = userptr->DesiredOutputs;
	Lhid = tptr->HiddenLayers;

	Hout = arrayptr->LayerOutputs;


	Outputs = arrayptr->Outputs;

	OutLayer = Lhid+1;

	for(i =0; i < Nout; i++)
	{
		error = Outputs[i] - Hout[OutLayer][i];
		arrayptr->MSError[i] += error*error;
		if (userptr->SaveOutput)
		{
			fprintf(userptr->fpout,"%lf ",Hout[OutLayer][i]);
		}
	}
	if (userptr->SaveOutput)
	{
		fprintf(userptr->fpout,"\n ");
	}

	for(i = 0; i < Nout; i++)
	{
		TotalMse += arrayptr->MSError[i];
	}

	return(TotalMse);
}
void make2Dfrom1D(double *in,double **out,int M,int N)
{
	int i,j;

	for (i=0;i<M;i++)
		for (j=0;j<N;j++)
			out[i][j]=in[i+j*M];
}
void flatten(int N,double **in,double *out)
{
	int i,j;

	for (i=0;i<N;i++)
		for (j=0;j<N;j++)
			out[i*N+j]=in[i][j];
}
void dump_matrix(int M,int N,const char *fname,double **mat)
/* Dumps a matrix to a file */
{
	int i,j;
	FILE *fp;

	if ((fp=fopen(fname,"w"))==NULL)
	   printf("Could not open file...\n");
	for (i=0;i<M;i++)
	{
		for (j=0;j<N;j++)
			fprintf(fp,"%f ",mat[i][j]);
		fprintf(fp,"\n");
	}
	fclose(fp);
}

void dump_vector(int N,char *fname,double *mat)
/* Dumps a matrix to a file */
{
        int i;
        FILE *fp;

        if ((fp=fopen(fname,"a"))==NULL)
           printf("Could not open file...\n");
        for (i=0;i<N;i++)
          fprintf(fp,"%f ",mat[i]);
	fprintf(fp,"\n");
        fclose(fp);
}


