/******************************************************************************

                     Split a data file to two new files

*******************************************************************************/
# include <stdlib.h>
# include <string.h>
# include <stdio.h>
# include <math.h>

void main(int argc,char *argv[])
{
   int IX=3;
   int IY=4009;
   int IZ=234;

   char  Infile[60];
   FILE *infile,*outfile1,*outfile2;
   char data[100];//, outfile1[100],outfile2[100];
   float p,pb;          //  p: The probability of a pattern belonging to file new1.txt;
   int i,j,N,M;             //  The column numbers for each pattern
   float *x,*t;
	float rand1(int *, int *, int *);

	if(argc != 7) 
	{
		printf("Usage: split training_file output_file_1 output_file_2 proportion num_inputs num_outputs\n");
		exit(1);
	}
	p=atof(argv[4]);
    N=atoi(argv[5]);
    M=atoi(argv[6]);
	
	x= (float *)malloc(sizeof(float)*N);
    t= (float *)malloc(sizeof(float)*M);
	/*  Open the files to read or write data  */
     infile=fopen(argv[1],"r");
	 if (infile==NULL)
	 {
		 printf("cannot open training file!\n");
		 exit(1);
	 }

	 outfile1=fopen(argv[2],"w");
	 if (outfile1==NULL){
		 printf("cannot create %s ! \n",argv[2]);
		 exit(1);
	 }
       
     outfile2=fopen(argv[3],"w");
	 if (outfile2==NULL)
         {
		 printf("cannot open %s ! \n",argv[3]);
		 exit(1);
	 }

	 while(!feof(infile))
	 {
         pb=rand1(&IX,&IY,&IZ);
         for(i=0;i<N && !feof(infile);i++) 
         {
	         fscanf(infile,"%f",&x[i]);
		 if (feof(infile))
                      break;	
		     if(pb<=p) 
		        fprintf(outfile1," %f ",x[i]);
		     else
		        fprintf(outfile2," %f ",x[i]); 
         }

		 for (i=0;i<M && !feof(infile);i++)
		 {
		    fscanf(infile,"%f%",&t[i]);
                      if (feof(infile))
                         break;
		    if(pb<=p) 
			  fprintf(outfile1," %f ",t[i]);
		    else 
			 fprintf(outfile2," %f ",t[i]);
		 }
		 if(pb<=p) 
			  fprintf(outfile1,"\n",t[i]);
	     else 
			 fprintf(outfile2,"\n",t[i]);
		 
	 }                /*  End of while loop  */

	 fclose(infile);
	 fclose(outfile1);
	 fclose(outfile2);

}                      /* End of main  */
     

        
/******************************************************************************

	Subroutine : Random no. Generator
		-the random nos. generated are uniformly distributed
		 between 0 and 1.
		-IX, IY, IZ are the SEEDS

*******************************************************************************/


float rand1(int *ix, int *iy, int *iz)
{
	int ixx, iyy, izz;
	float itemp;
	float temp;


	ixx=(*ix)/177;
	*ix=171*(*ix%177)-2*ixx;

	if(*ix < 0)
		*ix+=30269;

	iyy=(*iy)/176;
	*iy=176*(*iy%176)-2*iyy;

	if(*iy < 0)
		*iy+=30307;

	izz=(*iz)/178;
	*iz=170*(*iz%178)-2*izz;

	if(*iz < 0)
		*iz+=30323;

	temp=(float)(*ix)/30629.0+(float)(*iy)/30307.0+(float)(*iz)/30323.0;
	itemp=floor(temp);
	return (temp-itemp);

}

