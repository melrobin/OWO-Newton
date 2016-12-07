// University of Texas at Arlington
// Image Processing and Neural Networks Lab
// Director: Dr M. T. Manry
// Conjugate Gradient with OLF
// Authors: Rohit Rawat, Babu Hemanth Aswathappa, Praveen Jesudhas

// Updated: 07/09/2012 - Rohit Rawat
//(1) Modified main() to accept command line arguments.
//(2) Produces weights file in format compatible with NuMap.

// Updated: 5/15/2011 - Rohit Rawat
//(1) Corrected the OLS routine.
//(2) Replaced gets with a safer version.


// Updated: 5/6/2011 - Rohit Rawat
//(1) Fixed random initialization of input weights.
//(2) Added symmetry to R matrix computation
//(3) Fixed net control and added net control verification
//(4) Initialized saved weights to OLS weights
//(5) Changed type of Ep to double from int
//(6) Changed N+k to N+1+k in delta_po computation
//(7) Added a debug mode, moved net control verification and display
//		of non-essential information to debug mode.
//(8) Routed all memory allocations through wrapper functions
//(9) Added calculation of error on the validation file.

// Updated: 5/20/2011 - Rohit Rawat
//(1) Correction to OLS.

/* DEBUG_MODE: if defined, prints additional debugging information.
	(net control verification during initialization, inner product of
	gradients and direction	vectors (R), sum of gradient energies (Xn),
	number of successive reverts to previous weights (Ierr) during training.

	if commented, the program displays only relevant information.
*/
//#define DEBUG_MODE

// Names of files written to
#define RESULTS_FILE "results.txt"
#define WEIGHTS_FILE "mlp_weights.wts"
#define MSE_FILE "mse_values.txt"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <assert.h>

/* Utility functions */
double* createArray(int size);
double** createMatrix(int rows, int cols);
void freeArray(double *ptr);
void freeMatrix(double **ptr, int rows);
double innerProd(double **a, double **b, int rows, int cols);
void safe_gets(char *str, int size);

int main(int argc, char *argv[])
{
	double slete(double , double );
	double OLS1(double **, int, int, int, double **, double **, double *);
	
	time_t rawtime;
	struct tm * timeinfo;

	char fname_trg[200];
	char fname_tst[200];
	int N, M, Nh, L, Nv;
	FILE *infile, *binfile;
	FILE *result, *mse_file, *weights_file;


	int n, k, i, it, Nit,j;
	int Ierr, reset_Ierr;

	double *x, *xv,*ty, *mm, *vv,*Error, *f, *y, **goi_wts,**goh_wts,**p,**poh,**poi;
	double *mx;
	double **w, *net, *O,*Fgx,*Dyp;
	double *m_hu, *sd_hu, m_d, sd_d,Err,Xd,B1,Xn,zp;

	double **r, **c, *Et;

	double **Wout;

	double **saved_wo, **saved_w, **d;

	double MSE, *E, dez,dez2, Ep;

	double *delta_po, *delta_p;

	double **g_wts;

	double temp, R;

	double z;

	double gradient_energy = 0.0;
	double gradient_energy1 = 0.0;
	double gradient_energy2 = 0.0;

	// name of temporary biinary file
	char *binfile_name = "cgtemp.tmp";


	void exit_handler();

	printf("Current working directory is:\n");
	system("cd");
	printf("\n");
	
	printf("Multi-Layer Perceptron (CG with OLF)\n");

	if(argc < 7)
	{
		char *exe_name = strrchr(argv[0], '\\');
		printf("Insufficient number of arguments. \n\nUsage: \n%s <training_file> <validation_file> <number_of_inputs> <number_of_outputs> <number_of_hidden_units> <number_of_iterations>",
			(exe_name)?exe_name+1:argv[0]);
		printf("\n\nPlease enter this information now:\n");

		printf("Enter training file name: ");
		safe_gets(fname_trg, 200);
	
		printf("Enter testing file name: ");
		safe_gets(fname_tst, 200);

		printf("Enter the number of inputs: ");
		scanf("%d", &N);

		printf("Enter the number of outputs: ");
		scanf("%d", &M);

		printf("Enter the number of hidden units: ");
		scanf("%d", &Nh);

		printf("Enter the number of BP iterations: ");
		scanf("%d", &Nit);

		atexit(exit_handler);
	}
	else
	{
		strcpy(fname_trg, argv[1]);
		strcpy(fname_tst, argv[2]);
		N = atoi(argv[3]);
		M = atoi(argv[4]);
		Nh = atoi(argv[5]);
		Nit = atoi(argv[6]);
	}

	result=fopen(RESULTS_FILE, "w");
	if(result == NULL)
	{
		printf("Unable to open file: %s\n", RESULTS_FILE);
		exit(0);
	}
	mse_file = fopen(MSE_FILE,"w");
	if(mse_file == NULL)
	{
		printf("Unable to open file: %s\n", MSE_FILE);
		exit(0);
	}
    
	time ( &rawtime );
	timeinfo = localtime ( &rawtime );
	fprintf (result, "The current date/time is: %s\n", asctime (timeinfo) );

	fprintf(result,"Name of the input file: %s\n\n", fname_trg);
	fprintf(result,"Name of the validation file: %s\n\n", fname_tst);
	fprintf(result,"Number of inputs, N = %d\n\n", N);
	fprintf(result,"Number of outputs, M = %d\n\n", M);
	fprintf(result,"Number of hidden units, Nh = %d\n\n", Nh);
	fprintf(result,"Number of iterations (for backpropagation), Nit = %d\n\n", Nit);

	L = N+1+Nh;

	x = createArray(L);
	mx = createArray(N);
	xv = createArray(N);
	ty = createArray(M);
	y = createArray(M);
	mm = createArray(Nh);
	vv = createArray(Nh);
	m_hu = createArray(Nh);
	sd_hu = createArray(Nh);
	net = createArray(Nh);
	O = createArray(Nh);
	Error = createArray(M);
	E = createArray(M);
	f = createArray(Nh);
	Fgx = createArray(Nh);
	Dyp = createArray(M);
	Et = createArray(M);
	delta_po = createArray(M);
	delta_p = createArray(Nh);

	r = createMatrix(L, L);
	c = createMatrix(M, L);

	w = createMatrix(Nh, N+1);
	g_wts = createMatrix(Nh, N+1);
	p = createMatrix(Nh, N+1);
	d = createMatrix(Nh, N+1);
	saved_w = createMatrix(Nh, N+1);

	Wout = createMatrix(M, L);
	goi_wts = createMatrix(M, N+1);
	goh_wts = createMatrix(M, Nh);
	poi = createMatrix(M, N+1);
	poh = createMatrix(M, Nh);
	saved_wo = createMatrix(M, L);

    
	infile = fopen(fname_trg, "r"); // input pattern file
	if(infile==NULL)
	{
		printf("Unable to open pattern file: %s\n", fname_trg);
		exit(0);
	}

	binfile = fopen("tempfile.txt", "w+");        // binary file used subsequently
	if(binfile==NULL)
	{
		printf("Unable to create temporary file:");
		exit(0);
	}

	// 1. Calculation of input means and variance.

	for(n=0; n<N; n++)
	{
		mx[n] = 0.0;
		xv[n] = 0.0;
	}
    
	Nv = 0;
	while(!feof(infile))
	{
		for(n=0; n<N; n++){
			fscanf(infile, "%lf ", &x[n]);}
			fscanf(infile, "\n");
		for(i=0; i<M; i++){
			fscanf(infile, "%lf ", &ty[i]);}
			fscanf(infile, "\n");
			
		if(feof(infile))
			break;

		for(n=0; n<N; n++)
		{
			mx[n] += x[n];
			xv[n] += x[n]*x[n];
		}

		Nv++;
	}
	printf("Nv = %d\n", Nv);

	for(n=0; n<N; n++)
	{
		mx[n] /= Nv;    // normalize means
		xv[n] /= Nv;    // normalize means
		xv[n] -= mx[n]*mx[n];
	}

	// 1b. Saving zero mean inputs to the binary file

	rewind(infile);
	Nv = 0;
	while(!feof(infile))
	{
		for(n=0; n<N; n++){
			fscanf(infile, "%lf ", &x[n]);}
			fscanf(infile, "\n");
		for(i=0; i<M; i++){
			fscanf(infile, "%lf ", &ty[i]);}
			fscanf(infile, "\n");
			
		if(feof(infile))
			break;

		for(n=0; n<N; n++)
			x[n] -= mx[n];


		for(n=0; n<N; n++){
			fprintf(binfile, "%lf ", x[n]);}
			fprintf(binfile,"\n");	
		for(i=0; i<M; i++){
			fprintf(binfile, "%lf ", ty[i]);}
			fprintf(binfile,"\n");

		Nv++;
	}
	printf("Wrote Nv = %d zero mean patterns to binary file.\n", Nv);

	fclose(infile);

	// 2. Random initialization of input weights.

	for(k=0; k<Nh; k++)
	{
		for(n=0; n<N; n++)
		{
			w[k][n] = slete(1.0, 0.0);  
			w[k][n] = w[k][n]/sqrt(xv[n]);
		}

		w[k][N]=slete(1.0, 0.0);

	}

	for(i=0; i<Nh; i++)
	{
		mm[i] = 0.0;
		vv[i] = 0.0;
	}

	for(i=0;i<Nh;i++)
	{
		for(k=0;k<N+1;k++)
		{
			mm[i]+=w[i][k];
			vv[i]+=w[i][k]*w[i][k];
		}
	}
	for(i=0; i<Nh; i++)
	{
		mm[i] /= (N+1);	// rohit: (1) changed to N+1
		vv[i] /= (N+1);
		vv[i] -= mm[i]*mm[i];
	}

	for (i=0;i<Nh;i++)
	{
		for(k=0;k<N;k++)
			w[i][k]=(w[i][k]-mm[i])/vv[i];
	}

	// 3. Calculation of means of hidden unit net functions' means and std-dev

	for(k=0; k<Nh; k++)
	{
		m_hu[k] = 0.0;
		sd_hu[k] = 0.0;
	}

	fflush(binfile);

	rewind(binfile);
	Nv = 0;
	x[N] = 1.0;
	while(!feof(binfile))
	{
		for(n=0; n<N; n++){
			fscanf(binfile, "%lf ", &x[n]);}
			fscanf(binfile, "\n");
		for(i=0; i<M; i++){
			fscanf(binfile, "%lf ", &ty[i]);}
			fscanf(binfile, "\n");

		if(feof(binfile))
			break;

		for(k=0; k<Nh; k++)
		{
			net[k] = 0.0;
			for(n=0; n<N+1; n++)
				net[k] += w[k][n]*x[n];
			m_hu[k] += net[k];
			sd_hu[k] += net[k]*net[k];
		}

		Nv++;
	}

	for(k=0; k<Nh; k++)
	{
		m_hu[k] /= Nv;
		sd_hu[k] /= Nv;
		sd_hu[k] -= m_hu[k]*m_hu[k];
		sd_hu[k] = sqrt(sd_hu[k]);
	}

#ifdef DEBUG_MODE
	printf("Nv = %d\n", Nv);
	for(k=0; k<Nh; k++)
		printf("%f,%f\n", m_hu[k], sd_hu[k]);
#endif

	// 4. Perform net control

	m_d = 0.5;               // Desired Hidden-unit means
	sd_d = 1;                // Desired Hidden-unit std-devs

	for(k=0; k<Nh; k++)
	{
		for(n=0; n<N+1; n++)
		{
			w[k][n] *= sd_d/sd_hu[k];
		}
		w[k][N] = w[k][N]-(m_hu[k]*sd_d/sd_hu[k])+m_d;	// rohit: (3) Fixed net control
	}

#ifdef DEBUG_MODE
	// 5. Verify net control

	for(k=0; k<Nh; k++)
	{
		m_hu[k] = 0.0;
		sd_hu[k] = 0.0;
	}

	rewind(binfile);
	Nv = 0;
	x[N] = 1.0;
	while(!feof(binfile))
	{
		for(n=0; n<N; n++){
			fscanf(binfile, "%lf ", &x[n]);}
			fscanf(binfile, "\n");
		for(i=0; i<M; i++){
			fscanf(binfile, "%lf ", &ty[i]);}
			fscanf(binfile, "\n");

		if(feof(binfile))
			break;

		for(k=0; k<Nh; k++)
		{
			net[k] = 0.0;
			for(n=0; n<N+1; n++)
				net[k] += w[k][n]*x[n];
			m_hu[k] += net[k];
			sd_hu[k] += net[k]*net[k];
		}

		Nv++;
	}

	for(k=0; k<Nh; k++)
	{
		m_hu[k] /= Nv;
		sd_hu[k] /= Nv;
		sd_hu[k] -= m_hu[k]*m_hu[k];
		sd_hu[k] = sqrt(sd_hu[k]);
	}


	printf("Nv = %d\n", Nv);
	for(k=0; k<Nh; k++)
		printf("%f,%f\n", m_hu[k], sd_hu[k]);
#endif

	// 6. Calculation of R and C matrices
	for(i=0; i<M; i++)
		E[i] = 0.0;

	for(n=0; n<L; n++)
	{
		for(i=0; i<L; i++)
			r[n][i] = 0;
	}

	for(i=0; i<M; i++)
	{
		for(n=0; n<L; n++)
			c[i][n] = 0;
		Et[i] = 0;
	}



	rewind(binfile);
	Nv = 0;
	x[N] = 1.0;
	while(!feof(binfile))
	{
		for(n=0; n<N; n++){
			fscanf(binfile, "%lf ", &x[n]);}
			fscanf(binfile, "\n");
		for(i=0; i<M; i++){
			fscanf(binfile, "%lf ", &ty[i]);}
			fscanf(binfile, "\n");

		if(feof(binfile))
			break;

		for(i=0; i<M; i++)
			Et[i] += ty[i]*ty[i];

		for(k=0; k<Nh; k++)
		{
			net[k] = 0.0;
			for(n=0; n<N+1; n++)
				net[k] += w[k][n]*x[n];
			O[k] = 1/(1+exp(-net[k]));  // Activation function
			x[N+1+k] = O[k];        // New basis vector
		}
		for(n=0; n<L; n++)
		{
			for(i=0; i<=n; i++)	// rohit: (2) added symmetry
				r[n][i] += x[n]*x[i];
			for(i=0; i<M; i++)
				c[i][n] += ty[i]*x[n];
		}

		Nv++;
	}
	//printf("b. Nv = %d\n", Nv);

	for(n=0; n<L; n++)
	{
		for(i=0; i<=n; i++)
		{
			r[n][i] /= Nv;
			r[i][n] = r[n][i];
		}
	}

	for(i=0; i<M; i++)
	{
		for(n=0; n<L; n++)
			c[i][n] /= Nv;
		Et[i] /= Nv;
	}

	// 7. Perform OWO to calculate output weights
	for (i=0;i<M;i++)
		for(j=0;j<L;j++)
			Wout[i][j]=0.0;
	
	MSE =  OLS1(Wout, N,M,L,r,c,Et);

	fprintf(result, "Initial network error = %lf\n", MSE);
	printf("Initial network error = %lf\n",MSE);

	// Initialize old weights 
	// rohit: (3) moved saving weight to after OLS call
	for(k=0; k<Nh; k++)
	{
		for(n=0; n<N+1; n++)
			saved_w[k][n] = w[k][n];
	}

	for(i=0; i<M; i++)
	{
		for(n=0; n<L; n++)
			saved_wo[i][n] = Wout[i][n];	// rohit: (4) initialize to OLS weight instead of 0
	}



#ifdef DEBUG_MODE
	printf("\nR = [ < G , P > + < Goi , Poi > + < Goh , Poh >]/Xn\n");
	printf("Xn = E(g_oi) + E(g_oh) + E(g_hi)\n");

	fprintf(result,"\n\nIt#\tZ\t\tMSE\t\tR\t\tXn\t\tIerr\n\n");
	printf("\n\nIt#\tZ\t\tMSE\t\tR\t\tXn\t\tIerr\n\n");
#else
	fprintf(result,"\n\nIt#\tZ\t\tMSE\n\n");
	printf("\n\nIt#\tZ\t\tMSE\n\n");
#endif

	for(k=0; k<Nh; k++)
	{
		for(n=0; n<N+1; n++)
		{
			p[k][n]=0.0;
		}
	}

	for(i=0; i<M; i++)
	{
		for(n=0; n<N+1; n++)        
		{
			poi[i][n]=0;
		}

		for(k=0; k<Nh; k++)        
		{
			poh[i][k] =0;
		}
	}

	Xd=1;
	B1=0;
	z = 0.05;	// rohit: apparantly it is independent of initialization
	zp = z;		// rohit: zp was not intialized, but initialization is not required


	Ep = 1e20;	// rohit: (5) changed type of Ep to double from int
	reset_Ierr = 1;

	for(it = 0; it<=Nit; it++)
	{
		//if the error increased in the last iteration, Ierr is allowed to increment
		if(reset_Ierr == 0)
			reset_Ierr = 1;
		else
			Ierr = 0;

		// Initializing gradients to zero for this iteration

		for(k=0; k<Nh; k++)
		{
			for(n=0; n<N+1; n++)
			{
				g_wts[k][n] = 0.0;
			}
		}

		for(i=0; i<M; i++)
		{
			for(n=0; n<N+1; n++)        
			{
				goi_wts[i][n] =0;
			}

			for(k=0; k<Nh; k++)        
			{
				goh_wts[i][k] =0;
			}
		}

		// A pass through the data file
		rewind(binfile);
		x[N] = 1.0;
		for (i=0;i<M;i++)
		{
			Error[i]=0;
			y[i]=0;
		}


		Nv = 0;
		while(!feof(binfile))
		{
			for(n=0; n<N; n++){
			fscanf(binfile, "%lf ", &x[n]);}
			fscanf(binfile, "\n");
		for(i=0; i<M; i++){
			fscanf(binfile, "%lf ", &ty[i]);}
			fscanf(binfile, "\n");

			if(feof(binfile))
				break;

			for(k=0; k<Nh; k++)
			{
				net[k] = 0.0;
				for(n=0; n<N+1; n++)
					net[k] += w[k][n]*x[n];
				O[k] = 1/(1+exp(-net[k]));  // Activation function
				f[k] =  O[k]*(1-O[k]);
				x[N+1+k] = O[k];        // New basis vector
			}



			// calculate output deltas;

			for(i=0; i<M; i++)
			{
				temp = 0.0;         
				for(n=0; n<L; n++)
					temp += Wout[i][n] * x[n];
				y[i]= temp;
				delta_po[i] = 2*(ty[i] - temp);
			}

			for(i=0;i<M;i++)
				Error[i]+=(y[i]-ty[i])*(y[i]-ty[i]);

			// calculate input deltas

			for(k=0; k<Nh; k++)
			{
				temp = 0.0;
				for(i=0; i<M; i++)
					temp += delta_po[i]*Wout[i][N+k+1];
				delta_p[k] = temp * O[k]*(1-O[k]); 
			}

			// calculate gradients

			for(k=0; k<Nh; k++)
			{
				for(n=0; n<N+1; n++)        
				{
					g_wts[k][n] += delta_p[k]*x[n];
				}

			}
			for(i=0; i<M; i++)
			{
				for(n=0; n<N+1; n++)        
					goi_wts[i][n] += delta_po[i]*x[n];

				for(k=0; k<Nh; k++)        
					goh_wts[i][k] += delta_po[i]*x[N+1+k];	// rohit: (6) changed N+k to N+1+k
			}


			Nv++;	
		} // while(!feof(binfile))



		Err=0;
		for ( i=0;i<M;i++)
		{
			Error[i]= Error[i]/Nv;
			Err+=Error[i];
			//	printf("%d\t%lf\t%lf   \n", it+1,z,Error[i]);
		}

		fprintf(mse_file, "%lf\n", Err);

		if(Err <Ep)
		{

			for (i=0;i<M;i++)
				for(j=0;j<L;j++)
					saved_wo[i][j]= Wout[i][j];

			for(k=0; k<Nh; k++)
			{
				for(n=0; n<N+1; n++)
					saved_w[k][n] = w[k][n];
			}

		}
		else
		{
			Ierr = Ierr +1;
			printf("Reading the previous weights back again\n");

			for (i=0;i<M;i++)
				for(j=0;j<L;j++)
					Wout[i][j]=saved_wo[i][j];

			for(k=0; k<Nh; k++)
			{
				for(n=0; n<N+1; n++)
					w[k][n] = saved_w[k][n];
			}
			Err = Ep;
			z = zp/2;
			// Change hidden weights

			for(k=0; k<Nh; k++)
			{
				for(n=0; n<N+1; n++)
					w[k][n] += z * p[k][n];

			}

			for(i=0; i<M; i++)
			{
				for(n=0; n<N+1; n++)        
					Wout[i][n] += z*poi[i][n];

				for(k=0; k<Nh; k++)        
					Wout[i][N+k+1] += z*poh[i][k];               
			}

			reset_Ierr = 0;	// rohit: (7) Ierr was reset to zero at the beginning of the iteration, 
			// it should not be done when the error has increased and we are counting
			//printf("It = %d, Ierr = %d\n", it, Ierr);

			continue;

		}




		// Normalize the gradients
		for(k=0; k<Nh; k++)
		{
			for(n=0; n<N+1; n++)
			{
				g_wts[k][n] /= Nv;
			}
		}
		for(i=0; i<M; i++)
		{
			for(n=0; n<N+1; n++)        
				goi_wts[i][n] /=Nv;

			for(k=0; k<Nh; k++)        
				goh_wts[i][k] /=Nv;


		}

		gradient_energy = 0.0;
		gradient_energy1 = 0.0;
		gradient_energy2 = 0.0;
		for(k=0; k<Nh; k++)
		{
			for(n=0; n<N+1; n++)       
			{
				gradient_energy += g_wts[k][n]*g_wts[k][n];
			}

		}

		for(i=0; i<M; i++)
		{
			for(n=0; n<N+1; n++)        
			{
				gradient_energy1 += goi_wts[i][n]*goi_wts[i][n];
			}

		}

		for(k=0; k<Nh; k++)
		{
			for(i=0; i<M; i++)        
			{
				gradient_energy2 += goh_wts[i][k]*goh_wts[i][k];
			}

		}

		Xn = gradient_energy + gradient_energy1 + gradient_energy2;
		B1=Xn/Xd;
		Xd=Xn;

		// rohit: find R for printing R = [ < G , P > + < Goi , Poi > + < Goh , Poh >]/Xn

		R = 0.0;
#ifdef DEBUG_MODE
		R = innerProd(g_wts, p, Nh, N+1) + innerProd(goi_wts, poi, M, N+1) + innerProd(goh_wts, poh, M, Nh);
#endif

		if (Ierr >0)
		{
			for(k=0; k<Nh; k++)
			{
				for(n=0; n<N+1; n++)
				{
					p[k][n]=0.0;
				}
			}

			for(i=0; i<M; i++)
			{
				for(n=0; n<N+1; n++)        
				{
					poi[i][n]=0;
				}

				for(k=0; k<Nh; k++)        
				{
					poh[i][k] =0;
				}
			}
		}

		for(k=0; k<Nh; k++)
		{
			for(n=0; n<N+1; n++)
				p[k][n]= B1*p[k][n]+g_wts[k][n];
		}

		for(i=0; i<M; i++)
		{
			for(n=0; n<N+1; n++)        
				poi[i][n]= B1*poi[i][n]+ goi_wts[i][n];

			for(k=0; k<Nh; k++)        
				poh[i][k] =B1*poh[i][k]+goh_wts[i][k];
		}

		if(it!=0)
		{
#ifdef DEBUG_MODE
			fprintf(result, "%d\t%lf\t%lf\t%lf\t%lf\t%d \n", it,z,Err, R, Xn, Ierr);
			printf("%d\t%lf\t%lf\t%lf\t%lf\t%d \n", it,z,Err, R, Xn, Ierr);
#else
			fprintf(result, "%d\t%lf\t%lf \n", it,z,Err);
			printf("%d\t%lf\t%lf \n", it,z,Err);
#endif
		}

		if (it!= Nit)
		{

			// Calculation of invisible learning factor
			x[N]=1.0;

			rewind(binfile);


			dez=0;
			dez2 =0;

			while(!feof(binfile))
			{
            
            for(n=0; n<N; n++){
			fscanf(binfile, "%lf ", &x[n]);}
			fscanf(binfile, "\n");
		    for(i=0; i<M; i++){
			fscanf(binfile, "%lf ", &ty[i]);}
			fscanf(binfile, "\n");

				if(feof(binfile))
					break;

				for(k=0; k<Nh; k++)
				{
					net[k] = 0.0;
					for(n=0; n<N+1; n++)
						net[k] += w[k][n]*x[n];
					O[k] = 1/(1+exp(-net[k]));  // Activation function
					x[N+1+k] = O[k];        // New basis vector
				}
				for(i=0; i<M; i++)
				{
					y[i] = 0.0;         
					for(n=0; n<L; n++)
						y[i] += Wout[i][n] * x[n];
				}

				for (k=0;k<Nh;k++)
				{
					Fgx[k]=0;
					for (n=0;n<N+1;n++)
						Fgx[k]+=p[k][n]*x[n];
					Fgx[k]=Fgx[k]*O[k]*(1-O[k]);
				}

				for(i=0;i<M;i++)
				{
					Dyp[i]=0;
					for(k=0;k<Nh;k++)
						Dyp[i]+=poh[i][k]*O[k]+Wout[i][N+k+1]*Fgx[k];

					for (n=0;n<N+1;n++)
						Dyp[i]+=poi[i][n]*x[n];

					dez+= (ty[i]-y[i])*Dyp[i];
					dez2+=Dyp[i]*Dyp[i];
				}
			}	// end of data pass

			z = dez/dez2;

			// Change hidden weights

			for(k=0; k<Nh; k++)
			{
				for(n=0; n<N+1; n++)
					w[k][n] += z * p[k][n];
			}

			for(i=0; i<M; i++)
			{
				for(n=0; n<N+1; n++)        
					Wout[i][n] += z*poi[i][n];

				for(k=0; k<Nh; k++)        
					Wout[i][N+k+1] += z*poh[i][k];   // was N+k in CGpraveen code
			}

			zp=z;
		} // end of if(last iteration)

	}
	printf("End of CG iterations.\n");

	fclose(binfile);
	remove("tempfile.txt");

	fclose(mse_file);

	// SAVING THE WEIGHTS TO A NETWORK FILE
	weights_file=fopen(WEIGHTS_FILE,"w");
	if(weights_file == NULL)
	{
		printf("Unable to write weights to file %s\n", WEIGHTS_FILE);
		exit(0);
	}
	fprintf(weights_file, "500\n");
	fprintf(weights_file, fname_trg);
	fprintf(weights_file, "\n");
	fprintf(weights_file, "%d\n", N);
	fprintf(weights_file, "%d\n", M);
	fprintf(weights_file, "%d\n", Nh);
	fprintf(weights_file, "The_Training_Algorithm_is_CG");

	// adjust weights for input means
	for(i=0; i<M; i++)
	{
		temp = 0.0;
		for(j=0; j<N; j++)
			temp += Wout[i][j] * mx[j];
		Wout[i][N] -= temp;
	}

	for(i=0; i<Nh; i++)
	{
		temp = 0.0;
		for(j=0; j<N; j++)
			temp += w[i][j] * mx[j];
		w[i][N] -= temp;
	}

	// saving weights from input & hidden units to output units
	for(i=0;i<L;i++)
	{
		fprintf(weights_file, "\n");
		for(k=0;k<M;k++)
			fprintf(weights_file, "%lf\t", Wout[k][i]);
	}

	// saving weights from input units to hidden units 
	for(j=0;j<Nh;j++)
	{
		fprintf(weights_file, "\n");
		for(i=0;i<N+1;i++)
			fprintf(weights_file, "%lf\t", w[j][i]);
	}   
	fclose(weights_file);
	// END OF SAVING WEIGHTS

	printf("\n!  Training Complete for MLP with %d Hidden Layers over %d iterations  !\n", Nh, Nit );


	//8. Recalculate error from real input (validation)
	printf("Performing validation on %s\n", fname_tst);

	for(i=0; i<M; i++)
		E[i] = 0.0;

	infile = fopen(fname_tst, "r");
	if(infile==NULL)
	{
		printf("Unable to open test file %s.\n", fname_tst);
		exit(0);
	}

	Nv = 0;
	x[N] = 1.0;
	while(!feof(infile))
	{
		for(n=0; n<N; n++)
		{
			fscanf(infile, "%lf", &x[n]);
			// x[n] -= mx[n]; // input means have been adjusted for in the weights
		}
		for(i=0; i<M; i++)
			fscanf(infile, "%lf", &ty[i]);
		if(feof(infile))
			break;

		for(k=0; k<Nh; k++)
		{
			net[k] = 0.0;
			for(n=0; n<N+1; n++)
				net[k] += w[k][n]*x[n];
			O[k] = 1/(1+exp(-net[k]));  // Activation function
			x[N+1+k] = O[k];        // New basis vector
		}
		for(i=0; i<M; i++)
		{
			temp = ty[i];
			for(n=0; n<L; n++)
				temp -= Wout[i][n]*x[n];
			E[i] += temp*temp;
		}
		Nv++;
	}
	printf("Nv (TEST file) = %d\n", Nv);
	fprintf(result, "Nv (TEST file) = %d\n", Nv);

	MSE = 0.0;
	for(i=0; i<M; i++)
	{
		MSE += E[i];
	}
	MSE /= Nv;

	printf("Testing error: %lf\n", MSE);
	fprintf(result, "Testing error: %lf\n", MSE);
	fclose(result);

	freeArray(x);
	freeArray(mx);
	freeArray(xv);
	freeArray(ty);
	freeArray(y);
	freeArray(mm);
	freeArray(vv);
	freeArray(m_hu);
	freeArray(sd_hu);
	freeArray(net);
	freeArray(O);
	freeArray(Error);
	freeArray(E);
	freeArray(f);
	freeArray(Fgx);
	freeArray(Dyp);
	freeArray(Et);
	freeArray(delta_po);
	freeArray(delta_p);

	freeMatrix(r, L);
	freeMatrix(c, M);

	freeMatrix(w, Nh);
	freeMatrix(g_wts, Nh);
	freeMatrix(p, Nh);
	freeMatrix(d, Nh);
	freeMatrix(saved_w, Nh);

	freeMatrix(Wout, M);
	freeMatrix(goi_wts, M);
	freeMatrix(goh_wts, M);
	freeMatrix(poi, M);
	freeMatrix(poh, M);
	freeMatrix(saved_wo, M);

	return 0;
}

void exit_handler()
{
	putchar('\n');
	system("PAUSE");
}

/*******************************************/


double OLS1(double **Wout, int N, int M, int L, double **R, double **C,double *Et)
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
		{
			A[i][j]=0.0;
		}

		for (j = 0; j < M; j++)
		{
			Wt[j][i] = 0.0;
		}
	}

	for (i = 0; i < M; i++ )
	{
		E[i] = 0.0;		
	}

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



/******************************************************************************
Subroutine : Random no. Generator
-the random nos. generated are uniformly distributed
between 0 and 1.
-IX, IY, IZ are the SEEDS
*******************************************************************************/

double rand1(int *ix, int *iy, int *iz)
{
	int ixx, iyy, izz;
	double itemp;
	double temp;


	ixx = (*ix) / 177;
	*ix = 171 * (*ix % 177) - 2 * ixx;

	if(*ix < 0)
		*ix += 30269;

	iyy = (*iy) / 176;
	*iy = 176 * (*iy % 176) - 2 * iyy;

	if(*iy < 0)
		*iy += 30307;

	izz = (*iz) / 178;
	*iz = 170 * (*iz % 178) - 2 * izz;

	if(*iz < 0)
		*iz += 30323;

	temp = (double)(*ix) / 30629.0 + (double)(*iy) / 30307.0 + (double)(*iz) / 30323.0;
	itemp = floor(temp);
	return (temp - itemp);
}

/********* Gaussian Random No. generator **************************************/

double slete(double std, double xmean)
{
	double rand1(int *, int *, int *);
	// Seeds for the random number generator function
#define PI 3.1415926    
	static int IX = 3;
	static int IY = 4009;
	static int IZ = 234;	
	return (xmean + std * cos(2 * PI * rand1(&IX, &IY, &IZ)) * sqrt(-2.0 * log(rand1(&IX, &IY, &IZ))));
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
