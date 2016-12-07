#include <string>
#include <iomanip>
#include "utils.h"
#include <valarray>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include "mlp.h"
#include <cassert>
#ifndef MYDEFS_H
#define MYDEFS_H
#include "mydefs.h"
#endif
#include "matops.h"
#include <string.h>
extern double slete(double, double);
extern double ols(int,matrix &,matrix &,matrix &,int,const valarray<double> &,double *,matrix &);
using namespace std;
int ols_linear_solve(valarray<double> &,matrix &,valarray<double>,vector<int> & );
matrix compute_HWO_direction(const matrix&,const matrix&, size_t&);
valarray<double> compute_newton_direction(const matrix &,const valarray<double> &,size_t& );
std::vector<double> MLP::train_CG(int Nit,size_t N,size_t M,size_t Nh,string fname)
{
	
	double OLS1(double **, int, int, int, double **, double **, double *);
	extern double* createArray(int);
        extern double** createMatrix(int, int);
        extern void freeArray(double *);
        extern void freeMatrix(double **, int);
         extern void safe_gets(char *, int);
        extern double innerProd(double **, double **, int, int);
	//time_t rawtime;
	//struct tm * timeinfo;

	char fname_trg[200];
	char fname_tst[200];
	int L, Nv;
	FILE *infile, *binfile;
	FILE *result, *mse_file, *weights_file;

        vector<double> results;
	int n, k, i, it,j;
	int Ierr, reset_Ierr;

	//double *x,*ty;
       double  **goi_wts,**goh_wts,**p,**poh,**poi;
	//double *mx,*xv,*net,*mm,*vv;
	double **w,*Fgx,*Dyp;
	//double *m_hu, *sd_hu,
        double m_d, sd_d,Err,Xd,B1,Xn,zp;

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

	printf("\n");
	
	printf("Multi-Layer Perceptron (CG with OLF)\n");


		strcpy(fname_trg, fname.c_str());
		strcpy(fname_tst, fname.c_str());


	L = N+1+Nh;

        valarray<double> x(L),mx(N),xv(N),net(Nh),mm(Nh),vv(Nh),m_hu(Nh),sd_hu(Nh),y(M),O(Nh),f(Nh);

        valarray<double> ty(M),Error(M);

	//Error = createArray(M);
	E = createArray(M);
	//f = createArray(Nh);
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


	Nv = 0;
	while(!feof(infile))
	{
		for(n=0; n<N; n++)
			fscanf(infile, "%lf ", &x[n]);

		if(feof(infile))
                    break;

		for(i=0; i<M; i++)
			fscanf(infile, "%lf ", &ty[i]);

                mx += x;
                xv += x*x;
		Nv++;
	}
	printf("Nv = %d\n", Nv);
        mx /= Nv;
        xv /= Nv;
        xv -= mx*mx;


	// 1b. Saving zero mean inputs to the binary file

	rewind(infile);
	Nv = 0;
	while(!feof(infile))
	{
		for(n=0; n<N; n++)
		  fscanf(infile, "%lf ", &x[n]);
                if(feof(infile))
                    break;

		for(i=0; i<M; i++)
			fscanf(infile, "%lf ", &ty[i]);

			x -= mx;

		for(n=0; n<N; n++){
			fprintf(binfile, "%lf ", x[n]);}
			fprintf(binfile,"\n");	
		for(i=0; i<M; i++){
			fprintf(binfile, "%lf ", ty[i]);}
			fprintf(binfile,"\n");

		Nv++;
	}
	printf("Wrote Nv = %d zero mean patterns to binary file.\n", Nv);

	//fclose(infile);

	// 2. Random initialization of input weights.

	for(k=0; k<Nh; k++)
	{ 
              w[k][N]=slete(0,1);
		for(n=0; n<N; n++)
		{
                  w[k][n] = slete(0.0, 1.0);  
		w[k][n] = w[k][n]/sqrt(xv[n]);
		}
	}

	for(i=0;i<Nh;i++)
	{
		for(k=0;k<N+1;k++)
		{
			mm[i]+=w[i][k];
			vv[i]+=w[i][k]*w[i][k];
		}
	}


//        mm /= (N+1);
  //      vv /= (N+1);
    //    vv =- mm*mm;
//	for (i=0;i<Nh;i++)
//	{
//		for(k=0;k<N;k++);
//			w[i][k]=(w[i][k]-mm[i])/vv[i];
//	}

	// 3. Calculation of means of hidden unit net functions' means and std-dev

	//fflush(binfile);

	rewind(infile);
	Nv = 0;
	
	while(!feof(infile))
	{
		for(n=0; n<N; n++){
			fscanf(infile, "%lf ", &x[n]);}
			//fscanf(binfile, "\n");
                   if(feof(infile))
			break;
		for(i=0; i<M; i++)
			fscanf(infile, "%lf ", &ty[i]);
			//fscanf(binfile, "\n");

		x-=mx;
               x[N] = 1.0;
		for(k=0; k<Nh; k++)
		{
			net[k] = 0.0;
			for(n=0; n<N+1; n++)
				net[k] += w[k][n]*x[n];
		}
                m_hu += net;
                sd_hu += net*net;
		Nv++;
	}
        m_hu /= Nv;
        sd_hu /= Nv;
        sd_hu -= m_hu*m_hu;
        sd_hu=sqrt(sd_hu);

	// 4. Perform net control

	m_d = 0.5;               // Desired Hidden-unit means
	sd_d = 1;                // Desired Hidden-unit std-devs

	for(k=0; k<Nh; k++)
	{
		for(n=0; n<N+1; n++)
			w[k][n] *= sd_d/sd_hu[k];
	
		w[k][N] = w[k][N]-(m_hu[k]*sd_d/sd_hu[k])+m_d;	// rohit: (3) Fixed net control
	}


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

	rewind(infile);
	Nv = 0;
	x[N] = 1.0;
	while(!feof(infile))
	{
		for(n=0; n<N; n++)
			fscanf(infile, "%lf ", &x[n]);
//			fscanf(binfile, "\n");
                if(feof(infile))
                        break;

		for(i=0; i<M; i++)
			fscanf(infile, "%lf ", &ty[i]);
		//	fscanf(binfile, "\n");

//		if(feof(binfile))
//			break;
                x-=mx;
 x[N] = 1.0;
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
			//for(i=0; i<=n; i++)	// rohit: (2) added symmetry
                          for (i=0;i<L;i++)
				r[n][i] += x[n]*x[i];
			for(i=0; i<M; i++)
				c[i][n] += ty[i]*x[n];
		}

		Nv++;
	}
	printf("b. Nv = %d\n", Nv);
       // ofstream rfile("R.txt");
	for(n=0; n<L; n++)
	{
		//for(i=0; i<=n; i++)
                for (i=0;i<L;i++)
		{
			r[n][i] /= Nv;
         //         rfile << r[n][i] << " ";
			//r[i][n] = r[n][i];
		}
           // rfile << endl;
	}

	for(i=0; i<M; i++)
	{
		for(n=0; n<L; n++)
			c[i][n] /= Nv;
		Et[i] /= Nv;
	}
       // cout << r[0][0] << " " << r[0][1] << " " << r[0][2] << endl;
	// 7. Perform OWO to calculate output weights
	for (i=0;i<M;i++)
		for(j=0;j<L;j++)
			Wout[i][j]=0.0;
	
	MSE =OLS1(Wout,N,M,L,r,c,Et);

	//fprintf(result, "Initial network error = %lf\n", MSE);
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

	//fprintf(result,"\n\nIt#\tZ\t\tMSE\t\tR\t\tXn\t\tIerr\n\n");
	printf("\n\nIt#\tZ\t\tMSE\t\tR\t\tXn\t\tIerr\n\n");
#else
	//fprintf(result,"\n\nIt#\tZ\t\tMSE\n\n");
	printf("\n\nIt#\tZ\t\tMSE\n\n");
#endif

	for(k=0; k<Nh; k++)
		for(n=0; n<N+1; n++)
			p[k][n]=0.0;
	

	for(i=0; i<M; i++)
	{
		for(n=0; n<N+1; n++)        
			poi[i][n]=0;	
		for(k=0; k<Nh; k++)        
			poh[i][k] =0;
	
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
			for(n=0; n<N+1; n++)
				g_wts[k][n] = 0.0;
		
		for(i=0; i<M; i++)
		{
			for(n=0; n<N+1; n++)        
				goi_wts[i][n] =0;

			for(k=0; k<Nh; k++)        
				goh_wts[i][k] =0;
		
		}

		// A pass through the data file
		rewind(infile);
		
		Nv = 0;
		while(!feof(infile))
		{
			for(n=0; n<N; n++)
			fscanf(infile, "%lf ", &x[n]);
		
                    if(feof(infile))
				break;
		for(i=0; i<M; i++)
			fscanf(infile, "%lf ", &ty[i]);
			x-=mx;
               x[N] = 1.0;

			for(k=0; k<Nh; k++)
			{
				net[k] = 0.0;
				for(n=0; n<N+1; n++)
					net[k] += w[k][n]*x[n];
				O[k] = 1/(1+exp(-net[k]));  // Activation function
				
				x[N+1+k] = O[k];        // New basis vector
			}

f =  O*(1.-O);

			// calculate output deltas;

			for(i=0; i<M; i++)
			{
				temp = 0.0;         
				for(n=0; n<L; n++)
					temp += Wout[i][n] * x[n];
				y[i]= temp;
				delta_po[i] = 2*(ty[i] - temp);
			}

                         Error += (y-ty)*(y-ty);
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
				for(n=0; n<N+1; n++)        
					g_wts[k][n] += delta_p[k]*x[n];
			
			for(i=0; i<M; i++)
			{
				for(n=0; n<N+1; n++)        
					goi_wts[i][n] += delta_po[i]*x[n];

				for(k=0; k<Nh; k++)        
					goh_wts[i][k] += delta_po[i]*x[N+1+k];	// rohit: (6) changed N+k to N+1+k
			}

			Nv++;	
		} // while(!feof(binfile))

			Error/= Nv;
			Err=Error.sum();

		if(Err <Ep)
		{
			for (i=0;i<M;i++)
				for(j=0;j<L;j++)
					saved_wo[i][j]= Wout[i][j];

			for(k=0; k<Nh; k++)
				for(n=0; n<N+1; n++)
					saved_w[k][n] = w[k][n];
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
				for(n=0; n<N+1; n++)
					w[k][n] += z * p[k][n];

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
			for(n=0; n<N+1; n++)
				g_wts[k][n] /= Nv;
		
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
			for(n=0; n<N+1; n++)       
				gradient_energy += g_wts[k][n]*g_wts[k][n];
		
		for(i=0; i<M; i++)
			for(n=0; n<N+1; n++)        	
				gradient_energy1 += goi_wts[i][n]*goi_wts[i][n];

		for(k=0; k<Nh; k++)
			for(i=0; i<M; i++)        
				gradient_energy2 += goh_wts[i][k]*goh_wts[i][k];
	
		Xn = gradient_energy + gradient_energy1 + gradient_energy2;
		B1=Xn/Xd;
		Xd=Xn;

		if (Ierr >0)
		{
			for(k=0; k<Nh; k++)
				for(n=0; n<N+1; n++)
					p[k][n]=0.0;
			
			for(i=0; i<M; i++)
			{
				for(n=0; n<N+1; n++)        
					poi[i][n]=0;
				
				for(k=0; k<Nh; k++)        
					poh[i][k] =0;
				
			}
		}

		for(k=0; k<Nh; k++)
			for(n=0; n<N+1; n++)
				p[k][n]= B1*p[k][n]+g_wts[k][n];
	
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
	//		fprintf(result, "%d\t%lf\t%lf\t%lf\t%lf\t%d \n", it,z,Err, R, Xn, Ierr);
			printf("%d\t%lf\t%lf\t%lf\t%lf\t%d \n", it,z,Err, R, Xn, Ierr);
#else
	//		fprintf(result, "%d\t%lf\t%lf \n", it,z,Err);
			printf("%d\t%lf\t%lf \n", it,z,Err);
#endif
		}

		if (it!= Nit)
		{

			// Calculation of invisible learning factor
		
			rewind(infile);
			dez=0;
			dez2 =0;

			while(!feof(infile))
			{
            
            		for(n=0; n<N; n++)
			  fscanf(infile, "%lf ", &x[n]);
		
     			if(feof(infile))
				break;
		    for(i=0; i<M; i++)
			fscanf(infile, "%lf ", &ty[i]);
		
			x-=mx;
				x[N]=1.0;

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

	freeArray(E);
//	freeArray(f);
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

	return results;

}

vector<double> MLP::train_OWO_HWO(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int k=0;
   size_t nlin=0;
   double Eg;
   double alpha,MSE,d1,d2;
   int Nu;
   const double dmean = 0.5;
   const double dstd = 1.;
   vector<double> return_vals;
   int M=outputs,N=inputs,Nh=hidden;
   Nu=N+Nh+1;
   matrix Wa(Nu,M),G_hwo(Nh,N+1);
 
   int Niw=(N+1)*Nh;
   valarray<double> g(Niw),d(Niw);
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   while(k < Nit)
    {
      MSE=owo();
      return_vals.push_back(MSE);
      compute_derivs(false,Eg);
      G_hwo=compute_HWO_direction(Ri,G,nlin);
      g=unroll(G_hwo);
      alpha=linesearch(d1,d2,G_hwo);
      update_weights(alpha,G_hwo);
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << alpha ;
      cout << setw(13) << d1 ;
      cout << setw(13) << d2;
      cout << setw(13) << nlin << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << alpha ;
      report << setw(13) << d1;
      report << setw(13) << d2;
      report << setw(13) << nlin << endl;

      k++;
     } //while
                
      report.close();
      return(return_vals);
}

vector<double> MLP::train_OWO_BP(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int k=0;
   size_t nlin=0;
   double Eg;
   double alpha,MSE,d1,d2;
   const double dmean = 0.5;
   const double dstd = 1.;
   vector<double> result_vals;
   init_weights();

   calculate_net_stats();

   net_control(dmean,dstd);
  ofstream linefile("owo-bp.txt");
         linefile << W << endl << endl;
   for (size_t i=0;i< net_mean.size();i++)
      linefile << net_mean[i] << " " << net_std[i] << endl;
      linefile.close();
   while(k < Nit)
   {
      MSE=owo();
      result_vals.push_back(MSE);
      compute_derivs(false,Eg);
      alpha=linesearch(d1,d2,G);
      update_weights(alpha,G);
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << alpha ;
      cout << setw(13) << d1;
      cout << setw(13) << d2;
      cout << setw(13) << nlin << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << alpha ;
      report << setw(13) << d1;
      report << setw(13) << d2;
      report << setw(13) << nlin << endl;

      k++;
     } //while
                
      report.close();
      return(result_vals);
}
vector<double> MLP::train_LM(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
  vector<double> return_vals;
  ofstream report;
  int k=0,i,j,iters;
  size_t nlin;
   double Eg;
   double MSE,lambda,alpha;
   int Nu;
   double MSE_after_weight_change; 
   const double dmean = 0.5;
   const double dstd = 1.;
   int M=outputs,N=inputs,Nh=hidden;
   Nu=N+Nh+1;
   matrix Wa(Nu,M);
 
   int Niw=(N+1)*Nh, Now=M*(N+Nh+1),Nw=Niw+Now;
   valarray<double> g(Niw),d(Niw),dt(Nw),gt(Nw),tempW(Niw),gsave(Nw);
   valarray<double> tempWoi(M*(N+1)),tempWoh(M*Nh);
   matrix DW(Nh,N+1),DWoh(M,Nh),DWoi(M,N+1),L(Now,Now),Hsave(Nw,Nw);
   valarray<double> diags(Ht.getRow());
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
 //  for (i=0;i<M;i++)
   // {
     // for (j=0;j<N+1;j++)
       //  Woi[i][j]=slete(0,1);
      // for (int j1=0;j1<Nh;j1++)
       //  Woh[i][j1]=slete(0,1);
   // }
   MSE=owo();
   //MSE=1e20;
   lambda=1e-2;
   while(k < Nit)
   {
       return_vals.push_back(MSE);

       compute_derivs(true,Eg);
       for (i=0;i<Ht.getRow();i++)
          Ht[i][i] += lambda;
       Hsave=Ht; //save the Hessian...
       gsave=Gt; //and the gradient
       dt=compute_newton_direction(Ht,Gt,nlin);
       L=chol(Ht);
       tempW=dt[slice(0,Niw,1)];
       tempWoi=dt[slice(Niw,M*(N+1),1)];
       tempWoh=dt[slice(Niw+M*(N+1),M*Nh,1)];
       DW=reshape(tempW,Nh,N+1);
       DWoi=reshape(tempWoi,M,N+1);
       DWoh=reshape(tempWoh,M,Nh);
   //    alpha=linesearch(0,1,DW,1e-3,iters);
       W += DW;
       Woi += DWoi;
       Woh += DWoh;
       MSE_after_weight_change=J2(W);
       if (MSE > MSE_after_weight_change) 
        {
          lambda/=10;
          cout << MSE << endl;
          MSE=MSE_after_weight_change;
     //     Hsave=Ht;
     //     gsave=Gt;
        }
       else
       {
         lambda *=10;
         Ht=Hsave;
         Gt=gsave;
         W -= DW;
         Woi -= DWoi;
         Woh -= DWoh;
         k--;
       }
      cout << setw(4)  << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << lambda ;
      cout << setw(13) << rcond(Ht);
      cout << setw(13) << MSE_after_weight_change ;
      cout << setw(13) << nlin << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << lambda ;
      report << setw(13) << MSE_after_weight_change;
      report << setw(13) << (d*(H*d)).sum() ;
      report << setw(13) << nlin << endl;
      k++;
     } //while
      report.close();
     return(return_vals);
  return(return_vals);
}
vector<double> MLP::train_OWO_Newton(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int k=0;
   size_t nlin;
   double Eg;
   double alpha,MSE;
   int Nu;
   double MSE_after_weight_change; 
   const double dmean = 0.5;
   const double dstd = 1.;
   vector<double> return_vals;
   int M=outputs,N=inputs,Nh=hidden;
   Nu=N+Nh+1;
   matrix Wa(Nu,M);
 
   int Niw=(N+1)*Nh, Now=M*(N+Nh+1),Nw=Niw+Now,iters;
   valarray<double> g(Niw),d(Niw),dt(Nw),gt(Nw);
   matrix D(Nh,N+1),L(Now,Now);

   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);

   while(k < Nit)
   {
       MSE=owo();
       return_vals.push_back(MSE);
       //  Wa=get_output_weights();
       //double alpha2=J(Wa);
       compute_derivs(true,Eg);
       g=unroll(G);
       d=compute_newton_direction(H,g,nlin);
       
       //L=chol(Ho);
      // for (i=0;i<Ht.getRow();i++)
        //  Ht[i][i] += t;
       //dt=compute_newton_direction(Ht,gt,nlin);
      // L=chol(Ht);

       //d_o=compute_newton_direction(Ho,g_out,nlin);
       D=reshape(d,Nh,N+1);
       alpha=linesearch(0,1,D,1e-3,iters);
      // alpha=linesearch(d);
       //alpha = linesearch(D);
//       if (k==0)
//       {
//         SVDRECORD T;
//         T=svd(H);
//         ofstream linefile("singulars.txt"),searchfile("linesearch.txt");
//         linefile << R << endl << endl;
//         for (int z =0;z<T.S.size();z++)
//             linefile << T.S[z] << endl;
//         for (double z=0;z<1.5;z+=.01)
//            searchfile << z << " " << J2(W+z*D) << endl;
//         linefile.close();
//         searchfile.close();
//       }
 //      alpha1=linesearch(d_o);
      MSE_after_weight_change=J2(W+alpha*D);
       //if (MSE_after_weight_change > MSE ) 
     //for some reason the objective
     //function is not decrementing.  Perhaps we have gone too far in the
     //line search, so we cut alpha down to see
      // {
      //   alpha/=2;
      //   cout << J2(W+alpha*D) << endl;
      // }
      update_weights(alpha,D);
      //MSE=J2(W);
      cout << setw(4)  << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << alpha ;
      //cout << setw(13) << MSE_after_weight_change;
      cout << setw(13) << rcond(H);
   //   cout << setw(13) << MSE_after_weight_change;
      cout << setw(13) << (d*(H*d)).sum() ;
      cout << setw(13) << nlin << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << alpha ;
      //report << setw(13) << MSE;
      report << setw(13) << MSE_after_weight_change;
      report << setw(13) << (d*(H*d)).sum() ;
      report << setw(13) << nlin << endl;
      k++;
     } //while
      report.close();
     return(return_vals);
}
vector<double> MLP::train_OIG_BP(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int i,k=0;
   size_t nlin=0;
   double Eg;
   double alpha=0,MSE;
   double alpha1=0; 
   const double dmean = 0.5;
   const double dstd = 1.;
   vector<double> return_vals;
   int N=inputs,Nh=hidden;//,num_dependent_rows;
   matrix D(Nh,N+1);
 
   valarray<double> r(N+1);
   vector<int> dependent_rows;
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   calculate_net_stats();
   while(k < Nit)
   {
     MSE=owo();
     return_vals.push_back(MSE);
      compute_derivs(true,Eg);
      D=G;
      compute_derivs_OIG(*this,D);
      ols_linear_solve(r,H_ig,g_ig,dependent_rows);
      for (i=0;i<Nh;i++)
          D.scale_row(i,r);
      W +=  D;
   
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << alpha ;
      cout << setw(13) << alpha1;
      cout << setw(13) << dependent_rows.size() << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << alpha ;

      report << setw(13) << alpha1;
      report << setw(13) << nlin << endl;

      k++;
     } //while
                
     report.close();
     return(return_vals);
}

vector<double> MLP::train_OIG_HWO(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int i,k=0;
   size_t nlin=0;
   double Eg;
   double alpha=0,MSE;
   double alpha1=0; 
   const double dmean = 0.5;
   const double dstd = 1.;

   int N=inputs,Nh=hidden;//,num_dependent_rows;
   matrix G_hwo(Nh,N+1),D(Nh,N+1);
 
   valarray<double> r(N+1);
   vector<double>return_vals;
   vector<int> dependent_rows;
     init_weights();
     calculate_net_stats();
     net_control(dmean,dstd);
     calculate_net_stats();
    while(k < Nit)
    {
      MSE=owo();
      return_vals.push_back(MSE);
      compute_derivs(true,Eg);
      G_hwo=compute_HWO_direction(Ri,G,nlin);
      D=G_hwo;
      compute_derivs_OIG(*this,D);
      ols_linear_solve(r,H_ig,g_ig,dependent_rows);
      for (i=0;i<Nh;i++)
          D.scale_row(i,r);
        W +=  D;
   
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << alpha ;
      cout << setw(13) << alpha1;
      cout << setw(13) << dependent_rows.size() << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << alpha ;

      report << setw(13) << alpha1;
      report << setw(13) << nlin << endl;

      k++;
     } //while
                
     report.close();
     return(return_vals);
}


