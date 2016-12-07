#include <valarray>
#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string>
#ifdef _WIN32
  #include <mkl.h>
  #define F77_FUNC(func) (func)
#else
  #include "fortran.h"
  #define F77_FUNC(func) func##_
#endif

#include "random.h"
#include "mlp.h"
#include "matrix.h"
using namespace std;

extern double h(double);
void MLP::net_control(double mean,double std)
{
   int i,j;
   for (i=0;i<Nh;i++)
     for (j=0;j<=N;j++)
        W[i][j] *= (std/net_std[i]);
                
   for (i=0;i<Nh;i++)
     W[i][N]= W[i][N] - net_mean[i]*std/net_std[i] + mean;   
   
}

size_t MLP::calculate_stats()
/* Calculates the mean, variance and total energy in the datafile */
//Inputs:
//N Number of inputs
//Nh Number of hidden units
//Nc Number of outputs
//OUTPUT:
//Nv: number of patterns that are read

//Xmean:  The mean vector for the training file
//Xstd:  The standard deviation vector for the training file
{
  int i;
  ifstream infile;
  valarray<double> X(N+1),Xmean(N+1),ty(M);
  
  infile.open(trnFile.c_str());

  if (!infile.is_open())
  {
     cout << "Error reading file in stats routine" << endl;
     exit(0);
  }
  
   while (!infile.eof())
     {
        for (i = 0; i < N; i++)
           infile >> X[i];

        X[N]=1.0;

        if (infile.eof())
           break;
        else
           for (i=0;i<M;i++)
                infile >> ty[i];
        Nv++;
        Xmean += X;
        Xstd += X*X;

        Et += ty*ty;
       }
        Xmean /= Nv;
        Xstd = Xstd / (double)(Nv) - Xmean*Xmean;
        Xstd = sqrt(Xstd);
        Et /= double(Nv);
        valarray<double> sub1(Nh),sub2(M);
        sub1=W*Xmean;
        sub2=Woi*Xmean;
        for (i=0;i<Nh;i++)
           W[i][N] -= sub1[i];
        for (i=0;i<M;i++)
           Woi[i][N] -= sub2[i];
        infile.close();
    return(Nv);
}
void MLP::calculate_net_stats()
//Inputs: N, the number of patterns
//                Nh, the number of hidden units
//                      M, the numbere of outputs
//              weights, the input weight array
//Outputs:  net_mean, the mean of the net function
//                      net_std, standard deviation of the net function
{
   int i,NPlusOne,intone=1;
   FILE *ifs;
   valarray<double> X(N+1),ty(M),net(Nh);
   double ZERO=0.0,ONE=1.0;

   NPlusOne = N + 1;

//The following lines are critical to ensure that we initialize all techniques from the same spot and should not be deleted because declaring a valarray automatically initializes to zero.
   for (i=0;i<Nh;i++)
   {
        net_mean[i]=0.;
        net_std[i]=0.;
   }
   ifs = fopen(trnFile.c_str(),"r");
   rewind(ifs);
   while (!feof(ifs))
    {
       for (i = 0; i < N && !feof(ifs); i++)
           fscanf(ifs, "%lf", &X[i]);

       if (feof(ifs))
           break;
       else
           for (i=0;i<M;i++)
                fscanf(ifs,"%lf",&ty[i]);
       X[N] = 1.0;
       dgemv_("T",&NPlusOne,&Nh,&ONE,&W[0][0],&NPlusOne,&X[0],
               &intone,&ZERO,&net[0],&intone);

           net_mean += net;
           net_std += net*net;
     }

   net_mean /= Nv;
   net_std = net_std/(double)Nv-net_mean*net_mean;
   net_std=sqrt(net_std);

}

void MLP::calculate_corr()
{
    int i,Nu,intone=1,NPlusOne,num_read=0;
    valarray<double> ty(M);
    double ZERO=0.0,one=1.0;
    ifstream infile,trnfile;
    ofstream outfile;
    trnfile.open(trnFile.c_str());
    Nu = N + Nh + 1;
    valarray<double> net(Nh);
    valarray<double> Xa(Nu),x(N+1);
   // valarray<double> Xmean_local(Nu);
    NPlusOne = N + 1;
    zero(R);
    zero(Ri);
    zero(C); 

    infile.open("weights.wts");
        
  //  for (i=0;i<N;i++)
  //    infile >> Xmean_local[i];
  //  infile.close();
    while (!trnfile.eof())
      {
          for (i = 0; i < N; i++)
           {
               trnfile >> Xa[i];
               x[i]=Xa[i];
           }
          if (trnfile.eof())  
            break;
                
          for (i=0;i<M;i++)
            trnfile >> ty[i];
        
          num_read++;
          x[N]=1.;
          Xa[N]=1.;
          F77_FUNC(dgemv)("T",&NPlusOne,&Nh,&one,&W[0][0],&NPlusOne,&Xa[0],&intone,&ZERO,&net[0],&intone);
                 
          net=net.apply(h);

          F77_FUNC(dcopy)(&Nh,&net[0],&intone,&Xa[N],&intone);

          Xa[N+Nh]=1.0;
                 
          F77_FUNC(dger)(&Nu,&Nu,&one,&Xa[0],&intone,&Xa[0],&intone,&R[0][0],&Nu);
          F77_FUNC(dger)(&NPlusOne,&NPlusOne,&one,&x[0],&intone,&x[0],&intone,&Ri[0][0],&NPlusOne);
          F77_FUNC(dger)(&Nu,&M,&one,&Xa[0],&intone,&ty[0],&intone,&C[0][0],&Nu);
        }
     R = R/Nv;
     Ri = Ri/Nv;
     C = C/Nv;              
   trnfile.close();
}

