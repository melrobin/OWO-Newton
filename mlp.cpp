
#include <cstring>
#include <cstdlib>
#include <math.h>
#include <string>
#include <iomanip>
#include "utils.h"
#include <valarray>
#include <iostream>
#include <fstream>
#include "mlp.h"
#ifdef _WIN32
  #include <mkl.h>
  #define F77_FUNC(func) (func)
#else
  #include "fortran.h"
  #define F77_FUNC(func) func##_
#endif


#include "matops.h"
#ifndef MYDEFS_H
#define MYDEFS_H
#include "mydefs.h"
#endif

double h(double x)
{
        return(1./(1+exp(-x)));
}

using namespace std;
extern double slete(double, double);
extern double slete1(double, double );

void MLP::init_weights()
{
   for (int i=0;i<Nh;i++)
   {
      W[i][N]=slete(0,1);
      for (int j=0;j<N;j++)
          W[i][j]=slete(0,1)/Xstd[j];
   }

}

void MLP::set_W(const matrix &X)
{
  if ((X.getRow() != W.getRow())|| (X.getCol() != W.getCol()))
  {
     cout << "set_W:  wrong dimensions of matrix" << endl;
     exit(0);
  }
  W=X;
}

void MLP::set_Woh(const matrix &X)
{
  if ((X.getRow() != Woh.getRow())|| (X.getCol() != Woh.getCol()))
  {
     cout << "set_W:  wrong dimensions of matrix" << endl;
     exit(0);
  }
  Woh=X;
}

void MLP::set_Woi(const matrix &X)
{
  if ((X.getRow() != Woi.getRow())|| (X.getCol() != Woi.getCol()))
  {
     cout << "set_Woi:  wrong dimensions of matrix" << endl;
     exit(0);
  }
  Woi=X;
}

void MLP::write_weights(string output)
{
   ofstream outfile(output.c_str());
   outfile << "500" << endl;
   outfile << trnFile << endl;
   outfile << N << endl;
   outfile << M << endl;
   outfile << Nh << endl;
   outfile << "The_training_algorithm_is_OWO-Newton" << endl;
   outfile << transpose(Woi);
   outfile << transpose(Woh);
   outfile << W;
   outfile.close();
}
double MLP::J2(const matrix& Wa) 
//The objective function is computed by making one pass through the data with
//the current weights.
//ASSUMPTIONS:  The input mean has been incorporated into the thresholds 
{
   ifstream trnfile,infile;
   double ZERO=0.0,ONE=1.0;
   valarray<double> e(M);
   int i,NPlusOne,num_patterns=0;
   int INTONE=1;
   valarray<double> x(N+1),t(M),y(M),net(Nh),temp(M);
 
   NPlusOne = N + 1;
   trnfile.open(trnFile.c_str());
   while (!trnfile.eof())
   {
      for (i=0;i<N;i++)
        trnfile >> x[i];
 
      x[N]=1.0;
      if (!trnfile.eof())
         for (i=0;i<M;i++)
             trnfile >> t[i];
      else
        break;
      num_patterns++; 
      F77_FUNC(dgemv)("T",&NPlusOne,&Nh,&ONE,&Wa[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
      net=net.apply(h);
                
      F77_FUNC(dgemv)("T",&Nh,&M,&ONE,&Woh[0][0],&Nh,&net[0],&INTONE,&ZERO,&y[0],&INTONE);
      F77_FUNC(dgemv)("T",&NPlusOne,&M,&ONE,&Woi[0][0],&NPlusOne,&x[0],&INTONE,&ONE,&y[0],&INTONE);
                
      temp = t-y;
      e += (temp*temp);
   }
   trnfile.close();
   return(e.sum()/num_patterns);
}
double MLP::J(const string fname) 
//The objective function is computed by making one pass through the data with
//the current weights.
//ASSUMPTIONS:  The input mean has been incorporated into the thresholds 
{
   ifstream trnfile,infile;
   double ZERO=0.0,ONE=1.0;
   valarray<double> e(M);
   int i,NPlusOne,num_patterns=0;
   int INTONE=1;
   valarray<double> x(N+1),t(M),y(M),net(Nh),temp(M);
 
   NPlusOne = N + 1;
   infile.open(fname.c_str());
   while (!infile.eof())
   {
      for (i=0;i<N;i++)
        infile >> x[i];
 
      x[N]=1.0;
      if (!infile.eof())
         for (i=0;i<M;i++)
             infile >> t[i];
      else
        break;
      num_patterns++; 
      F77_FUNC(dgemv)("T",&NPlusOne,&Nh,&ONE,&W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
      net=net.apply(h);
                
      F77_FUNC(dgemv)("T",&Nh,&M,&ONE,&Woh[0][0],&Nh,&net[0],&INTONE,&ZERO,&y[0],&INTONE);
      F77_FUNC(dgemv)("T",&NPlusOne,&M,&ONE,&Woi[0][0],&NPlusOne,&x[0],&INTONE,&ONE,&y[0],&INTONE);
                
      temp = t-y;
      e += (temp*temp);
   }
   infile.close();
   return(e.sum()/num_patterns);
}
double MLP::J2(const matrix& Wa,const matrix& W1,const matrix& W2) 
//The objective function is computed by making one pass through the data with
//the current weights.
//ASSUMPTIONS:  The input mean has been incorporated into the thresholds 
{
   ifstream trnfile,infile;
   double ZERO=0.0,ONE=1.0;
   valarray<double> e(M);
   int i,NPlusOne,num_patterns=0;
   int INTONE=1;
   valarray<double> x(N+1),t(M),y(M),net(Nh),temp(M);
 
   NPlusOne = N + 1;
   trnfile.open(trnFile.c_str());
   while (!trnfile.eof())
   {
      for (i=0;i<N;i++)
        trnfile >> x[i];
 
      x[N]=1.0;
      if (!trnfile.eof())
         for (i=0;i<M;i++)
             trnfile >> t[i];
      else
        break;
      num_patterns++; 
      F77_FUNC(dgemv)("T",&NPlusOne,&Nh,&ONE,&Wa[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
      net=net.apply(h);
                
      F77_FUNC(dgemv)("T",&Nh,&M,&ONE,&W2[0][0],&Nh,&net[0],&INTONE,&ZERO,&y[0],&INTONE);
      F77_FUNC(dgemv)("T",&NPlusOne,&M,&ONE,&W1[0][0],&NPlusOne,&x[0],&INTONE,&ONE,&y[0],&INTONE);
                
      temp = t-y;
      e += (temp*temp);
   }
   trnfile.close();
   return(e.sum()/num_patterns);
}
double MLP::J2(const matrix& W1,const matrix& W2) 
//The objective function is computed by making one pass through the data with
//the current weights.
//ASSUMPTIONS:  The input mean has been incorporated into the thresholds 
{
   ifstream trnfile,infile;
   double ZERO=0.0,ONE=1.0;
   valarray<double> e(M);
   int i,NPlusOne,num_patterns=0;
   int INTONE=1;
   valarray<double> x(N+1),t(M),y(M),net(Nh),temp(M);
 
   NPlusOne = N + 1;
   trnfile.open(trnFile.c_str());
   while (!trnfile.eof())
   {
      for (i=0;i<N;i++)
        trnfile >> x[i];
 
      x[N]=1.0;
      if (!trnfile.eof())
         for (i=0;i<M;i++)
             trnfile >> t[i];
      else
        break;
      num_patterns++; 
      F77_FUNC(dgemv)("T",&NPlusOne,&Nh,&ONE,&W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
      net=net.apply(h);
                
      F77_FUNC(dgemv)("T",&Nh,&M,&ONE,&W2[0][0],&Nh,&net[0],&INTONE,&ZERO,&y[0],&INTONE);
      F77_FUNC(dgemv)("T",&NPlusOne,&M,&ONE,&W1[0][0],&NPlusOne,&x[0],&INTONE,&ONE,&y[0],&INTONE);
                
      temp = t-y;
      e += (temp*temp);
   }
   trnfile.close();
   return(e.sum()/num_patterns);
}
void MLP::spit_info()
{
   cout << "Class MethodTraining file: " << trnFile.c_str() << endl;
   cout << "Number of patterns: " << Nv << endl;
   cout << "Number of hidden units: " << Nh << endl;
   cout << "Number of inputs: " << N << endl;
   cout << "Number of outputs: " << M << endl << endl;
}

void MLP::update_weights(double alpha,const matrix& D)
{
    W += alpha*D;
}
void MLP::save_weights()
{
//    memcpy(&W_old[0][0],&W[0][0],(N+1)*Nh*sizeof(double));
}

void MLP::restore_weights()
{
  // memcpy(&W[0][0],&W_old[0][0],(N+1)*Nh*sizeof(double));

}

int MLP::getNumOutputs(void)
{
	return (N);
}
MLP::MLP(int N1,int N2,int N3,string fname):net_std(N3),Et(N2)
{
   N=N1;
   Nh=N3;
   M=N2;
   Nv=0;
   Nu = N + Nh + 1;
   Xstd=valarray<double>(N1+1);
   net_mean=valarray<double>(N3);
   //H=matrix(Nh*(N+1),Nh*(N+1));
   H_OIT=matrix((N+1)*(N+1),(N+1)*(N+1));
   //H_ONT=matrix(Nh*Nh,Nh*Nh);
   G=matrix(Nh,N+1);
   G1=matrix(Nh,N+1);
   R=matrix(Nu,Nu);
   Ri=matrix(N+1,N+1);
   C=matrix(M,Nu);
   W=matrix(Nh,N+1);
   Woi=matrix(M,N+1);
   Woh=matrix(M,Nh);
   H_ig=matrix(N+1,N+1);
   g_ig=valarray<double>(N+1);
   Gt=valarray<double>((N+1)*(Nh+M)+M*Nh);
   G_oit=matrix(N+1,N+1);
   Ht=matrix((N+1)*(Nh+M)+M*Nh,(N+1)*(Nh+M)+M*Nh);
   trnFile=fname;
}
MLP::~MLP()
{
}

size_t MLP::init_mlp(string trnFile)
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
  this->trnFile=trnFile;
  if (!infile.is_open())
  {
     cout << "Error reading file in stats routine" << endl;
     exit(0);
  }
  Nv=0;
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

matrix MLP::get_gradient(void)
{
  size_t M=G.getRow(),N=G.getCol();
  matrix temp(M,N);
  temp = G;
  return(temp);
}

matrix MLP::get_output_weights(void)
{
   matrix temp;
   temp=vertcat(true,true,Woi,Woh);
   return(temp);
}
matrix MLP::get_input_weights(void)
{
  matrix temp(Nh,N+1);
  temp=W;
  return(temp);
}
double MLP::J(const matrix &Wa)
{
  matrix weights(Nh,Nu);
  matrix Wb;
  valarray<double> E(M),temp1(M),temp2(M);
  Wb=vertcat(true,true,Woi,Woh);
  temp1=Wb.extract_row(N);
  for (int i=N+1;i<Nu;i++) 
   {
      temp2=Wb.extract_row(i);
      Wb.insert_row(i-1,temp2);
   }
  Wb.insert_row(N+Nh,temp1);
  E=Et+diag(transpose(Wb)*R*Wb-2.*C*Wb);
  return(E.sum());
}
double act(double x)
{
	return(1./(1+exp(-x)));
}
valarray<double> MLP::process_pattern(valarray<double>& x,valarray<double>& net,valarray<double>& O)//,const matrix &D)
{
   valarray<double> y(M);//net(Nh),y(M),O(Nh);
   int INTONE=1;
   double ZERO = 0.0,ONE=1.0;
   int NPlusOne=N+1;
   matrix tempW(Nh,N+1);

   tempW=W;
 //  for (int i=0;i<Nh;i++)
 //    net[i]=O[i]=0;

   dgemv_("T",&NPlusOne,&Nh,&ONE,&tempW[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
 
     O=net.apply(act);
 
                //Part 1 of compute the output variable y
     dgemv_("T",&Nh,&M,&ONE,&Woh[0][0],&Nh,&O[0],&INTONE,&ZERO,&y[0],&INTONE);
                //Part 2 which adds the bypass weights
     dgemv_("T",&NPlusOne,&M,&ONE,&Woi[0][0],&NPlusOne,&x[0],&INTONE,&ONE,&y[0], &INTONE);
 
   return y;

}
valarray<double> MLP::process_pattern(valarray<double>& x)
{
   valarray<double> net(Nh),y(M),O(Nh);
   int INTONE=1;
   double ZERO = 0.0,ONE=1.0;
   int NPlusOne=N+1;

   dgemv_("T",&NPlusOne,&Nh,&ONE,&W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
 
     O=net.apply(act);
 
                //Part 1 of compute the output variable y
     dgemv_("T",&Nh,&M,&ONE,&Woh[0][0],&Nh,&O[0],&INTONE,&ZERO,&y[0],&INTONE);
                //Part 2 which adds the bypass weights
     dgemv_("T",&NPlusOne,&M,&ONE,&Woi[0][0],&NPlusOne,&x[0],&INTONE,&ONE,&y[0], &INTONE);
 
   return y;

}

