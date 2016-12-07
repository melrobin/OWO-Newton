#include <string>
#include <iomanip>
#include "utils.h"
#include <valarray>
#include <ctime>
#include <iostream>
#include <fstream>
#include "mlp.h"
#include <cassert>
#include <vector>
#ifndef MYDEFS_H
#define MYDEFS_H
#include "mydefs.h"
#endif
#include "matops.h"
using namespace std;

extern double tic(void);
extern double toc(void);
extern void setup(int,char **,CONFIG& );

extern "C" void init_genrand(unsigned long );
valarray<double> compute_newton_direction(const matrix ,valarray<double>,size_t & );

void read_weights(int& N, int& M, int &Nh, string& trnfile,string fname,matrix& Woi, matrix& Woh, matrix& W)
{
   ifstream infile(fname.c_str());
   int r,i,j;
   string method,version;

   infile >> version;
   infile >> trnfile;
   infile >> N;
   infile >> M;
   infile >> Nh;
   infile >> method;

//   W=matrix(Nh,N+1);
  // Woh=matrix(M,Nh);
  // Woi=matrix(M,N+1);

  // matrix tempW(Nh,N+1),tempWoi(M,N+1),tempWoh(M,Nh);
   valarray<double> woi(M*(N+1)),woh(M*Nh),w(Nh*(N+1));
   for (r=0; r < (N+1)*M; r++)
        infile >> woi[r];

   for (r=0; r < Nh*M; r++)
        infile >> woh[r];

   for (r=0; r < Nh*(N+1); r++)
        infile >> w[r];
   
   Woi=reshape(woi,M,N+1);
   Woh=reshape(woh,M,Nh);
   r=0;
   W=matrix(Nh,N+1);
   for (i=0;i<Nh;i++)
     for (j=0;j<=N;j++)
       W[i][j]=w[r++];
   //W=reshape(w,Nh,N+1); 
   infile.close();
   //cout << W << endl;
}
int main(int argc,char *argv[])
{
   int Nh,N,M,i,num_patterns=0;
  
   if (argc < 3)
   {
      cout <<"Usage: process weight_file datafile"<< endl;
      exit(0);
   }
   //ifstream infile(config.trainingfile.c_str());
   //if (!infile.is_open()) //bad training file
   //{
     // cout << "Bad training file...try again!" << endl;
     // exit(0);
  // }
  // infile.close();
  
   string fname(argv[1]),trnfile,tstfile(argv[2]);
   ifstream infile(tstfile.c_str());
   ofstream outfile("outvals.txt");
   matrix W,Woh,Woi;
   double e;
  // W=matrix(Nh,N+1);
  // Woh=matrix(M,Nh);
  // Woi=matrix(M,N+1);
   //tstfile=argv[2].c_str();
   read_weights(N,M,Nh,trnfile,fname,Woi,Woh,W);
   //W(Nh,N+1),Woh(M,Nh),Woi(M,N+1);
   valarray<double> x(N+1),t(M),y(M),temp(M);
   MLP *testMLP = new MLP(N,M,Nh,trnfile);
   testMLP->set_W(W);
   testMLP->set_Woh(Woh);
   testMLP->set_Woi(Woi);
   //cout << testMLP->J(tstfile) << endl;
   
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
      y=testMLP->process_pattern(x);
                
      for (i=0;i<M;i++)
         outfile << t[i] << " ";
      for (i=0;i<M;i++)
         outfile << y[i] << " ";
      outfile << endl;
      temp=t-y;
      e += (temp*temp).sum();
   }
   cout << "Error: " << e/num_patterns << endl;
   //report.close();
   return(0);
}

