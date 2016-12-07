#include <string>
#include <iomanip>
#include "utils.h"
#include <valarray>
#include <iostream>
#include <fstream>
#include "mlp.h"
#ifndef MYDEFS_H
#define MYDEFS_H
#include "mydefs.h"
#endif
#include "matops.h"
using namespace std;

extern "C" void init_genrand(unsigned long );

int main(int argc,char *argv[])
{
   
   int Nit,M,N,Nh,seed;
   double MSE;

   if (argc != 7)
   {
      cout <<"usage: gen_init_wts training_file outputfile random_seed inputs outputs hidden_units"<< endl;
      exit(0);
   }

   seed=atoi(argv[3]);
   init_genrand(seed);
   N=atoi(argv[4]);
   M=atoi(argv[5]);
   Nh=atoi(argv[6]);
   cout << "Generating initial weight file" << endl;  
   
   MLP *myMLP = new MLP(N,M,Nh,argv[1]);
           myMLP->calculate_stats();
   myMLP->init_weights();
   myMLP->calculate_net_stats();
   myMLP->net_control(.5,1);
  
       MSE=myMLP->owo();
 
  myMLP->write_weights(argv[2]);
  cout << MSE << endl;
        return(0);
}

