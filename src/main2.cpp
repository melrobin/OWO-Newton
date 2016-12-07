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
#include <QtGui/QtGui>
#ifndef MYDEFS_H
#define MYDEFS_H
#include "mydefs.h"
#endif
#include "matops.h"
#include "numap_ui.h"
using namespace std;

extern double tic(void);
extern double toc(void);
extern void setup(int,char **,CONFIG& );

extern "C" void init_genrand(unsigned long );
valarray<double> compute_newton_direction(const matrix ,valarray<double>,size_t & );

void read_weights(string fname,matrix& Woi, matrix& Woh, matrix& W)
{
   ifstream infile(fname.c_str());
   size_t N,M,Nh,r;
   string method,version,trnfile;

   infile >> version;
   infile >> trnfile;
   infile >> N;
   infile >> M;
   infile >> Nh;
   infile >> method;
   matrix tempW(Nh,N+1),tempWoi(M,N+1),tempWoh(M,Nh);
   valarray<double> woi(M*(N+1)),woh(M*Nh),w(Nh*(N+1));
   for (r=0; r < (N+1)*M; r++)
        infile >> woi[r];

   for (r=0; r < Nh*M; r++)
        infile >> woh[r];

   for (r=0; r < Nh*(N+1); r++)
        infile >> w[r];
   
   Woi=reshape(woi,M,N+1);
   Woh=reshape(woh,M,Nh);
   W=reshape(w,Nh,N+1); 
   infile.close();
}
int main(int argc,char *argv[])
{
   CONFIG config;
   ofstream report;
   double tstart,tstop;
   vector<double> results, comparison,results2,comparison2;
   vector< vector<double> > total_results;
   int Nit,M,N,Nh,num_patterns;

   if (argc < 2)
   {
      cout <<"Not enough arguments.  You must specify a configuration file...aborting! "<< endl;
      exit(0);
   }
   setup(argc,argv,config);
   ifstream infile(config.trainingfile.c_str());
   if (!infile.is_open()) //bad training file
   {
      cout << "Bad training file...try again!" << endl;
      exit(0);
   }
   infile.close();
   Nit = config.iterations;
   N = config.inputs;
   report.open(config.report_file.c_str());
   cout << "Starting MLP Training" << endl;  
   M = config.outputs;
   
  // mlp2->calculate_stats();
for (Nh=1; Nh<15;Nh++)
{
   MLP *mlp2 = new MLP(N,M,Nh,config.trainingfile);
   MLP *myMLP = new MLP(N,M,Nh,config.trainingfile);
   cout << "MLP Training Report" << endl;
   cout << "Training file: " << config.trainingfile << endl;
   cout << "Inputs: " << config.inputs << endl;
   cout << "Outputs: " << config.outputs << endl;
   cout << "Hidden Units: " << Nh << endl;
   cout << "Number of Patterns: " << num_patterns << endl;
   cout << "Number of iterations: " << Nit << endl;
   cout<<setfill('_')<<setw(80)<<"_"<<endl;
   cout<<setfill(' ');
                
   cout <<setw(9) << "Iteration";
   cout << setw(11) << "MSE";
   cout << setw(9) << "Z" ;
   cout << setw(15) << "Metric 1" ;
   cout << setw(13) << "Metric 2";
   cout << setw(16) << "Rank" << endl;

    init_genrand(448247);
//    myMLP->set_Nh(Nh);
   num_patterns = myMLP->calculate_stats();
  
   results=myMLP->train_ONT(Nit,config.inputs,M,Nh,config.trainingfile);
init_genrand(448247);
 mlp2->calculate_stats();
   results2=mlp2->train_OIT_HWO(Nit,config.inputs,M,Nh,config.trainingfile);
   //myMLP->set_Nh(Nh);
  // cout << results.back();
   comparison.push_back(results.back());
   comparison2.push_back(results2.back());
  //total_results.push_back(results); 
 //if (config.optimal_learning_factor)
  //{
  // cout << "ONT Validation error: "<< validate(myMLP,config.val_file) << endl;
    //report << "ONT Validation error: "<< validate(myMLP,config.val_file) << endl;
  }
cout << endl;
 for (int i=0;i<comparison.size();i++)
   cout << comparison[i] << " " << comparison2[i] << endl;

 //  init_genrand(448247);
//MLP *mlp3=new MLP(N,M,config.hidden_units,config.trainingfile);
//   tstart=tic();
//mlp3->calculate_stats();
  // results =  mlp3->train_OIT_HWO(Nit,N,M,config.hidden_units,config.trainingfile);     
  // tstop=toc();
  // cout << "Time elapsed: " << tstop-tstart << " seconds" << endl;
  // cout << "Percentage error reduction: " << (results.front()-results.back())/results.front() * 100. << endl << endl;
  // report << "Time elapsed: " << tstop-tstart << " seconds" << endl;
  // total_results.push_back(results);
 //   if (config.optimal_learning_factor)
  //  {
  //      cout << "OIT Validation error: "<< validate(mlp3,config.val_file) << endl;
  //      report << "OIT Validation error: "<< validate(mlp3,config.val_file) << endl;
  //  }

//   init_genrand(448247);
//   tstart=tic();
//   results=mlp2->train_CG(Nit,config.inputs,M,config.hidden_units,config.trainingfile);
//   tstop=toc();
//   cout << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   cout << "Percentage error reduction: " << (results.front()-results.back())/results.front() * 100. << endl;
//   report << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   total_results.push_back(results);
//   if (config.optimal_learning_factor)
//   {
//       cout << "CG Validation error: "<< validate(mlp2,config.val_file) << endl;
//       report << "CG Validation error: "<< validate(mlp2,config.val_file) << endl;
//   }

 //  int rows=total_results.size();
 //  int cols=total_results[0].size();
 //  ofstream outfile1(config.chart_file.c_str());
 //  for (int i=0;i < cols; i++)
 //    {
  //    for (int j=0; j < rows;j++)
    //     outfile1 << total_results[j][i] << "   ";
    //  outfile1 << endl;
    // }
 //  outfile1.close();  
 //  report.close();
        return(0);
}

