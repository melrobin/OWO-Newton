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
#include <sstream>
#include <QtGui/QtGui>
#include <unistd.h>
#ifndef MYDEFS_H
#define MYDEFS_H
#include "mydefs.h"
#endif
#include "matops.h"
#include "numap_ui.h"
using namespace std;
extern "C" double genrand_real3(void);
extern vector<int> randperm(int);
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
   vector<double> results;
   vector< vector<double> > total_results;
   int Nit,M,N;
   int c;
   opterr = 0;

   if (argc < 2)
   {
      cout <<"Not enough arguments.  You must specify a configuration file...aborting! "<< endl;
      exit(0);
   }

//   QApplication app(argc, argv);
//   QMainWindow *window=new QMainWindow;
//   Ui_NuMap numap;
//   numap.setupUi(window);
  
   //  window->resize(320, 240);
//     window->show();
    // window.show();
    // window.setWindowTitle(
      //   QApplication::translate("toplevel", "Top-level widget"));
//   app.exec();
   setup(argc,argv,config);
   while ((c=getopt(argc,argv,"in:"))!= -1)
      switch(c)
        {
           case 'i':
              cout << optarg << endl;
              break;
           case 'n':
              cout << "Inputs: " << optarg;
              N = atoi(optarg);
           break;
           default:
              abort();
              break;
        }     
   ifstream infile(config.trainingfile.c_str());
   if (!infile.is_open()) //bad training file
   {
      cout << "Bad training file...try again!" << endl;
      exit(0);
   }
   infile.close();
   Nit = config.iterations;
   //N = config.inputs;
   report.open(config.report_file.c_str());
   cout << "Starting MLP Training" << endl;  
   M = config.outputs;
   
   MLP *myMLP = new MLP(N,M,config.hidden_units,config.trainingfile);
   MLP *mlp2 = new MLP(N,M,config.hidden_units,config.trainingfile);

 //  size_t num_patterns = myMLP->calculate_stats();
   size_t num_patterns = mlp2->calculate_stats();
   cout << "MLP Training Report" << endl;
   cout << "Training file: " << config.trainingfile << endl;
   cout << "Inputs: " << config.inputs << endl;
   cout << "Outputs: " << config.outputs << endl;
   cout << "Hidden Units: " << config.hidden_units << endl;
   cout << "Number of Patterns: " << num_patterns << endl;
   cout << "Number of iterations: " << Nit << endl;
   report << "MLP Training Report" << endl;
   report << "Training file: " << config.trainingfile << endl;
   report << "Inputs: " << config.inputs << endl;
   report << "Outputs: " << config.outputs << endl;
   report << "Hidden Units: " << config.hidden_units << endl;
   report << "Number of Patterns: " << num_patterns << endl;
   report << "Number of iterations: "<< Nit << endl;
   cout<<setfill('_')<<setw(80)<<"_"<<endl;
   cout<<setfill(' ');
                
   cout <<setw(9) << "Iteration";
   cout << setw(11) << "MSE";
   cout << setw(9) << "Z" ;
   cout << setw(15) << "Metric 1" ;
   cout << setw(13) << "Metric 2";
   cout << setw(16) << "Rank" << endl;
   //report <<setw(9) << "Iteration";
   //report << setw(11) << "MSE";
   //report << setw(9) << "Z" ;
   //report << setw(15) << "Metric 3" ;
   //report << setw(13) << "d1";
   //report << setw(16) << "Rank" << endl;

//    init_genrand(448247);
//   tstart=tic();
//   results=myMLP->train_ONT(Nit,config.inputs,M,config.hidden_units,config.trainingfile);
//    tstop=toc();
//    cout << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   cout << "Percentage error reduction: " << (results.front()-results.back())/results.front() * 100. << endl;
//   report << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//  total_results.push_back(results); 
// if (config.optimal_learning_factor)
//{
//   cout << "ONT Validation error: "<< validate(myMLP,config.val_file) << endl;
//report << "ONT Validation error: "<< validate(myMLP,config.val_file) << endl;
//}
 //   init_genrand(448247);
 //  tstart=tic();
 //  results=myMLP->train_LM(Nit,config.inputs,M,config.hidden_units,config.trainingfile);
 //  total_results.push_back(results);
 //  tstop=toc(); 

/*    init_genrand(448247);
//   tstart=tic();
//   results=myMLP->train_OIT_HWO(Nit,config.inputs,M,config.hidden_units,config.trainingfile);
//total_results.push_back(results);
   tstop=toc(); */
   results.clear();;
   init_genrand(448247);          
   tstart=tic();
   results=mlp2->train_OWO_Newton(Nit,config.inputs,M,config.hidden_units,config.trainingfile);
   tstop=toc();
   cout << "Time elapsed: " << tstop-tstart << " seconds" << endl;
   cout << "Percentage error reduction: " << (results.front()-results.back())/results.front() * 100. << endl;
   report << "Time elapsed: " << tstop-tstart << " seconds" << endl;
 //  total_results.push_back(results);
 //   vector<int> trains;
 //  valarray<double> myResults(Nit);
 //  if (config.optimal_learning_factor)
 //  {
  //    trains=randperm(config.k);
   //   for (int i=0;i<config.k;i++)
    //  {
     //    ostringstream ss;
     //    ss << "file" << trains[i]<<".txt";
     //    mlp2->init_mlp(ss.str());
     //    results = mlp2->train_OWO_Newton(Nit,config.inputs,M,config.hidden_units,ss.str());
      //   for (int ii=0;ii<Nit;ii++)
       //     myResults[ii] += results[ii];
     // }
   //   myResults /= config.k;
  //    for (int ii=0;ii<Nit;ii++)
   //        cout << ii+1 << " " <<myResults[ii]<<endl;
   //  total_results.push_back(myResults/config.k);
 //     cout << "OWO-Newton Validation error: "<< validate(mlp2,config.val_file) << endl;
   //   report << "OWO-Newton Validation error: "<< validate(mlp2,config.val_file) << endl;
  // }
 mlp2->write_weights(config.weight_file); 
 /*  init_genrand(448247);
   tstart=tic();
 //  results = myMLP->train_OIG_BP(Nit,N,M,config.hidden_units,config.trainingfile); 
//   tstop=toc();
//   cout << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   cout << "Percentage error reduction: " << (results.front()-results.back())/results.front() * 100. << endl;
//   report << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   total_results.push_back(results);
*/
//   init_genrand(448247);
//  MLP *mlp3=new MLP(N,M,config.hidden_units,config.trainingfile);
//   tstart=tic();
//   mlp3->calculate_stats();
//   results =  mlp3->train_OIT_HWO(Nit,N,M,config.hidden_units,config.trainingfile);     
//   tstop=toc();
//   cout << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   cout << "Percentage error reduction: " << (results.front()-results.back())/results.front() * 100. << endl << endl;
//   report << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   total_results.push_back(results);
//    if (config.optimal_learning_factor)
//    {
//        cout << "OIT Validation error: "<< validate(mlp3,config.val_file) << endl;
//        report << "OIT Validation error: "<< validate(mlp3,config.val_file) << endl;
//    }
 //  MLP *mlp3=new MLP(N,M,config.hidden_units,config.trainingfile);//,*mlp2 = new MLP(N,M,config.hidden_units,config.trainingfile);;
 //  init_genrand(448247);
 //  mlp3->calculate_stats();
 //  tstart=tic();
 //  results = mlp3->train_BP(Nit,N,M,config.hidden_units,config.trainingfile);
//   tstop=toc();
//   cout << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   cout << "Percentage error reduction: " << (results.front()-results.back())/results.front() * 100. << endl ;
//   report << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   total_results.push_back(results);
//   results.clear(); 
//    if (config.optimal_learning_factor)
//    {
//        cout << "BP Validation error: "<< validate(mlp3,config.val_file) << endl;
//        report << "BP Validation error: "<< validate(mlp3,config.val_file) << endl//;
  //  }

//   init_genrand(448247);
//   tstart=tic();
//   results=myMLP->train_OWO_HWO(Nit,config.inputs,M,config.hidden_units,config.trainingfile);
//   tstop=toc();
//   cout << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   cout << "Percentage error reduction: " << (results.front()-results.back())/results.front() * 100. << endl;
//   report << "Time elapsed: " << tstop-tstart << " seconds" << endl;
//   total_results.push_back(results);

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
//   results.clear();

   //myMLP->write_weights();
   //string fname("weights.txt");
 //       read_weights(fname,Woi,Woh,W);
   //     MLP *testMLP = new MLP(config.inputs,config.outputs,config.hidden_units,config.trainingfile);
     //   testMLP->set_W(W);
      //  testMLP->set_Woh(Woh);
       // testMLP->set_Woi(Woi);
        //cout << testMLP->J2(W) << endl;
//   int rows=total_results.size();
 //  int cols=total_results[0].size();
 //  ofstream outfile1(config.chart_file.c_str());
 //  for (int i=0;i < cols; i++)
 //    {
 //     for (int j=0; j < rows;j++)
 //        outfile1 << total_results[j][i] << "   ";
 //     outfile1 << endl;
 //    }
//   outfile1.close();  
   report.close();
        return(0);
}

