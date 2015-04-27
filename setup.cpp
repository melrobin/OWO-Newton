#include <iostream>
#include <fstream>
#include <string>
#include "mydefs.h"
using namespace std;
void setup(int argc,char *argv[],CONFIG &config)
{
  ifstream infile(argv[1]); //open the configuration file
  infile >> config.trainingfile;
  infile >> config.inputs;
  infile >> config.outputs;
  infile >> config.hidden_units;
  infile >> config.iterations;
  infile >> config.optimal_learning_factor;
  infile >> config.report_file;
  infile >> config.weight_file;
  infile >> config.k;
  infile >> config.chart_file;
  infile >> config.val_file;
  infile.close();
}
  
