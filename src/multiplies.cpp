#include <iostream>
#include <fstream>
#include <sstream>
using namespace std;
extern void setup(int,string,CONFIG& );
int main(int argc,char *argv[])
{
   int N,M,Nh,Nv,Nu,Niw,Nw,N_OWO,N_OLS,N_Newton,N_OWO_NEWTON,N_OWO_BP,N_LM,Nit;
   CONFIG config;
   stringstream ss;
   string arg;
 
   if (argc < 3)
   {
      cout << "Usage: multiplies config_file output_file" << endl;
      exit(0);
   }
   ss.str(argv[1]);
   ss >> arg;
   //Pull information from configuration file and set up variables
   setup(argc,arg,config);
   Nit = config.iterations;
   N = config.inputs;
   Nv = config.num_patterns;  
   M = config.outputs;
   //Open up the output file for writing
   ofstream outfile(argv[2]);
   
   // Let's compute some preliminary variables to get our computations started 
   // These will be the basis for our higher order calculations
   Nu=N+Nh+1;
   Niw=(N+1)*Nh;
   Nw=Niw+M*Nh;
   // Now we start computations on the algorithms
   N_OWO=Nv*(Nu+1)*(M+Nu/2);
   N_OLS=Nu*(Nu+1)*(M+Nu*(2*Nu+1)/6+3/2);
   N_Newton=(1+M/2)*Nv*(N+1)*(N+2)+Nh*(2*M+1)*(N+1)*Nv;
   N_OWO_NEWTON=N_OWO+N_OLS+N_Newton;
   N_OWO_BP=Nv*(2*Nh*(N+2)+M*(Nu+1)+Nu*(Nu+1)/2+M*(N+6*Nh+4))+N_OLS+Nh*(N+1);
   N_LM=N_OWO_BP+Nv*(M*Nu*(Nu+3*Niw)+4*Niw*Niw)+Nw*Nw*Nw+Nw*Nw;
   outfile << "OWO_BP " << "OWO_Newton " << "LM " << endl;
   for (int i=1;i<=Nit;i++)
     outfile << i*N_OWO_BP << " " << i*N_OWO_NEWTON << " " << i*N_LM << endl;
   outfile.close();
   return(0);
}

