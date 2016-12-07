# include <cstdlib>
#include <iostream>
#include <valarray>
#include <string>
#include <sstream>
#include <vector>
#include <fstream>
using namespace std;
extern int random_int(int,int);
extern "C" double genrand_real3(void);
typedef std::pair<double,int> mypair;
bool comparator (const mypair& left, const mypair& right)
{
  return (left.first < right.first);
}

int main(int argc,char *argv[])
{
   ifstream infile;
   int i,j,N,M,k,Nv=0;             //  The column numbers for each pattern
  vector<string> fnames;

	if(argc != 5) 
	{
	    cout << "Usage: splitk training_file num_inputs num_outputs k"<< endl;
		exit(1);
	}
	k=atof(argv[4]);
    N=atoi(argv[2]);
    M=atoi(argv[3]);
	infile.open(argv[1]);
      valarray<double> x(N),t(M);
      for (i=0;i<k;i++)
     {
        ostringstream ss;
        ss  << "file"<<i<<".txt";
        fnames.push_back(ss.str());
     }
     ofstream outfiles[k];
     for (i=0;i<fnames.size();i++)
     {
        outfiles[i].open(fnames[i].c_str());
    //    cout << fnames[i] << endl;
     }
     while (!infile.eof())
     {
         j = random_int(0,k-1);
           for (i=0;i<N;i++)
              infile >> x[i];
 
           //x[N]=1.;
           if (infile.eof())
             break;
           for (i=0;i<M;i++)
                 infile >> t[i];

          for (i=0;i<N;i++)
             outfiles[j] << x[i] << " ";
          for (i=0;i<M;i++)
             outfiles[j] << t[i] << " ";
          outfiles[j] << endl;   
           Nv++;
      }
  cout << "Number of patterns processed: "<< Nv << endl;
  for (i=0;i<k;i++)
     outfiles[i].close();
}                      /* End of main  */
     
