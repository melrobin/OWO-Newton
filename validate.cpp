#include <string>
//#include <iomanip>
#include "utils.h"
#include <valarray>
//#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include "mlp.h"
//#include <cassert>
#ifndef MYDEFS_H
#define MYDEFS_H
#include "mydefs.h"
#endif
#include "matops.h"
#include "fortran.h"

double validate(MLP* X,string val_fname)
{
   double E=0;
   int N=X->N,Nh=X->Nh,M=X->M,i,num_patterns=0;
   valarray<double> x(N+1),O(Nh),net(Nh),y(M),e(M),t(M);
   ifstream tstfile(val_fname.c_str());
   while (!tstfile.eof())
   {
      for (i=0;i<N;i++)
        tstfile >> x[i];

      x[N]=1.;
      if (!tstfile.eof())
        for (i = 0; i < M; i++)
            tstfile >> t[i];
      else
        break;
      num_patterns++;
      y=X->process_pattern(x,net,O);
        e=t-y;

     E += (e*e).sum();
    }
   tstfile.close();
   cout << "Validation patterns processed: " << num_patterns << endl;
   return(E/num_patterns);
}

