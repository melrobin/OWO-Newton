#include <iostream>
#include "mlp.h"

using namespace std;
double  **MLP::AllocateMatrix(int m, int n)
{
        double  **dptr;
        int i;

        dptr = (double **)calloc(m, sizeof(double *));
        //dptr = new double[m*n];
        if(dptr == NULL)
        {
                cout << "Not enough memory to allocate buffer." << endl;
                cout << "Aborting......\n"<< endl;
                /* terminate program if out of memory   */
                exit(1);
        }
       // dptr[0]=(double *)calloc(m*n,sizeof(double));
        dptr[0]=new double[m*n];
        for(i = 1; i < m ; i++)
        {
                dptr[i]=dptr[0]+i*n;
                if(dptr[i] == NULL)
                {
                   cout <<  "Not enough memory to allocate buffer. \n" << endl;
                   cout <<   "Aborting......\n" << endl;
                        /* terminate program if out of memory   */
                   exit(1);
                }
        }
        return(dptr);
}

