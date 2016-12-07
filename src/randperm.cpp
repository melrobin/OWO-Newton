#include <iostream>
#include <vector>
using namespace std;
extern vector<int> randperm(int);
extern "C" double genrand_real3(void);

int main()
{
   vector<int> values;
   int i,N=10;

   values=randperm(N);
   for (i=0;i<N;i++)
       cout << values[i] << endl;
   return(0);
}

