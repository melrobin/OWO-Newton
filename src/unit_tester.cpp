#include <valarray>
#include <vector>
#include "matops.h"
#include "matrix.h"

using namespace std;
double g(double x)
{
  return(cos(x)-sin(x));
}
double golden_search(double a,double a1,double eps,int& Nit)
{
   const double tau=(3-sqrt(5))/2;
   double b,b1;
   
   Nit=0;
   double F1,F2,F3,F4;
   while ((a1-a) > eps)
   {
       
       b=a+tau*(a1-a);
       b1=a1-tau*(a1-a);
           F1=g(a),F2=g(a1),F3=g(b),F4=g(b1);

       if (F3 < F4)
         if (F1 < F3)
           a1=b;
         else
           a1=b1;
       else if (F3 > F4)
         if (F4>F2)
           a=b1;
         else
           a=b;
       else
       {
         a=b;
         a1=b1;
       }
       Nit++;
    }
   return(a);       
}      
   
int main()
{
  double arr[]={1.8,2.88,2.05,-.89,5.25,-2.95,-.95,-3.8,1.58,-2.69,-2.9,-1.04,-1.11,-.66,-.59,.8};
  double a[30] = {
            8.79,  6.11, -9.15,  9.57, -3.49,  9.84,
            9.93,  6.91, -7.93,  1.64,  4.02,  0.15,
            9.83,  5.04,  4.86,  8.83,  9.80, -8.99,
            5.45, -0.27,  4.85,  0.74, 10.00, -6.02,
            3.16,  7.98,  3.01,  5.80,  4.27, -5.31
        };
   double val1[]={1,2,3,4},val2[]={5,6,7,8,0,-1},val3[2]={.5,-.7};
//  double alpha=1.;
  valarray<double> mat(a,30);
  valarray<double> val(val1,4),vil(val2,6),vil3(val3,2);
//  vector<int> remove;
    matrix A,B,C,D(10,10);
    SVDRECORD t;
    A=matrix(mat,6,5);
    B=matrix(val,2,2),C=matrix(vil,2,3);
//  remove.push_back(0);
//  remove.push_back(2);
//  remove.push_back(4);
 // T=shrink_matrix(5,5,remove);
  //cout << A << endl << endl;
  //cout << T*A*transpose(T) << endl;
//  cout << B+alpha*C << endl;

//  cout << B.scale_row(vil3) << endl;
//  vil3=B.sum_rows();
//  cout << vil3[0] << " " << vil3[1] << endl;
 //int iters;
  //cout << golden_search(0,5,1e-5,iters) << endl;
  //cout << iters << " iterations of Golden Search" << endl;
 // cout << 1./rcond(A) << endl;
 // t=svd(A);
  D.set_submatrix(4,5,A);
  cout << A << endl;
  cout << D << endl;
  //cout << horzcat(B,C);
  //cout << horzcat(C,B);
  return(0);
}
 
