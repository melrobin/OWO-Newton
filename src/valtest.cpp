#include <iostream>
#include <valarray>
using namespace std;

int main()
{
  double val1[]={1,2,3,4},val2[]={5,6,7,8,9,10};
  size_t len[]={2,2},strides[]={5,1},len2[]={2,3},strides2[]={5,1};
  valarray<double> m1(val1,4),m2(val2,6),mat(10);
  gslice mySlice1(0,valarray<size_t>(len,2),valarray<size_t>(strides,2));
  gslice mySlice2(2,valarray<size_t>(len2,2),valarray<size_t>(strides2,2));
  mat[mySlice1]=10.*m1;
  mat[mySlice2]=2*2.*m2;
  for (int i =0;i<10;i++)
    cout << mat[i] << endl;
  return(0);
} 
