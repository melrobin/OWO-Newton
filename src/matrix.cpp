#include <iostream>
#include <valarray>
#include <cassert>
#include <algorithm>
#include <stdexcept>
#ifdef _WIN32
  #include <mkl.h>
  #define F77_FUNC(func) (func)
#else
  #include "fortran.h"
  #define F77_FUNC(func) func##_
#endif

//#ifdef _WIN32
//#include <mkl.h>
//#define F77_FUNC(func) (func)
//#else
//#include "fortran.h"
//#define F77_FUNC(func) (func_)
//#endif
#include "matrix.h" 
using namespace std;

matrix::matrix()
{
    // init(1,1);
      nrows=1;
   ncols=1;
   //assert((m>0)&&(n>0));
   try
   {
      data = new valarray<double>(1);
   }
   catch(exception& ba)
   {
      cout << "something caught: " << ba.what() << endl;
   }

}

matrix::matrix( size_t m, size_t n)
{
   // init(m,n);
      nrows=m;
   ncols=n;
   assert((m>0)&&(n>0));
   try
   {
      data = new valarray<double>(m*n);
   }
   catch(exception& ba)
   {
      cout << "something caught: " << ba.what() << endl;
   }

}
matrix::matrix(const matrix& X)
{
   nrows=X.nrows;
   ncols=X.ncols;
   data=new valarray<double>(nrows*ncols);
   *data = *(X.data);
}
matrix::matrix(size_t M,size_t N,double *A)
{
  this->nrows=M;
  this->ncols=N;
  data=new valarray<double>(M*N);
  for (size_t i = 0;i<M*N;i++)
    (*data)[i]=A[i];
}

matrix::matrix(const valarray<double> &X, size_t M,size_t N)
{
   if (M*N != X.size())
     {
        cout << "Incorrect number of data elements for dimensions chosen!" << endl;
        exit(0);
     }
   nrows=M;
   ncols=N;
   data = new valarray<double>(M*N);
   *data=X;
}

matrix matrix::operator-() const
{
  matrix temp(this->nrows,this->ncols);
  (*temp.data) = -*(this->data);
  return(temp);
}
matrix::~matrix()
{
  delete data;
}
void matrix::init(size_t m, size_t n)
{
   nrows=m;
   ncols=n;
   assert((m>0)&&(n>0));
   try
   {
      data = new valarray<double>(m*n);
   }
   catch(exception& ba)
   {
      cout << "something caught: " << ba.what() << endl;
   }
}

int matrix::getRow() const
{
  return nrows;
}

int matrix::getCol() const
{
   return ncols;
}

matrix& matrix::operator=(const matrix & X)
{
   if (&X != this)
    {
//     init(X.getRow(),X.getCol());
    nrows=X.getRow();
   ncols=X.getCol();
//   assert((m>0)&&(n>0));
//   try
 //  {
 //     data = new valarray<double>(m*n);
 //  }
 //  catch(exception& ba)
 //  {
 //     cout << "something caught: " << ba.what() << endl;
 //  }
   delete data;
   data = new valarray<double>(X.getRow()*X.getCol());
   *data=*(X.data);
    }
   return(*this);
}

matrix& matrix::make_from_array(double **X,int rows,int cols)
{
  *data=valarray<double>(&X[0][0],rows*cols);
   nrows = rows;
   ncols = cols;
   return(*this);
}

matrix matrix::operator*(const matrix& B)
{
  matrix temp(nrows,B.ncols);
  matrix b;
  double zero=0.0,one=1.0;

  b=B;
  if (ncols != b.nrows)
     throw 1;

  F77_FUNC(dgemm)("N","N",&b.ncols,&nrows,&b.nrows,&one,&b.data[0][0],&b.ncols,&data[0][0],&ncols,&zero,&temp.data[0][0],&b.ncols);
  return(temp);
}

matrix matrix::operator*(const double k) const
 
{
  matrix temp(nrows,ncols);
  *(temp.data) = k*(*data);;
  return temp;
}
void matrix::scale_column(int i,double k)
{
  int m=this->getRow(),n=this->getCol();
  slice_array<double> t = (*data)[slice(i,m,n)];
  valarray<double> k1(k,m);
  t *= k1; 
}

void matrix::scale_row(int i,double k)
{
  int n=this->getCol();
  slice_array<double> t = (*data)[slice(i*n,n,1)];
  valarray<double> k1(k,n);
  t *= k1;

}
void matrix::set_row(size_t i,const std::valarray<double> &x)
{
  int n=this->getCol();
  if (n != x.size())
  {
     cout << "set_row:  Check your dimensions" << endl;
     throw 1;
  }
  //slice_array<double> t = (*data)[slice(i*n,n,1)];
  //valarray<double> k1(k,n);
  //t *= k1;
  (*data)[slice(i*n,n,1)]=x;

}
void matrix::insert_row(size_t i,const std::valarray<double> &x)
{
  size_t m=this->getRow(),n=this->getCol();
  slice insert_slice(i*n,n,1);
  if (i>=m)
  {
    cout << "insert_row:  Improper row number" << endl;
    exit(0);
  }
  if (x.size() != n)
  {
    cout << "insert_row:  Improper number of elements to insert into matrix" << endl;
    exit(0);
  }
  (*data)[insert_slice]=x;
  
}
std::valarray<double> matrix::extract_column(size_t i)
{
  size_t m=getRow(),n=getCol();
  
  if (i >= n)
  {
     cout << "extract_column:  Improper column number" << endl;
     exit(0);
  }
  valarray<double> temp(m);
  temp = (*data)[slice(i,m,n)];
  return(temp);
}
std::valarray<double> matrix::extract_row(size_t i) const
{
   size_t m=this->getRow(),n=this->getCol();
   if (i >= m)
   {
      cout << "extract_row:  Improper row number" << endl;
      exit(0);
   }

   valarray<double> temp(n);
   temp = (*data)[slice(i*n,n,1)];
   return(temp);

}

std::valarray<double> matrix::extract_row(size_t i) 
{
   size_t m=this->getRow(),n=this->getCol();
   if (i >= m)
   {
      cout << "extract_row:  Improper row number" << endl;
      exit(0);
   }

   valarray<double> temp(n);
   temp = (*data)[slice(i*n,n,1)];
   return(temp);

}
void matrix::set_submatrix(size_t start_row, size_t start_column, const matrix &X)
//Sets the submatrix starting at (start_row, start_column) to the matrix X
{
  size_t M=this->getRow(),N=this->getCol(),M1=X.getRow(),N1=X.getCol();
  size_t start=start_column+start_row*N;
  size_t sizes[]={M1,N1},strides[]={N,1};
  gslice gs(start,valarray<size_t> (sizes,2),valarray<size_t> (strides,2));
  if ((start_row+M1 > M) || (start_column+N1 > N))
    {
       cout << "set_submatrix:  not enough space to implement" << endl;
       exit(0);
    }
  this->data[0][gs]=X.data[0];  
}
matrix matrix::scale_row(const valarray<double> &k) const
{
   size_t m=getRow(),n=getCol();
   int i;
   
   if (n != k.size())
   {
      cout << "scale_row:  Improper matrix size" << endl;
      exit(0);
   }
   matrix temp(m,n);
  
   valarray<double> temp1(n);
   for (i=0;i<m;i++)
    {
       temp1=(*data)[slice(i*n,n,1)];
       temp1 *= k;
       temp.data[0][slice(i*n,n,1)] = temp1;
    }
    return(temp);

}

void matrix::scale_row(int i,const valarray<double> &k)
{
   size_t m=getRow(),n=getCol();
   //int i;

   if (n != k.size())
   {
      cout << "scale_row:  Improper matrix size" << endl;
      exit(0);
   }
   matrix temp(m,n);

   valarray<double> temp1(n);
       temp1=(*data)[slice(i*n,n,1)];
       temp1 *= k;
       data[0][slice(i*n,n,1)] = temp1;

}

valarray<double> matrix::sum_rows() const
{
   size_t m=getRow(),n=getCol();
   int i;

   valarray<double> temp1(n);
   for (i=0;i<m;i++)
       temp1 += (*data)[slice(i*n,n,1)];
    return(temp1);

}

matrix matrix::scale_column(const valarray<double> &k) const
{
   size_t m=getRow(),n=getCol();
   int i;

   if (m != k.size())
   {
      cout << "scale_column:  Improper matrix size" << endl;
      exit(0);
   }
   matrix temp(m,n);
  
      valarray<double> temp1(m);
         for (i=0;i<n;i++)
          {
             temp1 = (*data)[slice(i,m,n)];
             temp1 *= k;
             temp.data[0][slice(i,m,n)] = temp1;
          }
         return(temp);
   
 }
   
ostream& operator<<(ostream& os, const matrix &m)
{
   int mval=m.getRow();
   int nval=m.getCol();
   for(int i=0; i < mval; i++) 
   {
      for(int j=0; j < nval; j++)
        os << m(i,j) << " ";
      flush(cout);
      os << endl;
   }
   return os;
}
matrix& matrix::operator+=(const matrix& X)
{
   *data += *(X.data);
  return (*this);
}
matrix& matrix::operator-=(const matrix& X)
{
   *data -= *(X.data);
  return (*this);
}
void matrix::clear()
{
   for (size_t i=0;i<(this->data)->size();i++)
	(*(this->data))[i]=0;
}
matrix& matrix:: operator/(const double x)
{
	*(this->data) /= x;
	return(*this);
}

matrix matrix::operator+(const matrix &b)
{

   if ((nrows != b.nrows) || (ncols != b.ncols))
      throw 1;
   matrix temp(nrows,ncols);
   temp = *this;
   *(temp.data) += *(b.data);
   return(temp);
}
matrix matrix::operator-(const matrix &b)
{
   if ((nrows != b.nrows) || (ncols != b.ncols))
      throw 1;
   matrix temp(nrows,ncols);
   temp = *this;
   *(temp.data) -= *(b.data);
   return(temp);

}
valarray<double> matrix::operator*(const valarray<double> &x)
{
  double zero=0.0,one=1.0; 
  int intone=1;
  assert(ncols == x.size());
  valarray<double> temp(nrows);
  valarray<double> b=x;

  dgemv_("T",&ncols,&nrows,&one,&((*data)[0]),&ncols,&b[0],&intone,&zero,&temp[0],&intone);
  return(temp);
}
matrix operator*(double k,const matrix& m)
{
  return m*k;
}
double matrix::operator()( size_t r,size_t c) const
{
   return row(r)[c];
}

