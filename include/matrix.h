#pragma once
#include <iostream>
#include <valarray>
#include <vector>
#include "slice_iter.h"

class matrix
{
/*! Private function to initialize the class.  I am unsure why it is still here.
Probably is safe to remove. */
   void init(size_t,size_t);
/*! Number of rows of the matrix */
   int nrows;
/*! Number of columns of the matrix */
   int ncols;
/*! The data in this class is represented as a pointer to a valarray.  There was a good
reason for doing it this way instead of a regular valarray, but I cannot remember. 
It is probably a good idea for it to remain this way because it would break a 
substantial amount of code...probably */
   valarray<double> *data;

   public:
/*! The constructor of the class.  Probably unneeded here */
      matrix();
/*! Initializes a matrix of a particular number of rows and columns */
      matrix(size_t ,size_t );
/*!  I believe that this is the copy constructor, but I am unsure */
      matrix(const matrix &);
/*! Constructs a matrix from a valarray.  Very useful.  Does this in the
MATLAB reshape sense, meaning columns first as we traverse the valarray */
      matrix(const valarray<double> &,size_t,size_t);
/*! Useful for constructing a matrix from an array.  The parameters should probably
be changed to be more consistent with the previous member function */
      matrix(size_t,size_t,double *);
/*! Class destructor.  For the most part unused */
      ~matrix();
      int getRow() const;
      int getCol() const;
      matrix(double **,int ,int);
      slice_iter<double> row(size_t);
      slice_iter<double> column(size_t);
      slice_iter<double> row(size_t) const;
      slice_iter<double> column(size_t) const;
      slice_iter<double> operator[](size_t i) const { return row(i); }
      slice_iter<double> operator[](size_t i) { return row(i); }  
      double& operator()(size_t,size_t);
      double operator() (size_t,size_t) const;
      friend std::ostream& operator<<(std::ostream &os, const matrix& m);
      friend matrix kron(const matrix &,const matrix & );
      friend std::valarray<double> unroll(const matrix &);
      friend matrix reshape(const std::valarray<double> &,size_t, size_t);
      friend matrix vertcat(bool, bool, const matrix &, const matrix &);
      friend matrix horzcat(const matrix &, const matrix &);
      friend valarray<double> diag(const matrix &);
      friend void zero(matrix &);
      friend matrix element_mult(const matrix &,const matrix &);
      friend double rcond(const matrix &);
      matrix& operator+=(const matrix& );
      matrix& operator-=(const matrix&);
      matrix operator+(const matrix &);
      matrix operator-(const matrix &);
      matrix operator*(const matrix &);
      matrix operator*(const double) const;
      matrix operator-() const;
      friend matrix operator*(double, const matrix&);
      matrix& operator=(const matrix &);
      void insert_row(size_t,const valarray<double> &);
      void set_row(size_t,const valarray<double> &);
      void set_submatrix(size_t, size_t, const matrix &);
      valarray<double> sum_rows(void) const;
      matrix& make_from_array(double **,int,int);
      void clear();
      matrix& operator/(const double);
      void scale_column(int,double);
      void scale_row(int,double);
      void scale_row(int ,const valarray<double> &);
      matrix scale_row(const valarray<double> &) const;
      matrix scale_column(const valarray<double> &) const;
      std::valarray<double> extract_column(size_t);
      std::valarray<double> extract_row(size_t) const;
      std::valarray<double> extract_row(size_t);
      std::valarray<double> operator* (const std::valarray<double> &);
      
};
inline slice_iter<double> matrix::column(size_t i) const
{
  return slice_iter<double>(data,slice(i,nrows,ncols));
}
inline slice_iter<double> matrix::row(size_t i) const
{
  return slice_iter<double>(data,slice(i*ncols,ncols,1));
}

inline slice_iter<double> matrix::column(size_t i)
{
  return slice_iter<double>(data,slice(i,nrows,ncols));
}
inline slice_iter<double> matrix::row(size_t i)
{
  return slice_iter<double>(data,slice(i*ncols,ncols,1));
}
inline double& matrix::operator()( size_t r,size_t c)
{
   return row(r)[c];
} 
/*! This is the structure that the svd() function returns */
struct SVDRECORD
{
/*! valarray of the singular values */
	valarray<double> S;
/*! We return the transpose of V in the singular value decomposition.
 This should probably be changed so that it returns a matrix V properly
 as is well known V is a square matrix */
	matrix V;
/*! We return matrix U */
	matrix U;
} ;




