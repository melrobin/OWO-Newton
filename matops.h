#include <valarray>
#include "matrix.h"
using namespace std;
extern matrix reshape(const valarray<double> &,size_t, size_t);
extern valarray<double> unroll(const matrix &);
extern matrix vertcat(bool, bool,const matrix &, const matrix &);
extern matrix transpose(const matrix &);
extern valarray<double> diag(const matrix &);
extern void zero(matrix &);
extern matrix shrink_matrix(size_t, size_t, const vector<int> );
extern matrix element_mult(const matrix&,const matrix&);
extern matrix qr(const matrix &);
extern matrix chol(const matrix &);
extern double rcond(const matrix &);
extern SVDRECORD svd(const matrix &);
