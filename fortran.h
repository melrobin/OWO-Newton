extern "C"
{
extern void dgemv_(const char *,int *,int *,const double *,double *,int *,double *,const int *,const double *,double *, const int *);
extern double ddot_(int *,double *,const int *,double *,const int *);
extern void dgesv_(int *,int *,double *,int *, int *,double *,int *,int *);
extern void dgels_(char *,int *,int *,int *,double *,int *,double *,int *,double *,int *,int *);
extern void daxpy_(int *,double *,double *,const int *,double *,const int *);
extern void dscal_(int *,const double *,double *,const int *);
extern void dcopy_(int *,double *,const int *,double *,const int *);
extern double dnrm2_(int *,double *,const int *);
extern void dgemm_(const char *,const char *,int *,int *,int *,double *,
		  double *,int *,double *,int *,double *,double *,int *);
extern double dposv_(char *, int *,int *,double *,int *,double *,int *,int *);
extern void dsysv_(char *,int *,int *,double *,int *,int *,double *, int *,double *,int *,int * );
extern void dpotrf_(char *,int *,double *,int *,int *);
extern void dgeqrf_(int *,int *,double *,int *,double *,double *,int *,int *);
extern void dger_(int *,int *,double *,double *,int *,double *,int *,double *,int *);
extern void dgelsy_(int *,int *,int *,double *,int *,double *,int *,int *,double *,int *,double *,int *,int *);
extern void dgelsd_(int *,int *,int *,double *,int *,double *,int *,double *, double *,int *,double *,int *,int *,int *);
extern void dgetrf_(int *,int *,double *,int *,int *,int *);
extern void dgecon_(const char *,int *,double *,int *,double *,double *,double *,int *,int *);
extern double dlange_(const char *,int *,int *,double *,int *,double *);
extern void dgesvd_(const char *,const char *,int *,int *,double *,int *,double *,double *,int *,double *,int *,double *,int *,int *);
}
