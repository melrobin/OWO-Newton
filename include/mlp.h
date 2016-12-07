#include <string>
#include <valarray>
#ifndef MYDEFS_H
#define MYDEFS_H
#include "mydefs.h"
#include "matrix.h"
#include "ann.h"
//#include "lmdefs.h"
#endif
using namespace std;
class MLP : public ANN
{
   int Nh, M,Nu;
   int N;
   size_t Nv ;
   matrix G,R,Ri,C,Woh,Woi,W,H_ig,Ht,H_OIT,H_ONT,G_oit,G1;
   valarray<double> net_mean,net_std,Xstd;
   valarray<double> Et,g_ig,Gt;
   string trnFile;
   public:
   void setNumInputs(int num);
   void setNumOutputs(int num);
   void set_Woh(const matrix &);
   void set_Woi(const matrix &);
   void set_W(const matrix &);
   int getNumInputs() const;
   int getNumOutputs(void) ;
   MLP(int,int,int,string);
   ~MLP(void);
   size_t get_num_patterns(){return(Nv);}
   vector<double> train_BP(int,size_t,size_t,size_t,string);
   vector<double> train_CG(int,size_t,size_t,size_t,string);
   vector<double> train_OWO_BP(int,size_t,size_t,size_t,string);
   vector<double> train_OWO_Newton(int ,size_t,size_t ,size_t,string);
   vector<double> train_OIG_HWO(int ,size_t,size_t ,size_t,string);
   vector<double> train_OIT_HWO(int ,size_t,size_t ,size_t,string);
   vector<double> train_OIG_BP(int ,size_t,size_t ,size_t,string);  
   vector<double> train_OIT_BP(int ,size_t,size_t ,size_t,string);
   vector<double> train_OWO_HWO(int,size_t,size_t,size_t,string);
   vector<double> train_LM(int,size_t,size_t,size_t,string);
   vector<double> train_OIT(int,size_t,size_t,size_t,string);
   vector<double> train_ONT(int,size_t,size_t,size_t,string);
   vector<double> train_MOLF(int,size_t,size_t,size_t,string);
   friend void compute_derivs_OIG(MLP&,const matrix& );
   friend void compute_derivs_OIT(MLP&,const matrix& );
   friend void compute_derivs_MOLF(MLP &,const matrix& );
   friend void compute_derivs_ONT(MLP &,const matrix& );
   friend double validate(MLP *, string);
   size_t calculate_stats();
   size_t init_mlp(string);
   void spit_info();
   void calculate_net_stats();
   void init_weights();
   void net_control(double,double);
   double owo();
   void write_weights(string);
   matrix get_output_weights(void);
   matrix get_input_weights(void);
   matrix get_Hessian(void);
   matrix get_gradient(void);
   void set_Nh(size_t n){Nh=n;}
   void calculate_corr();
   void update_weights(double,const matrix&);
   void save_weights();
   void restore_weights();
   void compute_derivs(bool,double &);
   matrix compute_gradient(void);
   double J(const string );
   double J(const matrix &);
   double J2(const matrix &);
   double J2(const matrix&,const matrix&,const matrix&);
   double J2(const matrix& ,const matrix& ) ;
   double linesearch(const matrix &);
   double linesearch(const matrix &,int&);
 //  valarray<double> linesearch(matrix& ,matrix&,matrix& );
   double linesearch(matrix &,matrix &,matrix &);
   double linesearch(double ,double,const matrix&,double,int&);
   double linesearch(double ,double ,const matrix&,const matrix &, const matrix &,double ,int& );
   double linesearch(double a,double a1,const matrix &, const matrix &,double ,int& );
   valarray<double> process_pattern(valarray<double>& );
   valarray<double> process_pattern(valarray<double>& ,valarray<double>&,valarray<double>&);
};

