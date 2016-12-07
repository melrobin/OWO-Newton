#include <iostream>
#include <valarray>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include "matrix.h"
#ifdef _WIN32
  #include <mkl.h>
  #define F77_FUNC(func) (func)
#else
  #include "fortran.h"
  #define F77_FUNC(func) func##_
#endif

#include "mlp.h"
#include "matops.h" 
extern double h(double);

using namespace std;

void MLP::compute_derivs(bool doHessian,double& grad_energy)
{
//PURPOSE:  Batch calculation of the gradient for single hidden unit layer
//          MLP training.  This gradient is the same as used in BP.
//INPUT:    trnfname  training file name                   STRING
//          N         Number of inputs                     INTEGER
//          M         Number of outputs                    INTEGER
//          Nh        Number of hidden units               INTEGER
//          Woh       hidden unit to output weight array   DOUBLE    M x Nh
//          Woi       input to output array size           DOUBLE   M x (N+1)
//          W         input to hidden units layer array    DOUBLE   Nh x  N+1
//OUTPUT:
//          G Gradient of the error function               DOUBLE   Nh x  N+1
//                                                 
//RETURNS:  Energy in the gradient--the trace of G^T*G
  double ZERO=0.0,ONE=1.0;
  int NPlusOne;
  size_t i;
  int INTONE=1;//,k;
  valarray<double> x(N+1),t(M),y(M),net(Nh),e(M),temp(Nh),xa(N+Nh+1);
  valarray<double> delta_o(M),delta_p(Nh),ones(1.,Nh);
  valarray<double> f_prime_net(Nh),O(Nh),temp_row(Nh);
  valarray<double> g((N+1)*Nh);
  ifstream trnfile;
  int Niw=(N+1)*Nh,Now=M*(N+Nh+1);
  int Nu=N+Nh+1;
  valarray<double> dy_dw(Niw),temp_g_row(N+1);
  matrix temp_mat(Nh,N+1),Hoi(Now,Niw),tempHoi(Now,Niw),tempHoi1(Nu,Niw);
  matrix dy_dWoi(M,N+1),dy_dWoh(M,Nh),H_R(Niw,Niw);
  NPlusOne = N + 1;
  valarray<double> wo(Nu); 

  //Because we store the derivatives that we calculate with this method in the
  //class itself, we must be sure to zero them out for each iteration.
//  zero(H);
  zero(G);
  zero(G1);
  zero(Ht);
  //zero(Hoi);

  for (i=0;i<Gt.size();i++)
   Gt[i]=0.;

  trnfile.open(trnFile.c_str());

  while (!trnfile.eof())
  {
     for (i=0;i<N;i++)
       trnfile >> x[i];
 
     x[N]=1.;
     if (!trnfile.eof())
         for (i=0;i<M;i++)
            trnfile >> t[i];
     else
       break;
        
     F77_FUNC(dgemv)("T",&NPlusOne,&Nh,&ONE,&W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
     O=net.apply(h);
     y = process_pattern(x); 
     xa[slice(0,N,1)]=x[slice(0,N,1)];
     xa[slice(N,Nh,1)]=O;
     xa[N+Nh]=1.;
     
    //Calculate the output delta.  In BP, this will be used to calculate 
    //the delta for the hidden units.
     e=t-y;
     delta_o=2.*e;
    // F77_FUNC(dger)(&M,&Nu,&ONE,&e[0],&INTONE,&xa[0],&INTONE,&dE_dWo[0][0],&M);   
     //Calculate the derivative of the net function, assumimg a sigmoidal
     //activation
     f_prime_net=(ones-O)*O;
      //Calculate the delta for the hidden units using the BP algorithm
     F77_FUNC(dgemv)("N",&Nh,&M,&ONE,&Woh[0][0],&Nh,&delta_o[0],&INTONE,&ZERO,
                &delta_p[0],&INTONE);
     delta_p *= f_prime_net;
        //Building the incremental part of the gradient which depends on 
        //the hidden unit delta and the current pattern.  This is a rank 1 
        //outer product
     F77_FUNC(dger)(&NPlusOne,&Nh,&ONE,&x[0],&INTONE,&delta_p[0],&INTONE,&G[0][0],&NPlusOne);
     if (doHessian)
     {  

        for (i=0;i<M;i++)
        {
         //This implements the Hessian equations for the input weights
          zero(temp_mat);
          //zero(tempHoi1);
          temp_row = Woh.extract_row(i);
          temp_row *= f_prime_net;
          F77_FUNC(dger)(&NPlusOne,&Nh,&ONE,&x[0],&INTONE,&temp_row[0],&INTONE,&temp_mat[0][0],&NPlusOne);
         // F77_FUNC(dger)(&NPlusOne,&Nh,&e[i],&x[0],&INTONE,&temp_row[0],&INTONE,&G1[0][0],&NPlusOne);
          dy_dw=unroll(temp_mat);
      //    G1 += e[i]*temp_mat;
          //Form the Gauss-Newton approximation of the input weight Hessian
          F77_FUNC(dger)(&Niw,&Niw,&ONE,&dy_dw[0],&INTONE,&dy_dw[0],&INTONE,&H_R[0][0],&Niw);
          F77_FUNC(dger)(&Niw,&Nu,&ONE,&dy_dw[0],&INTONE,&xa[0],&INTONE,&Hoi[i*Nu][0],&Niw);
          //Hoi.set_submatrix(i*Nu,0,tempHoi1); //Build the temporary Hoi matrix one step at a time.
        }
      }//if
     // Hoi += 2.*tempHoi/Nv; //add the accumulated result to Hoi
   }//while
   calculate_corr();
   for (i=0;i<M;i++)
    {
      wo[slice(0,N+1,1)]=Woi.extract_row(i);
      wo[slice(N+1,Nh,1)]=Woh.extract_row(i);
      Gt[slice(Niw+i*Nu,Nu,1)]=-R*wo+C.extract_row(i);
      Ht.set_submatrix(Niw+i*Nu,Niw+i*Nu,2.*R);
    }
   G = G / Nv;
 //  G1 = 2.*G1/Nv;

   H_R = 2.*H_R / Nv;
   Hoi = 2.*Hoi/Nv/Nv;
   Ht.set_submatrix(0,0,H_R);
 //  Ht.set_submatrix(Niw,0,Hoi);
 //  Ht.set_submatrix(0,Niw,transpose(Hoi));
   Gt[slice(0,Niw,1)]=unroll(G);
  // Gt[slice(Niw,M*Nu,1)]=-unroll(dE_dWo/Nv);
   grad_energy = (g*g).sum(); 
   trnfile.close();
}

matrix MLP::compute_gradient(void)
{
//PURPOSE:  Batch calculation of the gradient for single hidden unit layer
//          MLP training.  This gradient is the same as used in BP.
//INPUT:    trnfname  training file name                   STRING
//          N         Number of inputs                     INTEGER
//          M         Number of outputs                    INTEGER
//          Nh        Number of hidden units               INTEGER
//          Woh       hidden unit to output weight array   DOUBLE    M x Nh
//          Woi       input to output array size           DOUBLE   M x (N+1)
//          W         input to hidden units layer array    DOUBLE   Nh x  N+1
//OUTPUT:
//          G Gradient of the error function               DOUBLE   Nh x  N+1
//                                                 
//RETURNS:  Energy in the gradient--the trace of G^T*G
  double ZERO=0.0,ONE=1.0;
  int NPlusOne;
  size_t i;
  int INTONE=1;//,k;
  valarray<double> x(N+1),t(M),y(M),net(Nh),e(M),temp(Nh),xa(N+Nh+1);
  valarray<double> delta_o(M),delta_p(Nh),ones(1.,Nh);
  valarray<double> f_prime_net(Nh),O(Nh),temp_row(Nh);
  valarray<double> g((N+1)*Nh);
  ifstream trnfile;
  int Niw=(N+1)*Nh,Now=M*(N+Nh+1);
  int Nu=N+Nh+1;
 // valarray<double> dy_dw(Niw),temp_g_row(N+1);
//  matrix temp_mat(Nh,N+1),Hoi(Now,Niw),tempHoi(Now,Niw),tempHoi1(Nu,Niw);
//  matrix dy_dWoi(M,N+1),dy_dWoh(M,Nh);
  NPlusOne = N + 1;
//  valarray<double> wo(Nu); 

  //Because we store the derivatives that we calculate with this method in the
  //class itself, we must be sure to zero them out for each iteration.
  zero(G);

 // for (i=0;i<Gt.size();i++)
 //  Gt[i]=0.;

  trnfile.open(trnFile.c_str());

  while (!trnfile.eof())
  {
     for (i=0;i<N;i++)
       trnfile >> x[i];
 
     x[N]=1.;
     if (!trnfile.eof())
         for (i=0;i<M;i++)
            trnfile >> t[i];
     else
       break;
        
     F77_FUNC(dgemv)("T",&NPlusOne,&Nh,&ONE,&W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
     O=net.apply(h);
     y = process_pattern(x); 
    // xa[slice(0,N,1)]=x[slice(0,N,1)];
    // xa[slice(N,Nh,1)]=O;
    // xa[N+Nh]=1.;
     
    //Calculate the output delta.  In BP, this will be used to calculate 
    //the delta for the hidden units.
     e=t-y;
     delta_o=2.*e;
    // F77_FUNC(dger)(&M,&Nu,&ONE,&e[0],&INTONE,&xa[0],&INTONE,&dE_dWo[0][0],&M);   
     //Calculate the derivative of the net function, assumimg a sigmoidal
     //activation
     f_prime_net=(ones-O)*O;
      //Calculate the delta for the hidden units using the BP algorithm
     F77_FUNC(dgemv)("N",&Nh,&M,&ONE,&Woh[0][0],&Nh,&delta_o[0],&INTONE,&ZERO,
                &delta_p[0],&INTONE);
     delta_p *= f_prime_net;
        //Building the incremental part of the gradient which depends on 
        //the hidden unit delta and the current pattern.  This is a rank 1 
        //outer product
     F77_FUNC(dger)(&NPlusOne,&Nh,&ONE,&x[0],&INTONE,&delta_p[0],&INTONE,&G[0][0],&NPlusOne);
    } //while
   G=G/Nv;
   return(G); 
     
}
void compute_derivs_OIG(MLP& X,const matrix& D)
{
    double ZERO=0.0,ONE=1.0;
   int i,NPlusOne;
   int INTONE=1;
   int M=X.M,Nh=X.Nh,N=X.N;
   valarray<double> x(N+1),t(M),y(M),net(Nh),e(M);
   valarray<double> ones(1.,Nh);
   valarray<double> f_prime_net(Nh),O(Nh),temp_row(N+1);
   matrix V(M,N+1),dyp_dr(M,N+1);
   ifstream trnfile;
   NPlusOne = N + 1;

   zero(X.H_ig);
   for (i=0;i<N+1;i++)
     X.g_ig[i]=0;

   trnfile.open(X.trnFile.c_str());

  while (!trnfile.eof())
  {
     for (i=0;i<X.N;i++)
       trnfile >> x[i];
 
     x[N]=1.;
     if (!trnfile.eof())
         for (i = 0; i < X.M; i++)
            trnfile >> t[i];
     else
       break;
        
     dgemv_("T",&NPlusOne,&Nh,&ONE,&X.W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
     O=net.apply(h);
     y = X.process_pattern(x); 
 
    //Calculate the output delta.  In BP, this will be used to calculate 
    //the delta for the hidden units.
     e=t-y;
 
        //Calculate the derivative of the net function, assumimg a sigmoidal
        //activation
     f_prime_net=(ones-O)*O;
 
     V = X.Woh.scale_row(f_prime_net)*D;
     dyp_dr=V.scale_row(x);
     X.g_ig += transpose(dyp_dr)*e;

     for (i = 0;i < X.M; i++)
     {
        temp_row = dyp_dr.extract_row(i);
        dger_(&NPlusOne,&NPlusOne,&ONE,&temp_row[0],&INTONE,&temp_row[0],&INTONE,&X.H_ig[0][0],&NPlusOne);
     }
   
   }
   X.g_ig = 2.*X.g_ig/X.Nv;
   X.H_ig = 2.*X.H_ig/X.Nv; 
   trnfile.close();
}
void compute_derivs_OIT(MLP& X,const matrix& D)
{
   double ZERO=0.0,ONE=1.0;
   int i,NPlusOne;
   int INTONE=1;
   int M=X.M,Nh=X.Nh,N=X.N;
   valarray<double> x(N+1),t(M),y(M),net(Nh),e(M);
   valarray<double> ones(1.,Nh);
   valarray<double> f_prime_net(Nh),O(Nh),temp_row(N+1);
   ifstream trnfile;
   NPlusOne = N + 1;
   int NPlusOne2=NPlusOne*NPlusOne; // (N+1)^2
   valarray<double> temp_dy_dr(NPlusOne2),temp_g_row(N+1);
   matrix dy_droit(N+1,N+1),K(Nh,N+1);
   zero(X.H_OIT);
   zero(X.G_oit);

   K=D;
   trnfile.open(X.trnFile.c_str());

  while (!trnfile.eof())
  {
     for (i=0;i<X.N;i++)
       trnfile >> x[i];
 
     x[N]=1.;
     if (!trnfile.eof())
         for (i = 0; i < X.M; i++)
            trnfile >> t[i];
     else
       break;
        
     F77_FUNC(dgemv)("T",&NPlusOne,&Nh,&ONE,&X.W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
     O=net.apply(h);
     y = X.process_pattern(x); 
     e=t-y;

     f_prime_net=(ones-O)*O;

     for (i = 0;i < X.M; i++)
     {         
        zero(dy_droit);  
        temp_row = X.Woh.extract_row(i);
        temp_row *= f_prime_net;
        for (int k=0;k<Nh;k++)
        {
            temp_g_row=K.extract_row(k);
            F77_FUNC(dger)(&NPlusOne,&NPlusOne,&temp_row[k],&x[0],&INTONE,&temp_g_row[0],&INTONE,&dy_droit[0][0],&NPlusOne);
        } //for k=0 to Nh
        X.G_oit+=e[i]*dy_droit;
        temp_dy_dr=unroll(dy_droit);
        F77_FUNC(dger)(&NPlusOne2,&NPlusOne2,&ONE,&temp_dy_dr[0],&INTONE,&temp_dy_dr[0],&INTONE,&X.H_OIT[0][0],&NPlusOne2);
      }//for i=0 to M
    }//while
   X.H_OIT=2.*X.H_OIT/X.Nv;
   X.G_oit = -2.*X.G_oit/X.Nv;
   trnfile.close();
}

void compute_derivs_MOLF(MLP& X,const matrix& D)
{
   double ZERO=0.0,ONE=1.0;
   int i,NPlusOne;
   int INTONE=1;
   int M=X.M,Nh=X.Nh,N=X.N;
   valarray<double> x(N+1),t(M),y(M),net(Nh),e(M);
   valarray<double> ones(1.,Nh);
   valarray<double> f_prime_net(Nh),O(Nh),temp_row(N+1);
   ifstream trnfile;
   NPlusOne = N + 1;
   //int Nh2=X.Nh*X.Nh; // (N+1)^2
   valarray<double> g_molf(Nh),v(Nh),dy_dz(Nh);
   matrix K(Nh,N+1),H_molf(Nh,Nh);
  // zero(X.H_molf);
  // zero(X.G_molf);

   K=D;
   trnfile.open(X.trnFile.c_str());

  while (!trnfile.eof())
  {
     for (i=0;i<X.N;i++)
       trnfile >> x[i];
 
     x[N]=1.;
     if (!trnfile.eof())
         for (i = 0; i < X.M; i++)
            trnfile >> t[i];
     else
       break;
        
     F77_FUNC(dgemv)("T",&NPlusOne,&Nh,&ONE,&X.W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
     O=net.apply(h);
     y = X.process_pattern(x); 
     e=t-y;

     f_prime_net=(ones-O)*O;

     for (i = 0;i < X.M; i++)
     {          
        temp_row = X.Woh.extract_row(i);
        temp_row *= f_prime_net;
        v=K*x;
        dy_dz = temp_row*v;
        g_molf += e[i]*dy_dz; //Compute the MOLF gradient
      //  for (int m=0;m<Nh;k++)
      //  {
      //      temp_g_row=K.extract_row(m);
      //      F77_FUNC(dger)(&NPlusOne,&NPlusOne,&temp_row[k],&x[0],&INTONE,&temp_g_row[0],&INTONE,&dy_droit[0][0],&NPlusOne);
       // } //for k=0 to Nh
       // X.G_oit+=e[i]*dy_droit;
       // temp_dy_dr=unroll(dy_dr);
        F77_FUNC(dger)(&Nh,&Nh,&ONE,&dy_dz[0],&INTONE,&dy_dz[0],&INTONE,&H_molf[0][0],&Nh);
      }//for i=0 to M
    }//while
    H_molf = 2.*H_molf/X.Nv;
    g_molf = 2.*g_molf/X.Nv;
   trnfile.close();
}
