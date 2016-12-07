#include <string>
#include <iomanip>
#include "utils.h"
#include <valarray>
#include <ctime>
#include <iostream>
#include <fstream>
#include <vector>
#include "mlp.h"
#include <cassert>
#ifndef MYDEFS_H
#define MYDEFS_H
#include "mydefs.h"
#endif
#include "matops.h"
#include "fortran.h"
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/io.hpp>
extern double slete(double, double);
extern double ols(int,matrix &,matrix &,matrix &,int,const valarray<double> &,double *,matrix &);
using namespace std;
int ols_linear_solve(valarray<double> &,matrix &,valarray<double>,vector<int> & );
matrix compute_HWO_direction(const matrix&,const matrix&, size_t&);
valarray<double> compute_newton_direction(const matrix &,const valarray<double> &,size_t& );
extern double h(double);
std::vector<double> MLP::train_BP(int Nit,size_t N,size_t M,size_t Nh,string fname)
{
   vector<double> return_vals;
   ifstream trnfile;
extern double h(double);
  double ZERO=0.0,ONE=1.0,MSE=0.;
  int NPlusOne=N+1,INTONE=1;
	int i,j,k,u,v,it;
	double dmean=0.5, dstd=1.0,z3;
	int Nu,Nv;
	valarray<double> x(N+1);
valarray<double> f_prime_net(Nh),O(Nh);
  valarray<double> delta_o(M),delta_p(Nh),ones(1.,Nh);
	int nh=Nh,n;
       matrix H(3,3);
      matrix G(Nh,N+1),Goi(M,N+1),Goh(M,Nh);
	Nu=N+Nh+1;
      valarray<double> net(Nh);
       valarray<double> gg(3),t(M),tt(M),y(M),e(M);
     matrix R(Nu,Nu),C(M,Nu),woo(M,Nu);
       trnfile.open(fname.c_str());
 int mm=M;

    init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   MSE = owo();
  
  for(it=0;it<Nit;it++)
  { 
       Nv=0;
         return_vals.push_back(MSE);
       trnfile.clear();
            trnfile.seekg(ios::beg);
       while (!trnfile.eof())
       {
           for (i=0;i<N;i++)
              trnfile >> x[i];
 
           x[N]=1.;
           if (trnfile.eof())
             break;
           for (i=0;i<M;i++)
                 trnfile >> t[i];
           //else
             // break;
           Nv++;
           dgemv_("T",&NPlusOne,&nh,&ONE,&W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
           O=net.apply(h);
           y = process_pattern(x,net,O); 
     
    //Calculate the output delta.  In BP, this will be used to calculate 
    //the delta for the hidden units.
     e=t-y;
     delta_o=2.*e;
    // F77_FUNC(dger)(&M,&Nu,&ONE,&e[0],&INTONE,&xa[0],&INTONE,&dE_dWo[0][0],&M);   
     //Calculate the derivative of the net function, assumimg a sigmoidal
     //activation
     f_prime_net=(ones-O)*O;
      //Calculate the delta for the hidden units using the BP algorithm
     dgemv_("N",&nh,&mm,&ONE,&Woh[0][0],&nh,&delta_o[0],&INTONE,&ZERO,&delta_p[0],&INTONE);
     delta_p *= f_prime_net;
        //Building the incremental part of the gradient which depends on 
        //the hidden unit delta and the current pattern.  This is a rank 1 
        //outer product
     int nn=N;
     dger_(&NPlusOne,&nh,&ONE,&x[0],&INTONE,&delta_p[0],&INTONE,&G[0][0],&NPlusOne);
     dger_(&NPlusOne,&mm,&ONE,&x[0],&INTONE,&delta_o[0],&INTONE,&Goi[0][0],&NPlusOne);
     dger_(&nh,&mm,&ONE,&O[0],&INTONE,&delta_o[0],&INTONE,&Goh[0][0],&nh);
       }//end of while loop.  Makes one pass through the data to calculate initial gradients I guess.      
        G =G/ Nv;
        Goi=Goi/Nv;
        Goh =Goh/Nv;
          z3=linesearch(G,Goi,Goh);
          Woi += z3*Goi;
          Woh += z3*Goh;
          W += z3*G;
        MSE=J2(W);
        cout << it+1 << setw(15) << MSE << setw(12)<< z3 << endl;
     }
   trnfile.close();

	return (return_vals);
}
std::vector<double> MLP::train_ONT(int Nit,size_t N,size_t M,size_t Nh,string fname)
{
   double ZERO=0.0,ONE=1.0;
   int i,j,NPlusOne;
   size_t nlin,k=0;
   int INTONE=1;
   valarray<double> x(N+1),t(M),y(M),net(Nh),e(M);
   valarray<double> ones(1.,Nh);
   valarray<double> f_prime_net(Nh),O(Nh),temp_row(Nh);
   ifstream trnfile;
   NPlusOne = N + 1;
   int Nh2=Nh*Nh;
   valarray<double> g_ont(Nh2);
   matrix H_ont(Nh2,Nh2),R_ont(Nh,Nh),G(Nh,N+1);
   double MSE,Eg,alpha; 
   vector<double> return_vals;
   const double dmean = 0.5;
   const double dstd = 1.;
   //int pattern_count=0;
   valarray<double> d(Nh2),unrolled_dy_dr(Nh2),dn_dr(Nh);
  matrix G_hwo(Nh,N+1),D(Nh,N+1),tempD(Nh,N+1);
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   int n=Nh,iters;
   matrix dy_dr(Nh,Nh),temp_g_mat(Nh,Nh);
   valarray<double> g(Nh*(N+1));
   trnfile.open(fname.c_str());
   while(k < Nit)
   {
      trnfile.clear();
      trnfile.seekg(0,ios::beg); 
  
     MSE=owo();
     zero(temp_g_mat);
     zero(H_ont);
     return_vals.push_back(MSE);
     G=compute_gradient();
     G_hwo=compute_HWO_direction(Ri,G,nlin);
     D=G_hwo;
   //  for (int ii=0;ii < Nh;ii++)
     //  for (int jj=0;jj< N+1;jj++)
       //   D[ii][jj]=slete(-10,10);
     while (!trnfile.eof())
     {
        
        for (i=0;i<N;i++)
          trnfile >> x[i];

        x[N]=1.;
        if (!trnfile.eof())
          for (i = 0; i < M; i++)
            trnfile >> t[i];
        else
          break;
        y=process_pattern(x,net,O);//,tempD);
        e=t-y;
       
        f_prime_net=(ones-O)*O;

        dn_dr=D*x;
        for (i = 0;i < M; i++)
        {          
           zero(dy_dr);
           temp_row = Woh.extract_row(i);
           temp_row *= f_prime_net;
           dger_(&n,&n,&ONE,&dn_dr[0],&INTONE,&temp_row[0],&INTONE,&dy_dr[0][0],&n);
           temp_g_mat += e[i]*dy_dr; //Compute the ONT gradient
           unrolled_dy_dr=unroll(dy_dr);
           dger_(&Nh2,&Nh2,&ONE,&unrolled_dy_dr[0],&INTONE,&unrolled_dy_dr[0],&INTONE,&H_ont[0][0],&Nh2);
        }//for i=0 to M ;
      }//while 1

      g_ont=unroll(temp_g_mat);
      H_ont = 2.*H_ont/Nv;
      g_ont = 2.*g_ont/Nv;
      
      d=compute_newton_direction(H_ont,g_ont,nlin);
      R_ont=reshape(d,Nh,Nh);
     
      if (k==30)
       {
         ofstream outfile("linetest.txt");
         for (double zz=-.3;zz<1;zz+=.001)
           outfile << zz<< " " << J2(W+zz*D) << endl;
         outfile.close();
       }
 //alpha=linesearch(R_ont*D);
       // alpha /=2;
     // double alpha2;
      alpha=linesearch(-1,1,R_ont*D,1e-4,iters);
      //alpha2=linesearch(R_ont*D,iters);
    //  alpha=linesearch(R_ont*D,iters);
      update_weights(alpha,R_ont*D);
      
      cout <<setw(4) << k+1 << " ";

      cout << setw(15) << MSE;
      cout << setw(15) << J2(W);
      cout << setw(15) << alpha;
      cout << setw(15) << iters;
     // cout << setw(15) << (unroll(G)*(unroll(R_ont*D))).sum();
       cout << setw(15) << nlin << endl;
      k++;

     } //while k
      trnfile.close();          
    
     return(return_vals);
}

std::vector<double> MLP::train_MOLF(int Nit,size_t N,size_t M,size_t Nh,string fname)
{
   double ZERO=0.0,ONE=1.0;
   int i,NPlusOne;
   size_t nlin,k=0;
   int INTONE=1;
   valarray<double> x(N+1),t(M),y(M),net(Nh),e(M);
   valarray<double> ones(1.,Nh);
   valarray<double> f_prime_net(Nh),O(Nh),temp_row(N+1);
   ifstream trnfile;
   NPlusOne = N + 1;
   valarray<double> g_molf(Nh),v(Nh),dy_dz(Nh);
   matrix K(Nh,N+1),H_molf(Nh,Nh),R_molf(Nh,Nh);
   double MSE,Eg; 
   vector<double> return_vals;
   const double dmean = 0.5;
   const double dstd = 1.;
   matrix D(Nh,N+1);
   valarray<double> d(Nh);
  cout << "MOLF" << endl;
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   D=G;
   int n=Nh;
   trnfile.open(fname.c_str());
   while(k < Nit)
   {
     MSE=owo();
     return_vals.push_back(MSE);
     compute_derivs(false,Eg);
     while (!trnfile.eof())
     {
        for (i=0;i<N;i++)
          trnfile >> x[i];
 
        x[N]=1.;
        if (!trnfile.eof())
          for (i = 0; i < M; i++)
            trnfile >> t[i];
        else
          break;
        
        dgemv_("T",&NPlusOne,&n,&ONE,&W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
        O=net.apply(h);
        y = process_pattern(x); 
        e=t-y;

        f_prime_net=(ones-O)*O;

        for (i = 0;i < M; i++)
        {          
           temp_row = Woh.extract_row(i);
           temp_row *= f_prime_net;
           v=G*x;
           dy_dz = temp_row*v;
           g_molf += e[i]*dy_dz; //Compute the MOLF gradient
     
        dger_(&n,&n,&ONE,&dy_dz[0],&INTONE,&dy_dz[0],&INTONE,&H_molf[0][0],&n);
        }//for i=0 to M
      }//while 1
      H_molf = 2.*H_molf/Nv;
      g_molf = 2.*g_molf/Nv;
      
      d=compute_newton_direction(H_molf,g_molf,nlin);
      for (i=0;i<Nh;i++)
        R_molf[i][i]=d[i];
    
      update_weights(1,R_molf*G);
   
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE  << endl;
      k++;
     } //while 2
                
    trnfile.close(); 
     return(return_vals);
}

std::vector<double> MLP::train_CG(int Nit,size_t N,size_t M,size_t Nh,string fname)
{
   int INTONE=1;
   ifstream pinfile(fname.c_str());
   double ONE=1.;
   vector<double> return_vals;
   int n, k, i, it,iters;
   double Xd,B1,Xn;
   double MSE;
   double z;
   double gradient_energy = 0.0;
   double gradient_energy1 = 0.0;
   double gradient_energy2 = 0.0;
   const double dmean=0.5,dstd=1.;
   int jj;
   int NPlusOne=N+1;
   valarray<double> x(N+1),net(Nh);
   valarray<double> y(M),O(Nh),f(Nh),e(M);
   valarray<double> delta_po(M),delta_p(Nh),Fgx(Nh),Dyp(M);
   valarray<double> ty(M),Error(M),Et(M),cc(3);
   matrix g_wts(Nh,N+1),p(Nh,N+1);
   matrix goi_wts(M,N+1),goh_wts(M,Nh),poi(M,N+1),poh(M,Nh);
 
   init_weights();
   calculate_net_stats(); 
        net_control(dmean,dstd);
        MSE=owo();

	Xd=1;
	B1=0;
	for(it = 0; it<Nit; it++)
	{
            zero(g_wts);
            zero(goi_wts);
            zero(goh_wts);		
            MSE=0;
            pinfile.clear();
            pinfile.seekg(ios::beg);
	    while(!pinfile.eof())
	    {
		for(n=0; n<N; n++)
		  pinfile >> x[n];
		      
                if(pinfile.eof())
		  break;
		for(i=0; i<M; i++)
		  pinfile >> ty[i];
                  
                x[N] =1.0;
                net = W*x;
                O=net.apply(h);
                f =  O*(1.-O);

                y=Woi*x+Woh*O;
                e=ty-y;
                delta_po=2.*e;
                MSE += (e*e).sum();
                         
                delta_p=f*(transpose(Woh)*delta_po);

                n=Nh;
                jj=M;
		dger_(&NPlusOne,&n,&ONE,&x[0],&INTONE,&delta_p[0],&INTONE,&g_wts[0][0],&NPlusOne);
                dger_(&NPlusOne,&jj,&ONE,&x[0],&INTONE,&delta_po[0],&INTONE,&goi_wts[0][0],&NPlusOne);
                dger_(&n,&jj,&ONE,&O[0],&INTONE,&delta_po[0],&INTONE,&goh_wts[0][0],&n);	
	
	     } // while(!feof(binfile))

	     MSE/=Nv;
             return_vals.push_back(MSE);
   
	     g_wts = g_wts / Nv;	
             goi_wts=goi_wts/Nv;
             goh_wts=goh_wts/Nv;
//p=B1*p+g_wts;
  //           poi=B1*poi+goi_wts;
    //         poh=B1*poh+goh_wts;
             gradient_energy = 0.0;
	     gradient_energy1 = 0.0;
	     gradient_energy2 = 0.0;

	     gradient_energy=diag(transpose(g_wts)*g_wts).sum();
             gradient_energy1=diag(transpose(goi_wts)*goi_wts).sum();
             gradient_energy2=diag(transpose(goh_wts)*goh_wts).sum();
	
	     Xn = gradient_energy + gradient_energy1 + gradient_energy2;
	     B1=Xn/Xd;
	     Xd=Xn;

             p=B1*p+g_wts;
             poi=B1*poi+goi_wts;
             poh=B1*poh+goh_wts;	
             double z=linesearch(p,poi,poh);
        //     cc[0]=linesearch(W,p);
            //z=linesearch(0,2,p,poi,poh,1e-3,iters);
             cout << setw(7) << it << setw(15) << MSE << setw(13) << z  << endl;
 
             W+=z*p;
             Woi += z*poi;
             Woh += z*poh;
	}
        pinfile.close();
	return return_vals;

}

vector<double> MLP::train_OWO_HWO(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int k=0,iters;
   size_t nlin=0;
   double Eg;
   double alpha,MSE,d1,d2;
   const double dmean = 0.5;
   const double dstd = 1.;
   vector<double> return_vals;
   int N=inputs,Nh=hidden;
   matrix G_hwo(Nh,N+1);
 
   int Niw=(N+1)*Nh;
   valarray<double> g(Niw),d(Niw);
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   while(k < Nit)
    {
      MSE=owo();
      return_vals.push_back(MSE);
      compute_derivs(false,Eg);
      G_hwo=compute_HWO_direction(Ri,G,nlin);
      g=unroll(G_hwo);
      alpha=linesearch(G_hwo,iters);
      update_weights(alpha,G_hwo);
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << alpha ;
      cout << setw(13) << d1 ;
      cout << setw(13) << d2;
      cout << setw(13) << nlin << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << alpha ;
      report << setw(13) << d1;
      report << setw(13) << d2;
      report << setw(13) << nlin << endl;

      k++;
     } //while
                
      report.close();
      return(return_vals);
}

vector<double> MLP::train_OWO_BP(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int k=0,iters;
   size_t nlin=0;
  // double Eg;
   double alpha,MSE,d1,d2;
   const double dmean = 0.5;
   const double dstd = 1.;
   matrix G_bp(Nh,N+1);
   vector<double> result_vals;
   init_weights();

   calculate_net_stats();

   net_control(dmean,dstd);

   while(k < Nit)
   {
      MSE=owo();
      result_vals.push_back(MSE);
      G_bp=compute_gradient();
     // compute_derivs(false,Eg);
     // alpha=linesearch(G,iters);
      //alpha=linesearch(W,G_bp);
      alpha=linesearch(G_bp,iters);
      update_weights(alpha,G_bp);
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << alpha << endl ;
      //cout << setw(13) << iters << endl;
      //cout << setw(13) << d2;
      //cout << setw(13) << nlin << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << alpha ;
   //   report << setw(13) << d1;
   //   report << setw(13) << d2;
      report << setw(13) << nlin << endl;

      k++;
     } //while
                
      report.close();
      return(result_vals);
}
vector<double> MLP::train_LM(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
  vector<double> return_vals;
  ofstream report;
  int k=0,i;
  size_t nlin;
   double Eg;
   double MSE,alpha1;
 //  int Nu;
 //  double MSE_after_weight_change; 
   const double dmean = 0.5;
   const double dstd = 1.;
   int M=outputs,N=inputs,Nh=hidden;
 
   int Niw=(N+1)*Nh, Now=M*(N+Nh+1),Nw=Niw+Now,iters;
   valarray<double> g(Niw),d(Niw),dt(Nw),tempW(Niw);//,gt(Nw),tempW(Niw),gsave(Nw);
   valarray<double> tempWoi(M*(N+1)),tempWoh(M*Nh);
   matrix DW(Nh,N+1),DWoh(M,Nh),DWoi(M,N+1),Hsave(Nw,Nw);
   valarray<double> diags(Ht.getRow()),cc(3);
   init_weights();
   matrix zeroDW(Nh,N+1),zeroDWoi(M,N+1),zeroDWoh(M,Nh);
 /*  for (i=0;i<M;i++)
    {
      for (int j=0;j<N+1;j++)
         Woi[i][j]=slete(0.5,1);
       for (int j1=0;j1<Nh;j1++)
         Woh[i][j1]=slete(0.5,1);
    } */
calculate_net_stats();
net_control(dmean,dstd);
   MSE=owo();
   double alpha,alpha2;
   while(k < Nit)
   {
     MSE=J2(W);
       return_vals.push_back(MSE);
       compute_derivs(true,Eg);
     //  for (i=0;i<Ht.getRow();i++)
      //    Ht[i][i]+= lambda;
       dt=compute_newton_direction(Ht,Gt,nlin);
       tempW=dt[slice(0,Niw,1)];
       tempWoi=dt[slice(Niw,M*(N+1),1)];
       tempWoh=dt[slice(Niw+M*(N+1),M*Nh,1)];
       DW=reshape(tempW,Nh,N+1);
       DWoi=reshape(tempWoi,M,N+1);
       DWoh=reshape(tempWoh,M,Nh);
//double aa,bb;
//alpha=linesearch(aa,bb,DW);
//alpha=linesearch(0,1,DW,1e-3,iters);
// alpha1=linesearch(0,1,zeroDW,DWoi,zeroDWoh,1e-3,iters);
//alpha2=linesearch(0,1,zeroDW,zeroDWoi,DWoh,1e-3,iters);
cc=linesearch(DW,DWoh,DWoi);
       W += cc[0]*DW;
       Woi += cc[1]*DWoi;
       Woh += cc[2]*DWoh;
 
      cout << setw(4)  << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(13) << cc[0];
 cout << setw(12) << cc[1];
cout << setw(12) << cc[2];
      cout << setw(9) << alpha;
      cout << setw(11) << nlin << endl;

      k++;
     } //while
      report.close();
     return(return_vals);
  return(return_vals);
}
vector<double> MLP::train_OWO_Newton(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream weight_file;

   int i,k=0;
   size_t nlin;
   double alpha,MSE;
   const double dmean = 0.5;
   const double dstd = 1.;
   double ONE=1.;
   int INTONE=1;
   vector<double> return_vals;
   int N=inputs,Nh=hidden,M=outputs;
   ifstream trnfile;
   valarray<double> ones(1.,Nh);
   int NPlusOne= N+1;
   int Niw=(N+1)*Nh,iters;
   valarray<double> g(Niw),d(Niw),f_prime_net(Nh),temp_row(Nh),dy_dw(Niw),e(M),y(M);
   valarray<double> x(N+1),t(M),net(Nh),O(Nh);
   matrix D(Nh,N+1),temp_mat(Nh,N+1),G(Nh,N+1);
   matrix H_R(Niw,Niw),tempW(Nh,N+1);
   weight_file.open("initial_weights.txt");
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   trnfile.open(fname.c_str());

   while(k < Nit)
   {
       zero(H_R);
       zero(G);
       trnfile.clear();
       trnfile.seekg(ios::beg);
       MSE=owo();
       return_vals.push_back(MSE);
       while (!trnfile.eof())
       {
        for (i=0;i<N;i++)
          trnfile >> x[i];
 
        x[N]=1.;
        if (!trnfile.eof())
          for (i = 0; i < M; i++)
            trnfile >> t[i];
        else
          break;
 
        y=process_pattern(x,net,O);
        e=t-y;
       
        f_prime_net=(ones-O)*O;
        for (i=0;i<M;i++)
        {
         //This implements the Hessian equations for the input weights
          zero(temp_mat);
          temp_row = Woh.extract_row(i);
          temp_row *= f_prime_net;
          dger_(&NPlusOne,&Nh,&ONE,&x[0],&INTONE,&temp_row[0],&INTONE,&temp_mat[0][0],&NPlusOne);
          dger_(&NPlusOne,&Nh,&e[i],&x[0],&INTONE,&temp_row[0],&INTONE,&G[0][0],&NPlusOne);
          dy_dw=unroll(temp_mat);
          //Form the Gauss-Newton approximation of the input weight Hessian
          dger_(&Niw,&Niw,&ONE,&dy_dw[0],&INTONE,&dy_dw[0],&INTONE,&H_R[0][0],&Niw);
        } // for i
   }//while NOT EOF
       H_R = 2.*H_R/Nv;
       G = 2.*G/Nv;
       g=unroll(G);
       d=compute_newton_direction(H_R,g,nlin);
       D=reshape(d,Nh,N+1);
 
      // if (k==0)
      //     for (l=0;l<2;l+=.005)  
      //       lineout << l << " " <<J2(W+l*D) << endl;
       //alpha=linesearch(0,2,D,1e-3,iters);
        alpha=linesearch(D,iters);
      // tempW=W;
      //alpha=linesearch(D);
     //  tempW+= alpha*D;
     //  alpha+=linesearch(tempW,D);
//MSE=owo();  
//return_vals.push_back(MSE);
      update_weights(alpha,D);
      cout << setw(4)  << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(13) << J2(W);
   //   MSE=J2(W);
      cout << setw(12) << alpha ;
    //  cout << setw(13) << iters;
      //cout << setw(13) << (d*(H*d)).sum() ;
      cout << setw(13) << nlin << endl;
      k++;
     } //while k
      //report.close();
     //  lineout.close();
       trnfile.close();
     return(return_vals);
}
vector<double> MLP::train_OIT(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report,weight_file;
   int k=0;
   size_t nlin;
   double Eg;
   double MSE; 
   const double dmean = 0.5;
   const double dstd = 1.;
   vector<double> return_vals;
   int N=inputs,Nh=hidden;

   valarray<double> d((N+1)*(N+1)),g_oit((N+1)*(N+1));
   matrix D(Nh,N+1),R_oit(N+1,N+1);
   vector<int> dependent_rows;
   init_weights();

   calculate_net_stats();

   net_control(dmean,dstd);

   while(k < Nit)
   {
       MSE=owo();

       return_vals.push_back(MSE);
       compute_derivs(true,Eg);
  
       g_oit=unroll(G_oit);
       d=compute_newton_direction(H_OIT,-g_oit,nlin);
      // ols_linear_solve(d,H_OIT,-g_oit,dependent_rows);
       R_oit=reshape(d,N+1,N+1);

      update_weights(1,G*R_oit);
      cout << setw(4)  << k+1 << " ";
      cout << setw(15) << MSE;

      cout << setw(13) << rcond(H_OIT) << endl;

      k++;
     } //while
      report.close();
     return(return_vals);
}
vector<double> MLP::train_OIG_BP(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int i,k=0;
   size_t nlin=0;
   double Eg;
   double alpha=0,MSE;
   double alpha1=0; 
   const double dmean = 0.5;
   const double dstd = 1.;
   vector<double> return_vals;
   int N=inputs,Nh=hidden;
   matrix D(Nh,N+1);
 
   valarray<double> r(N+1);
   vector<int> dependent_rows;
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   while(k < Nit)
   {
     MSE=owo();
     return_vals.push_back(MSE);
      compute_derivs(true,Eg);
      D=G;
      compute_derivs_OIG(*this,D);
      ols_linear_solve(r,H_ig,g_ig,dependent_rows);
      for (i=0;i<Nh;i++)
          D.scale_row(i,r);
      W +=  D;
   
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << alpha ;
      cout << setw(13) << alpha1;
      cout << setw(13) << dependent_rows.size() << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << alpha ;

      report << setw(13) << alpha1;
      report << setw(13) << nlin << endl;

      k++;
     } //while
                
     report.close();
     return(return_vals);
}
vector<double> MLP::train_OIT_BP(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int k=0;
   size_t nlin=0;
   double Eg;
   double MSE; 
   const double dmean = 0.5;
   const double dstd = 1.;
   vector<double> return_vals;
   int N=inputs,Nh=hidden;
   matrix D(Nh,N+1),R_oit(N+1,N+1);
 
   valarray<double> d((N+1)*(N+1)),g_oit((N+1)*(N+1));
   vector<int> dependent_rows;
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   while(k < Nit)
   {
     MSE=owo();
     return_vals.push_back(MSE);
      compute_derivs(true,Eg);
      D=G;

     //G_hwo=compute_HWO_direction(Ri,G,nlin);
      compute_derivs_OIT(*this,D);
      g_oit=unroll(G_oit);
      d=compute_newton_direction(H_OIT,-g_oit,nlin);
       R_oit=reshape(d,N+1,N+1);
    
      update_weights(1,D*R_oit);
   
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE  << endl;
     // cout << setw(12) << alpha ;
     // cout << setw(13) << alpha1;
     // cout << setw(13) << dependent_rows.size() << endl;
     // report <<setw(4) << k+1 << " ";
    //  report << setw(15) << MSE;
     // report << setw(12) << alpha ;

     // report << setw(13) << alpha1;
     // report << setw(13) << nlin << endl;

      k++;
     } //while
                
     report.close();
     return(return_vals);
}
vector<double> MLP::train_OIG_HWO(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int i,k=0;
   size_t nlin=0;
   double Eg;
   double alpha=0,MSE;
   double alpha1=0; 
   const double dmean = 0.5;
   const double dstd = 1.;

   int N=inputs,Nh=hidden;//,num_dependent_rows;
   matrix G_hwo(Nh,N+1),D(Nh,N+1);
 
   valarray<double> r(N+1);
   vector<double>return_vals;
   vector<int> dependent_rows;
     init_weights();
     calculate_net_stats();
     net_control(dmean,dstd);
     calculate_net_stats();
    while(k < Nit)
    {
      MSE=owo();
      return_vals.push_back(MSE);
      compute_derivs(true,Eg);
      G_hwo=compute_HWO_direction(Ri,G,nlin);
      D=G_hwo;
      compute_derivs_OIG(*this,D);
      ols_linear_solve(r,H_ig,g_ig,dependent_rows);
      for (i=0;i<Nh;i++)
          D.scale_row(i,r);
        W +=  D;
   
      cout <<setw(4) << k+1 << " ";
      cout << setw(15) << MSE;
      cout << setw(12) << alpha ;
      cout << setw(13) << alpha1;
      cout << setw(13) << dependent_rows.size() << endl;
      report <<setw(4) << k+1 << " ";
      report << setw(15) << MSE;
      report << setw(12) << alpha ;

      report << setw(13) << alpha1;
      report << setw(13) << nlin << endl;

      k++;
     } //while
                
     report.close();
     return(return_vals);
}
vector<double> MLP::train_OIT_HWO(int Nit,size_t inputs,size_t outputs,size_t hidden,string fname)
{
   ofstream report;
   int k=0;
   size_t nlin=0;
   double Eg;
   double MSE;
   double alpha; 
   const double dmean = 0.5;
   const double dstd = 1.;
   vector<double> return_vals;
   int N=inputs,Nh=hidden;
   matrix D(Nh,N+1),R_oit(N+1,N+1),G_hwo(Nh,N+1);
 
   valarray<double> d((N+1)*(N+1)),g_oit((N+1)*(N+1));
   //vector<int> dependent_rows;
   init_weights();
   calculate_net_stats();
   net_control(dmean,dstd);
   while(k < Nit)
   {
     MSE=owo();
     return_vals.push_back(MSE);
      compute_derivs(true,Eg);
     
     G_hwo=compute_HWO_direction(Ri,G,nlin); 
     D=G_hwo;
      compute_derivs_OIT(*this,D);
      g_oit=unroll(G_oit);
      d=compute_newton_direction(H_OIT,-g_oit,nlin);
       R_oit=reshape(d,N+1,N+1);
      alpha=linesearch(D*R_oit);
      update_weights(alpha,D*R_oit);
   
      cout <<setw(4) << k+1;
      cout << setw(15) << MSE;
      cout << setw(12) << alpha << endl;
     // cout << setw(13) << alpha1;
     // cout << setw(13) << dependent_rows.size() << endl;
     // report <<setw(4) << k+1 << " ";
    //  report << setw(15) << MSE;
     // report << setw(12) << alpha ;

     // report << setw(13) << alpha1;
     // report << setw(13) << nlin << endl;

      k++;
     } //while
                
     report.close();
     return(return_vals);
}

