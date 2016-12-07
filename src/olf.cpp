#include <iostream>
#include <valarray>
#include <fstream>
#include <vector>
#include <string>
#include <cassert>
#include <cmath>
#include "matrix.h"
#include "fortran.h"
#include "matops.h"
#include "mlp.h"

using namespace std;
extern double h(double);
valarray<double> compute_newton_direction(const matrix &,const valarray<double> &, size_t& );
/*valarray<double> MLP::linesearch(matrix& D,matrix& Doh, matrix & Doi)
{

   double ZERO=0.0,ONE=1.0,H;
   int i,NPlusOne, Nv=0;
   int INTONE=1,INTTHREE=3,ipiv[3],info;
   size_t nlin;
   valarray<double> x(N+1),t(M),y(M),net(Nh),e(M);
   valarray<double> ones(1.,Nh),g(3),d(.5,3),tempg(3);
   valarray<double> temp(Nh),dy_dz1(M),dy_dz2(M),dy_dz3(M),g1(M);
   valarray<double> f_prime_net(Nh),O(Nh),delta_o(M);
   matrix tempW(Nh,N+1),tempWoi(M,N+1),tempWoh(M,Nh);
   matrix H3(3,3);
   ifstream trnfile;
   NPlusOne = N + 1;
   trnfile.open(trnFile.c_str());
   d[1]=d[2]=0;
   d[0]=.1;
   cout << "Here" << endl;
   for (int j=0;j<7;j++)
   {
       trnfile.clear();
      trnfile.seekg(0,ios::beg);
     zero(H3);
     H=0;
     Nv=0;
     for (int kk=0;kk<3;kk++)
          g[kk]=0.;
           for (int kk=0;kk<M;kk++)
          g1[kk]=0.;
     tempW=W+d[0]*D;
     tempWoh=Woh+d[1]*Doh;
     tempWoi=Woi + d[2]*Doi;
   while (!trnfile.eof())
   {
      for (i=0;i<N;i++)
        trnfile >> x[i];

      x[N]=1.0;
		
      if (!trnfile.eof())
         for (i=0;i<M;i++)
	   trnfile >> t[i];
      else
         break;

      Nv++;
      dgemv_("T",&NPlusOne,&Nh,&ONE,&tempW[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);		
      O=net.apply(h);
      dgemv_("T",&Nh,&M,&ONE,&tempWoh[0][0],&Nh,&O[0],&INTONE,&ZERO,&y[0],&INTONE);

		//Add the bypass weight contributions
      dgemv_("T",&NPlusOne,&M,&ONE,&tempWoi[0][0],&NPlusOne,&x[0],&INTONE,&ONE,&y[0],&INTONE);

      e = t-y;
     f_prime_net=(1.-O)*O;

        dgemv_("T",&NPlusOne,&Nh,&ONE,&D[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&temp[0],&INTONE);
	temp *=f_prime_net;

	dy_dz1=tempWoh*temp;
        dy_dz2=Doh*O;
        dy_dz3=Doi*x;
 
	 g1 += e*dy_dz1;
   //      g[1] += (e*dy_dz2).sum();
   //      g[2] += (e*dy_dz3).sum();
        H+=(dy_dz1*dy_dz1).sum();
        // H3[0][0] += (dy_dz1*dy_dz1).sum();
   //      H3[1][1] =1;//+= (dy_dz2*dy_dz2).sum();
    //     H3[2][2] =1;//+= (dy_dz3*dy_dz3).sum();
     //   H3[0][1] += (dy_dz1*dy_dz2).sum();
     //    H3[0][2] += (dy_dz1*dy_dz3).sum();
     //    H3[1][2] += (dy_dz2*dy_dz3).sum();

	} //While loop
       // g[0]=g1.sum();
//	H3[2][1]=H3[1][2];
//        H3[2][0]=H3[0][2];
//        H3[1][0]=H3[0][1];
     //   chol(H3);
   // d+=compute_newton_direction(H3,g,nlin);
      d[0] += g1.sum()/H;
     //  dgesv_(&INTTHREE,&INTONE,&H3[0][0],&INTTHREE,ipiv,&g[0],&INTTHREE,&info);
     //  assert(info==0);
      // d += g;
       cout << d[0] << " " << d[1]<< " " << d[2]<< " " << H << " " << g1.sum()<< " " << j <<endl;
  } //for j
    trnfile.close();
   return(d);
} */
double MLP::linesearch(const matrix &D)
//PURPOSE:  Implements the Newton-Raphson method for a line search
//RETURNS:  THE LEARNING FACTOR
{
   double ZERO=0.0,ONE=1.0,z=1;
   const int NITERS=5;
   int i,NPlusOne, Nv=0,INTONE=1;
   size_t nlin;
   valarray<double> x(N+1),t(M),y(M),net(Nh),e(M);
   valarray<double> ones(1.,Nh);
   valarray<double> temp(Nh),dy_dz(M),g(M);
   valarray<double> f_prime_net(Nh),O(Nh),delta_o(M);;
   valarray<double> H(M);
   matrix tempW(Nh,N+1),tempA(Nh,N+1);
   ifstream trnfile;
   double F0,F1;
   NPlusOne = N + 1;
   trnfile.open(trnFile.c_str());
   tempA=W;
  // F0=J2(tempA);
  // F1=J2(tempA+D);
  // z=F0/(2*F1); 
  // cout << z << endl;
  for (int uu=0;uu<NITERS;uu++)
  {
    trnfile.clear();
      trnfile.seekg(0,ios::beg);
   for (int vv=0;vv<M;vv++)
    g[vv]=H[vv]=0;
    //g.fill(0.);
   tempW=W+z*D;
   H=0;
   while (!trnfile.eof())
   {
      for (i=0;i<N;i++)
        trnfile >> x[i];

      x[N]=1.0;
		
      if (!trnfile.eof())
         for (i=0;i<M;i++)
	   trnfile >> t[i];
      else
         break;

      Nv++;
     // y=process_pattern(x,net,O);
      dgemv_("T",&NPlusOne,&Nh,&ONE,&tempW[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);		
      O=net.apply(h);
f_prime_net=(1.-O)*O;
      dgemv_("T",&Nh,&M,&ONE,&Woh[0][0],&Nh,&O[0],&INTONE,&ZERO,&y[0],&INTONE);
      dgemv_("T",&NPlusOne,&Nh,&ONE,&W[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&net[0],&INTONE);
      O=net.apply(h);
		//Add the bypass weight contributions
      dgemv_("T",&NPlusOne,&M,&ONE,&Woi[0][0],&NPlusOne,&x[0],&INTONE,&ONE,&y[0],&INTONE);

      e = t-y;
     

        dgemv_("T",&NPlusOne,&Nh,&ONE,&D[0][0],&NPlusOne,&x[0],&INTONE,&ZERO,&temp[0],&INTONE);
	temp *=f_prime_net;

	dy_dz=Woh*temp; 
	 g += e*dy_dz;
         H += (dy_dz*dy_dz);
  //       H += (dy_dz*dy_dz).sum();

	} //While loop
   g/=Nv;
   H/=Nv;
 //  double a=2.*J2(tempW+D)+2.*J2(tempW)-4.*J2(tempW+.5*D);
 //  double b=J2(tempW+D)-J2(tempW)-a;
   //double c=-g.sum();
 //  z=-b/2./a;
   //z = g.sum()/a/2;
    z += g.sum()/H.sum();
 //   cout << z << " " << g.sum() << endl;
  } //end for loop Newton iterations
    trnfile.close();
   return(z);
}
double MLP::linesearch(const matrix& D,int& iters)
{
   iters=0;
   double t=1.,alpha = 0.1, beta=0.7,K=J2(W);
   
   while (J2(W+t*D) > K+t*alpha*(unroll(G)*unroll(W)).sum()) 
    {
       t *= beta;
       iters++;
    }
   return(t);
}
/*double MLP::linesearch(double a,double a1,const matrix& D,double eps,int& Nit)
{
   const double tau=(3-sqrt(5))/2;
   double b,b1;

   Nit=0;
   double F3,F4;
   while ((a1-a) > eps)
   {

       b=a+tau*(a1-a);
       b1=a1-tau*(a1-a);
      F3=J2(W+b*D);
      F4=J2(W+b1*D);

       if (F3 < F4)
 
           a1=b1;
       else 
           a=b;

       Nit++;
    }
   return(a);
} */
double MLP::linesearch(double a1,double b1,const matrix & D,double eps,int& iters)
{
   const double c=(sqrt(5)-1)/2;
   double x1,x2,fx1,fx2,a,b;
   iters=0;
   a=a1;
   b=b1;
   
   x1=c*a+(1-c)*b;
   fx1=J2(W+x1*D);
   x2=(1-c)*a+c*b;
   fx2=J2(W+x2*D);
   while (b-a > eps)
  // for(i=0;i<N;i++)
   {
       iters++;
       if (fx1<fx2)
       {
          b=x2;
          x2=x1;
          fx2=fx1;
          x1=c*a+(1-c)*b;
          fx1=J2(W+x1*D);
       }
       else
       {
          a=x1;
          x1=x2;
          fx1=fx2;
          x2=(1-c)*a+c*b;
          fx2=J2(W+x2*D);
       };
       printf("%d %f %f %f %f %f\n",iters,x1,x2,fx1,fx2,b-a);
    //   if (fabs(b-a)<eps)
    //   {
    //      printf("Done after %d steps\n",i);
    //      return(a);
    //   }//if
    }//while
 //  if (i<= N)
 //    printf("Golden search succeeded!  Sweet\n");
 //  else
 //    printf("Really our golden search failed!  Please investigate the output.\n");
   return(a);
}
double MLP::linesearch(matrix &p, matrix &poi,matrix &poh)
{
   double z,dez=0,dez2=0;
   int i,n;
   ifstream infile;
   infile.open(trnFile.c_str());
   valarray<double> x(N+1),t(M),net(Nh),O(Nh),Dyp(M),Fgx(Nh),y(M);
   while(!infile.eof())
  {
            
           for(n=0; n<N; n++)
	      infile >> x[n];
		
     			if(infile.eof())
				break;
		    for(i=0; i<M; i++)
			infile >> t[i];	
	
				x[N]=1.0;
                           net = W*x;
                
                           O=net.apply(h);
                                y=Woi*x+Woh*O;
                           Fgx = p*x;
                           Fgx *= (O*(1.-O));

                             Dyp = poh*O+Woh*Fgx + poi*x;
                        dez += ((t-y)*Dyp).sum();
                             dez2 += (Dyp*Dyp).sum();
			}	// end of data pass
                        
			z = dez/dez2;
   infile.close();
   return(z);
} 
double MLP::linesearch(double a,double a1,const matrix& D,const matrix & D1, const matrix &D2,double eps,int& Nit)
{
   const double tau=(3-sqrt(5))/2;//,c=(sqrt(5)-1)/2.;
   double b,b1;//,x1,x2;

   Nit=0;
   double F3,F4;
   while ((a1-a) > eps)
   {

       b=a+tau*(a1-a);
       b1=a1-tau*(a1-a);

       F3=J2(W+b*D,Woi+b*D1,Woh+b*D2);
      F4=J2(W+b1*D,Woi+b1*D1,Woh+b1*D2);
if (F3 < F4)
           a1=b1;
       else
       Nit++;
    }
   return(a);
}
double MLP::linesearch(double a,double a1,const matrix & D1, const matrix &D2,double eps,int& Nit)
{
   const double tau=(3-sqrt(5))/2;
   double b,b1;

   Nit=0;
   double F1,F2,F3,F4;
   while ((a1-a) > eps)
   {

       b=a+tau*(a1-a);
       b1=a1-tau*(a1-a);
  
       F3=J2(Woi+b*D1,Woh+b*D2);
      F4=J2(Woi+b1*D1,Woh+b1*D2);

       if (F3 < F4)
  
           a1=b1;
       else
       Nit++;
    }
   return(a);
}

