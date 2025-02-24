#include <cmath>
#include <iostream>
#include <fstream>
#include <ctime>
#include "Random.h"
#include "Orbital.h"
using namespace std;

double wavefunction(double n, double exponent, double r, double x, double y, double z, int l, int m)
{
//note that the normalization factor of spherical harmonics is not important 
//because we only use the ratio of two wavefunctions to determine whether to accept or reject new point.
   double wf;
   if(l==0) wf = pow(r, n)*exp(-1*exponent*r);
   else if (l==1&&m==0) wf = pow(r, n)* exp(-1*exponent*r)*z/r;
   else if (l==1&&m==1) wf = pow(r, n)* exp(-1*exponent*r)*x/r;
   else if (l==1&&m==2) wf = pow(r, n)* exp(-1*exponent*r)*y/r;
   else if (l==2&&m==0) wf = pow(r, n)* exp(-1*exponent*r)*(x*y)/(r*r);
   else if (l==2&&m==1) wf = pow(r, n)* exp(-1*exponent*r)*(y*z)/(r*r);
   else if (l==2&&m==2) wf = pow(r, n)* exp(-1*exponent*r)*(z*x)/(r*r);
   else if (l==2&&m==3) wf = pow(r, n)* exp(-1*exponent*r)*(x*x-y*y)/(r*r);
   else if (l==2&&m==4) wf = pow(r, n)* exp(-1*exponent*r)*(3*z*z-r*r)/(r*r);
   return wf;
} 

int main(int argc, char** argv)
{
  if(argc<2) { 
    cout<<"Usage: executable orbital_input_file_name"<<endl;
    exit(1);
  }

//read orbital information
  Orbital* orblist = new Orbital[2]; 
  ifstream fin(argv[1]);
  int l,m;
  double power,exponent;
  for(int o1=0;o1<2;o1++){
    fin>>l>>m>>power>>exponent;
    orblist[o1].initialize(l,m,power,exponent);
  }   

//set integration parameters
  const int bmax=100000; // number of iterations in a sampling loop
  int iteration=50; // number of iteration in an averaging loop
  double range=2.0;  // range of variation in unit of 1/exponent for each wave function

//set random number seed
  int t=time(0);
  cout<<"Random Seed chosen by current time = "<<t<<endl;
  Random rnd;
  rnd.Init(t);

//initialize monte carlo integration
  double factor1=orblist[0].get_exponent();
  double factor2=orblist[1].get_exponent();
  double n1=orblist[0].get_power();
  double n2=orblist[1].get_power();
  int l1=orblist[0].get_l();
  int l2=orblist[1].get_l();
  int m1=orblist[0].get_m();
  int m2=orblist[1].get_m();
  double average=0.0, average2=0.0;
  double x1=0.5,y1=0.5,z1=0.5,x2=1.0,y2=1.0, z2=1.0;
  double x1_tmp,y1_tmp,z1_tmp,x2_tmp,y2_tmp,z2_tmp;
  double wf1,wf2,wf3,wf4,wf1_tmp,wf2_tmp,wf3_tmp,wf4_tmp;
  double delta1=1.0/factor1*range;
  double delta2=1.0/factor2*range; 
  double r1,r2,r3,r4,r1_tmp,r2_tmp,r3_tmp,r4_tmp; 
  double density, density_tmp,value;
  double convertor=27.2;
  r1=sqrt(x1*x1+y1*y1+z1*z1);
  r2=sqrt(x2*x2+y2*y2+z2*z2);
  wf1=wavefunction(n1, factor1, r1, x1, y1, z1, l1, m1);
  wf2=wavefunction(n1, factor1, r2, x2, y2, z2, l1, m1);
  wf3=wavefunction(n2, factor2, r1, x1, y1, z1, l2, m2);
  wf4=wavefunction(n2, factor2, r2, x2, y2, z2, l2, m2);
  density=wf1*wf4*wf1*wf4; 
  value=wf1*wf2*wf3*wf4/density/sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));

//start monte carlo integration
  for(int i=0;i<iteration;i++){ // averaging loop
    double integral=0.0;
    int count=0; 
    int count2=0;
    for(int b=0;b<bmax;b++){ // sampling loop
      x1_tmp=x1+delta1*(rnd.Uniform()-0.5);
      y1_tmp=y1+delta1*(rnd.Uniform()-0.5);
      z1_tmp=z1+delta1*(rnd.Uniform()-0.5);
      x2_tmp=x2+delta2*(rnd.Uniform()-0.5);
      y2_tmp=y2+delta2*(rnd.Uniform()-0.5);
      z2_tmp=z2+delta2*(rnd.Uniform()-0.5);
      r1_tmp=sqrt(x1_tmp*x1_tmp+y1_tmp*y1_tmp+z1_tmp*z1_tmp);
      r2_tmp=sqrt(x2_tmp*x2_tmp+y2_tmp*y2_tmp+z2_tmp*z2_tmp);
      wf1_tmp=wavefunction(n1, factor1, r1_tmp, x1_tmp, y1_tmp, z1_tmp, l1, m1);
      wf2_tmp=wavefunction(n1, factor1, r2_tmp, x2_tmp, y2_tmp, z2_tmp, l1, m1);
      wf3_tmp=wavefunction(n2, factor2, r1_tmp, x1_tmp, y1_tmp, z1_tmp, l2, m2);
      wf4_tmp=wavefunction(n2, factor2, r2_tmp, x2_tmp, y2_tmp, z2_tmp, l2, m2);
      density_tmp=wf1_tmp*wf4_tmp*wf1_tmp*wf4_tmp;
      double prob=fabs(density_tmp/density);
      double pick=rnd.Uniform();
      if(pick<prob) { 
         x1=x1_tmp; y1=y1_tmp; z1=z1_tmp;
         x2=x2_tmp; y2=y2_tmp; z2=z2_tmp;
         density=density_tmp;
         double  distance=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+(z1-z2)*(z1-z2));
         value=wf1_tmp*wf2_tmp*wf3_tmp*wf4_tmp/density/distance;
         count++;
      }
      integral+=value;
    } 
    double rate=double(count)/bmax;
    integral/=bmax;
    cout<<"integral = "<<integral*convertor<<endl;
    cout<<"acceptance rate = "<< rate <<endl;
    average+=integral;
    average2+=integral*integral;
  }
  average/=iteration;
  average2/=iteration;
  double uncertainty=sqrt((average2-average*average)/(iteration-1));
  cout<<"==================================================="<<endl;
  cout<<"The average of integral in unit of eV = "<<average*convertor<<endl;
  cout<<"The uncertainty of integral = "<<uncertainty*convertor<<endl;
  cout<<"==================================================="<<endl;
  return 0;
}

