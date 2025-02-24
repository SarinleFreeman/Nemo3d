# include <cmath>
# include <iostream>
# include <fstream>

using namespace std;

FILE *fp1;
FILE *fp2;

int main(){

//Note: need to set this to the number of atoms

int N=876160;

//Just another way to define the dot
//double hbar=6.582e-16; //eV-s
//double m_eff=0.19;
//double m0=0.511e6/(9e16*1e18); // eV-s^2/nm^2
//double m_star=m0*m_eff;

//E0=1e-4 %eV
//lambda=hbar/sqrt(m_star*E0)
//lambda=40;
//double dot_decay_constant=40.0;
//double lambda=dot_decay_constant/sqrt(2);
//double E0=(hbar*hbar)/(lambda*lambda*m_star);
//double omega=E0/hbar;
//double curve=m_star*0.5*omega*omega;

//cout<<"E0="<<E0<<endl;
//cout<<"Curve="<<curve<<endl;

//Note: Reading the ascii file here
//It can be changed to binary file read also, but need to modify the fscanf format
//Note: In the ascii file, the header information needs to be deleted i.e. from the first line, there has to have x y z coordinates
fp1=fopen("siO2_si_dot.nd_rAtom_0_ascii", "r");
  if (fp1==NULL)
    {
      printf("Can't open input file 1\n");
      exit(-1);
    }

  fp2=fopen("SQD_V", "w");
  if (fp2==NULL)
    {
      printf("Can't open output file 2\n");
      exit(-1);
    }

  double a0=0.543095;

  double xx=0.0;
  double yy=0.0;
  double zz=0.0;
  double pot=0.0;
  double f3=0.0;

  double x1=20.0;
  double y1=20.0;

  //curve=1e-4; 
  double a=1e-4; //eV/nm^2

  //define these to be the field in Si
  double Ex=0.0;  //eV/nm
  double Ez=0.01; //eV/nm

  double sio2_length=7*a0;
  double dielectric_ratio=11.9/3.9;
  //field in SiO2
  double Ez_sio2=dielectric_ratio*Ez;

for(int i=0;i<N;i++){

        fscanf(fp1,"   %lf   %lf   %lf\n", &xx, &yy, &zz);
        if((i<4)||(i==N-1))
                printf("%12.10le %12.10le %12.10le\n", xx, yy, zz);

          pot=a*(xx-x1)*(xx-x1)+a*(yy-y1)*(yy-y1);
	 
	  if(zz<=sio2_length)	
	  f3=Ez_sio2*zz+Ex*xx;
	  else
	  f3=Ez*(zz+(dielectric_ratio-1)*sio2_length)+Ex*xx;	 

          pot=-pot-f3;
          //fprintf(fp3,"%lf %lf %lf\n",xx, yy,zz);

        fprintf(fp2,"%12.10le\n",pot);
        //if((i<4)||(i==N-1))
        //      printf("sum=%12.10le\n", pot);

  }

fclose(fp1);
fclose(fp2);

return 0;
}
