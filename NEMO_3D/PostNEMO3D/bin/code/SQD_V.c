# include <math.h>
# include <stdio.h>
# include <stdlib.h>

//# include "nrutil.c"


FILE *fp2;
FILE *fp3;

main()
{
  int i,j,k,m;
  double xx,yy,zz;
  double Six[8],Siy[8],Siz[8];
  double asi,dump,pot;
  
  fp2=fopen("SQD_V", "w");
  if (fp2==NULL)
    {
      printf("Can't open input file\n");
      exit(-1);
    }

  fp3=fopen("r_atom", "w");
  if (fp3==NULL)
    {
      printf("Can't open input file\n");
      exit(-1);
    }
 
  asi = 5.430950000000e-01;

  Six[0] = asi;
  Siy[0] = asi;
  Siz[0] = asi;
  
  Six[1] = asi;
  Siy[1] = asi/2.0;
  Siz[1] = asi/2.0;
  
  Six[2] = asi/2.0;
  Siy[2] = asi/2.0;
  Siz[2] = asi;

  Six[3] = asi/2.0;
  Siy[3] = asi;
  Siz[3] = asi/2.0;
  
  Six[4] = 3.0*asi/4.0;
  Siy[4] = 3.0*asi/4.0;
  Siz[4] = 3.0*asi/4.0;
  
  Six[5] = 3.0*asi/4.0;
  Siy[5] = asi/4.0;
  Siz[5] = asi/4.0;
  
  Six[6] = asi/4.0;
  Siy[6] = asi/4.0;
  Siz[6] = 3.0*asi/4.0;

  Six[7] = asi/4.0;
  Siy[7] = 3.0*asi/4.0;
  Siz[7] = asi/4.0;


  // define the region of the potential (in Angstroms)

  double Lx=11.0;		//32.8;	 //37.2;
  double Ly=11.0;
  double Lz=11.0;		//32.8;		//22.0;
  double Nx=(int)(Lx/asi);
  double Ny=(int)(Ly/asi);
  double Nz=(int)(Lz/asi);		

  double x1=5.4;
  double z1=5.4;
  double y1=5.4;	
  double a=0.0005; //eV/nm^2
  double Ex=0.0;  //eV/nm
  double Ez=0.0; //eV/nm
  
  double f1=0.0;
  double f2=0.0;
  double f3=0.0;			

  for(i=0;i<Nx;i++)
    for(j=0;j<Ny;j++)
      for(k=0;k<Nz;k++)
	for(m=0;m<8;m++)
	  {xx = Six[m]+(double)i*asi;
	  yy = Siy[m]+(double)j*asi;
	  zz = Siz[m]+(double)k*asi;
	
	  f1=a*((xx-x1)*(xx-x1)+(yy-y1)*(yy-y1)+(zz-z1)*(zz-z1));
	  pot=f1; 
	
	 f3=Ez*zz+Ex*xx;
	 pot=-pot-f3;

	  fprintf(fp3,"%lf %lf %lf\n",xx, yy,zz);
	  fprintf(fp2,"%lf\n", pot);

	  }
 
  fclose(fp3); 
  fclose(fp2);
 
}

