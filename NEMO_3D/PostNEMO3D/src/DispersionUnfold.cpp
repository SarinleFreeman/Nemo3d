#include <stdio.h>
#include <stdlib.h>
#include<string.h>
// #include <math.h>
#include "r2tensor.h"
#include "rvector.h"
#include "ivector.h"
#include "cmatrix.h"
#include "cvector.h"

#define system_error printf

/*
#define	Es	0
#define	Ep	1
#define	Vss_s	2
#define	Vsp_s	3
#define	Vpp_s	4
#define	Vpp_p	5
#define	Sorb	0
#define	Xorb	1
#define	Yorb	2
#define	Zorb	3 */
// #define NORB	4 
#define SMALL	1.0e-6
#define SMALL1  5.0e-5
#define	A_SM_SC	1.0
#define TEENY	1.0e-14
#define TINY	1.0e-14

#define INTSIZE 10           /* maximum characters integer (used in itoa)*/



/* Modification of piksrt in Numerical Recipes in C by Press, et al.	*/
/* Sort array arr[1..n] in ascending numerical order by straight 	*/
/* insertion.  On output arr is replaced by the sorted array and 	*/
/* the array index[1..n] contains the indices of the corresponding	*/
/* UNSORTED array.  That is index[j] = (index in unsorted array of the	*/
/* j-th element in the sorted array).  This is used to match sorted	*/
/* eigenvalues with eigenvectors.  Although this is an N^2 algorithm,	*/
/* there are only 20 bands in the model, thus the arrays will be 20	*/
/* elements, and, as Press, et al note, for N < 50 or so this method	*/
/* which has lower overhead than quicksort or heapsort will be good.	*/
void piksrtix(int n,double arr[],int index[])
{
	int i,j,ix;
	double a;

	/* Initialize index array */
	for (i=1;i<=n;i++)
		index[i] = i;

	/* Pick out an element	*/
	for (j=2;j<=n;j++)
	{
		a = arr[j];
		ix = index[j];
		i = j - 1;

		/* Look for a place to insert it	*/
		while ((i > 0) && (arr[i] > a))
		{
			arr[i+1] = arr[i];
			index[i+1] = index[i];
			i--;
		}  /* while */

		arr[i+1] = a;
		index[i+1] = ix;
	}  /* for j */
}  /* piksrtix */



/* Normalize vectors stored as columns of Wr, Wi [1..n][1..n] */
void normalize_all(double **Wr, double **Wi, int n, int Nval)
{
	double sum_sq, norm;
	int i, col;

	/* Loop over vectors (columns) */
	for (col=1;col<=Nval;col++)
	{
		sum_sq = 0.0;
		for (i=1;i<=n;i++)
			sum_sq += (Wr[i][col]*Wr[i][col] + Wi[i][col]*Wi[i][col]);
		norm = sqrt(sum_sq);

		for (i=1;i<=n;i++)
		{
			Wr[i][col] /= norm;
			Wi[i][col] /= norm;
		}  /* for i */
	}  /* for col */
}  /* normalize_all */



/* polar form Wm, Wp[1..n][1..n], m=mag, p=phase.						*/
void rect_to_polar_mat(double **Wr, double **Wi, double **Wm, double **Wp, int n, int Nval)
{
	double norm, phase;
	int i, j;
	/* Loop through columns */
	for(j=1;j<=Nval;j++)
	{
		/* Loop through rows */
		for(i=1;i<=n;i++)
		{
			norm = sqrt(Wr[i][j]*Wr[i][j]+Wi[i][j]*Wi[i][j]);
			if (norm > TINY)
			{
				/* In quadrants I & II if Im part pos.  Use arcos */
				/* in QIII, QIV, Im part neg:  get 2PI-acos */
				if (Wi[i][j] > 0.0)
					phase = acos(Wr[i][j]/norm);
				else
					phase = 2.0*M_PI - acos(Wr[i][j]/norm);
			}  /* if (norm > TINY) */
			else
			{
				norm = 0.0;
				phase = 0.0;
			}  /* else */


			Wm[i][j] = norm;	
			Wp[i][j] = phase;
		}  /* for i */
	}  /* for j */
}  /* rect_to_polar_mat */



/* Initialize matrix of exponentials Ur,Ui[1..n][1..n], n=Nc[0]*Nc[1]*Nc[2] */
/* Determines minima and maxima for each direction internally.				*/
void init_exp_mat(double **Ur, double **Ui, int Nc[])
{
	double Bvec[3][3],Gv[3],GRdot,norm;
	int i,j,ix,iy,iz,jx,jy,jz,iix,iiy,iiz,Rv[3],Nmax[3],Nmin[3],row,col;

	/* Compute limits */
	for(i=0;i<=2;i++)
	{
		if (Nc[i] % 2)
		{
			Nmax[i] = (Nc[i]-1)/2;
			Nmin[i] = -Nmax[i];
		}
		else
		{
			Nmin[i] = -(Nc[i]-2)/2;
			Nmax[i] = Nc[i]/2;
		}
	}  /* for i */	

	/* Compute norm */
	norm = sqrt(((double)Nc[0])*((double)Nc[1])*((double)Nc[2]));

	/* Assign Supercell Recip Lattice Vectors; for SC they	*/
	/* are one-component each and orthogonal.				*/
	/* Each vector is stored as a COLUMN of Bvec.			*/
	for(i=0;i<=2;i++)
	{
		for(j=0;j<=2;j++)
			Bvec[j][i] = 0.0;

		/* Units are inverse small cell cube edge */
		Bvec[i][i] = 2.0*M_PI/((double)Nc[i]);
	}  /* for i */

	/* Compute matrix; Rows are same Rv, cols same Gv */
	for(jx=0;jx<Nc[0];jx++)
	{
		Rv[0] = jx;
		for(jy=0;jy<Nc[1];jy++)
		{
			Rv[1] = jy;
			for(jz=0;jz<Nc[2];jz++)
			{
				Rv[2] = jz;
				row = 1 + Rv[2] + Nc[2]*Rv[1] + Nc[2]*Nc[1]*Rv[0];

				/* Now compute G */
				for(ix=Nmin[0];ix<=Nmax[0];ix++)
				{
					iix = ix - Nmin[0];
					for(iy=Nmin[1];iy<=Nmax[1];iy++)
					{
						iiy = iy - Nmin[1];
						for(iz=Nmin[2];iz<=Nmax[2];iz++)
						{
							iiz = iz - Nmin[2];
							col =1 + iiz + Nc[2]*iiy + Nc[2]*Nc[1]*iix;

							/* For SC dot prod is trivial so do it manually */
							Gv[0] = ((double)ix)*Bvec[0][0];	/* + iy*Bvec[1][0] + iz*Bvec[2][0] */
							Gv[1] = ((double)iy)*Bvec[1][1];	/* + ix*Bvec[0][1] + iz*Bvec[2][1] */
							Gv[2] = ((double)iz)*Bvec[2][2];	/* + ix*Bvec[0][2] + iy*Bvec[1][2] */

							/* Dot product: R.G */
							GRdot = 0.0;
							for(i=0;i<=2;i++)
								GRdot += ((double)Rv[i])*Gv[i];

							/* Compute & store complex exp */
							Ur[row][col] = cos(GRdot)/norm;
							Ui[row][col] = sin(GRdot)/norm;
//printf("Ur[%d][%d] = %lf, Ui[%d][%d] = %lf\n",row,col,Ur[row][col],row,col,Ui[row][col]);
						}  /* for ix */
					} /* for iy */
				}  /* for iz */
			}  /* for jz */
		}  /* for jy */
	}  /* for jx */
}  /* init_exp_mat */



/* Take supercell orb coeffs stored as {s,x,y,z} for each atom for each	*/
/* state -- state is a col of Wr, Wi[1..n][1..n], n=NORB*ncell and		*/
/* project out with unitary matrix Ur, Ui[1..ncell][1..ncell].  Store	*/
/* result as Bcr, Bci[1..n][1..n], again with cols being states, rows	*/
/* representing coeff at a given supercell recip lattice vector G_n.	*/
/* Indexing on G_n is by (iz,iy,ix), iz=z-index, etc.  For each G_n, 	*/
/* order of orbs is still {s,x,y,z}.  Nc[0..2] gives # cells 0->x etc	*/
/* Normalize each state for each G_n; zero out states which are < TINY	*/
/* k[0..2] is supercell wavevector used to first dephase states 		*/
void coeffs_super2bulk_nm(double **Wr, double **Wi, double **Ur, double **Ui,
	double **Bcr, double **Bci, double k[], double **norm_mat, int Nc[], int NORB, int ncell, int NBASIS, int Nval)
{
	double norm, kRdot, exr, exi, tr, ti;
	int n, i, j, ii, jj, orb, state, Gix;
	int ix, iy, iz, Rv[3], row, row_rv, atom;
	
	/* First dephase the matrix W for this state all orbs */
	for(ix=0;ix<Nc[0];ix++)
	{
		Rv[0] = ix;
		for(iy=0;iy<Nc[1];iy++)
		{
			Rv[1] = iy;
			for(iz=0;iz<Nc[2];iz++)
			{
				Rv[2] = iz;
				row_rv = Rv[2] + Nc[2]*Rv[1] + Nc[2]*Nc[1]*Rv[0];

				/* Compute k.R */
				kRdot = 0.0;
				for(i=0;i<=2;i++)
					kRdot += ((double)Rv[i])*k[i];

				/* Compute exp(-i*k.R); remember cosine is EVEN! */
				exr = cos(kRdot);
				exi = sin(-kRdot);

				/* Dephase; replace entry in matrix by dephased entry */
				for(atom=0;atom<NBASIS;atom++)
				{ 		
					for(orb=0;orb<NORB;orb++)
					{
//                                      	row = 1 + NORB*row_rv + orb;
						row = 1 + NORB*NBASIS*row_rv + NORB*atom + orb;

						for(state=1;state<=Nval;state++)  
						{
							/* Cols of W are normalized; TINY elts are num errors--zero */
							if (fabs(Wr[row][state]) < TINY)
								Wr[row][state] = 0.0;
							if (fabs(Wi[row][state]) < TINY)
								Wi[row][state] = 0.0;

							tr = exr*Wr[row][state] - exi*Wi[row][state];
							ti = exr*Wi[row][state] + exi*Wr[row][state];
							Wr[row][state] = tr;
							Wi[row][state] = ti;
						}  /* for state */
					}  /* for orb */
				} /* for atom */
			}  /* for ix */
		}  /* for iy */
	}  /* for iz */
	
	/* Do the projection and normalization one state at a time */
	for(state=1;state<=Nval;state++)  
	{
		/* Proceed by orbital type on projections */
		for(atom=0;atom<NBASIS;atom++)
		{ 
			for(orb=0;orb<NORB;orb++)
			{
				/* Project out with herm conj of U */
				for(i=1;i<=ncell;i++)
				{
					/* ii is index in W, Bc */
//					ii = 1 + NORB*(i-1) + orb;
					ii = 1 + NORB*NBASIS*(i-1) + NORB*atom + orb; 

					Bcr[ii][state] = 0.0;
					Bci[ii][state] = 0.0;
	
					for(j=1;j<=ncell;j++)
					{
						/* jj is index in W */
//						jj = 1 + NORB*(j-1) + orb;
						jj = 1 + NORB*NBASIS*(j-1) + NORB*atom + orb; 

						Bcr[ii][state] += (Ur[j][i]*Wr[jj][state] + Ui[j][i]*Wi[jj][state]);
						Bci[ii][state] += (Ur[j][i]*Wi[jj][state] - Ui[j][i]*Wr[jj][state]);
					}  /* for j */
				}  /* for i */
			}  /* for orb */
		} /* for atom */ 

		/* Now normalize bulk coeffs for each G_n for this state */
		for(i=0;i<ncell;i++)
		{
			/* Compute index of G_n for s-orb */
//			Gix = i*NORB + 1;
			Gix = i*NORB*NBASIS + 1; /* index of G_n for s-orb of atom 1 */
			norm = 0.0;
			for(atom=0;atom<NBASIS;atom++) 
			{ 
				for(orb=0;orb<NORB;orb++)
				{
//					ii = Gix + orb;
					ii = 1 + NORB*NBASIS*i + NORB*atom + orb; 

					/* Zero out tiny coeffs which might be left over */
					if (fabs(Bcr[ii][state]) < TINY)
						Bcr[ii][state] = 0.0;
				
					if (fabs(Bci[ii][state]) < TINY)
						Bci[ii][state] = 0.0;

					norm += Bcr[ii][state]*Bcr[ii][state] + Bci[ii][state]*Bci[ii][state];
				}  /* for orb */
			} /* for atom */ 

			norm = sqrt(norm);
			norm_mat[i+1][state] = norm;

//			printf("\na[%d][%d] = %le ",i+1,state,norm_mat[i+1][state]);

			/* Normalize if > TINY */
			if (norm > TINY)
			{
				for(atom=0;atom<NBASIS;atom++)
				{ 
					for(orb=0;orb<NORB;orb++)
					{
//						ii = Gix + orb;
						ii = 1 + NORB*NBASIS*i + NORB*atom + orb; 
						Bcr[ii][state] /= norm;
						Bci[ii][state] /= norm;
					}  /* for orb */
				} /* for atom */ 
			}  /* if */
		} /* for i */
	}  /* for state */
}  /* coeffs_super2bulk_nm */



/* This function takes an integer as input and converts it to a string. */
/* itoa is used to generate the names of files which contain eigenvectors */
/* from NEMO3D */

/*--------------------------------------------------------------------+
| itoa() - manage the sign, compute the string equivalent, and call   |
| memcpy().                                                           |
+--------------------------------------------------------------------*/
char *itoa(int value)
{
int count,                   /* number of characters in string       */
    i,                       /* loop control variable                */
    sign;                    /* determine if the value is negative   */
char *ptr,                   /* temporary pointer, index into string */
     *string,                /* return value                         */
     *temp;                  /* temporary string array               */

count = 0;
if ((sign = value) < 0)      /* assign value to sign, if negative    */
   {                         /* keep track and invert value          */
   value = -value;
   count++;                  /* increment count                      */
   }

/* allocate INTSIZE plus 2 bytes (sign and NULL)                     */
temp = (char *) malloc(INTSIZE + 2);
if (temp == NULL)
   {
   return(NULL);
   }
memset(temp,'\0', INTSIZE + 2);

string = (char *) malloc(INTSIZE + 2);
if (string == NULL)
   {
   return(NULL);
   }
memset(string,'\0', INTSIZE + 2);
ptr = string;                /* set temporary ptr to string          */

/*--------------------------------------------------------------------+
| NOTE: This process reverses the order of an integer, ie:            |
|       value = -1234 equates to: char [4321-]                        |
|       Reorder the values using for {} loop below                    |
+--------------------------------------------------------------------*/
do {
   *temp++ = value % 10 + '0';   /* obtain modulus and or with '0'   */
   count++;                      /* increment count, track iterations*/
   }  while (( value /= 10) >0);

if (sign < 0)                /* add '-' when sign is negative        */
   *temp++ = '-';

*temp-- = '\0';              /* ensure null terminated and point     */
                             /* to last char in array                */

/*--------------------------------------------------------------------+
| reorder the resulting char *string:                                 |
| temp - points to the last char in the temporary array               |
| ptr  - points to the first element in the string array              |
+--------------------------------------------------------------------*/
for (i = 0; i < count; i++, temp--, ptr++)
   {
   memcpy(ptr,temp,sizeof(char));
   }

return(string);
} /* itoa */



int unfold_dispersion(void)
{
	double **Hr, **Hi, **Wr, **Wi, *val, a_sm, k[3], params[6];
	double **Ur, **Ui, **Bcr, **Bci, **norm_mat, bz_max;
	double **H4r, **H4i, **W4r, **W4i, *val4, em[6], **Hwork, **Hwork4;
	double **Bvec, kpG[3], **Wm, **Wp, **W4m, **W4p, Gv[3];
	double kx, ky, kz;
        int err, status, Nc[3], i, j, ii, jj, ncell, n, *ixv, *ixv4, n1;
	int orb, iz, iy, ix, iix, iiy, iiz, ic, Nmin[3], Nmax[3], Eig_Count, iter_k, bandmodel, NORB, NBASIS, Nval;
	int periodic, periodic_y, periodic_z, nbulk, Ncc[3]; 
	char iname[81],res, s[80000], iname1[81]; 
	FILE *fp, *fp1;

	cmatrix Hc, Cvec, Hc4, Cvec4;
	cvectr Cval, Cval4;
	
	/* For simple cubic sm cell size is conv unit cell cube edge */
	a_sm = A_SM_SC;

	//printf("Enter the lattice parameter of a unit cell>");
	//scanf("%lf",&a_sm);

	bz_max = M_PI/a_sm;	

	printf("Enter number of small cells in x, y, z> ");
        scanf("%d%d%d",&Nc[0],&Nc[1],&Nc[2]);

	printf("Enter the number of basis atoms (8 for 1 unit cell of Si)>");
	scanf("%d",&NBASIS);
	//NBASIS=8; /* 8 Si atoms per unit cell */
	printf("Enter the number of orbitals per atom>");
	scanf("%d",&NORB);
	//NORB=10; /* for Bands_10_sp3d5ss_nospin */

	printf("\nPeriodic in y? (y/n) > ");
        scanf("%s",&res);
        if(res == 'y' || res == 'Y')  periodic_y=1;
        else periodic_y=0;

	printf("Periodic in z? (y/n) > ");
	scanf("%s",&res);
	if(res == 'y' || res == 'Y')  periodic_z=1; 
	else periodic_z=0; 

        /* Adjust the unit cell size and basis size according to the periodicity of supercell */
        /* Ncc[i]'s are the number of unit cells in direction i, these are defined taking into */
        /* account the periodicity of the lattice */
	if(periodic_y==1 && periodic_z==1) /* Bulk */
	{ 
		Ncc[0] = Nc[0]; Ncc[1] = Nc[1]; Ncc[2] = Nc[2];
		NBASIS = NBASIS; 
	}

	if(periodic_y==1 && periodic_z==0) /* Quantum well confined in z direction */
	{ 
		Ncc[0] = Nc[0]; Ncc[1] = Nc[1]; Ncc[2] = 1;
		NBASIS = NBASIS*Nc[2];
	}

	if(periodic_y==0 && periodic_z==1) /* Quantum well confined in y direction */
	{ 
		Ncc[0] = Nc[0]; Ncc[1] = 1; Ncc[2] = Nc[2];
		NBASIS = NBASIS*Nc[1];
	}

	if(periodic_y==0 && periodic_z==0) /* Nanowire confined in y and z directions */
	{ 
		Ncc[0] = Nc[0]; Ncc[1] = 1; Ncc[2] = 1;
		NBASIS = NBASIS*Nc[1]*Nc[2];
	}


	/* Compute number of cells and matrix dimension */
	ncell = Ncc[0]*Ncc[1]*Ncc[2]; /* Number of cells in a supercell */
	n = NBASIS*NORB*ncell; /* Size of eigenvector */

	/* Get the k-vector of superlattice */
	printf("\nEnter k-vector of superlattice, kx, ky, kz.\n");
        printf("Units [Pi/N[i]a[i]] a=small cell conv unit cell cube> ");
        scanf("%lf %lf %lf",&k[0],&k[1],&k[2]);

        for(i=0;i<=2;i++)  k[i] *= M_PI/Nc[i]; 

	printf("Enter the total number of eigenvectors (or eigenvalues)>");
	scanf("%d",&Nval);

	/* Allocate matrices and vectors */
//	Hc = Cmatrix(n,n);
//	Cvec = Cmatrix(n,n);
//	Cval = Cvectr(n);
// 	Hc4 = Cmatrix(NORB*nbulk,NORB*nbulk);
// 	Cvec4 = Cmatrix(NORB*nbulk,NORB*nbulk);
// 	Cval4 = Cvectr(NORB*nbulk);
//	Hr = R2tensor(1,n,1,n);
//	Hi = R2tensor(1,n,1,n);

//	Wr = R2tensor(1,n,1,n);
        Wr = R2tensor(1,n,1,Nval); 
//	Wi = R2tensor(1,n,1,n);
        Wi = R2tensor(1,n,1,Nval);
//	Bcr = R2tensor(1,n,1,n);
        Bcr = R2tensor(1,n,1,Nval);
//	Bci = R2tensor(1,n,1,n);
        Bci = R2tensor(1,n,1,Nval); 
//	Wm = R2tensor(1,n,1,n);
        Wm = R2tensor(1,n,1,Nval); 
//      Wp = R2tensor(1,n,1,n);
        Wp = R2tensor(1,n,1,Nval); 
//      Hwork = R2tensor(1,n,1,n);
        Hwork = R2tensor(1,n,1,Nval); 
//	val = Rvectr(n+1);
	val = Rvectr(Nval+1);
//	norm_mat = R2tensor(1,Nc[0]*Nc[1]*Nc[2],1,n);
	norm_mat = R2tensor(1,ncell,1,Nval);
	Ur = R2tensor(1,ncell,1,ncell);
	Ui = R2tensor(1,ncell,1,ncell);
//	H4r = R2tensor(1,NORB*nbulk,1,NORB*nbulk);
//	H4i = R2tensor(1,NORB*nbulk,1,NORB*nbulk);
//	W4r = R2tensor(1,NORB*nbulk,1,NORB*nbulk);
//	W4i = R2tensor(1,NORB*nbulk,1,NORB*nbulk);
//	W4m = R2tensor(1,NORB*nbulk,1,NORB*nbulk);
//	W4p = R2tensor(1,NORB*nbulk,1,NORB*nbulk);
//	Hwork4 = R2tensor(1,NORB*nbulk,1,NORB*nbulk);
	Bvec = R2tensor(0,2,0,2);
//	val4 = Rvectr(NORB*nbulk+1);
//	ixv = Ivectr(n+1);
	ixv = Ivectr(Nval+1);
//	ixv4 = Ivectr(NORB*nbulk+1);


       /* Read the eigenvalues and eigenvectors from NEMO3D */

	printf("Enter input file name containing eigenvectors and eigenvalues (Si_wire, Si_well or Si_bulk)> ");
	scanf("%s",iname1); 

        /* Read eigenvalues */
        strcpy(iname,iname1);
        strcat(iname,".nd_Ek_ascii");
        printf("%s\n",iname);

        fp = fopen(iname,"rt");
        if (fp == NULL)
                system_error("No such input file.");

        for(ii=1;ii<=5;ii++) fgets(s,8000,fp); /* Skip first 5 lines */

        for(i=1;i<=Nval;i++)
        {
                fgets(s,8000,fp);
                sscanf(s,"%lf %lf %lf %lf",&kx,&ky,&kz,&val[i]); 
				/* Only eigenvalue is stored. kx, ky, kz are entered by the user separately */
        }
        fclose(fp);

	/* Read eigenvectors */
	for(j=1;j<=Nval;j++)
	{
		strcpy(iname,iname1);
		strcat(iname,".nd_evec_");
		strcat(iname,itoa(j));
		strcat(iname,"_ascii");
		printf("%s\n",iname);

	        fp = fopen(iname,"rt");
        	if (fp == NULL)
               		system_error("No such input file.");

		for(ii=1;ii<=7;ii++) fgets(s,8000,fp); /* Skip first 7 lines */

		i=1;
		while (fgets(s,8000,fp)!=NULL)
		{
			sscanf(s,"%lf %lf",&Wr[i][j], &Wi[i][j]);
			i++;
		}
	        fclose(fp);
	}


	/* Assign Supercell Recip Lattice Vectors; for SC they	*/
	/* are one-component each and orthogonal.				*/
	/* Each vector is stored as a COLUMN of Bvec.			*/
	for(i=0;i<=2;i++)
	{
		for(j=0;j<=2;j++)
			Bvec[j][i] = 0.0;

		/* Units are inverse small cell cube edge */
		Bvec[i][i] = 2.0*M_PI/((double)Ncc[i]);
	}  /* for i */

	/* Compute min and max integers for Supercell Recip vects */
	for(i=0;i<=2;i++)
	{
		if (Ncc[i] % 2)
		{
			Nmax[i] = (Ncc[i]-1)/2;
			Nmin[i] = -Nmax[i];
		}
		else
		{
			Nmin[i] = -(Ncc[i]-2)/2;
			Nmax[i] = Ncc[i]/2;
		}
	}  /* for i */	

//	piksrtix(n,val,ixv);
	piksrtix(Nval,val,ixv); 

	/* Normalize */
	normalize_all(Wr,Wi,n,Nval);

	/* Save polar forms of supercell states */
	rect_to_polar_mat(Wr,Wi,Wm,Wp,n,Nval);

	/* Compute Unitary Marrix */
	init_exp_mat(Ur,Ui,Ncc);

	/* Dephase supercell states, project out bulk states */
	coeffs_super2bulk_nm(Wr,Wi,Ur,Ui,Bcr,Bci,k,norm_mat,Ncc,NORB,ncell,NBASIS,Nval);

	/* We are done with Wr, Wi -- supercell states are in polar form	*/
	/* in Wp, Wm, so reuse for polar forms Bcr, Bci, (mag,ph)=(Wr,Wi)	*/
	rect_to_polar_mat(Bcr,Bci,Wr,Wi,n,Nval);

//	fp = fopen("Orb_coeffs.out","w");
//	if (fp == NULL)
//		system_error("Cannot open orbital file for output.");

	fp1 = fopen("E-k.out","w");
	if (fp1 == NULL)
                system_error("Cannot open E-k file for output.");

//	fprintf(fp,"Nc[0] = %3d   Nc[1] = %3d   Nc[2] = %3d\n",Nc[0],Nc[1],Nc[2]);
//	fprintf(fp,"k = (%17.10le, %17.10le, %17.10le)[Pi/a_sm]\n\n",k[0]/M_PI,k[1]/M_PI,k[2]/M_PI);
//	fprintf(fp,"  n0  n1  n2    Orbital Coefficient\t\t\tWeight Coefficient\n\n");

	for(j=1;j<=Nval;j++) 
	{
		jj = ixv[j];
		for(ix=0;ix<Ncc[0];ix++)
		{
		iix = ix + Nmin[0];
			for(iy=0;iy<Ncc[1];iy++)
			{
			iiy = iy + Nmin[1];

				for(iz=0;iz<Ncc[2];iz++)
				{
 					iiz = iz + Nmin[2];
					/* Row index of s-orb for this G_n */
//					ic = ix + Nc[0]*iy + Nc[0]*Nc[1]*iz;
					ic = iz + Ncc[2]*iy + Ncc[2]*Ncc[1]*ix;
					i = 1 + NORB*ic;

					if (norm_mat[ic+1][jj]>SMALL1)// && iiy==0 && iiz==0)

					{
						kx = k[0]/M_PI + (double(iix)/M_PI)*Bvec[0][0];
						if (kx<-1) kx+=2;
						if (kx>1) kx-=2;

						ky = k[1]/M_PI + (double(iiy)/M_PI)*Bvec[1][1];
						if (ky<-1) ky+=2;
						if (ky>1) ky-=2;

						kz = k[2]/M_PI + (double(iiz)/M_PI)*Bvec[2][2];
						if (kz<-1) kz+=2;
						if (kz>1) kz-=2;

						if (ky==0 && kz==0) /* to plot E-k in [1 0 0] direction */
						//if (ky==kx && kz==0) /* to plot E-k in [1 1 0] direction */
						//if (ky==kx && kz==kx) /* to plot E-k in [1 1 1] direction */
						{
							fprintf(fp1,"%17.10le %17.10le %17.10le ",kx,ky,kz); 
							fprintf(fp1,"  %f \n",val[j]);
						}
					}

//					fprintf(fp,"%3d %3d %3d     %17.10le  %17.10le rad  %17.10le\n",
//					iix,iiy,iiz,Wr[i][jj],Wi[i][jj],norm_mat[ic+1][jj]);

					/* Print the other orb coeffs without G_n for ease of reading */
					for(orb=1;orb<NORB;orb++)
					{
						ii = i + orb;
//						fprintf(fp,"\t\t%17.10le  %17.10le rad\n",Wr[ii][jj], Wi[ii][jj]);
					}  /* for orb */
				}  /* for iz */
			}  /* for iy */
		}  /* for ix */
//		fprintf(fp,"\n");
	}  /* for j */
//	fclose(fp);
fclose(fp1); 

	/* Free up the allocated storage */	
//	rm_cmatrix(&Hc);
//	rm_cmatrix(&Cvec);
//	rm_cvectr(&Cval);
//	rm_cmatrix(&Hc4);
//	rm_cmatrix(&Cvec4);
//	rm_cvectr(&Cval4);

//	rm_r2tensor(&Hr,1);
//	rm_r2tensor(&Hi,1);
	rm_r2tensor(&Wr,1);
	rm_r2tensor(&Wi,1);
	rm_r2tensor(&Bcr,1);
	rm_r2tensor(&Bci,1);
	rm_r2tensor(&Wm,1);
	rm_r2tensor(&Wp,1);
	rm_r2tensor(&Hwork,1);
	rm_r2tensor(&norm_mat,1);
	rm_r2tensor(&Ur,1);
	rm_r2tensor(&Ui,1);
//	rm_r2tensor(&H4r,1);
//	rm_r2tensor(&H4i,1);
//	rm_r2tensor(&W4r,1);
//	rm_r2tensor(&W4i,1);
//	rm_r2tensor(&W4m,1);
//	rm_r2tensor(&W4p,1);
//	rm_r2tensor(&Hwork4,1);
	rm_r2tensor(&Bvec,0);
	rm_rvectr(&val);
	rm_ivectr(&ixv);
//	rm_rvectr(&val4);
//	rm_ivectr(&ixv4);

return(0);
}  /* main */
