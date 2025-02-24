#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>

#include "r2tensor.h"
#include "rvector.h"
#include "ivector.h"

//#include "/home/boykin/lib/utility.c"

#define system_error printf

#define	FCC_P_CUBE	4		/* # atoms in NEMO-3D cube */
#define NEMO_HEAD	7		/* # lines to toss in NEMO3D ASCII header */
#define EPS_UNITARY	1.0e-10	/* effective 0 for unitary test */
#define EPS		1.0e-10
#define TR 2 /* All energies are doubly degenerate due to time reversal symmetry */


/* Convert positive integer n to string s */
/* Modified itoa from Kernighan & Ritchie, 2ed, p. 64 */
/* Builds in the reverse function from p. 62 */
void pos_itoa(int n, char s[])
{
	int i, j, nn;
	char c;
	
	/* Generate digits in reverse order */
	nn = n;
	i = 0;
	do
	{
		s[i++] = nn % 10 + '0';
	}
	while ((nn /= 10) > 0);
	
	s[i] = '\0';
	
	/* Now reverse the string in place */
	for(i=0, j=strlen(s)-1; i<j; i++, j--)
	{
		c = s[i];
		s[i] = s[j];
		s[j] = c;
	}  /* for */
}  /* pos_itoa */



/* Construct complete file name from template + ASCII(funm) + suffix	*/
/* return as fname.														*/
void full_filename(char template1[], int num, char suffix[], char fname[])
{
	char fnum[10];

	/* ASCII conversion */
	pos_itoa(num,fnum);

	/* Put together the name */
	strcpy(fname,template1);
	strcat(fname,fnum);
	strcat(fname,suffix);
}  /* full_filename */





/* Clean up input stream pointed to by fp */
void hoover(FILE *fp)
{
	char ch;

	do
		fscanf(fp,"%c",&ch);
	while (ch != '\n');
}  /* hoover */


/* Clean out input till = found */
void clean_to_eq(FILE *fp)
{
	char toss;

	do
		fscanf(fp,"%c",&toss);
	while (toss != '=');
}  /* clean_to_eq */




/* Get the header of the nemo3d ASCII file pointed to by fp	*/
/* return the number of orbitals, Norb and atoms Natom		*/
void get_header(FILE *fp, int *Norb, int *Natom)
{
	int j, n;

	/* Toss 2 lines, look for "=", get # orbs */
	for(j=1;j<=2;j++)
		hoover(fp);

	clean_to_eq(fp);
	fscanf(fp,"%d",&n);
	*Norb = n;
	hoover(fp);

	/* Look for "=", get # atoms */
	clean_to_eq(fp);
	fscanf(fp,"%d",&n);
	*Natom = n;
}  /* get_header */


/* Get nemo3d energy header, return max # states */
int maximum_states(FILE *fp)
{
	int i, maximillian;

	/* Discard first 2 lines */
	for(i=1;i<=2;i++)
		hoover(fp);

	/* Look for equals sign */
	clean_to_eq(fp);

	/* Get max # of states */
	fscanf(fp,"%d",&maximillian);
	
	return(maximillian);
}  /* maximum_states */

	




/* Get energies from nemo3d Ek ASCII file pointed to by fp	*/
/* Get only those energies GREATER THAN emin and LESS THAN	*/
/* emax. Use all_energies[0..max_states-1] return SUBVECTOR	*/ 
/* energies[0..*negy-1] setting a pointer to start energy	*/
/* in list of all also return negy = #energy states found &	*/
/* min_ix and max_ix, indices of minimum and maximum vects	*/
/* NOTE THAT VECTORS ARE UNIT-INDEXED.						*/
void get_nemo3d_energies(FILE *fp, double emin, double emax, double *all_energies,
				int max_states, int *negy, int *min_ix, int *max_ix)
{
	int i;
	double kx, ky, kz, egy;

	/* Discard first 5 lines */
	for(i=1;i<=5;i++)
		hoover(fp);

	/* Read in all states assuming ascending-sorted list */
	*min_ix = -1;
	*max_ix = -1;
	for(i=0;i<max_states;i++)
	{
		/* toss k since remap works only for k=0	*/
		fscanf(fp,"%lf%lf%lf%lf",&kx,&ky,&kx,&egy);
		all_energies[i] = egy;
	
		/* Get index of minimum vector:  UNIT indexed 	*/
		/* so vector (i+1) corresponds to 1st energy i 	*/
		/* Also attach convenience pointer energies to	*/
		/* this element of all_energies.				*/
		if ((egy > emin) && (*min_ix == -1))
			*min_ix = i + 1;

		/* Get index of maximum vector:  UNIT indexed 	*/
		/* so vector i corresponds to last energy (i-1)	*/
		if ((egy > emax) && (*max_ix == -1))
			*max_ix = i;
	}  /* for i */

	/* Now compute negy */
	*negy = (*max_ix - *min_ix + 1);

}  /* get_nemo3d_energies */
	
		
		
	




/* Read complex nemo3d vector from nemo3d ASCII file pointed to by fp	*/
/* Reorder so that the anion of a cell follows the cation of its cell	*/
/* Thus the new ordering is (by atom #):  (0,4); (1,5); (2,6); (3,7)	*/
/* Norb= #orbitals, Ncube= #NEMO3D cubes; nrec=2*FCC_P_CUBE*Ncube*Norb	*/
/* Return in vector vr, vi[0..nrec-1].  Discard header (1st 7 lines)	*/
void get_nemo3d_vector(FILE *fp, double vr[], double vi[], int Ncube, int Norb)
{
	int i, j, caoffst, cube_len, orb, ix, t;
	double re_pt, im_pt;

	/* Toss header -- we already looked at it to determine nrec */
	for(j=1;j<=NEMO_HEAD;j++)
		hoover(fp);

	/* Cation-anion offset and cube spacing */
	caoffst = Norb*FCC_P_CUBE;
	cube_len = 2*caoffst;

	/* i:  loop over cubes; t:  loop over atom type (c=0, a=1)	*/
	/* j:  loop over fcc within a cube.							*/
	for(i=0;i<Ncube;i++)
		for(t=0;t<2;t++)
			for(j=0;j<FCC_P_CUBE;j++)
				for(orb=0;orb<Norb;orb++)
				{
					/* Compute index of this orbital */
					/* for its atom in return vector */
					ix = orb + (2*j + t)*Norb + i*cube_len;
					fscanf(fp,"%lf%lf",&re_pt,&im_pt);
					vr[ix] = re_pt;
					vi[ix] = im_pt;
				}  /* for orb */
}  /* get_nemo3d_vector */




/* Index for calculation of time reversal state */
long get_index(int g,int s,int orb, int Norb)
{
        return (g*2+s)*(Norb/2)+orb+1;
} /* get_index */


/* Compute time reversal state vr_tr[1...nrec], vi_tr[1...nrec] corresponding to the 
   original nemo3d state vr[1...nrec], vi[1...nrec] */
void compute_tr_state(double *vr,double *vi,double *vr_tr,double *vi_tr,int Norb,int nrec)
{
        int o_max, orbital, g, ix0, ix1;

        o_max=Norb/2;
        for(int g=0;g<(nrec/Norb);g++)
           {
             for(int orbital=0;orbital<o_max;orbital++)
                {
		  ix0 = get_index(g,0,orbital,Norb);
	          ix1 = get_index(g,1,orbital,Norb);
                  vr_tr[ix1] = vr[ix0];
                  vi_tr[ix1] = -vr[ix0];
                  vr_tr[ix0] = -vr[ix1];
                  vi_tr[ix0] = vr[ix1];
                } /* for orbital */
           } /*for g */
} /* compute_tr_state */



/* Compute real space coordinates in terms of FCC prim dir vects av[0..2][0..2]	*/
/* where row is vector, col coords, 0->x etc.  Integer n[j] (of n[0..2]) gives	*/
/* # of units of av[j][0..2].  Store the components in the row-th row of		*/
/* Rv[0..N_FCC-1][0..2], again cols are components, 0->x etc.					*/
void real_space_coords(double **av, int n[], double **Rv, int row)
{
	int i, j;

	/* Loop over components */
	for(i=0;i<3;i++)
	{
		/* Start out with anion or cation as appropriate */
		Rv[row][i] = 0.0; 

		/* Loop over FCC direct lattice vectors */
		for(j=0;j<3;j++)
			Rv[row][i] += ((double)n[j])*av[j][i];
	}  /* for i */
}  /* real_space_coords */





void gen_Rvecs(double **R_a, int Nrect[])
{
	int i, ii, j, ip1, ip2, nemo[3], nn_a[3], n_a[3];
	int iix[3], nemo_ix, aix, err;
	double **av;

	/* Initialize FCC primitive direct lattice vectors 	*/
	/* row index is vector, col coord, 0->x etc			*/
	/* This routine is called ONCE only so ok to allocate & free locally */
	av = R2tensor(0,2,0,2);//allocate_real_matrix(0,2,0,2);	

	for(i=0;i<=2;i++)
	{
		for(j=0;j<=2;j++)
			av[i][j] = 0.5;

		av[i][i] = 0.0;
	}  /* for i */

	/* Initialize index for NEMO3D order */
	iix[0] = 0;
	iix[1] = 2;
	iix[2] = 1;

	/* Now loop through the NEMO cubes */
	for(nemo[0]=0;nemo[0]<Nrect[0];nemo[0]++)
		for(nemo[1]=0;nemo[1]<Nrect[1];nemo[1]++)
			for(nemo[2]=0;nemo[2]<Nrect[2];nemo[2]++)
			{
				/* Compute NEMO index */
				nemo_ix = ((nemo[0]*Nrect[1] + nemo[1])*Nrect[2] + nemo[2])*FCC_P_CUBE;
				aix = nemo_ix;
				
				/* Components in terms of av[i] for atom 0 of NEMO cube */
				for(i=0;i<3;i++)
				{
					ip1 = (i+1)%3;
					ip2 = (i+2)%3;
					nn_a[i] = -nemo[i] + nemo[ip1] + nemo[ip2];
				}  /* for i */

				/* Origin cell */
				real_space_coords(av,nn_a,R_a,aix);

				/* Map prim cells at a1-a3 = av[0]-av[2] */
				for (ii=0;ii<3;ii++)
				{
					i = iix[ii];
					aix = nemo_ix + ii + 1;

					ip1 = (i+1)%3;
					ip2 = (i+2)%3;
					
					n_a[i] = nn_a[i] - 1;
					n_a[ip1] = nn_a[ip1];
					n_a[ip2] = nn_a[ip2];

					real_space_coords(av,n_a,R_a,aix);
				}  /* for i */	
		}  /* for nemo[2] */

		/* Free our local storage */
		rm_r2tensor(&av,0);//free_real_matrix(av,0,2,0);
}  /* gen_Rvecs */




/* Generate small cell k vectors, ie supercell Gs.  Return as Gv[0..Nc-1][0..2]	*/
/* x->0 etc, Nc=4*Nrect[0]*Nrect[1]*Nrect[2], Nrect[0..2] = #cubes in x-z		*/
void gen_SmCell_k(double **Gv, int Nrect[])
{
	int Nx, Ny, Nz, tNx, tNy, tNz, count, NN[3], nx, ny, nz, n[3];
	int i, j, ii, ip1, ip2, nrec, total, Q_l[3], Q_u[3];
	double nv[4][3], d[3], dd[3], bv[3][3], twoPI, **mat, sum;
	double isqrt3, sqrt3PI;

	/* Constants used here */
	sum = sqrt(3.0);
	isqrt3 = 1.0/sum;
	sqrt3PI = M_PI*sum;
	
	/* Get # of 4-FCC cubes along each Cartesian Axis */
	Nx = Nrect[0];
	Ny = Nrect[1];
	Nz = Nrect[2];
	
	tNx = 2*Nx;
	tNy = 2*Ny;
	tNz = 2*Nz;
	
	NN[0] = Ny*Nz;
	NN[1] = Nx*Nz;
	NN[2] = Nx*Ny;
	
	/* # FCC prim cells in the supercell */
	nrec = 4*Nx*Ny*Nz;
	
	/* Set nv[0..2][] to be unit vects corresponding to recip lattice 			*/
	/* vects bv[0..2][]; set nv[3][] to be redundant unit vect (1,1,1)/sqrt(3)	*/
	for(i=0;i<3;i++)
	{
		for(j=0;j<3;j++)
			nv[i][j] = isqrt3;

		nv[i][i] = -isqrt3;

		/* For this one only we use i to index components not vector # */
		nv[3][i] = isqrt3;
	}  /* for i */

	/* Allocate storage */
	mat = R2tensor(0,nrec-1,0,2);//allocate_real_matrix(0,nrec-1,0,2);
	
	/* FCC prim recip lattice vectors */
	twoPI = 2.0*M_PI;

	/* Compute limits on Supercell part of K-vector:  has simple-cubic limits */
	for(i=0;i<3;i++)
		if (Nrect[i] % 2)
		{
			Q_l[i] = -(Nrect[i]-1)/2;
			Q_u[i] = (Nrect[i]-1)/2;
		}
		else
		{
			Q_l[i] = -(Nrect[i]-2)/2;
			Q_u[i] = Nrect[i]/2;
		}
	
	/* Generate FCC prim recip lattice vectors */
	for(i=0;i<3;i++)
	{
		ip1 = (i+1)%3;
		ip2 = (i+2)%3;
		bv[i][i] = -twoPI;
		bv[i][ip1] = twoPI;
		bv[i][ip2] = twoPI;
	}  /* for i */
	
	count = 0;
	for(nx=Q_l[0];nx<=Q_u[0];nx++)
		for(ny=Q_l[1];ny<=Q_u[1];ny++)
			for(nz=Q_l[2];nz<=Q_u[2];nz++)
			{
				n[0] = ny*Nz + Ny*nz;
				n[1] = nx*Nz + Nx*nz;
				n[2] = nx*Ny + Nx*ny;
				
				/* k for 0 additional wavevector:  mat is in the 	*/
				/* basis of FCC prim recip lattice vects.  That	is	*/
				/* mat[i][j] = coefficient of bv[j][0..2]			*/
				for(i=0;i<=2;i++)
				{
					dd[i] = 0.5*((double)n[i])/((double)NN[i]);
					mat[count][i] = dd[i];
				}  /* for i */
				count++;

				/* The 3 other wavevectors due to 4-FCC cube are 	*/
				/* (2*PI/a ){(1,0,0); (0,1,0); (0,0,1)}.  Note that */ 
				/* (2*Pi/a)(1,0,0) = 0.5*(bv[1]+bv[2]), etc.		*/
				for(i=0;i<3;i++)
				{
					ip1 = (i+1) % 3;
					ip2 = (i+2) % 3;

					for(j=0;j<3;j++)
						d[j] = dd[j];
					
					d[ip1] += 0.5;
					d[ip2] += 0.5;

					for(j=0;j<3;j++)
						if (d[j] > 0.5)
							d[j] -= 1.0;
						else if (d[j] <= -0.5)
							d[j] += 1.0;

					for(j=0;j<3;j++)
						mat[count][j] = d[j];
					
					count++;
				}  /* for i */
			}  /* for nz */
			
	total = count;
	
	/* Now loop through states; compute sm cell k	*/
	/* We convert from recip lattice vector basis	*/
	/* to Cartesian basis here.						*/
	for(i=0;i<total;i++)
	{
		for(j=0;j<=2;j++)
		{
			sum = 0.0;
			for(ii=0;ii<=2;ii++)
				sum += bv[ii][j]*mat[i][ii];
			Gv[i][j] = sum;
		}  /* for j */

		/* Shift Gv back into 1st zone if it falls out.				*/
		/* Geometry:  n*(r-r_0)=0 is eq of plane w/outward norm n	*/
		/* on which r_0 falls.  n*(r-r_0) < 0 if r is "inside"		*/
		/* First do the 8 [111] faces as +/- pairs.					*/
		for(ii=0;ii<4;ii++)
		{
			/* Dot product */
			sum = 0.0;
			for(j=0;j<3;j++)
				sum += Gv[i][j]*nv[ii][j];

			/* For these faces n*r_0=sqrt3PI */
			if (sum > sqrt3PI)
				for(j=0;j<3;j++)
					Gv[i][j] -= 2.0*sqrt3PI*nv[ii][j];
			else if (sum <= - sqrt3PI)
				for(j=0;j<3;j++)
					Gv[i][j] += 2.0*sqrt3PI*nv[ii][j];
		}  /* for ii */

		/* Now do the 6 [100] faces as +/- pairs */
		for(ii=0;ii<3;ii++)
		{
			/* Dot prod is trivial:  Gv[i][ii] */
			if (Gv[i][ii] > twoPI)
				Gv[i][ii] -= 2.0*twoPI;
			else if (Gv[i][ii] <= -twoPI)
				Gv[i][ii] += 2.0*twoPI;
		}  /* for ii */
	}  /* for i */
		
	/* Free storage */
	rm_r2tensor(&mat,0);//free_real_matrix(mat,0,nrec-1,0);
	
}  /* gen_SmCell_k */




/* Construct unitary projection matrix Ur, Ui[0..Nc-1][0..Nc-1] using G vectors	*/
/* Gv[0..Nc-1][0..2] and FCC prim cell vectors Rv[0..Nc-1][0..2] where row is	*/
/* vector and cols are components, 0->x etc.  									*/
void projection_matrix(double **Gv, double **Rv, double **Ur, double **Ui, int Nc)
{
	int r, c, j;
	double norm, RGdot;

	/* Overall normalization */
	norm = sqrt((double)Nc);

	/* Loop over rows, then cols of Ur, Ui */
	for(r=0;r<Nc;r++)
		for(c=0;c<Nc;c++)
		{
			/* Compute Rv[r].Gv[c]	*/
			RGdot = 0.0;
			for(j=0;j<3;j++)
				RGdot += Rv[r][j]*Gv[c][j];

			/* Compute exp(I*Rv[r].Gv[c])/sqrt(Nc), store in Ur, Ui	*/
			Ur[r][c] = cos(RGdot)/norm;
			Ui[r][c] = sin(RGdot)/norm;
		}  /* for c */
}  /* projection_matrix */





/* Compute hc(U)*v chosing every nth element of v, starting with element start	*/
/* Ur, Ui[0..nn-1][0..nn-1], vr, vi[0..nn-1]; return as uvr,uvi[0..nn-1].		*/
void U_dagger_v(double **Ur, double **Ui, double vr[], double vi[], 
				double uvr[], double uvi[], int start, int nth, int nn)
{
	int i, j, jj;
	double sumr, sumi;
	
	/* Loop over rows of hc(U), columns of cc(U)	*/
	for(i=0;i<nn;i++)
	{
		sumr = 0.0;
		sumi = 0.0;
		
		/* Summation loop */
		for(j=0;j<nn;j++)
		{
			jj = start + j*nth;
			
			sumr += (Ur[j][i]*vr[jj] + Ui[j][i]*vi[jj]);
			sumi += (Ur[j][i]*vi[jj] - Ui[j][i]*vr[jj]);
		}  /* for j */
		
		uvr[i] = sumr;
		uvi[i] = sumi;
	}  /* for i */
}  /* U_dagger_v */



/* Check unitarity of matrix Ur, Ui[0..n-1][0..n-1].  eps is effectively 0	*/
/* Return 1 for unitary, 0 for not unitary (to within eps, of course).		*/
int check_unitarity(double **Ur, double **Ui, int n, double eps)
{
	int i, j, m, pass;
	double sumr, sumi;

	/* Set pass test to 1 and compute matrix (U+)*U:  should be 1 */
	i = 0;
	pass = 1;
	while ((i<n) && (pass))
	{
		j = 0;
		while ((j<n) && (pass))
		{
			sumr = 0.0;
			sumi = 0.0;

			/* Compute (U+)*U[i][j] */
			for(m=0;m<n;m++)
			{
				sumr += (Ur[m][i]*Ur[m][j]+Ui[m][i]*Ui[m][j]);
				sumi += (Ur[m][i]*Ui[m][j]-Ui[m][i]*Ur[m][j]);
			}  /* for m */

			/* Different checks for diagonal and off-diag elements */
			if (j == i)
				pass = ((fabs(sumi) < eps) && (fabs(1.0-sumr) < eps));
			else
				pass = ((fabs(sumi) < eps) && (fabs(sumr) < eps));
			j++;
		}  /* while j */
//printf("\ni = %d",i);
		i++;
	}  /* while i */

	return(pass);
}  /* check_unitarity */



/* Project out probabilities for wavevectors (as ordered by gen_Gvec) from vector	*/
/* vr, vi[0..nrec], where nrec = 2*FCC_P_CUBE*Ncube*Norb.  Ncube = #nemo3d cubes,	*/
/* Norb = #orbitals, Nc = #remapped fcc prim cells = FCC_P_CUBE*Ncube.  Work-space	*/
/* is the vector wvr, wvi[0..Nc-1].  Store result (after projecting out and summing	*/
/* over all orbitals) as row eix of matrix Ekmat[0..Negy-1][0..Nc-1] where Negy =	*/
/* #energy eigenvectors we read in from the nemo3d output.  Cols of Ekmat are 		*/
/* prim cell k's (ie supercell G's) with order as determined by gen_Gvec.  The		*/
/* input vector MUST be the re-ordered vector read in by get_nemo3d_vector, which 	*/
/* is different from the ordering in the nemo3d output ASCII file.					*/
/* Use unitary projection matrix Ur, Ui[0..Nc-1][0..Nc-1] to do projection.			*/
void project_probs(double **Ekmat, double **Ur, double **Ui, double vr[], double vi[],
double wvr[], double wvi[], int Ncube, int Norb, int Nc, int Negy, int eix)
{
	int j, cube, Norbfcc, orb;

	/* Zero row eix of Ekmat since we deal with eix-th energy only */
	for(j=0;j<Nc;j++)
		Ekmat[eix][j] = 0.0;

	/* Orbitals per FCC is the skip factor */
	Norbfcc = 2*Norb;
	
	/* Loop over orbitals PER FCC.  This is 2*orbs/atom!	*/
	/* U_dagger_v does the multiply over fccs.  Note that	*/
	/* we resolve down to fcc (prim cell) vectors, NOT atom	*/
	/* vectors.  Thus we use total # orbs per fcc not atom	*/
	for(orb=0;orb<Norbfcc;orb++)
	{
		/* Start in 0-th FCC, orb-th orbital	*/
		U_dagger_v(Ur,Ui,vr,vi,wvr,wvi,orb,Norbfcc,Nc); 

		/* Add in probabilities for this orbital	*/
		/* Note: Orbital is cation-s, anion-pz, etc	*/
		for(j=0;j<Nc;j++)
			Ekmat[eix][j] += (wvr[j]*wvr[j]+wvi[j]*wvi[j]);
	}  /* for orb */
}  /* project_probs */

				
	

int fccproj(void)
{
	int min_ix, max_ix, negy, max_states, err, i, j, Ncube, Norb, Natom;
	int Nc, num, maxbands, eix, kix, nrec, is_unitary, n_nemo[3], tr, eix_tr;
	double emin, emax, *energies, *all_energies, **Gv, **Rv, **Ur, **Ui;
	double **Ekmat, tmp_vec[3], *vr, *vi, *wvr, *wvi, *vr_tr, *vi_tr, *sel_energies;
	char ename[81], iname[81], suffix[20], template1[81];
	FILE *fp;

	/* Get Ek input file */
	printf("Enter name of Ek ASCII file> ");
	scanf("%s",ename);
	
	/* Get min and max energies */
	printf("Enter min and max energies> ");
	scanf("%lf%lf",&emin,&emax);

	/* Attempt open */
	fp = fopen(ename,"r");
	if (fp == NULL)
		system_error("No such energy input file.");

	/* Get energies from nemo3d ASCII file */
	max_states = maximum_states(fp);

	/* Allocate storage */
	all_energies = Rvectr(max_states/*-1*/);//allocate_real_vector(0,max_states-1);

	/* Rewind file then get input */
	rewind(fp);
	get_nemo3d_energies(fp,emin,emax,all_energies,max_states,&negy,
					&min_ix,&max_ix);

	fclose(fp);

//	printf("min = %3d  max = %3d\n",min_ix,max_ix); // removed to be able to submit in queue

	/* Assign energies */
	energies = &(all_energies[min_ix-1]);

	/* Test print */
	for(i=0;i<negy;i++)
	{
		j = i + min_ix;
//		printf("%3d   %19.12le\n",j,energies[i]);  // removed to be able to submit in queue
	}  /* for i */

	/* Get nemo3d header from 1st vector file */
	printf("Enter file name template> ");
	scanf("%s",template1);
	
	/* Set suffix for all */
	strcpy(suffix,"_ascii");
	
	/* Set full file name to 1st file */
	full_filename(template1,min_ix,suffix,iname);

	printf("%s\n",iname);

	/* Attempt open & header read */
	fp = fopen(iname,"r");
	if (fp == NULL)
		system_error("No such vector file.");
	
	get_header(fp,&Norb,&Natom);

	fclose(fp);

	/* Compute # cubes, # prim cells */
	Ncube = Natom/(2*FCC_P_CUBE);
	Nc = Natom/2;

	/* Test print */
	printf("Norb = %3d  Natom = %3d  Ncube = %3d  Nc = %3d\n",Norb,Natom,Ncube,Nc);

	/* Compute # bands and records in nemo3d eigenvector */
	maxbands = 2*Norb;
	nrec = 2*Nc*Norb;

	/* Get nemo3d structure dimensions */
	printf("Enter # nemo3d cubes in x, y, z> ");
	scanf("%d%d%d",&n_nemo[0],&n_nemo[1],&n_nemo[2]);
	
	/* Run checks on these dimensions */
	if (Nc != FCC_P_CUBE*n_nemo[0]*n_nemo[1]*n_nemo[2])
		printf("# nemo3d cubes disagrees with output vector header!\n");

	/* Now allocate storage */
	Gv = R2tensor(0,Nc-1,0,2);//allocate_real_matrix(0,Nc-1,0,2);
	Rv = R2tensor(0,Nc-1,0,2);//allocate_real_matrix(0,Nc-1,0,2);
	Ur = R2tensor(0,Nc-1,0,Nc-1);//allocate_real_matrix(0,Nc-1,0,Nc-1);
	Ui = R2tensor(0,Nc-1,0,Nc-1);//allocate_real_matrix(0,Nc-1,0,Nc-1);
	vr = Rvectr(nrec/*-1*/);//allocate_real_vector(0,nrec-1);
	vi = Rvectr(nrec/*-1*/);//allocate_real_vector(0,nrec-1);
        vr_tr = Rvectr(nrec/*-1*/);//allocate_real_vector(0,nrec-1);
        vi_tr = Rvectr(nrec/*-1*/);//allocate_real_vector(0,nrec-1);
	wvr = Rvectr(Nc/*-1*/);//allocate_real_vector(0,Nc-1);
	wvi = Rvectr(Nc/*-1*/);//allocate_real_vector(0,Nc-1);
	sel_energies = Rvectr(TR*negy);
	Ekmat = R2tensor(0,TR*negy-1,0,Nc-1);//allocate_real_matrix(0,negy-1,0,Nc-1);

        /* copy all_energies to sel_energies and duplicate energies for time reversal states */
        j=0;
        for(i=min_ix-1;i<=max_ix-1;i++)
        { sel_energies[j]=all_energies[i];
          j++;
          sel_energies[j]=all_energies[i];
          j++;
        }

	/* Compute small-cell k */
	gen_SmCell_k(Gv,n_nemo);

	/* Print available small cell k to file */
	fp = fopen("SmCell_k.txt","w");
	if (fp == NULL)
		system_error("Cannot open k file for write.");

	for(i=0;i<Nc;i++)
	{
		for(j=0;j<3;j++)
			tmp_vec[j] = Gv[i][j]/M_PI;
		fprintf(fp,"%6d   %12.5le   %12.5le   %12.5le\n",
			i,tmp_vec[0],tmp_vec[1],tmp_vec[2]);
	}  /* for i */

	fclose(fp);

	/* Generate R vectors */
	gen_Rvecs(Rv,n_nemo);

	/* Construct projection matrix */
	projection_matrix(Gv,Rv,Ur,Ui,Nc);

printf("\nTEST: before unitarity check\n");

	/* Check unitarity and report */
	is_unitary = check_unitarity(Ur,Ui,Nc,EPS_UNITARY);

printf("\nTEST: after unitarity check\n");

	if (is_unitary)
		printf("Unitarity check passed.\n");
	else
		printf("Unitarity check failed.\n");

	/* The main event:  Read in nemo3d ASCII eigenvector	*/
	/* Project out probabs, store in Ekmat, get another		*/
	for(i=min_ix;i<=max_ix;i++)
	{
	   for(tr=0;tr<TR;tr++)
	      {
		/* Construct file-name of relevant eigenvector */
		full_filename(template1,i,suffix,iname);
		fp = fopen(iname,"r");
		if (fp == NULL)
			system_error("No such input file.");

		/* nemo3d vector */
		get_nemo3d_vector(fp,vr,vi,Ncube,Norb);
	
		/* Energy index */
		eix = i - min_ix;
		eix_tr = TR*(eix)+tr;

                   /* Compute Time Reversal State */
                   /* Note that this time reversal state is eigenvector at k=-k.
                      This change of sign of k will not matter for k=0 */

                   if (tr==1)
                      {
                       /* compute time reversal state */
                       compute_tr_state(vr,vi,vr_tr,vi_tr,Norb,nrec);

                       /* copy Wr_tr, Wi_tr to Wr, Wi respectively */
                       for (int ii=1;ii<=nrec;ii++)
                           {
                             vr[ii]=vr_tr[ii];
                             vi[ii]=vi_tr[ii];
                           } /* for i */
                      } /* if tr */


		/* Project out this energy eigenvector */
		project_probs(Ekmat,Ur,Ui,vr,vi,wvr,wvi,Ncube,Norb,Nc,TR*negy,eix_tr);

		/* We are done with this file */
		fclose(fp);
	      } /* for tr */
	   printf("done with eigenvector %d (i.e. %d out of %d eigenvectors)\n",i,i-min_ix+1,max_ix-min_ix+1);
	}  /* for i */

	printf("projected OK\n");

	/* Save Ekmatrix, k-vectors, and energies to binary file */
	fp = fopen("Projected_states.bin","wb");
	if (fp == NULL)
		system_error("Cannot open output binary file for write");
		
	fwrite(&Nc,sizeof(int),1,fp);
	fwrite(&negy,sizeof(int),1,fp);
	fwrite(&maxbands,sizeof(int),1,fp);
	
	/* We can do the entire list of energies in one fell swoop */
	fwrite(sel_energies,sizeof(double),TR*negy,fp); //NOTE NOTE

	/* Write out matrix of k-vectors */
	for(i=0;i<Nc;i++)
		for(j=0;j<3;j++)
			fwrite(&(Gv[i][j]),sizeof(double),1,fp);
	
	/* Now write out Ekmatrix */
	for(i=0;i<TR*negy;i++)
		for(j=0;j<Nc;j++)
			fwrite(&(Ekmat[i][j]),sizeof(double),1,fp);

	/* We are done finally! */
	fclose(fp);

	/* Free up storage */
	rm_r2tensor(&Gv,0);//free_real_matrix(Gv,0,Nc-1,0);
	rm_r2tensor(&Rv,0);//free_real_matrix(Rv,0,Nc-1,0);
	rm_r2tensor(&Ur,0);//free_real_matrix(Ur,0,Nc-1,0);
	rm_r2tensor(&Ui,0);//free_real_matrix(Ui,0,Nc-1,0);
	rm_rvectr(&all_energies);//free_real_vector(all_energies,0);
	rm_rvectr(&vr);//free_real_vector(vr,0);
	rm_rvectr(&vi);//free_real_vector(vi,0);
	rm_rvectr(&vr_tr);//free_real_vector(vr,0);
        rm_rvectr(&vi_tr);//free_real_vector(vi,0);
	rm_rvectr(&wvr);//free_real_vector(wvr,0);
	rm_rvectr(&wvi);//free_real_vector(wvi,0);
	rm_rvectr(&sel_energies);
	rm_r2tensor(&Ekmat,0);//free_real_matrix(Ekmat,0,negy-1,0);

	return(0);
}  /* fccproj */

