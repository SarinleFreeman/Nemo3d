#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>
#include <time.h>

#include "r2tensor.h"
#include "rvector.h"
#include "ivector.h"

#include "FCCproj.h"

//#include "/home/boykin/lib/utility.c"

#define system_error printf

#define	FCC_P_CUBE	4		/* # atoms in NEMO-3D cube */
#define NEMO_HEAD	7		/* # lines to toss in NEMO3D ASCII header */
#define EPS_UNITARY	1.0e-10	/* effective 0 for unitary test */
#define EPS		1.0e-10
#define TR 1 /* All energies are doubly degenerate due to time reversal symmetry */
#define A_SM_SC 1.0

#define NEMO_HEAD_rAtom 6

/* Get energies from nemo3d Ek ASCII file pointed to by fp	*/
/* Get only those energies GREATER THAN emin and LESS THAN	*/
/* emax. Use all_energies[0..max_states-1] return SUBVECTOR	*/ 
/* energies[0..*negy-1] setting a pointer to start energy	*/
/* in list of all also return negy = #energy states found &	*/
/* min_ix and max_ix, indices of minimum and maximum vects	*/
/* NOTE THAT VECTORS ARE UNIT-INDEXED.						*/
void get_nemo3d_energies_cubic(FILE *fp, double emin, double emax, double *all_energies,
				int max_states, int *negy, int *min_ix, int *max_ix, 
                                double *k)                                           /* RENAMED AND MODIFIED TO READ K*/
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
		fscanf(fp,"%lf%lf%lf%lf",&kx,&ky,&kz,&egy);
                k[0]=kx; k[1]=ky; k[2]=kz;
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

	if (*max_ix == -1)
		*max_ix = max_states;

	/* Now compute negy */
	*negy = (*max_ix - *min_ix + 1);

}  /* get_nemo3d_energies */
	
		
/* Read complex nemo3d vector from nemo3d ASCII file pointed to by fp   */
/* Return in vector Wr, Wi[0..nrec-1].  Discard header (1st 7 lines)    */
void get_nemo3d_vector_cubic(FILE *fp, double vr[], double vi[])       
{
        int i,j;
        double re_pt, im_pt;
        char s[8000];

        /* Toss header -- we already looked at it to determine nrec */
        for(j=1;j<=NEMO_HEAD;j++)
            hoover(fp);
        i=1;
        while (fgets(s,8000,fp)!=NULL)
              {
                sscanf(s,"%lf%lf", &re_pt, &im_pt);
                vr[i] = re_pt;
                vi[i] = im_pt;
                i++;
              } /* while */
}  /* get_nemo3d_vector_cubic */



void gen_Rvecs_cubic(double **R_a, int Nrect[])  /* renamed */
                                                 /* Note that first index of R_a starts from 0 */
{
	int i, nemo[3];
	int nemo_ix, aix;

	/* Now loop through the NEMO cubes */
	for(nemo[0]=0;nemo[0]<Nrect[0];nemo[0]++)
		for(nemo[1]=0;nemo[1]<Nrect[1];nemo[1]++)
			for(nemo[2]=0;nemo[2]<Nrect[2];nemo[2]++)
			{
				/* Compute NEMO index */
				nemo_ix = ((nemo[0]*Nrect[1] + nemo[1])*Nrect[2] + nemo[2]);//*FCC_P_CUBE; // NO NEED TO MULTIPLY BY FCC_P_CUBE
				aix = nemo_ix;
 
                                for(i=0;i<3;i++)
                                	R_a[aix][i] = nemo[i];
				
			}  /* for nemo[2] */

}  /* gen_Rvecs */




void gen_SmCell_kpGvecs_cubic(double **kpGvecs, double *k, int *Nrect, double bz_max) /* renamed */
                                                   /* Note that first index of kpGvecs starts from 0 */
{
        int i,j,Nmin[3],Nmax[3],ix,iy,iz,iix,iiy,iiz,vindx;
        double Gv[3],kpG[3];
        double **Bvec;

        Bvec = R2tensor(0,2,0,2);

        /* Assign Supercell Recip Lattice Vectors; for SC they  */
        /* are one-component each and orthogonal.                               */
        /* Each vector is stored as a COLUMN of Bvec.                   */
        double twoPI = 2.0*M_PI;
        for(i=0;i<3;i++)
           {
             for(j=0;j<3;j++)
                Bvec[j][i] = 0.0;

                /* Units are inverse small cell cube edge */
                Bvec[i][i] = twoPI/((double)Nrect[i]);
           }  /* for i */

        /* Compute min and max integers for Supercell Recip vects */
        for(i=0;i<3;i++)
           {
                if (Nrect[i] % 2)
                {
                        Nmax[i] = (Nrect[i]-1)/2;
                        Nmin[i] = -Nmax[i];
                }
                else
                {
                        Nmin[i] = -(Nrect[i]-2)/2;
                        Nmax[i] = Nrect[i]/2;
                }
           }  /* for i */

        /* Now compute bulk vectors for each recip latt vect */
        for(ix=Nmin[0];ix<=Nmax[0];ix++)
            for(iy=Nmin[1];iy<=Nmax[1];iy++)
                for(iz=Nmin[2];iz<=Nmax[2];iz++)
                   {
                        iix = ix - Nmin[0];
                        iiy = iy - Nmin[1];
                        iiz = iz - Nmin[2];

                        /* For SC dot prod is trivial so do it manually */
                        Gv[0] = ((double)ix)*Bvec[0][0];        /* + iy*Bvec[1][0] + iz*Bvec[2][0] */
                        Gv[1] = ((double)iy)*Bvec[1][1];        /* + ix*Bvec[0][1] + iz*Bvec[2][1] */
                        Gv[2] = ((double)iz)*Bvec[2][2];        /* + ix*Bvec[0][2] + iy*Bvec[1][2] */

                        for(i=0;i<3;i++)
                           {
                               kpG[i] = k[i] + Gv[i];

                               /* Test if out of1st BZ; shift back if it is */
                               if (kpG[i] > bz_max)
                                   kpG[i] -= 2.0*bz_max;
                                   else if (kpG[i] < -bz_max)
                                   kpG[i] += 2.0*bz_max;
                           }  /* for i */

			vindx = ((iix*Nrect[1] + iiy)*Nrect[2] + iiz);

                        for(i=0;i<3;i++)
                            kpGvecs[vindx][i] = kpG[i];///M_PI;
                     }  /* for iz*/
        rm_r2tensor(&Bvec,0);
} /* gen_SmCell_kpGvecs */



/* Project out probabilities for wavevectors (as ordered by gen_kpGvec) from vector	*/
/* vr, vi[0..nrec], where nrec = 2*Nc*FCC_P_CUBE*Norb. Work-space               	*/
/* is the vector wvr, wvi[0..nrec-1].  Store result (after projecting out and summing	*/
/* over all orbitals) as row eix of matrix Ekmat[0..Negy-1][0..Nc-1] where Negy =	*/
/* #energy eigenvectors we read in from the nemo3d output.  Cols of Ekmat are 		*/
/* prim cell k's (ie supercell kpGvecs's) with order as determined by gen_kpGvecs.      */
/* Use unitary projection matrix Ur, Ui[0..Nc-1][0..Nc-1] to do projection.		*/
void project_probs_cubic(double **Ekmat, double **Ur, double **Ui, double vr[], double vi[],
double wvr[], double wvi[], int NBASIS, int Norb, int Nc_nano, int Negy, int eix, int *NBASIS_nano,
double k[], int *n_nanostruct) /* int NBASIS - not required */
{
	int i, j, orb;
        double norm, kRdot, exr, exi, tr, ti;
        int n, ii, jj, state, Gix;
        int ix, iy, iz, Rv[3], row, row_rv, atom;
	int *Norbnano, *CNorbnano;

	/* Zero row eix of Ekmat since we deal with eix-th energy only */
	for(j=0;j<Nc_nano;j++)
		Ekmat[eix][j] = 0.0;

	/* Compute orbitals per unit cell of nanowire also compute cumulative orbitals for indexing */
	Norbnano = Ivectr(Nc_nano);
	CNorbnano = Ivectr(Nc_nano+1);

	Norbnano[0] = NBASIS_nano[0]*Norb;
	CNorbnano[0] = 0;
	for(j=1;j<Nc_nano;j++)
	{
		Norbnano[j] = NBASIS_nano[j]*Norb;
		CNorbnano[j] = CNorbnano[j-1] + Norbnano[j-1];
	}
	CNorbnano[Nc_nano] = CNorbnano[Nc_nano-1] + Norbnano[Nc_nano-1];

        /* First dephase the matrix W for this state all orbs */
        for(ix=0;ix<n_nanostruct[0];ix++)
        {
                Rv[0] = ix;
                for(iy=0;iy<n_nanostruct[1];iy++)
                {
                        Rv[1] = iy;
                        for(iz=0;iz<n_nanostruct[2];iz++)
                        {
                                Rv[2] = iz;
                                row_rv = Rv[2] + n_nanostruct[2]*Rv[1] + n_nanostruct[2]*n_nanostruct[1]*Rv[0];

                                /* Compute k.R */
                                kRdot = 0.0;
                                for(i=0;i<=2;i++)
                                        kRdot += ((double)Rv[i])*k[i];

                                /* Compute exp(-i*k.R); remember cosine is EVEN! */
                                exr = cos(kRdot);
                                exi = sin(-kRdot);

                                /* Dephase; replace entry in matrix by dephased entry */
				for(int atomorb=0;atomorb<Norbnano[i];atomorb++)
				{
					row = CNorbnano[ix] + atomorb;
   	                                tr = exr*vr[row] - exi*vi[row];
                                        ti = exr*vi[row] + exi*vr[row];
                                        vr[row] = tr;
                                        vi[row] = ti;
                                } /* for atomorb */
                        }  /* for ix */
                }  /* for iy */
        }  /* for iz */

        /* Do the projection and normalization one state at a time */
        /* Proceed by orbital type on projections */
        /* Project out with herm conj of U */
        for(i=0;i<Nc_nano;i++)
        {
		for(int atomorb=0;atomorb<Norbnano[i];atomorb++)
		{
                	/* ii is index in v, wv */
			ii = CNorbnano[i] + atomorb;
                        wvr[ii] = 0.0;
                        wvi[ii] = 0.0;

                        for(j=0;j<Nc_nano;j++)
                        {
                             	/* jj is index in v */
				jj = CNorbnano[j] + atomorb;
                                wvr[ii] += (Ur[j][i]*vr[jj] + Ui[j][i]*vi[jj]);
                                wvi[ii] += (Ur[j][i]*vi[jj] - Ui[j][i]*vr[jj]);
                        }  /* for j */
                }  /* for atomorb */
	}  /* for i */

        /* Now normalize bulk coeffs for each G_n for this state */
        for(i=0;i<Nc_nano;i++)
        {
                norm = 0.0;
		for(int atomorb=0;atomorb<Norbnano[i];atomorb++)
		{
			ii = CNorbnano[i] + atomorb;
                        norm += wvr[ii]*wvr[ii] + wvi[ii]*wvi[ii];
                } /* for atomorb */

		Ekmat[eix][i] = norm;
	} /* for i */

	rm_ivectr(&Norbnano);
	rm_ivectr(&CNorbnano);

}  /* project_probs_cubic */


void define_small_cell_of_nanostructure(int *n_nemo, int *n_nanostruct, int *NBASIS, int *periodic_y, int *periodic_z)
{
        if(*periodic_y==1 && *periodic_z==1) /* Bulk */
        { 
		n_nanostruct[0] = n_nemo[0]; n_nanostruct[1] = n_nemo[1]; n_nanostruct[2] = n_nemo[2]; 
          	*NBASIS = 2*FCC_P_CUBE; 
	}

        if(*periodic_y==1 && *periodic_z==0) /* Quantum well confined in z direction */
        { 
		n_nanostruct[0] = n_nemo[0]; n_nanostruct[1] = n_nemo[1]; n_nanostruct[2] = 1;           
          	*NBASIS = 2*FCC_P_CUBE*n_nemo[2]; 
	}

        if(*periodic_y==0 && *periodic_z==1) /* Quantum well confined in y direction */
        { 
		n_nanostruct[0] = n_nemo[0]; n_nanostruct[1] = 1; n_nanostruct[2] = n_nemo[2];           
	        *NBASIS = 2*FCC_P_CUBE*n_nemo[1]; 
	}

        if(*periodic_y==0 && *periodic_z==0) /* Nanowire confined in y and z directions */
        { 
		n_nanostruct[0] = n_nemo[0]; n_nanostruct[1] = 1; n_nanostruct[2] = 1;           
          	*NBASIS = 2*FCC_P_CUBE*n_nemo[1]*n_nemo[2]; 
	}
} /* define_small_cell_of_nanostructure */


void read_nemo3d_atom_pos_gen_NBASIS_nano(FILE *fp,int *NBASIS_nano)
{
        int i;
        double xo,yo,zo,xn,yn,zn;
        char s[8000];

        for(i=1;i<=NEMO_HEAD_rAtom;i++)
           hoover(fp);

        fgets(s,8000,fp);
        sscanf(s,"%lf %lf %lf",&xo,&yo,&zo);

        i=0;
        NBASIS_nano[0]=1;
        while (fgets(s,8000,fp)!= NULL)
        { sscanf(s,"%lf %lf %lf",&xn,&yn,&zn);
          if (xn>xo)
             { i++;
               xo=xn; yo=yn; zo=zn;
               NBASIS_nano[i]=1;
             } /* if */
          else
             NBASIS_nano[i]++;
        } /* while */
} /* read_nemo3d_atom_pos_gen_NBASIS_nano */
					

int cubicproj(void) 
{
	int min_ix, max_ix, negy, max_states, err, i, j, Ncube, Norb, Natom;
	int Nc, num, maxbands, eix, kix, nrec, is_unitary, n_nemo[3], tr, eix_tr;
	double emin, emax, *energies, *all_energies, k[3], **kpGvecs, **Rv, **Ur, **Ui;
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
	get_nemo3d_energies_cubic(fp,emin,emax,all_energies,max_states,&negy,
					&min_ix,&max_ix,k);

	fclose(fp);

	printf("min = %3d  max = %3d\n",min_ix,max_ix);

	/* Assign energies */
	energies = &(all_energies[min_ix-1]);

	/* Test print */
	for(i=0;i<negy;i++)
	{
		j = i + min_ix;
		printf("%3d   %19.12le\n",j,energies[i]); 
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

	/* Compute # cubes, # fcc (cubic) unit cells */
	Ncube = Natom/(2*FCC_P_CUBE);
	//Nc = Natom/2;
        Nc = Ncube;

	/* Test print */
	printf("Norb = %3d  Natom = %3d  Ncube = %3d  Nc = %3d\n",Norb,Natom,Ncube,Nc);

	/* compute records in nemo3d eigenvector */
	//nrec = 2*Nc*FCC_P_CUBE*Norb; /* MODIFIED */

	/* Get nemo3d structure dimensions */
	printf("Enter # nemo3d cubes in x, y, z> ");
	scanf("%d%d%d",&n_nemo[0],&n_nemo[1],&n_nemo[2]);
	
	/* Run checks on these dimensions */
	if (Nc != FCC_P_CUBE*n_nemo[0]*n_nemo[1]*n_nemo[2])
		printf("# nemo3d cubes disagrees with output vector header!\n");

        /* Get periodicity info */
        char res;
        int periodic_y=0, periodic_z=0;
 
        printf("\nPeriodic in y? (y/n) > ");
        scanf("%s",&res);
        if(res == 'y' || res == 'Y')  periodic_y=1;
        else periodic_y=0;

        printf("Periodic in z? (y/n) > ");
        scanf("%s",&res);
        if(res == 'y' || res == 'Y')  periodic_z=1;
        else periodic_z=0;

        /* Adjust the unit cell size and basis size according to the periodicity of supercell */
        /* n_nanostruct[i]'s are the number of unit cells in direction i, these are defined taking into */
        /* account the periodicity of the lattice */
	int NBASIS=0, n_nanostruct[3], Nc_nano;
        define_small_cell_of_nanostructure(n_nemo,n_nanostruct,&NBASIS,&periodic_y,&periodic_z);
        Nc_nano = n_nanostruct[0]*n_nanostruct[1]*n_nanostruct[2];

        /* read nemo3d rAtom file and generate NBASIS_nano[1...Nc_nano] */
	int *NBASIS_nano;
	NBASIS_nano = Ivectr(Nc_nano);
	char rname[81];
	/* Get rAtom input file */
	printf("Enter name of rAtom ASCII file> ");
	scanf("%s",rname);
	/* Attempt open */
	fp = fopen(rname,"r");
	if (fp==NULL)
		system_error("No such atom position input file.");

	/* Read rAtom file and generate atom position data into unit_cells_nanostrcture */
	read_nemo3d_atom_pos_gen_NBASIS_nano(fp,NBASIS_nano);

	fclose(fp);

	/* Compute length of eigenvector */
	long natoms=0;
	for(i=0;i<Nc_nano;i++) 
	natoms += NBASIS_nano[i]; /* natoms = total number of atoms */
	nrec=natoms*Norb; /* length of eigenvector */


	/* Now allocate storage */
	kpGvecs = R2tensor(0,Nc_nano-1,0,2);//allocate_real_matrix(0,Nc-1,0,2);
	Rv = R2tensor(0,Nc_nano-1,0,2);//allocate_real_matrix(0,Nc-1,0,2);
	Ur = R2tensor(0,Nc_nano-1,0,Nc_nano-1);//allocate_real_matrix(0,Nc-1,0,Nc-1);
	Ui = R2tensor(0,Nc_nano-1,0,Nc_nano-1);//allocate_real_matrix(0,Nc-1,0,Nc-1);
	vr = Rvectr(nrec/*-1*/);//allocate_real_vector(0,nrec-1);
	vi = Rvectr(nrec/*-1*/);//allocate_real_vector(0,nrec-1);
        vr_tr = Rvectr(nrec/*-1*/);//allocate_real_vector(0,nrec-1);
        vi_tr = Rvectr(nrec/*-1*/);//allocate_real_vector(0,nrec-1);
//	wvr = Rvectr(Nc_nano/*-1*/);//allocate_real_vector(0,Nc-1);
//	wvi = Rvectr(Nc_nano/*-1*/);//allocate_real_vector(0,Nc-1);
	wvr = Rvectr(nrec/*-1*/);//allocate_real_vector(0,Nc-1);
	wvi = Rvectr(nrec/*-1*/);//allocate_real_vector(0,Nc-1);
	sel_energies = Rvectr(TR*negy);
	Ekmat = R2tensor(0,TR*negy-1,0,Nc-1);//allocate_real_matrix(0,negy-1,0,Nc-1);

        /* copy all_energies to sel_energies and duplicate energies for time reversal states */
        j=0;
        for(i=min_ix-1;i<=max_ix-1;i++)
        { sel_energies[j]=all_energies[i];
          j++;
//          sel_energies[j]=all_energies[i];
//        j++;
        }

	/* Compute small-cell k */
        /* For simple cubic sm cell size is conv unit cell cube edge */
        double a_sm = A_SM_SC;
        double bz_max = M_PI/a_sm;
        /* scale k according to Nc in x,y,z */
         for(i=0;i<3;i++) k[i] *= M_PI/n_nanostruct[i];

        gen_SmCell_kpGvecs_cubic(kpGvecs,k,n_nanostruct,bz_max); /* renamed */

	/* Print available small cell k to file */
	fp = fopen("SmCell_k.txt","w");
	if (fp == NULL)
		system_error("Cannot open k file for write.");

	for(i=0;i<Nc_nano;i++)
	{
		for(j=0;j<3;j++)
			tmp_vec[j] = kpGvecs[i][j]/M_PI;
		fprintf(fp,"%6d   %12.5le   %12.5le   %12.5le\n",
			i,tmp_vec[0],tmp_vec[1],tmp_vec[2]);
	}  /* for i */

	fclose(fp);

	/* Generate R vectors */
	gen_Rvecs_cubic(Rv,n_nanostruct);

	/* Construct projection matrix */
	projection_matrix(kpGvecs,Rv,Ur,Ui,Nc_nano);

	/* Check unitarity and report */
	is_unitary = check_unitarity(Ur,Ui,Nc_nano,EPS_UNITARY);

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
		get_nemo3d_vector_cubic(fp,vr,vi); /* an argument specifying the size of eigenvector should be passed */
	
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
		project_probs_cubic(Ekmat,Ur,Ui,vr,vi,wvr,wvi,NBASIS,Norb,Nc_nano,TR*negy,eix_tr,NBASIS_nano,k,n_nanostruct); 

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
		
	fwrite(&Nc_nano,sizeof(int),1,fp);
	fwrite(&negy,sizeof(int),1,fp);

        /* compute # of bands, actual number of bands found from projection algorithm will always be less than maxbands */
        maxbands = NBASIS*Norb; 

	fwrite(&maxbands,sizeof(int),1,fp);  
	
	/* We can do the entire list of energies in one fell swoop */
	fwrite(sel_energies,sizeof(double),TR*negy,fp);

	/* Write out matrix of k-vectors */
	for(i=0;i<Nc_nano;i++)
		for(j=0;j<3;j++)
			fwrite(&(kpGvecs[i][j]),sizeof(double),1,fp);
	
	/* Now write out Ekmatrix */
	for(i=0;i<TR*negy;i++)
		for(j=0;j<Nc_nano;j++)
			fwrite(&(Ekmat[i][j]),sizeof(double),1,fp);

	/* We are done finally! */
	fclose(fp);

	/* Free up storage */
	rm_r2tensor(&kpGvecs,0);//free_real_matrix(Gv,0,Nc-1,0);
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
	rm_ivectr(&NBASIS_nano);

	return(0);
}  /* cubicproj */

