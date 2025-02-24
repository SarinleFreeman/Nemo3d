#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>

#include "r2tensor.h"
#include "rvector.h"
#include "ivector.h"

#include "FCCguess.h"

#define system_error printf

//#include "/home/boykin/lib/utility.c"

#define	FCC_P_CUBE	4		/* # atoms in NEMO-3D cube */
#define NEMO_HEAD	7		/* # lines to toss in NEMO3D ASCII header */
#define EPS_UNITARY	1.0e-10	/* effective 0 for unitary test */
#define EPS		1.0e-10
#define EPS_K		1.0e-6
#define TR 1 /* All energies are doubly degenerate due to time reversal symmetry */



/* Determine if vector k[0..2] and k_unit[0..2], 0 -> x, 1 -> y, 2 -> z */
/* lie along the same direction or an equivalent direction for simple   */
/* cubic.  (For example, the bands along [110] and [-110] must be the   */
/* same in simple cubic.)  k_unit is of unit length; k has arbitrary    */
/* length.  Use eps (a small value) as tolerance on the comparison.             */
/* Return +1 if they are along the same direction, 0 if not.                    */
/* Uses fact that if k and k_unit are parallel the projection of k onto */
/* k_unit is the same as the magnitude of k.                                                    */
int comp_dir_cubic(double k[], double k_unit[], double eps)
{
        int i, same;
        double k_test[3], kmag, proj, temp;

        /* First attempt to align by choosing same signs on unit vector */
        for(i=0;i<=2;i++)
        {
                temp = fabs(k_unit[i]);
                if (k[i] < 0.0)
                        k_test[i] = -temp;
                else
                        k_test[i] = temp;
        }  /* for i */

        /* Now compute projection of k onto k_test and compare to mag k */
        kmag = 0.0;
        proj = 0.0;

        for(i=0;i<=2;i++)
        {
                kmag += k[i]*k[i];
                proj += k[i]*k_test[i];
        }  /* for i */

        kmag = sqrt(kmag);

        /* Compare directions to within eps */
        if (fabs(proj - kmag) <= eps)
                same = 1;
        else
                same = 0;

        return(same);
}  /* comp_dir_cubic */




int cubicguess(void)
{
	int i, j, Nc, same, negy, maxbands, kix, *nb;
	double **Ekmat, **begymat, **bsdevmat, **bprobmat, *energies, **kpGvecs;
	double kmag, control[3], k_unit[3], kv[3];
	FILE *fp, *fp1, *fp2, *fp3, *fp4;

	/* Get binary input file */
	fp = fopen("Projected_states.bin","rb");
	if (fp == NULL)
		system_error("No such binary input file.");

	/* Get Nc, #energies, max # bands */
	fread(&Nc,sizeof(int),1,fp);
	fread(&negy,sizeof(int),1,fp);
	fread(&maxbands,sizeof(int),1,fp);

	/* Allocate storage */
	kpGvecs = R2tensor(0,Nc-1,0,2);//allocate_real_matrix(0,Nc-1,0,2);
	energies = Rvectr(TR*negy/*-1*/);//allocate_real_vector(0,negy-1);
	begymat = R2tensor(0,maxbands-1,0,Nc-1);//allocate_real_matrix(0,maxbands-1,0,Nc-1);
	bsdevmat = R2tensor(0,maxbands-1,0,Nc-1);//allocate_real_matrix(0,maxbands-1,0,Nc-1);
	bprobmat = R2tensor(0,maxbands-1,0,Nc-1);//allocate_real_matrix(0,maxbands-1,0,Nc-1);
	Ekmat = R2tensor(0,TR*negy-1,0,Nc-1);//allocate_real_matrix(0,negy-1,0,Nc-1);
	nb = Ivectr(Nc/*-1*/);//allocate_integer_vector(0,Nc-1);
	
	/* Get energy vecor */
	fread(energies,sizeof(double),TR*negy,fp);

	/* Get k-vectors */
	for(i=0;i<Nc;i++)
		for(j=0;j<3;j++)
			fread(&(kpGvecs[i][j]),sizeof(double),1,fp);
	
	/* Now read Ekmatrix */
	for(i=0;i<TR*negy;i++)
		for(j=0;j<Nc;j++)
			fread(&(Ekmat[i][j]),sizeof(double),1,fp);
	fclose(fp);
 
        /* Write ProbvsE data at chosen k_small_cell to a file */
        fp = fopen("ProbvsE.txt","w");
        for(j=0;j<Nc;j++) {
            if(j==22) // choose == # for diff small cell ks
               for(i=0;i<TR*negy;i++)
                   fprintf(fp,"%12.10le %12.10le\n",energies[i],Ekmat[i][j]); 
           }
        fclose(fp);

	/* Get control parameters for band guessing */
	printf("Enter min gap, tiny prob, min band prob> ");
	scanf("%lf%lf%lf",&control[0],&control[1],&control[2]);

	/* Now guess bands; here loop over small cell k, indexed by kix */
	for(kix=0;kix<Nc;kix++)
		band_guess(energies,Ekmat,control,TR*negy,
					begymat,bsdevmat,bprobmat,nb,kix);

	/* Finally print out these results */
	fp = fopen("band_guesses.txt","w");
	if (fp == NULL)
		system_error("Cannot open output for writing.");

	for(kix=0;kix<Nc;kix++)
	{
		fprintf(fp,"kix = %3d:  k = (%12.5le, %12.5le, %12.5le)(Pi/a):  nb = %3d\n\n",
				kix,kpGvecs[kix][0]/M_PI,kpGvecs[kix][1]/M_PI,kpGvecs[kix][2]/M_PI,nb[kix]);
		/* Loop thorough bands found at each sm cell k */
		for(j=0;j<nb[kix];j++)
			fprintf(fp,"%19.12le     %12.5le     %7.5lf\n",begymat[j][kix],
				bsdevmat[j][kix],bprobmat[j][kix]);

		fprintf(fp,"\n\n");
	}  /* for kix */

	/* We are done */
	fclose(fp);

	/* Now write out energies, sdev, prob for [100], [110], [111] */
	fp1 = fopen("Ek_all.txt","w");
	fp2 = fopen("Sdev_all.txt","w");
	fp3 = fopen("Prob_all.txt","w");
	fp4 = fopen("Ek_plot.txt","w");
	
	for(i=0;i<3;i++)
	{
		/* Set unit vector for comparison */
		for(j=0;j<=i;j++)
			k_unit[j] = 1.0/sqrt((double)i+1.0);

		for(j=i+1;j<3;j++)
			k_unit[j] = 0.0;

		/* Print out unit vectors to each file */
		fprintf(fp1,"k_unit = (%12.5le, %12.5le, %12.5le)\n\n",
			k_unit[0],k_unit[1],k_unit[2]);
		fprintf(fp2,"k_unit = (%12.5le, %12.5le, %12.5le)\n\n",
			k_unit[0],k_unit[1],k_unit[2]);		
		fprintf(fp3,"k_unit = (%12.5le, %12.5le, %12.5le)\n\n",
			k_unit[0],k_unit[1],k_unit[2]);
		fprintf(fp4,"k_unit = (%12.5le, %12.5le, %12.5le)\n\n",
                        k_unit[0],k_unit[1],k_unit[2]);
                fprintf(fp4,"   kmag          Energy        Esdev         bprob\n");
               
		for(kix=0;kix<Nc;kix++)
		{
			/* Copy a vector, get mag, and test */
			kmag = 0.0;
			for(j=0;j<3;j++)
			{
				kv[j] = kpGvecs[kix][j];
				kmag += kv[j]*kv[j];
			}  /* for j */

			/* We print out magnitude of k normalized to PI */
			kmag = sqrt(kmag)/M_PI;

			/* Always count (0,0,0) as a vaild state */
			if (kmag < EPS_K)
				same = 1;
			else
				same = comp_dir_cubic(kv,k_unit,EPS_K);
			if (same)
			{
				fprintf(fp1,"%12.5le ",kmag);
				fprintf(fp2,"%12.5le ",kmag);
				fprintf(fp3,"%12.5le ",kmag);

				for(j=0;j<nb[kix];j++)
				{
					fprintf(fp1,"  %12.5le",begymat[j][kix]);
					fprintf(fp2,"  %12.5le",bsdevmat[j][kix]);
					fprintf(fp3,"  %12.5le",bprobmat[j][kix]);
					fprintf(fp4,"  %12.5le  %12.5le  %12.5le  %12.5le\n",
                                                kmag,begymat[j][kix],bsdevmat[j][kix],bprobmat[j][kix]);
				}  /* for j */
				
				fprintf(fp1,"\n");
				fprintf(fp2,"\n");
				fprintf(fp3,"\n");
			}  /* if (same) */
		}  /* for kix */

		fprintf(fp1,"\n");
		fprintf(fp2,"\n");
		fprintf(fp3,"\n");
	}  /* for i */

	/* Close files */
	fclose(fp1);
	fclose(fp2);
	fclose(fp3);
	fclose(fp4);

	/* Free up storage */
	rm_r2tensor(&kpGvecs,0);//free_real_matrix(Gv,0,Nc-1,0);
	rm_rvectr(&energies);//free_real_vector(energies,0);
	rm_r2tensor(&begymat,0);//free_real_matrix(begymat,0,maxbands-1,0);
	rm_r2tensor(&bprobmat,0);//free_real_matrix(bprobmat,0,maxbands-1,0);
	rm_r2tensor(&bsdevmat,0);//free_real_matrix(bsdevmat,0,maxbands-1,0);
	rm_r2tensor(&Ekmat,0);//free_real_matrix(Ekmat,0,negy-1,0);
	rm_ivectr(&nb);//free_integer_vector(nb,0);

	return(0);
}  /* fccguess */
	
