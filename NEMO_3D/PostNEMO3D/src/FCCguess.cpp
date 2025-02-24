#include <stdio.h>
#include <stdlib.h>
//#include <math.h>
#include <string.h>

#include "r2tensor.h"
#include "rvector.h"
#include "ivector.h"

#define system_error printf

//#include "/home/boykin/lib/utility.c"

#define	FCC_P_CUBE	4		/* # atoms in NEMO-3D cube */
#define NEMO_HEAD	7		/* # lines to toss in NEMO3D ASCII header */
#define EPS_UNITARY	1.0e-10	/* effective 0 for unitary test */
#define EPS		1.0e-10
#define EPS_K		1.0e-6
#define TR 2 /* All energies are doubly degenerate due to time reversal symmetry */


/* Determine if vector k[0..2] and k_unit[0..2], 0 -> x, 1 -> y, 2 -> z	*/
/* lie along the same direction or an equivalent direction for FCC		*/
/* (For example, the bands along [110] and [-110] must be the			*/
/* same in FCC.)  k_unit is of unit length; k has arbitrary				*/
/* length.  Use eps (a small value) as tolerance on the comparison.		*/
/* Return +1 if they are along the same direction, 0 if not.			*/
/* Uses fact that if k and k_unit are parallel the projection of k onto	*/
/* k_unit is the same as the magnitude of k.							*/
int comp_dirFCC(double k[], double k_unit[], double eps)
{
	int i, j, same;
	double k_test[3], kmag, proj, temp;

	/* First, copy magnitude of k components and compute length of k */
	kmag = 0.0;
	for(i=0;i<3;i++)
	{
		k_test[i] = fabs(k[i]);
		kmag += k[i]*k[i];
	}  /* for i */
	
	kmag = sqrt(kmag);
	
	/* k_unit is assumed positive; try projecting for 3 permutations of k_test */
	i = 0;
	same = 0;

	while ((i<3) && (!same))
	{
		/* Compute projection of k_test onto k_unit and compare to mag k */
		proj = 0.0;

		for(j=0;j<=2;j++)
			proj += k_unit[j]*k_test[j];

		/* Compare directions to within eps */
		if (fabs(proj - kmag) <= eps)
			same = 1;
		else
		{ 
			/* Permute components of k_test */
			temp = k_test[2];
			k_test[2] = k_test[1];
			k_test[1] = k_test[0];
			k_test[0] = temp;
			i++;
		}  /* else */
	}  /* while */
	
	return(same);
}  /* comp_dirFCC */




/* Guess band energies and std deviations as well as number of degenerate	*/
/* or near-degenerate bands at sm cell k indexed by kix, ordered as with 	*/
/* gen_Gvec.  INPUTS:  energies[0..ns-1],all available energies; 			*/
/* Ekmat[0..ns-1][0..Nc-1]; Ekmat[i][j] = prob for i-th	energy at j-th 		*/
/* small-cell wavevector, ordered as with gen_Gvec; control[0..2] ([0]= 	*/
/* min resolvable gap; [1] = negligible probability,[2] = min cum prob for 	*/
/* a band) which controls band definition.  OUTPUT:  nb[0..nc-1]; nb[kix] = */
/* number of resolvable bands at kix; bprobmat[0..nb[kix]-1][0..Nc-1]; 		*/
/* bprobmat[j][kix] = tot prob for j-th resolvable band at kix; 			*/
/* begymat[0..nb[kix]-1][0..Nc-1]; begymat[j][kix] = mean energy of j-th	*/ 	
/* resolvable band at kix; bsdevmat[0..nb[kix]-1][0..Nc-1];bsdevmat[j][kix]=*/ 
/* std dev of j-th resolvable band at kix.									*/
void band_guess(double energies[], double **Ekmat, double control[], int ns,
		double **begymat, double **bsdevmat, double **bprobmat, 
		int nb[], int kix)
{
	int i, j, nband, in_gap, try_gap, start, i_start;
	double min_gap, tot_prob, tiny_prob, e, e2, gap_start, curr_prob;
	double esum, e2sum, temp, temp2, min_band, e_prev;

	/* Save control parameters */
	min_gap = control[0];
	tiny_prob = control[1];
	min_band = control[2];

	/* Clear out any large samples.  We need to begin within a gap	*/
	nb[kix] = 0;
	i = 0;
	start = 1;
	try_gap = 0;
	in_gap = 0;

	/* OK we have now discarded any partial bands at start of search region */
	/* start a real band accumulation at i									*/
	i_start = i;
	esum = 0.0;
	e2sum = 0.0;
	tot_prob = 0.0;
	e_prev = energies[i_start];
	for(i=i_start;i<ns;i++)
	{
		curr_prob = Ekmat[i][kix];
		e = energies[i];
		e2 = e*e;

		/* Accumulate probs and averages if we are in a band -- not gap	*/
		if (!in_gap)
		{

			/* Only add to current band if samples are closely spaced */
			if (e - e_prev < min_gap)
			{
				tot_prob += curr_prob;
				esum += e*curr_prob;
				e2sum += e2*curr_prob;
			}  /* if (e-e_prev < min_gap) */


			/* Now check if we are about to go into next gap */
			if (curr_prob < tiny_prob)			
			{


				/* If this is 1st tiny prob for this gap start recording */
				if (!try_gap)
				{
					try_gap = 1;
					gap_start = e;
				}  /* if (!try_gap) */
				else
				{
					/* Is this gap bigger than min resolvable gap? */
					if (e-gap_start > min_gap)
					{
						in_gap = 1;

						/* Save statistics -- we are done with previous band */
						/* Only  save if cum prob > min_band */
						if (tot_prob > min_band)
						{
							j = nb[kix];
							temp = esum/tot_prob;
							begymat[j][kix] = temp;
							temp2 = e2sum/tot_prob;
							bsdevmat[j][kix] = sqrt(temp2 - temp*temp);
							bprobmat[j][kix] = tot_prob;
							nb[kix] = j+1;
						}
		
						/* Reset sums for next band */
						esum = 0.0;
						e2sum = 0.0;
						tot_prob = 0.0;
						try_gap = 0;
					}  /* if (e-gap_start > min_gap) */
				}  /* else */
			}  /* if (curr_prob < tiny_prob) */
			else if (e - e_prev > min_gap)
			{
				/* Case of a big gap between big samples: 	*/
				/* Declare a band and move on with it.		*/
				j = nb[kix];
				temp = esum/tot_prob;
				begymat[j][kix] = temp;
				temp2 = e2sum/tot_prob;
				bsdevmat[j][kix] = sqrt(temp2 - temp*temp);
				bprobmat[j][kix] = tot_prob;
				nb[kix] = j+1;

				/* Start new band accumulation */
				in_gap = 0;
				tot_prob = curr_prob;
				esum = e*curr_prob;
				e2sum = e2*curr_prob;
			}  /* else if (e-e_prev > min_gap) */
		}  /* if (!in_gap) */
		else
		{
			/* We are in a gap; out of gap on large sample only */
			if (curr_prob > tiny_prob)
			{
				in_gap = 0;
				try_gap = 0;
				tot_prob = curr_prob;
				esum = e*curr_prob;
				e2sum = e2*curr_prob;
			}  /* if (curr_prob > tiny_prob) */
		}  /* else */
		
		/* Update previous energy */
		e_prev = e;
	}  /* for i */	

	/* Now record last band -- whatever is left over if we were NOT in gap */
	if (!in_gap)
	{
		j = nb[kix];
		temp = esum/tot_prob;
		begymat[j][kix] = temp;
		temp2 = e2sum/tot_prob;
		bsdevmat[j][kix] = sqrt(temp2 - temp*temp);
		bprobmat[j][kix] = tot_prob;
		nb[kix] = j+1;
	}  /* if (!in_gap) */
}  /* band_guess */

	
				
	

int fccguess(void)
{
	int i, j, Nc, same, negy, maxbands, kix, *nb;
	double **Ekmat, **begymat, **bsdevmat, **bprobmat, *energies, **Gv;
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
	Gv = R2tensor(0,Nc-1,0,2);//allocate_real_matrix(0,Nc-1,0,2);
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
			fread(&(Gv[i][j]),sizeof(double),1,fp);
	
	/* Now read Ekmatrix */
	for(i=0;i<TR*negy;i++)
		for(j=0;j<Nc;j++)
			fread(&(Ekmat[i][j]),sizeof(double),1,fp);
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
				kix,Gv[kix][0]/M_PI,Gv[kix][1]/M_PI,Gv[kix][2]/M_PI,nb[kix]);

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
		fprintf(fp4,"k_unit = (%12.5le, %12.5le, %12.5le)\n\n",k_unit[0],k_unit[1],k_unit[2]);

		for(kix=0;kix<Nc;kix++)
		{
			/* Copy a vector, get mag, and test */
			kmag = 0.0;
			for(j=0;j<3;j++)
			{
				kv[j] = Gv[kix][j];
				kmag += kv[j]*kv[j];
			}  /* for j */

			/* We print out magnitude of k normalized to PI */
			kmag = sqrt(kmag)/M_PI;

			/* Always count (0,0,0) as a vaild state */
			if (kmag < EPS_K)
				same = 1;
			else
				same = comp_dirFCC(kv,k_unit,EPS_K);

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
					fprintf(fp4,"  %12.5le  %12.5le  %12.5le  %12.5le\n",kmag,begymat[j][kix],bsdevmat[j][kix],bprobmat[j][kix]);
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
	rm_r2tensor(&Gv,0);//free_real_matrix(Gv,0,Nc-1,0);
	rm_rvectr(&energies);//free_real_vector(energies,0);
	rm_r2tensor(&begymat,0);//free_real_matrix(begymat,0,maxbands-1,0);
	rm_r2tensor(&bprobmat,0);//free_real_matrix(bprobmat,0,maxbands-1,0);
	rm_r2tensor(&bsdevmat,0);//free_real_matrix(bsdevmat,0,maxbands-1,0);
	rm_r2tensor(&Ekmat,0);//free_real_matrix(Ekmat,0,negy-1,0);
	rm_ivectr(&nb);//free_integer_vector(nb,0);

	return(0);
}  /* fccguess */
	
