/*****************************************************************************
The JPL Parallel Eigensolvers package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This file is part of
The Jet Propulsion Laboratory (JPL) Parallel Eigensolvers package.

This library is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty
of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)
  E. Robert Tisdale (E.Robert.Tisdale@jpl.nasa.gov)

Written by: Chris Bowen
            Gerhard Klimeck
            Fabiano Oyafuso
            E. Robert Tisdale

*****************************************************************************
$Source: /repo/eigen/src/lanczos_blk_old.c,v $
*****************************************************************************/
#include "lanczos_blk.h"
#ifdef INTEL_MODIFIED
	#include <memory.h>
	#include <emmintrin.h>
#endif //INTEL_MODIFIED


void single_blk_lanczos_iteration(
    void (*mul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const
    int* argument[],	/* in    matrix-vector multipier argument list	*/
    nml_dcscalar *r1_n,	/* inout r1_{n-1} = q1_n*beta_{n-1}		*/
    nml_dcscalar *r2_n,	/* inout r2_{n-1} = q2_n*beta_{n-1}		*/
    nml_dcscalar *q1_n,	/*   out current  complex Lanczos vector q1_{n}	*/
    nml_dcscalar *q2_n,	/*   out current  complex Lanczos vector q2_{n}	*/
    nml_dcscalar *q3_n,	/*   out current  complex Lanczos vector q1_{n}	*/
    nml_dcscalar *q4_n,	/*   out current  complex Lanczos vector q2_{n}	*/
    nml_extent	 m,	/* in    extent of vectors r_n, q_n and q_{n-1}	*/
    nml_extent	 iteration,/* in    the current iteration 	*/
    nml_dscalar	 *alpha,/*   out real  on-diagonal element alpha_{n}	*/
    nml_dscalar	 *alpha2,/*   out real off-diagonal element beta_{n}	*/
    nml_dcscalar *beta1,/* in    real off-diagonal element beta_{n-1}	*/
    nml_dcscalar *beta2,/* in    real off-diagonal element beta_{n-1}	*/
    nml_dscalar	 *gamma1,/*   out real off-diagonal element gamma_{n}	*/
    nml_dscalar	 *gamma2,/*   out real off-diagonal element gamma_{n}	*/
    nml_dscalar	 *ga1,/*   inout norm of r1_n 	*/
    nml_dscalar	 *ga2/*   inout norm of r2_n	*/
    ) {			/* Lanczos algorithm (iteration n)		*/
  /* q_{0} = 0, r_{0} = arbitrary and beta_{0} = ||r_{0}||.
	for n = 1, 2, 3, . . .
	      q_{n} <-- r_{n-1}/beta_{n-1}
	      r_{n} <-- A*q_{n}		// matrix-vector product
	  alpha_{n} <-- <q_{n} | r_{n}>	// vector-vector product
	      r_{n} <-- r_{n} - beta_{n-1}*q_{n-1} - alpha_{n}*q_{n}
	   beta_{n} <-- ||r_{n}||
  */
  nml_offset	j = 0;
  nml_dcscalar	alc = {0.0, 0.0};
  if(iteration==0){
  nml_dscalar g1_n1=nml_dcv_ssqr(r1_n);//g1=norm(r1) 
  nml_dscalar g1=0.0;
  Allreduce_MPI_sp(&g1_n1, &g1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  g1 = sqrt(g1);
  *ga1=g1;
  nml_dscalar g2_n1=nml_dcv_ssqr(r2_n);//g2=norm(r2) 
  nml_dscalar g2=0.0;
  Allreduce_MPI_sp(&g2_n1, &g2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  g2 = sqrt(g2);
  *ga2=g2;
  for (j = 0; j < m; ++j) {/* q1_{n} <-- r1_{n-1}/norm of r1_{n-1}	*/
    nml_dscalar	inv_g_n1 = 1.0/g1;
    /* first do the real part ... */
    nml_dscalar	*q1_n_j = &(q1_n[j].r);
    nml_dscalar	*r1_n_j = &(r1_n[j].r);
    *q1_n_j = *r1_n_j*inv_g_n1;
    /* now do the imaginary part */
    ++q1_n_j;
    ++r1_n_j;
    *q1_n_j = *r1_n_j*inv_g_n1;
    }

  //b2=r2'*q1 b2=q1'*r2
  for (j = 0; j < m; ++j) {
    /* alc = cadd(alc, cmul(cconj(q_n[j]), r_n[j]));			*/
    nml_dscalar *r2_n_j = &(r2_n[j].r);
    nml_dscalar *q1_n_j = &(q1_n[j].r);

    alc.r += (r2_n_j[0])*(q1_n_j[0]);
    alc.i += r2_n_j[1]*(q1_n_j[0]);
    alc.i -= (r2_n_j[0])*(q1_n_j[1]);
    alc.r += (r2_n_j[1])*(q1_n_j[1]);
    }
    nml_dscalar b2_n1[2] ;
    b2_n1[0]=alc.r;
    b2_n1[1]=alc.i;
    nml_dcscalar b2={0.0,0.0};
    { nml_dscalar       bl[2];
      bl[0]=0.0;
      bl[1]=0.0;
      Allreduce_MPI_sp(&b2_n1[0], &bl[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      b2.r = bl[0];
      b2.i = bl[1];
      }
    //printf("b2_n1[0]= %1.15e, [1]=%1.15e\n",b2_n1[0],b2_n1[1]);
    //printf("b2[0]= %1.15e, [1]=%1.15e\n",b2.r,b2.i);

  for (j = 0; j < m; ++j) { //q2=r2-b2*q1
    nml_dscalar *r2_n_j = &(r2_n[j].r);
    nml_dscalar *q1_n_j = &(q1_n[j].r);
    nml_dscalar *q2_n_j = &(q2_n[j].r);
    *q2_n_j=(*r2_n_j)-b2.r*(*q1_n_j)+b2.i*(q1_n_j[1]);
    ++q2_n_j;
    ++r2_n_j;
    *q2_n_j=(*r2_n_j)-b2.r*(q1_n_j[1])-b2.i*(*q1_n_j);
    }
    
     g2_n1=nml_dcv_ssqr(q2_n);//g2=norm(q2) 
     g2 = 0.0;
    {  nml_dscalar       g2_2 = 0.0;
      Allreduce_MPI_sp(&g2_n1, &g2_2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      g2 = sqrt(g2_2);
      }
  //printf("g2_n1 %1.15e,g2 %1.15e\n",g2_n1,g2);
  for (j = 0; j < m; ++j) {	/* q2_{n} <-- q2_{n-1}/g2		*/
    /*
    q_n[j].r = r_n[j].r/beta_n1;
    q_n[j].i = r_n[j].i/beta_n1;
    */

    nml_dscalar	inv_beta_n2 = 1.0/g2;
    /* first do the real part ... */
    nml_dscalar	*q2_n_j = &(q2_n[j].r);
    *q2_n_j = *q2_n_j*inv_beta_n2;
    /* now do the imaginary part */
    ++q2_n_j;
    *q2_n_j = *q2_n_j*inv_beta_n2;
    }
  //r1=A*q1
  mul(argument, r1_n, q1_n);	/* do matrix-vector multiply		*/
  //r2=A*q2
  mul(argument, r2_n, q2_n);	/* do matrix-vector multiply		*/
  //a1=q1'*r1
  alc.r=0.0; alc.i=0.0;
    for (j = 0; j < m; ++j) {
    nml_dscalar *r1_n_j = &(r1_n[j].r);
    nml_dscalar *q1_n_j = &(q1_n[j].r);

    alc.r += (*r1_n_j)*(*q1_n_j);
    alc.i += r1_n_j[1]*(*q1_n_j);
    ++q1_n_j;
    alc.i -= (*r1_n_j)*(*q1_n_j);
    ++r1_n_j;
    alc.r += (*r1_n_j)*(*q1_n_j);
       }
  nml_dscalar a1_n1=alc.r;
    { nml_dscalar       al1 = 0.0;
      Allreduce_MPI_sp(&a1_n1, &al1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      *alpha = al1;
      }
  //alpha2=q2'*r2
  alc.r=0.0; alc.i=0.0;
    for (j = 0; j < m; ++j) {
    nml_dscalar *r2_n_j = &(r2_n[j].r);
    nml_dscalar *q2_n_j = &(q2_n[j].r);

    alc.r += (*r2_n_j)*(*q2_n_j);
    alc.i += r2_n_j[1]*(*q2_n_j);
    ++q2_n_j;
    alc.i -= (*r2_n_j)*(*q2_n_j);
    ++r2_n_j;
    alc.r += (*r2_n_j)*(*q2_n_j);
       }
  nml_dscalar a2_n1=alc.r;
    { nml_dscalar       al2 = 0.0;
      Allreduce_MPI_sp(&a2_n1, &al2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      *alpha2 = al2;
      }
  
  //beta1=q2'*r1
  alc.r=0.0; alc.i=0.0;
    for (j = 0; j < m; ++j) {
    nml_dscalar *r1_n_j = &(r1_n[j].r);
    nml_dscalar *q2_n_j = &(q2_n[j].r);

    alc.r += (*r1_n_j)*(*q2_n_j);
    alc.i += r1_n_j[1]*(*q2_n_j);
    ++q2_n_j;
    alc.i -= (*r1_n_j)*(*q2_n_j);
    ++r1_n_j;
    alc.r += (*r1_n_j)*(*q2_n_j);
       }
  nml_dscalar b1_n1[2];
    b1_n1[0]=alc.r;
    b1_n1[1]=alc.i;
    { nml_dscalar       bl1[2];
	    bl1[0] = 0.0;
	    bl1[1] = 0.0;
      Allreduce_MPI_sp(&b1_n1[0], &bl1[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      (*beta2).r = bl1[0];
      (*beta2).i = bl1[1];
      }
  
    for (j = 0; j < m; ++j) {
      q3_n[j].r = 0.0 ;
      q3_n[j].i = 0.0;
      }
    for (j = 0; j < m; ++j) {
      q4_n[j].r = 0.0; 
      q4_n[j].i = 0.0;
      }
  *gamma1=0.0;
  *gamma2=0.0;
  (*beta1).r=0.0; 
  (*beta1).i=0.0; 
    //printf("alpha %1.15e, beta2 %1.15e %1.15e alphas2 %1.15e\n",*alpha,(*beta2).r,(*beta2).i,*alpha2);
  }
  else{

#ifdef INTEL_MODIFIED
#include "lanczos_blk_SSE.c"
#else
      nml_dscalar a1=*alpha;
      nml_dcscalar b2;
      b2.r=(*beta1).r;
      b2.i=(*beta1).i;
      nml_dscalar a2=*alpha2;
      nml_dcscalar b1;
      b1.r=(*beta2).r;
      b1.i=(*beta2).i;
      nml_dscalar g1=*gamma1;
      nml_dscalar g2=*gamma2;
      
    for (j = 0; j < m; ++j) {
      r1_n[j].r -= a1*q1_n[j].r + b1.r*q2_n[j].r - b1.i*q2_n[j].i + g1*q3_n[j].r + b2.r*q4_n[j].r + b2.i*q4_n[j].i;
      r1_n[j].i -= a1*q1_n[j].i + b1.r*q2_n[j].i + b1.i*q2_n[j].r + g1*q3_n[j].i + b2.r*q4_n[j].i - b2.i*q4_n[j].r;
      }
    for (j = 0; j < m; ++j) {
      r2_n[j].r -= b1.r*q1_n[j].r+b1.i*q1_n[j].i + a2*q2_n[j].r+g2*q4_n[j].r;
      r2_n[j].i -= b1.r*q1_n[j].i-b1.i*q1_n[j].r + a2*q2_n[j].i+g2*q4_n[j].i;
      }

	for (j = 0; j < m; ++j) {
		q3_n[j].r = q1_n[j].r ;
		q3_n[j].i = q1_n[j].i;
	}
    for (j = 0; j < m; ++j) {
		q4_n[j].r = q2_n[j].r ;
		q4_n[j].i = q2_n[j].i;
	}
	
	
#endif       

  nml_dscalar g1_n1=nml_dcv_ssqr(r1_n);//g1=norm(r1)
  nml_dscalar g1_l=0.0;
  Allreduce_MPI_sp(&g1_n1, &g1_l, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *gamma1 = sqrt(g1_l);
  *ga1=*gamma1;
  nml_dscalar g2_n2=nml_dcv_ssqr(r2_n);//g1=norm(r1)
  nml_dscalar g2_2=0.0;
  Allreduce_MPI_sp(&g2_n2, &g2_2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
  *ga2 = sqrt(g2_2);
  for (j = 0; j < m; ++j) {     /* q1_{n} <-- r1_{n-1}/norm of r1_{n-1}         */
    nml_dscalar inv_g_n1 = 1.0/(*gamma1);
    /* first do the real part ... */
    nml_dscalar *q1_n_j = &(q1_n[j].r);
    nml_dscalar *r1_n_j = &(r1_n[j].r);
    *q1_n_j = *r1_n_j*inv_g_n1;
    /* now do the imaginary part */
    ++q1_n_j;
    ++r1_n_j;
    *q1_n_j = *r1_n_j*inv_g_n1;
    }

  //b2=r2'*q1 b2=q1'*r2
  for (j = 0; j < m; ++j) {
    /* alc = cadd(alc, cmul(cconj(q_n[j]), r_n[j]));                    */
    nml_dscalar *r2_n_j = &(r2_n[j].r);
    nml_dscalar *q1_n_j = &(q1_n[j].r);

    alc.r += (*r2_n_j)*(*q1_n_j);
    alc.i += r2_n_j[1]*(*q1_n_j);
    ++q1_n_j;
    alc.i -= (*r2_n_j)*(*q1_n_j);
    ++r2_n_j;
    alc.r += (*r2_n_j)*(*q1_n_j);
    }
    nml_dscalar b2_n1[2];
    b2_n1[0] = alc.r;
    b2_n1[1] = alc.i;
    { nml_dscalar       bl_l[2];
	   bl_l[0] = 0.0;
	   bl_l[1] = 0.0;
      Allreduce_MPI_sp(&b2_n1[0], &bl_l[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      (*beta1).r = bl_l[0];
      (*beta1).i = bl_l[1];
      }
  for (j = 0; j < m; ++j) { //q2=r2-b2*q1
    nml_dscalar *r2_n_j = &(r2_n[j].r);
    nml_dscalar *q1_n_j = &(q1_n[j].r);
    nml_dscalar *q2_n_j = &(q2_n[j].r);
    *q2_n_j=(*r2_n_j)-(*beta1).r*(*q1_n_j)+(*beta1).i*(q1_n_j[1]);
    //if(j==0){ printf("r2 = %1.15e q1= %1.15e\n",*r2_n_j,*q1_n_j);}
    ++q2_n_j;
    ++r2_n_j;
    *q2_n_j=(*r2_n_j)-(*beta1).i*(*q1_n_j)-(*beta1).r*(q1_n_j[1]);
    //if(j==0){ printf("r2 = %1.15e q1= %1.15e\n",*r2_n_j,*q1_n_j);}
    }
  nml_dscalar g2_n1=nml_dcv_ssqr(q2_n);//g2=norm(q2)
    { nml_dscalar       gl = 0.0;
      Allreduce_MPI_sp(&g2_n1, &gl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      *gamma2 = sqrt(gl);
      }
  //printf("g2 %1.15e\n",g2);
  for (j = 0; j < m; ++j) {     /* q2_{n} <-- q2_{n-1}/g2               */
    /*
    q_n[j].r = r_n[j].r/beta_n1;
    q_n[j].i = r_n[j].i/beta_n1;
    */

    nml_dscalar inv_beta_n2 = 1.0/(*gamma2);
    /* first do the real part ... */
    nml_dscalar *q2_n_j = &(q2_n[j].r);
    *q2_n_j = *q2_n_j*inv_beta_n2;
    /* now do the imaginary part */
    ++q2_n_j;
    *q2_n_j = *q2_n_j*inv_beta_n2;
    }
  //r1=A*q1
  mul(argument, r1_n, q1_n);    /* do matrix-vector multiply            */
  //r2=A*q2
  mul(argument, r2_n, q2_n);    /* do matrix-vector multiply            */
  //a1=q1'*r1
  alc.r=0.0; alc.i=0.0;
    for (j = 0; j < m; ++j) {
    nml_dscalar *r1_n_j = &(r1_n[j].r);
    nml_dscalar *q1_n_j = &(q1_n[j].r);

    alc.r += (*r1_n_j)*(*q1_n_j);
    alc.i += r1_n_j[1]*(*q1_n_j);
    ++q1_n_j;
    alc.i -= (*r1_n_j)*(*q1_n_j);
    ++r1_n_j;
    alc.r += (*r1_n_j)*(*q1_n_j);
       }
  nml_dscalar a1_n1=alc.r;
    { nml_dscalar       al1 = 0.0;
      Allreduce_MPI_sp(&a1_n1, &al1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      *alpha = al1;
      }
  //alpha2=q2'*r2
  alc.r=0.0; alc.i=0.0;
    for (j = 0; j < m; ++j) {
    nml_dscalar *r2_n_j = &(r2_n[j].r);
    nml_dscalar *q2_n_j = &(q2_n[j].r);

    alc.r += (*r2_n_j)*(*q2_n_j);
    alc.i += r2_n_j[1]*(*q2_n_j);
    ++q2_n_j;
    alc.i -= (*r2_n_j)*(*q2_n_j);
    ++r2_n_j;
    alc.r += (*r2_n_j)*(*q2_n_j);
       }
  nml_dscalar a2_n1=alc.r;
    { nml_dscalar       al2 = 0.0;
      Allreduce_MPI_sp(&a2_n1, &al2, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      *alpha2 = al2;
      }

  //beta2=q2'*r1
  alc.r=0.0; alc.i=0.0;
    for (j = 0; j < m; ++j) {
    nml_dscalar *r1_n_j = &(r1_n[j].r);
    nml_dscalar *q2_n_j = &(q2_n[j].r);

    alc.r += (*r1_n_j)*(*q2_n_j);
    alc.i += r1_n_j[1]*(*q2_n_j);
    ++q2_n_j;
    alc.i -= (*r1_n_j)*(*q2_n_j);
    ++r1_n_j;
    alc.r += (*r1_n_j)*(*q2_n_j);
       }
  nml_dscalar b1_n1[2];
    b1_n1[0]=alc.r;
    b1_n1[1]=alc.i;
    { nml_dscalar       bl1[2];
	   bl1[0] = 0.0;
	   bl1[1] = 0.0;
      Allreduce_MPI_sp(&b1_n1[0], &bl1[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
      (*beta2).r = bl1[0];
      (*beta2).i = bl1[1];
      }

  }
}/*block lanczos*/
static int convMethodValue=2;
void find_block_converged_eigval(
    int		iterations,	/* in Lanczos iteration			*/
    nml_dvector *ALPHA,
    nml_dcvector *BETA,
    nml_dvector *GAMMA,
    nml_dvector	*eigv_conv,	/* inout actual converged eigenvalues	*/
    int		*pEigenvalues,	/*   out number converged eigenvalues	*/
    int		*n_nosuccess,	/*   out number no success		*/
    int		*n_old,		/* inout number old eigenvalues		*/
    nml_dvector	*eigv_old,	/* inout actual old eigenvalues		*/
    FILE	*fp_trace,	/* in trace   log text file pointer	*/
    int		verbose,	/* in verbose diagnostic messages	*/
    nml_dscalar	emin,		/* in eigenvalue search lower bound	*/
    nml_dscalar	emax,		/* in eigenvalue search upper bound	*/
    nml_dscalar	tolerance,	/* in eigenvalue convergence tolerance	*/
    nml_dscalar	resolution,	/* in eigenvalue separation minimum	*/
    FILE	*fp_log,	/* in message log text file pointer	*/
    nml_dscalar	ga1,	/* in norm of r1_n	*/
    nml_dscalar	ga2	/* in norm of r2_n	*/
    ) {
  const nml_dscalar	bias = 1e-80;	/* PLEASE EXPLAIN ME			*/
  const nml_dscalar	gain = 1e-8;	/* PLEASE EXPLAIN ME			*/
  nml_dscalar	deigv_min = DBL_MAX; /*double	deigv_min = 1.0e100;	*/
  nml_dscalar	 eigv_min = DBL_MAX; /*double	 eigv_min = 1.0e100;	*/
  int		Neval_in_range = 0;	/* number eigenvalues found	*/
  nml_dvector	*eval = nml_dv_new(2*iterations);
					/* actual eigenvalues found	*/
  *n_nosuccess = 0;
   int i=0;
   Neval_in_range=*pEigenvalues;  // *pEigenvalues inputs the number of eigenvalues 
                                  //solved by solve_pentamatrix_eigenvalue function
				  //and outputs the number of eigenvalues that has no
				  //ghost eigenvalues
				  
   for (i=0;i<Neval_in_range;++i) {
     eval[i]=eigv_conv[i];
   }
  *pEigenvalues = 0;

  /* convergence Method for eigenvalue
   * there are two Methods are implemented:
   *      1. check if the eigenvalue converge or not (old)
   *      2. check if the final component of the eigenvector
   *         is less than tolerance (new)
   */
   /*Method 1*/
   if(convMethodValue==1)
     { 
       for (i = 0; i < Neval_in_range; ++i) {
         if ((0 < i) && (fabs(eval[i-1] - eval[i]) < (bias + eval[i])*gain))
           continue;               // ignore other degenerate eigenvalues
         if (NULL != fp_trace) {
           fprintf(fp_trace, "%d %1.15e\n", iterations, eval[i]);
           fflush(fp_trace);
           }
       { int convskip = 0;
         int j = 0;
         for (j = 0; j < *pEigenvalues; ++j)
           if (fabs(eval[i] - eigv_conv[j]) < resolution){
             convskip = 1;
             //printf("%1.15e skipped\n",eval[i]);
	   }
	   if (0 == convskip) {
	     int j = 0;
	     for (j = 0; j < *n_old; ++j) {
	       nml_dscalar deigv = fabs(eigv_old[j] - eval[i]);
	       if (deigv < deigv_min) {
	         deigv_min = deigv;
	         eigv_min = eval[i];
	       }
             if (deigv < tolerance) {//check if the eigenvalue converges
               eigv_conv[*pEigenvalues] = eval[i];
               if (verbose && (NULL != fp_log)) {
                 fprintf(fp_log, "   eigenvalue candidate #%d:  %18.15f\n",
                 *pEigenvalues, eigv_conv[*pEigenvalues]*MTRX_SCALE);
                 fflush(fp_log);
               }
               ++(*pEigenvalues);
               *n_nosuccess = 0;
             }
           }
         }
       }
       }
     }
  /*Method 2*/
  else if(convMethodValue==2)
  {
    nml_dcvector        *xx = nml_dcv_new(2*iterations);
    nml_dcvector	*xx_prev=nml_dcv_new(2*iterations);
    nml_dcvector        *yy = nml_dcv_new(2*iterations);
    int      iteration = 0;
    int k;
    int flag;
     
    for (i = 0;i < Neval_in_range; ++i) {
      flag=1;
      for (iteration = 0; iteration < 2*iterations; ++iteration) {
        yy[iteration].r = 1.0*rand()/RAND_MAX;
        yy[iteration].i = 0.0;
      }
      solve_pentamatrix_vector(iterations, ALPHA,BETA,GAMMA,xx,yy,eval[i]);
      		

#ifdef DEBUGCONVERGE
      for(iteration=2*iterations-20;iteration<2*iterations;iteration++){
        printf("xx[%d].r=%1.15e xx[%d].i=%1.15e\n",iteration,xx[iteration].r,iteration,xx[iteration].i);
      }
      printf("ga1 = %1.15e, ga2 = %1.15e\n",ga1,ga2);
      printf("current %d eigenvalue = %1.15e tol=%1.15e\n",2*iterations-2,eval[i],tolerance);
#endif

      nml_dscalar mod=xx[2*iterations-2].r*xx[2*iterations-2].r+xx[2*iterations-2].i*xx[2*iterations-2].i;
      mod=sqrt(mod);
      nml_dscalar mod2=xx[2*iterations-1].r*xx[2*iterations-1].r+xx[2*iterations-1].i*xx[2*iterations-1].i;
      mod2=sqrt(mod2);
#ifdef DEBUGCONVERGE
      printf("%1.15e %1.15e ga1*mod+ga2*mod2 = %1.15e\n",mod, mod2,ga1*mod+ga2*mod2);
#endif
      if((ga1*mod+ga2*mod2)>10*tolerance/*&&fabs(xx[iteration].i)>tolerance*/)
          flag=0;



      //use Gram-schmidt algorithm to orthorgonalize two eigenvectors of a duplicate eigenvalue, and check if the resulted eigenvectors are orthogonal
      if(flag==1&&i>1&&(*pEigenvalues)>1&&fabs(eval[i-1]-eval[i])<1e-14&&fabs(eval[i-1]-eigv_conv[(*pEigenvalues)-1])<1e-14)
      {
        nml_dscalar invert=0.0;
	invert=nml_dcv_ssqr(xx);
	invert=1.0/sqrt(invert);
	for(k=0;k<2*iterations;k++){
	  xx_prev[k].r=xx_prev[k].r*invert;
	  xx_prev[k].i=xx_prev[k].i*invert;
	}
	invert=nml_dcv_ssqr(xx);
	invert=1.0/sqrt(invert);
	for(k=0;k<2*iterations;k++){
	  xx[k].r=xx[k].r*invert;
	  xx[k].i=xx[k].i*invert;
	}
        nml_dcscalar dotproduct={0.0, 0.0};
	//  dotproduct=xx_prev'*xx;
        for(k=0;k<2*iterations;k++)
        {
          nml_dscalar *xx_prev_n_j = &(xx_prev[k].r);
	  nml_dscalar *xx_n_j = &(xx[k].r);

          dotproduct.r += (*xx_prev_n_j)*(*xx_n_j);
          dotproduct.i -= xx_prev_n_j[1]*(*xx_n_j);
	  ++xx_n_j;
          dotproduct.i += (*xx_prev_n_j)*(*xx_n_j);
	  ++xx_prev_n_j;
          dotproduct.r += (*xx_prev_n_j)*(*xx_n_j);
	}
	if(fabs(dotproduct.r)<1e-14&&fabs(dotproduct.i)<1e-14){
	  continue;
	}
	//xx=xx-dotproduct*xx_prev=xx-(xx_prev'*xx)*xx_prev
        for(k=0;k<2*iterations;k++)
        {
	  xx[k].r=xx[k].r-dotproduct.r*xx_prev[k].r+dotproduct.i*xx_prev[k].i;
	  xx[k].i=xx[k].i-dotproduct.r*xx_prev[k].i-dotproduct.i*xx_prev[k].r;
	}
	invert=nml_dcv_ssqr(xx);
	invert=sqrt(invert);
	invert=1/invert;
	for(k=0;k<2*iterations;k++){
	  xx[k].r=xx[k].r*invert;
	  xx[k].i=xx[k].i*invert;
	}
	dotproduct.r=0.0;
	dotproduct.i=0.0;
	//dotproduct= xx_prev'*xx check if two eigenvectos after Gram-schmdit are orthorgonal
	for(k=0;k<2*iterations;k++)
	{
          nml_dscalar *xx_prev_n_j = &(xx_prev[k].r);
	  nml_dscalar *xx_n_j = &(xx[k].r);

          dotproduct.r += (*xx_prev_n_j)*(*xx_n_j);
          dotproduct.i -= xx_prev_n_j[1]*(*xx_n_j);
	  ++xx_n_j;
          dotproduct.i += (*xx_prev_n_j)*(*xx_n_j);
	  ++xx_prev_n_j;
          dotproduct.r += (*xx_prev_n_j)*(*xx_n_j);
	}
	if(dotproduct.r>1e-14||dotproduct.i>1e-14)
	  flag=0;
      }
      if(flag==1){
        int j = 0;
        for (j = 0; j < *n_old; ++j) {
          nml_dscalar deigv = fabs(eigv_old[j] - eval[i]);
          if (deigv < deigv_min) {
            deigv_min = deigv;
            eigv_min = eval[i];
          }
        }
        for(k=0;k<2*iterations;k++)
        {
	   xx_prev[k].r=xx[k].r;
	   xx_prev[k].i=xx[k].i;
        }
        {
        eigv_conv[*pEigenvalues]=eval[i];
        if (verbose && (NULL != fp_log)) {
          fprintf(fp_log, "   eigenvalue candidate #%d:  %18.15f\n",
          *pEigenvalues, eigv_conv[*pEigenvalues]*MTRX_SCALE);
          fflush(fp_log);
        }
        ++(*pEigenvalues);
        }
      }
//#define Timing_detail
#ifdef Timing_detail
      fprintf(fp_log,"%d = %1.15f\n",i,eval[i]);
#endif
    }
    nml_dcv_delete(xx);
    nml_dcv_delete(xx_prev);
    nml_dcv_delete(yy);
  }

  if (verbose) {
    if (deigv_min < DBL_MAX)
      if (NULL != fp_log)
	fprintf(fp_log, "    E = %e  dE_min = %e\n",
	  eigv_min*MTRX_SCALE, deigv_min*MTRX_SCALE);
  }

  *n_old = 0;
  { int	iteration = 0;
    for (iteration = 0; iteration < iterations; ++iteration)
      if ((emin/MTRX_SCALE < eval[iteration])
       && (eval[iteration] < emax/MTRX_SCALE)
       && ((0 == iteration) || (resolution < fabs(
	(eval[iteration-1] - eval[iteration])/(bias + eval[iteration]))))) {
	eigv_old[*n_old] = eval[iteration];
	++(*n_old);
        }
    }
    
  ++(*n_nosuccess);

  nml_dv_delete(eval);
  } /* find_block_converged_eigval */

int gramSchmidt(
     nml_dcvector* xx_prev,
     nml_dcvector* xx,
     nml_extent length
     ) {
      int k;
      //use Gram-schmidt algorithm to orthorgonalize two eigenvectors of a duplicate eigenvalue, and check if the resulted eigenvectors are orthogonal
      {
	nml_dscalar invert=0.0;
        for (k = 0; k < length; ++k)
          invert += xx[k].r*xx[k].r+xx[k].i*xx[k].i;
	nml_dscalar invertsum=0.0;
        Allreduce_MPI_sp(&invert, &invertsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        invertsum=sqrt(invertsum);
        invert=1/invertsum;
	for(k=0;k<length;k++){
	  xx[k].r=xx[k].r*invert;
          xx[k].i=xx[k].i*invert;
	}
	invert=0.0;
        for (k = 0; k < length; ++k)
          invert += xx[k].r*xx[k].r+xx[k].i*xx[k].i;
	invertsum=0.0;
        Allreduce_MPI_sp(&invert, &invertsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        invertsum=sqrt(invertsum);
        invert=1/invertsum;
	for(k=0;k<length;k++){
	  xx_prev[k].r=xx_prev[k].r*invert;
          xx_prev[k].i=xx_prev[k].i*invert;
	}
        nml_dcscalar dotproduct={0.0, 0.0};
        //  dotproduct=xx_prev'*xx;
        for(k=0;k<length;k++)
        {
          nml_dscalar *xx_prev_n_j = &(xx_prev[k].r);
          nml_dscalar *xx_n_j = &(xx[k].r);

          dotproduct.r += (*xx_prev_n_j)*(*xx_n_j);
          dotproduct.i -= xx_prev_n_j[1]*(*xx_n_j);
          ++xx_n_j;
          dotproduct.i += (*xx_prev_n_j)*(*xx_n_j);
          ++xx_prev_n_j;
          dotproduct.r += (*xx_prev_n_j)*(*xx_n_j);
        }
        nml_dscalar dotbeforesum[2];
        dotbeforesum[0]=dotproduct.r;
        dotbeforesum[1]=dotproduct.i;
	nml_dscalar dotsum[2];
	dotsum[0]=0.0;
	dotsum[1]=0.0;
	Allreduce_MPI_sp(&dotbeforesum[0], &dotsum[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
	dotproduct.r=dotsum[0];
	dotproduct.i=dotsum[1];
	
        if(fabs(dotproduct.r)<1e-14&&fabs(dotproduct.i)<1e-14){
          return 1;
        }
        //project xx to xx_prev , xx=xx-dotproduct*xx_prev=xx-(xx_prev'*xx)*xx_prev
        for(k=0;k<length;k++)
        {
          xx[k].r=xx[k].r-dotproduct.r*xx_prev[k].r+dotproduct.i*xx_prev[k].i;
          xx[k].i=xx[k].i-dotproduct.r*xx_prev[k].i-dotproduct.i*xx_prev[k].r;
        }
	invert=0.0;
        for (k = 0; k < length; ++k)
          invert += xx[k].r*xx[k].r+xx[k].i*xx[k].i;
	invertsum=0.0;
        Allreduce_MPI_sp(&invert, &invertsum, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        invertsum=sqrt(invertsum);
        invert=1/invertsum;
        //normalize the projected vector xx
        for(k=0;k<length;k++)
        {
          xx[k].r=xx[k].r*invert;
          xx[k].i=xx[k].i*invert;
        }
        dotproduct.r=0.0;
        dotproduct.i=0.0;
        //dotproduct= xx_prev'*xx check if two eigenvectors after Gram-schmdit are orthorgonal
        for(k=0;k<length;k++)
        {
          nml_dscalar *xx_prev_n_j = &(xx_prev[k].r);
          nml_dscalar *xx_n_j = &(xx[k].r);

          dotproduct.r += (*xx_prev_n_j)*(*xx_n_j);
          dotproduct.i -= xx_prev_n_j[1]*(*xx_n_j);
          ++xx_n_j;
          dotproduct.i += (*xx_prev_n_j)*(*xx_n_j);
          ++xx_prev_n_j;
          dotproduct.r += (*xx_prev_n_j)*(*xx_n_j);
        }
        dotbeforesum[0]=dotproduct.r;
        dotbeforesum[1]=dotproduct.i;
        dotsum[0]=0.0;
        dotsum[1]=0.0;
        Allreduce_MPI_sp(&dotbeforesum[0], &dotsum[0], 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        dotproduct.r=dotsum[0];
        dotproduct.i=dotsum[1];

        if(dotproduct.r<1e-14||dotproduct.i<1e-14)
	  return 1;
	else
	  return 0;
      }
    }/* gramSchmidt */

int penta_LDLT_f(
  int iteration,
  nml_dvector *aa,
  nml_dcscalar *b,
  nml_dvector *g,
  nml_dvector *D,
  nml_dcscalar *p,
  nml_dvector *q,
  nml_dscalar shift
  ) {
   int N=2*iteration;
   int j=0;
   D[0]=aa[0]-shift;
   p[0].r=b[1].r/D[0];
   p[0].i=b[1].i/D[0];
   D[1]=aa[1]-shift-(p[0].r*p[0].r+p[0].i*p[0].i)*D[0];
   for (j=2;j<N;j++){
     q[j-2]=g[j]/D[j-2];
     p[j-1].r=(b[j].r-q[j-2]*p[j-2].r*D[j-2])/D[j-1];
     p[j-1].i=(b[j].i+q[j-2]*p[j-2].i*D[j-2])/D[j-1];
     D[j]=aa[j]-shift-p[j-1].r*p[j-1].r*D[j-1]-p[j-1].i*p[j-1].i*D[j-1]
          -q[j-2]*q[j-2]*D[j-2];
	     
   }
   return 0;
   }/* penta_LDLT_f */
int penta_LDLT_s(
  int iteration,
  nml_dvector *D,
  nml_dcvector *p,
  nml_dvector *q,
  nml_dcvector *x,
  nml_dcvector *rhs
  ) {
    int N=2*iteration;
    int i=0;
    nml_dcvector *y =  nml_dcv_new(N);
    y[0].r=rhs[0].r;
    y[0].i=rhs[0].i;
    y[1].r=rhs[1].r-p[0].r*y[0].r+p[0].i*y[0].i; 
    y[1].i=rhs[1].i-p[0].i*y[0].r-p[0].r*y[0].i;
    for(i=2;i<N;i++) {
       y[i].r=rhs[i].r-p[i-1].r*y[i-1].r+p[i-1].i*y[i-1].i-q[i-2]*y[i-2].r;
       y[i].i=rhs[i].i-p[i-1].i*y[i-1].r-p[i-1].r*y[i-1].i-q[i-2]*y[i-2].i;
    } 
    x[N-1].r=y[N-1].r/D[N-1];
    x[N-1].i=y[N-1].i/D[N-1];
    x[N-2].r=y[N-2].r/D[N-2]-p[N-2].r*x[N-1].r-p[N-2].i*x[N-1].i;
    x[N-2].i=y[N-2].i/D[N-2]+p[N-2].i*x[N-1].r-p[N-2].r*x[N-1].i;
    for(i=N-3;i>=0;i--) {
       x[i].r=y[i].r/D[i]-p[i].r*x[i+1].r-p[i].i*x[i+1].i-q[i]*x[i+2].r;
       x[i].i=y[i].i/D[i]+p[i].i*x[i+1].r-p[i].r*x[i+1].i-q[i]*x[i+2].i;
    }  
    nml_dcv_delete(y);
    return 0;
  }
int penta_sturm(
  int iteration,
  double ref_lam,
  nml_dvector *aa,
  nml_dcscalar *b,
  nml_dvector *g,
  int *n_eigs
  ) {
    int flag=1;
    int N=2*iteration;
    nml_dvector   *D = nml_dv_new(2*iteration);
    nml_dcvector   *p = nml_dcv_new(2*iteration);
    nml_dvector   *q = nml_dv_new(2*iteration);
    *n_eigs=0;
    D[0]=aa[0]-ref_lam;
    if(D[0]<0)
      *n_eigs=*n_eigs+1;
    if(D[0]>-1e-10&&D[0]<1e-10)
      flag=0; 
    if (flag==0){
      nml_dv_delete(D);
      nml_dcv_delete(p);
      nml_dv_delete(q);
      return flag;
    }
    p[0].r=b[1].r/D[0];
    p[0].i=b[1].i/D[0];
    D[1]=aa[1]-ref_lam-(p[0].r*p[0].r+p[0].i*p[0].i)*D[0];
    if(D[1]<0)
      *n_eigs=*n_eigs+1;
    if(D[1]>-1e-10&&D[1]<1e-10)
      flag=0; 
    if(flag==0){
      nml_dv_delete(D);
      nml_dcv_delete(p);
      nml_dv_delete(q);
      return flag;
    }
    int j; 
    for (j=2;j<N;j++){
      q[j-2]=g[j]/D[j-2];
      
      p[j-1].r=(b[j].r-q[j-2]*p[j-2].r*D[j-2])/D[j-1];
      p[j-1].i=(b[j].i+q[j-2]*p[j-2].i*D[j-2])/D[j-1];
      D[j]=aa[j]-ref_lam-p[j-1].r*p[j-1].r*D[j-1]-p[j-1].i*p[j-1].i*D[j-1]
          -q[j-2]*q[j-2]*D[j-2];
      if(D[j]<0){
        *n_eigs=*n_eigs+1;
      }
      if(D[j]>-1e-10&&D[j]<1e-10)
        flag=0; 
      if(flag==0){
        nml_dv_delete(D);
        nml_dcv_delete(p);
        nml_dv_delete(q);
        return flag;
      }
    }
    nml_dv_delete(D);
    nml_dcv_delete(p);
    nml_dv_delete(q);
    return flag;
  }
int solve_pentamatrix_vector(
  int iteration,
  nml_dvector *aa,
  nml_dcvector *b,
  nml_dvector *g,
  nml_dcvector *x,    
  nml_dcvector *rhs,   
  nml_dscalar shift
  ) {
    int N=2*iteration;
    nml_dvector   *D = nml_dv_new(2*iteration);
    nml_dcvector   *p = nml_dcv_new(2*iteration);
    nml_dvector   *q = nml_dv_new(2*iteration);
    int i,j;
    nml_dscalar norm=0.0;
    for(j=0;j<N;j++){
       norm+=rhs[j].r*rhs[j].r+rhs[j].i*rhs[j].i;
    }
    nml_dscalar inv=1/sqrt(norm);
    for(j=0;j<N;j++){
       rhs[j].r=rhs[j].r*inv;
       rhs[j].i=rhs[j].i*inv;;
    }
    
    penta_LDLT_f(iteration,aa,b,g,D,p,q,shift);
    for (i=0;i<10;i++){
       penta_LDLT_s(iteration,D,p,q,x,rhs);
       norm=0.0;
       for(j=0;j<N;j++){
          norm+=x[j].r*x[j].r+x[j].i*x[j].i;
       }
       inv=1/sqrt(norm);
       for(j=0;j<N;j++){
	  rhs[j].r=x[j].r*inv;
	  rhs[j].i=x[j].i*inv;
       }
    }
    for(j=0;j<N;j++){
       x[j]=rhs[j];
    }
    nml_dv_delete(D);
    nml_dcv_delete(p);
    nml_dv_delete(q);
    return 0;
  }
int solve_pentamatrix_eigenvalue(
  int iteration,
  nml_dvector *aa,
  nml_dcvector *b,
  nml_dvector *g,
  nml_dvector *val,
  int *val_len,
  nml_dscalar  emin,
  nml_dscalar  emax,
  nml_extent requested1
  ) {
    int steps=80;
    double lower=emin;
    double upper=emax;
    double ref_lam;
    double step=(upper-lower)/steps;
    int n_eigs;
    int j=0;
    int flag=1;
    int i=0;
    nml_dscalar* bound = (nml_dscalar*)malloc(sizeof(nml_dscalar)*(iteration/2));
    int* num_eigval = (int*)malloc(sizeof(int)*(iteration/2));
    for(j=0;j<steps;j++){ 
       ref_lam=lower+j*step; 
       flag=penta_sturm(iteration,ref_lam,aa,b,g,&n_eigs); 
#ifdef DEBUGPENTA_STURM
       if(flag==0)
         printf("flag = 0 ref_lam = %1.7f\n",ref_lam);
       int looptime=0; 
       /*while(flag==0&&looptime<10){
         ref_lam=ref_lam+0.1*step;
         flag=penta_sturm(iteration,ref_lam,aa,b,g,&n_eigs);
	 looptime++;
       }*/
#endif
       if(i>0&&(n_eigs>num_eigval[i-1])){
          bound[i]=ref_lam-step;
          num_eigval[i]=num_eigval[i-1];
#ifdef DEBUGPENTA_STURM
	  printf("%d first bound = %1.15f num_eigval = %d\n",i,bound[i],num_eigval[i]);
#endif
          i++;
          bound[i]=ref_lam;
          num_eigval[i]=n_eigs;
#ifdef DEBUGPENTA_STURM
	  printf("%d first bound = %1.15f num_eigval = %d\n",i,bound[i],num_eigval[i]);
#endif
	  i++;
#ifdef DEBUGPENTA_STURM
	  if(n_eigs-num_eigval[0]>4*requested1)
	  	 break;
#endif
       }
       if(i==0){
          bound[i]=ref_lam;
          num_eigval[i]=n_eigs;
#ifdef DEBUGPENTA_STURM
	  printf("%d first bound = %1.15f num_eigval = %d\n",i,ref_lam,n_eigs);
#endif
	  i++;
       }

    }
    int real_step=i;
    int start_step=0;
    int k=real_step;
#ifdef DEBUGPENTA_STURM
    printf("\nreal_step = %d\n",real_step);
    for(i=0;i<real_step;i++)
      printf("%d first bound = %1.15f num_eigval = %d\n",i,bound[i],num_eigval[i]);
#endif
    double sub_step;
    int loop=0;
    while(loop<4){
      for(i=start_step+1;i<real_step;i++){
        if(num_eigval[i]>num_eigval[i-1]+1){
	   sub_step=(bound[i]-bound[i-1])/steps;    
	   for(j=0;j<steps;j++){
	      ref_lam=bound[i-1]+j*sub_step;
	      flag=penta_sturm(iteration,ref_lam,aa,b,g,&n_eigs);
	      if(k>real_step&&n_eigs>num_eigval[k-1]){
	         bound[k]=ref_lam-sub_step;
	         num_eigval[k]=num_eigval[k-1];
#ifdef DEBUGPENTA_STURM
                 printf("%d bound %d = %1.15f num_eigval = %d\n",loop,k,bound[k],num_eigval[k]);
#endif
	         k++;
	         bound[k]=ref_lam;
	         num_eigval[k]=n_eigs;
#ifdef DEBUGPENTA_STURM
                 printf("%d bound %d = %1.15f num_eigval = %d\n",loop,k,bound[k],num_eigval[k]);
#endif
	         k++;
	      }
	      if(k==real_step){
	         bound[k]=ref_lam;
	         num_eigval[k]=n_eigs;
#ifdef DEBUGPENTA_STURM
                 printf("%d bound %d = %1.15f num_eigval = %d\n",loop,k,bound[k],num_eigval[k]);
#endif
	         k++;
	      }
           }
        }
      }
      start_step=real_step;
      real_step=k;
      loop++;
    }

#ifdef DEBUGPENTA_STURM
    for(i=0;i<real_step;i++)
      printf("%d first bound = %1.15f num_eigval = %d\n",i,bound[i],num_eigval[i]);
#endif
    k=0;
    double ref_lam2=0.0;
    for(j=1;j<real_step;j++){
      if(num_eigval[j]==(num_eigval[j-1]+1)){
        double a1=bound[j-1];
#ifdef DEBUGPENTA_STURM
        printf("lower = %1.15f upper = %1.15f num_eigval = %d\n",bound[j-1],bound[j],num_eigval[j-1],num_eigval[j]);
#endif
        double b1=bound[j];
        i=0;
        while(b1-a1>1e-14&&i<60){
          i++;
          ref_lam=(a1+b1)/2;
          flag=penta_sturm(iteration,ref_lam,aa,b,g,&n_eigs);
          ref_lam2=ref_lam;
          if(n_eigs==num_eigval[j-1])
            a1=ref_lam2;
          else
            b1=ref_lam2;
	} 
	   
        val[k]=ref_lam2;
#ifdef DEBUGPENTA_STURM
	printf(" %d, %1.15e ref_lam2 = %1.15f\n",k,val[k],ref_lam2);
#endif
        k++;
      }
      else if(num_eigval[j]>(num_eigval[j-1]+1)&&j>=start_step){
        double a1=bound[j-1];
#ifdef DEBUGPENTA_STURM
        printf("lower = %1.15f upper = %1.15f num_eigval = %d\n",bound[j-1],bound[j],num_eigval[j-1],num_eigval[j]);
#endif
        double b1=bound[j];
        i=0;
        while(b1-a1>1e-14&&i<60){
          i++;
          ref_lam=(a1+b1)/2;
          flag=penta_sturm(iteration,ref_lam,aa,b,g,&n_eigs);
          ref_lam2=ref_lam;
            if(n_eigs==num_eigval[j-1])
              a1=ref_lam2;
            else
       	      b1=ref_lam2;
            }
          val[k]=ref_lam2;
#ifdef DEBUGPENTA_STURM
	  printf(" %d, %1.15e ref_lam2 = %1.15f\n",k,val[k],ref_lam2);
          printf(" %d, %1.15e ref_lam2 = %1.15f\n",k-1,val[k-1],ref_lam2);
#endif
          k++;
          val[k]=ref_lam2;
#ifdef DEBUGPENTA_STURM
	  printf(" %d, %1.15e ref_lam2 = %1.15f\n",k,val[k],ref_lam2);
#endif
          k++;
      }	   
    }
    *val_len=k;
#ifdef DEBUGPENTA_STURM
    for(i=0;i<k;i++){
      printf("eigenvalue = %1.15e \n",val[i]);
    }
#endif
   

    free(num_eigval);
    free(bound);
    return 0;
  }

int block_lanczos_pentadiagonal(	/*   out number of converged eigenvalues*/
    nml_dvector	 *eigv_conv,	/*   out eigenvalues			*/
    nml_dvector   *ALPHA, 
    nml_dvector   *GAMMA,
    nml_dcvector   *BETA,
    nml_extent	 *pIterations,	/*   out actual number of iterations	*/
    nml_dcscalar *r_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
    nml_dcscalar *r2_n,		/* inout r2_{n} = q_{n+1}*beta_{n}	*/
    nml_dcscalar *q_n,		/* inout current  complex Lanczos vector*/
    nml_dcscalar *q_n1,		/* inout previous complex Lanczos vector*/
    nml_extent	 length,	/* in extent of vectors r_n, q_n & q_n1	*/
    int		 requested1,	/* in requested number of eigenvalues	*/
    int		 requested2,	/* in requested number of eigenvalues	*/
    int	ConvCheckStartIter,	/* in convergence check start iteration	*/
    int ConvCheckSkipRate,	/* in convergence check skip rate	*/
    int		 imax,		/* in maximum number of iterations + 1	*/
    nml_dscalar	 emin1,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax1,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 emin2,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax2,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 tolerance,	/* in eigenvalue convergence tolerance	*/
    nml_dscalar	 resolution,	/* in eigenvalue separation minimum	*/
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const int*	 argument[],	/* matrix-vector multiply argument list	*/
    FILE	 *fp_trace,	/* in   trace log text file pointer	*/
    FILE	 *fp_log,	/* in message log text file pointer	*/
    int		 verbose	/* in verbose diagnostic messages	*/
    ) {
  /* Use Lanczos' algorithm to tridiagonalize complex hermitian matrix.	*/
  int	no_convergence = 1;	/* no convergence yet			*/
  int	n_old1 = 0;		/* number of old eigenvalues		*/
  int	n_old2 = 0;		/* number of old eigenvalues		*/
  int	iterations = 0;		/* Lanczos iteration number		*/
  int	skip = 0;		/* number of convergence checks skipped	*/
  int	skip_write = 0;		/* not actually used in this version	*/
  int	n_nosuccess1 = 0;	/* number of no success			*/
  int	n_nosuccess2 = 0;	/* number of no success			*/
  int	eigenvalues = 0;	/* number of converged eigenvalues	*/
  int	eigenvalues1 = 0;	/* number of converged eigenvalues in range1	*/
  int	eigenvalues2 = 0;	/* number of converged eigenvalues in range2	*/
  nml_dscalar alpha = 0.0;	/*		diagonal element	*/
  nml_dscalar alpha2 = 0.0;	/*		diagonal element	*/
  nml_dcscalar  beta1 ;	/*	    off diagonal element	*/
  nml_dcscalar  beta2 ;	/* previous off diagonal element	*/
  nml_dscalar  gamma1 = 0.0;	/* previous off diagonal element	*/
  nml_dscalar  gamma2 = 0.0;	/* previous off diagonal element	*/
  nml_dscalar* eigv_old1 = (nml_dscalar*)nml_malloc(sizeof(nml_dscalar)*(2*imax+2));
  nml_dscalar* eigv_old2 = (nml_dscalar*)nml_malloc(sizeof(nml_dscalar)*(2*imax+2));
  				/* variables to keep track of convergence */

  nml_dscalar  ga1 = 0.0;	/* previous off diagonal element	*/
  nml_dscalar  ga2 = 0.0;	/* previous off diagonal element	*/
  beta1.r=0.0;beta1.i=0.0;beta2.r=0.0;beta2.i=0.0;
  nml_dcscalar *q2_n=q_n1;
  nml_dcvector *q3_n = nml_dcv_new(length);
  nml_dcvector *q4_n = nml_dcv_new(length);

  if (mpi_n3d_id != mpi_n3d_masterid) {
#ifdef MPI3d
    do {
       int sndrcv[2];

      MPI_TIC(start_time2);

      single_blk_lanczos_iteration(
	matmul, argument, r_n,(nml_dcscalar*)r2_n, q_n,q2_n,(nml_dcscalar*)q3_n,(nml_dcscalar*)q4_n, length,iterations, 
        &alpha, &alpha2,&beta1, &beta2,&gamma1,&gamma2,&ga1,&ga2);
      MPI_TOC(time_lanczos_iter, start_time2);
      ALPHA[2*iterations]=alpha;
      BETA[2*iterations].r=beta1.r;
      BETA[2*iterations].i=beta1.i;
      GAMMA[2*iterations]=gamma1;
      ++iterations;
      ALPHA[2*iterations-1]=alpha2;
      BETA[2*iterations-1].r=beta2.r;
      BETA[2*iterations-1].i=beta2.i;
      GAMMA[2*iterations-1]=gamma2;
      ++skip;
      ++skip_write;

      MPI_TIC(start_time);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_TOC(time_barrier, start_time);
      MPI_Bcast(&sndrcv[0],  2, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);
      no_convergence = sndrcv[0];
      eigenvalues = sndrcv[1];
      MPI_Bcast(&eigv_conv[0], eigenvalues, MPI_DOUBLE,
                mpi_n3d_masterid, MPI_COMM_WORLD);
    } while (no_convergence);
#else
    error("Non-MPI code should not reach this statements.\n");
#endif
  }
  else {		/* (mpi_n3d_id == mpi_n3d_masterid)		*/
#ifdef WriteOutPentaMatrix
    FILE *fp=fopen("penmatrix.out","w");
#endif
    do {

      MPI_TIC(start_time2);
      single_blk_lanczos_iteration(
	matmul, argument, r_n,(nml_dcscalar *)r2_n, q_n,q2_n,(nml_dcscalar *)q3_n,(nml_dcscalar *)q4_n, length,iterations, 
        &alpha, &alpha2,&beta1, &beta2,&gamma1,&gamma2,&ga1,&ga2);
      //if(iterations==0){
      //  FILE *p1=fopen("q11_vec.m","w");
      //  FILE *p2=fopen("q21_vec.m","w");
      //  for(int i=0;i<length;i++) {
      //  fprintf(p1,"%1.15e %1.15e\n",q_n[i].r,q_n[i].i);
      //  }
      //  for(int i=0;i<length;i++) {
      //  fprintf(p2,"%1.15e %1.15e\n",q2_n[i].r,q2_n[i].i);
      //  }
      //  fclose(p1);
      //  fclose(p2);
      //}
      MPI_TOC(time_lanczos_iter, start_time2);
#ifdef WriteOutPentaMatrix
      fprintf(fp,"%1.15e %1.15e %1.15e %1.15e\n",alpha,beta1.r, beta1.i,gamma1);
      fprintf(fp,"%1.15e %1.15e %1.15e %1.15e\n",alpha2,beta2.r,beta2.i,gamma2);
#endif
      ALPHA[2*iterations]=alpha;
      BETA[2*iterations].r=beta1.r;
      BETA[2*iterations].i=beta1.i;
      GAMMA[2*iterations]=gamma1;
      ++iterations;
      ALPHA[2*iterations-1]=alpha2;
      BETA[2*iterations-1].r=beta2.r;
      BETA[2*iterations-1].i=beta2.i;
      GAMMA[2*iterations-1]=gamma2;
      
      ++skip;
      ++skip_write;

      if (verbose > 1 && (NULL != fp_log)) {
        fprintf(fp_log, "it=%d al=%9.5f be=%9.5f\n", iterations, alpha, beta1.r);
        fflush(fp_log);
        }
      else
        if (0 == iterations%100 && (NULL != fp_log)) {
          fprintf(fp_log, "Iteration count = %d\n", iterations);
          fflush(fp_log);
      }
      

      if ((ConvCheckStartIter < iterations) && (ConvCheckSkipRate <= skip)) {
        skip = 0;
        solve_pentamatrix_eigenvalue(iterations,ALPHA,BETA,GAMMA,eigv_conv,&eigenvalues1,emin1,emax1,requested1);

	find_block_converged_eigval(iterations, ALPHA,BETA,GAMMA, eigv_conv, &eigenvalues1,
	    &n_nosuccess1, &n_old1, eigv_old1, fp_trace, verbose, 
	    emin1, emax1, tolerance, resolution, fp_log,ga1,ga2);
        if (requested2 > 0) {
           solve_pentamatrix_eigenvalue(iterations,ALPHA,BETA,GAMMA,eigv_conv+eigenvalues1,&eigenvalues2,emin2,emax2,requested2);
           find_block_converged_eigval(iterations, ALPHA,BETA,GAMMA, eigv_conv+eigenvalues1 , &eigenvalues2,
                                 &n_nosuccess2, &n_old2, eigv_old2, fp_trace, verbose, 
                                 emin2, emax2, tolerance, resolution, fp_log,ga1,ga2);
        }
      }

      eigenvalues = eigenvalues1 + eigenvalues2;
      no_convergence = ((iterations + 1 < imax) 
                        && (   eigenvalues1 < requested1
                            || eigenvalues2 < requested2 )
                        );

#ifdef	MPI3d
      if (1 < mpi_n3d_numprocs) {
         int sndrcv[2];
         sndrcv[0] = no_convergence;
         sndrcv[1] = eigenvalues;

	MPI_TIC(start_time);
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_TOC(time_barrier, start_time);
	MPI_TIC(start_time);
	MPI_Bcast(&sndrcv[0], 2, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);
	MPI_Bcast(&eigv_conv[0], eigenvalues, MPI_DOUBLE, 
                  mpi_n3d_masterid, MPI_COMM_WORLD);
	MPI_TOC(time_bcast, start_time);
	if (NULL != fp_log)
	  fflush(fp_log);
	}
#endif/*MPI3d		*/
      } while (no_convergence);
#ifdef WriteOutPentaMatrix
     fclose(fp);
#endif
    }
  *pIterations = iterations;

  nml_free(eigv_old1);
  nml_free(eigv_old2);
  nml_dcv_delete(q3_n);
  nml_dcv_delete(q4_n);
  return eigenvalues;
  } /* block_lanczos_pentadiagonal */

/*  Find a few eigenvalues in a complex hermitian matrix
    using Lanczos' algorithm.						*/
nml_extent eigenvaluesBlockLanczos(	/* number of eigenvalues actually found	*/
    nml_dvector	 *value,	/*   out real eigenvalue vector		*/
    nml_dvector  *ALPHA,
    nml_dvector  *GAMMA,
    nml_dcvector *BETA,
    nml_extent	 *pIterations,	/*   out actual number of iterations	*/
    nml_dcvector *r_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
    nml_dcvector *r2_n,		/* inout r2_{n} = q_{n+1}*beta_{n}	*/
    nml_dcvector *q_n,		/* inout current  complex Lanczos vector*/
    nml_dcvector *q_n1,		/* inout previous complex Lanczos vector*/
    nml_extent	 length,	/* in extent of vectors r_n, q_n & q_n1	*/
    nml_extent	 requested1,	/* in requested number of eigenvalues	*/
    nml_extent	 requested2,	/* in requested number of eigenvalues	*/
    int	ConvCheckStartIter,	/* in convergence check start iteration	*/
    int	ConvCheckSkipRate,	/* in convergence check skip rate	*/
    nml_extent	 imax,		/* in maximum number of iterations + 1	*/
    nml_dscalar	 emin1,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax1,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 emin2,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax2,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 tolerance,	/* in eigenvalue convergence tolerance	*/
    nml_dscalar	 resolution,	/* in eigenvalue separation minimum	*/
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const int*	 argument[],	/* matrix-vector multiply argument list	*/
    FILE	 *fp_trace,	/* in   trace log text file pointer	*/
    FILE	 *fp_log,	/* in message log text file pointer	*/
    int		 verbose	/* in verbose diagnostic messages	*/
    ) {

  MPI_TIME_TIC(start_time_lanc_all);
  int eigenvalues = block_lanczos_pentadiagonal(
      value, ALPHA,GAMMA,BETA, pIterations, r_n,r2_n, q_n, q_n1, length,
      requested1, requested2, ConvCheckStartIter, ConvCheckSkipRate, imax,
      emin1, emax1, emin2, emax2, tolerance, resolution, matmul, argument,
      fp_trace, fp_log, verbose);
  MPI_TOC(time_lanczos_iter_all,
    (start_time_lanc_all + time_lanczos_setup_once));
#ifdef MPI_TIMING
  if (NULL != fp_log) {
    fprintf(fp_log, "\ntime_lanczos_iter_all = %g\n", time_lanczos_iter_all);
    fflush(fp_log);
    }
#endif /* MPI_TIMING */
  if (NULL != fp_log) {
    fprintf(fp_log, "CPU %d expected Broadcast of eigenvalues\n", mpi_n3d_id);
    fflush(fp_log);
    }
  MPI_TIC(start_time);
  MPI_Bcast(&eigenvalues, 1, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);
  MPI_TOC(time_bcast, start_time);
  if (NULL != fp_log) {
    fprintf(fp_log, "CPU %d received Broadcast of eigenvalues\n", mpi_n3d_id);
    fflush(fp_log);
    }

  /* Sort the eigenvalues in value in ascending order.			*/
  rsort(value, eigenvalues);

  return eigenvalues;
  } /* eigenvaluesLanczos	*/

static int convMethodVector=2;
/*  Compute the eigenvectors corresponding to the eigenvalues found
    in a complex hermitian matrix using Lanczos' algorithm.		*/
nml_extent eigenvectorsBlockLanczos(	/* number of eigenvalues actually found	*/
    nml_dcvector *cvalue,	/*   out complex eigenvalues		*/
    nml_dcmatrix *vector,	/*   out complex eigenvectors		*/
    const
    nml_dvector  *value,	/* in    real    eigenvalues		*/
    nml_extent	 eigenvalues,	/* in	 number of eigenvalues		*/
    nml_extent	 iterations,	/* in    number of iterations		*/
    nml_dvector   *ALPHA,
    nml_dvector   *GAMMA,
    nml_dcvector   *BETA,
    nml_dcscalar *r_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
    nml_dcvector *r2_n,		/* inout r2_{n} = q_{n+1}*beta_{n}	*/
    nml_dcscalar *q_n,		/* inout current  complex Lanczos vector*/
    nml_dcscalar *q_n1,		/* inout previous complex Lanczos vector*/
    nml_extent	 length,	/* in extent of vectors r_n, q_n & q_n1	*/
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const int*	 argument[],	/* matrix-vector multiply argument list	*/
    nml_dscalar  tolerance,     /* in eigenvalue convergence tolerance  */
    FILE	 *fp_trace,	/* in   trace log text file pointer	*/
    FILE	 *fp_log,	/* in message log text file pointer	*/
    int		 verbose	/* in verbose diagnostic messages	*/
    ) {

  //MPI_TIME_TIC(start_time_lanc_all);

  nml_dmatrix	*xm = nml_dm_new(eigenvalues, 2*iterations);
  nml_dmatrix	*xm_imag = nml_dm_new(eigenvalues, 2*iterations);
				/* real eigenvector matrix		*/
  nml_dscalar alpha = 0.0;      /*              diagonal element        */
  nml_dscalar alpha2 = 0.0;     /*              diagonal element        */
  nml_dcscalar  beta1 ; /*          off diagonal element        */
  nml_dcscalar  beta2 ; /* previous off diagonal element        */
  nml_dscalar  gamma1 = 0.0;    /* previous off diagonal element        */
  nml_dscalar  gamma2 = 0.0;    /* previous off diagonal element        */
  nml_dscalar  ga1 = 0.0;    /* previous off diagonal element        */
  nml_dscalar  ga2 = 0.0;    /* previous off diagonal element        */
  beta1.r=0.0;beta1.i=0.0;beta2.r=0.0;beta2.i=0.0;
  nml_dcscalar *q2_n=q_n1;
  nml_dcvector *q3_n = nml_dcv_new(length);
  nml_dcvector *q4_n = nml_dcv_new(length);
  int* num_iteration=NULL;
  int flag;
  /* two Method are implemented to operate on solving Ritz-vector
   *    1. all iterate at the same number of iterations
   *       for each Tridiagonal vector (old)
   *    2. selective to iterate for each Tridiagonal vector(new)
   */
   if(convMethodVector==2){
     num_iteration = (int*)nml_malloc(sizeof(int)*(eigenvalues));
   }
  if (mpi_n3d_masterid == mpi_n3d_id) {
    /* Compute the subspace eigenvectors on the master CPU only.
       Then distribute the vectors to the other CPUs			*/
    nml_dcvector	*xx = nml_dcv_new(2*iterations);
    nml_dcvector	*yy = nml_dcv_new(2*iterations);
    /*nml_dvector *oALPHA = nml_dv_new(2*iterations);
    nml_dvector *oGAMMA = nml_dv_new(2*iterations);
    nml_dcvector *oBETA = nml_dcv_new(2*iterations);*/
    /* copy:  tt <-- T */
	nml_offset	iteration = 0;
	int i;
	/*for(i=0;(unsigned int)i<2*iterations;i++) 
		oALPHA[i]=ALPHA[i];
	for(i=0;(unsigned int)i<2*iterations;i++) 
		oGAMMA[i]=GAMMA[i];
	for(i=0;(unsigned int)i<2*iterations;i++) {
		oBETA[i]=BETA[i];
	}*/
    /* Calculate eigenvectors of T for all converged eigenvzzzalues	*/
    { nml_offset	eigenvalue = 0;
      for (eigenvalue = 0; eigenvalue < eigenvalues; ++eigenvalue) {
		iteration = 0;
	cvalue[eigenvalue].r = value[eigenvalue];
	for (iteration = 0; iteration < 2*iterations; ++iteration) {
	  yy[iteration].r = 1.0*rand()/RAND_MAX;
	  yy[iteration].i = 0.0;
        }
	solve_pentamatrix_vector(iterations, ALPHA,BETA,GAMMA,xx,yy,value[eigenvalue]);

	/* Method 1 */
	if(convMethodVector==1){
	  
	  for (iteration = 0; iteration < 2*iterations; ++iteration) {
	    xm[eigenvalue][iteration] = xx[iteration].r;
	    xm_imag[eigenvalue][iteration] = xx[iteration].i;
          }
	}

	/* Method 2 */
	else if(convMethodVector==2){
	  num_iteration[eigenvalue]=iterations;
	  flag=0;
	  int num_con_checking=5;
	  for (iteration =0; iteration<2*iterations;++iteration){
	    if(iteration>0/*&&fabs(xx[iteration-1].r-xx[iteration].r)<tolerance
			    &&fabs(xx[iteration-1].i-xx[iteration].i)<tolerance*/
			    &&fabs(xx[iteration].r)<10*tolerance
			    /*&&fabs(xx[iteration].i)<10*tolerance*/){
	      if(flag<num_con_checking){
		if(iteration==(unsigned int)(num_iteration[eigenvalue]+1)&&flag<num_con_checking){
		  flag++;
		}
		else 
	          flag=0;
		num_iteration[eigenvalue]=iteration;
	      }
	    }
	    xm[eigenvalue][iteration] = xx[iteration].r;
	    xm_imag[eigenvalue][iteration] = xx[iteration].i;
	  }
	  //fprintf(fp_log,"eigenvalue = %d number of iteration = %d\n",eigenvalue,num_iteration[eigenvalue]);
	  //fflush(fp_log);
	}

	/*for(i=0;(unsigned int)i<2*iterations;i++) 
		ALPHA[i]=oALPHA[i];
	for(i=0;(unsigned int)i<2*iterations;i++) 
		GAMMA[i]=oGAMMA[i];
	for(i=0;(unsigned int)i<2*iterations;i++) 
		BETA[i]=oBETA[i];                */
      }
    }
    /* Need to distribute the subspace eigenvectors to the various CPU's
       so that they can assemble the eigenvectors locally.
       The data is given in a real matrix xm
       with eigenvalues times iterations elements.			*/
    if (NULL != fp_log) {
      fprintf(fp_log, "CPU %d starting Broadcast of eigenvectors\n", mpi_n3d_id);
      fflush(fp_log);
    }
    MPI_TIC(start_time);
    MPI_Bcast(&xm[0][0], eigenvalues*2*iterations, MPI_DOUBLE,
      mpi_n3d_masterid, MPI_COMM_WORLD);
    MPI_Bcast(&xm_imag[0][0], eigenvalues*2*iterations, MPI_DOUBLE,
      mpi_n3d_masterid, MPI_COMM_WORLD);
    MPI_TOC(time_bcast, start_time);
    if(convMethodVector==2)
      MPI_Bcast(&num_iteration[0], eigenvalues, MPI_INT,
      mpi_n3d_masterid, MPI_COMM_WORLD);
    if (NULL != fp_log) {
      fprintf(fp_log, "CPU %d finished Broadcast of eigenvectors\n", mpi_n3d_id);
      fflush(fp_log);
    }
    //fclose(fin);
    /*nml_dv_delete(oALPHA);
    nml_dv_delete(oGAMMA);
    nml_dcv_delete(oBETA);*/
    nml_dcv_delete(xx);
    nml_dcv_delete(yy);
  }

  else {			/* (mpi_n3d_masterid != mpi_n3d_id)		*/
#ifdef	MPI3d
    /* Receive the subspace eigenvectors from the master CPU.		*/
    if (NULL != fp_log) {
      fprintf(fp_log, "CPU %d expected Broadcast of eigenvectors\n", mpi_n3d_id);
      fflush(fp_log);
    }
    MPI_TIC(start_time);
    MPI_Bcast(&xm[0][0], eigenvalues*2*iterations, MPI_DOUBLE,
      mpi_n3d_masterid, MPI_COMM_WORLD);
    MPI_Bcast(&xm_imag[0][0], eigenvalues*2*iterations, MPI_DOUBLE,
      mpi_n3d_masterid, MPI_COMM_WORLD);
    MPI_TOC(time_bcast, start_time);
    if(convMethodVector==2)
      MPI_Bcast(&num_iteration[0], eigenvalues, MPI_INT,
      mpi_n3d_masterid, MPI_COMM_WORLD);
    if (NULL != fp_log) {	/* PLEASE EXPLAIN ME			*/
      fprintf(fp_log, "CPU %d received Broadcast of eigenvectors\n", mpi_n3d_id);
      fflush(fp_log);
    }
#else /*MPI3d		*/
    error("Non-MPI code should not reach this statements.\n");
#endif/*MPI3d		*/
  }

  /* Once the subspace eigenvectors are distributed,
     we can compute the overall eigenvectors.				*/

  if (mpi_n3d_id == mpi_n3d_masterid) {
    if (NULL != fp_log) {		/* PLEASE EXPLAIN ME		*/
      fprintf(fp_log,
	"Repeat Lanczos Iterations and construct eigenvectors\n"
	"iteration:\n");
      fflush(fp_log);
    }
  }
  { nml_offset iteration = 0;
    nml_dscalar dummy=0.0;
    for (iteration = 0; iteration < iterations; ++iteration) {
      if ((mpi_n3d_id == mpi_n3d_masterid) && (0 == iteration%100)) {
	if (NULL != fp_log) {
	  fprintf(fp_log, "%d..\n", iteration);
	  fflush(fp_log);
	}
      }

      MPI_TIC(start_time2);

      single_blk_lanczos_iteration(
        matmul, argument, r_n,(nml_dcscalar*)r2_n, q_n,q2_n,(nml_dcscalar*)q3_n,
	(nml_dcscalar*)q4_n, length,iteration, &alpha, &alpha2,&beta1, &beta2,
	&gamma1,&gamma2,&ga1,&ga2);

      MPI_TOC(time_lanczos_iter, start_time2);
      dummy=alpha; 
      {
/*#ifdef WriteOutVector
	FILE* fvec1=fopen("evec_1","w"); 
        FILE* fvec2=fopen("evec_2","w"); 
#endif*/
	nml_offset eigenvalue = 0;
	for (eigenvalue = 0; eigenvalue < eigenvalues; ++eigenvalue) {
	  nml_offset	j = 0;
	  if(convMethodVector==2){
	     if(iteration*2>(unsigned int)(num_iteration[eigenvalue]))
	       continue;
	  }
	  for (j = 0; j < length; ++j) {
	    vector[eigenvalue][j].r += xm[eigenvalue][2*iteration]*q_n[j].r-xm_imag[eigenvalue][2*iteration]*q_n[j].i+xm[eigenvalue][2*iteration+1]*q2_n[j].r-xm_imag[eigenvalue][2*iteration+1]*q2_n[j].i;
	    vector[eigenvalue][j].i += xm[eigenvalue][2*iteration]*q_n[j].i+xm_imag[eigenvalue][2*iteration]*q_n[j].r+xm[eigenvalue][2*iteration+1]*q2_n[j].i+xm_imag[eigenvalue][2*iteration+1]*q2_n[j].r;
/*#ifdef WriteOutVector
          if(eigenvalue==0)
	    fprintf(fvec1,"%1.15e %1.15e\n",vector[eigenvalue][j].r,vector[eigenvalue][j].i);
          if(eigenvalue==1)
	    fprintf(fvec2,"%1.15e %1.15e\n",vector[eigenvalue][j].r,vector[eigenvalue][j].i);
#endif*/
	    }
	  
	  }
/*#ifdef WriteOutVector
	fclose(fvec1);
	fclose(fvec2);
#endif*/
	}
      MPI_TIC(start_time);
      MPI_Barrier(MPI_COMM_WORLD);
      MPI_TOC(time_barrier, start_time);
      }
      nml_offset eigenvalue_n;
      nml_offset eigenvalue_n1;
      int rel=0;
      for (eigenvalue_n = 0; eigenvalue_n < eigenvalues; ++eigenvalue_n) {
	for(eigenvalue_n1=0;eigenvalue_n1<eigenvalues;++eigenvalue_n1){
        if(eigenvalue_n1>eigenvalue_n&&fabs(value[eigenvalue_n1]-value[eigenvalue_n])<1e-14){
          rel=gramSchmidt((nml_dcvector*) &(vector[eigenvalue_n][0]),(nml_dcvector*) &(vector[eigenvalue_n1][0]),length);
	  if(rel==1&&mpi_n3d_id==mpi_n3d_masterid)
	    printf("eigenvectors of %1.15f and %1.15f are orthorgonal\n",value[eigenvalue_n1],value[eigenvalue_n]);
	}
        MPI_Barrier(MPI_COMM_WORLD);
	}
      }
    }
  if (mpi_n3d_id == mpi_n3d_masterid) {
    if (NULL != fp_log) {
      fprintf(fp_log,
	"Done with Iterations and construction of eigenvectors\n");
      fflush(fp_log);
      }
    }
  if(convMethodVector==2)
    nml_free(num_iteration);
  nml_dm_delete(xm);
  nml_dm_delete(xm_imag);
  nml_dcv_delete(q3_n);
  nml_dcv_delete(q4_n);
  return eigenvalues;
  } /* eigenvectorsLanczos */

  void savePendiagonalMatrix(
      nml_dvector *oALPHA ,
      nml_dvector *oGAMMA ,
      nml_dcvector *oBETA ,
      int iterations,	/* in    Lanczos iterations			*/
      int processors,	/* in    number of MPI processors		*/
      FILE *fp_pendiag	/* inout pendiagonal matrix binary file pointer	*/
      ) {			/*    save real symmetric tridiagonal matrix T	*/
  nml_dvector	*T_write = nml_dv_new(8*iterations);
  int	iteration = 0;
      printf("%1.15f %1.15f %1.15f %1.15f\n",oALPHA[1],oBETA[1].r,oBETA[1].i,oGAMMA[1]);
  for (iteration = 0; iteration < 2*iterations; ++iteration)
  {  T_write[4*iteration  ] = oALPHA[iteration];
  }
  for (iteration = 0; iteration < 2*iterations; ++iteration)
  {  T_write[4*iteration+1] = oBETA[iteration].r;
     T_write[4*iteration+2] = oBETA[iteration].i;
  }
  for (iteration = 0; iteration < 2*iterations; ++iteration)
  {  T_write[4*iteration+3] = oGAMMA[iteration];
  }
  fwrite(&iterations, sizeof(int), 1, fp_pendiag);
  fwrite(&processors, sizeof(int), 1, fp_pendiag);
  fwrite(&T_write[0], sizeof(nml_dscalar), 8*iterations, fp_pendiag);

  nml_dv_delete(T_write);
  }	/* saveTridiagonalMatrix	*/

void restorePendiagonalMatrix(
    nml_dvector *oALPHA ,
    nml_dvector *oGAMMA,
    nml_dcvector *oBETA,
    nml_extent
	*pIterations,	/*   out Lanczos iterations pointer		*/
    nml_extent
	*pProcessors,	/* inout number of MPI processors		*/
    FILE *fp_pendiag	/* inout tridiagonal matrix binary file pointer	*/
    ) {			/* restore real symmetric tridiagonal matrix *pT*/
  fread(pIterations, sizeof(nml_extent), 1, fp_pendiag);
  fread(pProcessors, sizeof(nml_extent), 1, fp_pendiag);
  { nml_extent	iterations = *pIterations;
    iterations = *pIterations=iterations-50;
    nml_dvector	*T_read = nml_dv_new(8*iterations);
    fread(&T_read[0], sizeof(nml_dscalar), 8*iterations, fp_pendiag);
    oALPHA = nml_dv_new(2*iterations);
    oGAMMA = nml_dv_new(2*iterations);
    oBETA = nml_dcv_new(2*iterations);
    { nml_offset iteration = 0;
      for (iteration = 0; iteration < 2*iterations; ++iteration)
	oALPHA[iteration] = T_read[4*iteration];
      for (iteration = 0; iteration < 2*iterations; ++iteration)
	oBETA[iteration].r = T_read[4*iteration+1];
      for (iteration = 0; iteration < 2*iterations; ++iteration)
	oBETA[iteration].i = T_read[4*iteration+2];
      for (iteration = 0; iteration < 2*iterations; ++iteration)
	oGAMMA[iteration] = T_read[4*iteration+3];

      }

    nml_dv_delete(T_read);
    }
  }	/* restorePendiagonalMatrix	*/


/* Deprecated -- Calculates requested eigenvalues and eigenvectors
   of the Hamiltonian in the energy range [emin:emax].
   imax is the maximum number of iterations + 1.			*/
int eigensolve_blk_lanczos(		/* number of eigenvalues actually found	*/
    nml_dvector	 **pValue,	/* out real eigenvalue vector		*/
    int	ConvCheckStartIter,	/* in convergence check start iteration	*/
    int	ConvCheckSkipRate,	/* in convergence check skip rate	*/
    int		 imax,		/* in maximum number of iterations + 1	*/
    int		 requested1,	/* in requested number of eigenvalues	*/
    int		 requested2,	/* in requested number of eigenvalues	*/
    nml_extent	 length,	/* in extent of vectors r_n, q_n & q_n1	*/
    int*	 argument[],	/* matrix-vector multiply argument list	*/
    nml_dscalar	 emin1,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax1,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 emin2,		/* in eigenvalue search lower bound	*/
    nml_dscalar	 emax2,		/* in eigenvalue search upper bound	*/
    nml_dscalar	 tolerance,	/* in eigenvalue convergence tolerance	*/
    int		 verbose,	/* in verbose diagnostic messages	*/
    const char	 *fn_trace,	/* in   trace log text file name	*/
    const char	 *fn_pendiag,	/* in tridiagonal matrix binary file name */
    const char	 *fn_eigval,	/* in eigenvalue  vector binary file name */
    nml_dcvector *r_n,          /* inout r_{n} = q_{n+1}*beta_{n}       */
    const nml_dcvector *r_0,    /* in    initial guess vector           */
    nml_dcvector *q_n,          /* inout current  complex Lanczos vector*/
    nml_dcvector *q_n1,         /* inout previous complex Lanczos vector*/
    nml_dcvector **pCValue,     /*   out complex eigenvalue  vector     */
    nml_dcmatrix **pCVector,    /*   out complex eigenvector matrix     */
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    int		 restoreState	/* in read tridiagonal matrix from file	*/
    ) {
  unsigned int i = 0;
  nml_extent	iterations = 0;
  nml_extent	eigenvalues = 0;
  nml_d3bands	*T = NULL;
  nml_dvector   *ALPHA = NULL;
  nml_dvector   *GAMMA = NULL;
  nml_dcvector   *BETA = NULL;
  nml_dcvector *r2_n = NULL;
  nml_dcvector *r02_n = NULL;
  FILE		*fp_log = stdout;
  FILE		*fp_trace = (NULL == fn_trace)? NULL: fopen(fn_trace, "a");
  if ((NULL != fn_trace) && (NULL == fp_trace)) {
    error("Could not open the file \"%s\"", fn_trace);
  }

  //r02_n =  nml_dcv_new(length);
  //for(int i=0;(unsigned int)i<length;i++) {
  //   r02_n[i].r = 1.0*rand()/RAND_MAX;
  //   r02_n[i].i = 0.0;
  //}
  /* Compute the real symmetric tridiagonal matrix and eigenvalues
     or restore them from binary files saved in a previous computation.	*/
  if (restoreState) {
    /* Restore state = {tridiagonal matrix, eigenvalues}.		*/
    { /* restore real symmetric tridiagonal matrix T from file		*/
      nml_extent	processors;

      FILE *fp_pendiag = fopen(fn_pendiag, "r");
      if (NULL == fp_pendiag) {
	fprintf(stderr, "Could not open file %s\n", fn_pendiag);
	exit(EXIT_FAILURE);
	MPI_Finalize();		/* PLEASE EXPLAIN ME			*/
	}
      restorePendiagonalMatrix(ALPHA,GAMMA,BETA, &iterations, &processors, fp_pendiag);
      printf("%1.15f %1.15f %1.15f %1.15f\n",ALPHA[1],BETA[1].r,BETA[1].i,GAMMA[1]);
      if (processors != (nml_extent)mpi_n3d_numprocs) {
	if (NULL != fp_log) {
	  fprintf(fp_log,
	    "ACHTUNG!!!:  you are running on a different number of processors\n"
	    "from the run that generated your tridiagonal matrix\n");
        }
      }
      fclose(fp_pendiag);
    }

    { /* restore real eigenvalue vector *pValue	from file		*/
      FILE* fp_eigval = fopen(fn_eigval, "r");
      if (NULL == fp_eigval) {
	fprintf(stderr, "Could not open file %s\n", fn_eigval);
	exit(EXIT_FAILURE);
	MPI_Finalize();		/* PLEASE EXPLAIN ME			*/
	}
      restoreEigenvalues(pValue, &eigenvalues, fp_eigval);
      if (eigenvalues != iterations)
	if (NULL != fp_log)
	  fprintf(fp_log, "ACHTUNG!!!:  eigenvalues != iterations\n");

      fclose(fp_eigval);
    }
  }

  else { /* (!restoreState) */
    /* Compute the tridiagonal matrix and eigenvalues from scratch.	*/
    r2_n =  nml_dcv_new(length);
    r02_n =  nml_dcv_new(length);
    for(i=0;i<length;i++) {
	r2_n[i].r=r02_n[i].r= 1.0*rand()/RAND_MAX;
	r2_n[i].i=r02_n[i].i=0.0 ;
    }
    ALPHA = nml_dv_new(2*imax+2);
    GAMMA = nml_dv_new(2*imax+2);
    BETA = nml_dcv_new(2*imax+2);
    *pValue	= nml_dv_new(imax);	/* actual eigenvalues found	*/
    eigenvalues  = eigenvaluesBlockLanczos(*pValue, ALPHA,GAMMA,BETA, &iterations,
	r_n,r2_n, q_n, q_n1, length, requested1, requested2,
	ConvCheckStartIter, ConvCheckSkipRate, imax,
	emin1, emax1, emin2, emax2, tolerance, 100.0*tolerance, matmul, 
        (const int**)argument, fp_trace, fp_log, verbose);
    if (mpi_n3d_masterid == mpi_n3d_id) {
      /* Only the master processor
       * saves the pendiagonal matrix to a machine readable file and
       * saves the eigenvalues to a machine readable file.		*/

      /* Save state = {pendiagonal matrix, eigenvalues}.		*/
      { /* save real symmetric pendiagonal to file		*/
	FILE*	fp_pendiag = fopen(fn_pendiag, "w");
	if (NULL != fp_pendiag) {
           savePendiagonalMatrix(ALPHA,GAMMA,BETA, iterations, mpi_n3d_numprocs, fp_pendiag);
           fclose(fp_pendiag);
        }
      }


      { /* save real eigenvalue vector *pValue to file			*/
	FILE*	fp_eigval = fopen(fn_eigval, "w");
	if (NULL != fp_eigval) {
           saveEigenvalues(*pValue, eigenvalues, fp_eigval);
           fclose(fp_eigval);
        }
      }
    }
    iterations=iterations-50;
  }

  /* At this point,
     the real symmetric tridiagonal matrix and eigenvalues have been
     computed and saved to files or restored from files so that
     it is now possible to compute the corresponding eigenvectors.	*/
  if (NULL != r_0) {
    /* Set up the initial guess for the rerun of the Lanczos algorithm.	*/
    nml_dcscalar zero = nml_dcmplx(0.0, 0.0);
    nml_offset	j = 0;
    for (j = 0; j < length; ++j) {
      r_n[j] = r_0[j];
      r2_n[j]=r02_n[j];
      q_n1[j] = zero;
    }

    *pCValue  = nml_dcv_new(eigenvalues);
    *pCVector = nml_dcm_new(eigenvalues, length);
    eigenvectorsBlockLanczos(*pCValue, *pCVector, *pValue, eigenvalues,
	iterations, ALPHA,GAMMA,BETA, r_n,r2_n, q_n, q_n1, length, matmul,
	(const int**)argument, tolerance,fp_trace, fp_log, verbose);
    }

  if (fp_trace)   fclose(fp_trace);
  nml_dcv_delete(r2_n);
  nml_dcv_delete(r02_n);
  nml_dv_delete(ALPHA);
  nml_dv_delete(GAMMA);
  nml_dcv_delete(BETA);

  return eigenvalues;
} /* eigensolve_blk_lanczos	*/
