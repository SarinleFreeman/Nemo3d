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

THE NEWEST FILE
*****************************************************************************
$Source: /repo/eigen/src/lanczos.c,v $
*****************************************************************************/
//#define REDUCE_OPT

#include "lanczos.h"
#include "MPI_Timing.h"
#ifdef INTEL_MODIFIED
    #include <emmintrin.h>
    #include <pmmintrin.h>
    #include <memory.h>
#endif

#ifdef MPI_TIMING
    double time_lanczos_first;
    double time_lanczos_rest;
    double time_lanczos_allreduce;
	double time_lanczos_allreduce_tot = 0.0;
#endif

#ifdef NoGnuStatementExpre
void Allreduce_MPI_sp(double *pa,double *pb, int count, MPI_Datatype datatype, MPI_Op Operation, MPI_Comm comm)
{
    int i;
    for (i=0;i<count;i++){
        pb[i] = pa[i];
    }

    return;
}
#endif


/* Sort the first n elements of vector v in ascending order. */
void rsort(
    nml_dvector	*v,	// inout eigenvalue array			
    nml_extent n	// in    eigenvalue array extent		
) 
{	// bubble sort algorithm		
	//printf("*v =%x\n",v);
	//printf("n= %d\n",n);
    //for(unsigned int k=0;k<n;k++)
	//printf("at the begining rsort v[%d]=%1.15e\n",k,v[k]);

    nml_extent j = 1;
    if (n<2) return;

    for (j = 1; j <= n-1; ++j) {	// step down through the array		
        nml_dscalar	t = v[j];	// save the [possibly misplaced] element
        nml_extent i = j - 1;

        while ((0 <= i) && (t < v[i])) {	// lighter elements (bubbles)	
            v[i+1] = v[i];		// percolate back up through the array	
            --i;
        }

        v[i+1] = t;		// and put the element where it belongs	
    }

    //   for(unsigned int k=0;k<n;k++)
    //       printf("rsort v[%d]=%1.15e\n",k,v[k]);
}  // rsort() 


/* Sort the first n elements of vector v in ascending order. 
And in the same manner rearrange terms in a second vector */
void rsort_twovec(
    double	*v,	// inout eigenvalue array		
    int *second_vec ,
    int n	// in    eigenvalue array extent		
) 
{   // bubble sort algorithm		
    int j = 1;
    if (n<2) return;

    //  for (j = 0; j < n; ++j) 
    //  {
    //    printf("Inside rsort_twovec, before all: j=%d, ene=%1.15e, iter=%d\n",
    //	     j, v[j], second_vec[j]) ;
    //  }

    double t ;
    int t_from_secnd_vec ;

    for (j = 1; j <= n-1; ++j) {	// step down through the array		
        t = v[j];	// save the [possibly misplaced] element
        t_from_secnd_vec = second_vec[ j ] ;
        //  printf("At present: j=%d, t=%1.15e, tfsv=%d\n", j, t, t_from_secnd_vec) ;

        int i = j - 1;
		 
        while ((0 <= i) && (t < v[i])) {	// lighter elements (bubbles)	
            v[i+1] = v[i];		// percolate back up through the array
            second_vec[ i+1 ] = second_vec[ i ] ;
            --i;
        }

        v[i+1] = t;		// and put the element where it belongs	
        second_vec[ i+1 ] = t_from_secnd_vec ;
    }

    //  for (j = 0; j < n; ++j) 
    //  {
    //    printf("Inside rsort_twovec, after all: j=%d, ene=%1.15e, iter=%d\n",
    //     j, v[j], second_vec[j]) ;
    //  }
}  


/* Copies nml_dvector from source w to destination v.
 * It iterates copyLength times through v beginning at element vOffsetStart
 * and copying it to w beginning at element wOffsetStart.
 * Bounds checking enabled.						*/
nml_dvector* copyNmlDvector(nml_dvector* pv, nml_extent vOffsetStart,
    const nml_dvector* pw, nml_extent wOffsetStart, nml_extent copyLength)  
{
	
    nml_extent extentW = nml_dv_extent(pw);
    nml_extent extentV = nml_dv_extent(pv);
    nml_extent copyLengthW = extentW;
    nml_extent copyLengthV = extentV;

    /* Check if offset position exceeds size of vector			*/
    if (wOffsetStart >= extentW) {
        nml_message("In copyNmlDvector():: Error - start offset exceeds size of source vector!\n");
    }

    if (vOffsetStart >= extentV) {
        nml_message("In copyNmlDvector():: Error - start offset exceeds size of destination vector!\n");
    }

    /* Check if copy will exceed length of vector				*/
    if ((wOffsetStart+copyLength) > extentW) {
        nml_message("In copyNmlDvector():: copy length exceeds source dvector size! Adjusting length...\n");
        copyLengthW = extentW - wOffsetStart;
    }

    if ((vOffsetStart+copyLength) > extentV) {
        nml_message("In copyNmlDvector():: copy length exceeds source dvector size! Adjusting length...\n");
        copyLengthV = extentV - vOffsetStart;
    }

    /* Copy the largest size possible
     * which is the smaller of the two vector copyLengths			*/
    if (copyLengthW>copyLengthV) {
        copyLength = copyLengthV;
    }
    else {
        copyLength = copyLengthW;
    }

    /* Copy the vector elements						*/
    while (copyLength>0) {
        pv[vOffsetStart++] = pw[wOffsetStart++];
        copyLength--;
    }
	
    return pv;
} /* copyNmlDvector() */


/* Copies nml_d3bands from source U to destination T.
 * It iterates copyLength times through U beginning at element uOffsetStart
 * and copying it to T beginning at element tOffsetStart.
 * Bounds checking enabled.						*/
nml_d3bands* copyNmlD3bands(nml_d3bands* pT, nml_extent tOffsetStart,
    const nml_d3bands* pU, nml_extent uOffsetStart, nml_extent copyLength) 
{
    nml_extent extentT = nml_d3_extent(pT);
    nml_extent extentU = nml_d3_extent(pU);
    nml_extent copyLengthT = extentT;
    nml_extent copyLengthU = extentU;

    /* Check if offset position exceeds size of D3bands			*/
    if (tOffsetStart >= extentT) {
        nml_message("In copyNmlDvector():: Error - start offset exceeds size of source D3bands!\n");
        exit(1);
    }

    if (uOffsetStart >= extentU) {
        nml_message("In copyNmlDvector():: Error - start offset exceeds size of destination D3bands!\n");
        exit(1);
    }

    /* Check if copy will exceed length of D3bands			*/
    if ((tOffsetStart+copyLength) > extentT) {
        nml_message("In copyNmlDvector():: copy length exceeds source D3bands size! Adjusting length...\n");
        copyLengthT = extentT - tOffsetStart;
    }

    if ((uOffsetStart+copyLength) > extentU) {
        nml_message("In copyNmlDvector():: copy length exceeds source D3bands size! Adjusting length...\n");
        copyLengthU = extentU - uOffsetStart;
    }

    /* Copy the largest size possible
     * which is the smaller of the two vector copyLengths			*/
    if (copyLengthT>copyLengthU) {
        copyLength = copyLengthU;
    }
    else {
        copyLength = copyLengthT;
    }

    copyNmlDvector(pT[+1], tOffsetStart, pU[+1], uOffsetStart, copyLength);
    copyNmlDvector(pT[ 0], tOffsetStart, pU[ 0], uOffsetStart, copyLength);
    copyNmlDvector(pT[-1], tOffsetStart, pU[-1], uOffsetStart, copyLength);

    return pT;
} /* copyNmlD3bands() */



/* real symmetric [nonsingular] tridiagonal system solver		*/
int solve_rtridiag_ii(	/* Solve xT^{T} = y for x then return 0;	*/
    nml_dvector	*x,	/*   out real vector				*/
    const
    nml_d3bands	*T,	/* in    real symmetric tridiagonal matrix	*/
    const
    nml_dvector	*y	/* in    real vector				*/
) 
{
    /* solve Tx = y for x.				*/
    /* The super diagonal band is ignored since T[i+1][i] = T[i][i+1].	*/
    /* T = LU where L lower and U is upper triangular banded.		*/
    /* U[0][0] = T[0][0], U[i][i] = T[i][i] - L[i][i-1]*T[i-1][i],	*/
    /* U[i][i+1] = T[i][i+1] otherwise T[i][j] = 0			*/
    /* L[i][i] = 1, L[i+1][i] = U[i][i]/T[i+1][i], otherwise L[i][j] = 0 	*/
    /* Solve Lw = y for w then solve Ux = w for x.			*/

    const nml_extent	m = nml_d3_extent(T);
    nml_dvector	*g = nml_dv_new(m);		/* L[i+1][i] workspace	*/
    nml_dscalar	d = T[0][0];			/* U[i][i]   temporary	*/

    if ((m != nml_dv_extent(y)) || (m != nml_dv_extent(x)))
        error("In function solve_rtridiag_ii(\n"
              "  nml_dvector*, const nml_d3bands*, const nml_dvector*),\n"
              "vector and tridiagonal matrix dimensions are not equal.");

    if (0.0 == d) {
        error("In function solve_rtridiag_ii(\n"
              "  nml_dvector*, const nml_d3bands, const nml_dvector*),\n"
              "diagonal element U[0][0] = 0.\n");
    }

    /* Combine decomposition and forward substitution and divide w[j]	*/
    /* by d = U[j][j] so that isn't necessary to store all of the U[j][j]	*/
    /* until the backward substitution step.				*/
    x[0] = y[0]/d;				/* w[0]/U[0][0]		*/

    { 
	    nml_offset	j = 1;

        for (j = 1; j < m; ++j) {
            g[j] = T[-1][j]/d;	/* L[j][j-1] = U[j-1][j]/U[j-1][j-1]	*/
            d = T[0][j] - T[-1][j]*g[j];		/* U[j][j]		*/
            if (0.0 == d) {
	            /* return 1; */
                error ("In function solve_rtridiag_ii(\n"
                       "  nml_dvector*, const nml_d3bands*, const nml_dvector*),\n"
                       "diagonal element U[0][%d] = 0.\n", j);
            }

            x[j] = (y[j] - T[-1][j]*x[j-1])/d;	/* w[j]/U[j][j]		*/
        }

        for (j = m-1; 0 < j; --j) {			/* back substitution	*/
            x[j-1] -= g[j]*x[j];	/* (w[j-1] - x[j]*U[j-1][j])/U[j-1][j-1]*/
        }
    }

    nml_dv_delete(g);
    return 0;
}	/* solve_rtridiag_ii	*/


nml_dscalar InverseIteration(	/* return refined eigenvalue estimate	*/
    nml_d3bands	*h,	/* inout real symmetric tridiagonal matrix	*/
    nml_dvector	*x,	/* inout real workspace vector			*/
    nml_dvector	*y,	/* inout real eigenvector			*/
    nml_dscalar	lambda,	/* in    real eigenvalue estimate		*/
    FILE	*fp_log	/* in    message log file pointer		*/
) 
{
    nml_dscalar	de = lambda;
    //nml_dscalar	lambda0 = lambda;
    const nml_dscalar	minimum = 1e-4;		/* PLEASE EXPLAIN ME		*/
    nml_dscalar	maximum = DBL_MAX;	/* double	maximum = 1e10;	*/
    const int	iterations = 10;	/* PLEASE EXPLAIN ME		*/
    int		iteration = 0;
    const int	m = nml_dv_extent(x);
    nml_dvector	*y_old = nml_dv_new(m);
    int i;
    nml_d3bands	*hs = nml_d3_clone(h);	/* copy:  hs <-- h		*/
	
    do {
        nml_dscalar	nrm2 = 0.0;
        /*  nml_dv_ssub(hs[0], de);						*/
        for (i = 0; i < m; ++i) hs[0][i] -= de;
        solve_rtridiag_ii(x, hs, y);

        /*  nml_dv_vcpy(y_old, y);						*/
        for (i = 0; i < m; ++i) y_old[i] = y[i];

        /* y <-- x/||x|| */
        /*  nml_dv_sdiv(nml_dv_vcpy(y, x), nml_dv_sabs(x));	*/
        for (i = 0; i < m; ++i) nrm2 += x[i]*x[i];
        for (i = 0; i < m; ++i) y[i] = x[i]/sqrt(nrm2);

        /* maximum <-- infinity norm of |y| - |y_old|			*/
        maximum = 0.0;
        for (i = 0; i < m; ++i) {		/* PLEASE EXPLAIN ME		*/
            nml_dscalar	t = fabs(fabs(y[i])-fabs(y_old[i]));
            maximum = (maximum < t)? t: maximum;
        }

        #if 0
            lambda0 = lambda;
            lambda = 0.0;
            for (i = 0; i < m; ++i) lambda += y[i]*y[i]*h[0][i];
            for (i = 1; i < m; ++i) lambda += 2*y[i-1]*y[i]*h[-1][i];	// assume symmetric matrix
            de = lambda - lambda0;

            if (NULL != fp_log)
                fprintf(fp_log, "%25.18e   %25.18e  %e\n", lambda, de, maximum);
        #endif

        ++iteration;
    } while((minimum < maximum) && (iteration < iterations));

    for (i = 0; i < m; ++i) x[i] = y[i];

    if (iteration >= iterations)
        if (NULL != fp_log)
            fprintf(fp_log, "Bailed out on inverse iteration\n");

    nml_dm_delete(hs);
    nml_dv_delete(y_old);

    return lambda;
}	/* InverseIteration		*/


int originalInverseIteration(	/* returns 0				*/
    nml_d3bands	*T,	/* inout real symmetric tridiagonal matrix	*/
    nml_dvector	*x,	/* inout real workspace vector			*/
    nml_dvector	*y,	/* inout real eigenvector			*/
    nml_dscalar	lambda	/* in    real eigenvalue estimate		*/
) 
{
    /* The inverse iteration algorithm is used to find			*/
    /* the eigenvector y corresponding to the eigenvalue lambda.		*/
  
    const nml_dscalar de = 1e-8;		/* PLEASE EXPLAIN ME		*/

    const int m = nml_dv_extent(x);
    int i,j;
    int iterations=20;
    nml_dscalar nrm2=0.0;

    nml_dv_ssub(T[0], lambda+de);

    for (i = 0; i < m; ++i) nrm2 += y[i]*y[i];
    for (i = 0; i < m; ++i) y[i] = y[i]/sqrt(nrm2);

    for(j=0;j<iterations;j++){
        if (solve_rtridiag_ii(x, T, y)==0) {
            // nml_dv_ssub(T[0], lambda);
            nrm2=0.0;

            for (i = 0; i < m; ++i) nrm2 += x[i]*x[i];
            for (i = 0; i < m; ++i) y[i] = x[i]/sqrt(nrm2);
        }
        else
            printf("original inverse iteration bailed out\n");
    }

	nrm2=0.0;
    for (i = 0; i < m; ++i) x[i] = y[i];
	
    //for (i = 0; i < m; ++i) nrm2 += x[i]*x[i];

    return 0;
}	/* originalInverseIteration	*/


void single_lanczos_iteration(
    void (*mul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const int* argument[],	/* in    matrix-vector multipier argument list	*/
    
    #ifdef REDUCE_OPT
        int iteration,
    #endif // REDUCE_OPT

    nml_dcscalar *r_n,	/* inout r_{n-1} = q_n*beta_{n-1}		*/
    nml_dcscalar *q_n,	/*   out current  complex Lanczos vector q_{n}	*/
    nml_dcscalar *q_n1,	/* inout previous complex Lanczos vector q_{n-1}*/
    nml_extent	 m,	/* in    extent of vectors r_n, q_n and q_{n-1}	*/
    nml_dscalar	 *alpha,/*   out real  on-diagonal element alpha_{n}	*/
    nml_dscalar	 *beta,	/*   out real off-diagonal element beta_{n}	*/
    nml_dscalar beta_n1/* in    real off-diagonal element beta_{n-1}	*/
) 
{
    /* Lanczos algorithm (iteration n)		*/
    
    #ifdef REDUCE_OPT
    /* q_{0} = 0, beta_{0} = 0
       if iteration = 0
            r_{n} = A * q_{n};
        alpha_{n} = q_{n}' * r_{n}; 
	 
	   else
	    for n = 1, 2, 3, ...
	        q_{n} = r_{n} - alpha_{n} * q_{n} - beta_{n-1} * q_{n-1};
	        r_{n} = A * q_{n};
	     beta_{n} = norm(q_n);         // ALLREDUCE
	    alpha_{n} = q_{n}' * r_{n};   // ALLREDUCE
	        q_{n} = q_{n} / beta_{n};
	        r_{n} = r_{n} / beta_{n};
	    alpha_{n} = alpha_{n} / beta_{n}^2;
	*/
    #else  
    /* q_{0} = 0, r_{0} = arbitrary and beta_{0} = ||r_{0}||.
	   for n = 1, 2, 3, . . .
	      q_{n} <-- r_{n-1}/beta_{n-1}
	      r_{n} <-- A*q_{n}		// matrix-vector product
	  alpha_{n} <-- <q_{n} | r_{n}>	// vector-vector product
	      r_{n} <-- r_{n} - beta_{n-1}*q_{n-1} - alpha_{n}*q_{n}
	   beta_{n} <-- ||r_{n}||
    */
    #endif // REDUCE_OPT

    #ifdef REDUCE_OPT
    // First iteration
    if(iteration == 0)
    {
	    // r_{1} = A * q_{1}
	    mul(argument, r_n, q_n);	/* do matrix-vector multiply		*/
	
		/* alpha_{1} = <q_{1} | r_{1}>					*/
	    nml_dcscalar alc = { 0.0, 0.0 };
	    nml_dscalar	alphap = 0.0;
	    nml_offset	j = 0;

	    for (j = 0; j < m; ++j) {
	        /* alc = cadd(alc, cmul(cconj(q_n[j]), r_n[j]));			*/
	        nml_dscalar *r_n_j = &(r_n[j].r);
	        nml_dscalar *q_n_j = &(q_n[j].r);
	        nml_dscalar *q_n1_j = &(q_n1[j].r);

	        alc.r += (*r_n_j)*(*q_n_j);
	        alc.i += r_n_j[1]*(*q_n_j);

	        ++q_n_j;
	        alc.i -= (*r_n_j)*(*q_n_j);
	        ++r_n_j;
	        alc.r += (*r_n_j)*(*q_n_j);

			q_n1_j = 0;
			q_n1_j++;
			q_n1_j = 0;
			
	    }
	    
	    alphap = alc.r;

        /* Now we need to assemble one alpha. Just sum them up.		*/ 
		nml_dscalar	al = 0.0;

	    Allreduce_MPI_sp(&alphap, &al, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	    *alpha = al;
    }
    else
    #endif // REDUCE_OPT
    {
	    #ifdef REDUCE_OPT
        nml_offset j = 0;

        /* q_{n+1} = r_{n} - alpha_{n} * q_{n} - beta_{n-1} * q_{n-1}; */
		#ifdef INTEL_MODIFIED
			#pragma message ("INTEL_MODIFED: Lanczos loop vectorization is ON")
		    __m128d __r_nj, __q_nj, __q_n1j, __al, __beta_n1;

	    	__al = _mm_load_sd((double *)alpha);
			__al = _mm_unpacklo_pd(__al, __al);

			__beta_n1 = _mm_load_sd((double *)&beta_n1);
			__beta_n1 = _mm_unpacklo_pd(__beta_n1, __beta_n1);

			for( j = 0; j < m; j++)
			{
                __q_nj = _mm_load_pd((double *)&q_n[j]);
                __r_nj = _mm_load_pd((double *)&r_n[j]);
                __q_n1j = _mm_load_pd((double *)&q_n1[j]);

				_mm_store_pd((double *)&q_n1[j], __q_nj);

                __q_nj = _mm_mul_pd(__al, __q_nj);
				__q_n1j = _mm_mul_pd(__beta_n1, __q_n1j);
				__r_nj = _mm_sub_pd(__r_nj, __q_nj);
				__q_nj = _mm_sub_pd(__r_nj, __q_n1j);

				_mm_store_pd((double *)&q_n[j], __q_nj);
			}
		#else
		    for( j = 0; j < m; j++ )
		    {
                nml_dscalar al = *alpha;

		        register nml_dscalar REq = q_n[j].r;
			    register nml_dscalar IMq = q_n[j].i;
			
			    nml_dscalar *pr_n = &(r_n[j].r);
			    nml_dscalar *pq_n1 = &(q_n1[j].r);
			    nml_dscalar *pq_n = &(q_n[j].r);
			
			    *pq_n = *pr_n - al * REq - (beta_n1) * (*pq_n1);
				*pq_n1 = REq;
				
			    ++pr_n;
				++pq_n;
			    ++pq_n1;
			
		        *pq_n = *pr_n - al * IMq - (beta_n1) * (*pq_n1);    
				*pq_n1 = IMq;
            }		
		#endif // INTEL_MODIFIED
		
		/* r_{n+1} = A * q_{n+1}; */
	    mul(argument, r_n, q_n);
	
	    /*
	           beta_{n} = norm(q_{n+1});         // ALLREDUCE
            alpha_{n+1} = q_{n+1}' * r_{n+1};    // ALLREDUCE
        */
	    { 
		    nml_dcscalar alc = { 0.0, 0.0 };
			nml_dscalar alphap = 0.0;
			nml_dscalar betap = 0.0;
		
		    for (j = 0; j < m; ++j) {
		        /* alc = cadd(alc, cmul(cconj(q_n[j]), r_n[j]));			*/
		        nml_dscalar *r_n_j = &(r_n[j].r);
		        nml_dscalar *q_n_j = &(q_n[j].r);

		        alc.r += (*r_n_j)*(*q_n_j);
		        alc.i += r_n_j[1]*(*q_n_j);
		        
				betap += (*q_n_j) * (*q_n_j);
				
		        ++q_n_j;
		        alc.i -= (*r_n_j)*(*q_n_j);
		        ++r_n_j;
		        alc.r += (*r_n_j)*(*q_n_j);
		        
				betap += (*q_n_j) * (*q_n_j);
		    }
		
		    alphap = alc.r; 

		    /* Now we need to assemble one alpha. Just sum them up.		*/
		    { 
				nml_dscalar input[2] = { 0.0, 0.0 };
				nml_dscalar output[2] = { 0.0, 0.0 };
				input[0] = alphap; input[1] = betap;

		        Allreduce_MPI_sp(input, output, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		        *alpha = output[0];
		        *beta = sqrt(output[1]);
		    }
		}
		    
		/*    q_{n+1} = q_{n+1} / beta_{n};
	          r_{n+1} = r_{n+1} / beta_{n};
	      alpha_{n+1} = alpha_{n+1} / beta_{n}^2; 
	    */
	    {
		    #ifdef INTEL_MODIFIED
			    __m128d __inv_beta, __r_nj, __q_nj;
				nml_dscalar inv_beta = 1.0/(*beta);

				__inv_beta = _mm_load_sd((double *)&inv_beta);
				__inv_beta = _mm_unpacklo_pd(__inv_beta, __inv_beta);

				for(j = 0; j < m; j++)
				{
				    __r_nj = _mm_load_pd((double *)&r_n[j]);
				    __q_nj = _mm_load_pd((double *)&q_n[j]);

					__r_nj = _mm_mul_pd(__r_nj, __inv_beta);
					__q_nj = _mm_mul_pd(__q_nj, __inv_beta);

					_mm_store_pd((double *)&r_n[j], __r_nj);
					_mm_store_pd((double *)&q_n[j], __q_nj);
				}
			#else
			    for(j = 0; j < m; j++)
			    {
				    nml_dscalar *r_n_j = &(r_n[j].r);
		            nml_dscalar *q_n_j = &(q_n[j].r);
		
				    *r_n_j /= (*beta);
				    *q_n_j /= (*beta);

				    r_n_j++;
				    q_n_j++;

				    *r_n_j /= (*beta);
				    *q_n_j /= (*beta);
			    }
			#endif // INTEL_MODIFIED

		    *alpha /= ((*beta)*(*beta));
	    }    
        #else

        // original code comes here
		nml_offset	j = 0;
		
		#ifdef INTEL_MODIFIED
		    #include "lanczos_SSE_1.c"
		#else
		    for (j = 0; j < m; ++j) {	/* q_{n} <-- r_{n-1}/beta_{n-1}		*/
		        /*
		        q_n[j].r = r_n[j].r/beta_n1;
		        q_n[j].i = r_n[j].i/beta_n1;
		        */

		        nml_dscalar	inv_beta_n1 = 1.0/beta_n1;
		        /* first do the real part ... */
		        nml_dscalar	*q_n_j = &(q_n[j].r);
		        nml_dscalar	*r_n_j = &(r_n[j].r);
		        *q_n_j = *r_n_j*inv_beta_n1;
		        /* now do the imaginary part */
		        ++q_n_j;
		        ++r_n_j;
		        *q_n_j = *r_n_j*inv_beta_n1;
		    }
		#endif // INTEL_MODIFIED
		
		mul(argument, r_n, q_n);	/* do matrix-vector multiply		*/

		{ /* alpha_{n} <-- <q_{n} | r_{n}>					*/
		    nml_dcscalar	alc = {0.0, 0.0};
		    nml_dscalar	alphap = 0.0;
		    nml_offset	j = 0;
		    for (j = 0; j < m; ++j) {
		        /* alc = cadd(alc, cmul(cconj(q_n[j]), r_n[j]));			*/
		        nml_dscalar *r_n_j = &(r_n[j].r);
		        nml_dscalar *q_n_j = &(q_n[j].r);

		        alc.r += (*r_n_j)*(*q_n_j);
		        alc.i += r_n_j[1]*(*q_n_j);
		        ++q_n_j;
		        alc.i -= (*r_n_j)*(*q_n_j);
		        ++r_n_j;
		        alc.r += (*r_n_j)*(*q_n_j);
		    }

		    alphap = alc.r;

		    /* Now we need to assemble one alpha. Just sum them up.		*/
		    { 
			    nml_dscalar	al = 0.0;


		        Allreduce_MPI_sp(&alphap, &al, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		        *alpha = al;
		    }
		}

		{ /* compute r_{n} = q_{n+1}*beta_{n} and beta_{n} = ||r_{n}||	*/
		    nml_offset	j = 0;
		
		    #ifdef INTEL_MODIFIED
		    #include "lanczos_SSE_2.c"
		    #else
		    nml_dscalar	al = *alpha;
		    nml_dscalar	betap = 0.0;
		
		    for (j = 0; j < m; ++j) {
		        
		        /*
		        r_n[j].r -= al*q_n[j].r + beta_n1*q_n1[j].r;
		        r_n[j].i -= al*q_n[j].i + beta_n1*q_n1[j].i;
		        q_n1[j] = q_n[j];
		        betap += r_n[j].r*r[j].r + r_n[j].i*r_n[j].i;
		        */

		        register nml_dscalar REq = q_n[j].r;
		        register nml_dscalar IMq = q_n[j].i;

		        nml_dscalar *pr_n  = &(r_n[j].r);
		        nml_dscalar *pq_n1 = &(q_n1[j].r);

		        *pr_n -= al*REq + beta_n1*(*pq_n1);
		        betap += (*pr_n)*(*pr_n);

		        ++pr_n;
		        ++pq_n1;

		        *pr_n -= al*IMq + beta_n1*(*pq_n1);
		        betap += (*pr_n)*(*pr_n);

		        pq_n1[ 0] = IMq;
		        pq_n1[-1] = REq;
		    }
		    #endif // INTEL_MODIFIED
		    
		    { /* Assemble one beta. Sum them and then take sq root.		*/
		        nml_dscalar	bl = 0.0;

		        Allreduce_MPI_sp(&betap, &bl, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		        *beta = sqrt(bl);
		    }
		}
        #endif // REDUCE_OPT
	}
}	/* single_lanczos_iteration	*/


static int convMethodValue=2;

/* MODIFIED BY MAREK KORKUSINSKI, APRIL 11 2005 */

void find_converged_eigval(
    int		iterations,	/* in Lanczos iteration			*/
    const
    nml_d3bands	*T,		/* in real symmetric tridiagonal matrix	*/
    nml_dvector	*eigv_conv,	/* inout actual converged eigenvalues	*/
    nml_ivector * convIterationNumbers ,
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
    FILE	*fp_log		/* in message log text file pointer	*/
) 
{
    const nml_dscalar	bias = 1e-80;	/* PLEASE EXPLAIN ME			*/
    const nml_dscalar	gain = 1e-8;	/* PLEASE EXPLAIN ME			*/
    nml_dscalar	        deigv_min = DBL_MAX; /*double	deigv_min = 1.0e100;	*/
    nml_dscalar	        eigv_min = DBL_MAX; /*double	 eigv_min = 1.0e100;	*/
    int		        Neval_in_range = 0;	/* number eigenvalues found	*/
    nml_dvector	        *eval = nml_dv_new(iterations); /* actual eigenvalues found	*/

    /*  BEGIN BLOCK DIAGONALIZING THE TRIDIAGONAL MATRIX
        USING BLAS SUBROUTINE */

    const char	range[2] = "V";	/* find all eigenvalues in (vl, vu]	*/
    const char	order[2] = "E";	/* order eigenvalues for entire matrix	*/
    nml_dscalar	vl = emin/MTRX_SCALE;	/* search region lower bound	*/
    nml_dscalar	vu = emax/MTRX_SCALE;	/* search region upper bound	*/
    int		il = 0;		/* ignored search index lower bound	*/
    int		iu = 0;		/* ignored search index upperbound	*/
    nml_dscalar	abstol = 1.0e-13;	/* eigenvalue absolute tolerance*/
    int		info = 0;	/* eigenvalue convergence information	*/
    int		nsplit = 0;	/* number of diagonal blocks in matrix	*/
    nml_ivector	*iblock = nml_iv_new(iterations); /* eigenvalue block membership	*/
    nml_ivector	*isplit = nml_iv_new(iterations); /* index of split between blocks*/
    nml_dvector	*work = nml_dv_new(4*iterations); /* Fortran subroutine workspace	*/
    nml_ivector	*iwork = nml_iv_new(3*iterations);/* Fortran subroutine workspace	*/
    nml_dvector	*diag = nml_dv_clone(T[0]);	      /* Diagonal elements*/
    nml_dvector	*offd = nml_dv_clone(T[1]);	      /* Off diagonal elements*/

    /* BLAS: Compute the eigenvalues of a symmetric tridiagonal matrix.	*/
    dstebz(range, order, &iterations, &vl, &vu, &il, &iu, &abstol, diag, offd,
           &Neval_in_range, &nsplit, eval, iblock, isplit, work, iwork, &info);
  
    nml_dv_delete(offd);
    nml_dv_delete(diag);
    nml_iv_delete(iwork);
    nml_dv_delete(work);
    nml_iv_delete(isplit);
    nml_iv_delete(iblock);

    *n_nosuccess = 0;

    /* convergence Method for eigenvalue
     * there are two Methods are implemented:
     * 1. check if the eigenvalue converge or not (old)
     * 2. check if the final component of the eigenvector 
     * is less than tolerance (new)
     */

    /*Method 1*/ 
    if(convMethodValue==1)
    { 
        *pEigenvalues = 0;
        int	i = 0;

        for (i = 0; i < Neval_in_range; ++i) {
            if ((0 < i) && (fabs(eval[i-1] - eval[i]) < fabs(bias + eval[i])*gain))
	            continue;		// ignore other degenerate eigenvalues

            if (NULL != fp_trace) {
                fprintf(fp_trace, "%d %1.15e\n", iterations, eval[i]);
                fflush(fp_trace);
            }

            { 
                int convskip = 0;
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

                        if (deigv < tolerance) {
						    //check if the eigenvalue converges
                            eigv_conv[*pEigenvalues] = eval[i];
                            convIterationNumbers[ *pEigenvalues ] = iterations ;

                            if (verbose && (NULL != fp_log)) {
                                fprintf(fp_log, "   eigenvalue candidate #%d:  %18.15f\n", *pEigenvalues, eigv_conv[*pEigenvalues]*MTRX_SCALE);
                                fflush(fp_log);
                            }

                            ++(*pEigenvalues);
                            *n_nosuccess = 0;
    }}}}}}

    /*Method 2*/
    else if(convMethodValue==2) 
    {
        //    iterations=iterations-50;
        nml_dvector *xx = nml_dv_new(iterations);
        nml_dvector *yy = nml_dv_new(iterations);
        nml_d3bands *tt = nml_d3_new(iterations);
        int iteration;
        int flag;
        int	i = 0;

        for (i = 0; i < Neval_in_range; ++i) {
            /* SKIP DEGENERATE EIGENVALUES */
            if ((0 < i) && (fabs(eval[i-1] - eval[i]) < fabs(bias + eval[i])*resolution))
            { 
                // cout << "Degenerate! " << eval[i-1] << "   " << eval[i] << endl ;
                // ignore other degenerate eigenvalues
            }
            else
            {  
                /* begin block 'if not degenerate' */

                // Begin Srikant addition for tracing eigen values - Sept 06, 2006
                if (NULL != fp_trace) {
                    fprintf(fp_trace, "%d %1.15e\n", iterations, eval[i]);
                    fflush(fp_trace);
                }

                //End trace of Eigen value by Srikant
                /* check if this eigenvalue has already converged */

                int alreadyConvergedFlag = 0 ;

                for ( int indexj = 0 ; indexj < *pEigenvalues ; indexj++ )
                {
                    if (fabs(eigv_conv[indexj] - eval[i]) < resolution )
                        alreadyConvergedFlag = 1 ;
                }

                if ( 0 == alreadyConvergedFlag )
                {   
                    /* this eigenvalue has not yet been converged! */
                    flag=1;
					
                    /* copy:  tt <-- T */
                    copyNmlD3bands(tt, 0, T, 0, iterations);
					
                    for (iteration = 0; iteration < iterations; ++iteration) yy[iteration] = 1.0*rand()/RAND_MAX;
                          
                    originalInverseIteration(tt, xx, yy, eval[i]);

                    for (iteration = iterations-2; iteration < iterations; ++iteration) {
                        if(fabs(xx[iteration])>tolerance)
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
						
                        {
                            eigv_conv[*pEigenvalues]=eval[i];
                            convIterationNumbers[ *pEigenvalues ] = iterations ;

                            if (verbose && (NULL != fp_log)) {
                                fprintf(fp_log, "   eigenvalue candidate #%d:  %18.15f\n",
                                        *pEigenvalues, eigv_conv[*pEigenvalues]*MTRX_SCALE);
                                fflush(fp_log);
                            }
							
                            ++(*pEigenvalues);
                        }
                    }

                    #ifdef Timing_detail
                        fprintf(fp_log,"%d = %1.15f\n",i,eval[i]);
                    #endif
                } /* end block "this eigenvalue has not yet been converged" */
            } /* end block 'if not degenerate' */
        }   /* end loop over eigenvalues */

        nml_dv_delete(xx);
        nml_dv_delete(yy);
        nml_d3_delete(tt);
    }
	
    if (verbose) {
        /*if (deigv_min != 1e100)*/
        if (deigv_min < DBL_MAX)
            if (NULL != fp_log)
                fprintf(fp_log, "    E = %e  dE_min = %e\n", eigv_min*MTRX_SCALE, deigv_min*MTRX_SCALE);
    }

    *n_old = 0;

    { 
        int	iteration = 0;
        for (iteration = 0; iteration < iterations; ++iteration)
            if ((emin/MTRX_SCALE < eval[iteration])
             && (eval[iteration] < emax/MTRX_SCALE)
             && ((0 == iteration) || (resolution < fabs((eval[iteration-1] - eval[iteration])/(bias + eval[iteration]))))) 
            {
                eigv_old[*n_old] = eval[iteration];
                ++(*n_old);
            }
    }
    
    ++(*n_nosuccess);

    nml_dv_delete(eval);
} /* find_converged_eigval */



int lanczos_tridiagonal(	/*   out number of converged eigenvalues*/
    nml_dvector	 *eigv_conv,	/*   out eigenvalues			*/
    nml_ivector * convIterationNumbers ,
    nml_d3bands	 *T,		/* inout symmetric tridiagonal matrix	*/
    nml_extent	 *pIterations,	/*   out actual number of iterations	*/
    nml_dcscalar *r_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
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
) 
{
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
    nml_dscalar  beta = 0.0;	/*	    off diagonal element	*/
    nml_dscalar  beta_n1 = 0.0;	/* previous off diagonal element	*/
    nml_dscalar* eigv_old1 = (nml_dscalar*)malloc(sizeof(nml_dscalar)*imax);
    nml_dscalar* eigv_old2 = (nml_dscalar*)malloc(sizeof(nml_dscalar)*imax);
	
    /* variables to keep track of convergence */
    nml_dscalar  beta2 = nml_dcv_ssqr(r_n);	/* <r_n | r_n>		*/
    nml_dvector * eigv_conv1 = nml_dv_new( imax ) ;
    nml_dvector * eigv_conv2 = nml_dv_new( imax ) ;
    nml_ivector * convIterationNumbers1 = nml_iv_new( imax ) ;
    nml_ivector * convIterationNumbers2 = nml_iv_new( imax ) ;

    #ifndef REDUCE_OPT
        Allreduce_MPI_sp(&beta2, &beta_n1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta_n1 = sqrt(beta_n1);
    #endif // REDUCE_OPT

    if (mpi_n3d_id != mpi_n3d_masterid) {

        #ifdef MPI3d
        do {
            int sndrcv[2];

            #ifdef Timing_detail
                MPI_TIC(start_time2);
            #endif

            #ifdef MPI_TIMING
                MPI_TIC(start_time2);
            #endif

            single_lanczos_iteration(matmul, argument, 
	            #ifdef REDUCE_OPT
	                iterations,
	            #endif // REDUCE_OPT
	            r_n, q_n, q_n1, length, &alpha, &beta, beta_n1);

            #ifdef Timing_detail
            MPI_TOC(time_lanczos_iter, start_time2);
            #endif

            #ifdef MPI_TIMING
            if ( iterations == 0 )
                MPI_TOC(time_lanczos_first, start_time2);
            else
                MPI_TOC(time_lanczos_rest, start_time2);
            #endif
      
            beta_n1 = beta;
            ++iterations;
            ++skip;
            ++skip_write;

            MPI_Barrier(MPI_COMM_WORLD);
            MPI_Bcast(&sndrcv[0],  2, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

            no_convergence = sndrcv[0];
            eigenvalues = sndrcv[1];
            MPI_Bcast(&eigv_conv[0], eigenvalues, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);
            MPI_Bcast(&convIterationNumbers[0] , eigenvalues , MPI_INT , mpi_n3d_masterid, MPI_COMM_WORLD);

            /*
            if (NULL != fp_log) {
	            fprintf(fp_log, "cpuid=%d Received no_convernce=%d\n", mpi_n3d_id, no_convergence);
	            fflush(fp_log);
            }
            */
			
        } while (no_convergence);
        #else
        error("Non-MPI code should not reach this statements.\n");
        #endif
    }
    else {			/* (mpi_n3d_id == mpi_n3d_masterid)		*/
        #ifdef Timing_detail
        time_lanczos_iter=0.0;
        #endif

        #ifdef MPI_TIMING
        time_lanczos_first=0.0;
        time_lanczos_rest=0.0;
        #endif

        do {
            //int sndrcv[2];
            #ifdef Timing_detail
            MPI_TIC(start_time2);
            #endif

            single_lanczos_iteration( matmul, argument, 
	            #ifdef REDUCE_OPT
	            iterations,
                #endif // REDUCE_OPT
                r_n, q_n, q_n1, length, &alpha, &beta, beta_n1);

            #ifdef Timing_detail
            MPI_TOC(time_lanczos_iter, start_time2);
            #endif

            #ifdef MPI_TIMING
            if ( iterations == 0 )
                MPI_TOC(time_lanczos_first, start_time2);
            else
                MPI_TOC(time_lanczos_rest, start_time2);
            #endif

            #ifdef REDUCE_OPT
            T[0][iterations] = alpha;

            if(iterations != 0)
            {
                T[1][iterations-1] = beta;
                T[-1][iterations] = beta;
            }
            #else
            T[0][iterations] = alpha;
            T[1][iterations] = T[-1][iterations+1] = beta;
            #endif // REDUCE_OPT

			//printf("#### alpha=%f beta=%f\n",alpha,beta);

            beta_n1 = beta;
            ++iterations;
            ++skip;
            ++skip_write;

            if (verbose > 1 && (NULL != fp_log)) {
                fprintf(fp_log, "it=%d al=%9.5f be=%9.5f\n", iterations, alpha, beta);
                fflush(fp_log);
            }
            else
                if (0 == iterations%100 && (NULL != fp_log)) {
                    fprintf(fp_log, "Iteration count = %d of %d maximum\n", iterations,imax);
                    fflush(fp_log);
                }

            #ifdef Timing_detail
            int rate=0.0; 
            double find_converge=0.0;
            MPI_TIC(start_time3);
            #endif

            if ((ConvCheckStartIter < iterations) && (ConvCheckSkipRate <= skip)) {
                skip = 0;

                #ifdef Timing_detail
                rate++;
                #endif

                if ( ( requested1 > 0 ) && ( eigenvalues1 <= requested1 ) ) 
                {
                    find_converged_eigval(iterations, T, eigv_conv1, convIterationNumbers1, &eigenvalues1,
                                          &n_nosuccess1, &n_old1, eigv_old1, fp_trace, verbose, 
                                          emin1, emax1, tolerance, resolution, fp_log);
                }
                else
                {
                    if ( requested1 == 0) 
                        eigenvalues1 = 0 ;
                }

                if ( (requested2 > 0) && ( eigenvalues2 <= requested2 ) ) 
                {
                    find_converged_eigval(iterations, T, eigv_conv2, 
                                 convIterationNumbers2 ,
                                 &eigenvalues2,
                                 &n_nosuccess2, &n_old2, eigv_old2, fp_trace, verbose, 
                                 emin2, emax2, tolerance, resolution, fp_log);
                }
                else
                {
                    if ( requested2 == 0) {
                        eigenvalues2 = 0 ;
                }}

                for (int indexj = 0 ; indexj < eigenvalues1 ; indexj++)
                {
                    eigv_conv[ indexj ] = eigv_conv1[ indexj ] ;
                    convIterationNumbers[ indexj ] = convIterationNumbers1[ indexj ] ;
                }
	
                for (int indexj = 0 ; indexj < eigenvalues2 ; indexj++)
                {
                    eigv_conv[ indexj + eigenvalues1 ] = eigv_conv2[ indexj ] ;
                    convIterationNumbers[ indexj + eigenvalues1 ] = convIterationNumbers2[ indexj ] ;
                }
            }

            #ifdef Timing_detail
            MPI_TOC(find_converge,start_time3);
            printf("total find_converge time = %.2f total find_converge called %d times\n",find_converge,rate);
            #endif
			
            eigenvalues = eigenvalues1 + eigenvalues2;

            no_convergence = ((iterations + 1 < imax) 
                        && ( eigenvalues1 < requested1 || eigenvalues2 < requested2 ));

            #ifdef MPI3d
            if (1 < mpi_n3d_numprocs) {
                int sndrcv[2];
                sndrcv[0] = no_convergence;
                sndrcv[1] = eigenvalues;

                MPI_Barrier(MPI_COMM_WORLD);

                MPI_Bcast(&sndrcv[0], 2, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);
                MPI_Bcast(&eigv_conv[0], eigenvalues, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);
                MPI_Bcast(&convIterationNumbers[0] , eigenvalues , MPI_INT , mpi_n3d_masterid, MPI_COMM_WORLD);

                if (NULL != fp_log)
                    fflush(fp_log);
            }
            #endif/*MPI3d		*/

        } while (no_convergence);

        #ifdef Timing_detail
        printf("time of lanczos iteration in first sweep = %.2f\n",time_lanczos_iter);
        #endif

        #ifdef MPI_TIMING
        printf("Lanczos_1st = %.2f Lanczos_rest = %.2f\n", time_lanczos_first, time_lanczos_rest);
        #endif
    }

    *pIterations = iterations;

    free(eigv_old1);
    free(eigv_old2);

    nml_iv_delete( convIterationNumbers2 ) ;
    nml_iv_delete( convIterationNumbers1 ) ;

    nml_dv_delete( eigv_conv1 ) ;
    nml_dv_delete( eigv_conv2 ) ;

    return eigenvalues;
} /* lanczos_tridiagonal */

/*  Find a few eigenvalues in a complex hermitian matrix
    using Lanczos' algorithm.						*/
nml_extent eigenvaluesLanczos(	/* number of eigenvalues actually found	*/
    nml_dvector	 *value,	/*   out real eigenvalue vector		*/
    nml_ivector * convIterationNumbers ,
    nml_d3bands	 *T,		/*   out symmetric tridiagonal matrix	*/
    nml_extent	 *pIterations,	/*   out actual number of iterations	*/
    nml_dcvector *r_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
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
) 
{

    #ifdef MPI_TIMING
    MPI_TIME_TIC(start_time_lanc_all);
	#endif

    int eigenvalues = lanczos_tridiagonal(
        value, convIterationNumbers ,
        T, pIterations, r_n, q_n, q_n1, length, requested1, requested2,
        ConvCheckStartIter, ConvCheckSkipRate, imax,
        emin1, emax1, emin2, emax2, tolerance, resolution, matmul, argument,
        fp_trace, fp_log, verbose);

    #ifdef MPI_TIMING
    MPI_TOC(time_lanczos_iter_all, (start_time_lanc_all ));
    #endif

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

    MPI_Bcast(&eigenvalues, 1, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

    if (NULL != fp_log) {
        fprintf(fp_log, "CPU %d received Broadcast of eigenvalues\n", mpi_n3d_id);
        fflush(fp_log);
    }

    for(int i=0;i<eigenvalues;i++)
        printf("before rsort value[%d]=%1.15e  in iterations %d\n",i,value[i],convIterationNumbers[i]);

    /* Sort the eigenvalues in value in ascending order.			*/
    // rsort(value, eigenvalues);
    rsort_twovec(value, convIterationNumbers , eigenvalues) ;

    printf("after rsort\n");

    for(int i=0;i<eigenvalues;i++)
        printf("after rsort value[%d]=%1.15e  in iterations %d\n",i,value[i],convIterationNumbers[i]);

    return eigenvalues;
} /* eigenvaluesLanczos	*/


static int convMethodVector=3;

/*  Compute the eigenvectors corresponding to the eigenvalues found
    in a complex hermitian matrix using Lanczos' algorithm.		*/
nml_extent eigenvectorsLanczos(	/* number of eigenvalues actually found	*/
    nml_dcvector *cvalue,	/*   out complex eigenvalues		*/
    nml_dcmatrix *vector,	/*   out complex eigenvectors		*/
    const
    nml_dvector  *value,	/* in    real    eigenvalues		*/
    const
    nml_ivector * convIterationNumbers ,
    nml_extent	 eigenvalues,	/* in	 number of eigenvalues		*/
    nml_extent	 iterations,	/* in    number of iterations		*/
    const
    nml_d3bands	 *T,		/* in real symmetric tridiagonal matrix	*/
    nml_dcscalar *r_n,		/* inout r_{n} = q_{n+1}*beta_{n}	*/
    nml_dcscalar *q_n,		/* inout current  complex Lanczos vector*/
    nml_dcscalar *q_n1,		/* inout previous complex Lanczos vector*/
    nml_extent	 length,	/* in extent of vectors r_n, q_n & q_n1	*/
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    const int*	 argument[],	/* matrix-vector multiply argument list	*/
    nml_dscalar	 tolerance,	/* in eigenvalue convergence tolerance	*/
    FILE	 *fp_trace,	/* in   trace log text file pointer	*/
    FILE	 *fp_log,	/* in message log text file pointer	*/
    int		 verbose	/* in verbose diagnostic messages	*/
) 
{
    //MPI_TIME_TIC(start_time_lanc_all);

    nml_dmatrix	*xm = nml_dm_new(eigenvalues, iterations); /* real eigenvector matrix		*/
    nml_dscalar	beta_n1 = 0.0;	/* previous off diagonal element	*/
	
    #ifndef REDUCE_OPT
        nml_dscalar	beta2 = nml_dcv_ssqr(r_n);	/* <r_n | r_n>		*/

        Allreduce_MPI_sp(&beta2, &beta_n1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        beta_n1 = sqrt(beta_n1);
    #endif // REDUCE_OPT

    int* num_iteration;
    int flag;

    /* three Methods are implemented to operate on solving Ritz-vector
     * 	1. all iterate at the same number of iterations 
     * 	   for each Tridiagonal vector (old)
     * 	2. selective to iterate for each Tridiagonal vector(new)
     *  3. a method similar to 2, but using the information about the
     *     number of iterations used to get a converged eigenvalue
     *     This information is obtained in the eigenvalue computation
     *     and passed to this routine in the vector convIterationNumbers
     */

    if ( (convMethodVector==2) || (convMethodVector == 3) ){
        num_iteration = (int*)malloc(sizeof(int)*(eigenvalues));
    }

    if (mpi_n3d_masterid == mpi_n3d_id) {
        /* Compute the subspace eigenvectors on the master CPU only.
           Then distribute the vectors to the other CPUs			*/
        nml_dvector	*xx = nml_dv_new(iterations);
        nml_dvector	*yy = nml_dv_new(iterations);
        nml_d3bands	*tt = nml_d3_new(iterations);

        /* copy:  tt <-- T */
        /* nml_d3_3cpy(tt, T);
         * Since size of tt is less than size of T (imax),
         * nml_d3_3cpy() will correctly produce the warning.
         * Instead this specific case will use copyNmlD3bands()
         * which copies the two D3bands given offset starts and copy length. */

        nml_offset	iteration = 0;
        copyNmlD3bands(tt, 0, T, 0, iterations);

        /* Calculate eigenvectors of T for all converged eigenvalues	*/
        { 
            nml_offset	eigenvalue = 0;

            #ifdef MPI_TIMING
            MPI_TIC(start_time);
            double start_time_eigenvector=0.0;
            #endif

            for (eigenvalue = 0; eigenvalue < eigenvalues; ++eigenvalue) {
                iteration = 0;
                cvalue[eigenvalue].r = value[eigenvalue];

                if ( (convMethodVector == 1) || (convMethodVector == 2) )
                {
                    for (iteration = 0; iteration < iterations; ++iteration) 
                        yy[iteration] = 1.0*rand()/RAND_MAX;

                    originalInverseIteration(tt, xx, yy, value[eigenvalue]);
	    
                    /* Method 1*/
                    if(convMethodVector==1){
                        for (iteration = 0; iteration < iterations; ++iteration) {
                            xm[eigenvalue][iteration] = xx[iteration];
                    }}

                    /* Method 2*/
                    else if(convMethodVector==2){
                        num_iteration[eigenvalue]=iterations;
                        flag=0;
                        //nml_dscalar tolerance=1e-7;
                        int num_con_checking=5;

                        for (iteration = 0; iteration < iterations; ++iteration) {
                            if(iteration>0/*&&fabs(xx[iteration-1]-xx[iteration])<tolerance*/
                               && fabs(xx[iteration])<tolerance){
                                if(flag<num_con_checking){
                                    if(iteration==(unsigned int)(num_iteration[eigenvalue]+1)&&flag<num_con_checking){
                                        flag++;
                                    }
                                    else
	                                   flag=0;

                                    num_iteration[eigenvalue]=iteration;
                                }
                            }

                            xm[eigenvalue][iteration] = xx[iteration];
                        }

                        fprintf(fp_log,"eigenvalue= %d number of iteration =%d\n",eigenvalue,num_iteration[eigenvalue]); 
                        fflush(fp_log);
                    }  /* end method 2 */

                    copyNmlD3bands(tt, 0, T, 0, iterations);
                } /* end if method == 1 or method == 2 */
                else
                {   /* convMethodVector == 3 */
                    // printf("In eigenvectors: assuming method 3\n") ;
                    int requiredIterations ;
                    // save the number of iterations after which this eigenvalue converged
                    requiredIterations = convIterationNumbers[ eigenvalue ] ;
                    num_iteration[ eigenvalue ] = requiredIterations ;

                    // printf("For eigenvalue %d do %d iterations\n", eigenvalue , requiredIterations) ;
                    // zero the tridiagonal matrix tt and vectors
                    for (int indexj = 0 ; indexj < iterations ; indexj++)
                    {
                        tt[ -1 ][ indexj ] = 0.0 ;
                        tt[  0 ][ indexj ] = 0.0 ;
                        tt[  1 ][ indexj ] = 0.0 ;
	                    xx[ indexj ] = 0.0 ;
                        yy[ indexj ] = 0.0 ;
                    }

                    // copy the tridiagonal T into the tridiagonal tt, but only as many elements
                    // as the number of iterations necessary to converge this eigvalue
                    copyNmlD3bands( tt , 0 , T , 0 , requiredIterations ) ;
                    // now find the eigenvector corresponding to the eigenvalue

                    for (iteration = 0; iteration < requiredIterations; ++iteration) 
                        yy[iteration] = 1.0*rand()/RAND_MAX;

                    originalInverseIteration(tt, xx, yy, value[eigenvalue]);

                    // and copy it into the eigenvector matrix
	    
                    for (iteration = 0; iteration < requiredIterations; ++iteration) 
                        xm[eigenvalue][iteration] = xx[iteration];

                } /* end if method == 3 */
            }  /* end loop over eigenvalues */

            #ifdef Timing_detail
            MPI_TOC(start_time_eigenvector,start_time);
            printf("time for calculate eigenvectors of T = %.2f\n",start_time_eigenvector);
            #endif
        }

        /* Need to distribute the subspace eigenvectors to the various CPU's
           so that they can assemble the eigenvectors locally.
           The data is given in a real matrix xm
           with eigenvalues times iterations elements.			*/
        if (NULL != fp_log) {
            fprintf(fp_log, "CPU %d starting Broadcast of eigenvectors\n", mpi_n3d_id);
            fflush(fp_log);
        }

        MPI_Bcast(&xm[0][0], eigenvalues*iterations, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);

        if ((convMethodVector==2) || (convMethodVector == 3 ) )
            MPI_Bcast(&num_iteration[0], eigenvalues, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

        #ifdef Timing_detail
        printf("time for broadcast eigenvectors of T = %.2f\n",time_bcast);
        #endif

        if (NULL != fp_log) {
            fprintf(fp_log, "CPU %d finished Broadcast of eigenvectors\n", mpi_n3d_id);
            fflush(fp_log);
        }

        nml_dv_delete(xx);
        nml_dv_delete(yy);
        nml_d3_delete(tt);
    }
    else {			/* (mpi_n3d_masterid != mpi_n3d_id)		*/
        #ifdef	MPI3d
        /* Receive the subspace eigenvectors from the master CPU.		*/
        if (NULL != fp_log) {
            fprintf(fp_log, "CPU %d expected Broadcast of eigenvectors\n", mpi_n3d_id);
            fflush(fp_log);
        }

        MPI_Bcast(&xm[0][0], eigenvalues*iterations, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);

        if ((convMethodVector==2) || (convMethodVector == 3) )
            MPI_Bcast(&num_iteration[0], eigenvalues, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

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
            fprintf(fp_log, "Repeat Lanczos Iterations and construct eigenvectors\n" "iteration:\n");
            fflush(fp_log);
        }
    }

    { 
        nml_offset iteration = 0;

        #ifdef MPI_TIMING
        double start_time_iteration_eigenvector=0.0;              
        time_lanczos_iter=0.0;
        MPI_TIC(start_time3); 
        #endif

        #ifdef MPI_TIMING
        time_lanczos_first=0.0;
        time_lanczos_rest=0.0;
        #endif
		
        for (iteration = 0; iteration < iterations; ++iteration) {

            #ifdef REDUCE_OPT
            static nml_dscalar alpha = 0.0;
            static nml_dscalar beta  = 0.0;
		    #else
            nml_dscalar alpha = 0.0;
            nml_dscalar beta  = 0.0;
            #endif // REDUCE_OPT

            if ((mpi_n3d_id == mpi_n3d_masterid) && (0 == iteration%100)) {
                if (NULL != fp_log) {
                    fprintf(fp_log, "%d of %d max\n", iteration,iterations);
                    fflush(fp_log);
            }}

            #ifdef MPI_TIMING
            MPI_TIC(start_time2);
            #endif

            #ifdef MPI_TIMING
            MPI_TIC(start_time2);
            #endif

            single_lanczos_iteration(matmul, argument, 
                #ifdef REDUCE_OPT
                iteration,
                #endif // REDUCE_OPT
                r_n, q_n, q_n1, length, &alpha, &beta, beta_n1);

            #ifdef MPI_TIMING
            MPI_TOC(time_lanczos_iter, start_time2);
            #endif

            #ifdef MPI_TIMING
            if ( iterations == 0 )
                MPI_TOC(time_lanczos_first, start_time2);
            else
                MPI_TOC(time_lanczos_rest, start_time2);
            #endif

            beta_n1 = beta;
      
            { 
                nml_offset eigenvalue = 0;

                for (eigenvalue = 0; eigenvalue < eigenvalues; ++eigenvalue) {
                    nml_offset	j = 0;
                    if ((convMethodVector==2) || (convMethodVector == 3)) {
                        if(iteration>(unsigned int)(num_iteration[eigenvalue]))
                            continue;
                    }

                    for (j = 0; j < length; ++j) {
                        vector[eigenvalue][j].r += xm[eigenvalue][iteration]*q_n[j].r;
                        vector[eigenvalue][j].i += xm[eigenvalue][iteration]*q_n[j].i;
                    }
                }
            }

            MPI_Barrier(MPI_COMM_WORLD);
        }

        #ifdef Timing_detail
        MPI_TOC(start_time_iteration_eigenvector,start_time3);
        fprintf(fp_log,"time for calculate full eigenvectors = %.2f\n",start_time_iteration_eigenvector);
        fprintf(fp_log,"time for lanczos iteration = %.2f\n",time_lanczos_iter);
        fflush(fp_log);
        #endif
    }

    if (mpi_n3d_id == mpi_n3d_masterid) {
        if (NULL != fp_log) {
            fprintf(fp_log, "Done with Iterations and construction of eigenvectors\n");
            fflush(fp_log);
    }}

    if ( (convMethodVector==2) || (convMethodVector==3))
        free(num_iteration);

    nml_dm_delete(xm);

    return eigenvalues;
} /* eigenvectorsLanczos */


void saveTridiagonalMatrix(
    const nml_d3bands	*T,	/* in    real symmetric tridiagonal matrix	*/
    int iterations,	/* in    Lanczos iterations			*/
    int processors,	/* in    number of MPI processors		*/
    FILE *fp_tridiag	/* inout tridiagonal matrix binary file pointer	*/
) 
{			/*    save real symmetric tridiagonal matrix T	*/
    nml_dvector	*T_write = nml_dv_new(2*iterations);
    int	iteration = 0;

    for (iteration = 0; iteration < iterations; ++iteration)
        T_write[2*iteration  ] = T[0][iteration  ];
  
    for (iteration = 1; iteration < iterations; ++iteration)
        T_write[2*iteration-1] = T[1][iteration-1];
  
    fwrite(&iterations, sizeof(int), 1, fp_tridiag);
    fwrite(&processors, sizeof(int), 1, fp_tridiag);
    fwrite(&T_write[0], sizeof(nml_dscalar), 2*iterations, fp_tridiag);

    nml_dv_delete(T_write);
}	/* saveTridiagonalMatrix	*/


/* write out the Tridiagonal matrix in ascii with file name 'tridiag.m',
 * use X=[array] format, the array has two columns, the first column
 * is alpha value, the second column is beta value.
 * It can be used in matlab to do further analysis. It was written on 09/21/2004 
 * The old function saveTridiagonalMatrix only write out the matrix 
 * in binary format, and can be restored by restoreTridiagonalMatrix in NEMO 3-D
 */
void saveTridiagonalMatrixInAscii(
    const
    nml_d3bands	*T,	/* in    real symmetric tridiagonal matrix	*/
    int iterations,	/* in    Lanczos iterations			*/
    int processors,	/* in    number of MPI processors		*/
    FILE *fp_tridiag	/* inout tridiagonal matrix binary file pointer	*/
) 
{			/*    save real symmetric tridiagonal matrix T	*/
    int	iteration = 0;
    double zero=0.0;
    FILE* temp_tridiag=fopen("tridiag.m","w");

    fprintf(temp_tridiag,"N=%d\n",iterations);
    fprintf(temp_tridiag,"X=[\n");

    for (iteration =0;iteration<iterations; ++iteration)
    {
       if(iteration==0)
            fprintf(temp_tridiag, "%25.18e   %25.18e\n",zero,T[0][iteration]);
        else 
            fprintf(temp_tridiag, "%25.18e   %25.18e\n",T[1][iteration-1],T[0][iteration]);
    }

    fprintf(temp_tridiag,"];\n");
    fclose(temp_tridiag);
}	/* saveTridiagonalMatrixInAscii	*/


void restoreTridiagonalMatrix(
    nml_d3bands	**pT,	/*   out real symmetric tridiagonal matrix	*/
    nml_extent
	*pIterations,	/*   out Lanczos iterations pointer		*/
    nml_extent
	*pProcessors,	/* inout number of MPI processors		*/
    FILE *fp_tridiag	/* inout tridiagonal matrix binary file pointer	*/
) 
{			/* restore real symmetric tridiagonal matrix *pT*/
    fread(pIterations, sizeof(nml_extent), 1, fp_tridiag);
    fread(pProcessors, sizeof(nml_extent), 1, fp_tridiag);

    { 
        nml_extent	iterations = *pIterations;
        iterations = *pIterations=iterations-200;
        nml_dvector	*T_read = nml_dv_new(2*iterations);

        fread(&T_read[0], sizeof(nml_dscalar), 2*iterations, fp_tridiag);
        *pT = nml_d3_new(iterations);

        { 
            nml_offset iteration = 0;
            for (iteration = 0; iteration < iterations; ++iteration)
                (*pT)[0][iteration] = T_read[2*iteration];
 
            for (iteration = 1; iteration < iterations; ++iteration)
                (*pT)[1][iteration-1] = (*pT)[-1][iteration] = T_read[2*iteration-1];
        }

        nml_dv_delete(T_read);
    }
}	/* restoreTridiagonalMatrix	*/


void saveEigenvalues(
    const
    nml_dvector	*value,	/* in    real eigenvalue vector			*/
    int eigenvalues,	/* in    number of eigenvalues			*/
    FILE *fp_eigval	/* inout real eigenvalue binary file pointer	*/
) 
{			
    /* save real eigenvalue vector value		*/
    fwrite(&eigenvalues, sizeof(int), 1, fp_eigval);
    fwrite(value, sizeof(nml_dscalar), eigenvalues, fp_eigval);
}	/* saveEigenvalues		*/


void restoreEigenvalues(
    nml_dvector	**pValue,/*  out real eigenvalue vector			*/
    nml_extent
	*pEigenvalues,	/*   out number of eigenvalues			*/
    FILE *fp_eigval	/* inout real eigenvalue binary file pointer	*/
) 
{			
    /* restore real eigenvalue vector *pValue	*/
    fread(pEigenvalues, sizeof(nml_extent), 1, fp_eigval);
    *pValue = nml_dv_new(*pEigenvalues);
    fread(*pValue, sizeof(nml_dscalar), *pEigenvalues, fp_eigval);
}	/* restoreEigenvalues		*/


/* Deprecated -- Calculates requested eigenvalues and eigenvectors
 * of the Hamiltonian in the energy range [emin:emax].
 * imax is the maximum number of iterations + 1.			
 */
/* The control parameter that controls Tridiagonal matrix is written out
 * in ascii format
 */
static int writeOutAscii=0;

/*********************************************************************************************/

int eigensolve_lanczos(		/* number of eigenvalues actually found	*/
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
    const char	 *fn_tridiag,	/* in tridiagonal matrix binary file name */
    const char	 *fn_eigval,	/* in eigenvalue  vector binary file name */
    nml_dcvector *r_n,          /* inout r_{n} = q_{n+1}*beta_{n}       */
    const nml_dcvector *r_0,    /* in    initial guess vector           */
    nml_dcvector *q_n,          /* inout current  complex Lanczos vector*/
    nml_dcvector *q_n1,         /* inout previous complex Lanczos vector*/
    nml_dcvector **pCValue,     /*   out complex eigenvalue  vector     */
    nml_dcmatrix **pCVector,    /*   out complex eigenvector matrix     */
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
    int		 restoreState	/* in read tridiagonal matrix from file	*/
) 
{
    nml_extent	iterations = 0;
    nml_extent	eigenvalues = 0;
    nml_d3bands	*T = NULL;
    FILE		*fp_log = stdout;
    FILE		*fp_trace = (NULL == fn_trace)? NULL: fopen(fn_trace, "a");
    nml_ivector * convIterationNumbers ;

    convIterationNumbers = nml_iv_new( imax ) ;

    if (convIterationNumbers == NULL ) printf("ATTENTION!!!! ZERO POINTER!!!!") ;
    //  nml_iv_fill( convIterationNumbers , 0 ) ;
  
    #ifdef REDUCE_OPT
    nml_dcvector *q_n_tmp = nml_dcv_new(length);

    for(int i = 0; i < length; i++)
    {
		q_n[i].r = 1.0 * rand()/RAND_MAX;
		q_n[i].i = 0.0;
    }

	// q_1 <-- q_1 / norm(q_1);
	{
	    nml_dscalar q_norm = 0.0;

		for(int i = 0; i < length; i++)
		{
			nml_dscalar *pq_n = &(q_n[i].r);
			
			q_norm += (*pq_n)*(*pq_n);
		    ++pq_n;
			q_norm += (*pq_n)*(*pq_n);
		}

		{
		    /* Assemble one beta. Sum them and then take sq root.		*/
		    nml_dscalar	qn = 0.0;

		    Allreduce_MPI_sp(&q_norm, &qn, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

		    qn = sqrt(qn);

			for(int i = 0; i < length; i++)
			{
			    nml_dscalar *pq_n = &(q_n[i].r);

				*pq_n /= qn;
				q_n_tmp[i].r = *pq_n;
				pq_n++;
				*pq_n /= qn;
				q_n_tmp[i].i = *pq_n;
			}
		}
	}
    #endif // REDUCE_OPT

    if ((NULL != fn_trace) && (NULL == fp_trace)) 
        error("Could not open the file \"%s\"", fn_trace);

    /* Compute the real symmetric tridiagonal matrix and eigenvalues
     * or restore them from binary files saved in a previous computation.	
	 */

    if (restoreState) 
    { 
        /* restore real symmetric tridiagonal matrix T from file		*/
        nml_extent	processors;
        FILE *fp_tridiag = fopen(fn_tridiag, "r");
      
        if (NULL == fp_tridiag) 
        {
            fprintf(stderr, "Could not open file %s\n", fn_tridiag);
            exit(EXIT_FAILURE);
            MPI_Finalize();		/* PLEASE EXPLAIN ME			*/
        }

        restoreTridiagonalMatrix(&T, &iterations, &processors, fp_tridiag);
      
        if (processors != (nml_extent)mpi_n3d_numprocs) 
            if (NULL != fp_log) 
                fprintf(fp_log, "ACHTUNG!!!:  you are running on a different number of processors\n"
                                "from the run that generated your tridiagonal matrix\n");

        fclose(fp_tridiag);

        /* restore real eigenvalue vector *pValue	from file		*/
        FILE* fp_eigval = fopen(fn_eigval, "r");
        if (NULL == fp_eigval) 
        {
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
    else 
    { 
        /* (!restoreState) */
        /* Compute the tridiagonal matrix and eigenvalues from scratch.	*/
 
        *pValue	= nml_dv_new(imax);	/* actual eigenvalues found	*/
        T = nml_d3_new(imax);	/* real symmetric tridiagonal matrix	*/

        #ifdef Timing_detail
            double eigenvalueTime=0.0;
            MPI_TIC(time_lanczos_setup_once);
        #endif

        eigenvalues = eigenvaluesLanczos(
                      *pValue,  convIterationNumbers, T, &iterations,
                      r_n, q_n, q_n1, length, requested1, requested2,
                      ConvCheckStartIter, ConvCheckSkipRate, imax,
                      emin1, emax1, emin2, emax2, tolerance, 100.0*tolerance, matmul, 
                      (const int**)argument, fp_trace, fp_log, verbose);

        if (mpi_n3d_masterid == mpi_n3d_id) 
        {
            for (int indexj = 0 ; indexj < eigenvalues ; indexj++ )
            {
                printf("Eigenvalue: %.8e  in %d iterations\n",
                (*pValue)[ indexj ] , convIterationNumbers[ indexj ] ) ;
            }
        }

        #ifdef Timing_detail
            MPI_TOC(eigenvalueTime,time_lanczos_setup_once);
            printf("eigenvalueLanczos running time = %.2f\n",eigenvalueTime);
        #endif
      
        if (mpi_n3d_masterid == mpi_n3d_id) {
            /* Only the master processor
             * saves the tridiagonal matrix to a machine readable file and
             * saves the eigenvalues to a machine readable file.		
			 */

            /* Save state = {tridiagonal matrix, eigenvalues}.		*/
            { /* save real symmetric tridiagonal matrix T to file		*/
                FILE*	fp_tridiag = fopen(fn_tridiag, "w");

                if (NULL != fp_tridiag) {
                    saveTridiagonalMatrix(T, iterations, mpi_n3d_numprocs, fp_tridiag);

                    if(writeOutAscii==1)
                        saveTridiagonalMatrixInAscii(T, iterations, mpi_n3d_numprocs, fp_tridiag);

                    fclose(fp_tridiag);
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

        //    iterations=iterations-200;
    
    }   /* end the if(restoreState) condition  */

    printf("before eigenvectr %d\n",eigenvalues);

    /* At this point,
       the real symmetric tridiagonal matrix and eigenvalues have been
       computed and saved to files or restored from files so that
       it is now possible to compute the corresponding eigenvectors.	
    */
  
    if (NULL != r_0) 
    {
        nml_memory_report("before the nml_dcm_new"); 

        /* Set up the initial guess for the rerun of the Lanczos algorithm.	*/

        /* The two operations below are not necessary; explicit loops later on... */
        //      nml_dcv_vcpy(r_n, r_0);			
        //  nml_dcv_fill(q_n1, nml_dcmplx(0.0, 0.0));	

        nml_dcscalar zero = nml_dcmplx(0.0, 0.0);
        nml_offset	j = 0;

        #ifdef REDUCE_OPT
	        for( j = 0; j < length; j++)
		    {
		        q_n[j].r = q_n_tmp[j].r;
			    q_n[j].i = q_n_tmp[j].i;
			    r_n[j] = zero;
			    q_n1[j] = zero;
            }

            nml_dcv_delete(q_n_tmp);
        #else
            for (j = 0; j < length; ++j) 
	        {
                r_n[j] = r_0[j];
                q_n1[j] = zero;
            }
        #endif // REDUCE_OPT

        *pCValue  = nml_dcv_new(eigenvalues);
        *pCVector = nml_dcm_new(eigenvalues, length);

        nml_memory_report("after the nml_dcm_new"); 

        #ifdef Timing_detail
            MPI_TIC(time_lanczos_setup_once);
        #endif

        eigenvectorsLanczos(*pCValue, *pCVector, *pValue, convIterationNumbers,  eigenvalues,
                            iterations, T, r_n, q_n, q_n1, length, matmul, (const int**)argument,
                            tolerance,fp_trace, fp_log, verbose);
    }   // end condition if (NULL != r_0)

    #ifdef Timing_detail
        double eigenvectorTime=0.0;
        MPI_TOC(eigenvectorTime,time_lanczos_setup_once);
        printf("eigenvectorLanczos running time = %.2f\n",eigenvectorTime);
    #endif

    if (fp_trace)   fclose(fp_trace);

    nml_d3_delete(T);
    nml_iv_delete( convIterationNumbers ) ; 

    #ifdef MPI_TIMING
	    printf("[LANCZOS] TOTAL TIME SPENT ON ALLREDUCE = %.2f\n", time_lanczos_allreduce_tot);
	#endif 

    return eigenvalues;
}     /* eigensolve_lanczos	*/

