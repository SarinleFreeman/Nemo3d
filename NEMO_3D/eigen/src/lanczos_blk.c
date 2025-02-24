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
$Source: /repo/eigen/src/lanczos_blk.c,v $
*****************************************************************************/

//Written by Rajib Rahman, Faisal Saied, Gerhard Klimeck


#include "lanczos_blk.h"
#ifdef INTEL_MODIFIED
	#include <memory.h>
	#include <emmintrin.h>
#endif //INTEL_MODIFIED
#include <vector>
using std::vector;

#include <iostream>
using namespace std;

#define ARPACK_EIGENVECTOR_SOLVER

#ifdef  __cplusplus
extern "C" {
#endif

    #ifdef FORTRAN_UNDERSCORE
    #define zhbevx		  zhbevx_

    #ifdef ARPACK_EIGENVECTOR_SOLVER
	#define zgbtrf        zgbtrf_
    #define znaupd        znaupd_
    #define zneupd       zneupd_
    #define zhesv         zhesv_
	#define zgbtrs        zgbtrs_
    #endif

    #endif

    #ifdef FORTRAN_UNDERSCORE2
    #define zhbevx            zhbevx__

    #ifdef ARPACK_EIGENVECTOR_SOLVER
	#define zgbtrf        zgbtrf__
    #define znaupd        znaupd__
    #define zneupd        zneupd__
    #define zhesv         zhesv__
	#define zgbtrs        zgbtrs__
    #endif

    #endif

    #ifdef FORTRAN_ALLCAPS
    #define zhbevx            ZHBEVX

    #ifdef ARPACK_EIGENVECTOR_SOLVER
	#define zgbtrf        ZGBTRF
    #define znaupd       ZNAUPD
    #define zneupd       ZNEUPD
    #define zhesv         ZHESV
	#define zgbtrs        ZGBTRS
    #endif

    #endif

    extern void zhbevx(char *JOBZ, char *RANGE, char *UPLO, int *N, int *KD, struct complex *AB, int *LDAB, struct complex *Q, int *LDQ,
    double *VL, double *VU, int *IL, int *IU, double *ABSTOL, int *M, double *W, struct complex *Z, int *LDZ, struct complex *WORK,
    double *RWORK, int *IWORK, int *IFAIL, int *INFO);

    #ifdef ARPACK_EIGENVECTOR_SOLVER
    extern void znaupd(int*, char*, int*, char*, int*, double*, complex*, int*, complex*,
                int*, int*, int*, complex*, complex*, int*, double*, int*);
    extern void zneupd( int*, char*, int*, complex*, complex*, int*, complex*,
                 complex*, char*, int*, char*, int*, double*, complex*, int*,
                 complex*, int*, int*, int*, complex*, complex*, int*, double*, int*);
	extern void zhesv(char*, int*, int*, complex*,int*, int*, complex*, int*, complex *, int *, int *);
	extern void zgbtrf(int *, int *, int *, int *, complex *, int *, int *, int *); 
	extern void zgbtrs(char *, int *, int *, int *, int *, complex *, int *, int *, complex *, int *, int *);
    #endif

#ifdef  __cplusplus
}
#endif

struct complex {double re; double im;};

#ifdef MPI_TIMING
MPI_TIME_INIT(time_bl_QBqr);
MPI_TIME_INIT(time_bl_RAQ);
MPI_TIME_INIT(time_bl_Anext);
MPI_TIME_INIT(time_bl_Rnext);
MPI_TIME_INIT(time_bl_Qnext);
MPI_TIME_INIT(time_bl_convec);
MPI_TIME_INIT(time_bl_conval);
MPI_TIME_INIT(Tstart);
#endif
//test printing
void print_vector(vector< nml_dcscalar* > *q, int m, int p){
    nml_offset j=0;
    for(int i=0;i<p;i++){
        nml_dcvector *q1= (*q)[i];
        cout<<"\nvec:"<<i<<endl;
        for(j=0;j<m;j++){ 
            nml_dscalar *r1_n_j = &(q1[j].r);
            cout<<*r1_n_j<<"   "<<*(++r1_n_j)<<endl;
}}}

//test printing
void print_Matrix( nml_dcmatrix *a_ptr, int p, int q){
    //if (mpi_n3d_id == mpi_n3d_masterid) {
    for(int i=0;i<p;i++){
        for(int j=0;j<q;j++){
            cout<<a_ptr[i][j].r<<" "<<a_ptr[i][j].i<<"      ";
        }
        cout<<endl;
}}

//compute A'*B, for two vectors A and B
void dot_product(
    nml_dcvector *A,
    nml_dcvector *B,
    nml_dcscalar *C_ptr,
    nml_extent m
)
{
	(*C_ptr).r=0.0;
	(*C_ptr).i=0.0;

	nml_offset    j = 0;
    nml_dcscalar  alc = {0.0, 0.0};

	//needs to be improved
    for (j = 0; j < m; ++j) {
     	nml_dscalar *a = &(A[j].r);
        nml_dscalar *b = &(B[j].r);

        alc.r += a[0]*b[0] + a[1]*b[1];
        alc.i += a[0]*b[1] - a[1]*b[0];
    }
               
    Allreduce_MPI_sp(&alc, C_ptr, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(&alc, C_ptr, 2, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
}


//compute A'*B, for two vectors A and B
void mult_dot_product_hermitian(
    nml_dcmatrix *A,
    nml_dcmatrix *B,
    nml_dcvector *C,
	nml_extent block_size,
    nml_extent length
)
{
    static nml_dcvector *alc = NULL;
    int indx = 0;

    if(alc == NULL)
	{
	    alc = nml_dcv_new(block_size*(block_size+1)/2);
	}

	//needs to be improved
	for(int i = 0; i < block_size; i++)
	{
	    nml_dcvector *a = (nml_dcvector *)&A[i][0];

	    for(int j = 0; j <= i; j++)
		{
	        alc[indx].r = 0; alc[indx].i = 0;

		    nml_dcvector *b = (nml_dcvector *)&B[j][0];

		    for(int k = 0; k < length; k++)
			{
                alc[indx].r += a[k].r*b[k].r + a[k].i*b[k].i;
                alc[indx].i += a[k].r*b[k].i - a[k].i*b[k].r;
		    }

			indx++;
		}
    }
               
    Allreduce_MPI_sp(alc, C, block_size*(block_size+1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    //MPI_Allreduce(alc, C, block_size*(block_size+1), MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);

	//nml_dcv_delete(alc);
}


//returns -1 if vectors are not to be deflated
//returns index of the vector to be deflated i<p

int deflation_test( 
	vector < nml_dcvector* > *Q,
	nml_dcmatrix *Beta,
	nml_extent m, 
	int *p,
	int iteration,
	double deflation_tolerance ){

	int deflate=-1;

	for(int i=0;i<*p;i++){
		if(Beta[i][i].r<deflation_tolerance){
			deflate=i;
			break;
		}
	}

	//cout<<"Checking deflation"<<endl;
	//print_Matrix(Beta, *p, *p);

	if((deflate>-1)&&(deflate<*p)){
		cout<<"Deflation test positive: Removing vector "<<deflate<<endl;
		nml_dcvector *q1= (*Q)[deflate];
		(*Q).erase(((*Q).begin())+deflate);
		nml_dcv_delete(q1);
		//cout<<"Size="<<(*Q).size()<<endl;
		for(int k=deflate;k<(*p)-1;k++){
			for(int j=k;j<*p;j++){
				Beta[k][j].r=Beta[k+1][j].r;
				Beta[k][j].i=Beta[k+1][j].i;
				Beta[k+1][j].r=0.0;
				Beta[k+1][j].i=0.0;
			}
		}
		if(deflate==(*p)-1){
		    	Beta[deflate][deflate].r=0.0;
                    	Beta[deflate][deflate].i=0.0;
		}
	//print_Matrix(Beta,*p,*p);
	*p=*p-1;
	}

    return deflate;
}



//[Q, B]=blqr(R)
void block_Lanczos_QR( 
	vector< nml_dcvector* > *R, 
	vector< nml_dcvector* > *Q, 
	nml_dcmatrix *Beta, 
	nml_extent m, 
	int p
)
{
    //r=R(:,1)
    nml_dcvector *r1= (*R)[0];

    //B(1,1)=norm(r1)
    nml_dscalar g1_n1=nml_dcv_ssqr(r1);
    nml_dscalar g1=0.0;
    Allreduce_MPI_sp(&g1_n1, &g1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    g1 = sqrt(g1);
    Beta[0][0].r=g1;

    //Q(:,1)=r/B(1,1)
    nml_offset    j = 0;
    nml_dscalar inv_g_n1 = 1.0/g1;
    nml_dcvector *q1= (*Q)[0];
    for (j = 0; j < m; j++) {	   
        nml_dscalar *r1_n_j = &(r1[j].r);
        nml_dscalar *q1_n_j = &(q1[j].r);
        *q1_n_j = *r1_n_j*inv_g_n1;
        ++q1_n_j;
        ++r1_n_j;
        *q1_n_j = *r1_n_j*inv_g_n1;
    }

    nml_dcscalar b1={0.0,0.0};
    nml_dcvector *r2=nml_dcv_new(m);

    for(int i=1;i<p;i++){
        //r=R(:,i)
        nml_dcvector *r1= (*R)[i];

        //nml_dcvector *r2=nml_dcv_new(m);
        for(int l=0;l<m;l++){
            r2[l].r=r1[l].r;
            r2[l].i=r1[l].i;
        }

        for(int k=0;k<i;k++){
	        nml_dcvector *q1= (*Q)[k];
	        dot_product(q1,r2,&b1,m);
	        Beta[k][i]=b1;

	        for (j = 0; j < m; ++j) { 
    		    nml_dscalar *r2_n_j = &(r2[j].r);
    		    nml_dscalar *q1_n_j = &(q1[j].r);

    		    *r2_n_j=(*r2_n_j)-b1.r*(*q1_n_j)+b1.i*(q1_n_j[1]);
    		    ++r2_n_j;
    		    *r2_n_j=(*r2_n_j)-b1.r*(q1_n_j[1])-b1.i*(*q1_n_j);
    	    }//end of j	
        }//end of k  

        g1_n1=nml_dcv_ssqr(r2);
        g1 = 0.0;

        Allreduce_MPI_sp(&g1_n1, &g1, 1, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
        g1 = sqrt(g1);
        Beta[i][i].r=g1;
        inv_g_n1 = 1.0/g1;

        nml_dcvector *q1= (*Q)[i];

        for (j = 0; j < m; ++j) {
	        // first do the real part ... 
    	    nml_dscalar *q1_n_j = &(q1[j].r);
	        nml_dscalar *r2_n_j = &(r2[j].r);

    	    *q1_n_j = *r2_n_j*inv_g_n1;

    	    /* now do the imaginary part */
    	    ++q1_n_j;
	        ++r2_n_j;
    	    *q1_n_j = *r2_n_j*inv_g_n1;
        }//end of j
    }//end of i	

    nml_dcv_delete(r2);
} //end of blqr


//indx=1 => orthogonal
//indx=0 => non-orthogonal
void orthogonality_test(
    vector< nml_dcvector* > *R,
    nml_extent m,
    int p,
    int *indx
)
{
    nml_dcscalar b1={0.0,0.0};

    for(int i=0;i<p;i++)
	    indx[i]=1;

    for(int i=1;i<p;i++){
	    nml_dcvector *r1= (*R)[i];
	    int flag=1;

	    for(int j=0;j<i;j++){
		    nml_dcvector *r2= (*R)[j];
		    dot_product(r1,r2,&b1,m);
		    if((b1.r>1e-14) || (b1.i>1e-14)){
			    flag=0;
			    break;
		    }
	    }//end of j

	    indx[i]=flag;
    }//end of i
}


void single_iteration(

    #ifndef BLANCZ_MATMULT 
    void (*mul)(const int**, nml_dcscalar*, const nml_dcscalar*),
	#else
    void (*mul)(const int**, nml_dcmatrix*, const nml_dcmatrix*),
	#endif

    const int* argument[],              // in    matrix-vector multipier argument list
    vector <nml_dcvector*> *R_n,        // inout r1_{n-1} = q1_n*beta_{n-1}
    vector <nml_dcvector*> *Q_n,        //   Q_current
    vector <nml_dcvector*> *Q_n_1,      //   Q_old
    nml_extent   m,                     // in    extent of vectors r_n, q_n and q_{n-1}
    nml_extent   iteration,             // in    the current iteration
    int p,                              // block size
    int block_size,	
    nml_dcmatrix  *Alpha,               //   out real  on-diagonal element alpha_{n}
    nml_dcmatrix  *Beta,                 //   out real off-diagonal element beta_{n}
    int *def,
    double deflation_tolerance,
    int deflation_on		
) 
{
    //cout<<"beginning single iteration:"<<iteration<<endl;

    nml_offset j = 0;
    int deflate=-1;
    int new_p=p;

    for(int i=0;i<block_size;i++)
        for(int k=0;k<block_size;k++){
            Beta[i][k].r=0.0;
            Beta[i][k].i=0.0;
		    Alpha[i][k].r=0.0;
            Alpha[i][k].i=0.0;
        }

    if(iteration!=0){
	    //[Q,B]=blqr(R)
	    block_Lanczos_QR(R_n,Q_n,Beta,m,p);


        //deflation test
        if((p>1)&&(deflation_on!=0)){
	        deflate=deflation_test(Q_n, Beta, m, &new_p, iteration, deflation_tolerance);
	        *def=deflate;
        }
    }

    #ifndef BLANCZ_MATMULT
        //R=A*Q
        for(int i=0;i<new_p;i++){
            nml_dcvector *r1_n= (*R_n)[i];
            nml_dcvector *q1_n= (*Q_n)[i];

            //r1=A*q1
            mul(argument, r1_n, q1_n);
        }
	#else
	    mul(argument, (nml_dcmatrix *)&R_n[0][0], (nml_dcmatrix *)&Q_n[0][0]);
	#endif

    if(iteration!=0){
        //R = R - Qold * Bj';
        for(int i=0;i<new_p;i++){
	        nml_dcvector *r1_n=(*R_n)[i];

            for(j=0;j<m;j++){
                nml_dcscalar t1={0.0, 0.0};
                for(int k=i;k<p;k++){
                    nml_dcscalar q1_n= (*Q_n_1)[k][j];
                    nml_dcscalar b_n=Beta[i][k];
                    t1.r+=(q1_n).r*b_n.r+(q1_n).i*b_n.i;
                    t1.i+=-(q1_n).r*b_n.i+(q1_n).i*b_n.r;
                }//end of k
                r1_n[j].r=r1_n[j].r-t1.r;
                r1_n[j].i=r1_n[j].i-t1.i;
            }//end of j
        }//end of i	
    } //end of iteration !=0

    // Aj = Q'*R;
	// Note: This place seems to be the only place to reduce MPI_Allreduce calls
	//#define BLK_LANCZ_ALLREDUCE
	#ifndef BLK_LANCZ_ALLREDUCE 
    for(int i=0;i<new_p;i++){
        nml_dcvector *q1_n= (*Q_n)[i];
        for(int j=0;j<=i;j++){
            nml_dcvector *r1_n= (*R_n)[j];
            nml_dcscalar t1={0.0, 0.0};
            dot_product(q1_n,r1_n,&t1,m);
            Alpha[i][j].r=Alpha[j][i].r=t1.r;
            if(i==j)
			    Alpha[i][j].i=0.0;
		    else{
			    Alpha[i][j].i=t1.i;
                Alpha[j][i].i=-t1.i;
		    }
        }// end of j
    }//end of i
	#else
	{
	    nml_dcmatrix *q1_n = (nml_dcmatrix *)&Q_n[0][0];
		nml_dcmatrix *r1_n = (nml_dcmatrix *)&R_n[0][0];
		static nml_dcvector *t1 = NULL; //= nml_dcv_new(new_p * (new_p + 1) / 2);
		int indx = 0;

		if(t1 == NULL)
		    t1 = nml_dcv_new(new_p * (new_p + 1) / 2);

		mult_dot_product_hermitian(q1_n, r1_n, t1, new_p, m);

		for(int i = 0; i < new_p; i++)
		{
		    for(int j = 0; j <= i; j++)
			{
			    Alpha[i][j].r = Alpha[j][i].r = t1[indx].r;
				if(i == j) Alpha[i][j].i = 0;
				else
				{
				    Alpha[i][j].i = t1[indx].i;
					Alpha[j][i].i = -t1[indx].i;
				}
				indx++;
			}
		}
	}
	#endif

    //R = R - Q*Aj;
    for(int i=0;i<new_p;i++){
	    nml_dcvector *r1_n=(*R_n)[i];
        for(j=0;j<m;j++){
            nml_dcscalar t1={0.0, 0.0};
            for(int k=0;k<new_p;k++){
                nml_dcscalar q1_n= (*Q_n)[k][j];
                nml_dcscalar a_n=Alpha[k][i];
                t1.r+=(q1_n).r*a_n.r-(q1_n).i*a_n.i;
                t1.i+=(q1_n).r*a_n.i+(q1_n).i*a_n.r;
            }//end of k
			
            r1_n[j].r=r1_n[j].r-t1.r;
            r1_n[j].i=r1_n[j].i-t1.i;
        }//end of j
    }//end of i

    //Q_old=Q
    for(int i=0;i<new_p;i++)
        for(j=0;j<m;j++)
            (*Q_n_1)[i][j]=(*Q_n)[i][j];

    if(deflate!=-1){
	    nml_dcvector *q1= (*Q_n_1)[p-1];
        (*Q_n_1).erase(((*Q_n_1).begin())+p-1);
        nml_dcv_delete(q1);
    }

    //cout<<"exiting single iteration"<<endl;
} //end of single iteration


void find_block_converged_eigvalues_Method4(
    nml_dvector *new_eig,     // in/out actual converged eigenvalues
    int 	*n_new,
    nml_dvector *old_eig,
    int 	*n_old,	
    nml_dvector *converged_eig,
    int         *found,             // out number of converged eigenvalues
    int 	requested,	
    int         *n_nosuccess,   // out number of no success
    FILE        *fp_trace,      // in trace log text file pointer
    int         verbose,        // in verbose diagnostic messages
    nml_dscalar tolerance,      // in eigenvalue convergence tolerance
    FILE        *fp_log         // in message log text file pointer
) 
{
    double resolution=100*tolerance;
    double *converge_temp=new double[requested];
    int num_converged=*found;

    cout<<"convergence checking M4"<<endl;

    #if 0
    cout<<"Old eigenvalues: "<<*n_old<<endl;
    for(int i=0;i<*n_old;i++){
        cout<<"eig:"<<i<<"       "<<old_eig[i]<<endl;
    }
    
    cout<<"New eigenvalues: "<<*n_new<<endl;
    for(int i=0;i<*n_new;i++){
        cout<<"eig:"<<i<<"       "<<new_eig[i]<<endl;
    }
	#endif // DEBUG
  
    (*found)=0;
    for(int i=0;i<*n_new;i++){
        for(int j=0;j<*n_old;j++){
            if(fabs(new_eig[i]-old_eig[j])<tolerance){
		        converge_temp[*found]=new_eig[i];
		        old_eig[j]=500;
		        (*found)++;
		        break;                
	    }}

  	    if((*found)==requested)
		    break;
    }
  
 
    for(int i=0;i<num_converged;i++){
        for(int j=0;j<*found;j++){
            if(fabs(converged_eig[i]-converge_temp[j])<resolution){
		        converged_eig[i]=converge_temp[j];
                converge_temp[j]=500;
                break;
    }}}

    for(int i=0;i<*found;i++){
        if(fabs(converge_temp[i])<490){
            converged_eig[num_converged]=converge_temp[i];
		    num_converged++;
		    *n_nosuccess=0;	

		    if(num_converged==requested)
                break;
    }}

    *found=num_converged;

    for(int i=0;i<*n_new;i++)
        old_eig[i]=new_eig[i];

    *n_old=*n_new;

    cout<<"Number of converged eigenvalues found: "<<*found<<endl;

    for(int i=0;i<*found;i++)
	    cout<<converged_eig[i]<<endl;

    delete [] converge_temp;
    cout<<"Exiting block_converged "<<endl;
} //end of block_converged_eigval


void eig_sort(nml_dvector *eig, int dim){
    int cur=0;
    double temp=0.0;

    for(int i=0;i<dim-1;i++){
        cur=i+1;
        while(cur<dim){
	        if(eig[i]>eig[cur]){
	            temp=eig[cur];
	            eig[cur]=eig[i];
	            eig[i]=temp;		
	        }

	        cur++;
	}}
}


void find_block_converged_eigval_Method2(
    int         tridiag_size,   // in Lanczos iteration
    int p,                      // block size 
    nml_dvector *eigv_conv,     // inout actual converged eigenvalues
    int         *pEigenvalues,  // out number converged eigenvalues
    int         *n_nosuccess,   // out number no success
    int         *n_old,         // inout number old eigenvalues
    nml_dvector *eigv_old,      // inout actual old eigenvalues
    FILE        *fp_trace,      // in trace   log text file pointer
    int         verbose,        // in verbose diagnostic messages
    nml_dscalar emin,           // in eigenvalue search lower bound
    nml_dscalar emax,           // in eigenvalue search upper bound
    nml_dscalar tolerance,      // in eigenvalue convergence tolerance
    nml_dscalar resolution,     // in eigenvalue separation minimum
    FILE        *fp_log,        // in message log text file pointer
    nml_dcmatrix *Beta,
    nml_dcmatrix *xx 	
) 
{
    cout<<"Entered Block converged successfully"<<endl; 
    const nml_dscalar bias = 1e-80;   // PLEASE EXPLAIN ME
    const nml_dscalar gain = 1e-8;    // PLEASE EXPLAIN ME
    nml_dscalar deigv_min = DBL_MAX;    //double        deigv_min = 1.0e100;
    nml_dscalar eigv_min = DBL_MAX;    //double         eigv_min = 1.0e100;
    int  Neval_in_range = 0;    // number eigenvalues found
    nml_dvector *eval = nml_dv_new(tridiag_size);
                                        // actual eigenvalues found
    *n_nosuccess = 0;
    int j=0;
    Neval_in_range=*pEigenvalues;

    for (j=0;j<Neval_in_range;++j) 
        eval[j]=eigv_conv[j];
   	
    *pEigenvalues = 0;
    int iteration = 0;
    int flag;
    int num=tridiag_size-p;
    int k=0;
 
    for (j = 0;j < Neval_in_range; ++j) {
	    nml_dscalar product=0.0;
        flag=1;
	    //int row=j*p*(iterations+1);
	    nml_dscalar mod=0;

	    for(k=0;k<p;k++){
	   	    mod=xx[num+k][j].r*xx[num+k][j].r+xx[num+k][j].i*xx[num+k][j].i;	
		    mod=sqrt(mod);
		    product+=mod*Beta[k][k].r;
	    }//end of for k

        if(product>10*tolerance)
            flag=0;

	    cout<<"flag="<<flag<<endl;
    	if(flag==1){
        	int l = 0;
        	for (l = 0; l < *n_old; ++l) {
          		nml_dscalar deigv = fabs(eigv_old[l] - eval[j]);
          		if (deigv < deigv_min) {
            		deigv_min = deigv;
            		eigv_min = eval[j];
          		}//end of if
        	}//end of for j

            { // label 1
        	    eigv_conv[*pEigenvalues]=eval[j];
        	    if (verbose && (NULL != fp_log)) {
          		    fprintf(fp_log, "   eigenvalue candidate #%d:  %18.15f\n",
          		    *pEigenvalues, eigv_conv[*pEigenvalues]*MTRX_SCALE);
          		    fflush(fp_log);
        	    }

        	    ++(*pEigenvalues);
            } // end of label 1
        } //end of flag==1

		//#define Timing_detail
        #ifdef Timing_detail
        fprintf(fp_log,"%d = %1.15f\n",j,eval[j]);
        #endif
    } // end of for j (eigenvalues)

    if (verbose) {
        if (deigv_min < DBL_MAX)
            if (NULL != fp_log)
                fprintf(fp_log, "    E = %e  dE_min = %e\n", eigv_min*MTRX_SCALE, deigv_min*MTRX_SCALE);
    } // end of if

    *n_old = 0;

    {   //label 2
	    int iteration = 0;
    	for (iteration = 0; iteration < tridiag_size; ++iteration)
      		if (  (emin/MTRX_SCALE < eval[iteration]) && (eval[iteration] < emax/MTRX_SCALE)
       		   && ((0 == iteration) || (resolution < fabs((eval[iteration-1] - eval[iteration])/(bias + eval[iteration]))))) {
        	    eigv_old[*n_old] = eval[iteration];
        		++(*n_old);
        	}
    } // end of label 2

    ++(*n_nosuccess);

    nml_dv_delete(eval);
    cout<<"Leaving block converged successfully"<<endl;

} // find_block_converged_eigval


int gramSchmidt(
    nml_dcvector* xx_prev,
    nml_dcvector* xx,
    nml_extent length
) 
{
    int k;
    // Use Gram-schmidt algorithm to orthorgonalize two eigenvectors of a duplicate eigenvalue, 
	// and check if the resulted eigenvectors are orthogonal.
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
        invert += xx_prev[k].r*xx_prev[k].r+xx_prev[k].i*xx_prev[k].i;

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
	
    if(fabs(dotproduct.r)<1e-14&&fabs(dotproduct.i)<1e-14)
        return 1;

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

}// gramSchmidt 


//find eigenvalues of the block tridiagonal matrix using LAPACK
void solve_block_tridiagonal_eigenvalues(
    int tridiag_size,
    int block_size,	 	
    nml_dcmatrix *triad,
    nml_dvector *val,
    int *val_len,
    nml_dscalar  emin,
    nml_dscalar  emax,
    int requested1,
    double tolerance,
	int is_cb
) 
{

    // find solution using LAPACK routine zhbevx, all the arguments have to 
    // be pointers and you have to add an underscore to the routine name 

    char JOBZ='N';                  // (N,V) Compute eigenvalues only
    char RANGE='V';                 // (A,V) all eigenvalues in the half-open interval (VL,VU]
    char UPLO='L';                  // (U,L) Lower triangle of A is stored
    int N=tridiag_size;             // The order of the matrix A
    int KD=block_size;              // the number of subdiagonals
    int LDAB=KD+1;                  // The leading dimension of the array AB

    // struct complex AB (complex)  On entry, the upper or lower triangle of the Hermi-tian band matrix A, 
	// stored in the first KD+1 rows of the array.

    struct complex *AB=new struct complex [N*LDAB];
    for(int i=0;i<N;i++)
	    for(int j=0;j<LDAB;j++){
		    AB[i*LDAB+j].re=triad[i][j].r;
		    AB[i*LDAB+j].im=triad[i][j].i;
 	    }

    int LDQ=1;                      // The leading dimension of the array Q (not referenced)
    struct complex *Q=new struct complex[LDQ]; // (complex)   the N-by-N unitary matrix used in the reduction to tridiagonal form
    double VL=emin;                 // If RANGE='V', the lower bound of the interval to be searched for eigenvalues
    double VU=emax;                 // the upper bound of the interval to be searched for eigenvalues
    int IL=2;                       // the index (from smallest to largest) of the smallest eigenvalue to be returned
    int IU=3;                       // the index (from smallest to largest) of the largest eigenvalue to be returned.
    double ABSTOL=tolerance;        // The absolute error tolerance for the eigenvalues
    int M=0;                        // The total number of eigenvalues found
    double *W=new double [N];       // The first M elements contain the selected eigen-values in ascending order
    int LDZ=1;                      // The leading dimension of the array Z
    struct complex *Z=new struct complex [LDZ]; // (complex)  the first M columns of Z contain the orthonormal eigenvectors of the matrix 
						                        // corresponding to the selected eigenvalues (not referenced).
    struct complex *WORK=new struct complex[N]; // (complex)  COMPLEX*16 array, dimension (N)
    double *RWORK=new double [7*N];             // DOUBLE PRECISION array, dimension (7*N)
    int *IWORK=new int [5*N];                   // INTEGER array, dimension (5*N)
    int *IFAIL=new int [1];                     // If INFO > 0, then IFAIL contains the indices of the eigenvectors 
	                                            // that failed to converge.
    int INFO=0;                     //0:  successful exit, < 0:  if INFO = -i, the i-th argument had an illegal

    zhbevx( &JOBZ, &RANGE, &UPLO, &N, &KD, AB, &LDAB, Q, &LDQ, &VL, &VU, 
	        &IL, &IU, &ABSTOL, &M, W, Z, &LDZ, WORK, RWORK, IWORK, IFAIL, &INFO );

    cout<<"Number of eigenvalues in the range "<<emin<<" to "<<emax<<":  "<<M<<endl;

#if 0 // bug fix that sometimes we cannot get the eigenvalues in the valence band
    if(requested1!=-1)
        if(requested1<M)
	        M=requested1;

    for(int i=0;i<M;i++) {
        //printf("%i:    %e\n", i, W[i]);
		if(is_cb)
            *val=W[i];
		else
		    *val=W[M-i];
        val++; 
    }

    *val_len=M;
#else
    int req_eig;
	if(requested1!=-1)
	    if(requested1<M)
		    req_eig=requested1;
		else req_eig=M;
	for(int i=0; i<req_eig; i++)
	{
	    if(is_cb) *val=W[i];
		else *val=W[M-i];
		val++;
	}

	*val_len=req_eig;
#endif

    cout<<"INFO="<<INFO<<endl;

    delete []  W;
    delete [] WORK;
    delete [] RWORK;
    delete [] IWORK;
    delete [] IFAIL;
    delete [] Q;
    delete [] Z;
  
    //return M;
}

#ifdef ARPACK_EIGENVECTOR_SOLVER

void solve_block_tridiagonal_eigenvectors(int tridiag_size, int block_size, nml_dcmatrix *triad, nml_dvector *val,
                                          int *val_len, nml_dscalar  emin, nml_dscalar  emax, int requested1,
										  double tolerance, nml_dcmatrix *xx,int is_cb) 
{
	// input parameters for ARPACK eigensolvers
    int IDO = 0; // Initially set to 0
	char BMAT[] = "I";
	int N = tridiag_size;
	char WHICH[] = "LR";

	// just for sure make 1.5 times requested eigenvalues
	int NEV = 2*(*val_len);
	double TOL = tolerance;
	nml_dcvector *RESID = nml_dcv_new(N);
	int NCV = 4*NEV; 
	nml_dcvector *V = nml_dcv_new(NCV * N);
	int LDV = tridiag_size;
	int IPARAM[11];
	int IPNTR[14] = {0,0,0,0,0,0,0,0,0,0,0,0,0,0};
	nml_dcvector *WORKD = nml_dcv_new(3 * N);
	int LWORKL = (3 * NCV + 5) * NCV;
	nml_dcvector *WORKL = nml_dcv_new(LWORKL);
	nml_dvector *RWORK = nml_dv_new(NCV);
	int INFO = 0;

	// input parameters for ZNEUPD
	int RVEC = 1;
	char HOWMNY = 'A';
//	int SELECT[NCV];
        int *SELECT = (int*) nml_malloc(sizeof(int)*NCV);
	nml_dcvector *D = nml_dcv_new(NEV + 1);
	nml_dcvector *Z = nml_dcv_new(NEV *  N);
	int LDZ = N;
	complex SIGMA;
	//nml_dcvector *WORKEV = nml_dcv_new(2 * NCV);
	nml_dcvector *WORKEV = nml_dcv_new(3 * NCV);

    // Shift and invert operation variables
	nml_dcvector *B = nml_dcv_new(N);

	printf("[BLKLANCZ] Solve block tridiagonal eigenvectors\n");
	printf("[BLKLANCZ] tridiag_size = %d, block_size = %d, val_len = %d\n", tridiag_size, block_size, *val_len);
	printf("[BLKLANCZ] NCV = %d, tolerance = %1.3e, emin = %f, emax = %f\n", NCV, tolerance, emin, emax);

	if(!is_cb) WHICH[0]='S';WHICH[1]='R'; // For valence band eigenvalues

    // Shift set
	SIGMA.re = 0;
	SIGMA.im = 0;

	// REPACK in LINPACK format
	int KD = block_size;
	int LDAB = 3 * KD + 1;
    nml_dcmatrix *AB = nml_dcm_new(N,LDAB);
//	int IPIV[N];
        int *IPIV = (int*) nml_malloc(sizeof(int)*N);
	{
	    for(int i = 0; i < N; i++) IPIV[i] = 0;

	    for(int i = 0; i < N; i++)
		    for(int j = 0; j < KD+1; j++)
			{
			    AB[i][2*KD+j].r = triad[i][j].r;
			    AB[i][2*KD+j].i = triad[i][j].i;
			}

		for(int i = 0; i < N; i++) {
		    for(int j = 1; j <= KD; j++) {
			    if((i+j) < N) {
                    AB[i+j][2*KD-j].r = AB[i][2*KD+j].r;
                    AB[i+j][2*KD-j].i = -AB[i][2*KD+j].i;
		}}}

        // Make (A-sigma*I)
		if(is_cb)
		{
		    for(int i = 0; i < N; i++)
                AB[i][2*KD].r -= emin;
		}
		else
		{
		    for(int i = 0; i < N; i++)
                AB[i][2*KD].r -= emax;
		}
	}

	zgbtrf(&N, &N, &KD, &KD, (complex *)&AB[0][0], &LDAB, &IPIV[0], &INFO);
	printf("[BLKLANCZ] ZGBTRF INFO = %d\n", INFO);

    // IPARAM set
	IPARAM[0] = 1; // 0
	IPARAM[1] = 0;
	IPARAM[2] = 1500;
	IPARAM[3] = 1;
	IPARAM[4] = 0; // This is number of converged ritz values
	IPARAM[5] = 0;
	IPARAM[6] = 1; // Mode selection, Ax = mu*x
	IPARAM[7] = 0;
	IPARAM[8] = 1;
	IPARAM[9] = 1;
	IPARAM[10] = 1;

	do
	{
        znaupd(&IDO, BMAT, &N, WHICH, &NEV, &TOL, (complex *)&RESID[0], &NCV, (complex *)&V[0], &LDV, 
		       &IPARAM[0], &IPNTR[0], (complex *)&WORKD[0], (complex *)&WORKL[0], &LWORKL, (double *)&RWORK[0], &INFO);
	    
		#if 0
		    printf("IDO = %d, INFO = %d\n", IDO, INFO);
		    printf("IPARAM=%d %d %d %d %d %d %d %d %d %d %d\n", IPARAM[0], IPARAM[1], IPARAM[2], IPARAM[3], IPARAM[4],
		            IPARAM[5], IPARAM[6], IPARAM[7], IPARAM[8], IPARAM[9], IPARAM[10]);
		    printf("IPNTR=%d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", IPNTR[0], IPNTR[1], IPNTR[2], IPNTR[3], IPNTR[4],
		            IPNTR[5],IPNTR[6],IPNTR[7],IPNTR[8],IPNTR[9],IPNTR[10],IPNTR[11],IPNTR[12],IPNTR[13]);
		#endif

		nml_dcvector*  in = &WORKD[IPNTR[0]-1];
		nml_dcvector* out = &WORKD[IPNTR[1]-1];

	    // Compute Y = OP * X, where OP = inv[TRIAD-sigma*I]
		// This is a linear problem of (TRIAD-sigma*I)*Y = X
		// It can be quite easy to solve using Gauss elimination
		// printf("[BLKLANCZ] IDO = %d\n", IDO);
	    if(IDO == -1 || IDO == 1)
		{
			// N = tridiag size
			int NRHS = 1; // Number of columns of Y and X 
			// B = X
			int LDB = N;
			int INFO1 = 0;
			char TRANS = 'N';

            // 0. Copy Input to B 
			for (int i = 0; i < N; i++)
			{
                B[i].r = in[i].r;
				B[i].i = in[i].i;
			}

			// 1. USE ZHESV function: This is a linear solver for hermitian matrix
			zgbtrs(&TRANS, &N, &KD, &KD, &NRHS, (complex *)&AB[0][0], &LDAB, &IPIV[0], (complex *)&B[0], &LDB, &INFO1);

            if(INFO1 != 0)
			{
			    printf("[BLKLANCZ] (TRIAD-sigma*I)*Y = X cannot be computed!\n");
				break;
			}

			// 2. Copy B to out
			for (int i = 0; i < N; i++)
			{
                out[i].r = B[i].r;
				out[i].i = B[i].i;
			}
		}
		// Compute Y = M * X, since M = I, just copy X to Y
		else if ( IDO == 2 )
		{
            for(int i = 0; i < N; i++)
			{
                out[i].r = in[i].r;
				out[i].i = in[i].i;
			}
		}
		else if ( IDO == 3 )
		{
		    if(IPARAM[0] == 0)
			{
                int NP = IPARAM[7];

                //printf("NP = %d NCV-NEV = %d\n", NP, NCV-NEV);
				// compute and return the shifts in the first NP locations of WORKL
				// What does this mean?? What to compute
				//for(int i = 0; i < NCV; i++)
			    //   printf("WORKL[IPNTR[5]+%d] = %1.10e + i %1.10e\n", i, WORKL[IPNTR[5]-1+i].re, WORKL[IPNTR[5]-1+i].im);

				// return the shifts in the first NP locations of WORKL
				for(int i = 0; i < NP; i++)
				{
                    WORKL[IPNTR[13]-1+i].r = SIGMA.re;
					WORKL[IPNTR[13]-1+i].i = SIGMA.im;
				}
			}
		}
		else 
		    break;
    } while ( IDO != 99);

	if (INFO < 0) printf("[BLKLANCZ] Error solving znaupd error = %d\n", INFO);
	else
	{
	    if(mpi_n3d_id == 0)
		{
		    cout << "===== Results of ZNAUPD =====" << endl;
			cout << "The Number of Ritz values requested is " << NEV << endl;
			cout << "The Number of Arnoldi vectors generated is " << NCV << endl;
			cout << "What Portion of the spectrum: " << WHICH << endl;
			cout << "The Number of Converged Ritz values is " << IPARAM[4] << endl;
			cout << "The Number of Implicit Arnoldi update iterations taken "<< IPARAM[2] << endl;
			cout << "The Number of OP*x is " << IPARAM[8] << endl;
			cout << "The Number of reorthogonalization steps is " << IPARAM[10] << endl;
			cout << "The Convergence criterion is " << tolerance << endl;

		    zneupd(&RVEC, &HOWMNY, &SELECT[0], (complex *)&D[0], (complex *)&Z[0], &LDZ, 
			       &SIGMA, (complex *)&WORKEV[0], BMAT, &N, WHICH, &NEV, &TOL, 
			       (complex *)&RESID[0], &NCV, (complex *)&V[0], &LDV, IPARAM, IPNTR, 
				   (complex *)&WORKD[0], (complex *)&WORKL[0], &LWORKL, (double *)&RWORK[0], &INFO);
		}
    }

	printf("[BLKLANCZ] ZNEUPD INFO = %d\n", INFO);

    if (INFO != 0) printf("[BLKLANCZ] Error solving zneupd\n");
	else if (IPARAM[4] > 0 && !mpi_n3d_id)
	{
	    printf("[BLKLANCZ] Found Eigenvalues!\n");
		for(int i = 0; i < IPARAM[4]; i++)
		    printf("           D[%d] = %f + i %f\n", i, D[i].r, D[i].i);

		printf("[BLKLANCZ] Reverting to original values\n");
		for(int i = 0; i < IPARAM[4]; i++)
		{
		    if(is_cb)
		        D[i].r = (1.0)/(D[i].r) + emin;
			else
		        D[i].r = (1.0)/(D[i].r) + emax;
		    printf("           D[%d] = %f\n", i, D[i].r);
	        //val[i] = D[i].r;		
		}

        // Filter out values and return corresponding vectors
		int indx_temp = 0;
		int min_indx;
//		double val_temp[*val_len];
                double *val_temp = (double*) nml_malloc(sizeof(double)*(*val_len));
		double min_temp;
//		char mark[IPARAM[4]];
                char *mark = (char*) nml_malloc(sizeof(char)*IPARAM[4]);
		memset(&mark[0], 0, sizeof(char)*IPARAM[4]);

		for(int i = 0; i < *val_len; i++)
		{
			min_indx = -1;
			min_temp = 500;
			int j;

		    for(j = 0; j < IPARAM[4]; j++)
			{
			    if((fabs(val[i]-D[j].r) < min_temp) && (mark[j] == 0))
				{
				    min_temp = fabs(D[j].r-val[i]);
					min_indx = j;
				}
			}

            // Find the nearest eigenvalues that matches tolerance criteria
			if(min_indx != -1)
			{
			    //if(min_temp <= TOL)
				{
			        val_temp[indx_temp] = D[min_indx].r;//val[i];
				    mark[min_indx] = 1;

				    for(int k = 0; k < N; k++)
				    {
                        xx[indx_temp][k].r = Z[min_indx * N + k].r;
                        xx[indx_temp][k].i = Z[min_indx * N + k].i;
				    }

				    indx_temp++;

				    printf("[BLKLANCZ] val[%d] = %1.10e, D[%d] = %1.10e\n", i, val[i], min_indx, D[min_indx].r);
				}
			}
			else
			    printf("[BLKLANCZ] Index matching error tolerance not found! val[%d] = %1.10e min_error = %1.10e\n", i, val[i], min_temp);
			   
		}

		*val_len = indx_temp;
		for(int i = 0; i < indx_temp; i++)
		    val[i] = val_temp[i];

  	        nml_free(mark);
  	        nml_free(val_temp);
	}

	#if 0
	if(INFO == 0)
    {
	    for(int j = 0; j < *val_len; j++) {
            for(int k = 0; k < N; k++) {
                xx[j][k].r=Z[j * N + k].r;
		        xx[j][k].i=Z[j * N + k].i;
    }}}
	#endif

	// Deleting
  	nml_free(SELECT);
  	nml_free(IPIV);
    nml_dv_delete(RWORK);
    nml_dcv_delete(WORKD);
    nml_dcv_delete(D);
    nml_dcv_delete(WORKL);
	nml_dcv_delete(B);
	nml_dcv_delete(RESID);
	nml_dcv_delete(WORKEV);
	nml_dcv_delete(Z);
	nml_dcm_delete(AB);
	nml_dcv_delete(V);

    cout<<"[BLKLANCZ] Exit from ARPACK eigensolver"<<endl;
}

#else // END OF ARPACK_EIGENVECTOR_SOLVER

void solve_block_tridiagonal_eigenvectors(
    int tridiag_size,
    int block_size,
    nml_dcmatrix *triad,
    nml_dvector *val,
    int *val_len,
    nml_dscalar  emin,
    nml_dscalar  emax,
    int requested1,
    double tolerance,
    nml_dcmatrix *xx	
) 
{
    // find solution using LAPACK routine zhbevx, all the arguments have to
    // be pointers and you have to add an underscore to the routine name
    char JOBZ='V';       // (N,V) Compute eigenvectors also
    char RANGE='V';      // (A,V) all eigenvalues in the half-open interval (VL,VU]
    char UPLO='L';       // (U,L) Lower triangle of A is stored
    int N=tridiag_size;  // The order of the matrix A
    int KD=block_size;   // the number of subdiagonals
    int LDAB=KD+1;       // The leading dimension of the array AB

    // struct complex AB (complex)  On entry, the upper or lower triangle of the Hermi-tian band matrix A, 
	// stored in the first KD+1 rows of the array.
    struct complex *AB= new struct complex [N*LDAB];
    for(int i=0;i<N;i++)
        for(int j=0;j<LDAB;j++){
            AB[i*LDAB+j].re=triad[i][j].r;
            AB[i*LDAB+j].im=triad[i][j].i;
        }

    int LDQ=N;                  // The leading dimension of the array Q
    struct complex *Q=new struct complex [LDQ*N];     
	                            // (complex)   the N-by-N unitary matrix used in the reduction to tridiagonal form
    double VL=emin;             // If RANGE='V', the lower bound of the interval to be searched for eigenvalues
    double VU=emax;             // the upper bound of the interval to be searched for eigenvalues
    int IL=2;                   // the index (from smallest to largest) of the smallest eigenvalue to be returned
    int IU=3;                   // the index (from smallest to largest) of the largest eigenvalue to be returned.
    double ABSTOL=tolerance;    // The absolute error tolerance for the eigenvalues
    int M=0;                    // The total number of eigenvalues found
    double *W=new double [N];   // The first M elements contain the selected eigen-values in ascending order
    int LDZ=N;                  // The leading dimension of the array Z
    struct complex *Z=new struct complex [LDZ*N]; 
	                            // (complex)  the first M columns of Z contain the orthonormal eigenvectors 
							    // matrix corresponding to the selected eigenvalues.
    struct complex *WORK=new struct complex [N]; // (complex)  COMPLEX*16 array, dimension (N)
    double *RWORK=new double [7*N];              // DOUBLE PRECISION array, dimension (7*N)
    int *IWORK=new int [5*N];                    // INTEGER array, dimension (5*N)
    int *IFAIL=new int [N];     //If INFO > 0, then IFAIL contains the indices of the eigenvectors that failed to converge.
    int INFO=0;                 //0:  successful exit
                                //< 0:  if INFO = -i, the i-th argument had an illegal

    zhbevx( &JOBZ, &RANGE, &UPLO, &N, &KD, AB, &LDAB, Q, &LDQ, &VL, &VU, &IL, 
	        &IU, &ABSTOL, &M, W, Z, &LDZ, WORK, RWORK, IWORK, IFAIL, &INFO );

    cout<<"Number of eigenvalues in the range "<<emin<<" to "<<emax<<":  "<<M<<endl;

    if(requested1!=-1){
        if(requested1<M)
            M=requested1;
    }  

    cout<<"Info="<<INFO<<endl;

    if(INFO>0){
	    cout<<"IFAIL"<<endl;
	    for(int i=0;i<M;i++)
		    cout<<IFAIL[i]<<"  "<<endl;
    }

    for(int i=0;i<M;i++){
        printf("%i:    %f\n", i, W[i]);
        *val=W[i];
        val++;
    }

    *val_len=M;

    cout<<"val len:"<<*val_len<<endl;

    for(int j=0;j<M;j++){
        //cout<<"Num="<<j<<endl;
        for(int k=0;k<N;k++){
            xx[j][k].r=Z[j*N+k].re;
		    xx[j][k].i=Z[j*N+k].im;
        }
    }

    delete [] AB;
    delete [] W;
    delete [] WORK;
    delete [] RWORK;
    delete [] IWORK;
    delete [] IFAIL;
    delete [] Q;
    delete [] Z;

    cout<<"Exit from LAPACK eigensolver"<<endl;
}
#endif // ARPACK_EIGENVECTOR_SOLVER


int block_lanczos_blocktridiagonal(     //   out number of converged eigenvalues
    nml_dvector  *eigv_conv,    //   out eigenvalues                    
    nml_dcmatrix *triad,
    nml_extent   *pIterations,  //   out actual number of iterations    
    int p,                      // block size                           
    vector <nml_dcvector*> *R_n,         // inout r_{n} = q_{n+1}*beta_{n}       
    vector <nml_dcvector*> *Q_n,         // inout current  complex Lanczos vector
    vector <nml_dcvector*> *Q_n1,        // inout previous complex Lanczos vector
    nml_extent   length,        // in extent of vectors r_n, q_n & q_n1 
    int          requested1,    // in requested number of eigenvalues   
    int          requested2,    // in requested number of eigenvalues   
    int ConvCheckStartIter,     // in convergence check start iteration 
    int ConvCheckSkipRate,      // in convergence check skip rate       
    int          imax,          // in maximum number of iterations + 1  
    nml_dscalar  emin1,         // in eigenvalue search lower bound     
    nml_dscalar  emax1,         // in eigenvalue search upper bound     
    nml_dscalar  emin2,         // in eigenvalue search lower bound     
    nml_dscalar  emax2,         // in eigenvalue search upper bound     
    nml_dscalar  tolerance,     // in eigenvalue convergence tolerance  
    nml_dscalar  resolution,    // in eigenvalue separation minimum     

	#ifndef BLANCZ_MATMULT
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
	#else
    void (*matmul)(const int**, nml_dcmatrix*, const nml_dcmatrix*),
	#endif

    const int*   argument[],    // matrix-vector multiply argument list 
    FILE         *fp_trace,     // in   trace log text file pointer     
    FILE         *fp_log,       // in message log text file pointer     
    int          verbose,        // in verbose diagnostic messages 
    int *tridiag_size,
    double deflation_tolerance,
    int max_no_success,
    int *found_vb,
    int *found_cb,
    int deflation_on,
    int conv_method			       
) 
{
    // Use Lanczos' algorithm to tridiagonalize complex hermitian matrix. 

    cout<<"Entering block lanczos blocktridiagonal"<<endl;	
    int converged = 1;     // no convergence yet                   
    int n_old1 = 0;             // number of old eigenvalues            
    int n_old2 = 0;             // number of old eigenvalues            
    int iterations = 0;         // Lanczos iteration number             
    int skip = 0;               // number of convergence checks skipped 
    int skip_write = 0;         // not actually used in this version    
    int n_nosuccess1 = 0;       // number of no success                 
    int n_nosuccess2 = 0;       // number of no success                 
    int eigenvalues = 0;        // number of converged eigenvalues      
    int eigenvalues1 = 0;       // number of converged eigenvalues in range1    
    int eigenvalues2 = 0;       // number of converged eigenvalues in range2    
    int added=1;			//solve for some extra eigenvalues as well
    int res2=0;
    int deflate=-1;
    int block_size=p;
    int found1=0;
    int found2=0;
    nml_dscalar* eigv_old1 = (nml_dscalar*)nml_malloc(sizeof(nml_dscalar)*(requested1+added));
    nml_dscalar* eigv_old2 = (nml_dscalar*)nml_malloc(sizeof(nml_dscalar)*(requested2+added));
    nml_dscalar* eigv_new1 = (nml_dscalar*)nml_malloc(sizeof(nml_dscalar)*(requested1+added));
    nml_dscalar* eigv_new2 = (nml_dscalar*)nml_malloc(sizeof(nml_dscalar)*(requested2+added));
    nml_dscalar* eigv_conv1=(nml_dscalar*)nml_malloc(sizeof(nml_dscalar)*(requested1));
    nml_dscalar* eigv_conv2=(nml_dscalar*)nml_malloc(sizeof(nml_dscalar)*(requested2));

    int indx=0; 
    nml_dcmatrix *Alpha = nml_dcm_new(block_size, block_size); 
    nml_dcmatrix *Beta =  nml_dcm_new(block_size, block_size);

    if (mpi_n3d_id != mpi_n3d_masterid) {
        #ifdef MPI3d
        do {
            int sndrcv[2];   //what's this?

            MPI_TIC(start_time2);
            single_iteration(matmul, argument, R_n, Q_n, Q_n1, length, iterations, p, 
			                 block_size, Alpha, Beta, &deflate, deflation_tolerance, deflation_on);

  	        if(deflate!=-1)
		        p--;

	        MPI_TOC(time_lanczos_iter,start_time2);

            ++iterations;
            ++skip;
            ++skip_write;

            MPI_TIC(start_time);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_TOC(time_barrier, start_time);
            MPI_Bcast(&sndrcv[0],  2, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

            converged = sndrcv[0];
            eigenvalues = sndrcv[1];

            MPI_Bcast(&eigv_conv[0], eigenvalues, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);
        } while ((converged==0)&&(iterations+1<imax));
        #else
        error("Non-MPI code should not reach this statements.\n");
        #endif
    }
    else {                // (mpi_n3d_id == mpi_n3d_masterid)             
        #ifdef WriteOutPentaMatrix
        FILE *fp=fopen("penmatrix.out","w");
        #endif

        do {
            MPI_TIC(start_time2);
            single_iteration(matmul, argument, R_n, Q_n, Q_n1,length,iterations, p, 
			                 block_size, Alpha, Beta, &deflate, deflation_tolerance, deflation_on);
            MPI_TOC(time_lanczos_iter, start_time2);

            #ifdef WriteOutPentaMatrix
            //fprintf(fp,"%1.15e %1.15e %1.15e %1.15e\n",alpha,beta1.r, beta1.i,gamma1);
            //fprintf(fp,"%1.15e %1.15e %1.15e %1.15e\n",alpha2,beta2.r,beta2.i,gamma2);
            #endif

	        //writing data to triadiagonal matrixi
	        if(iterations!=0){
 		        indx+=p;
		        for(int j=0;j<p;j++)
                	for(int i=0;i<=j;i++){
				        triad[indx-p+j][p-j+i]=Beta[i][j];
			}}
	
		    if(deflate!=-1){
			    cout<<"Value of deflate in block_tridiagonal:"<<deflate<<endl;
                p--;
			    cout<<"Value of p:"<<p<<endl;
		    }	

	        for(int j=0; j<p;j++){
		        for(int i=j;i<p;i++){
			        triad[indx+j][i-j]=Alpha[i][j];
		    }}

	        *tridiag_size+=p;
            //cout<<"Iteration="<<iterations<<"    tri_size:  "<<*tridiag_size<<endl;

            /*	cout<<"After Iteration="<<iterations<<endl;
	            cout<<"Alpha"<<endl;
	            print_Matrix(Alpha,p,p);
	            cout<<"Beta"<<endl;
	            print_Matrix(Beta,p,p);
 	            cout<<"Triad"<<endl;
	            print_Matrix(triad, (iterations+1)*p, p+1);
            */	

            ++skip;
            ++skip_write;
        
            if (verbose > 1 && (NULL != fp_log)) {
        	    //fprintf(fp_log, "it=%d al=%9.5f be=%9.5f\n", iterations, alpha, beta1.r);
        	    //fflush(fp_log);
            }
            else
               if (0 == iterations%100 && (NULL != fp_log)) {
                  fprintf(fp_log, "Iteration count = %d\n", iterations);
                  fflush(fp_log);
               }

            if (ConvCheckStartIter < iterations){	
      	        if(found1<requested1)	
		            n_nosuccess1++;

      	        if(found2<requested2)	
        	        n_nosuccess2++;
            }	

            if ((ConvCheckStartIter < iterations) && (ConvCheckSkipRate <= skip)) {
                skip = 0;
 
	            if((requested1>0)&&(found1<requested1)){
    		        solve_block_tridiagonal_eigenvalues(*tridiag_size, block_size, triad,
                        eigv_new1, &eigenvalues1, emin1, emax1, requested1+added, tolerance,0);

                    find_block_converged_eigvalues_Method4(eigv_new1, &eigenvalues1, eigv_old1, &n_old1, eigv_conv1, &found1,
                        requested1, &n_nosuccess1,  fp_trace, verbose, tolerance, fp_log);

		            *found_vb=found1;			
	            }//end of requested1

                if ((requested2 > 0)&&(found2<requested2)) {
		            solve_block_tridiagonal_eigenvalues(*tridiag_size, block_size, triad, 
			          eigv_new2, &eigenvalues2, emin2, emax2, requested2+added, tolerance,1);

                    find_block_converged_eigvalues_Method4(eigv_new2, &eigenvalues2, eigv_old2, &n_old2, eigv_conv2, &found2, 
       			       requested2, &n_nosuccess2,  fp_trace, verbose, tolerance, fp_log);
		
		            *found_cb=found2;
                }//end of requested2
            }//end of convergence check

            eigenvalues = found1 + found2;

	        if (((found1 >= requested1) && (found2 >=requested2))||(iterations+2==imax)){
		        for(int i=0;i<found1;i++)
			        eigv_conv[i]=eigv_conv1[i];

		        for(int i=0;i<found2;i++)
                    eigv_conv[i+found1]=eigv_conv2[i];

		        eig_sort(eigv_conv,found1+found2);

		        for(int i=0;i<eigenvalues;i++)
			        cout<<"Eig "<<i<<" :  "<<eigv_conv[i]<<endl;	    

                converged=1;
	        }//end of if
            else
                converged=0;

	        if((n_nosuccess1>=max_no_success)||(n_nosuccess2>=max_no_success)){
		        converged=1;	
		        cout<<"Maximum number of no_successes reached. Termination in progress"<<endl;
		        cout<<"no_success_vb: "<<n_nosuccess1<<"     no_success_cb:   "<<n_nosuccess2<<endl; 
	        }

            //cout<<"Convergence Result: "<<converged<<endl;

            #ifdef  MPI3d
            if (1 < mpi_n3d_numprocs) {
                int sndrcv[2];
                sndrcv[0] = converged;
                sndrcv[1] = eigenvalues;

                MPI_TIC(start_time);
                MPI_Barrier(MPI_COMM_WORLD);
                MPI_TOC(time_barrier, start_time);
                MPI_TIC(start_time);
                MPI_Bcast(&sndrcv[0], 2, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);
                MPI_Bcast(&eigv_conv[0], eigenvalues, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);
                MPI_TOC(time_bcast, start_time);

                if (NULL != fp_log)
                    fflush(fp_log);
            }
            #endif

            iterations++;

        //} while(iterations<imax);
        } while ((converged==0)&&(iterations+1<imax));

        #ifdef WriteOutPentaMatrix
        fclose(fp);
        #endif
    }

    *pIterations = iterations;

    nml_free(eigv_old1);
    nml_free(eigv_old2);
    nml_free(eigv_new1);
    nml_free(eigv_new2);
    nml_free(eigv_conv1);
    nml_free(eigv_conv2);
    nml_dcm_delete(Alpha);
    nml_dcm_delete(Beta);

    cout<<"end of block tridiagonal"<<endl;

    return eigenvalues;
} // block_lanczos_pentadiagonal 

                              
nml_extent eigenvaluesBlockLanczos(     // number of eigenvalues actually found 
    nml_dvector  *value,        //   out real eigenvalue vector        
    nml_dcmatrix *triad,
    nml_extent   *pIterations,  //   out actual number of iterations    
    int p,                      // block size                           
    vector <nml_dcvector*> *R_n,         // inout r_{n} = q_{n+1}*beta_{n}       
    vector <nml_dcvector*> *Q_n,         // inout current  complex Lanczos vector
    vector <nml_dcvector*> *Q_n1,        // inout previous complex Lanczos vector
    nml_extent   length,        // in extent of vectors r_n, q_n & q_n1 
    nml_extent   requested1,    // in requested number of eigenvalues   
    nml_extent   requested2,    // in requested number of eigenvalues   
    int ConvCheckStartIter,     // in convergence check start iteration 
    int ConvCheckSkipRate,      // in convergence check skip rate       
    nml_extent   imax,          // in maximum number of iterations + 1  
    nml_dscalar  emin1,         // in eigenvalue search lower bound     
    nml_dscalar  emax1,         // in eigenvalue search upper bound     
    nml_dscalar  emin2,         // in eigenvalue search lower bound     
    nml_dscalar  emax2,         // in eigenvalue search upper bound     
    nml_dscalar  tolerance,     // in eigenvalue convergence tolerance  
    nml_dscalar  resolution,    // in eigenvalue separation minimum     

	#ifndef BLANCZ_MATMULT
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
	#else
    void (*matmul)(const int**, nml_dcmatrix*, const nml_dcmatrix*),
	#endif 

    const int*   argument[],    // matrix-vector multiply argument list 
    FILE         *fp_trace,     // in   trace log text file pointer     
    FILE         *fp_log,       // in message log text file pointer     
    int          verbose,       // in verbose diagnostic messages
    int *tridiag_size,
    double deflation_tolerance,
    int max_no_success,
    int *found_vb,
    int *found_cb,
    int deflation_on,
    int conv_method				  
) 
{
    cout<<"Reached eigenvalues_block_lanczos"<<endl; 

    MPI_TIME_TIC(start_time_lanc_all);

    int eigenvalues = block_lanczos_blocktridiagonal(
      value, triad, pIterations, p, R_n, Q_n, Q_n1, length,
      requested1, requested2, ConvCheckStartIter, ConvCheckSkipRate, imax,
      emin1, emax1, emin2, emax2, tolerance, resolution, matmul, argument,
      fp_trace, fp_log, verbose, tridiag_size, deflation_tolerance, max_no_success,
      found_vb, found_cb, deflation_on, conv_method);
    
    MPI_TOC(time_lanczos_iter_all, (start_time_lanc_all + time_lanczos_setup_once));

    #ifdef MPI_TIMING
    if (NULL != fp_log) {
        fprintf(fp_log, "\ntime_lanczos_iter_all = %g\n", time_lanczos_iter_all);
        fflush(fp_log);
    }
    #endif // MPI_TIMING 

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

    cout<<"end of eigenvalues_block_lanczos"<<endl;
    // Sort the eigenvalues in value in ascending order.                  
    //rsort(value, eigenvalues);

    return eigenvalues;
}//end of method


//static int convMethodVector=2;

void saveBlockTridiagonal(
    nml_dcmatrix *triad ,
    int triad_size,   // in    Lanczos iterations 
    int p 
) 
{   // save real symmetric tridiagonal matrix T  
	FILE *fp_block=fopen("Block_tridiagonal.dat","w");

	if(NULL!=fp_block){
	    fprintf(fp_block,"%i %i\n",triad_size,4);
      	for(int i=0;i<triad_size;i++){
		    for(int j=0;j<p+1;j++){
		        fprintf(fp_block,"%i %i %e %e\n",i+1,j+1,triad[i][j].r,triad[i][j].i);
		}} 

   	    fclose(fp_block);
    }else
        cout<<"Error in opening file to dump block tridiagonal"<<endl;
}     // saveBlockTridiagonalMatrix  



void saveBlockTridiagonal_MatrixFormat(
    nml_dcmatrix *triad ,
    int iterations,   // in Lanczos iterations
    int p
) 
{  // save real symmetric tridiagonal matrix T

	int elements=0;
    FILE *fp_block=fopen("Block_tridiagonal_Matrix_Format.dat","w");
    if(NULL!=fp_block){
        fprintf(fp_block,"%i %i\n",(2*(iterations*p*(p+1))+(p*(p+1)))-(iterations+1)*p,4);
        for(int i=0;i<iterations*p;i++){
            for(int j=0;j<p+1;j++){
			elements++;
            fprintf(fp_block,"%i %i %e %e\n",i+j+1,i+1,triad[i][j].r,triad[i][j].i);

			if(i!=(i+j)){
				elements++;
				fprintf(fp_block,"%i %i %e %e\n",i+1,i+j+1,triad[i][j].r,(-triad[i][j].i));
		}}}

	    int count=p;
	    for(int i=iterations*p;i<(iterations+1)*p;i++){
            for(int j=0;j<count+1;j++){
			    count--;
			    elements++;
                fprintf(fp_block,"%i %i %e %e\n",i+j+1,i+1,triad[i][j].r,triad[i][j].i);

			    if(i!=(i+j)){
				    elements++;
                    fprintf(fp_block,"%i %i %e %e\n",i+1,i+j+1,triad[i][j].r,(-triad[i][j].i));
         }}}

	    fclose(fp_block);
    }
	else 
        cout<<"Error in opening file to dump block tridiagonal"<<endl;
	
    cout<<"Number of elements="<<elements<<endl;
    cout<<"Calculated Number of elements="<<(2*(iterations*p*(p+1))+(p*(p+1)))-(iterations+1)*p<<endl;
}     // saveBlockTridiagonalMatrix



void readBlockTridiagonal(nml_dcmatrix *triad){
    int row=0;
    int col=0; 
    int row_index=0;
    int col_index=0;
    double real_part=0.0;
    double image_part=0.0;

    FILE *fp_block=fopen("Block_tridiagonal.dat","r");
    if(NULL!=fp_block){
        fscanf(fp_block,"%i %i\n",&row,&col);
        for(int i=0;i<row;i++){
            fscanf(fp_block,"%i %i %e %e\n",&row_index,&col_index,&real_part, &image_part);
	        row_index--;
	        col_index--;
	        cout<<"row="<<row_index<<"   col="<<col_index<<"  real_p="<<real_part<<"  im_p="<<image_part<<endl;	   
            triad[row_index][col_index].r=real_part;
	        triad[row_index][col_index].i=image_part;
        }
        
	    fclose(fp_block);
    }
	else
        cout<<"Error in opening file to read block tridiagonal"<<endl;
}


nml_extent eigenvectorsBlockLanczos(	// number of eigenvalues actually found	
    nml_dcvector *cvalue,	//   out complex eigenvalues		
    nml_dcmatrix *vector,	//   out complex eigenvectors
    nml_dcmatrix *xx,		// In Ritz vectors
    nml_dcmatrix *xx2,	
    const
    nml_dvector  *value,	// in    real    eigenvalues		
    nml_extent	 eigenvalues,	// in	 number of eigenvalues
    nml_extent   eigvb,
    nml_extent   eigcb,			
    nml_extent	 iterations,	// in    number of iterations
    int 	triad_size,		
    nml_dcmatrix   *triad,
    std::vector <nml_dcvector*> *R_n, // inout r_{n} = q_{n+1}*beta_{n}	
    std::vector <nml_dcvector*> *Q_n, // Q_n
    std::vector <nml_dcvector*> *Q_n1, //Q_n_1	
    nml_extent	 length,	// in extent of vectors r_n, q_n & q_n1	
    int p,			// block size

	#ifndef BLANCZ_MATMULT
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
	#else
    void (*matmul)(const int**, nml_dcmatrix*, const nml_dcmatrix*),
	#endif

    const int*	 argument[],	// matrix-vector multiply argument list	
    nml_dscalar  tolerance,     // in eigenvalue convergence tolerance  
    FILE	 *fp_trace,	// in   trace log text file pointer	
    FILE	 *fp_log,	// in message log text file pointer	
    int		 verbose,	// in verbose diagnostic messages	
    double def_tol,		// deflation tolerance	
    int deflation_on,
    int conv_vec_method
) 
{
    //MPI_TIME_TIC(start_time_lanc_all);

	//nml_extent Eigenvalues = eigenvalues;
    //MPI_Bcast(&Eigenvalues, 1, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);
    
    //nml_dmatrix	*xm = nml_dm_new(eigenvalues, triad_size);
    //nml_dmatrix	*xm_imag = nml_dm_new(eigenvalues, triad_size);
	//if((xm == NULL) || (xm_imag == NULL))
	//    printf("[BLKLANCZ] Cannot allocate memory!!\n");

    nml_dcmatrix *Alpha=nml_dcm_new(p, p);	
    nml_dcmatrix *Beta=nml_dcm_new(p, p);   

    int* num_iteration=NULL;
    int flag;
  
    if(conv_vec_method==2)
        num_iteration = (int*)nml_malloc(sizeof(int)*(eigenvalues));

    /*if (mpi_n3d_masterid == mpi_n3d_id)*/ {
	    nml_offset iteration = 0;
	    int i;
    
   	    //label1
		// Just copy xx/xx2 -> xm/xm_img 
    	{ 
		    nml_offset eigenvalue = 0;

      		for (eigenvalue = 0; eigenvalue < eigenvalues; ++eigenvalue) {
		        iteration = 0;
		        cvalue[eigenvalue].r = value[eigenvalue];

		        // Method 1 
		        if(conv_vec_method==1){
	  		        for (iteration = 0; iteration < triad_size; iteration++) {
				        if(eigenvalue<eigvb){
	    				    //xm[eigenvalue][iteration] = xx[eigenvalue][iteration].r;
	    				    //xm_imag[eigenvalue][iteration] =xx[eigenvalue][iteration].i;
        			    }
						else{
					        //xm[eigenvalue][iteration] = xx2[eigenvalue-eigvb][iteration].r;
                            //xm_imag[eigenvalue][iteration] =xx2[eigenvalue-eigvb][iteration].i;
				        }  
			        }
		        }//end of method1

		        // Method 2 
		        else if(conv_vec_method==2){
	  		        num_iteration[eigenvalue]=iterations;
	  		        flag=0;
	  		        int num_con_checking=5;

	  		        for (iteration =0; iteration<triad_size;iteration++){ //for 2
			            if(eigenvalue<eigvb){	
	    		            if(iteration>0 && fabs(xx[eigenvalue][iteration].r)<10*tolerance){//if 3
	      		                if(flag<num_con_checking){//if 2
				                    if(iteration==(unsigned int)(num_iteration[eigenvalue]+1)&&flag<num_con_checking)
		  				                flag++;
				                    else 
	          				            flag=0;

			       	                num_iteration[eigenvalue]=iteration;
							}}

	    		    	    //xm[eigenvalue][iteration] = xx[eigenvalue][iteration].r;
	    		    	    //xm_imag[eigenvalue][iteration] = xx[eigenvalue][iteration].i;
			            }
						else{
				            if(iteration>0 && fabs(xx2[eigenvalue-eigvb][iteration].r)<10*tolerance){//if 3
                               	if(flag<num_con_checking){//if 2
                                    if(iteration==(unsigned int)(num_iteration[eigenvalue]+1)&&flag<num_con_checking)
                                        flag++;
                                   	else
                                        flag=0;
										
                                    num_iteration[eigenvalue]=iteration;
                            }}

                            //xm[eigenvalue][iteration] = xx2[eigenvalue-eigvb][iteration].r;
                            //xm_imag[eigenvalue][iteration] = xx2[eigenvalue-eigvb][iteration].i;
			   	        } //end of eigcb 		
	  		        }// end of for 2 (iteration)
		        }// end of method 2
       	    }//end of eigenvalue loop
        }// end of label1
#if 0
        printf("[BLKLANCZ] CPU %d starting Broadcast of eigenvectors\n", mpi_n3d_id);

        //cout<<"triad size:"<<triad_size<<endl;
        //cout<<"eigenvalues:"<<eigenvalues<<endl;

        MPI_TIC(start_time);

        MPI_Bcast(&xm[0][0], eigenvalues*triad_size, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);
        MPI_Bcast(&xm_imag[0][0], eigenvalues*triad_size, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);

        if(conv_vec_method==2)
            MPI_Bcast(&num_iteration[0], eigenvalues, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

        MPI_TOC(time_bcast, start_time);

        printf("[BLKLANCZ] CPU %d finished Broadcast of eigenvectors\n", mpi_n3d_id);
#endif

        //fclose(fin);
    }
#if 0
    else {			// (mpi_n3d_masterid != mpi_n3d_id)		
        //cout<<"triad size:"<<triad_size<<endl;
        //cout<<"eigenvalues:"<<eigenvalues<<endl;

        //#ifdef MPI3d
        // Receive the subspace eigenvectors from the master CPU.		
        printf("[BLKLANCZ] CPU %d expected Broadcast of eigenvectors\n", mpi_n3d_id);


        MPI_TIC(start_time);

        MPI_Bcast(&xm[0][0], eigenvalues*triad_size, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);
        MPI_Bcast(&xm_imag[0][0], eigenvalues*triad_size, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);

        if(conv_vec_method==2)
            MPI_Bcast(&num_iteration[0], eigenvalues, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

        MPI_TOC(time_bcast, start_time);


        printf("[BLKLANCZ] CPU %d received Broadcast of eigenvectors\n", mpi_n3d_id);
    }

    if (mpi_n3d_id == mpi_n3d_masterid) {
        if (NULL != fp_log) {		// PLEASE EXPLAIN ME		
            fprintf(fp_log, "Repeat Lanczos Iterations and construct eigenvectors\n" "iteration:\n");
            fflush(fp_log);
    }}
#endif

    //start of label2
    {
        nml_offset iteration = 0;
        nml_dscalar dummy=0.0;
        int pc=p;		
        int def=-1;
        int tri_size=0;
     
        //start of for iteration 
        for (iteration = 0; iteration < iterations; iteration++) {
	        if ((mpi_n3d_id == mpi_n3d_masterid) && (0 == iteration%100)) {
		        if (NULL != fp_log) {
	  		        fprintf(fp_log, "%d..\n", iteration);
	  		        fflush(fp_log);
		    }}

	        MPI_TIC(start_time2);
  	        single_iteration(matmul, argument, R_n, Q_n, Q_n1, length, iteration, pc, p, Alpha, Beta, &def, def_tol, deflation_on);	
	
	        if(def!=-1){
                cout<<"Value of deflate in block_tridiagonal:"<<def<<endl;
                pc--;
                cout<<"Value of pc:"<<pc<<endl;
            }	

            tri_size+=pc;

      	    MPI_TOC(time_lanczos_iter, start_time2);

      	    {  //label 3
	            nml_offset eigenvalue = 0;
	            for (eigenvalue = 0; eigenvalue < eigenvalues; ++eigenvalue) {
	  	            nml_offset	j = 0;
	  	            if(conv_vec_method==2){
	     		        if(tri_size>(unsigned int)(num_iteration[eigenvalue]) )
	       		            continue;
	  	            }

	  	            for(int k=0;k<pc;k++){
		                nml_dcvector *q_n= (*Q_n)[k];	

                        #if 0
	  	                for (j = 0; j < length; ++j) {
		                    vector[eigenvalue][j].r += 
						        xm[eigenvalue][tri_size-pc+k]*q_n[j].r-xm_imag[eigenvalue][tri_size-pc+k]*q_n[j].i;
		                    vector[eigenvalue][j].i += 
						        xm[eigenvalue][tri_size-pc+k]*q_n[j].i+xm_imag[eigenvalue][tri_size-pc+k]*q_n[j].r;
						}
						#else
						if(eigenvalue < eigvb)
						{
	  	                    for (j = 0; j < length; ++j) {
		                        vector[eigenvalue][j].r += 
						            xx[eigenvalue][tri_size-pc+k].r*q_n[j].r-xx[eigenvalue][tri_size-pc+k].i*q_n[j].i;
		                        vector[eigenvalue][j].i += 
						            xx[eigenvalue][tri_size-pc+k].r*q_n[j].i+xx[eigenvalue][tri_size-pc+k].i*q_n[j].r;
							}

						}
						else
						{
	  	                    for (j = 0; j < length; ++j) {
		                        vector[eigenvalue][j].r += 
						            xx2[eigenvalue-eigvb][tri_size-pc+k].r*q_n[j].r-xx2[eigenvalue-eigvb][tri_size-pc+k].i*q_n[j].i;
		                        vector[eigenvalue][j].i += 
						            xx2[eigenvalue-eigvb][tri_size-pc+k].r*q_n[j].i+xx2[eigenvalue-eigvb][tri_size-pc+k].i*q_n[j].r;
							}
						}
						#endif
	    	       }//end of k
	            } // end of for eigenvalue
	        }// end of label 3

            MPI_TIC(start_time);
            MPI_Barrier(MPI_COMM_WORLD);
            MPI_TOC(time_barrier, start_time);

            //cout<<"End of iteration:"<<iteration<<endl; 	
        }// end of for iteration

        nml_offset eigenvalue_n;
        nml_offset eigenvalue_n1;
        int rel=0;

        for (eigenvalue_n = 0; eigenvalue_n < eigenvalues; eigenvalue_n++) {
	        for(eigenvalue_n1=eigenvalue_n+1;eigenvalue_n1<eigenvalues;eigenvalue_n1++){
                if((fabs(value[eigenvalue_n1]-value[eigenvalue_n])<tolerance)){
	  	            rel=gramSchmidt((nml_dcvector*) &(vector[eigenvalue_n][0]),(nml_dcvector*) &(vector[eigenvalue_n1][0]),length);
	  	            if(rel==1&&mpi_n3d_id==mpi_n3d_masterid)
	    		        printf("eigenvectors of %1.15f and %1.15f are orthorgonal\n",value[eigenvalue_n1],value[eigenvalue_n]);
	            }// end of if

                MPI_Barrier(MPI_COMM_WORLD);
        }}
    } //end of label2

    if (mpi_n3d_id == mpi_n3d_masterid) {
        if (NULL != fp_log) {
            fprintf(fp_log, "Done with Iterations and construction of eigenvectors\n");
            fflush(fp_log);
    }}

    if(conv_vec_method==2)
        nml_free(num_iteration);
 
    //nml_dm_delete(xm);
    //nml_dm_delete(xm_imag);
    nml_dcm_delete(Alpha);
    nml_dcm_delete(Beta);

    cout<<"End of eigenvectors routine"<<endl;

    return eigenvalues;
} // eigenvectorsBlockLanczos


int eigensolve_blk_lanczos(             // number of eigenvalues actually found 
    nml_dvector  **pValue,      // out real eigenvalue vector           
    int ConvCheckStartIter,     // in convergence check start iteration 
    int ConvCheckSkipRate,      // in convergence check skip rate       
    int          imax,          // in maximum number of iterations + 1  
    int          requested1,    // in requested number of eigenvalues   
    int          requested2,    // in requested number of eigenvalues   
    nml_extent   length,        // in extent of vectors r_n, q_n & q_n1
    int*         argument[],    // matrix-vector multiply argument list 
    nml_dscalar  emin1,         // in eigenvalue search lower bound     
    nml_dscalar  emax1,         // in eigenvalue search upper bound     
    nml_dscalar  emin2,         // in eigenvalue search lower bound    
    nml_dscalar  emax2,         // in eigenvalue search upper bound     
    nml_dscalar  tolerance,     // in eigenvalue convergence tolerance  
    int          verbose,       // in verbose diagnostic messages       
    const char   *fn_trace,     // in   trace log text file name        
    const char   *fn_pendiag,   // in tridiagonal matrix binary file name 
    const char   *fn_eigval,    // in eigenvalue  vector binary file name 
    nml_dcvector *r_n,          // inout r_{n} = q_{n+1}*beta_{n}       
    const nml_dcvector *r_0,    // in    initial guess vector           
    nml_dcvector *q_n,          // inout current  complex Lanczos vector (to be discarded)
    nml_dcvector *q_n1,         // inout previous complex Lanczos vector (to be dsicarded)
    nml_dcvector **pCValue,     //   out complex eigenvalue  vector     
    nml_dcmatrix **pCVector,    //   out complex eigenvector matrix     

	#ifndef BLANCZ_MATMULT 
    void (*matmul)(const int**, nml_dcscalar*, const nml_dcscalar*),
	#else 
    void (*matmul)(const int**, nml_dcmatrix*, const nml_dcmatrix*),
	#endif

    int          restoreState,   // in read tridiagonal matrix from file  (to be discarded) 
    int p,
    double deflation_tolerance,
    int eig_dump,
    int num_vb,
    int num_cb,
    int max_no_success,
    int save_BT_band,
    int save_BT_SMF,
    int read_BT_band,
    int deflation_on,
    int conv_method,
    int conv_vec_method 	  
) 
{
    // sunnylee - test
	//tolerance /= 100;

    int tridiag_size=0; 		//Out tridiagonal size at convergence      
    nml_extent    iterations = 0;
    nml_extent    eigenvalues = 0;
    int found_vb=0;
    int found_cb=0;
   
    FILE          *fp_log = stdout;
    FILE          *fp_trace = (NULL == fn_trace)? NULL: fopen(fn_trace, "a");

    if ((NULL != fn_trace) && (NULL == fp_trace)) 
        error("Could not open the file \"%s\"", fn_trace);

    *pValue     = nml_dv_new(imax*p);                                                               

    //initialize random starting vectors
    vector< nml_dcvector* > R_n;
    vector< nml_dcvector* > Q_n;
    vector< nml_dcvector* > Q_n1;
    vector< nml_dcvector* > Q_n_temp;

    for(int i=0;i<p;i++){
        nml_dcvector* r=nml_dcv_new(length);
        nml_dcvector* q=nml_dcv_new(length);
	    nml_dcvector* q1=nml_dcv_new(length);
	    nml_dcvector* q2=nml_dcv_new(length);

        for(int j=0;j<length;j++){
            r[j].r=1.0*rand()/RAND_MAX;
            r[j].i=0.0;
        }  

        R_n.push_back(q);
        Q_n_temp.push_back(r);
	    Q_n.push_back(q2);
        Q_n1.push_back(q1);
    }

    //Orthonormalize the starting vectors 
    int *indx=new int[p];
    nml_dcmatrix *tri_diag = nml_dcm_new(imax*p, p+1);
    nml_dcmatrix *Beta=nml_dcm_new(p,p);

    //cout<<"Calling BL_QR"<<endl;
    block_Lanczos_QR(&Q_n_temp,&Q_n,Beta,length,p);

    //copy Q_n to Q_n_temp for future use when computing eigenvectors
    for(int i=0;i<p;i++){
	    for(int j=0;j<length;j++){
		    Q_n_temp[i][j].r=Q_n[i][j].r;
		    Q_n_temp[i][j].i=Q_n[i][j].i;
	}}

    //cout<<"Computing eigenvalues"<<endl;

    eigenvalues  = eigenvaluesBlockLanczos(*pValue, tri_diag, &iterations, p,
        &R_n, &Q_n, &Q_n1, length, requested1, requested2,
        ConvCheckStartIter, ConvCheckSkipRate, imax,
        emin1, emax1, emin2, emax2, tolerance, 100.0*tolerance, matmul,
        (const int**)argument, fp_trace, fp_log, verbose, &tridiag_size, 
	    deflation_tolerance, max_no_success, &found_vb, &found_cb, deflation_on, conv_method);

    if(mpi_n3d_id == mpi_n3d_masterid){
	    cout<<"[BLKLANCZ] Iterations :"<<iterations<<endl;
		cout<<"[BLKLANCZ] Eigenvalues:"<<eigenvalues<<endl;
	    cout<<"[BLKLANCZ] Found_vb   :"<<found_vb<<endl;
	    cout<<"[BLKLANCZ] Found_cb   :"<<found_cb<<endl;
    }

    //MPI_Bcast(&eigenvalues, 1, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

    iterations--;

    //print eigenvalues
    if(mpi_n3d_id == mpi_n3d_masterid)
        for(int i=0;i<eigenvalues;i++)
	        cout<<"Converged Eig: "<<i<<"   "<<(*pValue)[i]<<endl;

    //saving block tridiagonal in band format
    if((save_BT_band!=0)&&(mpi_n3d_id == mpi_n3d_masterid)){
        saveBlockTridiagonal(tri_diag, tridiag_size, p);	
        cout<<"Saved Matrix"<<endl;
    }

    saveBlockTridiagonal_MatrixFormat(tri_diag, iterations, p);
    //cout<<"tri_diag"<<endl;
    //print_Matrix( tri_diag, tridiag_size, p+1);

    /*nml_dcmatrix *tri_diag1 = nml_dcm_new(imax*p, p+1);
      readBlockTridiagonal(tri_diag1);
      cout<<"tri_diag1"<<endl;
      print_Matrix( tri_diag1, imax*p, p+1);
    */

    //dump eigenvalues of triad
    if((mpi_n3d_id == mpi_n3d_masterid)&&(eig_dump==1)){
	    nml_dvector *val=nml_dv_new(num_vb+num_cb);
	    int len_vb=0;
        int len_cb=0;
       
	    cout<<"Calling LAPACK eigensolver"<<endl;

	    if(num_vb>0)
	        solve_block_tridiagonal_eigenvalues(tridiag_size, p, tri_diag, val, &len_vb, emin1, emax1, num_vb, tolerance,0);
        if(num_cb>0)
            solve_block_tridiagonal_eigenvalues(tridiag_size, p, tri_diag, val+len_vb, &len_cb, emin2, emax2, num_cb, tolerance,1);

        if((len_vb+len_cb)>0){
	        cout<<"Number of eigenvalues found="<<len_vb+len_cb<<endl;

	        for(int i=0;i<len_vb+len_cb;i++)
        	    cout<<"Eig "<<i<<" :  "<<val[i]<<endl;
				
	        FILE *fp_eig=fopen("Eigenvalues.dat","w");
            if(NULL!=fp_eig){
                for(int i=0;i<len_vb+len_cb;i++)
                    fprintf(fp_eig, "%e\n",val[i]);        
                    fclose(fp_eig);
      	    }else
                cout<<"Error in opening file to dump LAPACK eigenvalues"<<endl;
        }	
		
        nml_dv_delete(val);	
    } // eig_dump
	
    MPI_Bcast(&tridiag_size, 1, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

    //compute eigenvectors
    if (NULL != r_0) {
 	    int num_vb1=found_vb;
        int num_cb1=found_cb;
		double delta = 0.0001;

        int req=-1;

	    double range1_min=emin1 - delta;
        double range1_max=emax1 + delta;
        double range2_min=emin2 - delta;
        double range2_max=emax2 + delta;
	
        #if 0
	    if(mpi_n3d_id == mpi_n3d_masterid){
	        double delta=0.0001;	
	        cout<<"Calling LAPACK eigenvector solver"<<endl;

       	    if(found_vb>0){
		        range1_min=(*pValue)[0]-delta;
		        range1_max=(*pValue)[found_vb-1]+delta;

		        cout<<"range1 min:"<<range1_min<<endl;
	            cout<<"range1 max:"<<range1_max<<endl;

	            solve_block_tridiagonal_eigenvalues(tridiag_size, p, tri_diag, val1, 
				                    &num_vb1, range1_min, range1_max, -1, tolerance);

                //MPI_Bcast(&num_vb1, 1, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

		        cout<<"Number found in VB="<<num_vb1<<endl;

		        for(int i=0;i<num_vb1;i++)
			        cout<<val1[i]<<endl;
	        }

	        if(found_cb>0){
                range2_min=(*pValue)[found_vb]-delta;
                range2_max=(*pValue)[found_vb+found_cb-1]+delta;

		        cout<<"range min:"<<range2_min<<endl;
	            cout<<"range max:"<<range2_max<<endl;

		        solve_block_tridiagonal_eigenvalues(tridiag_size, p, tri_diag, val2, &num_cb1,
					                range2_min, range2_max, -1, tolerance);

                //MPI_Bcast(&num_cb1, 1, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

		        cout<<"Number found in CB="<<num_cb1<<endl;

		        for(int i=0;i<num_cb1;i++)
                    cout<<val2[i]<<endl;
            }// end of cb
        } //mpi_n3d
        #endif

        nml_dcmatrix *xx;  //= nml_dcm_new(num_vb1, tridiag_size); // eig_vectors for VB
	    nml_dcmatrix *xx2; //= nml_dcm_new(num_cb1, tridiag_size); // eig_vectors for CB
	    nml_dvector  *val; //= nml_dv_new(num_vb1+num_cb1);        // eig_values for VB+CB

		if(mpi_n3d_id == mpi_n3d_masterid)
		{
		    nml_dvector *val1 = nml_dv_new(tridiag_size);
		    nml_dvector *val2 = nml_dv_new(tridiag_size);
            xx  = nml_dcm_new(num_vb1, tridiag_size); // eig_vectors for VB
	        xx2 = nml_dcm_new(num_cb1, tridiag_size); // eig_vectors for CB
	        val = nml_dv_new(num_vb1+num_cb1);        // eig_values for VB+CB

	        if(num_vb1>0){
		        for(int i=0;i<num_vb1;i++){
			        val1[i] = (*pValue)[i];
		    }}

	        if(num_cb1>0){
                for(int i=num_vb1;i<num_cb1+num_vb1;i++){
				    val2[i-num_vb1] = (*pValue)[i];
            }}
        
	        if(found_vb>0){
                //MPI_Barrier(MPI_COMM_WORLD);

                double  elapsed_time1=-MPI_Wtime();
                solve_block_tridiagonal_eigenvectors(tridiag_size,p,tri_diag, val1, &num_vb1,
                            range1_min, range1_max, -1 , tolerance,xx,0);
                elapsed_time1+=MPI_Wtime();
                cout<<"Time taken to solve for Ritz vectors in VB range:"<<elapsed_time1<<endl;
	        }

	        if(found_cb>0){
                //MPI_Barrier(MPI_COMM_WORLD);

                double  elapsed_time2=-MPI_Wtime();
                solve_block_tridiagonal_eigenvectors(tridiag_size,p,tri_diag, val2, &num_cb1,
                            range2_min, range2_max, -1, tolerance,xx2,1);
                elapsed_time2+=MPI_Wtime();
                cout<<"Time taken to solve for Ritz vectors in CB range:"<<elapsed_time2<<endl;
	        }

	        eigenvalues = num_vb1 + num_cb1;

	        if(num_vb1>0){
		        for(int i=0;i<num_vb1;i++){
			        val[i]=val1[i];
		    }}

	        if(num_cb1>0){
                for(int i=num_vb1;i<num_cb1+num_vb1;i++){
                    val[i]=val2[i-num_vb1];
            }}

            nml_dv_delete(val1);
            nml_dv_delete(val2);
        } // mpi_master...
		
        nml_dcm_delete(tri_diag);

        if(mpi_n3d_id == mpi_n3d_masterid){
            cout<<"Saving eigenvalue-vector map to file"<<endl;

            FILE *fp_eig1=fopen("Converged_Eig.dat","w");
            if(NULL!=fp_eig1){
                for(int i=0;i<found_vb+found_cb;i++)
                    fprintf(fp_eig1, "%e\n",(*pValue)[i]);

                fclose(fp_eig1);
            }else
                cout<<"Error in opening file to dump eigenvalue-vector map"<<endl;
        }


	    for(int i=0;i<p;i++){
            for(int j=0;j<length;j++){
		        R_n[i][j].r=0.0;
		        R_n[i][j].i=0.0;
                Q_n[i][j].r=Q_n_temp[i][j].r;
                Q_n[i][j].i=Q_n_temp[i][j].i;
		        Q_n1[i][j].r=0.0;
		        Q_n1[i][j].i=0.0;
        }}

	    // eigenvector solve
    	// Set up the initial guess for the rerun of the Lanczos algorithm. 
		// I need a barrier here to synchronize
		MPI_Barrier(MPI_COMM_WORLD);

        MPI_Bcast(&num_vb1, 1, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);
        MPI_Bcast(&num_cb1, 1, MPI_INT, mpi_n3d_masterid, MPI_COMM_WORLD);

        printf("[BLKLANCZ] CPU: %d, eigenvalues: %d (VB:%d, CB:%d) tridiag_size = %d\n", mpi_n3d_id, num_cb1+num_vb1, num_vb1, num_cb1, tridiag_size);

    	*pCValue  = nml_dcv_new(num_cb1+num_vb1);
    	*pCVector = nml_dcm_new(num_cb1+num_vb1, length);

		if(mpi_n3d_id != mpi_n3d_masterid)
		{
		    if(num_vb1 > 0)
                xx  = nml_dcm_new(num_vb1, tridiag_size); // eig_vectors for VB

			if(num_cb1 > 0)
	            xx2 = nml_dcm_new(num_cb1, tridiag_size); // eig_vectors for CB

	        val = nml_dv_new(num_cb1+num_vb1);        // eig_values for VB+CB

			if(xx == NULL || xx2 == NULL || val == NULL)
			    printf("[BLKLANCZ] Memory not alllocated!!!\n");
		}
        
		#if 0 
	    if(num_vb1>0){
		    for(int i=0;i<num_vb1;i++){
			    val[i]=val1[i];
		}}

	    if(num_cb1>0){
            for(int i=num_vb1;i<num_cb1+num_vb1;i++){
                val[i]=val2[i-num_vb1];
        }}
		#endif

		MPI_Bcast(&val[0], num_vb1+num_cb1, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);

        for(int i = 0 ; i < num_vb1+num_cb1 ; i++)
	        (*pValue)[i]=val[i];

        if(num_vb1 > 0)
            MPI_Bcast(&xx[0][0], 2*num_vb1*tridiag_size, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);

		if(num_cb1 > 0)
            MPI_Bcast(&xx2[0][0], 2*num_cb1*tridiag_size, MPI_DOUBLE, mpi_n3d_masterid, MPI_COMM_WORLD);

        //MPI_Barrier(MPI_COMM_WORLD);
		
		printf("[BLKLANCZ] CPU: %d Calling Eigenvector Solver\n", mpi_n3d_id);

    	eigenvectorsBlockLanczos(*pCValue, *pCVector, xx, xx2, val, num_vb1+num_cb1, num_vb1, num_cb1,
        	iterations, tridiag_size, tri_diag, &R_n, &Q_n, &Q_n1, length, p, matmul,
        	(const int**)argument, tolerance,fp_trace, fp_log, verbose, deflation_tolerance, deflation_on, conv_vec_method);

        cout<<"return from eigenvector solver: about to finsih"<<endl;

        nml_dv_delete(val);

        if(num_vb1>0)
            nml_dcm_delete(xx);

        if(num_cb1>0)
            nml_dcm_delete(xx2);

		eigenvalues = num_vb1 + num_cb1;

        cout<<"Exit from eigensolver"<<endl; 
    } // end of eigenvector calculation

	nml_dcm_delete(Beta);
	// sunnylee - test
	//nml_dcm_delete(tri_diag);

    cout<<"Successful Exit"<<endl;
 
    return eigenvalues;
} // eigensolve_blk_lanczos     
                                    

