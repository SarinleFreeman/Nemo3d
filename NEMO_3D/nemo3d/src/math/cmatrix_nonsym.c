/*****************************************************************************
The Jet Propulsion Laboratory (JPL) NanoElectronicMOdeling-3D package.
Copyright (C) 2002 California Institute of Technology (Caltech)

This application is free software, which you can redistribute and/or modify
under the terms of the GNU Lesser General Public License as published by the
Free Software Foundation; either version 2.1 of the License, or (at your
option) any later version.

This library is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this library; see the file COPYING. If not, write to the
Free Software Foundation, Inc.,
59 Temple Place, Suite 330,
Boston, MA  02111-1307  USA

For additional information, please contact
  Gerhard Klimeck (gekco@jpl.nasa.gov)
  Fabiano Oyafuso (fabiano@jpl.nasa.gov)

Written by:  Chris Bowen
             Gerhard Klimeck
             Fabiano Oyafuso
             Seungwon Lee
             Olga Lazarenkova
             Hook Hua

This product includes software developed by the Apache Software Foundation
(http://www.apache.org/).

*****************************************************************************
$Header: /repo/nemo3d/src/math/cmatrix_nonsym.c,v 1.9 2005/06/02 18:46:26 marek Exp $
*****************************************************************************/

/*
 * Contains eigen-routines for complex non-symmetric Hermitian matrices.
 */

#include "cmatrix_nonsym.h"



/**
 * Complex comparator used for passing into qsort().
 */
int complexCompare(const void *complexA, const void *complexB) {
    if      ( ((complex*)complexA)->r <  ((complex*)complexB)->r ) return -1;
    else if ( ((complex*)complexA)->r >  ((complex*)complexB)->r ) return  1;
    else return  0;
}



/**
 * Computes the eigenvalues of a complex non-symmetric cmatrix 'h'
 * and returns the eigenvalues in the cvector 'cvEigenValue'.
 */
void ceigval_ns_full( cvectr cvEigenValue, cmatrix h ) {
    int n, i, j;
    int ldh, ldz, ilo, ihi;
    cvectr tau, work;
    int lwork, info;
    char job, compz;
    cmatrix z;
    complex temp; /* Used in Swap() */


    /* --------------------------------------------------------------------------
     * Reduction to Hessenberg form by elimination method.
     * Cmatrix 'h' is a complex nonsymetric matrix, which gets replaced
     * by an upper Hessenberg matrix with the same eigenvalues.
     */
    //    n = rdim(h);
    //  ldh = cdim(h);

    n = (int) nml_dcm_extent1( h ) ;
    ldh = (int) nml_dcm_extent2( h ) ;

    tau = Cvectr(n-1);
    lwork = 5120; /* Optimal lwork size */
    work = Cvectr(lwork);

    ilo = 1;
    ihi = n;

    /* Transpose matrix */
    for (i=0; i < n; i++ ) {
        for (j=0; j < i; j++ ) {
            Swap( h[i][j], h[j][i] );
        }
    }


    /* Call to LAPACK
     * ZGEHRD reduces a complex general matrix to upper Hessenberg form by a
     * unitary similarity transformation.
     * Note: In EISPACK, function is COMHES() (non-unitary) and CORTH() (unitary).
     * http://www.netlib.org/lapack/complex16/zgehrd.f
     */
    zgehrd(&n, &ilo, &ihi, &h[0][0], &ldh, tau, work, &lwork, &info);


    if (info == 0) {
        /*
        printf("In ceigval_ns_full():: Hessenberg Elimination -> Successful exit from LAPACK's zgehrd(). Optimal LWORK is %f.\n", work[0].r );
        */
    }
    else if (info < 0) {
        n3d_warning("In ceigval_ns_full():: Hessenberg Elimination -> From call to LAPACK's zgehrd(), the %d-th argument had an illegal value.\n", abs(info) );
    }

    rm_cvectr(&work);
    rm_cvectr(&tau);


    /* --------------------------------------------------------------------------
     * All eigenvalues of a complex upper Hessenberg matrix by the QR algorithm.
     * Reduction to Hessenberg form by elimination method.
     * Cmatrix 'h' is a complex a complex upper Hessenberg matrix.
     * The contents of 'h' are unspecified on exit.
     */
    job = 'E';   /* compute eigenvalues only */
    compz = 'N'; /* no Schur vectors are computed; */
    //    n = rdim(h);

    n = (int) nml_dcm_extent1( h ) ;
    ilo = 1;
    ihi = n;

    //    ldh = cdim(h); /* trailing dimension in C/C++ */
    //    ldz = cdim(h); /* trailing dimension in C/C++ */

    ldh = (int) nml_dcm_extent2( h ) ;
    ldz = ldh ;
    z = Cmatrix(n,ldz); /* not referenced in LAPACK if compz = 'N' */
    lwork = 160; /* Optimal lwork size */
    work = Cvectr(lwork);

    /* Call to LAPACK
     * ZHSEQR computes the eigenvalues of a complex upper Hessenberg matrix H. 
     * Note: In EISPACK, function is COMQR()
     * http://www.netlib.org/lapack/complex16/zhseqr.f
     */
    zhseqr(&job, &compz, &n, &ilo, &ihi, &h[0][0], &ldh, cvEigenValue, &z[0][0], &ldz, work, &lwork, &info);

    if (info == 0) {
        /*
        printf("In ceigval_ns_full():: QR Algorithm -> Successful exit from LAPACK's zhseqr(). Optimal LWORK is %f.\n", work[0].r );
        */
    }
    else if (info < 0) {
        n3d_warning("In ceigval_ns_full():: QR Algorithm -> From call to LAPACK's zhseqr(), the %d-th argument had an illegal value.\n", abs(info) );
    }
    else if (info > 0) {
        n3d_warning("In ceigval_ns_full():: QR Algorithm -> From call to LAPACK's zhseqr(), ZHSEQR failed to compute all the eigenvalues in a total of 30*(IHI-ILO+1) iterations.\n");
    }

    rm_cvectr(&work);
    rm_cmatrix(&z);


    /* Sort the eigenvalues based on the real part */
    qsort( (void*)cvEigenValue, n, sizeof(complex), complexCompare );


    /* Transpose matrix
    for (i=0; i < n; i++) {
        for (j=0; j < i; j++) {
            Swap( h[i][j], h[j][i] );
        }
    }
    */

#ifdef WRITE_EIGENVALUES_TO_FILE
    { /* Write eigenvalues to file */
        FILE* fp = fopen("ham-eigenvalues","w");
        for (int i=0; i<n; i++) {
            fprintf(fp, "%d %f %f\n", i, cvEigenValue[i].r, cvEigenValue[i].i);
        }
    }
#endif
}




/**
 * Computes the eigenvalues and eigenvectors of input cmatrix 'h'. 
 */
void ceigvec_ns_full( cvectr cvEigenValue, cmatrix cmEigenVector, cmatrix h ) {
  // int n = rdim(h);
  int n = (int) nml_dcm_extent1( h ) ;

    int i, j;
    int lda;
    cvectr work;
    int lwork, info;
    char jobz, uplo;
#ifdef C99_AUTO
    double rwork[3*n-2], w[n];
#else
    double *rwork = (double *) malloc(sizeof(double)*3*n-2*sizeof(double));
    double *w = (double *) malloc(sizeof(double)*n);
#endif
    complex temp; /* Used in Swap() */


    /* Copy cmatrix h to cmEigenVector */
    for (i=0; i<n; i++) {
        for (j=0; j<n; j++) {
            cmEigenVector[i][j] = h[i][j];
        }
    }

    /* Transpose matrix */
    for (i=0; i<n; i++) {
        for (j=0; j<i; j++) {
            Swap( cmEigenVector[i][j], cmEigenVector[j][i] );
        }
    }


    /* Call to LAPACK
     * ZHEEV computes all eigenvalues and, optionally, eigenvectors of a
     * complex Hermitian matrix A.
     * http://www.netlib.org/lapack/complex16/zheev.f
     */
    jobz = 'V'; /* Compute eigenvalues and eigenvectors. */
    uplo = 'U'; /* Upper triangle of A is stored. */
    //    lda = cdim(h);

    lda = (int) nml_dcm_extent2( h ) ;

    lwork = 2 * lda - 1 ;
    /*    lwork = 5280;  Optimal lwork size */
    work = Cvectr(lwork);
    zheev(&jobz, &uplo, &n, &cmEigenVector[0][0], &lda, w, work, &lwork, rwork, &info);
    if (info == 0) {
        /*
        printf("In ceigvec_ns_full():: Successful exit from LAPACK's zheev(). Optimal LWORK is %f.\n", work[0].r );
        */
    }
    else if (info < 0) {
        n3d_warning("In ceigvec_ns_full():: From call to LAPACK's zheev(), the %d-th argument had an illegal value.\n", abs(info) );
    }
    else if (info > 0) {
        n3d_warning("In ceigvec_ns_full():: From call to LAPACK's zheev(), the algorithm failed to converge, %d off-diagonal elements of an intermediate tridiagonal form did not converge to zero.\n", info);
    }
    rm_cvectr(&work);


    /* Copy the real eigenvalues set in w[] to the real part of the cvector cvEigenValue. */
    for (i=0; i<n; i++) {
        cvEigenValue[i].r = w[i];
    }


    
    /* Transpose cmatrix cmEigenVector */
    for (i=0; i<n; i++) {
        for (j=0; j<i; j++) {
            Swap( cmEigenVector[i][j], cmEigenVector[j][i] );
        }
    }


    { /* Write eigenvectors to file */
        FILE* fp = fopen("ham-eigenvectors","w");
        for (i=0; i<n; i++) {
            for (j=0; j<n; j++) {
                if (cmEigenVector[i][j].r!=0.0 || cmEigenVector[i][j].i!=0.0) {
                    fprintf(fp, "%d %d %f %f\n", i, j, cmEigenVector[i][j].r, cmEigenVector[i][j].i);
                }
            }
        }
    }


    { /* Write eigenvalues to file */
        FILE* fp = fopen("ham-eigenvalues","w");
        for (i=0; i<n; i++) {
            fprintf(fp, "%d %f %f\n", i, cvEigenValue[i].r, cvEigenValue[i].i);
        }
    }
#ifndef C99_AUTO
  free(rwork);
  free(w);
#endif
}
