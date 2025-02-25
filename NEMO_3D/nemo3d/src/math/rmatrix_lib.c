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
$Header: /repo/nemo3d/src/math/rmatrix_lib.c,v 1.6 2005/03/04 21:07:51 dseaman
Exp $
*****************************************************************************/

#include "rmatrix_lib.h"

/*
 * Allocate memory for a square real operator matrix of dimension n.
 * Returns the pointer from Rmatrix(n,n)
 */
roperator Roperator(int n) { return Rmatrix(n, n); }

/*
 * 1) Computes a LU decomposition of a real square roperator (matrix) "a0"
 *    Store the LU decomposed rmatrix as roperator "a".
 * 2) Solves the system "Ax = b" for x. A is the LU decomposed square real
 * matrix. Returns the rvector x
 */
rvectr luDecompositionBackSubstitution(roperator a0, roperator a, rvectr x,
                                       rvectr b) {

  int i, j, m, n, nrhs, lda, ldb, info;
  char trans;
#ifdef C99_AUTO
  int ipiv[rdim(a0)];
#else
  int *ipiv = (int *)malloc(sizeof(int) * n);
#endif

  m = cdim(a0);   /* trailing dimension in C/C++ */
  n = rdim(a0);   /* leading  dimension in C/C++ */
  lda = rdim(a0); /* leading  dimension in C/C++ */

  /* Transpose-copy matrix */
  for (i = 0; i < n; i++) {
    for (j = 0; j < n; j++) {
      a[j][i] = a0[i][j];
    }
  }

  /* Call to LAPACK
   * DGETRF computes an LU factorization of a general M-by-N matrix A
   * using partial pivoting with row interchanges.
   * The factorization has the form
   *     A = P * L * U
   * where P is a permutation matrix, L is lower triangular with unit
   * diagonal elements (lower trapezoidal if m > n), and U is upper
   * triangular (upper trapezoidal if m < n).
   * http://www.netlib.org/lapack/double/dgetrf.f
   */
  dgetrf(&m, &n, &a[0][0], &lda, (int **)&ipiv, &info);

  if (info == 0) {
    /*
    printf("In luDecompositionBackSubstitution():: LU factorization ->
    Successful exit from LAPACK's dgetrf().\n");
    */
  } else if (info < 0) {
    n3d_warning(
        "In luDecompositionBackSubstitution():: LU factorization -> From call "
        "to LAPACK's dgetrf(), the %d-th argument had an illegal value.\n",
        abs(info));
  } else if (info > 0) {
    n3d_warning(
        "In luDecompositionBackSubstitution():: LU factorization -> From call "
        "to LAPACK's dgetrf(), U(i,i) is exactly zero. The factorization has "
        "been completed, but the factor U is exactly singular, and division by "
        "zero will occur if it is used to solve a system of equations.\n");
  }

  /* Copy rvector b to rvector x,
     operations occur on rvector x so as to preserve rvector b. */
  cp_rvectr(x, b);

  trans = 'N'; /*  A * X = B  (No transpose A') */
  nrhs = 1;    /* The number of right hand sides columns of the matrix B. (1 =>
                  column vector) */
  ldb = nml_dv_extent(b);

  /* Call to LAPACK
   * DGETRS solves a system of linear equations
   *     A * X = B  or  A' * X = B
   * with a general N-by-N matrix A using the LU factorization computed by
   * DGETRF. http://www.netlib.org/lapack/double/dgetrs.f
   */
  dgetrs(&trans, &n, &nrhs, &a[0][0], &lda, (int **)&ipiv, x, &ldb, &info);

  if (info == 0) {
    /*
    printf("In luDecompositionBackSubstitution():: BackSubstitution ->
    Successful exit from LAPACK's dgetrs().\n");
    */
  } else if (info < 0) {
    n3d_warning(
        "In luDecompositionBackSubstitution():: BackSubstitution -> From call "
        "to LAPACK's dgetrs(), the %d-th argument had an illegal value.\n",
        abs(info));
  }
#ifndef C99_AUTO
  free(ipiv);
#endif
  return x;
}
