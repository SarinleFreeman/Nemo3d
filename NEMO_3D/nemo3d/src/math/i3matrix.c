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
$Header: /repo/nemo3d/src/math/i3matrix.c,v 1.3 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

/*T_HEADER_FILE{i3matrix.c}
  This file contains creator and destructor functions for
  three-dimensional integer matrices.
  The mathematical objects are constructed such that they
  can be passed into FORTRAN and F90 by passing the address of the very first
  object element.  Then the multi-dimensional object can be dealt with
  in the FORTRAN or F90 code by the usual matrix addressing scheme.  However the
  ordering of the indices is inverted from the usual C notation. Within the C
  code the starting and ending point of each dimension is stored and retrieved
  with macros like \verb|vdim1| REFERENCE{vdim1,i3matrix.h}. T*/

#include "i3matrix.h"

/*T_HEADER_SHORT{I3matrix}
  Create the memory and pointer to a integer 3 dimensional matrix of dimensions
\\
  \verb|t[0...nrow-1][0...ncol-1][0...ndep-1]| and return the pointer to it.
  The dimensions of the matrix are stored using the following macros:
  \verb|vdim1| REFERENCE{vdim1,i3matrix.h},
  \verb|vdim2| REFERENCE{vdim2,i3matrix.h}, and
  \verb|vdim3| REFERENCE{vdim3,i3matrix.h}.
T*/
i3matrix I3matrix(int nrow, int ncol, int ndep) {
  int i, j, k, n1 = nrow - 1, n2 = ncol - 1, n3 = ndep - 1;
  i3matrix t;

  t = (int ***)nml_malloc((size_t)((nrow + 3) * sizeof(int **)));
  if (!t)
    fprintf(LOGFILE, "I3matrix unable to allocate memory on stage 1.\n");
  t += 3;
  vdim1(t) = nrow;
  vdim2(t) = ncol;
  vdim3(t) = ndep;

  t[0] = (int **)nml_malloc((size_t)((nrow * ncol) * sizeof(int *)));
  if (!t[0])
    fprintf(LOGFILE, "I3matrix unable to allocate memory on stage 2.\n");

  t[0][0] = (int *)nml_malloc((size_t)((nrow * ncol * ndep) * sizeof(int)));
  if (!t[0][0])
    fprintf(LOGFILE, "I3matrix unable to allocate memory on stage 3.\n");

  for (j = 0 + 1; j <= n2; j++)
    t[0][j] = t[0][j - 1] + ndep;
  for (i = 0 + 1; i <= n1; i++) {
    t[i] = t[i - 1] + ncol;
    t[i][0] = t[i - 1][0] + ncol * ndep;
    for (j = 0 + 1; j <= n2; j++)
      t[i][j] = t[i][j - 1] + ndep;
  }

  for (i = 0; i <= n1; i++) {
    for (j = 0; j <= n2; j++) {
      for (k = 0; k <= n3; k++) {
        t[i][j][k] = 0;
      }
    }
  }

  return t;
}

/*T_HEADER_SHORT{rm_i3matrix}
  Free up memory assigned to the I3matrix structure.
  The lower bound of the first dimension must be given as an input for the
  proper destruction of the object. T*/
void rm_i3matrix(i3matrix *t_ptr) {
  i3matrix t = *t_ptr;

  if (!t)
    return;

  nml_free(t[0][0]);
  nml_free(t[0]);
  nml_free(t - 3);

  *t_ptr = NULL;
  return;
}

/*T_HEADER_SHORT{check_i3matrix}
  Check whether or not the target \verb|i3matrix|
  exists and if it exists, whether is has the right
  dimension.  If the target does not exist, we create it here from scratch.
  If it exists, we need to check whether the existing dimensions are correct.
  If they are correct we just return the structure.  If they are incorrect,
  we remove the existing structure and make a new one of the correct dimension.
  \\ Note that there is a restriction to the reshaping of the i3matrix:  The
  first index cannot be reshaped!
  T*/
i3matrix check_i3matrix(i3matrix t, int nrow, int ncol, int ndep) {
  if (!t)
    return I3matrix(nrow, ncol, ndep);

  /* Check whether all the dimensions are correct. */
  if (nrow == vdim1(t) && ncol == vdim2(t) && ndep == vdim3(t)) {
  } else {
    rm_i3matrix(&t);
    t = I3matrix(nrow, ncol, ndep);
  }

  return t;
}

#ifdef I3MATRIX_TEST
int main(int argc, char *argv[]) {

  int i, ii, iii, n1 = 10, n2 = 20, n3 = 30;
  i3matrix I = I3matrix(n1, n2, n3);

  for (i = 0; i < n1; i++) {
    for (ii = 0; ii < n2; ii++) {
      for (iii = 0; iii < n3; iii++) {
        I[i][ii][iii] = (i + 1) * (ii + 1) + (iii + 1);
      }
    }
  }
  for (i = 0; i < n1; i++) {
    for (ii = 0; ii < n2; ii++) {
      for (iii = 0; iii < n3; iii++) {
        printf(" %d", I[i][ii][iii]);
      }
      printf("\n");
    }
    printf("\n");
  }

  memory_report(stdout);
  rm_i3matrix(&I);
  memory_report(stdout);

  return 0;
}
#endif
