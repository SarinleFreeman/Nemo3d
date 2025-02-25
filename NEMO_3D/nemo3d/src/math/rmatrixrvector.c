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
$Header: /repo/nemo3d/src/math/rmatrixrvector.c,v 1.3 2004/10/27 19:43:10 gekco
Exp $
*****************************************************************************/

/*T_HEADER_FILE{rmatrixrvector.c}
  This file contains creator and destructor functions for
  three dimensional real matrices. Sparse memory allocation is enabled
  by addressing though a double matrix/vector  scheme.  Only the lower lying
  rvectors can be handed down to fortran. T*/

#include "rmatrixrvector.h"

/*T_HEADER_SHORT{Rmatrixrvectr}
  Create the memory and pointers to a 3 dimensional
  matrix of vectors.  Do not allocate memory for the lower 1 dimension.
T*/

rmatrixrvectr Rmatrixrvectr(int d1, int d2) {
  int i, j;
  rmatrixrvectr t;

  /* Allocate memory for poiters to pointers to rows. */
  t = (real ***)nml_malloc((size_t)((d1 + 2) * sizeof(real **)));
  t += 2;

  vdim1(t) = d1;
  vdim2(t) = d2;

  /* Allocate memory of pointers to rows and set pointers to them. */
  t[0] = (real **)nml_malloc((size_t)((d1 * d2) * sizeof(real *)));

  for (i = 1; i < d1; i++)
    t[i] = t[i - 1] + d2;

  for (i = 0; i < d1; i++) {
    for (j = 0; j < d2; j++) {
      t[i][j] = NULL;
    }
  }

  /* Return the pointer. */

  return t;
}

/*T_HEADER_SHORT{rm_rmatrixrvectr}
  Free up memory assigned to the Rmatrixrvectr structure.
   T*/
void rm_rmatrixrvectr(rmatrixrvectr *t_ptr) {
  rmatrixrvectr t = *t_ptr;
  int d1, d2;
  int i1, i2;

  if (!t)
    return;

  d1 = vdim1(t);
  d2 = vdim2(t);

  for (i1 = 0; i1 < d1; i1++)
    for (i2 = 0; i2 < d2; i2++) {
      rm_rvectr(&t[i1][i2]);
      t[i1][i2] = NULL;
    }
  nml_free(t[0]);
  nml_free(t - 2);

  *t_ptr = NULL;
  return;
}

/*T_HEADER_SHORT{rm_rmatrixrvectr_shell_only}
  Free up memory assigned to the Rmatrixrvectr structure. Do not free the memory
  of the submatrices.
   T*/
void rm_rmatrixrvectr_shell_only(rmatrixrvectr *t_ptr) {
  rmatrixrvectr t = *t_ptr;

  if (!t)
    return;

  nml_free(t[0]);
  nml_free(t - 2);

  *t_ptr = NULL;
  return;
}

/* Attach a data cmatrix to the holding structure.  Do not create memory for d
   the data vector itself.  */
int add_elem_rmatrixrvectr(rvectr data_vector, int pos1, int pos2,
                           rmatrixrvectr hold_matrix) {
  if (pos1 > vdim1(hold_matrix)) {
    die("Error in aadd_elem_rmatrixrvectr\n"
        "Requested Position 1 in the matrix of vectors\n"
        "exceeds holding capacity.\n"
        "pos=%d dim=%d\n",
        pos1, vdim1(hold_matrix));
    return -1;
  }
  if (pos2 > vdim2(hold_matrix)) {
    die("Error in aadd_elem_rmatrixrvectr\n"
        "Requested Position 2 in the matrix of vectors\n"
        "exceeds holding capacity.\n"
        "pos=%d dim=%d\n",
        pos2, vdim2(hold_matrix));
    return -1;
  }

  hold_matrix[pos1][pos2] = data_vector;

  return 0;
}

#ifdef RMATRIXRVECTR_TEST
int main(int argc, char *argv[]) {

  int i1, i2, n1 = 2, n2 = 3, i3, i4, n3 = 4;
  rmatrixrvectr C = Rmatrixrvectr(n1, n2);
  rvectr sample = NULL;

  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {
      sample = Rvectr(n3);
      for (i3 = 0; i3 < n3; i3++) {
        printf("(i1=%d,i2=%d) (i3=%d)  ", i1, i2, i3);
        sample[i3] = i1 * 100.0 + i2 * 10.0 + i3 * 1.0;
        printf("%4d  ", (int)sample[i3]);
      }
      add_elem_rmatrixrvectr(sample, i1, i2, C);
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");

  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2; i2++) {
      for (i3 = 0; i3 < cdim(C[i1][i2]); i3++) {
        printf("(i1=%d,i2=%d) (i3=%d)  ", i1, i2, i3);
        printf("%4d  ", (int)C[i1][i2][i3]);
        printf("\n");
      }
      printf("\n");
    }
    printf("\n");
  }
  printf("\n");

  memory_report(stdout);
  rm_rmatrixrvectr(&C);
  memory_report(stdout);

  return 0;
}
#endif
