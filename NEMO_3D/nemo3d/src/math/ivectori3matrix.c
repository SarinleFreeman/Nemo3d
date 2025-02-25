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
$Header: /repo/nemo3d/src/math/ivectori3matrix.c,v 1.3 2004/10/27 19:43:10 gekco
Exp $
*****************************************************************************/

#include "ivectori3matrix.h"

/* Create memory to a vector of Ivector */
ivectri3matrix Ivectri3matrix(int n) {
  int i;
  ivectri3matrix I = (i3matrix *)nml_calloc((n + 1), sizeof(i3matrix));
  if (I == NULL) {
    die("unable to allocate memory in Ivectri3matrix.\n");
  }
  I++;
  vdim(I) = n;

  for (i = 0; i < n; i++) {
    I[i] = NULL;
  }

  return I;
}

/* Remove the memory associated with the vector to Ivectr. */
void rm_ivectri3matrix(ivectri3matrix *I_ptr) {
  ivectri3matrix I;
  int i;

  if (!I_ptr || !(I = *I_ptr))
    return;

  for (i = 0; i < vdim(I); i++) {
    rm_i3matrix(&I[i]);
    I[i] = NULL;
  }

  nml_free(--I);

  *I_ptr = NULL;

  return;
}

/* Remove the memory associated with the vector to Ivectr.
   Only remove the memory associated with the
   holding vector; not the datavectors. */
void rm_ivectri3matrix_shellonly(ivectri3matrix *I_ptr) {
  ivectri3matrix I;

  if (!I_ptr || !(I = *I_ptr))
    return;

  nml_free(--I);

  *I_ptr = NULL;

  return;
}

/* Attach a data vector to the holding structure.  Do not create memory for d
   the data vector itself.  */
int add_elem_ivectri3matrix(i3matrix data_matrix, int pos,
                            ivectri3matrix hold_vectr) {
  if (pos > vdim(hold_vectr)) {
    die("Error in add_elem_ivectri3matrix\n"
        "Requested Position in the vector of vectors\n"
        "exceeds holding capacity.\n"
        "pos=%d dim=%d\n",
        pos, vdim(hold_vectr));
    return -1;
  }

  if (hold_vectr[pos] != NULL) {
    rm_i3matrix(&hold_vectr[pos]);
  }

  hold_vectr[pos] = data_matrix;

  return 0;
}

#ifdef IVECTI3MATRIX_TEST
int main(int argc, char *argv[]) {

  int i1, i2, i3, i4, n1 = 3, n2 = 5, n3 = 7, n4 = 3;
  ivectri3matrix I = Ivectri3matrix(n1);
  i3matrix sample = NULL;

  for (i1 = 0; i1 < n1; i1++) {
    sample = I3matrix(n2 + i1, n3 + i1, n4 + i1);
    add_elem_ivectri3matrix(sample, i1, I);
    for (i2 = 0; i2 < n2 + i1; i2++)
      for (i3 = 0; i3 < n3 + i1; i3++)
        for (i4 = 0; i4 < n4 + i1; i4++) {
          sample[i2][i3][i4] = i1 + i2 + i3 + i4;
          printf("(%d %d %d %d) = %d, %d\n", i1, i2, i3, i4, sample[i2][i3][i4],
                 I[i1][i2][i3][i4]);
        }
  }
  printf("\n");

  for (i1 = 0; i1 < n1; i1++) {
    for (i2 = 0; i2 < n2 + i1; i2++)
      for (i3 = 0; i3 < n3 + i1; i3++)
        for (i4 = 0; i4 < n4 + i1; i4++) {
          printf("(%d %d %d %d) = %d\n", i1, i2, i3, i4, I[i1][i2][i3][i4]);
        }
  }
  printf("\n");

  memory_report(stdout);
  rm_ivectri3matrix(&I);
  memory_report(stdout);

  return 0;
}
#endif
