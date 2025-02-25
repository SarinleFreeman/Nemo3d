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
$Header: /repo/nemo3d/src/math/ivectorvector.c,v 1.4 2004/12/02 12:09:56 marek
Exp $
*****************************************************************************/

#include "ivectorvector.h"

/* Create memory to a vector of Ivector */
ivectrvectr Ivectrvectr(int n) {
  int i;
  ivectrvectr I = (ivectr *)nml_calloc((n + 1), sizeof(ivectr));
  if (I == NULL) {
    die("unable to allocate memory in Ivectrvectr.\n");
  }
  I++;
  vdim(I) = n;

  for (i = 0; i < n; i++) {
    I[i] = NULL;
  }

  return I;
}

/* Remove the memory associated with the vector to Ivectr. */
void rm_ivectrvectr(ivectrvectr *I_ptr) {
  ivectrvectr I;
  int i;

  if (!I_ptr || !(I = *I_ptr))
    return;

  for (i = 0; i < vdim(I); i++) {
    rm_ivectr(&I[i]);
    I[i] = NULL;
  }

  nml_free(--I);

  *I_ptr = NULL;

  return;
}

/* Remove the memory associated with the vector to Ivectr.
   Only remove the memory associated with the
   holding vector; not the datavectors. */
void rm_ivectrvectr_shellonly(ivectrvectr *I_ptr) {
  ivectrvectr I;

  if (!I_ptr || !(I = *I_ptr))
    return;

  nml_free(--I);

  *I_ptr = NULL;

  return;
}

/* Attach a data vector to the holding structure.  Do not create memory for d
   the data vector itself.  */
int add_elem_ivectrvectr(ivectr data_vectr, int pos, ivectrvectr hold_vectr) {
  if (pos > vdim(hold_vectr)) {
    die("Error in add_elem_ivectrvectr\n"
        "Requested Position in the vector of vectors\n"
        "exceeds holding capacity.\n"
        "pos=%d dim=%d\n",
        pos, vdim(hold_vectr));
    return -1;
  }

  if (hold_vectr[pos] != NULL) {
    rm_ivectr(&hold_vectr[pos]);
  }

  hold_vectr[pos] = data_vectr;

  return 0;
}

#ifdef IVECTVECTR_TEST
int main(int argc, char *argv[]) {

  int i, ii, n = 10;
  ivectrvectr I = Ivectrvectr(n);
  ivectr sample = NULL;

  for (i = 0; i < n; i++) {
    sample = Ivectr(i + 2);
    for (ii = 0; ii < i + 2; ii++) {
      sample[ii] = ii + 2;
    }
    add_elem_ivectrvectr(sample, i, I);
  }
  printf("\n");

  for (i = 0; i < n; i++) {
    for (ii = 0; ii < nml_iv_extent(I[i]); ii++) {
      printf(" %d", I[i][ii]);
    }
    printf("\n");
  }
  printf("\n");

  memory_report(stdout);
  rm_ivectrvectr(&I);
  memory_report(stdout);

  return 0;
}
#endif
