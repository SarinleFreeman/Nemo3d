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
$Header: /repo/nemo3d/src/math/cmatrix_lib.c,v 1.3 2004/12/02 12:09:56 marek Exp
$
*****************************************************************************/

#include "cmatrix_lib.h"

/*
 * Allocate memory for a square complex operator matrix of dimension n.
 * Returns the pointer from Cmatrix(n,n)
 */
coperator Coperator(int n) { return Cmatrix(n, n); }

/*
 * Solves the system "y = Ax" for cvector y. It multiplies cmatrix A with
 * cvector y. Returns the resultant cvector y.
 */
cvectr mul_cmatrix(cvectr y, cmatrix a, cvectr x) {
  int i, j;
  complex *vtrA;

  if (cdim(a) != nml_dcv_extent(x)) {
    n3d_warning("In mul_cmatrix():: Columns of A  !=  dim of x.\n");
  }

  /* Remove any current cvector of y */
  if ((y != NULL) && (nml_dcv_extent(y) != rdim(a))) {
    rm_cvectr(&y);
  }

  if (!y) {
    y = Cvectr(rdim(a));
  }

  for (i = 0; i < rdim(a); i++) {
    vtrA = a[i];
    y[i] = czero;

    for (j = 0; j < cdim(a); j++) {
      y[i].r += vtrA[j].r * x[j].r - vtrA[j].i * x[j].i;
      y[i].i += vtrA[j].i * x[j].r + vtrA[j].r * x[j].i;
    }
  }

  return y;
}

/*
 * Copies elements of cmatrix cmSource into cmatrix cmTarget.
 * cmTarget is returned.
 */
cmatrix cp_cmatrix(cmatrix cmTarget, cmatrix cmSource) {
  int i;

  /* If same cmatrix, nothing to do */
  if (cmTarget == cmSource) {
    return cmTarget;
  }

  /* If no cmatrix cmTarget exists, create one with the same dimension as
   * cmSource. */
  if (!cmTarget) {
    cmTarget = Cmatrix(rdim(cmSource), cdim(cmSource));
  } else { /* cmTarget exists */
    /* Make sure cmTarget has the same row and column dimensions as cmSource */
    if (rdim(cmTarget) != rdim(cmSource) || cdim(cmTarget) != cdim(cmSource)) {
      rm_cmatrix(&cmTarget);
      cmTarget = Cmatrix(rdim(cmSource), cdim(cmSource));
    }
  }

  for (i = 0; i < nml_dcv_extent(*cmSource); i++) {
    (*cmTarget)[i] = (*cmSource)[i];
  }

  return cmTarget;
}

/*
 * Performs a transpose on the cmatrix cm.
 */
cmatrix ctranspose_nocopy(cmatrix cm) {
  int i, j;
  complex tempComplex;

  for (i = 0; i < vdim(cm); i++) {
    for (j = i + 1; j < vdim(cm); j++) {
      tempComplex = cm[i][j];
      cm[i][j] = cm[j][i];
      cm[j][i] = tempComplex;
    }
  }

  return cm;
}
