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
$Header: /repo/nemo3d/src/math/r3tensor.c,v 1.3 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

#include "r3tensor.h"

/*
 * Allocate memory for a 3-dimensional tensor matrix of real element types.
 */
r3tensor R3tensor(int dimRowLow, int dimRowHigh, int dimColLow, int dimColHigh,
                  int dimDepthLow, int dimDepthHigh) {
  r3tensor r3t;
  int i, j, k;
  int dimRow = dimRowHigh - dimRowLow + 1;
  int dimCol = dimColHigh - dimColLow + 1;
  int dimDepth = dimDepthHigh - dimDepthLow + 1;

  r3t = (real ***)nml_malloc((size_t)((dimRow + NR_END + 6) * sizeof(real **)));
  if (!r3t)
    fprintf(LOGFILE,
            "In R3tensor():: Unable to allocate memory for dimRow+NR_END+6.\n");

  r3t += 6;
  d_1l(r3t) = dimRowLow;
  d_1h(r3t) = dimRowHigh;
  d_2l(r3t) = dimColLow;
  d_2h(r3t) = dimColHigh;
  d_3l(r3t) = dimDepthLow;
  d_3h(r3t) = dimDepthHigh;
  r3t += NR_END;
  r3t -= dimRowLow;

  r3t[dimRowLow] = (real **)nml_malloc(
      (size_t)((dimRow * dimCol + NR_END) * sizeof(real *)));
  if (!r3t[dimRowLow])
    fprintf(LOGFILE, "In R3tensor():: Unable to allocate memory for "
                     "dimRow*dimCol+NR_END.\n");
  r3t[dimRowLow] += NR_END;
  r3t[dimRowLow] -= dimColLow;

  r3t[dimRowLow][dimColLow] = (real *)nml_malloc(
      (size_t)((dimRow * dimCol * dimDepth + NR_END) * sizeof(real)));
  if (!r3t[dimRowLow][dimColLow])
    fprintf(LOGFILE, "In R3tensor():: Unable to allocate memory for "
                     "dimRow*dimCol*dimDepth+NR_END.\n");
  r3t[dimRowLow][dimColLow] += NR_END;
  r3t[dimRowLow][dimColLow] -= dimDepthLow;

  for (j = dimColLow + 1; j <= dimColHigh; j++) {
    r3t[dimRowLow][j] = r3t[dimRowLow][j - 1] + dimDepth;
  }
  for (i = dimRowLow + 1; i <= dimRowHigh; i++) {
    r3t[i] = r3t[i - 1] + dimCol;
    r3t[i][dimColLow] = r3t[i - 1][dimColLow] + dimCol * dimDepth;
    for (j = dimColLow + 1; j <= dimColHigh; j++) {
      r3t[i][j] = r3t[i][j - 1] + dimDepth;
    }
  }

  for (i = dimRowLow + NR_END; i <= dimRowHigh + NR_END; i++) {
    for (j = dimColLow + NR_END; j <= dimColHigh + NR_END; j++) {
      for (k = dimDepthLow + NR_END; k <= dimDepthHigh + NR_END; k++) {
        r3t[i][j][k] = 0.0;
      }
    }
  }

  return r3t;
}

/*
 * Deallocate the memory for a r3tensor structure.
 * Note that the dimension size of the row's lower bound is needed as an
 * argument to the destructor.
 */
void rm_r3tensor(r3tensor *ptrR3t, int dimRowLow) {
  r3tensor r3t = *ptrR3t;
  int dimRowHigh, dimColLow, dimColHigh, dimDepthLow, dimDepthHigh;

  if (r3t != NULL) {

    r3t += dimRowLow;
    if (d_1l(r3t) != dimRowLow) {
      n3d_warning("In rm_r3tensor():: Incorrect match of 1st dimension's lower "
                  "bound %d, with the passed in value of %d.\n",
                  dimRowLow, d_1l(r3t));
    }

    dimRowHigh = d_1h(r3t);
    dimColLow = d_2l(r3t);
    dimColHigh = d_2h(r3t);
    dimDepthLow = d_3l(r3t);
    dimDepthHigh = d_3h(r3t);
    r3t -= dimRowLow;

    nml_free(r3t[dimRowLow][dimColLow] + dimDepthLow - NR_END);
    nml_free(r3t[dimRowLow] + dimColLow - NR_END);
    nml_free(r3t + dimRowLow - NR_END - 6);

    *ptrR3t = NULL;
  }

  return;
}

/*
 * Checks the validity of the r3tensor.
 * If r3t does not exists, it creates a new one. But if r3t does exist, a check
 * on the dimension is performed. If the dimensions are valid, the r3tensor is
 * returned. Otherwise, it removes the invalid r3tensor and creates a new one
 * with the correct dimensions.
 */
r3tensor check_r3tensor(r3tensor r3t, int d1l, int d1h, int d2l, int d2h,
                        int d3l, int d3h) {

  if (!r3t) {
    return R3tensor(d1l, d1h, d2l, d2h, d3l, d3h);
  }

  else {

    r3t += d1l;

    if (d_1l(r3t) != d1l) {
      n3d_warning("In check_r3tensor():: Incorrect match of 1st dimension's "
                  "lower bound %d, with the passed in value of %d.\n",
                  d1l, d_1l(r3t));
    }

    /* Check dimensions */
    if ((d_1h(r3t) == d1h) && (d_2l(r3t) == d2l) && (d_2h(r3t) == d2h) &&
        (d_3l(r3t) == d3l) && (d_3h(r3t) == d3h)) {
      r3t -= d1l;
    }
    /* Else, incorrect dimensions. So remove the old r3tensor and create a new
       one. */
    else {
      r3t -= d1l;
      rm_r3tensor(&r3t, d1l);
      r3t = R3tensor(d1l, d1h, d2l, d2h, d3l, d3h);
    }

    return r3t;
  }
}
