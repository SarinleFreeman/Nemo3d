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
$Header: /repo/nemo3d/src/math/c3tensor.c,v 1.3 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

#include "c3tensor.h"

/*
 * Allocate memory for a 3-dimensional tensor matrix of complex element types.
 */
c3tensor C3tensor(int dimRowLow, int dimRowHigh, int dimColLow, int dimColHigh,
                  int dimDepthLow, int dimDepthHigh) {
  c3tensor c3t;
  int i, j, k;
  int dimRow = dimRowHigh - dimRowLow + 1;
  int dimCol = dimColHigh - dimColLow + 1;
  int dimDepth = dimDepthHigh - dimDepthLow + 1;

  c3t = (complex ***)nml_malloc(
      (size_t)((dimRow + NR_END + 6) * sizeof(complex **)));
  if (!c3t)
    fprintf(LOGFILE,
            "In C3tensor():: Unable to allocate memory for dimRow+NR_END+6.\n");

  c3t += 6;
  d_1l(c3t) = dimRowLow;
  d_1h(c3t) = dimRowHigh;
  d_2l(c3t) = dimColLow;
  d_2h(c3t) = dimColHigh;
  d_3l(c3t) = dimDepthLow;
  d_3h(c3t) = dimDepthHigh;
  c3t += NR_END;
  c3t -= dimRowLow;

  c3t[dimRowLow] = (complex **)nml_malloc(
      (size_t)((dimRow * dimCol + NR_END) * sizeof(complex *)));
  if (!c3t[dimRowLow])
    fprintf(LOGFILE, "In C3tensor():: Unable to allocate memory for "
                     "dimRow*dimCol+NR_END.\n");
  c3t[dimRowLow] += NR_END;
  c3t[dimRowLow] -= dimColLow;

  c3t[dimRowLow][dimColLow] = (complex *)nml_malloc(
      (size_t)((dimRow * dimCol * dimDepth + NR_END) * sizeof(complex)));
  if (!c3t[dimRowLow][dimColLow])
    fprintf(LOGFILE, "In C3tensor():: Unable to allocate memory for "
                     "dimRow*dimCol*dimDepth+NR_END.\n");
  c3t[dimRowLow][dimColLow] += NR_END;
  c3t[dimRowLow][dimColLow] -= dimDepthLow;

  for (j = dimColLow + 1; j <= dimColHigh; j++) {
    c3t[dimRowLow][j] = c3t[dimRowLow][j - 1] + dimDepth;
  }
  for (i = dimRowLow + 1; i <= dimRowHigh; i++) {
    c3t[i] = c3t[i - 1] + dimCol;
    c3t[i][dimColLow] = c3t[i - 1][dimColLow] + dimCol * dimDepth;
    for (j = dimColLow + 1; j <= dimColHigh; j++) {
      c3t[i][j] = c3t[i][j - 1] + dimDepth;
    }
  }

  for (i = dimRowLow + NR_END; i <= dimRowHigh + NR_END; i++) {
    for (j = dimColLow + NR_END; j <= dimColHigh + NR_END; j++) {
      for (k = dimDepthLow + NR_END; k <= dimDepthHigh + NR_END; k++) {
        c3t[i][j][k] = czero;
      }
    }
  }

  return c3t;
}

/*
 * Deallocate the memory for a c3tensor structure.
 * Note that the dimension size of the row's lower bound is needed as an
 * argument to the destructor.
 */
void rm_c3tensor(c3tensor *ptrC3t, int dimRowLow) {
  c3tensor c3t = *ptrC3t;
  int dimRowHigh, dimColLow, dimColHigh, dimDepthLow, dimDepthHigh;

  if (c3t != NULL) {

    c3t += dimRowLow;
    if (d_1l(c3t) != dimRowLow) {
      n3d_warning("In rm_c3tensor():: Incorrect match of 1st dimension's lower "
                  "bound %d, with the passed in value of %d.\n",
                  dimRowLow, d_1l(c3t));
    }

    dimRowHigh = d_1h(c3t);
    dimColLow = d_2l(c3t);
    dimColHigh = d_2h(c3t);
    dimDepthLow = d_3l(c3t);
    dimDepthHigh = d_3h(c3t);
    c3t -= dimRowLow;

    nml_free(c3t[dimRowLow][dimColLow] + dimDepthLow - NR_END);
    nml_free(c3t[dimRowLow] + dimColLow - NR_END);
    nml_free(c3t + dimRowLow - NR_END - 6);

    *ptrC3t = NULL;
  }

  return;
}

/*
 * Checks the validity of the c3tensor.
 * If c3t does not exists, it creates a new one. But if c3t does exist, a check
 * on the dimension is performed. If the dimensions are valid, the c3tensor is
 * returned. Otherwise, it removes the invalid c3tensor and creates a new one
 * with the correct dimensions.
 */
c3tensor check_c3tensor(c3tensor c3t, int d1l, int d1h, int d2l, int d2h,
                        int d3l, int d3h) {

  if (!c3t) {
    return C3tensor(d1l, d1h, d2l, d2h, d3l, d3h);
  }

  else {

    c3t += d1l;

    if (d_1l(c3t) != d1l) {
      n3d_warning("In check_c3tensor():: Incorrect match of 1st dimension's "
                  "lower bound %d, with the passed in value of %d.\n",
                  d1l, d_1l(c3t));
    }

    /* Check dimensions */
    if ((d_1h(c3t) == d1h) && (d_2l(c3t) == d2l) && (d_2h(c3t) == d2h) &&
        (d_3l(c3t) == d3l) && (d_3h(c3t) == d3h)) {
      c3t -= d1l;
    }
    /* Else, incorrect dimensions. So remove the old c3tensor and create a new
       one. */
    else {
      c3t -= d1l;
      rm_c3tensor(&c3t, d1l);
      c3t = C3tensor(d1l, d1h, d2l, d2h, d3l, d3h);
    }

    return c3t;
  }
}
