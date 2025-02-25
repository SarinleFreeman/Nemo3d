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
$Header: /repo/nemo3d/src/base/matmul.h,v 1.4 2006/12/18 21:40:04 gekco Exp $
*****************************************************************************/

#ifndef MATMUL_H
#define MATMUL_H

#include "cmatrix.h"
#include "complex.h"
#include "cvector.h"
#include "realtype.h"

int cmatmul_spds(cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim);
int cmatmul_spds_20(cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim);
int cmatmul_spds_10(cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim);
int cmatmul_spds_hc(cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim);
int cmatmul_spds_hc_20(cvectr y, cmatrix h, cvectr x, int isy, int isx,
                       int hdim);
int cmatmul_spds_hc_10(cvectr y, cmatrix h, cvectr x, int isy, int isx,
                       int hdim);

int cmatmul_spds_sc(cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim,
                    complex scalar);

#define CMATMUL_SPDS_1_MACRO

#ifdef CMATMUL_SPDS_1_MACRO
#define cmatmul_spds_1(y, h, x, isy, isx, hdim)                                \
  {                                                                            \
    y[isy].r += h[0][0].r * x[isx].r - h[0][0].i * x[isx].i;                   \
    y[isy].i += h[0][0].r * x[isx].i + h[0][0].i * x[isx].r;                   \
  }
#define cmatmul_spds_hc_1(y, h, x, isy, isx, hdim)                             \
  {                                                                            \
    y[isy].r += h[0][0].r * x[isx].r + h[0][0].i * x[isx].i;                   \
    y[isy].i += h[0][0].r * x[isx].i - h[0][0].i * x[isx].r;                   \
  }

// Hoon and sunnylee - this is for "call by value" function,
// changed h[0][0] to val (it looks trivial, however, it helps to speed up a
// little
#define cmatmul_sparse_spds_1(y, val, x, isy, isx)                             \
  {                                                                            \
    y[isy].r += val * x[isx].r;                                                \
    y[isy].i += val * x[isx].i;                                                \
  }
#define cmatmul_sparse_spds_hc_1(y, val, x, isy, isx)                          \
  {                                                                            \
    y[isy].r += val * x[isx].r;                                                \
    y[isy].i += val * x[isx].i;                                                \
  }
// Hoon and sunnylee

#else /* CMATMUL_SPDS_1_MACRO */

int cmatmul_spds_1(cvectr y, cmatrix h, cvectr x, int isy, int isx, int hdim);
int cmatmul_spds_hc_1(cvectr y, cmatrix h, cvectr x, int isy, int isx,
                      int hdim);

int cmatmul_sparse_spds_1(cvectr y, double h, cvectr x, int isy, int isx);
int cmatmul_sparse_spds_hc_1(cvectr y, double h, cvectr x, int isy, int isx);

#endif /* CMATMUL_SPDS_1_MACRO */

#endif /* MATMUL_H */
