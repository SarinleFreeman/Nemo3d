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
$Header: /repo/nemo3d/src/math/cmatrix_nonsym.h,v 1.3 2005/04/22 20:16:03
dseaman Exp $
*****************************************************************************/

#ifndef CMATRIX_NONSYM_H
#define CMATRIX_NONSYM_H

#include <math.h>
#include <stddef.h>
#include <stdlib.h>

#include "cmatrix_lib.h"
#include "util_memory.h"
#include "util_messages.h"
#include <cmatrix.h>
#include <complex.h>
#include <cvector.h>
#include <matrix.h> /* used for cdim() and vdim() definitions */
#include <rmatrix.h>
#include <rvector.h>
#include <vector.h>

#ifndef MACHINE_EPSILON
#define MACHINE_EPSILON
const double MachineEpsilon = 1.0E-15;
#endif /* MACHINE_EPSILON */

#define Swap(g, h)                                                             \
  {                                                                            \
    temp = (g);                                                                \
    (g) = (h);                                                                 \
    (h) = temp;                                                                \
  }

#ifdef __cplusplus
extern "C" {
#endif

#ifdef FORTRAN_UNDERSCORE
#define zheev zheev_
#define zgehrd zgehrd_
#define zhseqr zhseqr_
#endif

#ifdef FORTRAN_UNDERSCORE2
#define zheev zheev__
#define zgehrd zgehrd__
#define zhseqr zhseqr__
#endif

#ifdef FORTRAN_ALLCAPS
#define zheev ZHEEV
#define zgehrd ZGEHRD
#define zhseqr ZHSEQR
#endif

void zheev(char *, char *, int *, nml_dcscalar *, int *, double *, cvectr,
           int *, double *, int *);
void zgehrd(int *, int *, int *, nml_dcscalar *, int *, cvectr, cvectr, int *,
            int *);
void zhseqr(char *, char *, int *, int *, int *, nml_dcscalar *, int *, cvectr,
            nml_dcscalar *, int *, cvectr, int *, int *);

#ifdef __cplusplus
}
#endif

void ceigval_ns_full(cvectr cvEigenValue, cmatrix h);
void ceigvec_ns_full(cvectr, cmatrix, cmatrix);

#endif /* CMATRIX_NONSYM_H */
