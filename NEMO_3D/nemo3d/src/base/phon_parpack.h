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
$Header: /repo/nemo3d/src/base/phon_parpack.h,v 1.4 2004/02/27 04:09:05 lol Exp
$
*****************************************************************************/

#ifndef PHON_PARPACK_H
#define PHON_PARPACK_H

#ifndef NO_PARPACK

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "MPI_Timing.h"
#include "realtype.h"

#include "dmk_mult.h"
#include "qd_struct.h"

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef FORTRAN_UNDERSCORE
#define pznaupd pznaupd_
#define pzneupd pzneupd_
#define zaxpy zaxpy_
#define pdznorm2 pdznorm2_
#endif

#ifdef FORTRAN_UNDERSCORE2
#define pznaupd pznaupd__
#define pzneupd pzneupd__
#define zaxpy zaxpy__
#define pdznorm2 pdznorm2__
#endif

#ifdef FORTRAN_ALLCAPS
#define pznaupd PZNAUPD
#define pzneupd PZNEUPD
#define zaxpy ZAXPY
#define pdznorm2 PDZNORM2
#endif

void pznaupd(int *, int *, char *, int *, char *, int *, double *, complex *,
             int *, complex *, int *, int *, int *, complex *, complex *, int *,
             double *, int *);
void pzneupd(int *, int *, char *, int *, complex *, complex *, int *,
             complex *, complex *, char *, int *, char *, int *, double *,
             complex *, int *, complex *, int *, int *, int *, complex *,
             complex *, int *, double *, int *);
void zaxpy(int *, complex *, complex *, int *, complex *, int *);
double pdznorm2(int *, int *, complex *, int *);

#ifdef __cplusplus
}
#endif /* __cplusplus */

int phon_parpack(
    qd_struct d,
    int iq_ph /*number of the point in the reciprocal space
                (needed for writing eigenvalues to
                d->phon.E_ph[iq_ph][d->d->opt.ExecParam.PhonDisp.N_branches+3])*/
);            // diagonalizer ;)

#endif /* NO_PARPACK */

#endif /* PHON_PARPACK_H */
