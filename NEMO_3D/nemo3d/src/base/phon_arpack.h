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
$Header: /repo/nemo3d/src/base/phon_arpack.h,v 1.1 2004/02/27 04:09:05 lol Exp $
*****************************************************************************/

#ifndef PHON_ARPACK_H
#define PHON_ARPACK_H

#ifdef NO_PARPACK

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "realtype.h"
#include "MPI_Timing.h"


#include "qd_struct.h"
#include "dmk_mult.h"


#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#ifdef FORTRAN_UNDERSCORE
#define znaupd   znaupd_
#define zneupd   zneupd_
#define zaxpy    zaxpy_
#define dznorm2  dznorm2_
#endif

#ifdef FORTRAN_UNDERSCORE2
#define znaupd   znaupd__
#define zneupd   zneupd__
#define zaxpy    zaxpy__
#define dznorm2  dznorm2__
#endif

#ifdef FORTRAN_ALLCAPS
#define znaupd   ZNAUPD
#define zneupd   ZNEUPD
#define zaxpy    ZAXPY
#define dznorm2  DZNORM2
#endif


void znaupd(int*, char*, int*, char*, int*, double*, complex*, int*, complex*,
	      int*, int*, int*, complex*, complex*, int*, double*, int*);
void zneupd(int*, char*, int*, complex*, complex*, int*, complex*,
	      complex*, char*, int*, char*, int*, double*, complex*, int*,
	      complex*, int*, int*, int*, complex*, complex*, int*, double*,
	      int*);
void zaxpy(int*, complex*, complex*, int*, complex*, int*);
double dznorm2(int*, int*, complex*, int*);

#ifdef __cplusplus
}
#endif /* __cplusplus */


int phon_arpack( qd_struct d,
                 int iq_ph/*number of the point in the reciprocal space
                            (needed for writing eigenvalues to
                            d->phon.E_ph[iq_ph][d->d->opt.ExecParam.PhonDisp.N_branches+3])*/
                 );//diagonalizer ;)

#endif /* NO_PARPACK */

#endif /* PHON_PARPACK_H */
