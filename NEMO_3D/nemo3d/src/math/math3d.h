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
$Header: /repo/nemo3d/src/math/math3d.h,v 1.2 2003/10/08 16:18:03 hook Exp $ 
*****************************************************************************/

#ifndef MATH3D_H
#define MATH3D_H 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <rvector.h>
#include <cvector.h>
#include <rtridiag.h>
#include <cmatrix.h>

#include "io_utils.h"
#include "realtype.h"
#include "cvector_lib.h"
#include "rvector_lib.h"
#include "ctridiag.h"

#include "lanczos.h"

void eigvec_srt( cvectr vec, cmatrix cmat, int n );
real invit ( rtridiag h, rvectr x, rvectr y, real lambda );

/* Called by eigval_s_tridi_QD() below, which is nolonger needed. */
int td_ql_im_QD (rvectr d, rvectr e, roperator z, int n);

/*
 * Called from eigsys3d_par.[ch]::srch_converged_eigval_par(),
 * which is now made obsolete by LAPACK.
 */
rvectr eigval_s_tridi_QD (rvectr d, rtridiag t, int n);

#endif

