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
$Header: /repo/nemo3d/src/base/opt3d.h,v 1.3 2005/11/12 19:14:42 gekco Exp $
*****************************************************************************/

#ifndef OPT3D_H
#define OPT3D_H 1

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "io_utils.h"

#include "Boolean.h"
#include "cmatrix_nonsym.h"
#include "constants_nemo.h"
#include "i2tensor.h"
#include "i3tensor.h"
#include "qd_struct.h"
#include "realtype.h"
#include "rvector.h"

real calc_rate(qd_struct d, rvectr init, rvectr finl, real px, real py,
               real pz);
int calc_opt_rates(qd_struct d, real px, real py, real pz);
real calc_opt_rates_int(qd_struct d, real px, real py, real pz);
int calc_opt_mat(qd_struct d, real px, real py, real pz);
int calc_opt_mat_old(qd_struct d, real px, real py, real pz);
int calc_opt_rates_int_ramp(qd_struct d, real *target, real start, real end,
                            int n);

int calc_opt_mat_cmplx(qd_struct d, real ephi, real etheta);

int calc_rate_cmplx(rvectr optmat_vec, qd_struct d, cvectr init, cvectr finl,
                    real px, real py, real pz);

rvectr make_egrid_lineshape(rvectr trans, int n_trans, real Gamma, int *ne);

#endif
