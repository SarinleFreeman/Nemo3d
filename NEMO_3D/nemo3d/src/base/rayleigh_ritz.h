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
$Header: /repo/nemo3d/src/base/rayleigh_ritz.h,v 1.5 2004/10/27 19:43:10 gekco Exp $ 
*****************************************************************************/

#ifndef RAYLEIGH_RITZ_H
#define RAYLEIGH_RITZ_H 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <rvector.h>
#include <i2tensor.h>
#include <i3tensor.h>

#include "nml_global.h"
#include "realtype.h"
#include "constants_nemo.h"
#include "cmatrix_nonsym.h"
#include "qd_struct.h"
#include "run3d_mpi.h"
#include "eigsys3d_par.h"
#include "strain3d.h"
#include "mat3d.h"
#include "ivectorvector.h"
#include "matsite_init.h"
#define RESTART_GEKCO
#ifdef RESTART_GEKCO	
#include "cmacopt.h"
#endif /* RESTART_GEKCO */

#include "h_cvectr_mult.h"

void linSolv_CGNE(cvectr x, qd_struct d, real shift, cvectr r, real tol, 
                  real* xx, real* xHx);
void rayleighQuotientIter(qd_struct d);


int min_rayleigh(qd_struct d);
real calc_rayleigh( rvectr p, void *data );
real calc_rayleigh_sq( rvectr p, void *data );
void calc_rayleigh_d( rvectr p, rvectr xi, void *data );
complex calc_rayleigh_compl( cvectr x, int sz, qd_struct d );



#ifdef MPI3d
void calc_rayleigh_sq_macopt( rvectr p, int sz, rvectr xi, int, void *data, int itercount, real *onevalue );
void calc_rayleigh_sq_cmacopt( cvectr p, int sz, cvectr xi, int, void *data, int itercount, real *onevalue );
#else /* MPI3d */
void calc_rayleigh_sq_macopt( rvectr p, rvectr xi, void *data, int itercount, real *onevalue );
void calc_rayleigh_sq_cmacopt( cvectr p, cvectr xi, void *data, int itercount, real *onevalue );
#endif /* MPI3d */


#endif /* RAYLEIGH_RITZ_H */


