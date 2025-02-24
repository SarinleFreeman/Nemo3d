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
$Header: /repo/nemo3d/src/base/eigsys3d.h,v 1.2 2003/10/08 16:16:26 hook Exp $ 
*****************************************************************************/

#ifndef EIGSYS3D_H
#define EIGSYS3D_H 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <rvector.h>
#include <i2tensor.h>
#include <i3tensor.h>
#include <r3tensor.h>
#include <c3tensor.h>

#include "io_utils.h"
#include "Boolean.h"
#include "realtype.h"
#include "constants_nemo.h"
#include "cmatrix_nonsym.h"
#include "opt3d.h"
#include "mb_ham_spds_micro.h" 

#include "math3d.h"
#include "qd_struct.h"
#include "h_cvectr_mult.h"
#include "rayleigh_ritz.h"
#include "i_3dout_files.h"
#include "matmul.h"
#include "run3d_mpi.h"



int sym_lanc_it_c( qd_struct d, real shift, real *alpha, real *beta, real beta_m1, cvectr r, cvectr q, cvectr q_m1, cvectr q_m2 );

int h_eigval ( qd_struct d, real s );
int h_ceigval ( qd_struct d, real s );
int readstates_driver_c( qd_struct d );

int ns_lanc_driver( qd_struct d );
int sym_lanc_it( qd_struct d, real shift, real *alpha, real *beta, real beta_m1, rvectr r, rvectr q, rvectr q_m1, rvectr q_m2 );
int nonsym_lanc_it( qd_struct d, real shift, complex *alpha, complex *beta, complex *, cvectr r, cvectr s, cvectr q, cvectr q_m1, cvectr p, cvectr p_m1 );


#endif


