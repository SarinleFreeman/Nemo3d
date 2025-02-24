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
$Header: /repo/nemo3d/src/base/eigsys3d_par.h,v 1.15 2008/06/03 07:48:10 hoonryu Exp $ 
*****************************************************************************/

#ifndef EIGSYS3D_PAR_H
#define EIGSYS3D_PAR_H

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <complex.h>
#include <rvector.h>

#include "io_utils.h"
#include "nml_global.h"
#include "lanczos.h"
#include "lanczos_blk.h"
#include "qd_struct.h"
#include "run3d_mpi.h"
#include "realtype.h"
#include "cmatrix_nonsym.h"
#include "h_cvectr_mult.h"
#include "h_cmat_mult.h"
#include "eigsys3d.h"
#include "ham_local.h"
#include "rayleigh_ritz.h"

#ifdef MPI_TIMING
#include "MPI_Timing.h"
#endif

#include "arpack_driver.h"
#include "parpack_driver.h"

#if (defined USE_TRACEMIN_PAR || defined USE_TRACEMIN_SER)
#include "ctracemin_driver.h"
#include "qtracemin_driver.h"
#endif

#if (!defined MPI3d && defined FAKE_MPI)
#include "mpi_fake.h"
#endif /* (defined MPI3d && !defined FAKE_MPI) */

int lanczos_initialize(qd_struct d);
int lanczos_driver(qd_struct d);
int h_eigval_spds ( qd_struct d, real emin, real emax );
void compute_electronic_structure( qd_struct d );
void verify_ham( qd_struct d);
void verify_ham_col(qd_struct d, int col0, int colf);
void generate_wave_vectors( qd_struct d, rmatrix *kr, int* Nk );

#endif /* EIGSYS3D_PAR_H */
