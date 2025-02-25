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
$Header: /repo/nemo3d/src/base/ham_local.h,v 1.4 2004/10/27 19:43:10 gekco Exp $
*****************************************************************************/

#ifndef HAM_LOCAL_H
#define HAM_LOCAL_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "Boolean.h"
#include "io_utils.h"
#include "qd_struct.h"

#include "nml_global.h"
#include "realtype.h"
#include <rvector.h>

#include "MPI_Timing.h"
#include "cmatrix_nonsym.h"
#include "constants_nemo.h"
#include "h_cvectr_mult.h"
#include "mb_ham_spds_micro.h"
#include "run3d_mpi.h"

void local_Hoff_zincBlende(qd_struct d, cmatrix h, int l, int m,
                           int map_mm_to_n[8][8], int map_mb_to_m[8][4],
                           real a_lattice[3], int block, int nbr5, int nbr6,
                           int nbr7, real phase5, real phase6, real phase7);
void local_bandstruct_zincBlende_cubic(qd_struct d);

void local_bandstruct_cubic(qd_struct d);
#endif
