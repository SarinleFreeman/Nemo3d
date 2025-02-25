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
$Header: /repo/nemo3d/src/base/strain3d.h,v 1.8 2005/01/14 20:24:39 marek Exp $
*****************************************************************************/

#ifndef STRAIN3D_H
#define STRAIN3D_H 1

#include "Boolean.h"
#include "cmacopt.h"
#include "constants_nemo.h"
#include "i_3dout_files.h"
#include "io_utils.h"
#include "macopt.h"
#include "macopt_par.h"
#include "nml_global.h"
#include "qd_struct.h"
#include "realtype.h"
#include "run3d_mpi.h"
#include "util_stringops.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "MPI_Timing.h"

void min_estrain(qd_struct d);
real calc_estrain(qd_struct d, rvectr p, rvectr rcvL, rvectr rcvR);

void initStrainComm(qd_struct d);

void calc_estrain_macopt(rvectr p, int sz, rvectr xi, void *data, int itercount,
                         real *onevalue);

void calc_estrain_cmacopt(cvectr p, cvectr xi, void *data);

int entryInRcvWithThisAtom(int cell, int atom, ivectr rcvEntries,
                           int noAtomsPerCell);

#endif
