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
$Header: /repo/nemo3d/src/postprocessing/TimeReversal.h,v 1.1 2004/08/26
21:28:07 swlee Exp $
*****************************************************************************/

#ifndef __timereversal_h_
#define __timereversal_h_
#include "cvector_lib.h"
#include "i_3dout_files.h"
#include "io_utils.h"
#include "qd_struct.h"
#include "util_stringops.h"
#include <fstream>
#include <vector>

class TimeReversal {

public:
  TimeReversal(qd_struct din, int mpiID) : d(din), my_id(mpiID) {
    g_min = d->cell_s[my_id] * d->geo.AtomsPerCellMax();
    g_max = d->cell_ln[my_id] * d->geo.AtomsPerCellMax() + g_min;
    o_max = d->NBasisStates / 2;
    numAtoms = d->n_atom_tot;
    _cf_tr = Cvectr(numAtoms * o_max * 2);
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);
    if (mpi_numprocs == 1)
      g_max = numAtoms;
  }

  ~TimeReversal() {}

  void get_data(complex *);
  void apply_time_reversal_operator(char *);

  inline complex &cf(int g, int s, int o) {
    return _cf[(g * 2 + s) * o_max + o];
  }

  inline complex &cf_tr(int g, int s, int o) {
    return _cf_tr[(g * 2 + s) * o_max + o];
  }

private:
  qd_struct d;
  int g_min, g_max, o_max, numAtoms, my_id;
  complex *_cf;
  complex *_cf_tr;
  int mpi_numprocs;
};

#endif
