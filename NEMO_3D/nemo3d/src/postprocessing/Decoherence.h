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
$Header: /repo/nemo3d/src/postprocessing/Decoherence.h,v 1.1 2004/08/26 21:28:07
swlee Exp $
*****************************************************************************/

#ifndef __decoherence_h_
#define __decoherence_h_
#include "RandRandom.h"
#include "cvector_lib.h"
#include "io_utils.h"
#include "qd_struct.h"
#include "util_stringops.h"
#include <fstream>
#include <iomanip>
#include <math.h>
#include <vector>

class Decoherence {

public:
  Decoherence(qd_struct din, int e, int h, int mpiID)
      : d(din), e_max(e), h_max(h), my_id(mpiID) {
    g_min = d->cell_s[my_id] * 8;
    g_max = d->cell_ln[my_id] * 8 + g_min;
    g_num = g_max - g_min;
    o_max = d->NBasisStates / 2;
    numAtoms = d->n_atom_tot;
    decoherence_time = new double[e_max + h_max];
    dephasing_time = new double[e_max + h_max];
    alloy_dephasing_time = new double[e_max + h_max];
    interface_dephasing_time = new double[e_max + h_max];
    MagneticField = 0.0;
    RandomOrientation = true;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);
    if (mpi_numprocs == 1)
      g_max = numAtoms;
  }

  ~Decoherence() {
    delete[] decoherence_time;
    delete[] dephasing_time;
  }

  double MagneticField;
  bool RandomOrientation;

  void compute_decoherence_time();
  void compute_spin_correlator(const char &, int n);
  void get_data(complex *, complex *, int *, int *, int *, double *);
  void get_geometry(double, double, double, double, double, double, double);
  bool CationInsideDot(int);
  bool CationNearInterface(int);

  inline complex &cf_e(int e, int g, int s, int o) {
    return _cf_e[((e * numAtoms + g) * 2 + s) * o_max + o];
  }
  inline complex &cf_h(int h, int g, int s, int o) {
    return _cf_h[((h * numAtoms + g) * 2 + s) * o_max + o];
  }
  inline double &lattice(int iat, int dim) { return _lattice[iat * 3 + dim]; }
  inline int &neighbor(int ia, int in) { return _nn[ia * 4 + in]; }

private:
  qd_struct d;
  int e_max, h_max, g_min, g_max, g_num, o_max, numAtoms, my_id;
  complex *_cf_e;
  complex *_cf_h;
  int *_atom;
  int *_shape;
  int *_nn;
  double Xorigin, Yorigin, Zorigin, Xdim, Ydim, Zdim;
  double Xcenter, Ycenter, Zcenter, Xradius2, Yradius2, Zradius2;
  double Composition;
  double *_lattice;
  double *decoherence_time;
  double *dephasing_time;
  double *alloy_dephasing_time;
  double *interface_dephasing_time;
  int mpi_numprocs;
};

#endif
