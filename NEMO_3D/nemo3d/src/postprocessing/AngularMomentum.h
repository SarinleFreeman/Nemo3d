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
$Header: /repo/nemo3d/src/postprocessing/AngularMomentum.h,v 1.1 2004/08/26
21:28:07 swlee Exp $
*****************************************************************************/

#ifndef __angularmomentum_h_
#define __angularmomentum_h_
#include "cvector_lib.h"
#include "h_cvectr_mult.h"
#include "io_utils.h"
#include "qd_struct.h"
#include <fstream>

class AngularMomentum {

public:
  AngularMomentum(qd_struct din, int NumElectron, int NumHole, int mpiID)
      : d(din), e_max(NumElectron), h_max(NumHole), my_id(mpiID) {
    g_min = d->cell_s[my_id] * 8;
    g_max = d->cell_ln[my_id] * 8 + g_min;
    o_max = d->NBasisStates / 2;
    numAtoms = d->n_atom_tot;
    MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);
    if (mpi_numprocs == 1)
      g_max = numAtoms;
  }

  ~AngularMomentum() {
    rm_cvectr(&MatrixLx);
    rm_cvectr(&MatrixLy);
    rm_cvectr(&MatrixLz);
  }

  void compute_angular_momentum_spin(const char &);
  complex compute_global_angular_momentum(int, const char &, const char &);
  complex compute_global_angular_momentum_square(int, const char &,
                                                 const char &);
  complex compute_local_angular_momentum(int, int, const char &, const char &);
  complex compute_spin(int, int, const char &, const char &);
  void prepare_MatrixL(void);
  void set_origin(double, double, double);
  double Hamiltonian(int, int, int, int);
  void get_data(complex *, complex *, int *, int *, int *, double *);
  inline complex &cf_e(int e, int g, int s, int o) {
    return _cf_e[((e * numAtoms + g) * 2 + s) * o_max + o];
  }
  inline complex &cf_h(int h, int g, int s, int o) {
    return _cf_h[((h * numAtoms + g) * 2 + s) * o_max + o];
  }
  inline double &lattice(int iat, int dim) { return _lattice[iat * 3 + dim]; }
  inline int &nn(int iat, int nn) { return _nn[iat * 4 + nn]; }

private:
  qd_struct d;
  int e_max, h_max;
  int my_id, g_min, g_max, o_max, numAtoms;
  double X0, Y0, Z0;
  complex *_cf_e;
  complex *_cf_h;
  int *_atom;
  int *_shape;
  int *_nn;
  double *_lattice;
  complex *MatrixLx;
  complex *MatrixLy;
  complex *MatrixLz;
  int mpi_numprocs;
};

#endif
