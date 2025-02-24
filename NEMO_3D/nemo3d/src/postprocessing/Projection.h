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
$Header: /repo/nemo3d/src/postprocessing/Projection.h,v 1.1 2004/08/26 21:28:07 swlee Exp $
*****************************************************************************/

#ifndef __projection_h_
#define __projection_h_
#include <vector>
#include <fstream>
#include "cvector_lib.h"
#include "io_utils.h"
#include "qd_struct.h"
#include "util_stringops.h"

class Projection{

public:
  
  Projection(qd_struct din, int e, int h, int mpiID)
  : d(din), e_max(e), h_max(h), my_id(mpiID) 
  {
    g_min = d->cell_s[my_id]*8;
    g_max = d->cell_ln[my_id]*8+g_min;
    o_max = d->NBasisStates/2;
    numAtoms = d->n_atom_tot;   
    _coef       = Cvectr(2*o_max);

   MPI_Comm_size(MPI_COMM_WORLD, &mpi_numprocs);
   if(mpi_numprocs==1) g_max = numAtoms;

  }
  
  ~Projection()
  {
    rm_cvectr(&_coef);
  }

  void get_data(complex*, complex*, int*, int*, int*, double*);
  void compute_projection();
  void compute_weight_on_shape(int);
  void compute_qd_wl_extent(int);
  void compute_extent(int);
  void compute_overlap(int, int); 
  void compute_overlap_without_spin(int, int);
  void print_amplitude(int);
  void compute_electron_hole_alignment(int, int);

  inline complex& cf_e(int e, int g, int s, int o) 
    { return _cf_e[((e*numAtoms+g)*2+s)*o_max+o];}
  inline complex& cf_h(int h, int g, int s, int o)
    { return _cf_h[((h*numAtoms+g)*2+s)*o_max+o];}
  inline complex& coef(int s, int o)
    { return _coef[s*o_max+o];}
  inline double& lattice(int iat, int dim)
    { return _lattice[iat*3+dim];}	  
  inline int& neighbor(int ia, int in)
    { return _nn[ia*4+in];}

private:
 qd_struct d; 
 int e_max, h_max, g_min, g_max, o_max, numAtoms, my_id;
 complex* _cf_e;
 complex* _cf_h;
 complex* _coef;
 int* _atom;
 int* _shape;
 double* _lattice;
 int* _nn;
 int mpi_numprocs;
};  

#endif
