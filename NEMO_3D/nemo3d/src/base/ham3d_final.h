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
$Header: /repo/nemo3d/src/base/ham3d_final.h,v 1.8 2003/11/20 23:10:58 swlee Exp $ 
*****************************************************************************/

#ifndef HAM3D_FINAL_H
#define HAM3D_FINAL_H 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "io_utils.h"

#include "realtype.h"
#include "constants_nemo.h"
#include <rvector.h>
#include "i2tensor.h"
#include "i3tensor.h"
#include "cmatrix_nonsym.h"
#include "qd_struct.h"
#include "run3d_mpi.h"
#include "strain3d.h"
#include "mat3d.h"
#include "ivectorvector.h"

#include "MPI_Timing.h"

void basis_init(qd_struct d );
void nbr_init(qd_struct d);

void adjust_topology(qd_struct d);
int vecmap_init ( qd_struct d, int proc, ivectr vecmap_proc );
int mpi_segment_init( qd_struct d, int nproc );
int get_total_hamiltonian_dimension_and_offsets(qd_struct d);
void adjust_crystal_for_electronic_calculation(qd_struct d);

void setup_crystal_structure(qd_struct d);
void setup_hamiltonian(qd_struct d);
void setup_hamiltonian_for_postprocessing(qd_struct d);
void setup_crystal_structure_for_postprocessing(qd_struct d);

void ham_file_init_spds(qd_struct d, geo_struct* di);

#endif
