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
$Header: /repo/nemo3d/src/base/h_cvectr_mult.h,v 1.11 2005/12/24 03:37:14 gekco Exp $ 
*****************************************************************************/


#ifndef H_CVECTR_MULT_H
#define H_CVECTR_MULT_H 1

#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include <rvector.h>
#include <i2tensor.h>
#include <i3tensor.h>
#include <r3tensor.h>
#include <c3tensor.h>

#include "realtype.h"
#include "cvector_lib.h"
#include "cmatrix_nonsym.h"
#include "constants_nemo.h"
#include "math3d.h"

#include "opt3d.h"
#include "mb_ham_spds_micro.h"

#include "MPI_Timing.h"

#include "qd_struct.h"
#include "i_3dout_files.h"
#include "matmul.h"
#include "mat_def.h"
#include "Hmult_stored.h"


#include "run3d_mpi.h"
#include "blas_run3d_f.h"
#include "mb_ham_spds_micro2.h"


// prototype for hamiltonian storage base class
class HamStore {
 public:
   int Ham_use_hermiticity;
   cmatrixmatrix Ham;
   ivectr Ham_neighbor;
   imatrix Ham_offset_x, Ham_offset_y, Ham_yc;
   imatrix Ham_transp;
   imatrix Ham_atom_nbr_2_m;
   imatrix Ham_atom_nbr_2_l;
   imatrix l_m_2_atom;

   HamStore() {;}
};



void matvect_mult2(const int** argList, complex* y, const complex* x);
void matvect_mult(const int** argList, complex* y, const complex* x);

void Hmult_spds_complete( cvectr y, qd_struct d, real shift, cvectr x );

int Hmult_spds_col_fullstored_par ( cvectr y, cvectr yc, qd_struct d, real s, cvectr x, int proc );
int Hmult_spds_col ( cvectr y, cvectr yc, qd_struct d, real s, cvectr x, int proc );
int Hmult_spds_col_prep_store ( cvectr y, cvectr yc, qd_struct d, real s, cvectr x, int proc );

int Hmult ( rvectr y, qd_struct d, real s, rvectr x );
int Hmult_c ( cvectr y, qd_struct d, real s, cvectr x );

void H_DiagStrainCorr_OffDiag_Wrapper(qd_struct d, int l, int m, int AtomType_l_m, int n, int cindx, 
			    int aindx, rvectr nnv, rvectr nnv0, real a_lattice[3]) ;


void H_Diag_Wrapper(qd_struct, int l, int m, int i, int j,int k, int AtomType_l_m, ivectr Basis);

void Ham_to_Hzb_20(qd_struct d);

int Randomize_H ( qd_struct d, int proc );

void compute_ritz_residual_norm_deviation(qd_struct d);

void compute_ritz_value       ( qd_struct d );
void compute_residual         ( qd_struct d );
void compute_eigenvector_norm ( qd_struct d );

complex vect_dot_vect (complex*, complex*, qd_struct d);

void count_nbrs_for_sparse_zincblende_ham(int* Nout, int* Nin, qd_struct d);

#endif



