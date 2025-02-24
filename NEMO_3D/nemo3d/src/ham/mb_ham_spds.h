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
$Header: /repo/nemo3d/src/ham/mb_ham_spds.h,v 1.4 2003/10/08 16:17:22 hook Exp $
*****************************************************************************/

#ifndef MB_HAM_SPDS_H
#define MB_HAM_SPDS_H 1

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <rvector.h>
#include <complex.h>
#include <cmatrix.h>
#include "realtype.h"
#include "rmatrix_lib.h"
#include "mb_ham_spds_micro.h"


#define Make_sp3s_s_new_hamiltonian Make_spds_s_hamiltonian


/* These functions provide an interface to the spds* hamiltonian
   generator functions generated for 3D simulation.  These conform
   to the Nemo Hamiltonian generator format. */
int Make_spd_s_hamiltonian ( real *param, real *mb_strain, real Ev_offset, real kxa, real kya, cmatrix d, cmatrix s );
int Make_spds_s_hamiltonian( real *param, real *mb_strain, real Ev_offset, real kxa, real kya, cmatrix d, cmatrix s );
int Make_spdd_s_hamiltonian( real *param, real *mb_strain, real Ev_offset, real kxa, real kya, cmatrix d, cmatrix s );

void fold_t_k_sp3s_s_new ( real kxa, real kya, cmatrix d, cmatrix s, cmatrix d_fold, cmatrix s_fold );

void fold_t_k_spd_s ( real kxa, real kya, cmatrix d, cmatrix s, cmatrix d_fold, cmatrix s_fold );
void fold_t_k_spds_s ( real kxa, real kya, cmatrix d, cmatrix s, cmatrix d_fold, cmatrix s_fold );
void fold_t_k_spdd_s ( real kxa, real kya, cmatrix d, cmatrix s, cmatrix d_fold, cmatrix s_fold );
real mb_energy_spds( rvectr p, void *data1, void *data2, void *data3 );
real mb_mstar_spds( rvectr p, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, int xdir, int ydir, int zdir, rmatrix mtens );
real mb_mstar_spds_111( rvectr p, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length, rmatrix mtens );
real mb_newton_spds( rvectr p, rvectr param, rvectr mb_strain, real unstrnd_cubic_cell_length );
#endif /* MB_HAM_SPDS_H */

