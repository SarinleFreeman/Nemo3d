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
$Header: /repo/nemo3d/src/base/phon_Keating.h,v 1.4 2004/02/27 04:09:05 lol Exp
$
*****************************************************************************/

#ifndef phon_Keating_h
#define phon_Keating_h

#include "Boolean.h"
#include "cmacopt.h"
#include "constants_nemo.h"
#include "i_3dout_files.h"
#include "io_utils.h"
#include "macopt.h"
#include "macopt_par.h"
#include "qd_struct.h"
#include "realtype.h"
#include "run3d_mpi.h"
#include "util_memory.h"
#include "util_stringops.h"
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "MPI_Timing.h"

// #include <iostream.h>
#include <iostream>

#include "phon_arpack.h"
#include "phon_output.h"
#include "phon_parpack.h"
#include "phon_struct.h"

int N_DM(qd_struct d,
         int l, //-th unit cell
         int m, //-th atom
         int xi //-th component
);
// Returnes index of the element in the Dynamical Matrix

void Map_DM(qd_struct d);
/*Mapping of the indexies of the dynamical matrix in
  d->phon.DMmap[d->geo.N_Cell][d->geo.cellgeom.AtomsPerCellMax()]*/

void initial_assignments(
    qd_struct d,
    int n_qph,     // point in the reciprocal space
    int l,         // cell
    int m,         // atom
    cvectr phase,  //[n_max] phase associated with n-th neighbor
    rvectr alpha,  //[n_max] bond length VFF constant associated with n-th
                   // neighbor of m-th atom
    rvectr beta,   //[n_max] bond bending VFF constant associated with n-th
                   // neighbor of m-th atom
    rvectr sqrtMn, //[n_max] square root of the atomic mass associated with n-th
                   // neighbor
    rmatrix
        nnv, //[n_max] strained bond length between n-th neighbor and m-th atom
    rmatrix
        nnv0, //[n_max] realxed bond length between n-th neighbor and m-th atom
    rvectr dssq, //[n_max] strained bond length squared between n-th neighbor
                 // and m-th atom
    rvectr d0sq  //[n_max] realxed bond length squared between n-th neighbor and
                 // m-th atom
);
// Assign a set of arrays associated to the n-th neighbor of m-th atom in l-th
// unit cell n_max=d->geo.Neighbors(m)

void Add_element(qd_struct d, int n_row, int n_col, real tail_r, real tail_i);
/*Adds tail to the Dynamical matrix d->phon.DM[row][i-1] and
  the corresponding column indexes to d->phon.Ind_DM[row][i].
  Here 0<i<=d->phon.states_in_row */

void Add(qd_struct d, int row, int col, real tail_r, real tail_i);
/*Decide either add the element to DM on the current processor,
  what DM indexies n_row and n_col should be
  and call Add_element*/

void Add_Herm(qd_struct d, int row, int col, real tail_r, real tail_i);
/*Add tail and tail* to "hermitian" parts of DM*/

real DM_A1(real alpha, // bond length distorsion force constant array of
                       // N_Neighbor[m] elements
           real dssq,  // square of the length of the strained bond
           real d0sq,  // square of the length of the relaxed bond
           real nnv    // component of the strained bond
);
/*Computes A1:
  A1=d^2 E_alpha/du_x(n)du_x(n)=d^2 E_alpha/du_x(m)du_x(m)
  -A1=d^2 E_alpha/du_x(n)du_x(m)=d^2 E_alpha/du_x(m)du_x(n)
*/

real DM_A2(real alpha, // bond length distorsion force constant array of
                       // N_Neighbor[m] elements
           real d0sq,  // square of the length of the relaxed bond
           real x_ij, real y_ij // components of the strained bond
);
/*Computes A2:
  A2=d^2 E_alpha/du_x(n)du_y(n)=d^2 E_alpha/du_x(m)du_y(m)
  -A2=d^2 E_alpha/du_x(n)du_y(m)=d^2 E_alpha/du_x(m)du_y(n)
*/

real doublet(real beta1, real beta2, // bond bending distorsion force constant
             real d0sq1, real d0sq2  // square of the length of the relaxed bond
);
/*Computes coefficient
  0.375*sqrt(beta1*beta2)/sqrt(d0sq1*d0sq2)
  for the given V-shaped bond doublet*/

real DM_D(rvectr d1, rvectr d2, // vectors of the relaxed bonds of the V doublet
          rvectr r1, rvectr r2 // vectors of the strained bonds of the V doublet
);
/*Computes D:
  Term for bond bending distorsion
  to be added in different elements of DM
  with the same displacement projection
*/

void Add_A1(qd_struct d,
            real A1,      // element to be added
            int mDM,      // row
            int nDM,      // column
            real sqrtMm,  // square root of the atomic mass of m-th atom
            real sqrtMn,  // square root of the atomic mass ofthe n-th neighbor
                          // of m-th atom
            real phase_r, // real part of the phase
            real phase_i  // imaginary part of the phase
);
/*Add elements A1 to Dynamical Matrix*/

void Add_A2(
    qd_struct d,
    real A2,     // element to be added
    int l,       // cell
    int m,       // atom
    int xim,     // displacement component of the m-th atom
    int cindx,   // cell of the n-th neighbor of the m-th atom
    int aindx,   // atom of the n-th neighbor of the m-th atom
    int xin,     // displacement component of the n-th neighbor of the m-th atom
    real sqrtMm, // square root of the atomic mass of m-th atom
    real sqrtMn, // square root of the atomic mass ofthe n-th neighbor of m-th
                 // atom
    real phase_r, // real part of the phase
    real phase_i  // imaginary part of the phase
);
/*Add elements A2 to Dynamical Matrix*/

void Add_BXX(qd_struct d,
             real V,             // factor for the V doublet
             real D,             //=(r_mn1*r_mn2-d_mn1*d_mn2)
             real xn1, real xn2, // components of the bond length between m-th
                                 // atom and n1(n2) neighbor
             int cell_m, int atom_m, // cell and atom indexes of the vertex atom
             int cell_n1,
             int atom_n1, // cell and atom indexes of the n1-th neighbor atom
             int cell_n2,
             int atom_n2,  // cell and atom indexes of the n2-th neighbor atom
             int xi,       // projection of the displacement
             real sqrtMm,  // square root of the atomic mass of m-th atom
             real sqrtMn1, // square root of the atomic mass of the n1-th
                           // neighbor of m-th atom
             real sqrtMn2, // square root of the atomic mass of the n2-th
                           // neighbor of m-th atom
             real phase_n1_r, // real part of the phase of the n1-th neighbor of
                              // m-th atom
             real phase_n1_i, // imaginary part of the phase of the n1-th
                              // neighbor of m-th atom
             real phase_n2_r, // real part of the phase of the n2-th neighbor of
                              // m-th atom
             real phase_n2_i  // imaginary part of the phase of the n2-th
                              // neighbor of m-th atom
);
/*Add elements B to Dynamical Matrix
  when projections of the displacements are the same*/

void Add_BX1X2(
    qd_struct d,
    real V, // factor for the V doublet
    real x1n1, real x1n2, real x2n1,
    real x2n2, // components of the bond length between m-th atom and n1(n2)
               // neighbor
    int cell_m, int atom_m, // cell and atom indexes of the vertex atom
    int cell_n1,
    int atom_n1, // cell and atom indexes of the n1-th neighbor atom
    int cell_n2,
    int atom_n2,      // cell and atom indexes of the n2-th neighbor atom
    int xi1, int xi2, // projections of the displacement
    real sqrtMm,      // square root of the atomic mass of m-th atom
    real sqrtMn1,     // square root of the atomic mass of the n1-th neighbor of
                      // m-th atom
    real sqrtMn2,     // square root of the atomic mass of the n2-th neighbor of
                      // m-th atom
    real
        phase_n1_r, // real part of the phase of the n1-th neighbor of m-th atom
    real phase_n1_i, // imaginary part of the phase of the n1-th neighbor of
                     // m-th atom
    real
        phase_n2_r, // real part of the phase of the n2-th neighbor of m-th atom
    real phase_n2_i // imaginary part of the phase of the n2-th neighbor of m-th
                    // atom
);
/*Add elements B to Dynamical Matrix
  when projections of the displacements are different*/

void comb_DM(
    ivectr index, //[d->phon.Ind_DM[n_row][i]+1] vector d->phon.Ind_DM[n_row][i]
    cvectr DMrow, //[d->phon.Ind_DM[n_row][i]] vector d->phon.DM[n_row][i-1]
    int beg, int end // lower and upper bounds of the arrays
);
/*Modification of quick_sort to
  sort d->phon.Ind_DM[n_row][i]
  in ascending order.
  Also complex elements of DM
  d->phon.DM[n_row][i-1] are swapped.*/

void DM_zeros(qd_struct d);
// Assignes zeros to elements with value less then zero

void DM_Keating(qd_struct d,
                int n_qph, //-th point in the reciprocal space
                int cell0, // initial cell index
                int cellF  // final cell index
);
// Generates valence-field-force part of dynamical matrix

void phon_disp(qd_struct d);
/*Computes phonon dispersion for
  d->Nq_ph points in the reciprocal space starting from
  d->q_ph_min to d->q_ph_min */

////////////////////Printouts are in phon_output.[ch]//////////////////

#endif
