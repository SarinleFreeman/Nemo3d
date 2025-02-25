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
$Header: /repo/nemo3d/src/base/phon_struct.h,v 1.6 2004/02/16 03:36:56 gekco Exp
$
*****************************************************************************/

#ifndef PHON_STRUCT_H
#define PHON_STRUCT_H

#include <nml_dcmatrix.h>
#include <nml_dmatrix.h>
#include <nml_imatrix.h>
#include <stdlib.h>

#include "ivectori3matrix.h"
#include "realtype.h"
#include <cmatrix.h>
#include <cvector.h>
#include <i2tensor.h>
#include <i3tensor.h>
#include <imatrix.h>
#include <r3tensor.h>
#include <rmatrix.h>
#include <rvector.h>
/* #include "mat_struct.h" */
#include "cmatrixmatrix.h"
#include "ivectorvector.h"
#include "rmatrixrvector.h"

#include <sstream>
// #include <ostream>
#include <iomanip>
#include <iostream>

using std::ostringstream;
using std::setw;

/**
 * Structure to define all computed phonon stuff.
 */
class phon_struct {

public:
  int states_in_row; //!< maximum number of nonzero elements in every row in
                     //!< Keating model

  // The storage of the Dynamical Matrix is almost similar to that by
  // Harwell-Boeing
  cmatrix DM;     //!< dynamical matrix
  imatrix Ind_DM; /*!< indexies of nonzero elements in m-th row of DM.
                       [i][0] element contains the number of nonzero elements
                       in i-th row*/

  // Map
  i3matrix DMmap; /*!< [total # of unit cells][# of atoms per unit cell][3 proj]
                   Map of indexies in the dynamical matrix*/
  // Phonon spectrum
  rmatrix E_ph; /*!< [N_branches+3][Nq_ph] matrix of phonon energy
                  dispersion branches
                  The first 3 columns contain
                  q_ph_x,q_ph_y,q_ph_z components*/
  int nbr_uc;   /*!< =1 for the standart Keating model with second-neares
                  neighbors interaction limit.   Number of neighbors unit cells in
                  x   direction, whose   atoms interact with the boundary atoms of
                  the   given unit cell*/
  // Parallel stuff:
  ivectr seg_col_s,
      seg_col_ln; /*!< [d->nproc] Initial index and the number of the
                    DM elements on every processor.
                    The partition is by columns.*/
  ivectr min_row,
      max_row; /*!< [d->nproc] minimal and maximal indexies of the rows of DM
                 to be stored on the mpi_n3d_id processor*/
  ivectr szCommL, szCommR; /*!< Number of rows in the buffer
                             to be exchanged with the left/right processor*/
  ivectr N_cols,
      N_rows; /*!< Number of columns/rows of DM stored on every processor*/
  ivectr mpi_left,
      mpi_right; /*!< Left/right processors mpi_n3d_id  for communication*/
  // Parts of the eigenvectors for communication
  cvectr wL; //!< part of the vector to be transfered to Left processor
  cvectr wR; //!< part of the vector to be transfered to Right processor
  cvectr vL; //!< part of the vector to be received from Left processor
  cvectr vR; //!< part of the vector to be received from Right processor

  void ph_s_alloc(int N_uc,       // number of unit cells
                  int N_a,        // number of atoms per unit cell
                  int N_branches, // number of phonon branches
                  int Nq_ph,      // number of points in the reciprocal space]
                  int nproc       // number of processors
  ); // Allocates memory for phonon-related matrixes

  void comm_vectors_allocate(); // Allocates memory for parts of the eigenvector
                                // that participate in communication

  void
  comm_vectors_deallocate(); // Deallocates memory for parts of the eigenvector
                             // that participate in communication

  void DM_allocate(
      int N_rows,          // total number of rows ON THE PROCESSOR
      int DM_states_in_row // max number of elements per row on mpi_n3d_id proc
  );                       /*Allocates memory for the Dynamical matrix and
                             buffers for exchange with neighbor processors*/
  void DM_reallocate(
      int N_rows,          // total number of rows ON THE PROCESSOR
      int DM_states_in_row // max number of elements per row on mpi_n3d_id proc
  );                       /*Reallocates memory for the Dynamical matrix and
                           buffers for exchange with neighbor processors*/
  void ph_s_dealloc();     // Deallocates memory for phonon-related matrixes
  void ph_s_clear();       // Clears ALL elements of phon_struct
  void Print_DM();         // Prints out the Dynamical matrix
  void Print_parameters(); // Prints all the parameters
  void
  Picture_DM(int DM_total_size); /*Picture the Dynamical matrix stored on the
                               processor "0" - zero, "X" -- nonzero element,
                               every element is followed by (mpi_n3d_id) and
                               pre... by <row,col> in the whole DM*/
};

#endif
