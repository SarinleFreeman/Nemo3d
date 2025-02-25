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
$Header: /repo/nemo3d/src/base/Hmult_stored.h,v 1.7 2007/07/24 17:47:36 lee509
Exp $
*****************************************************************************/

#ifndef HMULT_STORED_H
#define HMULT_STORED_H

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "MPI_Timing.h"
#include "constants_nemo.h"
#include "qd_struct.h"
#include "realtype.h"

/*! \file Hmult_stored.h
 *  \brief Matrix-vector / matrix-matrix multiplier for sparse storage mode.
 */

int Hmult_spds_col_stored(cvectr y, cvectr yc, qd_struct d, real s, cvectr x,
                          int proc);
/*!< Matrix-vector multiplication routine for sparse storage mode y =
 * (H-shift*I)*x
 * @param y OUT: Output vector
 * @param yc OUT: Output vector for neighboring processor. It is sent to
 * neighbor processors via MPI_Sendrecv
 * @param d IN: QD_struct, see qd_struct.h
 * @param s IN: Energy shift
 * @param x IN: Input vector
 * @param proc IN: Processor number
 * @return length of vector that needs to be passed to neighbor processor
 */

int Hmatmult_spds_col_stored(cmatrix Y, cmatrix YC, qd_struct d, real s,
                             cmatrix X, int num_col, int proc);
/*!< Matrix-matrix multiplication routine for sparse storage mode Y =
 * (H-shift*I)*X
 * @param Y OUT: Output matrix
 * @param YC OUT: Output matrix for neighboring processor. It is sent to
 * neighbor processors via MPI_Sendrecv
 * @param d IN: QD_struct, see qd_struct.h
 * @param s IN: Energy shift
 * @param X IN: Input matrix
 * @param num_col IN: Number of columns of the matrix
 * @param proc IN: Processor number
 * @return length of vector that needs to be passed to neighbor processor
 */

void Hv_stored_20_Hzb(cvectr, cvectr, cvectr, qd_struct);
/*! Matrix-vector multiplication routine for sp3d5s* with spin orbital model.
 * Memory efficient mode. See HamZincBlende.c
 * @param y OUT: Output vector
 * @param yc OUT: Output vector for neighbor processor
 * @param x IN: Input vector
 * @param d IN: QD_struct, see qd_struct.h
 * @return void
 */

void Hv_stored_20_HzbComplex(cvectr, cvectr, cvectr, qd_struct);
/*! Matrix-vector multiplication routine for sp3d5s* with spin orbital model.
 * Speed efficient mode. See HamZincBlendeComplex.c
 * @param y OUT: Output vector
 * @param yc OUT: Output vector for neighbor processor
 * @param x IN: Input vector
 * @param d IN: QD_struct, see qd_struct.h
 * @return void
 */

void Hv_stored_10_sp3d5s_Hzb(cvectr, cvectr, cvectr, qd_struct);
/*! Matrix-vector multiplication routine for sp3d5s* without spin orbital model.
 * Memory efficientt mode. See HamZincBlende.c
 * @param y OUT: Output vector
 * @param yc OUT: Output vector for neighbor processor
 * @param x IN: Input vector
 * @param d IN: QD_struct, see qd_struct.h
 * @return void
 */

void Hv_stored_10_sp3d5s_HzbComplex(cvectr, cvectr, cvectr, qd_struct);
/*! Matrix-vector multiplication routine for sp3d5s* without spin orbital model.
 * Speed efficient mode. See HamZincBlendeComplex.c
 * @param y OUT: Output vector
 * @param yc OUT: Output vector for neighbor processor
 * @param x IN: Input vector
 * @param d IN: QD_struct, see qd_struct.h
 * @return void
 */

void Hv_stored_10_Hzb(cvectr, cvectr, cvectr, qd_struct);
/*! Matrix-vector multiplication routine for sp3s* with spin orbital model.
 * Memory efficient mode. See HamZincBlende.c
 * @param y OUT: Output vector
 * @param yc OUT: Output vector for neighbor processor
 * @param x IN: Input vector
 * @param d IN: QD_struct, see qd_struct.h
 * @return void
 */

void Hv_stored_10_HzbComplex(cvectr, cvectr, cvectr, qd_struct);
/*! Matrix-vector multiplication routine for sp3s* with spin orbital model.
 * Speed efficient mode. See HamZincBlendeComplex.c
 * @param y OUT: Output vector
 * @param yc OUT: Output vector for neighbor processor
 * @param x IN: Input vector
 * @param d IN: QD_struct, see qd_struct.h
 * @return void
 */

#endif
