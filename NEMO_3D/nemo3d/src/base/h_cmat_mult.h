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
$Header: /repo/nemo3d/src/base/h_cmat_mult.h,v 1.3 2007/07/24 16:53:15 lee509 Exp $ 
*****************************************************************************/


#ifndef H_CMAT_MULT_H
#define H_CMAT_MULT_H 

/*! \file h_cmat_mult.h
 *  \brief These routines support Matrix-Matrix Multiplication of Hamiltonian H
 *  Y = (H-sigma*I)*X. For external use as an API, use the function Hmatmult_spds_complete
 */

void matmat_mult(const int** argList, nml_dcmatrix* y, const nml_dcmatrix* x);
/*!< Provides external function pointer for Matrix-Matrix multiplication routine.
 * @param argList IN: Argument lists needed for matrix-matrix multiplication. 
 *                arg[0]: qd_struct, arg[1]: shift, arg[2]: number of columns for matrix multiplication, arg[3]: Nproj (Not used)
 * @param y OUT: Output complex matrix
 * @param x IN: Input complex matrix
 * @return void
 */

void Hmatmult_spds_complete( cmatrix Y, qd_struct d, real shift, cmatrix X, int num_col );
/*!< API for Matrix-Matrix multiplication. Y = (H-shift*I)*X
 * @param Y OUT: Output complex matrix
 * @param d IN: QD_struct, see qd_struct.h
 * @param shift IN: Shift value in energy
 * @param X IN: Input complex matrix
 * @param num_col IN: Number of columns of input matrix X
 * @return void
 */

int Hmatmult_spds_col_fullstored_par( cmatrix Y, cmatrix YC, qd_struct d, real s, cmatrix X, int num_col, int proc );
/*!< Matrix-Matrix multiplication routine for full storage mode. Y = (H-shift*I)*X
 * @param Y OUT: Output complex matrix
 * @param YC OUT: Output complex matrix for neighboring processor
 * @param d IN: QD_struct, see qd_struct.h
 * @param s IN: Shift value in energy
 * @param X IN: Input complex matrix
 * @param num_col IN: Number of columns of input matrix X
 * @param proc IN: Processor number
 * @return length of vector that needs to be passed to neighboring processor
 */

int Hmatmult_spds_col( cmatrix Y, cmatrix YC, qd_struct d, real s, cmatrix X, int num_col, int proc );
/*!< Matrix-Matrix multiplication routine for full storage mode. Y = (H-shift*I)*X
 * @param Y OUT: Output complex matrix
 * @param YC OUT: Output complex matrix for neighboring processor
 * @param d IN: QD_struct, see qd_struct.h
 * @param s IN: Shift value in energy
 * @param X IN: Input complex matrix
 * @param num_col IN: Number of columns of input matrix X
 * @param proc IN: Processor number
 * @return length of vector that needs to be passed to neighboring processor
 */

int Hmatmult_spds_col_prep_store( cmatrix Y, cmatrix YC, qd_struct d, real s, cmatrix X, int num_col, int proc );
/*!< Construction function for Hamiltonian (Y/YC/s/X/num_col are not used)
 * @param d IN: QD_struct, see qd_struct.h
 * @param proc IN: Processor number
 * @return 0
 */

#endif

